/*******************************************************************************
* Copyright 2014-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file SetupHalo.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include <mpi.h>
#include <map>
#include <set>
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "SetupHalo.hpp"
#include "SetupHalo_ref.hpp"

/*!
  Prepares system matrix data structure and creates data necessary necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/

void SetupHalo(SparseMatrix & A)
{
    if( A.geom->size > 1 )
    {
        local_int_t localNumberOfRows = A.localNumberOfRows;
        char  * nonzerosInRow = A.nonzerosInRow;
        global_int_t ** mtxIndG = A.mtxIndG;
        local_int_t ** mtxIndL = A.mtxIndL;
        double t1 = 0;

        int nproc = A.nproc;
        nproc = (nproc > A.numOfBoundaryRows) ? A.numOfBoundaryRows : nproc;

        local_int_t totalToBeSent = 0, totalToBeReceived = 0, cnt = 0;
        local_int_t number_of_neighbors = A.numberOfSendNeighbors;
        local_int_t *receiveLength = NULL, *sendLength = NULL, *map_neib_s = NULL;
        local_int_t *map_neib_r = A.work, *neighbors = NULL;
        global_int_t *map_send = NULL;

        double *sendBuffer = NULL;
        local_int_t  *elementsToSend   = NULL, *all_send = NULL, *all_recv = NULL;
        global_int_t *elementsToSend_G = NULL, *elementsToRecv_G = NULL;
        std::map< local_int_t, local_int_t > externalToLocalMap;
        local_int_t receiveEntryCount = 0, sendEntryCount = 0;

        map_send      = (global_int_t*) MKL_malloc( sizeof(global_int_t)*A.numOfBoundaryRows*number_of_neighbors, 512 );
        map_neib_s    = (local_int_t *) MKL_malloc( sizeof(local_int_t )*A.geom->size                           , 512 );
        neighbors     = (local_int_t *) MKL_malloc( sizeof(local_int_t )*number_of_neighbors                    , 512 );
        receiveLength = (local_int_t *) MKL_malloc( sizeof(local_int_t )*number_of_neighbors                    , 512 );
        sendLength    = (local_int_t *) MKL_malloc( sizeof(local_int_t )*number_of_neighbors                    , 512 );

        if ( map_neib_s == NULL || map_send == NULL || neighbors == NULL || receiveLength == NULL || sendLength == NULL ) return;

        for ( local_int_t i = 0; i < A.geom->size; i ++ ) map_neib_s[i] = 0;

        for ( local_int_t i = 0; i < A.geom->size; i ++ )
        {
            if ( map_neib_r[i] > 0 )
            {
                neighbors[cnt] = i;
                map_neib_s[i] = cnt++;
            }
        }
#ifndef HPCG_NO_OPENMP
        #pragma omp parallel for num_threads(nproc)
#endif
        for ( local_int_t i = 0; i < A.numOfBoundaryRows*number_of_neighbors;    i ++ ) map_send[i] = -1;

        for ( local_int_t i = 0; i < number_of_neighbors; i ++)
        {
            receiveLength[i] = 0;
            sendLength[i] = 0;
        }
#ifndef HPCG_NO_OPENMP        
        #pragma omp parallel for reduction(+:totalToBeSent)
#endif
        for ( local_int_t row = 0; row < A.numOfBoundaryRows; row++ )
        {
            local_int_t i = A.boundaryRows[row];
            for ( local_int_t j = 0; j < nonzerosInRow[i]; j++ )
            {
                if( mtxIndL[i][j] < 0 )
                {
                    global_int_t curIndex = mtxIndG[i][j];
                    int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
                    local_int_t jj = map_neib_s[rankIdOfColumnEntry];

                    if ( map_send[jj*A.numOfBoundaryRows + row] < 0 )
                    {
                        map_send[jj*A.numOfBoundaryRows + row] = A.localToGlobalMap[i];
                        totalToBeSent++;
                    }
                }
            }
        }
        sendBuffer       = (double      *) MKL_malloc(sizeof(double      )*totalToBeSent, 512);
        elementsToSend   = (local_int_t *) MKL_malloc(sizeof(local_int_t )*totalToBeSent, 512);
        elementsToSend_G = (global_int_t*) MKL_malloc(sizeof(global_int_t)*totalToBeSent, 512);

        if ( sendBuffer == NULL || elementsToSend == NULL || elementsToSend_G == NULL ) return;

        for ( local_int_t i = 0; i < number_of_neighbors; i ++ )
        {
            for ( local_int_t j = 0; j < A.numOfBoundaryRows; j ++)
            {
                global_int_t ind = map_send[i*A.numOfBoundaryRows + j];
                if ( ind >= 0 )
                {
                    elementsToSend_G[sendEntryCount] = ind;
                    sendLength[i] ++;
                    sendEntryCount++;
                }
            }
        }

        all_send = (local_int_t *) MKL_malloc(sizeof(local_int_t )*A.geom->size, 512);
        all_recv = (local_int_t *) MKL_malloc(sizeof(local_int_t )*A.geom->size, 512);

        if ( all_send == NULL || all_recv == NULL ) return;

        for( local_int_t i = 0 ; i < A.geom->size; i++ ) { all_send[i] = 0; all_recv[i] = 0; }
        for (int i = 0; i < number_of_neighbors; i++)
        {
            all_send[ neighbors[i] ] = sendLength[i];
        }
        MPI_Alltoall( all_send, 1, MPI_INT, all_recv, 1, MPI_INT, MPI_COMM_WORLD );
        local_int_t totalToBeRecv = 0;
        for (int i = 0; i < number_of_neighbors; i++)
        {
            receiveLength[ i ] = all_recv[ neighbors[i] ];
            totalToBeRecv += receiveLength[ i ];
        }
        elementsToRecv_G = (global_int_t*) MKL_malloc(sizeof(global_int_t)*totalToBeRecv, 512);

        if ( elementsToRecv_G == NULL ) return;

        int MPI_MY_TAG = 98;

        MPI_Request *request = (MPI_Request*) MKL_malloc(sizeof(MPI_Request)*number_of_neighbors, 512);

        if (request == NULL ) return;

        global_int_t *elementsToRecv_pt = elementsToRecv_G;

        for (int i = 0; i < number_of_neighbors; i++) {
            local_int_t n_recv = receiveLength[i];
            MPI_Irecv(elementsToRecv_pt, n_recv, MPI_LONG_LONG_INT, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
            elementsToRecv_pt += n_recv;
        }
        global_int_t *elementsToSend_pt = elementsToSend_G;
        for (int i = 0; i < number_of_neighbors; i++) {
            local_int_t n_send = sendLength[i];
            MPI_Send(elementsToSend_pt, n_send, MPI_LONG_LONG_INT, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
            elementsToSend_pt += n_send;
        }
        MPI_Status status;
        for (int i = 0; i < number_of_neighbors; i++) {
            if ( MPI_Wait(request+i, &status) ) {
                std::exit(-1); // TODO: have better error exit
            }
        }
#ifndef HPCG_NO_OPENMP
        #pragma omp parallel for
#endif
        for( local_int_t i = 0; i < totalToBeSent; i++ )
        {
            elementsToSend[ i ] = A.globalToLocalMap[ elementsToSend_G[ i ] ];
        }

        receiveEntryCount = 0;
        for (local_int_t i = 0; i < number_of_neighbors; i++)
        {
            for( local_int_t j = 0; j < receiveLength[i]; j ++, receiveEntryCount++ ) 
            {
                externalToLocalMap[elementsToRecv_G[receiveEntryCount]] = localNumberOfRows + receiveEntryCount; 
            }
        }
#ifndef HPCG_NO_OPENMP
        #pragma omp parallel for
#endif
        for ( local_int_t row = 0; row < A.numOfBoundaryRows; row++ )
        {
            local_int_t i = A.boundaryRows[row];
            for ( local_int_t j = 0; j < nonzerosInRow[i]; j++ )
            {
                if ( mtxIndL[i][j] < 0 ) // My column index, so convert to local index
                {
                    global_int_t curIndex = mtxIndG[i][j];
                    mtxIndL[i][j] = externalToLocalMap[curIndex];
                }
            }
        }

        local_int_t *sdispls = (local_int_t*) MKL_malloc( sizeof(local_int_t)*A.geom->size, 128 );
        local_int_t *rdispls = (local_int_t*) MKL_malloc( sizeof(local_int_t)*A.geom->size, 128 );
        local_int_t *scounts = (local_int_t*) MKL_malloc( sizeof(local_int_t)*A.geom->size, 128 );
        local_int_t *rcounts = (local_int_t*) MKL_malloc( sizeof(local_int_t)*A.geom->size, 128 );
        local_int_t tmp_s = 0, tmp_r = 0;

        if(sdispls == NULL || rdispls == NULL || scounts == NULL || rcounts == NULL) return;

        for( local_int_t i = 0; i < A.geom->size; i++ )
        {
            scounts[i] = 0;
            rcounts[i] = 0;
            sdispls[i] = 0;
            rdispls[i] = 0;
        }
        for( local_int_t i = 0; i < number_of_neighbors; i++ )
        {
            local_int_t root = neighbors[i];
            scounts[root] = sendLength[i];
            rcounts[root] = receiveLength[i];
            sdispls[root] = tmp_s; tmp_s+=sendLength[i];
            rdispls[root] = tmp_r; tmp_r+=receiveLength[i];
        }
        A.scounts = scounts;
        A.rcounts = rcounts;
        A.sdispls = sdispls;
        A.rdispls = rdispls;

        A.numberOfExternalValues = externalToLocalMap.size();
        A.localNumberOfColumns = A.localNumberOfRows + A.numberOfExternalValues;
        A.numberOfSendNeighbors = number_of_neighbors;
        A.totalToBeSent = totalToBeSent;
        A.elementsToSend = elementsToSend;
        A.neighbors = neighbors;
        A.receiveLength = receiveLength;
        A.sendLength = sendLength;
        A.sendBuffer = sendBuffer;

        MKL_free(map_send);
        MKL_free(map_neib_s);
        MKL_free(elementsToSend_G);
        MKL_free(all_send);
        MKL_free(all_recv);
        MKL_free(elementsToRecv_G);
        MKL_free(request);
    } else {
        A.numberOfExternalValues = 0;
        A.localNumberOfColumns = A.localNumberOfRows;
        A.numberOfSendNeighbors = 0;
        A.totalToBeSent = 0;
        A.elementsToSend = 0;
    }
}
