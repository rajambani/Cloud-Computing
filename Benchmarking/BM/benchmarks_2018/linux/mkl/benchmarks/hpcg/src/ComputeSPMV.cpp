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
 @file ComputeSPMV.cpp

 HPCG routine
 */

#include "ComputeSPMV.hpp"
#include "ComputeSPMV_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"

#include <mpi.h>
#include "Geometry.hpp"
#include <cstdlib>
#endif
/*!
  Routine to compute sparse matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This routine calls the reference SpMV implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV_ref
*/

int ComputeSPMV( const SparseMatrix & A, Vector & x, Vector & y)
{
    sparse_status_t status = SPARSE_STATUS_SUCCESS;
    struct optData *optData = (struct optData *)A.optimizationData;
    struct matrix_descr descr;
    sparse_matrix_t csrA = (sparse_matrix_t)optData->csrA;
    sparse_matrix_t csrB = (sparse_matrix_t)optData->csrB;
    descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descr.mode = SPARSE_FILL_MODE_FULL;
    descr.diag = SPARSE_DIAG_NON_UNIT;

    #ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
    #endif

    status = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descr, x.values, 0.0, y.values );

    if ( A.geom->size > 1 )
    {
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        status = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrB, descr, x.values, 0.0, optData->dtmp );
        #ifndef HPCG_NO_OPENMP
        #pragma omp parallel for
        #endif
        #pragma ivdep
        for (local_int_t i=0; i<optData->nrow_b; i++)
        {
            const MKL_INT ind = optData->bmap[i];
            y.values[ind] += optData->dtmp[i];
        }
    }

    return 0;
}

int ComputeSPMV_DOT( const SparseMatrix & A, Vector & x, Vector & y, double & pAp)
{
    sparse_status_t status = SPARSE_STATUS_SUCCESS;
    struct optData *optData = (struct optData *)A.optimizationData;
    struct matrix_descr descr;
    sparse_matrix_t csrA = (sparse_matrix_t)optData->csrA;
    sparse_matrix_t csrB = (sparse_matrix_t)optData->csrB;
    descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descr.mode = SPARSE_FILL_MODE_FULL;
    descr.diag = SPARSE_DIAG_NON_UNIT;

    #ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
    #endif

    status = mkl_sparse_d_dotmv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descr, x.values, 0.0, y.values, &pAp );

    if ( A.geom->size > 1 )
    {
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        status = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrB, descr, x.values, 0.0, optData->dtmp );
        double pAp_loc = 0, pAp_red = 0;
        #ifndef HPCG_NO_OPENMP
        #pragma omp parallel for reduction(+:pAp_loc,pAp_red)
        #endif
        #pragma ivdep
        for (local_int_t i=0; i<optData->nrow_b; i++)
        {
            const MKL_INT ind = optData->bmap[i];
            const double x_cur = x.values[ind];
            pAp_red += y.values[ind]*x_cur;
            y.values[ind] += optData->dtmp[i];
            pAp_loc += y.values[ind]*x_cur;
        }
        pAp += (pAp_loc - pAp_red);
    }

    return 0;
}
