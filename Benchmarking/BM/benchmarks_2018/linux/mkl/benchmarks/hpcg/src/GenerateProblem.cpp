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
 @file GenerateProblem.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
#include "GenerateProblem.hpp"
#include "GenerateProblem_ref.hpp"


/*!
  Routine to generate a sparse matrix, right hand side, initial guess, and exact solution.

  @param[in]  A        The generated system matrix
  @param[inout] b      The newly allocated and generated right hand side vector (if b!=0 on entry)
  @param[inout] x      The newly allocated solution vector with entries set to 0.0 (if x!=0 on entry)
  @param[inout] xexact The newly allocated solution vector with entries set to the exact solution (if the xexact!=0 non-zero on entry)

  @see GenerateGeometry
*/

void GenerateProblem(SparseMatrix & A, Vector * b, Vector * x, Vector * xexact) {

  // The call to this reference version of GenerateProblem can be replaced with custom code.
  // However, the data structures must remain unchanged such that the CheckProblem function is satisfied.
  // Furthermore, any code must work for general unstructured sparse matrices.  Special knowledge about the
  // specific nature of the sparsity pattern may not be explicitly used.
  // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  // below may result in global range values.
#if 0
  GenerateProblem_ref(A,b,x,xexact);
#else
  // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  // below may result in global range values.
//  double t_1 = dsecnd();
  global_int_t nx = A.geom->nx;
  global_int_t ny = A.geom->ny;
  global_int_t nz = A.geom->nz;
  global_int_t npx = A.geom->npx;
  global_int_t npy = A.geom->npy;
  global_int_t npz = A.geom->npz;
  global_int_t ipx = A.geom->ipx;
  global_int_t ipy = A.geom->ipy;
  global_int_t ipz = A.geom->ipz;
  global_int_t gnx = nx*npx;
  global_int_t gny = ny*npy;
  global_int_t gnz = nz*npz;
  double omp1,omp2;

  local_int_t localNumberOfRows = nx*ny*nz; // This is the size of our subblock
  // If this assert fails, it most likely means that the local_int_t is set to int and should be set to long long
  assert(localNumberOfRows>0); // Throw an exception of the number of rows is less than zero (can happen if int overflow)
  local_int_t numberOfNonzerosPerRow = 27; // We are approximating a 27-point finite element/volume/difference 3D stencil

  global_int_t totalNumberOfRows = ((global_int_t) localNumberOfRows)*((global_int_t) A.geom->size); // Total number of grid points in mesh
  // If this assert fails, it most likely means that the global_int_t is set to int and should be set to long long
  assert(totalNumberOfRows>0); // Throw an exception of the number of rows is less than zero (can happen if int overflow)
//  double t1 = dsecnd();
  // Allocate arrays that are of length localNumberOfRows
  char * nonzerosInRow = new char[localNumberOfRows];
  global_int_t ** mtxIndG = new global_int_t*[localNumberOfRows];
  local_int_t  ** mtxIndL = new local_int_t*[localNumberOfRows];
  double ** matrixValues = new double*[localNumberOfRows];
  double ** matrixDiagonal = new double*[localNumberOfRows];
  if (b!=0) InitializeVector(*b, localNumberOfRows);
  if (x!=0) InitializeVector(*x, localNumberOfRows);
  if (xexact!=0) InitializeVector(*xexact, localNumberOfRows);
  double * bv = 0;
  double * xv = 0;
  double * xexactv = 0;
  if (b!=0) bv = b->values; // Only compute exact solution if requested
  if (x!=0) xv = x->values; // Only compute exact solution if requested
  if (xexact!=0) xexactv = xexact->values; // Only compute exact solution if requested
  int nproc = A.nproc;
  A.localToGlobalMap.resize(localNumberOfRows);
  // Now allocate the arrays pointed to
  A.mtxL = (local_int_t*) MKL_malloc(sizeof(local_int_t )*numberOfNonzerosPerRow*localNumberOfRows, ALIGN );
  A.mtxG = (global_int_t*)MKL_malloc(sizeof(global_int_t)*numberOfNonzerosPerRow*localNumberOfRows, ALIGN );
  A.mtxA = (double*)      MKL_malloc(sizeof(double      )*numberOfNonzerosPerRow*localNumberOfRows, ALIGN );

  A.boundaryRows = new local_int_t[nx*ny*nz - (nx-2)*(ny-2)*(nz-2)];

  if ( A.mtxL == NULL || A.mtxG == NULL || A.mtxA == NULL || nonzerosInRow == NULL || A.boundaryRows == NULL
       || mtxIndG == NULL || mtxIndL == NULL || matrixValues == NULL || matrixDiagonal == NULL )
  {
      return;
  }

  A.numOfBoundaryRows = 0;
  local_int_t numOfBoundaryRows = 0;
#ifndef HPCG_NO_OPENMP
#pragma omp parallel reduction(+:numOfBoundaryRows) num_threads(nproc)
#endif
{
#ifndef HPCG_NO_OPENMP
 #pragma omp for nowait
#endif
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      A.boundaryRows[y*nx + x] = y*nx + x;
      numOfBoundaryRows++;
    }
  }
#ifndef HPCG_NO_OPENMP
  #pragma omp for nowait
#endif
  for (int z = 1; z < nz - 1; z++) {
    for (int x = 0; x < nx; x++) {
      A.boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + x ] = z*ny*nx + x;
      numOfBoundaryRows++;
    }
    for (int y = 1; y < ny - 1; y++) {
      A.boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + nx + 2*(y-1)] = (z*ny + y)*nx;
      numOfBoundaryRows++;
      A.boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + nx + 2*(y-1)+1] = (z*ny + y)*nx + nx - 1;
      numOfBoundaryRows++;
    }
    for (int x = 0; x < nx; x++) {
      A.boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + nx + 2*(ny-2) + x] = (z*ny + (ny - 1))*nx + x;
      numOfBoundaryRows++;
    }
  }
#ifndef HPCG_NO_OPENMP
  #pragma omp for nowait
#endif
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      A.boundaryRows[ny*nx + 2*(nz-2)*(nx+ny-2) + y*nx + x] = ((nz - 1)*ny + y)*nx + x;
      numOfBoundaryRows++;
    }
  }
}
  A.numOfBoundaryRows = numOfBoundaryRows;
  local_int_t localNumberOfNonzeros = 0;
  local_int_t *map_neib_r = (local_int_t*) MKL_malloc( sizeof(local_int_t)*A.geom->size, 512 );

  if ( map_neib_r == NULL ) return;

  for ( local_int_t i = 0; i < A.geom->size; i ++ ) map_neib_r[i] = 0;
#ifndef HPCG_NO_OPENMP
#pragma omp parallel reduction(+:localNumberOfNonzeros) num_threads(nproc)
#endif
{
#ifndef HPCG_NO_OPENMP
    int ithr = omp_get_thread_num();
    int nthr = omp_get_num_threads();
#else
    int ithr = 0;
    int nthr = 1;
#endif

    local_int_t works = (nz - 2)*(ny - 2);
    local_int_t begin = ((ithr  )*works)/nthr;
    local_int_t end   = ((ithr+1)*works)/nthr;
    for (local_int_t i = begin; i < end; i++)
    {
        local_int_t iz = i/(ny - 2) + 1;
        local_int_t iy = i%(ny - 2) + 1;
        
        for (local_int_t ix=1; ix<nx-1; ix++)
        {
            global_int_t giy = ipy*ny+iy;
            global_int_t giz = ipz*nz+iz;
            global_int_t gix = ipx*nx+ix;
            local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
            global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
            A.localToGlobalMap[currentLocalRow ] = currentGlobalRow;
            mtxIndG[currentLocalRow]      = A.mtxG + currentLocalRow*numberOfNonzerosPerRow;
            mtxIndL[currentLocalRow]      = A.mtxL + currentLocalRow*numberOfNonzerosPerRow;
            matrixValues[currentLocalRow] = A.mtxA + currentLocalRow*numberOfNonzerosPerRow;
            char numberOfNonzerosInRow = 0;
            double * currentValuePointer = matrixValues[currentLocalRow]; // Pointer to current value in current row
            global_int_t * currentIndexPointerG = mtxIndG[currentLocalRow]; // Pointer to current index in current row
            local_int_t  * currentIndexPointerL = mtxIndL[currentLocalRow]; // Pointer to current index in current row
            for (int sz=-1; sz<=1; sz++) {
                
                *(currentValuePointer + 0) = -1.0;
                *(currentValuePointer + 1) = -1.0;
                *(currentValuePointer + 2) = -1.0;
                *(currentValuePointer + 3) = -1.0;
                *(currentValuePointer + 4) = -1.0;
                *(currentValuePointer + 5) = -1.0;
                *(currentValuePointer + 6) = -1.0;
                *(currentValuePointer + 7) = -1.0;
                *(currentValuePointer + 8) = -1.0;
                
                local_int_t offset = currentLocalRow + sz*ny*nx;
                *(currentIndexPointerL + 0) = offset - nx - 1;
                *(currentIndexPointerL + 1) = offset - nx;
                *(currentIndexPointerL + 2) = offset - nx + 1;
                *(currentIndexPointerL + 3) = offset - 1;
                *(currentIndexPointerL + 4) = offset;
                *(currentIndexPointerL + 5) = offset + 1;
                *(currentIndexPointerL + 6) = offset + nx - 1;
                *(currentIndexPointerL + 7) = offset + nx;
                *(currentIndexPointerL + 8) = offset + nx + 1;
                currentValuePointer  += 9;
                currentIndexPointerL += 9;
            } // end sz loop
            *(currentValuePointer - 14) = 26.0;
            matrixDiagonal[currentLocalRow] = currentValuePointer - 14;
            numberOfNonzerosInRow += 27;
            nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
            localNumberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
            if (b!=0)      bv[currentLocalRow] = 26.0 - ((double) (numberOfNonzerosInRow-1));
            if (x!=0)      xv[currentLocalRow] = 0.0;
            if (xexact!=0) xexactv[currentLocalRow] = 1.0;
        } // end ix loop
    }
#ifndef HPCG_NO_OPENMP
#pragma omp for
#endif
    for (int i = 0; i < A.numOfBoundaryRows; i++) {
    local_int_t currentLocalRow = A.boundaryRows[i];

    local_int_t iz = currentLocalRow/(ny*nx);
    local_int_t iy = currentLocalRow/nx%ny;
    local_int_t ix = currentLocalRow%nx;

          global_int_t giz = ipz*nz+iz;
          global_int_t giy = ipy*ny+iy;
          global_int_t gix = ipx*nx+ix;
          
          global_int_t sz_begin = std::max<global_int_t>(-1, -giz);
          global_int_t sz_end = std::min<global_int_t>(1, gnz - giz - 1);
          
          global_int_t sy_begin = std::max<global_int_t>(-1, -giy);
          global_int_t sy_end = std::min<global_int_t>(1, gny - giy - 1);
          
          global_int_t sx_begin = std::max<global_int_t>(-1, -gix);
          global_int_t sx_end = std::min<global_int_t>(1, gnx - gix - 1);

//          local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;

          global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
          /*
          if( A.geom->size > 1 )
          {
              #pragma omp critical
              A.globalToLocalMap[currentGlobalRow] = currentLocalRow;
          }
          */
        mtxIndG[currentLocalRow]      = A.mtxG + currentLocalRow*numberOfNonzerosPerRow;
        mtxIndL[currentLocalRow]      = A.mtxL + currentLocalRow*numberOfNonzerosPerRow;
        matrixValues[currentLocalRow] = A.mtxA + currentLocalRow*numberOfNonzerosPerRow;
        A.localToGlobalMap[currentLocalRow] = currentGlobalRow;
        char numberOfNonzerosInRow = 0;
        double * currentValuePointer = matrixValues[currentLocalRow]; // Pointer to current value in current row
        global_int_t * currentIndexPointerG = mtxIndG[currentLocalRow]; // Pointer to current index in current row
        local_int_t  * currentIndexPointerL = mtxIndL[currentLocalRow];
        for (global_int_t sz=sz_begin; sz<=sz_end; sz++) {
            for (global_int_t sy=sy_begin; sy<=sy_end; sy++) {
                for (global_int_t sx=sx_begin; sx<=sx_end; sx++) {
                    global_int_t curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx;
                     local_int_t    col = currentLocalRow +sz*nx*ny+sy*nx+sx;
                    if (curcol==currentGlobalRow) {
                      matrixDiagonal[currentLocalRow] = currentValuePointer;
                      *currentValuePointer++ = 26.0;
                    } else {
                      *currentValuePointer++ = -1.0;
                    }
                    *currentIndexPointerG++ = curcol;
                    int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curcol);
                    if( A.geom->rank == rankIdOfColumnEntry )
                    {
                        *currentIndexPointerL++ = col;
                    } else {
                        map_neib_r[rankIdOfColumnEntry] ++;
                        *currentIndexPointerL++ = -1-rankIdOfColumnEntry;//(- col - 1);
                    }
//                      printf("%d %lf %d %d\n",currentLocalRow,currentValuePointer[-1],curcol,numberOfNonzerosInRow);
                    numberOfNonzerosInRow++;
                } // end sx loop
            } // end sy loop
        } // end sz loop
        nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
        localNumberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
        if( b!=0 ) bv[currentLocalRow] = 26.0 - ((double) (numberOfNonzerosInRow-1));
        if( x!=0 ) xv[currentLocalRow] = 0.0;
        if( xexact!=0 ) xexactv[currentLocalRow] = 1.0;
  }
}
    for (int i = 0; i < A.numOfBoundaryRows; i++) {
        local_int_t currentLocalRow = A.boundaryRows[i];
        
        local_int_t iz = currentLocalRow/(ny*nx);
        local_int_t iy = currentLocalRow/nx%ny;
        local_int_t ix = currentLocalRow%nx;
        
        global_int_t giz = ipz*nz+iz;
        global_int_t giy = ipy*ny+iy;
        global_int_t gix = ipx*nx+ix;
        global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
        A.globalToLocalMap[currentGlobalRow] = currentLocalRow;
    }

  local_int_t number_of_neighbors = 0;
  if( A.geom->size > 1 ) {
  for ( local_int_t i = 0; i < A.geom->size; i ++ ) number_of_neighbors += (map_neib_r[i] > 0); }
  A.work = map_neib_r;

  global_int_t totalNumberOfNonzeros = 0;
#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to sum all nonzeros
#ifdef HPCG_NO_LONG_LONG
  MPI_Allreduce(&localNumberOfNonzeros, &totalNumberOfNonzeros, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  long long lnnz = localNumberOfNonzeros, gnnz = 0; // convert to 64 bit for MPI call
  MPI_Allreduce(&lnnz, &gnnz, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  totalNumberOfNonzeros = gnnz; // Copy back
#endif
#else
  totalNumberOfNonzeros = localNumberOfNonzeros;
#endif
  // If this assert fails, it most likely means that the global_int_t is set to int and should be set to long long
  // This assert is usually the first to fail as problem size increases beyond the 32-bit integer range.
  assert(totalNumberOfNonzeros>0); // Throw an exception of the number of nonzeros is less than zero (can happen if int overflow)

  A.title = 0;
  A.numberOfSendNeighbors = number_of_neighbors;
  A.totalNumberOfRows = totalNumberOfRows;
  A.totalNumberOfNonzeros = totalNumberOfNonzeros;
  A.localNumberOfRows = localNumberOfRows;
  A.localNumberOfColumns = localNumberOfRows;
  A.localNumberOfNonzeros = localNumberOfNonzeros;
  A.nonzerosInRow = nonzerosInRow;
  A.mtxIndG = mtxIndG;
  A.mtxIndL = mtxIndL;
  A.matrixValues = matrixValues;
  A.matrixDiagonal = matrixDiagonal;
//  printf("GenerateProblem: ( %1.4f %1.4f %1.4f %1.4f ) %1.4f %1.4f\n",m1,m2,m3,m5,t4-t2,dsecnd()-t_1);fflush(0);
//  abort();
#endif
  return;
}

/*
  char *nonzerosInRow   = (char*)         MKL_malloc( sizeof(char        )*localNumberOfRows, ALIGN );
  global_int_t *mtxIndG = (global_int_t*) MKL_malloc( sizeof(global_int_t)*localNumberOfRows*numberOfNonzerosPerRow, ALIGN );
   local_int_t *mtxIndL = (local_int_t* ) MKL_malloc( sizeof( local_int_t)*localNumberOfRows*numberOfNonzerosPerRow, ALIGN );
  double *matrixValues  = (double*      ) MKL_malloc( sizeof(double      )*localNumberOfRows*numberOfNonzerosPerRow, ALIGN );
  double *matrixDiagonal= (double*      ) MKL_malloc( sizeof(double      )*localNumberOfRows, ALIGN );
*/
