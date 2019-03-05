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
 @file ComputeMG.cpp

 HPCG routine
 */

#include "ComputeMG.hpp"
#include "ComputeSYMGS.hpp"
#include "ComputeSPMV.hpp"
#include "ComputeRestriction_ref.hpp"
#include "ComputeProlongation_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#include <mpi.h>
#include "Geometry.hpp"
#include <cstdlib>
#endif

int ComputeMG_ref(const SparseMatrix & A, const Vector & r, Vector & x);

/*!
  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG_ref
*/
int ComputeMG(const SparseMatrix  & A, const Vector & r, Vector & x)
{
    int ierr = 0;

    if (A.mgData!=0) // Go to next coarse level if defined
    {
        const int numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps;
        const int numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps;

        if ( numberOfPresmootherSteps > 1 )
        {
            sparse_status_t status = SPARSE_STATUS_SUCCESS;
            struct optData *optData = (struct optData *)A.optimizationData;
            struct matrix_descr descr;
            sparse_matrix_t csrA = (sparse_matrix_t)optData->csrA;
            descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
            descr.mode = SPARSE_FILL_MODE_FULL;
            descr.diag = SPARSE_DIAG_NON_UNIT;

            if ( mkl_sparse_d_symgs(SPARSE_OPERATION_NON_TRANSPOSE, csrA, descr, 0.0, r.values, x.values) != SPARSE_STATUS_SUCCESS ) ierr ++;

            for ( int i = 1; i < numberOfPresmootherSteps; ++i ) ierr += ComputeSYMGS(A, r, x);
            ierr += ComputeSPMV(A, x, (*A.mgData->Axf));
        } else
        {
            ierr += ComputeSYMGS_MV(A, r, x, (*A.mgData->Axf));
        }

        ierr += ComputeRestriction_ref(A, r);
        ierr += ComputeMG(*A.Ac,*A.mgData->rc, *A.mgData->xc);
        ierr += ComputeProlongation_ref(A, x);

        for ( int i = 0; i < numberOfPostsmootherSteps; ++i ) ierr += ComputeSYMGS(A, r, x);

        if ( ierr != 0 ) return 1;
    } else
    {
        sparse_status_t status = SPARSE_STATUS_SUCCESS;
        struct optData *optData = (struct optData *)A.optimizationData;
        struct matrix_descr descr;
        sparse_matrix_t csrA = (sparse_matrix_t)optData->csrA;
        descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        descr.mode = SPARSE_FILL_MODE_FULL;
        descr.diag = SPARSE_DIAG_NON_UNIT;

        status = mkl_sparse_d_symgs(SPARSE_OPERATION_NON_TRANSPOSE, csrA, descr, 0.0, r.values, x.values);

        if ( status != SPARSE_STATUS_SUCCESS ) return 1;
    }

    return 0;
}
