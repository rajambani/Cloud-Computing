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
 @file ComputeResidual.cpp

 HPCG routine
 */
#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "Vector.hpp"

#ifdef HPCG_DETAILED_DEBUG
#include <fstream>
#include "hpcg.hpp"
#endif

#include <cmath>  // needed for fabs
#include "ComputeResidual.hpp"
#ifdef HPCG_DETAILED_DEBUG
#include <iostream>
#endif

/*!
  Routine to compute the inf-norm difference between two vectors where:

  @param[in]  n        number of vector elements (local to this processor)
  @param[in]  v1, v2   input vectors
  @param[out] residual pointer to scalar value; on exit, will contain result: inf-norm difference

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeResidual(const local_int_t n, const Vector & v1, const Vector & v2, double & residual) {

  double * v1v = v1.values;
  double * v2v = v2.values;
  double local_residual = 0.0;

#ifndef HPCG_NO_OPENMP
  #pragma omp parallel default(none) shared(local_residual, v1v, v2v)
  {
    double threadlocal_residual = 0.0;
    #pragma omp for
    for (local_int_t i=0; i<n; i++) {
      double diff = std::fabs(v1v[i] - v2v[i]);
      if (diff > threadlocal_residual) threadlocal_residual = diff;
    }
    #pragma omp critical
    {
      if (threadlocal_residual>local_residual) local_residual = threadlocal_residual;
    }
  }
#else // No threading
  for (local_int_t i=0; i<n; i++) {
    double diff = std::fabs(v1v[i] - v2v[i]);
    if (diff > local_residual) local_residual = diff;
#ifdef HPCG_DETAILED_DEBUG
    HPCG_fout << " Computed, exact, diff = " << v1v[i] << " " << v2v[i] << " " << diff << std::endl;
#endif
  }
#endif

#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to collect all partial sums
  double global_residual = 0;
  MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  residual = global_residual;
#else
  residual = local_residual;
#endif

  return 0;
}
