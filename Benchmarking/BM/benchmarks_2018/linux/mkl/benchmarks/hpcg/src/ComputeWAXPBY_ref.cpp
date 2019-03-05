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
 @file ComputeWAXPBY_ref.cpp

 HPCG routine
 */

#include "ComputeWAXPBY_ref.hpp"
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
/*!
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This is the reference WAXPBY impmentation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY
*/
int ComputeWAXPBY_ref(const local_int_t n, const double alpha, const Vector & x,
    const double beta, const Vector & y, Vector & w) {

    assert(x.localLength>=n); // Test vector lengths
    assert(y.localLength>=n);

    const double * const xv = x.values;
    const double * const yv = y.values;
    double * const wv = w.values;
  if (alpha==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = xv[i] + beta * yv[i];
  } else if (beta==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + yv[i];
  } else  {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + beta * yv[i];
  }

  return 0;
}
