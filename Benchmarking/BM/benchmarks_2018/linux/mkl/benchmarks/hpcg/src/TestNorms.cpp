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
 @file TestNorms.cpp

 HPCG routine
 */

#include <cmath>
#include "TestNorms.hpp"

/*!
  Computes the mean and standard deviation of the array of norm results.

  @param[in] testnorms_data data structure with the results of norm test

  @return Returns 0 upon success or non-zero otherwise
*/
int TestNorms(TestNormsData & testnorms_data) {
 double mean_delta = 0.0;
 for (int i= 0; i<testnorms_data.samples; ++i) mean_delta += (testnorms_data.values[i] - testnorms_data.values[0]);
 double mean = testnorms_data.values[0] + mean_delta/(double)testnorms_data.samples;
 testnorms_data.mean = mean;

 // Compute variance
 double sumdiff = 0.0;
 for (int i= 0; i<testnorms_data.samples; ++i) sumdiff += (testnorms_data.values[i] - mean) * (testnorms_data.values[i] - mean);
 testnorms_data.variance = sumdiff/(double)testnorms_data.samples;

 // Determine if variation is sufficiently small to declare success
 testnorms_data.pass = (testnorms_data.variance<1.0e-6);

 return 0;
}
