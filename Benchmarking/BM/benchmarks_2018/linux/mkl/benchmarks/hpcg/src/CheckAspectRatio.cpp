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
 @file CheckAspectRatio.cpp

 HPCG routine
 */

#include <algorithm>

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

#include "hpcg.hpp"

#include "CheckAspectRatio.hpp"

int
CheckAspectRatio(double smallest_ratio, int x, int y, int z, const char *what, bool DoIo) {
  double current_ratio = std::min(std::min(x, y), z) / double(std::max(std::max(x, y), z));

  if (current_ratio < smallest_ratio) { // ratio of the smallest to the largest
    if (DoIo) {
      HPCG_fout << "The " << what << " sizes (" << x << "," << y << "," << z <<
        ") are invalid because the ratio min(x,y,z)/max(x,y,z)=" << current_ratio <<
        " is too small (at least " << smallest_ratio << " is required)." << std::endl;
      HPCG_fout << "The shape should resemble a 3D cube. Please adjust and try again." << std::endl;
      HPCG_fout.flush();
    }

#ifndef HPCG_NO_MPI
    MPI_Abort(MPI_COMM_WORLD, 127);
#endif

    return 127;
  }

  return 0;
}
