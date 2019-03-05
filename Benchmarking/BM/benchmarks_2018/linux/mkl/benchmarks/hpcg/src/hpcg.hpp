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
 @file hpcg.hpp

 HPCG data structures and functions
 */

#ifndef HPCG_HPP
#define HPCG_HPP

#include <fstream>

//#define DBL_EPSILON 2.2204460492503131E-16
#define fabs_hpcg(a) ((a) >= 0.0) ? (a) : -(a)

extern std::ofstream HPCG_fout;

struct HPCG_Params_STRUCT {
  int comm_size; //!< Number of MPI processes in MPI_COMM_WORLD
  int comm_rank; //!< This process' MPI rank in the range [0 to comm_size - 1]
  int numThreads; //!< This process' number of threads
  int nx; //!< Number of x-direction grid points for each local subdomain
  int ny; //!< Number of y-direction grid points for each local subdomain
  int nz; //!< Number of z-direction grid points for each local subdomain
  int runningTime; //!< Number of seconds to run the timed portion of the benchmark
  int runRealRef;  // default true, turn on reference implementation
  char yamlFileName[1024];
 
};
/*!
  HPCG_Params is a shorthand for HPCG_Params_STRUCT
 */
typedef HPCG_Params_STRUCT HPCG_Params;

extern int HPCG_Init(int * argc_p, char ** *argv_p, HPCG_Params & params);
extern int HPCG_Finalize(void);

#endif // HPCG_HPP
