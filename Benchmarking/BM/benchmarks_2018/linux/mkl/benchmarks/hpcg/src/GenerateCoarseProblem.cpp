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

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include <cassert>
#include "GenerateCoarseProblem.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "SetupHalo.hpp"

/*!
  Routine to construct a prolongation/restriction operator for a given fine grid matrix
  solution (as computed by a direct solver).

  @param[inout]  Af - The known system matrix, on output its coarse operator, fine-to-coarse operator and auxiliary vectors will be defined.

  Note that the matrix Af is considered const because the attributes we are modifying are declared as mutable.

*/

void GenerateCoarseProblem(const SparseMatrix & Af) {

  // Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  // below may result in global range values.
  global_int_t nxf = Af.geom->nx;
  global_int_t nyf = Af.geom->ny;
  global_int_t nzf = Af.geom->nz;

  local_int_t nxc, nyc, nzc; //Coarse nx, ny, nz
  assert(nxf%2==0); assert(nyf%2==0); assert(nzf%2==0); // Need fine grid dimensions to be divisible by 2
  nxc = nxf/2; nyc = nyf/2; nzc = nzf/2;
  local_int_t * f2cOperator = new local_int_t[Af.localNumberOfRows];
  local_int_t localNumberOfRows = nxc*nyc*nzc; // This is the size of our subblock
  // If this assert fails, it most likely means that the local_int_t is set to int and should be set to long long
  assert(localNumberOfRows>0); // Throw an exception of the number of rows is less than zero (can happen if "int" overflows)

  // Use a parallel loop to do initial assignment:
  // distributes the physical placement of arrays of pointers across the memory system
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< localNumberOfRows; ++i) {
    f2cOperator[i] = 0;
  }


  // TODO:  This triply nested loop could be flattened or use nested parallelism
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t izc=0; izc<nzc; ++izc) {
	  local_int_t izf = 2*izc;
	  for (local_int_t iyc=0; iyc<nyc; ++iyc) {
		  local_int_t iyf = 2*iyc;
		  for (local_int_t ixc=0; ixc<nxc; ++ixc) {
			  local_int_t ixf = 2*ixc;
			  local_int_t currentCoarseRow = izc*nxc*nyc+iyc*nxc+ixc;
			  local_int_t currentFineRow = izf*nxf*nyf+iyf*nxf+ixf;
			  f2cOperator[currentCoarseRow] = currentFineRow;
		  } // end iy loop
	  } // end even iz if statement
  } // end iz loop

  // Construct the geometry and linear system
  Geometry * geomc = new Geometry;
  GenerateGeometry(Af.geom->size, Af.geom->rank, Af.geom->numThreads, nxc, nyc, nzc, geomc);

  SparseMatrix * Ac = new SparseMatrix;
  InitializeSparseMatrix(*Ac, geomc);
  Ac->nproc = Af.nproc;
  GenerateProblem(*Ac, 0, 0, 0);
  SetupHalo(*Ac);
  Vector *rc = new Vector;
  Vector *xc = new Vector;
  Vector * Axf = new Vector;
  InitializeVector(*rc, Ac->localNumberOfRows);
  ZeroVector(*rc);
  InitializeVector(*xc, Ac->localNumberOfColumns);
  ZeroVector(*xc);
  InitializeVector(*Axf, Af.localNumberOfColumns);
  ZeroVector(*Axf);
  Af.Ac = Ac;
  MGData * mgData = new MGData;
  InitializeMGData(f2cOperator, rc, xc, Axf, *mgData);
  //printf("%p, %d\n",mgData,localNumberOfRows);fflush(0);
  Af.mgData = mgData;

  return;
}
