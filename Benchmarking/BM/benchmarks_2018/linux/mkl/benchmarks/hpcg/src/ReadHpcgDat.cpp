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

#include <cstdio>

#include "ReadHpcgDat.hpp"

static int
SkipUntilEol(FILE *stream) {
  int chOrEof;
  bool finished;

  do {
    chOrEof = fgetc( stream );
    finished = (chOrEof == EOF) || (chOrEof == '\n') || (chOrEof == '\r');
  } while (! finished);

  if ('\r' == chOrEof) { // on Windows, \r might be followed by \n
    int chOrEofExtra = fgetc( stream );

    if ('\n' == chOrEofExtra || EOF == chOrEofExtra)
      chOrEof = chOrEofExtra;
    else
      ungetc(chOrEofExtra, stream);
  }

  return chOrEof;
}

int
ReadHpcgDat(int *localDimensions, int *secondsPerRun) {
  FILE * hpcgStream = fopen("hpcg.dat", "r");

  if (! hpcgStream)
    return -1;

  SkipUntilEol(hpcgStream); // skip the first line

  SkipUntilEol(hpcgStream); // skip the second line

  for (int i = 0; i < 3; ++i)
    if (fscanf(hpcgStream, "%d", localDimensions+i) != 1 || localDimensions[i] < 16)
      localDimensions[i] = 16;

  SkipUntilEol( hpcgStream ); // skip the rest of the second line

  if (secondsPerRun!=0) { // Only read number of seconds if the pointer is non-zero
	  if (fscanf(hpcgStream, "%d", secondsPerRun) != 1 || secondsPerRun[0] < 0)
		  secondsPerRun[0] = 30 * 60; // 30 minutes
  }

  fclose(hpcgStream);

  return 0;
}
