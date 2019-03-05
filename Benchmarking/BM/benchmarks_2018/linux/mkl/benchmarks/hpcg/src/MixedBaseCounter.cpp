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

#include <map>

#include "MixedBaseCounter.hpp"

MixedBaseCounter::MixedBaseCounter(int *counts, int length) {
  this->length = length;

  int i;

  for (i = 0; i < 32; ++i) {
    this->max_counts[i] = counts[i];
    this->cur_counts[i] = 0;
  }
  // terminate with 0's
  this->max_counts[i]      = this->cur_counts[i]      = 0;
  this->max_counts[length] = this->cur_counts[length] = 0;
}

MixedBaseCounter::MixedBaseCounter(MixedBaseCounter & left, MixedBaseCounter & right) {
  this->length = left.length;
  for (int i = 0; i < left.length; ++i) {
    this->max_counts[i] = left.max_counts[i] - right.cur_counts[i];
    this->cur_counts[i] = 0;
  }
}

void
MixedBaseCounter::next() {
  for (int i = 0; i < this->length; ++i) {
    this->cur_counts[i]++;
    if (this->cur_counts[i] > this->max_counts[i]) {
      this->cur_counts[i] = 0;
      continue;
    }
    break;
  }
}

int
MixedBaseCounter::is_zero() {
  for (int i = 0; i < this->length; ++i)
    if (this->cur_counts[i])
      return 0;
  return 1;
}

int
MixedBaseCounter::product(int * multipliers) {
  int k=0, x=1;

  for (int i = 0; i < this->length; ++i)
    for (int j = 0; j < this->cur_counts[i]; ++j) {
      k = 1;
      x *= multipliers[i];
    }

  return x * k;
}
