''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
       Intel(R) Optimized High Performance Conjugate Gradient Benchmark
................................................................................

Package Contents 
----------------

This package contains Intel(R) Optimized High Performance Conjugate Gradient
Benchmark (Intel(R) Optimized HPCG) for Linux* OS.
The package is based on the reference implementation of the HPCG benchmark v3.0
with the following additions:

    bin/lib/intel64/libihpcg_avx.a
      - Library with HPCG* subroutines optimized for Intel(R) Advanced
        Vector Extensions (Intel(R) AVX).

    bin/lib/intel64/libihpcg_avx2.a
      - Library with HPCG* subroutines optimized for Intel(R) Advanced
        Vector Extensions 2 (Intel(R) AVX2).

    bin/lib/intel64/libihpcg_knl.a
      - Library with HPCG* subroutines optimized for Intel(R) Advanced
        Vector Extensions 512 (Intel(R) AVX512)(formerly Knights Landing).

    bin/lib/intel64/libihpcg_skx.a
      - Library with HPCG* subroutines optimized for Intel(R) Advanced
        Vector Extensions 512 (Intel(R) AVX512)(codename Skylake Server).

    setup/Make.*
      - Setup files for building a benchmark against optimizations provided
        in shared libraries. The build steps are described in the HPCG
        QUICKSTART guide. Intel(R) C++ Compiler should be present on a host
        used for building.

    bin/xhpcg_avx
    bin/xhpcg_avx2
    bin/xhpcg_knl
    bin/xhpcg_skx
      - For user convenience the package is supplied with prebuilt
        HPCG* launchers, allowing out-of-the-box benchmark execution.

        Setup runtime environment for Intel(R) C/C++ Compiler 15.0.x or later
        version. Run benchmark on Intel(R) AVX, Intel(R) AVX2 or Intel(R) AVX-512
        enabled Intel processors as follows:
            #> export OMP_NUM_THREADS=<N_processors * N_cores>
            #> export KMP_AFFINITY=compact
            #> mpiexec.hydra -np $num_procs bin/xhpcg_avx -n$size -t$time_to_run
            #> mpiexec.hydra -np $num_procs bin/xhpcg_avx2 -n$size -t$time_to_run
            #> mpiexec.hydra -np $num_procs bin/xhpcg_knl -n$size -t$time_to_run
            #> mpiexec.hydra -np $num_procs bin/xhpcg_skx -n$size -t$time_to_run

        Note: Problem size should be multiple of 8 and at least 24.

Reference HPCG code in the ./src directory was modified in order
to bundle Intel(R) architecture optimizations with the shared libraries.

===============================================================================
  Intel and the Intel logo are trademarks of Intel Corporation in the U.S.
  and/or other countries.

* Other names and brands may be claimed as the property of others.

  Copyright(C) 2014 Intel Corporation. All rights reserved.

  This Intel(R) Optimized  Technology Preview  for High  Performance  Conjugate
  Gradient  Benchmark ("Software")  is furnished under license  and may only be 
  used or copied  in accordance  with the terms  of that  license.  No license, 
  express  or implied, by estoppel  or otherwise,  to any intellectual property 
  rights is granted by this document. The Software is subject to change without 
  notice,  and should not be construed  as a commitment by Intel Corporation to 
  market, license, sell  or support any product or technology. Unless otherwise 
  provided   for  in  the  license  under  which  this  Software  is  provided, 
  the Software  is provided  AS IS,  with no warranties of any kind, express or 
  implied. Except as expressly permitted by the Software license, neither Intel
  Corporation nor its suppliers assumes any responsibility or liability for any
  errors  or inaccuracies that may appear herein. Except as expressly permitted
  by the Software license,  no part  of the Software  may be reproduced, stored 
  in a retrieval system,  transmitted in any form,  or distributed by any means 
  without the express written consent of Intel Corporation.
======================================================================
 -- High Performance Conjugate Gradients (HPCG) Benchmark
    HPCG - 3.0 - November 11, 2015

    Michael A. Heroux
    Scalable Algorithms Group, Center for Computing Research
    Sandia National Laboratories, Albuquerque, NM

    Piotr Luszczek
    Jack Dongarra
    University of Tennessee, Knoxville
    Innovative Computing Laboratory
    (C) Copyright 2013 All Rights Reserved

 -- Copyright notice and Licensing terms:

 Redistribution  and  use in  source and binary forms, with or without
 modification, are  permitted provided  that the following  conditions
 are met:

 1. Redistributions  of  source  code  must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce  the above copyright
 notice, this list of conditions,  and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. The name of the  University,  the name of the  Laboratory,  or the
 names  of  its  contributors  may  not  be used to endorse or promote
 products  derived   from   this  software  without  specific  written
 permission.

 -- Disclaimer:

 THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
======================================================================
