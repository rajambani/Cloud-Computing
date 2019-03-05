#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

cd /exports/home/rambani/cs553-pa1/BM/benchmarks_2018/linux/mkl/benchmarks/linpack
export OMP_NUM_THREADS=1
./runme_xeon64

