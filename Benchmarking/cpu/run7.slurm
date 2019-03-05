#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

./MyCPUBench.exe cpu_SP_1thread.dat cpu_SP_1thread.out.dat

