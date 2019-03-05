#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

./MyCPUBench.exe cpu_HP_2thread.dat cpu_HP_2thread.out.dat

