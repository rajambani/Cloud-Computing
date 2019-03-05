#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

./MyCPUBench.exe cpu_HP_4thread.dat cpu_HP_4thread.out.dat

