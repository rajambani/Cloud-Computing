#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

./MyDiskBench.exe disk-RR-1-2thread.dat disk-RR-1-1thread.out.dat

