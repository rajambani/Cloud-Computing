#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

./MyDiskBench.exe disk-WS-10000-1thread.dat disk-RR-1-1thread.out.dat

