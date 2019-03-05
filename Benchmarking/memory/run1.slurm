#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

./MyRAMBench.exe memory-RWR-1-1thread.dat memory-RWR-1-1thread.out.dat
