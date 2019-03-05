#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

./MyRAMBench.exe memory-RWR-10000000-1thread.dat memory-RWR-10000000-1thread.out.dat

