#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

./MyRAMBench.exe memory-RWR-1000000-2thread.dat memory-RWR-1000000-2thread.out.dat

