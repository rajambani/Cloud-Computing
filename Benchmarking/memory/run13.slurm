#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

./MyRAMBench.exe memory-RWS-1-1thread.dat memory-RWS-1-1thread.out.dat

