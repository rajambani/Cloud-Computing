#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

./MyRAMBench.exe memory-RWS-1000-4thread.dat memory-RWS-1000-4thread.out.dat

