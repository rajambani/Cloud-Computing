#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

./MyRAMBench.exe memory-RWS-10000000-2thread.dat memory-RWS-10000000-2thread.out.dat

