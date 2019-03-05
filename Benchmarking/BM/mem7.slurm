#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute
#SBATCH --output=pmbw-10-1.out

pmbw -p 1 -P 1 -s 1048576 -S 1048576 -M 1G -f ScanWrite64PtrSimpleLoop
