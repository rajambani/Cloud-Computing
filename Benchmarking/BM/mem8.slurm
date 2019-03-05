#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute
#SBATCH --output=pmbw-10-2.out

pmbw -p 2 -P 2 -s 1048576 -S 1048576 -M 1G -f ScanWrite64PtrSimpleLoop
