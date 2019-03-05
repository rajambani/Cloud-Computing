#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute
#SBATCH --output=pmbw-10-4.out

pmbw -p 4 -P 4 -s 1048576 -S 1048576 -M 1G -f ScanWrite64PtrSimpleLoop
