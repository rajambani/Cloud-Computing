#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute
#SBATCH --output=pmbw-1-1random.out

pmbw -p 1 -P 1 -s 1024 -S 1024 -M 1G -f PermRead64SimpleLoop

