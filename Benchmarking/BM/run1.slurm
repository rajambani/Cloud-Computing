#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --wait-all-nodes=1

srun ./run1.sh $SLURM_JOB_NODELIST


