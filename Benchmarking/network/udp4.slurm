#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --wait-all-nodes=1

srun ./run1.sh $SLURM_JOB_NODELIST

java MyNETBenchUDP network-UDP-32000-1thread.dat 1

