#!/bin/bash

#SBATCH --job-name="linsort20GB"
#SBATCH --output="linsort20GB.log"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

START_TIME=$SECONDS
LC_ALL=c sort /input/data-20GB.in > /tmp/data-20GB.out
ELAPSED_TIME=$(($SECONDS - $START_TIME))

valsort /tmp/data-20GB.out

echo "Total time for linux sort: $ELAPSED_TIME seconds" >> /exports/home/rambani/cs553-pa2a/linsort20GB.log

