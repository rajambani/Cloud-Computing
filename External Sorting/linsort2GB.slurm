#!/bin/bash

#SBATCH --job-name="linsort2GB"
#SBATCH --output="linsort2GB.log"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

START_TIME=$SECONDS
LC_ALL=c sort /input/data-2GB.in > /tmp/data-2GB.out
ELAPSED_TIME=$(($SECONDS - $START_TIME))

valsort /tmp/data-2GB.out

echo "Total time for linux sort: $ELAPSED_TIME seconds" >> /exports/home/rambani/cs553-pa2a/linsort2GB.log

