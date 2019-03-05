#!/bin/bash

#SBATCH --job-name="mysort2GB"
#SBATCH --output="mysort2GB.log"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

java MySort /input/data-2GB.in

echo "-------------------------------Valsort Results-----------------------------------" >> ~/cs553-pa2a/mysort2GB.log
valsort /tmp/output/output