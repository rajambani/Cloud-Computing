#!/bin/bash

#SBATCH --job-name="mysort20GB"
#SBATCH --output="mysort20GB.log"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

java MySort /input/data-20GB.in

echo "-------------------------------Valsort Results-----------------------------------" >> ~/cs553-pa2a/mysort20GB.log
valsort /tmp/output/output
