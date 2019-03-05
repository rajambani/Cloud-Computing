#!/bin/bash

echo "Hello, $USER."

gcc -o MyDiskBench.exe MyDiskBench.c -lpthread

for i in `seq 31 40`;
        do
                echo $i
		sbatch run$i.slurm
        done

squeue -u rambani


