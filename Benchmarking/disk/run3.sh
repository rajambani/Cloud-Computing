#!/bin/bash

echo "Hello, $USER."

gcc -o MyDiskBench.exe MyDiskBench.c -lpthread

for i in `seq 21 30`;
        do
                echo $i
		sbatch run$i.slurm
        done

squeue -u rambani


