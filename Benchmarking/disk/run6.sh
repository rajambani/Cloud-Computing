#!/bin/bash

echo "Hello, $USER."

gcc -o MyDiskBench.exe MyDiskBench.c -lpthread -w

for i in `seq 51 52`;
        do
                echo $i
		sbatch run$i.slurm
        done

squeue -u rambani


