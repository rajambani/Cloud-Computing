#!/bin/bash

echo "Hello, $USER."

gcc -o MyRAMBench.exe MyRAMBench.c -lpthread

for i in `seq 20 24`;
        do
                echo $i
		sbatch run$i.slurm
        done

squeue -u rambani


