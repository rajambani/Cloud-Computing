#!/bin/bash

echo "hello, $USER."

gcc -o MyCPUBench.exe MyCPUBench.c -lpthread -w

for i in `seq 11 12`;
        do
                echo $i
		sbatch run$i.slurm
        done

squeue -u rambani


