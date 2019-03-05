#!/bin/bash

echo "hello, $USER."

gcc -o MyCPUBench.exe MyCPUBench.c -lpthread -w

for i in `seq 6 10`;
        do
                echo $i
		sbatch run$i.slurm
		sleep 5
		sbatch run4$(20+i).slurm
        done

squeue -u rambani


