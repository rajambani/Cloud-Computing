#!/bin/bash

echo "hello, $USER."
j=0
javac MyNETBenchTCP.java

for i in `seq 11 12`;
        do
                echo $i
		sbatch run$i.slurm
		sleep 5
		j=$(( $i + $20 ))
		sbatch run$j.slurm
		sleep 5
        done

squeue -u rambani

