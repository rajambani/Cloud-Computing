#!/bin/bash

echo "hello, $USER."
j=21
javac MyNETBenchTCP.java

for i in `seq 1 12`;
        do
                echo $i
		sbatch run$i.slurm
		sleep 5
		j=$(( $i + $20 ))
		echo "$j"
		sbatch run$j.slurm
		sleep 5
        done

squeue -u rambani


