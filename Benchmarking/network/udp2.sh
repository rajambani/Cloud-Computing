#!/bin/bash

echo "hello, $USER."
j=0
javac MyNETBenchTCP.java

for i in `seq 6 10`;
        do
                echo $i
		sbatch udp$i.slurm
		sleep 5
		j=$(( $i + $20 ))
		sbatch udp$j.slurm
		sleep 5
        done

squeue -u rambani

