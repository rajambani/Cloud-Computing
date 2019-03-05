#!/bin/bash

echo "hello, $USER."

for i in `seq 1 3`;
        do
                echo $i
		sbatch runLinpack$i.slurm
        done

squeue -u rambani




