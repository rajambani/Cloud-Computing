#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --wait-all-nodes=1

java MyNETBenchUDP network-UDP-32000-4thread.dat 2

