#!/bin/bash -l
#SBATCH --job-name="test"
#SBATCH --nodes=125
#SBATCH --time=00:05:00
#SBATCH --output="testoutput.txt"
#SBATCH --error="errorlog.txt"
echo "Nodes on which it executes"
echo $SLURM_JOB_NODELIST
aprun -n 125 -d 16 -N 1 ./DistGEMM.exe 40000 125 5

