#!/bin/bash -l
#SBATCH --job-name="test"
#SBATCH --nodes=8
#SBATCH --time=00:00:20
#SBATCH --output="testoutput.txt"
#SBATCH --error="errorlog.txt"
echo "Nodes on which it executes"
echo $SLURM_JOB_NODELIST
aprun -n 8 -d 16 -N 1 ./DistGEMM.exe $Nnodes 800 8 2

