#!/bin/bash -l
export OMP_NUM_THREADS 16
NNodes = 8;
#SBATCH --job-name="test"
#SBATCH --nodes = $NNodes
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00
#SBATCH --output="testoutput.txt"

echo "Nodes on which it executes"
echo $SLURM_JOB_NODELIST
aprun -n $NNodes -d 16 -N 1 ./DistGEMM.exe $Nnodes 8 2

