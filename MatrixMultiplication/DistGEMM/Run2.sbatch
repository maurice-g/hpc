#!/bin/bash -l
#SBATCH --job-name="benchmark_DistGEMM"
#SBATCH --nodes=8
#SBATCH --time=00:30:00
#SBATCH --output="outbench2x2x2.txt"
#SBATCH --error="outbench2x2x2.err"
numNodes=8
cubeSize=2
echo "Nodes on which it executes"
echo $SLURM_JOB_NODELIST
for i in `seq 91 100`
do
	matrixSize=$((i*512))
	echo "Starting benchmark for matrix size $matrixSize ..."
	aprun -n $numNodes -d 16 -N 1 ./DistGEMM.exe $matrixSize $numNodes $cubeSize
done

