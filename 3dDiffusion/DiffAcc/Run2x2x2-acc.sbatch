#!/bin/bash -l
#SBATCH --job-name="Diffusioncopy"
#SBATCH --nodes=8
#SBATCH --time=02:00:00
#SBATCH --output="3dDiffusion-8Node-acccopy.txt"
#SBATCH --error="errorlog2x2x2acc.txt"
export OMP_NUM_THREADS=16
for i in $(seq 1 30)
do
	meshsize=$((i*50))
 	aprun -n 8 -N 1 -d 16 ./Diffusion3d.exe 0.25 $meshsize $meshsize $meshsize 2 2 2
done 
