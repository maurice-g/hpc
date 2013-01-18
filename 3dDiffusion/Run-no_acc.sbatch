#!/bin/bash -l
#SBATCH --job-name="meassinglenode"
#SBATCH --time=00:59:00
#SBATCH --output="3dDiffusion-1Node-intel03-noacc-firsttouch.txt"
#SBATCH --error="errorlog2.txt"
export OMP_NUM_THREADS=16
for i in $(seq 1 50)
do
	meshsize=$((i*20))
 	aprun -n 1 -N 1 -d 16 ./Diffusion3d.exe 0.25 $meshsize $meshsize $meshsize 1 1 1
done 
