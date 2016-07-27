#!/bin/zsh

#SBATCH -n 24
#SBATCH -t 60
#SBATCH -A GCAM

module load openmpi
module load mkl/15.0.1
module load gcc/6.1.0

mpirun -np 24 ./mcpar-rosen1.exe

