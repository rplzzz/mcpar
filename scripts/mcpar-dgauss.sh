#!/bin/zsh

#SBATCH n 24
#SBATCH -t 60
#SBATCH -A GCAM

module load openmpi
module load mkl/15.0.1
module load gcc/6.1.0

mpirun -np 8 $wdir/mcpar-dgauss-mpi

