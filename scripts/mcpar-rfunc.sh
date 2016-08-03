#!/bin/zsh

#SBATCH -n 96
#SBATCH -t 30
#SBATCH -A GCAM

module load openmpi
module load mkl/15.0.1
module load gcc/6.1.0

date

time mpirun -np 96 ./mcpar-rfunc /people/link593/wrk/food-demand/src/R/food-demand.R  /people/link593/wrk/food-demand/data/obsdata-test-grand.csv ../runs/rfunc-test.dat 1 5

date

