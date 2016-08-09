#!/bin/zsh

#SBATCH -n 192
#SBATCH -t 24:00:00
#SBATCH -A GCAM

module load openmpi
module load mkl/15.0.1
module load gcc/6.1.0

date

time mpirun -np 192 ./mcpar-rfunc /people/link593/wrk/food-demand/src/R/food-demand-mc.R  /people/link593/wrk/food-demand/data/food-dmnd-price-hist.csv ../runs/rfunc-hist.dat 1000 2000

date

