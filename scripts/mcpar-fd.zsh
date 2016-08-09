#!/bin/zsh

#SBATCH -n 192
#SBATCH -t 14:00:00
#SBATCH -A GCAM

module load openmpi
module load mkl/15.0.1
module load gcc/6.1.0

RFILE="$HOME/wrk/food-demand/src/R/food-demand-mc.R"
DFILE="$HOME/wrk/food-demand/data/food-dmnd-price.$SLURM_ARRAY_TASK_ID.dat"
OFILE="mc-food-dmnd.$SLURM_ARRAY_TASK_ID.dat"

date
time mpirun -np 192 ./mcpar-rfunc $RFILE $DFILE $OFILE 1000 2000
date
