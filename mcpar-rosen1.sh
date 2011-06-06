#PBS -l nodes=8
#PBS -l walltime=30:00
#PBS -A GCCARBON

wdir=/homes/rpl/wrk/mcpar
cd $wdir

mpirun -np 8 $wdir/mcpar-rosen1-mpi.exe

