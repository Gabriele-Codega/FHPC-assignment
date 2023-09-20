#!/bin/sh
#SBATCH --output=ordered_test_%j.out
#SBATCH --partition=THIN
#SBATCH --job-name=Game_of_Life

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00

module load architecture/Intel
# module load openMPI/4.1.5/gnu/12.2.1
module load intelMPI/2021.7.1


# srun -n 1 make timeit

size=1000
nsteps=10

mpirun -np $SLURM_NTASKS --map-by socket ./GameOfLife -r -f imgs/omp_strong_init.pgm -n $nsteps -s 0 -e 0 -t $1
