#!/bin/sh
#SBATCH --output=gol_test_%j.out
#SBATCH --partition=THIN
#SBATCH --job-name=Game_of_Life

#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00

module load architecture/Intel
# module load openMPI/4.1.5/gnu/12.2.1
module load intelMPI/2021.7.1


srun -n 1 make timeit

size=1000
nsteps=4

mpirun -np $SLURM_NTASKS --map-by socket ./GameOfLife -i -k 1000,1000 -f imgs/test_init.pgm -t $1
