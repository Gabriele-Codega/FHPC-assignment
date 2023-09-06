#!/bin/sh
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo Running MPI strong scalability test.
mpirun -np 16 GameOfLife -i -f imgs/mpi_strong_init.pgm -k 1000,1000

touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
for n in $(seq 1 $SLURM_NTASKS_PER_NODE)
do
    echo Currently using $n tasks.
    mpirun -np $n --map-by core GameOfLife -r -f imgs/mpi_strong_init.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!