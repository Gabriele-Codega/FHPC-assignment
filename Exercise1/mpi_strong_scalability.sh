#!/bin/sh
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

echo Running MPI strong scalability test.
echo Size is 5k, for 2000 steps, by socket.
mpirun -np 8 GameOfLife -i -f imgs/mpi_strong_init.pgm -k 5000,5000

touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
for n in $(seq 1 $SLURM_NTASKS)
do
    echo Currently using $n tasks.
    for i in $(seq 1 5)
    do
        mpirun -np $n --map-by socket --report-bindings GameOfLife -r -f imgs/mpi_strong_init.pgm -n 2000 -s 0 -e 1 -t $1
    done
done
echo done!

squeue -j $SLURM_JOB_ID