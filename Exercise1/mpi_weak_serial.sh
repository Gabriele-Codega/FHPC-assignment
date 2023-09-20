#!/bin/sh
#SBATCH --output=mpi_weak_serial_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:40:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

echo Running MPI weak scalability test.
echo size is 10k x 1k per process. 1000 steps. Map by numa.
touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
for n in $(seq 1 32)
do
    echo Current size is $n.
    mpirun -np 1 --map-by numa GameOfLife -r -f imgs/weak_init_$n.pgm -n 1000 -s 0 -e 1 -t $1
done
echo done!