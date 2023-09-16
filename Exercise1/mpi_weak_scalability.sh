#!/bin/sh
#SBATCH --output=mpi_weak_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=4
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

echo Running MPI weak scalability test.
echo size is 10k x 1k per process. 1000 steps. Map by numa.
touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
for n in $(seq 14 $SLURM_NTASKS)
do
    echo Currently using $n tasks.
    echo Generating grid.
    mpirun -np $n GameOfLife -i -f imgs/weak_init_$n.pgm -k 10000,$((1000 * $n))
    echo Running evolution.
    for i in $(seq 1 5)
    do
        mpirun -np $n --map-by numa GameOfLife -r -f imgs/weak_init_$n.pgm -n 1000 -s 0 -e 1 -t $1
    done
done
echo done!