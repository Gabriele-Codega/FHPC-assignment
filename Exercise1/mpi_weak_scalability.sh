#!/bin/sh
#SBATCH --output=mpi_weak_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=4
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
echo Running MPI weak scalability test.
touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
for n in $(seq 1 $SLURM_NTASKS)
do
    echo Currently using $n tasks.
    echo Generating grid.
    mpirun -np $n GameOfLife -i -f imgs/weak_init_$n.pgm -k 1000,$((1000 * $n))
    echo Running evolution.
    mpirun -np $n --map-by socket GameOfLife -r -f imgs/weak_init_$n.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!