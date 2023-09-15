#!/bin/sh
#SBATCH --output=write_test_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Write_test
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --time=00:10:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

srun -n 1 make timeit

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

mpirun -np $SLURM_NTASKS GameOfLife -i -f imgs/omp_strong_init.pgm -k 10000,10000



touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
for n in $(seq 60 $SLURM_CPUS_PER_TASK)
do
    export OMP_NUM_THREADS=$n
    echo Currently using $n threads.
    for i in $(seq 1 5)
    do
        mpirun -np $SLURM_NTASKS --map-by socket GameOfLife -r -f imgs/omp_strong_init.pgm -n 100 -s 0 -e 1 -t $1
    done
done

squeue -j $SLURM_JOB_ID