#!/bin/sh
#SBATCH --output=omp_strong_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:01:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

size=5000
nsteps=50

echo Running OMP strong scalability test, size $size, $nsteps steps, $SLURM_NTASKS tasks.
#mpirun -np $SLURM_NTASKS GameOfLife -i -f imgs/omp_strong_init.pgm -k $size,$size


touch $1
echo "# ordered evolution: size = $size, steps = $nsteps" >> $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
# for n in $(seq 1 $SLURM_CPUS_PER_TASK)
# do
#     export OMP_NUM_THREADS=$n
#     echo Currently using $n threads.
#     for i in $(seq 1 5)
#     do
        mpirun -np $SLURM_NTASKS --map-by socket GameOfLife -r -f imgs/omp_strong_init.pgm -n $nsteps -s 0 -e 0 -t $1
#     done
# done

squeue -j $SLURM_JOB_ID