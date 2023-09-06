#!/bin/sh
#SBATCH --output=omp_strong_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

echo Running OMP strong scalability test.
mpirun -np $SLURM_NTASKS GameOfLife -i -f imgs/omp_strong_init.pgm -k 1000,1000

touch $1
echo "# nprocs, nthreads, total, comm, grid, idle, write" >> $1
for n in $(seq 1 $SLURM_CPUS_PER_TASK)
do
    export OMP_NUM_THREADS=$n
    echo Currently using $n threads.
    mpirun -np $SLURM_NTASKS --map-by socket GameOfLife -r -f imgs/omp_strong_init.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!