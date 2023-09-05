#!/bin/sh
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00

modoule load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

echo Running OMP strong scalability test.
mpirun -np 1 GameOfLife -i -f imgs/omp_strong_init.pgm -k 1000,1000

touch $1
echo "# nprocs, nthreads, total, procrank, thid, comm, grid, idle, write" >> $1
for n in {1..64}
do
    export OMP_NUM_THREADS=$n
    echo Currently using $n threads.
    mpirun -np 1 --map-by socket GameOfLife -r -f imgs/omp_strong_init.pgm -n 5000 -s 0 -e 1 -t $1 >> $1
done
echo done!