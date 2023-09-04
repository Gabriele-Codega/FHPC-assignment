#!/bin/sh
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00

echo Running OMP strong scalability test.
touch $1
echo "# nproc, nthreads, time" >> $1
for n in {1..8}
do
    export OMP_NUM_THREADS=$n
    echo Currently using $n threads.
    mpirun -np 1 --map-by socket GameOfLife -r -f imgs/init.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!