#!/bin/sh
#SBATCH --partition=EPYC
#SBATCH --job-name=Game_of_Life
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00

modoule load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

echo Running MPI strong scalability test.
mpirun -np 1 GameOfLife -r -f imgs/init.pgm -n 5000 -s 0 -e 1 -t $1

touch $1
echo "# nproc, nthreads, time" >> $1
export OMP_NUM_THREADS=1
for n in {1..128}
do
    echo Currently using $n tasks.
    mpirun -np $n GameOfLife -r -f imgs/init.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!