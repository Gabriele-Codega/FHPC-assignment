#!/bin/sh
echo Running MPI strong scalability test.
touch $1
echo "# nproc, nthreads, time" >> $1
export OMP_NUM_THREADS=1
for n in {1..4}
do
    echo Currently using $n tasks.
    mpirun -np $n GameOfLife -r -f imgs/init.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!