#!/bin/sh
echo Running MPI weak scalability test.
touch $1
echo "# nproc, nthreads, time" >> $1
export OMP_NUM_THREADS=2
for n in {1..4}
do
    echo Currently using $n tasks.
    echo Generating grid.
    mpirun -np $n GameOfLife -i -f imgs/weak_init_$n.pgm -k 500,$((100 * $n))
    echo Running evolution.
    mpirun -np $n GameOfLife -r -f imgs/weak_init_$n.pgm -n 5000 -s 0 -e 1 -t $1
done
echo done!