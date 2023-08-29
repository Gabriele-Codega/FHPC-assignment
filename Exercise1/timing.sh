#!/bin/sh
rm timings.csv
touch timings.csv
echo "# nproc, nthreads, time\n" >> timings.csv
for n in {1..8}
do
    export OMP_NUM_THREADS=$n
    mpirun -np 1 GameOfLife -r -f imgs/init.pgm -n 10000 -s 0 -e 1
done
echo done!