#!/bin/sh
#SBATCH --output=mpi_weak_datagen_%j.out
#SBATCH --partition=EPYC
#SBATCH --job-name=Weak_data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:10:00

module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

echo Generating data for weak scalability. xsize = 10000, ysize = 1000 x n.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
for n in $(seq 1 32)
do
    echo Current size is $n.
    mpirun -np 1 GameOfLife -i -f imgs/weak_init_$n.pgm -k 10000,$((1000 * $n))
done
echo done!