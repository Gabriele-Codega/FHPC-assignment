#!/bin/sh
#SBATCH --job-name=gemm_benchmark
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=30

module load architecture/AMD
module load mkl/latest
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/gcodeg00/myblis/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLIS_NUM_THREADS=$SLURM_CPUS_PER_TASK

touch data/$1_benchmark.csv
echo "# m, n, k, time, gflops">>$1_benchmark.csv

for n in $(seq 2000 1000 20000)
do
    srun -n 1 ./gemm_$1.x $n $n $n >> $1_benchmark.csv
done
