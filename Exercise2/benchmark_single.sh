#!/bin/sh
#SBATCH --job-name=gemm_benchmark
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:40:00

module load architecture/AMD
module load mkl/latest
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/gcodeg00/myblis/lib:$LD_LIBRARY_PATH

export BLIS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

declare -a places=("threads" "cores" "numa_domains" "sockets" "ll_caches")

export OMP_PROC_BIND=spread

for place in "${places[@]}"
do 
    export OMP_PLACES=$place

    fname=data/thin/$1_benchmark_single_${OMP_PLACES}_${OMP_PROC_BIND}_numa.csv

    touch $fname
    echo "# running on EPYC" >>$fname
    echo "# thread allocation policy: $OMP_PLACES">>$fname
    echo "# m, n, k, time, gflops">>$fname

    for n in $(seq 2000 1000 20000)
    do
        for i in $(seq 1 5)
        do
            srun -n 1 ./gemm_$1_single.x $n $n $n >> $fname
        done
    done
done

squeue -j $SLURM_JOB_ID