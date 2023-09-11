#!/bin/sh
#SBATCH --job-name=gemm_benchmark
#SBATCH --partition=EPYC

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00

module load architecture/AMD
module load mkl/latest
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/gcodeg00/myblis/lib:$LD_LIBRARY_PATH

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

fname=data/epyc/$1_scalability_single.csv

touch $fname
echo "# running on EPYC" >>$fname
echo "# thread allocation policy: $OMP_PLACES, $OMP_PROC_BIND">>$fname
echo "# max number of threads: $SLURM_CPUS_PER_TASK">>$fname
echo "# m, n, k, time, gflops">>$fname

size=10000

for n in $(seq 1 $SLURM_CPUS_PER_TASK)
do 
    export BLIS_NUM_THREADS=$n
    export OMP_NUM_THREADS=$n
    srun -n 1 gemm_$1_single.x $size $size $size >>$fname
done

squeue -j $SLURM_JOB_ID
