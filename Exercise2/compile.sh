#!/bin/sh
module load architecture/AMD
module load mkl/latest
module load openBLAS/0.3.23-omp

srun -n 1 -p EPYC make