# GEMM benchmark
This folder contains the code used to perform benchmarks on level 3 BLAS functions from different libraries and the data from said benchmarks.

All the benchmarks were run on ORFEO, a computing cluster at Area Science Park in Trieste, Italy. The available nodes were THIN nodes, which use Intel Xeon Gold processors, and EPYC nodes, which use AMD Epyc processors. Details about the results can be found in the report.

## Compilation
There are two makefiles, which are almost completely equivalent and were duplicated simply for time reasons. To compile the code you will need the required libraries, which are `MKL`, `BLIS` and `OpenBLAS`. If compiling on ORFEO, `MKL` and `OpenBLAS` should already be available, while you will need to setup `BLIS` on your own.

Before compiling you should modify the paths in the Makefile to match your installation of the libraries. If running on ORFEO you will have to load the modules for MKL and OpenBLAS (`module load mkl/latest openBLAS/0.3.23-omp`).

In the makefile, the flag `PFLAG` controls whether the executable will use single or double precision floating point numbers. You can set it to either `-DUSE_DOUBLE` or `-DUSE_SINGLE`. Note that MakefileSingle will produce executables with the `_single` suffix.

Once the makefile is setup, run `make` and compile on the target machine.

## Usage
The executable accepts either 0 or 3 arguments, which determine the size of the matrices to multiply. When providing no arguments, the multiplication will be between a $2000\times 200$ and a $200\times 1000$ matrix. Otherwise, the arguments should be three integers $m,k,n$ and the matrices will be of size $m\times k$ and $k\times n$.

## Scripts
The folder contains some shell scripts, which can be used to schedule slurm jobs.

Before running the scripts, you should make sure that all the parameters for slurm are set correctly and you should change the other parameters to your likings (i.e. size of the matrices, number of measurements per size, number of cores, filename).

The scripts require an argument, which can be `mkl`, `oblas` or `blis`. This determines which library will be used and effectively determines which executable will be ran and the name of the output.

Note that `benchmark.sh` and `benchmark_single.sh` serve the same exact purpose, with the exception that the latter runs executables with `_single` suffix. Both test the performance over variable matrix size.

`scalability.sh` test the preformance over an increasing number of OMP threads.

## Data
There are a number of data files, in csv, which contain the results of the benchmarks. The data are saved as `m,k,n,time,gflops`.

## Files in this folder
- `gemm.c`
    Source code for the benchmark. Initialises memory, calls the gemm function and measures time and flops.
- `Makefile`
    Default makefile, compiles in double precision by default.
- `MakefileSingle`
    Completely equivalent to Makefile, except that it compiles executables with different names and with single precision. It was created only so that I could run single and double benchmarks at the same time.
- `*.sh`
    Shell scripts that run the actual benchmarks. 
- `data/`
    A folder containing the data from benchmarks. The data are organised by node, then by test, then by precision.