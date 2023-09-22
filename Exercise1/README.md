# Game of Life
The goal is to implement Conway's Game of Life (or a modification of that) in a hybrid MPI - OpenMP code.

The code is written in C and uses OpenMPI. More details about the implementation can be found in the report.

## Compilation
To compile the code you will need OpenMP and OpenMPI. If you are using ORFEO you can `module load OpenMPI/4.1.5/gnu/12.2.1` (or whichever the current version is).

To compile you can simply use `make` or `make timeit`. The latter compiles the code including calls to OpenMP's timing routine `omp_get_wtime()`, which are used to measure the time to perform certain operations.

## Usage
To execute the code, use `mpirun` and supply some command-line arguments:
- `-i` : initialise a playground
- `-r` : run a given playground (read from file)
- `-k x,y` : specify the x and y dimensions of the playground, separated by a comma
- `-f FILE_NAME` : either read or write from/to FILE (should be a .pgm P5 file)
- `-n INT` : number of evolution steps
- `-s INT` : frequency of checkpoints (if 0, only produces a checkpoint on the last iteration)
- `-e INT` : type of evolution scheme to follow
    - 0 : ordered evolution
    - 1 : static evolution
- `-t FILE_NAME` : **when compiled with `timeit`**, specify the file, in csv, on which the measurements should be written.

## Scripts
The scripts are meant to be used when code is compiled in `timeit` mode and are used to perform the scalability tests.

Before running a script you should set the parameters (slurm parameters, problem size, number of iterations, kind of evolution, number of measurements etc.)

The scripts require an argument, which is the name of the file on which data should be saved.

## Data
The data are saved in csv format for easy use.

The header always begins with '#' because at an early stage I used gnuplot to visualise the results and lines that begin with '#' are treated as comments. Sometimes there is more than one line in the header, to take note of some other parameters used for the run.

Data are stored as `nprocs, nthreads, total, comm, grid, idle, write`, where
- `nprocs` is the number of MPI tasks
- `nthreads` is the number of OMP threads
- `total` is the total time for the evolution
- `comm` is the time spent in communications
- `grid` is the time spent updating the grid
- `idle` is the time spent doing nothing/waiting for other threads
- `write` is the time spent to write on file

Note that all these quantities refer to the maximum time among all the threads, which is in fact the limiting factor to the overall speed.


## Files in this directory
- `src/*.c`   
    Source files, including the main, the evolution routines and some i/o utilities
- `inc/*.h`     
    Headers with function declarations and definitions
- `imgs/*.pgm`      
    Initial conditions and checkpoints
- `Makefile`
- `*.sh`       
    Shell scripts to run some scalability tests.
- `obj/`        
    Directory where object files will be placed at compile time
- `data/`
    A folder containing the results from various scalability tests, ran with different parameters. 