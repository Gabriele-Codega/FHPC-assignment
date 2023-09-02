# Game of Life
The goal is to implement Conway's Game of Life (or a modification of that) in a hybrid MPI - OpenMP code.

The code is written in C and uses OpenMPI.

**Current version is not the final version**

To compile run `make` or `make timeit`. The latter compiles the code including calls to OpenMP's timing routine `omp_get_wtime()`, which are used to measure the time spent in communications, operations on the playground and idle for each thread in each MPI process and an overall measure of the time spent by each MPI process in the evolution routine.

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


## Files in this directory
- `src/*.c`   
    Source files, including the main, the evolution routines and some i/o utilities
- `inc/*.h`     
    Headers with function declarations and definitions
- `imgs/*.pgm`      
    Initial conditions and checkpoints
- `Makefile`
- `timing.sh`       
    A shell script to run an OMP strong scalability test
- `obj/`        
    Directory where object files will be placed at compile time