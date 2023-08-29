#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

#include <evolution.h>
#include <io_utils.h>

#define BOT_TO_TOP 1
#define TOP_TO_BOT 2


void 
ordered_evolution(char* full_grid, char* neigh, const int n, const int s, 
                const int maxval, const int xsize, const int ysize,
                const int procwork, const int procoffset,
                const int thwork, const int thoffset)
{
    /* figure out ranks and ids */
    int procrank;
    int numproc;
    #pragma omp single copyprivate(procrank, numproc)
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&procrank);
        MPI_Comm_size(MPI_COMM_WORLD,&numproc);
    }

    int thid = omp_get_thread_num();
    int numthreads = omp_get_num_threads();

    int prev_proc = (procrank-1+numproc)%numproc;
    int next_proc = (procrank+1)%numproc;

    /* repoint grid for convenience */
    char* top_boundary = full_grid;
    char* grid = full_grid + xsize;
    char* bot_boundary = grid+procwork;
    
    
    char* cp_fname = (char*)malloc(26*sizeof(char));
    if (cp_fname==NULL){
        printf("Unable to allocate a filename");
        exit(2);
    }
    
    MPI_Status status;
    MPI_Request snd_bot, rcv_bot;

    #ifdef TIMEIT
    double tstart_comm, total_time_comm=0;
    double tstart_grid, total_time_grid=0;
    #endif

    /* send initial boundary from last to first process */
    if (procrank == numproc - 1 && thid == numthreads - 1)
    {
        MPI_Isend(bot_boundary - xsize, xsize, MPI_BYTE, next_proc, BOT_TO_TOP, MPI_COMM_WORLD,&snd_bot);
    }

    /* ordered evolution loop */
    for (int step = 1; step <=n ; step++){
        
        #ifdef TIMEIT
        // #pragma omp single copyprivate(tstart_comm,total_time_comm)
        tstart_comm = omp_get_wtime();
        #endif
        if (thid == numthreads - 1)
        {
            MPI_Irecv(bot_boundary, xsize, MPI_BYTE, next_proc, TOP_TO_BOT, MPI_COMM_WORLD, &rcv_bot);
        }
        #pragma omp barrier
        if (thid == 0)
        {
            MPI_Send(grid,xsize, MPI_BYTE,prev_proc,TOP_TO_BOT,MPI_COMM_WORLD);
        }
        #pragma omp barrier
        if (thid == numthreads - 1)
            MPI_Wait(&rcv_bot,MPI_STATUS_IGNORE);
        #pragma omp barrier

        #pragma omp master
        {
            MPI_Recv(top_boundary,xsize,MPI_BYTE,prev_proc,BOT_TO_TOP,MPI_COMM_WORLD,&status);
        }
        #pragma omp barrier
        #ifdef TIMEIT
        // #pragma omp single copyprivate(tstart_comm,total_time_comm)
        total_time_comm += omp_get_wtime()-tstart_comm;
        #endif

        #ifdef TIMEIT
        // #pragma omp single copyprivate(tstart_grid,total_time_grid)
        tstart_grid = omp_get_wtime();
        #endif
        #pragma omp for schedule(static) ordered 
        for (int i = 0; i<procwork/xsize;i++){
            #pragma omp ordered
            {
            int irow = i *xsize; 
            for (int j = 0;j<xsize;j++){
                int nneigh = 0;

                for(int ii = i-1; ii < i+2; ii++){
                    int iirow = ii * xsize;
                    for (int jj = j+xsize-1; jj < j+xsize+2; jj++){
                        int jjcol = jj%xsize;
                        nneigh += grid[iirow+jjcol] * (!(iirow == irow && jjcol == j));
                    }
                }

                grid[irow + j] = grid[irow+j] * (nneigh == 2) + (nneigh == 3);
            }
            }
        }
        #ifdef TIMEIT
        // #pragma omp single copyprivate(tstart_grid,total_time_grid)
        total_time_grid += omp_get_wtime()-tstart_grid;
        #endif


        #ifdef TIMEIT
        // #pragma omp single copyprivate(tstart_comm,total_time_comm)
        tstart_comm = omp_get_wtime();
        #endif
        if (thid == numthreads - 1)
        {
            MPI_Isend(bot_boundary - xsize, xsize, MPI_BYTE, next_proc, BOT_TO_TOP, MPI_COMM_WORLD,&snd_bot);
            //MPI_Barrier(MPI_COMM_WORLD);
        }
        #pragma omp barrier
        #ifdef TIMEIT
        // #pragma omp single copyprivate(tstart_comm,total_time_comm)
        total_time_comm += omp_get_wtime()-tstart_comm;
        #endif


        /* write a checkpoint if required */
        if (step%s == 0){
            write_checkpoint(cp_fname,step,grid,procrank,procwork,thoffset,thwork,xsize,ysize,maxval);
        }

    } /* evolution step */
    free(cp_fname);
    if (thid == numthreads - 1)
        MPI_Cancel(&snd_bot);
    #ifdef TIMEIT
    printf("Comm: %f, Grid: %f\n", total_time_comm, total_time_grid);
    #endif
}