#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

#include <evolution.h>
#include <io_utils.h>

#define BOT_TO_TOP 1
#define TOP_TO_BOT 2

void 
static_evolution(char* full_grid, char* neigh, const int n, const int s, 
                const int maxval, const int xsize, const int ysize,
                const int procwork, const int procoffset,
                const int thwork, const int thoffset)
{
    /* get info about process and thread */
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

    /* repoint the grid for convenience */
    char* top_boundary = full_grid;
    char* grid = full_grid + xsize;
    char* bot_boundary = grid+procwork;
    
    /* name of file to write checkpoints */
    char* cp_fname = (char*)malloc(26*sizeof(char));

    MPI_Request rcv_top,rcv_bot;

    #ifdef TIMEIT
    double tstart_comm, total_time_comm=0;
    double tstart_grid, total_time_grid=0;
    double tstart_idle, total_time_idle=0;
    double tstart_write, total_time_write=0;
    double tstart, tend=0;
    #endif
    

    #ifdef TIMEIT
    tstart = omp_get_wtime();
    #endif
    /* static evolution loop */
    for (int step = 1; step <=n ; step++){

        /* communicating boundaries */
        if (procrank%2==0){
            if (thid==0){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Irecv(top_boundary, xsize, MPI_BYTE, prev_proc ,BOT_TO_TOP,MPI_COMM_WORLD, &rcv_top);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif

            if (thid == numthreads-1){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Irecv(bot_boundary, xsize, MPI_BYTE, next_proc,TOP_TO_BOT,MPI_COMM_WORLD, &rcv_bot);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif

            if (thid==0){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Send(grid, xsize, MPI_BYTE, prev_proc, TOP_TO_BOT, MPI_COMM_WORLD);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif
        
            if (thid == numthreads-1){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Send(bot_boundary-xsize, xsize, MPI_BYTE, next_proc, BOT_TO_TOP, MPI_COMM_WORLD);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif

        } else {
            if (thid==0){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Send(grid, xsize, MPI_BYTE, prev_proc, TOP_TO_BOT, MPI_COMM_WORLD);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif
        
            if (thid == numthreads-1){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Send(bot_boundary-xsize, xsize, MPI_BYTE, next_proc, BOT_TO_TOP, MPI_COMM_WORLD);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif

            if (thid==0){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Irecv(top_boundary, xsize, MPI_BYTE, prev_proc, BOT_TO_TOP,MPI_COMM_WORLD, &rcv_top);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif

            if (thid == numthreads-1){
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Irecv(bot_boundary, xsize, MPI_BYTE, next_proc,TOP_TO_BOT,MPI_COMM_WORLD, &rcv_bot);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp barrier
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif

        }

        /* wait for top boundary */
        if (thid == 0){
            #ifdef TIMEIT
            tstart_comm = omp_get_wtime();
            #endif
            MPI_Wait(&rcv_top,MPI_STATUS_IGNORE);
            #ifdef TIMEIT
            total_time_comm += omp_get_wtime()-tstart_comm;
            #endif
        }
        #ifdef TIMEIT
        tstart_idle = omp_get_wtime();
        #endif
        #pragma omp barrier
        #ifdef TIMEIT
        total_time_idle += omp_get_wtime()-tstart_idle;
        #endif

        /* wait for bottom boundary */
        if (thid == numthreads - 1){
            #ifdef TIMEIT
            tstart_comm = omp_get_wtime();
            #endif
            MPI_Wait(&rcv_bot,MPI_STATUS_IGNORE);
            #ifdef TIMEIT
            total_time_comm += omp_get_wtime()-tstart_comm;
            #endif
        }
        #ifdef TIMEIT
        tstart_idle = omp_get_wtime();
        #endif
        #pragma omp barrier
        #ifdef TIMEIT
        total_time_idle += omp_get_wtime()-tstart_idle;
        #endif


        
        #ifdef TIMEIT
        tstart_grid = omp_get_wtime();
        #endif
        /* compute number of live neighbours */
        #pragma omp for schedule(static) 
        for (int i = 0; i<procwork/xsize;i++){
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

                neigh[irow + j] = nneigh;
            }
        }

        /* evolve according to neighbours */
        #pragma omp for schedule(static)
        for (int i = 0; i<procwork/xsize;i++){
            int irow = i *xsize; 
            for (int j = 0;j<xsize;j++){
                grid[irow + j] = grid[irow+j] * (neigh[irow+j] == 2) + (neigh[irow + j] == 3);
            }
        }
        #ifdef TIMEIT
        total_time_grid += omp_get_wtime()-tstart_grid;
        #endif


        
        /* write a checkpoint if required */
        if (step%s == 0){
            #ifdef TIMEIT
            tstart_write = omp_get_wtime();
            #endif
            write_checkpoint(cp_fname,step,grid,procrank,procoffset,thoffset,thwork,xsize,ysize,maxval);     
            #ifdef TIMEIT
            total_time_write += omp_get_wtime()-tstart_write;
            #endif
        }

    }/* evolution step */
    #ifdef TIMEIT
    tend = omp_get_wtime()-tstart;
    #endif

    free(cp_fname);
    #ifdef TIMEIT
    // printf("# nprocs, nthreads, total, procrank, thid, comm, grid, idle, write\n")
    printf("%d,%d,%f,%d,%d,%f,%f,%f,%f\n",numproc,numthreads,tend, procrank, thid, total_time_comm, total_time_grid, total_time_idle, total_time_write);
    //printf("(p: %d, t: %d) Comm: %f, Grid: %f, Idle: %f, Write: %f\n", procrank, thid, total_time_comm, total_time_grid, total_time_idle, total_time_write);
    #endif

}