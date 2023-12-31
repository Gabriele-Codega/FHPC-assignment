#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

#include <evolution.h>
#include <io_utils.h>

#define BOT_TO_TOP 1
#define TOP_TO_BOT 2

#ifdef TIMEIT
extern double* time_array;
#endif

void 
ordered_evolution(char* full_grid, char* neigh, const int n, const int s, 
                const int maxval, const int xsize, const int ysize,
                const int procwork, const int procoffset)
{
    /* initialise timing variables */
    #ifdef TIMEIT
    double tstart_comm, total_time_comm=0;
    double tstart_grid, total_time_grid=0;
    double tstart_idle, total_time_idle=0;
    double tstart_write, total_time_write=0;
    double tstart, tend=0;

    tstart = omp_get_wtime();
    #endif

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
    
    MPI_Request snd_bot, rcv_bot, snd_top, rcv_top;
    snd_bot = rcv_bot = snd_top = rcv_top = MPI_REQUEST_NULL;



    /* send initial boundary from last to first process */
    if (procrank == numproc - 1 && thid == numthreads - 1)
    {
        MPI_Isend(bot_boundary - xsize, xsize, MPI_BYTE, next_proc, BOT_TO_TOP, MPI_COMM_WORLD,&snd_bot);
    }
    #ifdef TIMEIT
    tstart_idle = omp_get_wtime();
    #endif
    #pragma omp barrier
    #ifdef TIMEIT
    total_time_idle += omp_get_wtime()-tstart_idle;
    #endif

    if (procrank != 0 && thid == 0)
    {
        MPI_Isend(grid, xsize, MPI_BYTE, prev_proc, TOP_TO_BOT, MPI_COMM_WORLD, &snd_top);
    }
    #ifdef TIMEIT
    tstart_idle = omp_get_wtime();
    #endif
    #pragma omp barrier
    #ifdef TIMEIT
    total_time_idle += omp_get_wtime()-tstart_idle;
    #endif

    /* ordered evolution loop */
    for (int step = 1; step <=n ; step++){

        #pragma omp master
        {
            #ifdef TIMEIT
            tstart_comm = omp_get_wtime();
            #endif
            MPI_Irecv(top_boundary,xsize,MPI_BYTE,prev_proc,BOT_TO_TOP,MPI_COMM_WORLD,&rcv_top);
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

        if (thid == numthreads - 1)
        {
            #ifdef TIMEIT
            tstart_comm = omp_get_wtime();
            #endif
            MPI_Irecv(bot_boundary, xsize, MPI_BYTE, next_proc, TOP_TO_BOT, MPI_COMM_WORLD, &rcv_bot);
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

        #pragma omp master
        {
            #ifdef TIMEIT
            tstart_comm = omp_get_wtime();
            #endif
            MPI_Wait(&rcv_top,MPI_STATUS_IGNORE);
            if(snd_top != MPI_REQUEST_NULL)
                MPI_Wait(&snd_top,MPI_STATUS_IGNORE);
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
        

        #pragma omp for schedule(static) ordered 
        for (int i = 0; i<procwork/xsize;i++){
            #ifdef TIMEIT
            tstart_idle = omp_get_wtime();
            #endif
            #pragma omp ordered
            {
            #ifdef TIMEIT
            total_time_idle += omp_get_wtime()-tstart_idle;
            #endif
            int irow = i *xsize; 
            if (thid == numthreads - 1 && i >= 1)
            {
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Wait(&rcv_bot,MPI_STATUS_IGNORE);
                #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }

            #ifdef TIMEIT
            tstart_grid = omp_get_wtime();
            #endif
            for (int j = 0;j<xsize;j++){
                int nneigh = 0;

                for(int ii = i-1; ii < i+2; ii++){
                    int iirow = ii * xsize;
                    for (int jj = j+xsize-1; jj < j+xsize+2; jj++){
                        int jjcol = jj%xsize;
                        nneigh += grid[iirow+jjcol] * (!(iirow == irow && jjcol == j));
                    }
                }
                /* Uncomment the following line to use standard rules */
                // grid[irow + j] = grid[irow+j] * (nneigh == 2) + (nneigh == 3);
                grid[irow + j] = (nneigh == 2) + (nneigh == 3);
            }
            #ifdef TIMEIT
            total_time_grid += omp_get_wtime()-tstart_grid;
            #endif

            if(thid == 0 && i == 0)
            {
                #ifdef TIMEIT
                tstart_comm = omp_get_wtime();
                #endif
                MPI_Isend(grid,xsize,MPI_BYTE,prev_proc,TOP_TO_BOT,MPI_COMM_WORLD,&snd_top);
                 #ifdef TIMEIT
                total_time_comm += omp_get_wtime()-tstart_comm;
                #endif
            }
            }
        }

        if (thid == numthreads - 1){
            #ifdef TIMEIT
            tstart_comm = omp_get_wtime();
            #endif
            MPI_Isend(bot_boundary-xsize, xsize, MPI_BYTE, next_proc, BOT_TO_TOP, MPI_COMM_WORLD,&snd_bot);
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



        /* write a checkpoint if required */
        if (step%s == 0){    
            #ifdef TIMEIT
            tstart_write = omp_get_wtime();
            #endif        
            write_checkpoint(cp_fname,step,grid,procrank,procoffset,procwork,xsize,ysize,maxval);
            #ifdef TIMEIT
            total_time_write += omp_get_wtime()-tstart_write;
            #endif
        }

    } /* evolution step */
    free(cp_fname);
    if (thid == numthreads - 1)
    {
        if(snd_bot != MPI_REQUEST_NULL){
            MPI_Cancel(&snd_bot);
            MPI_Request_free(&snd_bot);
        }

        if(rcv_bot != MPI_REQUEST_NULL){
            MPI_Cancel(&rcv_bot); 
            MPI_Request_free(&rcv_bot);
        }

    }
    #pragma omp barrier
    #pragma omp master
    {
        if(snd_top != MPI_REQUEST_NULL){
            MPI_Cancel(&snd_top);
            MPI_Request_free(&snd_top);
        }
        if(rcv_top != MPI_REQUEST_NULL){
            MPI_Cancel(&rcv_top); 
            MPI_Request_free(&rcv_top);
        }
    }
    #pragma omp barrier
    /* stop timer and write measurements in array */
    #ifdef TIMEIT
    tend = omp_get_wtime()-tstart;

    time_array[thid] = tend;
    time_array[numthreads + thid] = total_time_comm;
    time_array[2*numthreads + thid] = total_time_grid;
    time_array[3*numthreads + thid] = total_time_idle;
    time_array[4*numthreads + thid] = total_time_write;
    #endif
}