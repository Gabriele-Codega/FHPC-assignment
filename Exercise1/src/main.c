#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>

#include <evolution.h>
#include <io_utils.h>

#define INIT 1
#define RUN  2
#define K_DFLT 100
#define ORDERED 0
#define STATIC  1

#ifdef TIMEIT
double* time_array = NULL;
#endif


int main(int argc, char** argv){

    char fname_deflt[] = "game_of_life.pgm";

    int   action = 0;
    int   xsize  = K_DFLT;
    int   ysize  = K_DFLT;
    int   e      = ORDERED;
    int   n      = 10000;
    int   s      = 1;
    char *fname  = NULL;

    #ifdef TIMEIT
    char *optstring = "irk:e:f:n:s:t:";
    char* timefile = NULL;
    #else
    char *optstring = "irk:e:f:n:s:";
    #endif

    /* read input from user */
    int c;
    char* xval = NULL;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        
        case 'i':
        action = INIT; break;
        
        case 'r':
        action = RUN; break;
        
        case 'k':
        xval = (char*)malloc( strlen(optarg)+1 );
        xval = strsep(&optarg,",");
        xsize = atoi(xval);
        ysize = atoi(optarg);
        break;

        case 'e':
        e = atoi(optarg); break;

        case 'f':
        fname = (char*)malloc( strlen(optarg)+1 );
        sprintf(fname, "%s", optarg );
        break;

        case 'n':
        n = atoi(optarg); break;

        case 's':
        s = atoi(optarg); break;

        #ifdef TIMEIT
        case 't':
        timefile = (char*)malloc( strlen(optarg)+1 );
        sprintf(timefile, "%s", optarg );
        break;
        #endif

        default :
        printf("Argument -%c not known\n", c ); break;
        }
    }
    if (s==0) s = n;

    


    /* choose whether to initialise a playground or run */
    if (action == INIT){
        int provided;
        MPI_Init_thread(NULL,NULL,MPI_THREAD_SERIALIZED,&provided);
        if(provided < MPI_THREAD_SERIALIZED){
            printf("Provided thread level lower than required\n");
            fflush(stdout);
            MPI_Finalize();
            exit(1);
        }

        char* grid = NULL;

        #pragma omp parallel shared(grid)
        {
            int numproc;
            int procrank;
            #pragma omp single copyprivate(numproc,procrank)
            {
                MPI_Comm_size(MPI_COMM_WORLD,&numproc);
                MPI_Comm_rank(MPI_COMM_WORLD,&procrank);
            }
            int thid = omp_get_thread_num();


            int procwork ;
            MPI_Offset procoffset;

            procwork = ysize/numproc + (procrank< (ysize%numproc));
            procoffset = ysize/numproc * procrank + (procrank >= (ysize%numproc)) * (ysize%numproc) + (procrank < (ysize%numproc)) * procrank;
            procwork *= xsize;
            procoffset *= xsize;

            #pragma omp single
            {
                grid = (char*) malloc(procwork*sizeof(char));
                if (grid == NULL )
                {
                    printf("Could not allocate a grid. Aborting.\n");
                    fflush(stdout);
                    MPI_Finalize();
                    exit(3);
                }
            }
        
            long int seed = time(NULL) + thid + procrank;
            srand48(seed);
            #pragma omp for
                for (int i = 0; i<procwork/xsize;i++){
                    int irow = i *xsize; 
                    for (int j = 0;j<xsize;j++){
                        grid[irow + j] = (int)(drand48()*2);
                    }
                }

            MPI_File fhout;
            int writeoffset;

            #pragma omp single copyprivate(fhout,writeoffset)
            {
                if(procrank == 0){
                    write_header(fname,xsize,ysize,1,&writeoffset);
                }
                MPI_Barrier(MPI_COMM_WORLD);

                MPI_Bcast(&writeoffset, 1, MPI_INT, 0, MPI_COMM_WORLD);

                MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY,MPI_INFO_NULL,&fhout);
                MPI_Offset start;
                MPI_File_get_position_shared(fhout,&start);
                MPI_Barrier(MPI_COMM_WORLD);

                writeoffset += start;
            }

            MPI_Status status;
            #pragma omp master
                MPI_File_write_at_all(fhout,writeoffset+procoffset, (void*)(grid),procwork,MPI_BYTE,&status);
            #pragma omp barrier
            #pragma omp single
                MPI_File_close(&fhout);
        }
        MPI_Finalize();
        free(grid);
    } else if (action == RUN) {

        /* select correct evolution routine */
        void (*evolution)(char* , char* , const int , const int , 
                        const int , const int , const int ,
                        const int , const int ) = NULL;
        switch (e)
        {
        case ORDERED:
            evolution = ordered_evolution;
            break;
        case STATIC:
            evolution = static_evolution;
            break;
        default:
            printf("Unknown evolution type '%d'.\n",e);
            printf("Valid values are:\n\t- 0 for ordered\n\t- 1 for static\n");
            return 1;
        }

        /* initialise MPI */
        int provided;
        MPI_Init_thread(NULL,NULL,MPI_THREAD_SERIALIZED,&provided);
        if(provided < MPI_THREAD_SERIALIZED){
            printf("Provided thread level lower than required\n");
            fflush(stdout);
            MPI_Finalize();
            exit(1);
        }

        /* get information about process */
        int numproc;
        int procrank;
        MPI_Comm_size(MPI_COMM_WORLD,&numproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&procrank);
        int params[4];

        /* read the header of input file */
        if (procrank==0){    
            read_header(fname, params);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(params, 4, MPI_INT, 0, MPI_COMM_WORLD);
        int xsize = params[0];
        int ysize = params[1];
        int maxval = params[2];
        int globaloffset = params[3];

        /* open a shared file and get the starting point for the data */
        MPI_File fh;
        MPI_Offset filestart;
        MPI_File_open(MPI_COMM_WORLD, fname,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
        MPI_File_get_position_shared(fh, &filestart);
        MPI_Barrier(MPI_COMM_WORLD);
        globaloffset += filestart;
        
        /* figure out workload and offset for each process (in rows) */

        static int procwork ;
        static MPI_Offset procoffset;
        static int thwork;
        static int thoffset;
        #pragma omp threadprivate(procwork,procoffset,thwork,thoffset)

        procwork = ysize/numproc + (procrank< (ysize%numproc));
        procoffset = ysize/numproc * procrank + (procrank >= (ysize%numproc)) * (ysize%numproc) + (procrank < (ysize%numproc)) * procrank;

        char* mygrid = (char*)malloc(procwork*xsize + 2*xsize);
        char* myneigh = NULL; 
        if(e == STATIC){
            myneigh = (char*)malloc(procwork*xsize);
            if(!myneigh){
                printf("Could not allcoate neighbours. Aborting.\n");
                MPI_Finalize();
                exit(3);
            }
        }
        if(!mygrid){
            printf("Could not allcoate grid. Aborting.\n");
            MPI_Finalize();
            exit(3);
        }

        #pragma omp parallel copyin(procwork, procoffset)
        {
            int numthreads = omp_get_num_threads();
            int thid = omp_get_thread_num();
            thwork = procwork/numthreads + (thid< (procwork%numthreads));
            thoffset = procwork/numthreads * thid + (thid >= (procwork%numthreads)) * (procwork%numthreads) + (thid < (procwork%numthreads)) * thid;
            thwork *= xsize;
            thoffset *= xsize;
            procwork *= xsize;
            procoffset *= xsize;

            MPI_Status status;
            #pragma omp critical
                MPI_File_read_at_all(fh,globaloffset+procoffset+thoffset, (void*)(mygrid+xsize+thoffset),thwork,MPI_BYTE,&status);
            
            #pragma omp barrier
        }
        MPI_File_close(&fh);

        #pragma omp parallel
        {
            #ifdef TIMEIT
            #pragma omp single
                time_array = (double*)malloc(5*omp_get_num_threads() * sizeof(double));
            #endif
            (*evolution)(mygrid,myneigh,n,s,
                        maxval,xsize,ysize,
                        procwork, procoffset);
        }

        /* process the measurements, extracting maximum times and writing on file */
        #ifdef TIMEIT
        double* max = (double*)calloc(5,sizeof(double));
        int numthreads;
        #pragma omp parallel
            #pragma omp master
            numthreads = omp_get_num_threads();

        for (int i = 0;i<5;i++){
            int ii = i * numthreads;
            for (int j = 0; j < numthreads; j++){
                if (time_array[ii + j] > max[i]){
                    max[i] = time_array[ii + j];
                }
            }
        }

        double* local_maxs = NULL;
        double* global_max = NULL;
        if(procrank == 0){
            local_maxs = (double*)malloc(5 * numproc * sizeof(double));
            global_max = (double*)calloc(5,sizeof(double));
        }
        MPI_Gather(max, 5, MPI_DOUBLE, local_maxs, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (procrank == 0){
            for (int j = 0; j < 5; j++){
                for (int i = 0;i<numproc;i++){
                    if (local_maxs[i * 5 + j] > global_max[j]){
                        global_max[j] = local_maxs[i * 5 + j];
                    }
                }
            }
            
            FILE* timefh;
            if( (timefh = fopen(timefile,"a")) )
            {
                fprintf(timefh,"%d,%d,%f,%f,%f,%f,%f\n",numproc,numthreads,global_max[0],global_max[1],global_max[2],global_max[3],global_max[4]);
                fclose(timefh);
            } else {
                printf("Couldn't write times.\n");
            }
        }

        free(time_array);
        free(max);
        free(local_maxs);
        free(global_max);
        #endif

        free(mygrid);
        free(myneigh);
        MPI_Finalize();
    }

    return 0;
}