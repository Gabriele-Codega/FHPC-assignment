#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include <evolution.h>
#include <pgm_utils.h>

#define TCPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec +	\
		   (double)ts.tv_nsec * 1e-9)

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1



int main(int argc, char** argv){

    char fname_deflt[] = "game_of_life.pgm";

    int   action = 0;
    int   xsize  = K_DFLT;
    int   ysize  = K_DFLT;
    int   e      = ORDERED;
    int   n      = 10000;
    int   s      = 1;
    char *fname  = NULL;
    char *optstring = "irk:e:f:n:s:";

    int c;
    char* xval = NULL;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        
        case 'i':
        action = INIT; break;
        
        case 'r':
        action = RUN; break;
        
        case 'k':
        xval = (char*)malloc( sizeof(optarg)+1 );
        xval = strsep(&optarg,",");
        xsize = atoi(xval);
        ysize = atoi(optarg);
        break;

        case 'e':
        e = atoi(optarg); break;

        case 'f':
        fname = (char*)malloc( sizeof(optarg)+1 );
        sprintf(fname, "%s", optarg );
        break;

        case 'n':
        n = atoi(optarg); break;

        case 's':
        s = atoi(optarg); break;

        default :
        printf("argument -%c not known\n", c ); break;
        }
    }
    



    if (action == INIT){
        char* grid = (char*) malloc(xsize*ysize*sizeof(char));

        for (int i = 0; i<ysize;i++){
            int irow = i *xsize; 
            for (int j = 0;j<xsize;j++){
                grid[irow + j] = 0;
            }
        }
        grid[5 * xsize + 2] = 1;
        grid[5 * xsize + 3] = 1;
        grid[5 * xsize + 4] = 1;

        write_pgm_image((void*)grid,1,xsize,ysize,fname);
    } else if (action == RUN) {

        /* read initial condition from file */

        void* ptr = NULL; // redundant??
        xsize = 0;
        ysize = 0;
        int maxval = 0;

        read_pgm_image(&ptr,&maxval,&xsize,&ysize,fname);
        char* grid = (char*) ptr;
        char* neigh = NULL;

        void (*evolution)(char*, char*, const int, const int, const int, const int, const int) = NULL;

        struct timespec ts;
        double tstart,tend;

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

        tstart = TCPU_TIME;
        (*evolution)(grid,neigh,n,s,maxval,xsize,ysize);
        tend = TCPU_TIME-tstart;
        printf("Elapsed time: %f s\n",tend);


    }

    return 0;
}