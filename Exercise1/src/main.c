#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <evolution.h>
#include <pgm_utils.h>

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1



int main(int argc, char** argv){

    char fname_deflt[] = "game_of_life.pgm";

    int   action = 0;
    int   k      = K_DFLT;
    int   e      = ORDERED;
    int   n      = 10000;
    int   s      = 1;
    char *fname  = NULL;
    char *optstring = "irk:e:f:n:s:";

    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        
        case 'i':
        action = INIT; break;
        
        case 'r':
        action = RUN; break;
        
        case 'k':
        k = atoi(optarg); break;

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
        char* grid;
        grid = (char*) malloc(k*k*sizeof(char));

        for (int i = 0; i<k;i++){
            int irow = i *k; 
            for (int j = 0;j<k;j++){
                grid[irow + j] = 0;
            }
        }
        grid[5 * k + 4] = 1;
        grid[5 * k + 5] = 1;
        grid[5 * k + 6] = 1;

        void* ptr = (void*) grid; 
        int written = write_pgm_image(ptr,1,k,k,fname);
    } else if (action == RUN) {

        /* read initial condition from file */

        void* ptr;
        int xsize = 0;
        int ysize = 0;
        int maxval = 0;

        read_pgm_image(&ptr,&maxval,&xsize,&ysize,fname);
        char* grid = (char*) ptr;
        char* neigh = NULL;

        register void (*evolution)(char*, char*, const int, const int);

        switch(e){
            case ORDERED:
                evolution = ordered_evolution;
                break;

            case STATIC:
                evolution = static_evolution;
                break;

            default:
                printf("Evolution type %d unknown.\nValid values are 0 (ordered), 1 (static).",e);
                break;

        }/* switch */


        evolution(grid,neigh,n,s);



    }

    return 0;
}