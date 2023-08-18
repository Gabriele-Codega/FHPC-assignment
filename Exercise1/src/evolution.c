#include <stdlib.h>
#include <stdio.h>

#include <evolution.h>
#include <pgm_utils.h>

void ordered_evolution(char* grid, char* neigh, const int n, const int s, const int maxval, const int xsize, const int ysize){
    char* cp_fname = (char*)malloc(21*sizeof(char));
    if (cp_fname==NULL){
        printf("unable to allocate a filename");
        exit(2);
    }
    
    /* ordered evolution loop */
    for (int step = 1; step <=n ; step++){
        
        for (int i = 0; i<ysize;i++){
            int irow = i *xsize; 
            for (int j = 0;j<xsize;j++){
                int nneigh = 0;

                for(int ii = i+ysize-1; ii < i+ysize+2; ii++){
                    int iirow = (ii%ysize) * xsize;
                    for (int jj = j+xsize-1; jj < j+xsize+2; jj++){
                        int jjcol = jj%xsize;
                        nneigh += grid[iirow+jjcol] * (!(iirow == irow && jjcol == j));
                    }
                }

                grid[irow + j] = grid[irow+j] * (nneigh == 2) + (nneigh == 3);
            }
        }

        /* write a checkpoint if required */

        if (step%s == 0){
            void* ptr = (void*) grid;
            snprintf(cp_fname,21,"checkpoint_%05d.pgm",step);
            write_pgm_image(ptr,maxval,xsize,ysize,cp_fname);
        }
    } /* evolution step */
    free(cp_fname);
}

void static_evolution(char* grid, char* neigh, const int n, const int s, const int maxval, const int xsize, const int ysize){
    char* cp_fname = (char*)malloc(21*sizeof(char));
    
    /* static evolution loop */
    neigh = (char*)malloc(xsize*ysize*sizeof(char));
    for (int step = 1; step <=n ; step++){
        
        for (int i = 0; i<ysize;i++){
            int irow = i *xsize; 
            for (int j = 0;j<xsize;j++){
                int nneigh = 0;

                for(int ii = i+ysize-1; ii < i+ysize+2; ii++){
                    int iirow = (ii%ysize) * xsize;
                    for (int jj = j+xsize-1; jj < j+xsize+2; jj++){
                        int jjcol = jj%xsize;
                        nneigh += grid[iirow+jjcol] * (!(iirow == irow && jjcol == j));
                    }
                }

                neigh[irow + j] = nneigh;
            }
        }

        for (int i = 0; i<ysize;i++){
            int irow = i *xsize; 
            for (int j = 0;j<xsize;j++){
                grid[irow + j] = grid[irow+j] * (neigh[irow+j] == 2) + (neigh[irow + j] == 3);
            }
        }


        
        /* write a checkpoint if required */

        if (step%s == 0){
            snprintf(cp_fname,21,"checkpoint_%05d.pgm",step);
            write_pgm_image((void*)grid,maxval,xsize,ysize,cp_fname);
        }
    }/* evolution step */
    free(cp_fname);
    free(neigh);

}