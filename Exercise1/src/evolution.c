#include <evolution.h>
#include <pgm_utils.h>

void ordered_evolution(char* grid, char* neigh, const int n, const int s){
    char* cp_fname = (char*)malloc(20*sizeof(char));
    
    /* ordered evolution loop */
    for (int step = 1; step <=n ; step++){
        
        for (int i = 0; i<k;i++){
            int irow = i *k; 
            for (int j = 0;j<k;j++){
                int nneigh = 0;

                for(int ii = i+k-1; ii < i+k+2; ii++){
                    int iirow = (ii%k) * k;
                    for (int jj = j+k-1; jj < j+k+2; jj++){
                        int jjcol = jj%k;
                        nneigh += grid[iirow+jjcol] * (!(iirow == irow && jjcol == j));
                    }
                }

                grid[irow + j] = grid[irow+j] * (nneigh == 2) + (nneigh == 3);
            }
        }

        
        /* write a checkpoint if required */

        if (step%s == 0){
            sprintf(cp_fname,"checkpoint_%05d.pgm",step);
            write_pgm_image((void*)grid,maxval,xsize,ysize,cp_fname);
        }
    } /* evolution step */
}

void static_evolution(char* grid, char* neigh, const int n, const int s){
    char* cp_fname = (char*)malloc(20*sizeof(char));
    
    /* static evolution loop */
    neigh = (char*)malloc(k*k*sizeof(char));
    for (int step = 1; step <=n ; step++){
        
        for (int i = 0; i<k;i++){
            int irow = i *k; 
            for (int j = 0;j<k;j++){
                int nneigh = 0;

                for(int ii = i+k-1; ii < i+k+2; ii++){
                    int iirow = (ii%k) * k;
                    for (int jj = j+k-1; jj < j+k+2; jj++){
                        int jjcol = jj%k;
                        nneigh += grid[iirow+jjcol] * (!(iirow == irow && jjcol == j));
                    }
                }

                neigh[irow + j] = nneigh;
            }
        }

        for (int i = 0; i<k;i++){
            int irow = i *k; 
            for (int j = 0;j<k;j++){
                grid[irow + j] = grid[irow+j] * (neigh[irow+j] == 2) + (neigh[irow + j] == 3);
            }
        }


        
        /* write a checkpoint if required */

        if (step%s == 0){
            sprintf(cp_fname,"checkpoint_%05d.pgm",step);
            write_pgm_image((void*)grid,maxval,xsize,ysize,cp_fname);
        }
    }/* evolution step */
}