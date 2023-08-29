#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include <mpi.h>

#include <io_utils.h>


void read_header(const char* fname, int* params){
    int xsize,ysize,maxval;
    FILE* image_file;
    image_file = fopen(fname,"r");
    if (!image_file){
        printf("Unable to find file '%s'. Aborting.",fname);
    }
    char   *line = NULL;
    char    MagicN[2];

    size_t  k, n;
    
    // get the Magic Number
    k = fscanf(image_file, "%2s%*c", MagicN );

    // skip all the comments
    k = getline( &line, &n, image_file);
    while ( (k > 0) && (line[0]=='#') )
        k = getline( &line, &n, image_file);

    if (k > 0)
        {
        k = sscanf(line, "%d%*c%d%*c%d%*c", &xsize, &ysize, &maxval);
        if ( k < 3 )
        fscanf(image_file, "%d%*c", &maxval);
        }
    else
        {
        maxval = -1;    // this is the signal that there was an I/O error
                        // while reading the image header
        free( line );
        }
    free( line );
    params[0] = xsize;
    params[1] = ysize;
    params[2] = maxval;
    params[3] = ftell(image_file);
    fclose(image_file);
}


void read_data(MPI_File fh, char* img,const int thwork, const MPI_Offset thoffset, const MPI_Offset procoffset, const int globaloffset){
    
    MPI_Status status;
    #pragma omp critical
        MPI_File_read_at_all(fh,globaloffset+procoffset+thoffset, (void*)(img+thoffset),thwork,MPI_BYTE,&status);
    
    #pragma omp barrier
}

extern void 
write_header(const char* fname, const int xsize, const int ysize, const int maxval, int* offset);

extern void 
write_checkpoint(char* fname, const int step, const char* grid,
                const int procrank, const int procoffset,
                const int thoffset, const int thwork,
                const int xsize, const int ysize, const int maxval);
