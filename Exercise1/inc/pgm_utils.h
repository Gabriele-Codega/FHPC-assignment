#ifndef PGM_UTILS_H
#define PGM_UTILS_H

int write_pgm_image( void *, int , int , int , const char *);
void read_pgm_image( void **, int *, int *, int *, const char *);
void swap_image( void *, int , int , int  );

#endif