#ifndef EVOLUTION_H
#define EVOLUTION_H

void 
ordered_evolution(char* full_grid, char* neigh, const int n, const int s, 
                const int maxval, const int xsize, const int ysize,
                const int procwork, const int procoffset);


void 
static_evolution(char* full_grid, char* neigh, const int n, const int s, 
                const int maxval, const int xsize, const int ysize,
                const int procwork, const int procoffset);


#endif