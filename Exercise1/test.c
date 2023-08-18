#include <stdio.h>

int main(){
    int step = 3;
    char fname[20];
    sprintf(fname,"checkpoint_%05d.pgm",step);
    puts(fname);
    return 0;
}