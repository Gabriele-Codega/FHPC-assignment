#include <stdio.h>
#include <string.h>

int main(){
    char* string;
    string = strdup("ba,cs");
    char* token = strsep(&string,",");
    printf("%s,%s",token,string);
    return 0;
}