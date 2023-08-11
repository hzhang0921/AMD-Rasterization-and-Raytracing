#include <stdio.h>
#define DEBUG

int main(int argc, char** argv){
    #ifdef DEBUG
        printf("hello\n");
    #else 
        printf("goodbye!\n");
    #endif
}