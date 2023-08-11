#include <stdio.h>
#define PI 3.14159

int mainA(void) {
    int x= 0;

    printf("%s\n", "Haoyang Zhang");
    printf("%f\n", PI);
    return 0;
}

int mainB(void) {
    int x = 15;
    if (x > 100) {
        printf("x is humongous! \n");
    }
    else if (x > 10) {
        printf("x is big \n");
    }
    else {
        printf("x is puny \n");
    }
    return 0;
}

int mainC(void) {
    int i;
    for (i = 0; i < 1000000; i++) {
        printf("%d, %d\n", i, i*i);
    }
    return 0;
}

float negate(float value) {
    return (-1.0 * value);
}

float square(float value) {
    return (value * value);
}

int main(void) {
    mainA();
    mainB();
    // mainC();
    float y;
    y = negate(PI);
    printf("negated y = %f \n", y);
    printf("y squared = %f ", y * y);
}