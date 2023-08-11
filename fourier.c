/*************************************************************************/
/* defines                                                               */
/*************************************************************************/

#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>

#define PI 3.1415926

int main(void)
{
    float array[360];
    for (int x = 0; x < 360; x++) 
    {
        float theta = x/359.0 * 2 * PI;
        array[x] = sin(theta) + sin(2 * theta) + sin( 5 * theta);
    }

    for (int k = 1; k < 360; k++) // discrete sin() transform
    {
        float sum = 0;
        for (int x = 0; x < 360; x++) 
        {
            float theta = x/359.0 * 2 * PI;
            float test = sin(k * theta);
            sum += test * array[x];
        }
        printf("average = %f, k = %d\n", sum / 360.0, k);
    }
    return 0;
}

