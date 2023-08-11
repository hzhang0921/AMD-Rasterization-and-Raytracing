#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#define PI 3.14
#endif

/*************************************************************************/
/* header files                                                          */
/*************************************************************************/
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>

#define IMAGE_WIDTH 800
#define IMAGE_HEIGHT 800
#define HALF_WINDOW_SIZE 400
#define R 0
#define G 1
#define B 2

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int Mojave_WorkAround = 1; // used for a bug in Mojave
int draw_one_frame = 1;
int test = 0;
int textureIndex = 0;
int current[IMAGE_HEIGHT][IMAGE_WIDTH];
float depth[IMAGE_HEIGHT][IMAGE_WIDTH];
float colorbuffer[IMAGE_HEIGHT][IMAGE_WIDTH][3];

/*************************************************************************/
/* buffer functions                                                      */
/*************************************************************************/
/*
 * clear_colorbuffer() resets all pixels to a default value of black (0, 0, 0, 0)
 */
void clear_colorbuffer() {
    for (int y = IMAGE_HEIGHT - 1; y >= 0; --y) {
        for (int x = 0; x < IMAGE_WIDTH; ++x) {
            colorbuffer[x][y][R] = 0;
            colorbuffer[x][y][G] = 0;
            colorbuffer[x][y][B] = 0;
        }
    }
}

/*
 * clear_depthbuffer() resets all pixel depths to a default value of 1000 (1000 pixels away into the screen)
 */

void clear_depthbuffer() {
    for (int y = IMAGE_HEIGHT - 1; y >= 0; --y) {
        for (int x = 0; x < IMAGE_WIDTH; ++x) {
            depth[x][y] = 1000;
        }
    }
}

/*
 * draw_buffer()
 */
void draw_colorbuffer() {
    for (int y = IMAGE_HEIGHT - 1; y >= 0; --y) {
        for (int x = 0; x < IMAGE_WIDTH; ++x) {
            glBegin(GL_POINTS);
            glColor4f(colorbuffer[x][y][R], colorbuffer[x][y][G], colorbuffer[x][y][B], 1);
            glVertex2f(x - HALF_WINDOW_SIZE, y - HALF_WINDOW_SIZE); // subtract 400 since color buffer is from 0 - 800
            glEnd();
        }
    }
}

/*
 *
 */
void firstImage() {
    for (int y = IMAGE_HEIGHT - 1; y >= 0; --y) {
        for (int x = 0; x < IMAGE_WIDTH; ++x) {
            float r = float(x) / (IMAGE_WIDTH - 1);
            float g = float(y) / (IMAGE_HEIGHT -1 );
            float b = 0.25;

            colorbuffer[y][x][R] = r;
            colorbuffer[y][x][G] = g;
            colorbuffer[y][x][B] = b;
        }
    }
}


/*************************************************************************/
/* GLUT functions                                                        */
/*************************************************************************/
/*
 * display routine
 */
void display(void)
{
    if (draw_one_frame == 0)
        return;

    /*
     * clear color buffer
     */
    glClear(GL_COLOR_BUFFER_BIT);
    clear_colorbuffer();
    clear_depthbuffer();

    if (test == 0) 
    {

    }
    else if (test == 1) 
    {


    }
    else if (test == 2) 
    {

    }
    else if (test == 3) 
    {
 
    }
    else if (test == 4)
    {

    }
    else if (test == 5)
    {

    }
    else if (test == 6)
    {

    }
    else if (test == 7)
    {

    }
    else if (test == 8)
    {

    }
    else if (test == 9)
    {

    }
    else if (test == 10)
    {

    }
    else if (test == 11)
    {

    }
    else if (test == 12)
    {

    }
    else if (test == 13)
    {

    }
    else if (test == 14)
    {

    }
    else if (test == 15)
    {

    } 
    else if (test == 16)
    {

    }


    /*
     * show results
     */
    draw_colorbuffer();
    glutSwapBuffers();

    draw_one_frame = 0;
}

/*
 * Key routine
 */
static void Key(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'a':
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 't':
        test++;
        if (test > 16)
            test = 0;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 'T':
        test--;
        if (test < 0)
            test = 16;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 'q':
        exit(0);
        break;
    case '\033':
        exit(0);
        break;
    }
}

/*
 * main function
 */
int main(int argc, char **argv)
{
    char filename[128];
    for (int i = 0; i < argc; i++) 
    {
        printf("arg %d: %s \n", i, argv[i]);
        if (strcmp(argv[i], "-scene") == 0) { // '' is for char, "" is for string
            if (i+1 >= argc) {
                printf("ERROR: no scene file specified\n");
                exit(1);
            } 
            strcpy(filename, argv[++i]);
            printf("scene file is %s\n", filename);
        }
    }
    glutInit(&argc, argv);
    srand(time(NULL));

    /*
     * create window
     */
    glutInitWindowSize(HALF_WINDOW_SIZE, HALF_WINDOW_SIZE);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutCreateWindow("My OpenGL Window");

    /*
     * setup display(), Key() funtions
     */
    glutDisplayFunc(display);
    glutKeyboardFunc(Key);

    /*
     * setup OpenGL state
     */
    glClearColor(0.0, 0.0, 0.0, 1.0);
    gluOrtho2D(-HALF_WINDOW_SIZE, HALF_WINDOW_SIZE, -HALF_WINDOW_SIZE, HALF_WINDOW_SIZE); // Original for -400 -> 400
    // glPointSize( 1.0 );
    glColor4f(1.0, 0.0, 0.0, 1.0);

    if (Mojave_WorkAround)
    {
        // Necessary for Mojave. Has to be different dimensions than in glutInitWindowSize();
        glutReshapeWindow(HALF_WINDOW_SIZE * 2, HALF_WINDOW_SIZE);
    }

    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}