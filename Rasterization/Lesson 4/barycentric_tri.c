#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
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
#include <strings.h> // library for strings

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/
#define PI 3.1415926
#define SWAP(a, b)      \
    {                   \
        int temp = (a); \
        (a) = (b);      \
        (b) = temp;     \
    }
#define SWAPPOINT(a, b)   \
    {                     \
        POINT temp = (a); \
        (a) = (b);        \
        (b) = temp;       \
    }
#define ipart(x) ((int)(x))
#define fpart(x) ((x)-ipart(x))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CLAMP(a, low, high) (((a) < (low)) ? (low) : (((a) > (high)) ? (high) : (a)))
#define X 0
#define Y 1
#define Z 2
#define W 3
#define R 0
#define G 1
#define B 2
#define A 3

typedef struct point
{
    float pos[4];
    float color[4];
    float tex_coords[4]; // vector routines need 4 in professional for 3D.
    float world[4]; // imaginary 3D world coordinates with x y and z
} POINT;

typedef struct image
{
    int width;
    int height;
    unsigned char data[4096][4096][4]; // unsigned char is 0-255
} IMAGE;

typedef struct triangle
{
    int index_list[3]; // list of indices into the vertex list
} TRIANGLE;

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int half_window_size = 400; // used for standarizing coordinate plane
int Mojave_WorkAround = 1; // used for a bug in Mojave
int draw_one_frame = 1;
int test = 0;
int textureIndex = 0;
int current[800][800];
float depth[800][800];
float colorbuffer[800][800][4];
// IMAGE Textures[20]; // allocates memory for 20 images of 4096 x 4096 RGBA pixels
POINT vertex_list[1000000]; // vertices of the triangles
TRIANGLE triangle_list[1000000]; // list of triangles

float y_angle = 0;
int num_triangles = 0;
int num_verts = 0;

/*************************************************************************/
/* vector functions                                                      */
/*************************************************************************/
void add_vector(float a[4], float b[4], float result[4])
{
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
    result[3] = a[3] + b[3];
}

void subtract_vector(float a[4], float b[4], float result[4])
{
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
    result[3] = a[3] - b[3];
}

void multiply_vector(float a[4], float b[4], float result[4])
{
    result[0] = a[0] * b[0];
    result[1] = a[1] * b[1];
    result[2] = a[2] * b[2];
    result[3] = a[3] * b[3];
}

void scalar_multiply_vector(float f, float a[4], float result[4])
{
    result[0] = a[0] * f;
    result[1] = a[1] * f;
    result[2] = a[2] * f;
    result[3] = a[3] * f;
}

void set_vector(float v[4], float r, float g, float b, float a)
{
    v[0] = r;
    v[1] = g;
    v[2] = b;
    v[3] = a;
}

void copy_vector(float input[4], float output[4])
{
    output[0] = input[0];
    output[1] = input[1];
    output[2] = input[2];
    output[3] = input[3];
}

void UCset_vector(unsigned char v[4], unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    v[0] = r;
    v[1] = g;
    v[2] = b;
    v[3] = a;
}

/*************************************************************************/
/* buffer functions                                                      */
/*************************************************************************/
/*
 * clear_colorbuffer() resets all pixels to a default value of black (0, 0, 0, 0)
 */
void clear_colorbuffer() {
    for (int x = 0; x < 800; x++) {
        for (int y = 0; y < 800; y++) {
            colorbuffer[x][y][R] = 0;
            colorbuffer[x][y][G] = 0;
            colorbuffer[x][y][B] = 0;
            colorbuffer[x][y][A] = 0;
        }
    }
}

/*
 * clear_depthbuffer() resets all pixel depths to a default value of 1000 (1000 pixels away into the screen)
 */

void clear_depthbuffer() {
    for (int x = 0; x < 800; x++) {
        for (int y = 0; y < 800; y++) {
            depth[x][y] = 1000;
        }
    }
}

/*
 * draw_buffer()
 */
void draw_colorbuffer() {
    for (int x = 0; x < 800; x++) {
        for (int y = 0; y < 800; y++) {
            glColor4f(colorbuffer[x][y][R], colorbuffer[x][y][G], colorbuffer[x][y][B], colorbuffer[x][y][A]);
            glBegin(GL_POINTS);
            glVertex2f(x - 400, y - 400); // subtract 400 since color buffer is from 0 - 800
            glEnd();
        }
    }
}

/*************************************************************************/
/* utility functions                                                     */
/*************************************************************************/

/*
 * random_float()
 */
float random_float(float low, float high)
{
    return ((float)(low + (rand() / (float)RAND_MAX) * (high - low)));
}

/*
 * draw_RGBApoint_struct()
 */
void draw_RGBApoint_struct(POINT p) {
    p.pos[X] += 400;
    p.pos[Y] += 400;
    float mix_factor = 0.5;
    if (p.pos[X] >= 0 && p.pos[X] < 800 && p.pos[Y] >= 0 && p.pos[Y] < 800) { // checks for out of bounds values
        if (p.pos[Z] < depth[(int)p.pos[X]][(int)p.pos[Y]]) { // checks if the point is closer than the previous point
            depth[(int)p.pos[X]][(int)p.pos[Y]] = p.pos[Z];
            set_vector(colorbuffer[(int)p.pos[X]][(int)p.pos[Y]], p.color[R], p.color[G], p.color[B], p.color[A]); 
        }
        else {
            ;
        }
    }
}

/*
 * draw_line_color()
 */
void draw_line_color( POINT p0, POINT p1 ) {
    float     pos_delta[4];
    float     color_delta[4];
    float     pos_inc[4];
    float     color_inc[4];
    POINT    p;
    
    if( (int)(p1.pos[X] - p0.pos[X]) == 0 && (int)(p1.pos[Y] - p0.pos[Y])== 0 )
    {
        /*
         * degenerate case - just draw point
         */
        draw_RGBApoint_struct( p0 );
    }
    else if( (int)p1.pos[Y] == (int)p0.pos[Y] )
    {
        /*
         * horizontal line
         */
        if( p0.pos[X] > p1.pos[X] )
        {
            SWAPPOINT(p0,p1);
        }

        p = p0;
        
        subtract_vector( p1.pos,   p0.pos,   pos_delta   );
        subtract_vector( p1.color, p0.color, color_delta );

        float inv_slope = 1.0/pos_delta[X];

        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );

        for( int x = p0.pos[X]; x <= p1.pos[X]; x++ )
        {
            draw_RGBApoint_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
        }
    }
    else if( (int)p1.pos[X] == (int)p0.pos[X] )
    {
        /*
         * vertical line
         */
        if( p0.pos[Y] > p1.pos[Y] )
        {
            SWAPPOINT(p0,p1);
        }

        p = p0;
        
        subtract_vector( p1.pos,   p0.pos,   pos_delta   );
        subtract_vector( p1.color, p0.color, color_delta );

        float inv_slope = 1.0/pos_delta[Y];
        
        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );

        for( int y = p0.pos[Y]; y <= p1.pos[Y]; y++ )
        {
            draw_RGBApoint_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
        }
    }
    else if( fabs(p1.pos[X] - p0.pos[X]) > fabs(p1.pos[Y] - p0.pos[Y]) )
    {
        /*
         * x-major line
         */
        if( p0.pos[X] > p1.pos[X] )
        {
            SWAPPOINT(p0,p1);
        }

        p = p0;

        subtract_vector( p1.pos,   p0.pos,   pos_delta   );
        subtract_vector( p1.color, p0.color, color_delta );

        float inv_slope = 1.0/pos_delta[X];

        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );

        for( int x = p0.pos[X]; x <= p1.pos[X]; x++ )
        {
            draw_RGBApoint_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
        }
    }
    else
    {
        /*
         * y-major line
         */
        if( p0.pos[Y] > p1.pos[Y] )
        {
            SWAPPOINT(p0,p1);
        }

        p = p0;

        subtract_vector( p1.pos,   p0.pos,   pos_delta   );
        subtract_vector( p1.color, p0.color, color_delta );

        float inv_slope = 1.0/pos_delta[Y];

        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );

        for( int y = p0.pos[Y]; y <= p1.pos[Y]; y++ )
        {
            draw_RGBApoint_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
        }
    }
}

/*
 * triangle_area
 */
float triangle_area(POINT p0, POINT p1, POINT p2) {
     return((p0.pos[X]*(p1.pos[Y] - p2.pos[Y]) + p1.pos[X]*(p2.pos[Y] - p0.pos[Y]) + p2.pos[X]*(p0.pos[Y] - p1.pos[Y]))/2.0);
}

/*
 * barycentric coordinates
 */
void calculate_barycentric(POINT p0, POINT p1, POINT p2, POINT p, float *alpha, float *beta, float *gamma) {
    float area = triangle_area(p0, p1, p2);
    if (area == 0) {
        *alpha = 0;
        *beta = 0;
        *gamma = 0;
        return;
    }
    *alpha = triangle_area(p, p1, p2) / area;
    *beta = triangle_area(p0, p, p2) / area;
    *gamma = triangle_area(p0, p1, p) / area;
}

/*
 * barycentric_triangle
 */
void barycentric_triangle(POINT p0, POINT p1, POINT p2) {
    float alpha, beta, gamma; // weights of p0, p1, and p2
    POINT p;
    int minX = MIN(MIN(p0.pos[X], p1.pos[X]), p2.pos[X]);
    int maxX = MAX(MAX(p0.pos[X], p1.pos[X]), p2.pos[X]);
    int minY = MIN(MIN(p0.pos[Y], p1.pos[Y]), p2.pos[Y]);
    int maxY = MAX(MAX(p0.pos[Y], p1.pos[Y]), p2.pos[Y]);
    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            set_vector(p.pos, x, y, 0, 0);
            calculate_barycentric(p0, p1, p2, p, &alpha, &beta, &gamma);
            if (alpha >= 0 && beta >= 0 && gamma >= 0) { // checks as long as alpha, beta, and gamma are >= 0, then you know that the point is within the Triangle
                set_vector(p.color, alpha * p0.color[R] + beta * p1.color[R] + gamma * p2.color[R], // sets the color at that point to the weighted average of the colors of the vertices
                                    alpha * p0.color[G] + beta * p1.color[G] + gamma * p2.color[G],
                                    alpha * p0.color[B] + beta * p1.color[B] + gamma * p2.color[B],
                                    1);
                draw_RGBApoint_struct(p);
            }
        }
    }
}  

/*
 * overlapping_triangleTest()
 */
void overlapping_triangleTest() {
    POINT p0, p1, p2;
    POINT p3, p4, p5;
    for (int i = 0; i < 5; i++) {
        set_vector(p0.pos, random_float(-400, 400), random_float(-400, 400), random_float(0, 100), 1);
        set_vector(p0.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p1.pos, random_float(-400, 400), random_float(-400, 400), random_float(0, 100), 1);
        set_vector(p1.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p2.pos, random_float(-400, 400), random_float(-400, 400), random_float(0, 100), 1);
        set_vector(p2.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p3.pos, random_float(-400, 400), random_float(-400, 400), random_float(0, 100), 1);
        set_vector(p3.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p4.pos, random_float(-400, 400), random_float(-400, 400), random_float(0, 100), 1);
        set_vector(p4.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p5.pos, random_float(-400, 400), random_float(-400, 400), random_float(0, 100), 1);
        set_vector(p5.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        barycentric_triangle(p0, p1, p2);
        barycentric_triangle(p3, p4, p5);
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
        POINT p0, p1;
        set_vector(p0.pos, -100, -100, 0, 1);
        set_vector(p0.color, 1, 0, 0, 1);
        set_vector(p1.pos, 100, 100, 0, 1);
        set_vector(p1.color, 0, 1, 0, 1);
        draw_line_color(p0, p1);
    }
    else if (test == 1)
    {
        POINT p0, p1, p2;
        set_vector(p0.pos, -400, -100, 0, 1);
        set_vector(p0.color, 1, 0, 0, 1);
        set_vector(p1.pos, 300, 100, 0, 1);
        set_vector(p1.color, 0, 1, 0, 1);
        set_vector(p2.pos, 0, 300, 0, 1);
        set_vector(p2.color, 0, 0, 1, 1);
        barycentric_triangle(p0, p1, p2);

    }
    else if (test == 2)
    {
        overlapping_triangleTest();
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
    case 'y':
        y_angle++;
        if (y_angle > 360) y_angle = 0;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 'Y':
        y_angle--;
        if (y_angle < 0) y_angle = 360;
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
    glutInit(&argc, argv);
    srand(time(NULL));

    /*
     * create window
     */
    glutInitWindowSize(half_window_size, half_window_size);
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
    gluOrtho2D(-half_window_size, half_window_size, -half_window_size, half_window_size); // Original for -400 -> 400
    // glPointSize( 1.0 );
    glColor4f(1.0, 0.0, 0.0, 1.0);

    if (Mojave_WorkAround)
    {
        // Necessary for Mojave. Has to be different dimensions than in glutInitWindowSize();
        glutReshapeWindow(half_window_size * 2, half_window_size * 2);
    }

    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}