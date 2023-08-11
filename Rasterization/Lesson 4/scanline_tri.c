/*
 *
 * lines.c - simple GLUT app that draws one frame with a single point at origin
 *
 * To build:  gcc -framework OpenGL -framework GLUT lines.c -o lines
 *
 */
/*************************************************************************/
/* defines                                                               */
/*************************************************************************/
#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#define PI 3.1415926
#define SWAP(a, b) {int temp = (a); (a) = (b); (b) = temp;}
#define SWAPPOINT(a, b) {POINT temp = (a); (a) = (b); (b) = temp;}
#define ipart(x) ((int)(x))
#define fpart(x) ((x) - ipart(x))
#define X 0
#define Y 1
#define Z 2
#define W 3
#define R 0
#define G 1
#define B 2
#define A 3
#endif

typedef struct point{
    float pos[4];
    float color[4];
    float tex_coords[4]; // vector routines need 4 in professional for 3D. S T R Q
} POINT;

typedef struct image{
    int width;
    int height;
    unsigned char data[4096][4096][4]; // just for making code easy and flexible. char is 0-255
} IMAGE;

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


/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int half_window_size = 400;
int Mojave_WorkAround = 1;
int draw_one_frame = 1;
int test = 0;
int current[800][800];


/*************************************************************************/
/* Texture functions                                                     */
/*************************************************************************/
void texture_sample(IMAGE *tex, float s, float t, float tex_color[4]) { // passing address is good for memory management. 4 bytes vs 16 MB
    int u = s * tex -> width; // -> = (*tex). 
    int v = t * tex -> height;
    tex_color[R] = tex -> data[v][u][R] / 255; //Since R is a unigned char
}


/*************************************************************************/
/* vector functions                                                      */
/*************************************************************************/
void add_vector(float a[4], float b[4], float result[4]) {
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
    result[3] = a[3] + b[3];
}

void subtract_vector(float a[4], float b[4], float result[4]) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
    result[3] = a[3] - b[3];
}

void multiply_vector(float a[4], float b[4], float result[4]) {
    result[0] = a[0] * b[0];
    result[1] = a[1] * b[1];
    result[2] = a[2] * b[2];
    result[3] = a[3] * b[3];
}

void scalar_multiply_vector(float f, float a[4], float result[4]) {
    result[0] = a[0] * f;
    result[1] = a[1] * f;
    result[2] = a[2] * f;
    result[3] = a[3] * f;
}

void set_vector(float v[4], float r, float g, float b, float a) {
    v[0] = r;
    v[1] = g;
    v[2] = b;
    v[3] = a;
}

/*************************************************************************/
/* utility functions                                                     */
/*************************************************************************/

/*
 * random_float()
 */
float random_float( float low, float high ) {
    return( (float)(low + (rand()/(float)RAND_MAX)*(high - low)) );
}

/*
 * draw_point()
 */
void draw_point( float x, float y ) {
    /*
     * draw points
     */
    glBegin(GL_POINTS);
    glVertex2f( x , y ); // + or - to translate x and y by pixels
    glEnd();
}

/*
 * draw_point_struct()
 */
void draw_point_struct(POINT p) {
    /*
     * draw points
     */
    glColor4f(p.color[R], p.color[G], p.color[B], 1.0);
    glBegin(GL_POINTS);
    glVertex2f(p.pos[X], p.pos[Y]); // + or - to translate x and y by pixels
    glEnd();
}

/*
 * set_color()
 */
void set_color( float r, float g, float b, float a ) {
    glColor4f( r, g, b, a );
}

/*
 * draw_line_structs()
 */
void draw_line_color(POINT p0, POINT p1) {
    int steep = fabs(p1.pos[Y] - p0.pos[Y]) > fabs(p1.pos[X] - p0.pos[X]); // check if dy > dx

    POINT one, two;
    if (steep) { // IF STEEP, SWAP X AND Y
        one.pos[X] = p0.pos[Y];
        one.pos[Y] = p0.pos[X];
        two.pos[X] = p1.pos[Y];
        two.pos[Y] = p1.pos[X];
    } else { // ELSE SET THE VALUES TO DEFAULT VALUES
        one.pos[X] = p0.pos[X];
        one.pos[Y] = p0.pos[Y];
        two.pos[X] = p1.pos[X];
        two.pos[Y] = p1.pos[Y];
    }

    if (one.pos[X] > two.pos[X]) { // guarantees one is on the right of two
        SWAPPOINT(one, two);
    }
    
    float dx =  one.pos[X] - two.pos[X]; // CHANGE IN X
    float dy = fabs(one.pos[Y] - two.pos[Y]); // CHANGE IN Y (since we already confirmed that x1 > x0)
    float error = 2 * dy - dx; // BRESENHAM ERROR

    float color_delta[4];
    float color_inc[4];
    subtract_vector(p1.color, p0.color, color_delta); // CALCULATES CHANGE IN RGBA
    scalar_multiply_vector(1/dx, color_delta, color_inc); // DIVIDE CHANGE IN RGBA BY STEPS

    int ystep = (one.pos[Y] < two.pos[Y]) ? 1 : -1; // CHECKS IF LINE IS TO RIGHT OR LEFT
    float y = one.pos[Y];
    float x = one.pos[X];

    while(x <= two.pos[X]) {
        if(steep) {
            p0.pos[X] = y;
            p0.pos[Y] = x;
        } else {
            p0.pos[X] = x;
            p0.pos[Y] = y;
        }
        draw_point_struct(p0);
        add_vector(p0.color, color_inc, p0.color);
        error -= 2 * dy; // BRESENHAM INTEGRATION
        if(error < 0) {
            y += ystep;
            error += 2 * dx;
        }
        x++;
    }
}

/* 
 * draw_line_texture()
 */
void draw_line_texture(POINT p0, POINT p1) {
    int steep = fabs(p1.pos[Y] - p0.pos[Y]) > fabs(p1.pos[X] - p0.pos[X]);

    float x0, y0, x1, y1;
    if (steep) {
        x0 = p0.pos[Y];
        y0 = p0.pos[X];
        x1 = p1.pos[Y];
        y1 = p1.pos[X];
    } else {
        x0 = p0.pos[X];
        y0 = p0.pos[Y];
        x1 = p1.pos[X];
        y1 = p1.pos[Y];
    }

    if (x0 > x1) {
        SWAP(x0, x1);
        SWAP(y0, y1);
    }
    
    float dx = x1 - x0;
    float dy = fabs(y1 - y0);
    float error = 2 * dy - dx;

    float tex_delta[4];
    float tex_inc[4];
    subtract_vector(p1.tex_coords, p0.tex_coords, tex_delta); // gets the actual delta between texture coordinates
    scalar_multiply_vector(1/dx, tex_delta, tex_inc);

    int ystep = (y0 < y1) ? 1 : -1;
    float y = y0;
    float x = x0;

    while(x <= x1) {
        if(steep) {
            p0.pos[X] = y;
            p0.pos[Y] = x;
        } else {
            p0.pos[X] = x;
            p0.pos[Y] = y;
        }
        draw_point_struct(p0);
        add_vector(p0.tex_coords, tex_inc, p0.tex_coords);
        error -= 2 * dy;
        if(error < 0) {
            y += ystep;
            error += 2 * dx;
        }
        x++;
    }
}


/*
 * 
 */
int isCollinear(POINT p0, POINT p1, POINT p2) { // checks if the area is == 0, if yes, then these points cannot form a triangle
    float area = fabs(p0.pos[X]*(p1.pos[Y] - p2.pos[Y]) + p1.pos[X]*(p2.pos[Y] - p0.pos[Y]) + p2.pos[X]*(p0.pos[Y] - p1.pos[Y]));
    return area == 0;
}

/*
 * draw_flatbottom_triangle() 
 */
void draw_flatbottom_triangle(POINT p0, POINT p1, POINT p2) {
    float left_delta[4], left_inc[4], right_delta[4], right_inc[4];
    // if (p0.pos[Y] < p1.pos[Y]) {SWAPPOINT(p0, p1);} // confirms p0 is the highest point
    // if (p1.pos[X] < p2.pos[X]) {SWAPPOINT(p1, p2);} // 
    if (isCollinear(p0, p1, p2) || p1.pos[Y] != p2.pos[Y]) {
        return;
    }
    float inv_left = (p0.pos[X] - p1.pos[X])/(p0.pos[Y] - p1.pos[Y]); // provides slope from middle to left
    float inv_right = (p0.pos[X] - p2.pos[X])/(p0.pos[Y] - p2.pos[Y]); // provides slope from mmiddle to right

    POINT x0, x1;
    x0 = p0; // sets starting x to p0 x
    x1 = p0; // sets starting x to p0 x
    subtract_vector(p0.color, p1.color, left_delta);
    subtract_vector(p0.color, p2.color, right_delta);

    scalar_multiply_vector(1.0/(p0.pos[Y] - p1.pos[Y]), left_delta, left_inc);
    scalar_multiply_vector(1.0/(p0.pos[Y] - p1.pos[Y]), right_delta, right_inc);

    for (int y = p0.pos[Y]; y >= p1.pos[Y]; y--) { // for every iteration down from middle point
        x0.pos[Y] = y;
        x1.pos[Y] = y;
        x0.pos[X] += inv_left;
        x1.pos[X] += inv_right;
        draw_line_color(x0, x1);
        subtract_vector(x0.color, left_inc, x0.color);
        subtract_vector(x1.color, right_inc, x1.color);
    }
}

/*
 * draw_flattopped_trangle()
 */
void draw_flattopped_triangle(POINT p0, POINT p1, POINT p2) { //p0 is lowest, p1 is top right, p2 is top left
    float left_delta[4], left_inc[4], right_delta[4], right_inc[4];
    if (p0.pos[Y] > p1.pos[Y]) {SWAPPOINT(p0, p1);} // confirms p0 is the lowest
    if (p1.pos[X] < p2.pos[X]) {SWAPPOINT(p1, p2);} // confirms p1 is to the right of p2
    if (isCollinear(p0, p1, p2) || p1.pos[Y] != p2.pos[Y]) {
        return;
    }
    float inv_left = (p0.pos[X] - p1.pos[X])/(p0.pos[Y] - p1.pos[Y]); // provides slope from middle to left
    float inv_right = (p0.pos[X] - p2.pos[X])/(p0.pos[Y] - p2.pos[Y]); // provides slope from middle to right

    POINT x0, x1;
    x0 = p0; // sets starting x to p0 x
    x1 = p0; // sets starting x to p0 x
    subtract_vector(p0.color, p1.color, left_delta);
    subtract_vector(p0.color, p2.color, right_delta);

    scalar_multiply_vector(1.0/(p0.pos[Y] - p1.pos[Y]), left_delta, left_inc);
    scalar_multiply_vector(1.0/(p0.pos[Y] - p1.pos[Y]), right_delta, right_inc);

    for (int y = p0.pos[Y]; y <= p2.pos[Y]; y++) { // for every iteration down from middle point
        x0.pos[Y] = y;
        x1.pos[Y] = y;
        x0.pos[X] += inv_left;
        x1.pos[X] += inv_right;
        draw_line_color(x0, x1);
        add_vector(x0.color, left_inc, x0.color);
        add_vector(x1.color, right_inc, x1.color);
    }
}

/*
 * draw_triangle()
 */
void draw_triangle(POINT p0, POINT p1, POINT p2) {
    if (isCollinear(p0, p1, p2)) {
        return;
    }
    if (p0.pos[Y] < p1.pos[Y]) {SWAPPOINT(p0, p1);} // confirms p0 is >= p1
    if (p0.pos[Y] < p2.pos[Y]) {SWAPPOINT(p0, p2);} // confirms p0 is >= p2
    if (p1.pos[Y] < p2.pos[Y]) {SWAPPOINT(p1, p2);} // confirms p1 is >= p2 || p0 >= p1 >= p2
    if (p1.pos[Y] == p2.pos[Y] || p1.pos[Y] == p0.pos[Y]) {  //checks p1 == p2 or p1 == p0
        if (p0.pos[Y] > p1.pos[Y]) {
            draw_flatbottom_triangle(p0, p1, p2);
        }
        else if (p2.pos[Y] < p1.pos[Y]) {
            draw_flattopped_triangle(p0, p1, p2);
        }
    }
    else {
        POINT p3; //for diving triangle into a flat bottom and flat top
        float rgb_delta[4], rgb_inc[4];
        float inv_rgb = (p0.pos[X] - p1.pos[X])/(p0.pos[Y] - p1.pos[Y]); // provides slope from middle to middle_y point
        subtract_vector(p0.color, p1.color, rgb_delta);
        scalar_multiply_vector((p0.pos[Y] - p1.pos[Y])/(p0.pos[Y] - p2.pos[Y]), rgb_delta, p3.color);
        p3.pos[X] = p0.pos[X] + (p1.pos[Y] - p0.pos[Y]) / (p2.pos[Y] - p0.pos[Y]) * (p2.pos[X] - p0.pos[X]);
        p3.pos[Y] = p1.pos[Y];
        draw_flatbottom_triangle(p0, p1, p3);
        draw_flattopped_triangle(p1, p3, p2);
    }
}

/*
 * draw_random_triangles()
 */
void draw_random_triangles() {
    for (int i = 0; i < 1; i++) {
        POINT p0, p1, p2;
        set_vector(p0.pos, (int)random_float(-400, 400), (int)random_float(-400, 400), 0, 0);
        set_vector(p1.pos, (int)random_float(-400, 400), (int)random_float(-400, 400), 0, 0);
        set_vector(p2.pos, (int)random_float(-400, 400), (int)random_float(-400, 400), 0, 0);
        set_vector(p0.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p1.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        set_vector(p2.color, random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        draw_triangle(p0, p1, p2);
    }

}


/*************************************************************************/
/* GLUT functions                                                        */
/*************************************************************************/
/*
 * display routine
 */
void display( void )
{
    if( draw_one_frame == 0 )
        return;
    
    /*
     * clear color buffer
     */
    glClear( GL_COLOR_BUFFER_BIT );
	
    if (test == 0) {
        POINT p0, p1, p2;
        set_vector(p0.pos, 0,    300, 0, 0);
        set_vector(p1.pos, -500, 0, 0, 0);
        set_vector(p2.pos, 300, 0, 0, 0);
        set_vector(p0.color, 1, 0, 0, 1);
        set_vector(p1.color, 0, 0, 1, 1);
        set_vector(p2.color, 0, 1, 0, 1);
        draw_flatbottom_triangle(p0, p1, p2);
    }
    else if (test == 1) {
        POINT p0, p1, p2;
        set_vector(p0.pos, 0, -300, 0, 0);
        set_vector(p1.pos, -400, 0, 0, 0);
        set_vector(p2.pos, 300, 0, 0, 0);
        set_vector(p0.color, 1, 0, 0, 1);
        set_vector(p1.color, 0, 0, 1, 1);
        set_vector(p2.color, 0, 1, 0, 1);
        draw_flattopped_triangle(p0, p2, p1);
    }
    else if (test == 2) {
        draw_random_triangles();
    }
    else if (test == 3) {
        
    }
    else if (test == 4) {
        
    }
    else if (test == 5) {
       
    }
    else if (test == 6) {
        
    }

    /*
     * show results
     */
    glutSwapBuffers();
    
    draw_one_frame = 0;
}

/*
 * Key routine
 */
static void Key( unsigned char key, int x, int y )
{
    switch( key )
    {     
        case 'a':       		                         draw_one_frame = 1;     glutPostRedisplay();    break;
		case 't':       test++; if( test > 6 ) test = 0; draw_one_frame = 1;     glutPostRedisplay();    break;
        case 'y':       test--; if (test < 0 ) test = 6; draw_one_frame = 1;     glutPostRedisplay();    break;
		case 'q':       exit( 0 );                             											 break;
        case '\033':    exit( 0 );																		 break;
    }
}

/*
 * main function
 */
int main( int argc, char **argv )
{
    glutInit( &argc, argv );
    srand( time(NULL) );

    /*
     * create window
     */
    glutInitWindowSize( half_window_size, half_window_size );
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
    glutCreateWindow( "My OpenGL Window" );

    /*
     * setup display(), Key() funtions
     */
    glutDisplayFunc( display );
    glutKeyboardFunc( Key );

    /*
     * setup OpenGL state
     */
    glClearColor( 0.0, 0.0, 0.0, 1.0 );
    gluOrtho2D( -half_window_size, half_window_size ,-half_window_size, half_window_size );
    // glPointSize( 1.0 );
    glColor4f( 1.0, 0.0, 0.0, 1.0 );

	if( Mojave_WorkAround )
	{
		// Necessary for Mojave. Has to be different dimensions than in glutInitWindowSize();
		glutReshapeWindow( half_window_size * 2, half_window_size * 2 );
	}
	
    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}