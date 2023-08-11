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
} POINT;


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
 * draw_horizontal()
 */
void draw_horizontal(int x0, int y0, int x1, int y1) {
    if (y0 != y1) {
        return;
    }
    if (x0 > x1) {
        SWAP (x0, x1);
    }
    for (int x = x0; x < x1; x++) {
        draw_point(x, y0);
    }
}

/*
 * draw_vertical()
 */
void draw_vertical(int x0, int y0, int x1, int y1) {
    if (x0 != x1) {
        return;
    }
    if (y0 > y1) {
        SWAP (y0, y1);
    }
    for (int y = y0; y < y1; y++) {
        draw_point(x0, y);
    }
}

/*
 * draw_diagonal()
 */
void draw_diagonal(int x0, int y0, int x1, int y1) {
    float dx = x1 - x0;
    float dy = y1 - y0;
    float slope = dy/dx;
    float inv_slope = dx/dy;
    float y = y0;
    float x = x0;
    if (fabs(dy) > fabs(dx)) {
        if (y0 > y1) {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        x = x0;
        for (int y = y0; y < y1; y++) {
            draw_point(x, y);
            x += inv_slope;
        }
    }
    else {
        if (x0 > x1) {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        y = y0;
        for (int x = x0; x < x1; x++) {
            draw_point(x, y);
            y += slope;
        }    
    }

}

/*
* draw_line()
*/
void draw_line(int x0, int y0, int x1, int y1) {
    if (y0 == y1) {
        draw_horizontal(x0, y0, x1, y1);
    }
    else if(x0 == x1) {
        draw_vertical(x0, y0, x1, y1);
    }
    else {
        draw_diagonal(x0, y0, x1, y1);
    }
}

/*
 * draw_random
 */
void draw_random() {
    for (int i = 0; i < 100; i++) {
        int x0 = random_float(-400, 400);
        int y0 = random_float(-400, 400);
        int x1 = random_float(-400, 400);
        int y1 = random_float(-400, 400);
        float r = random_float(0, 1);
        float g = random_float(0, 1);
        float b = random_float(0, 1);
        set_color(r, g, b, 1);
        draw_line(x0, y0, x1, y1);
    }
}

/*
 * draw_fan()
 */
void draw_fan() {
    for (int angle = 0; angle < 360; angle += 10) {
        float y1 = sin(angle / 360.0 * 2 * PI);
        float x1 = cos(angle / 360.0 * 2 * PI);
        int radius = 200;
        draw_line(0, 0, radius * x1, radius * y1);
    }
}

/* 
 * draw_coordinate()
 */
void draw_coordinate() {
    int width = 10;
    for (int x = -400; x < 400; x += width) {
        draw_line(x, -400, x, 400);
    }
    for (int y = -400; y < 400; y += width) {
        draw_line(-400, y, 400, y);
    }
}

/* 
 * draw_Bresenham()
 */
void draw_bresenham(int x0, int y0, int x1, int y1) {
    int y = y0;
    int dx = x1 - x0;
    int dy = y1 - y0;
    int error = 2 * dy - dx;
    for (int x = x0; x < x1; x++) {
        draw_point(x, y);
        error += 2 * dy;
        if (error > 0) {
            y++;
            error -= 2 * dx;
        }
    }
}

/*
 * random_Bresenham()
 */
void random_bresenham() {
        for (int i = 0; i < 100; i++) {
        int x0 = random_float(-400, 400);
        int y0 = random_float(-400, 400);
        int x1 = random_float(-400, 400);
        int y1 = random_float(-400, 400);
        float r = random_float(0, 1);
        float g = random_float(0, 1);
        float b = random_float(0, 1);
        set_color(r, g, b, 1);
        draw_bresenham(x0, y0, x1, y1);
    }
}

/*
 * draw_wu()
 */
void draw_wu(int x0, int y0, int x1, int y1) {
    float dx = x1 - x0;
    float dy = y1 - y0;
    float slope = dy/dx;
    float inv_slope = dx/dy;
    int neg_slope = 1;
    if (slope < 0) {
        neg_slope = -1;
    }
    if (fabs(dy) > fabs(dx)) { //tests for steepness
        if (y0 > y1) {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        float x = x0;
        for (int y = y0; y < y1; y++) {
            set_color(fpart(x), fpart(x), fpart(x), 1);
            draw_point(ipart(x), y);
            set_color(1-fpart(x), 1-fpart(x), 1-fpart(x), 1);
            draw_point(ipart(x)+neg_slope, y);
            x += inv_slope;
        }
    }
    else {
        if (x0 > x1) {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        float y = y0;
        for (int x = x0; x < x1; x++) {
        set_color(fpart(y), fpart(y), fpart(y), 1);
        draw_point(x, ipart(y));
        set_color(1-fpart(y), 1-fpart(y), 1-fpart(y), 1);
        draw_point(x, ipart(y)+neg_slope);
        y += slope;
        }  
    }
}

/*
 * draw_antialiasing ()
 */
void draw_antialising(int x0, int y0, int x1, int y1) {
    if (y0 == y1) {
        draw_horizontal(x0, y0, x1, y1);
    }
    else if(x0 == x1) {
        draw_vertical(x0, y0, x1, y1);
    }
    else {
        draw_wu(x0, y0, x1, y1);
    }
}


/*
 * random_wu()
 */
void random_wu() {
    for (int i = 0; i < 100; i++) {
        int x0 = random_float(-400, 400);
        int y0 = random_float(-400, 400);
        int x1 = random_float(-400, 400);
        int y1 = random_float(-400, 400);
        float r = random_float(0, 1);
        float g = random_float(0, 1);
        float b = random_float(0, 1);
        set_color(r, g, b, 1);
        draw_antialising(x0, y0, x1, y1);
    }
}

/*
 * draw_line_structs()
 */
void draw_line_color(POINT p0, POINT p1) {
    float dx = p1.pos[X] - p0.pos[X];
    float dy = p1.color[Y] - p0.pos[Y];
    float dr = p1.color[R] - p0.color[R];
    float dg = p1.color[G] - p0.color[G];
    float db = p1.color[B] - p0.color[B];
    int steep = 0;
    if (fabs(dy) > fabs(dx)) {
        steep = 1;
    }
    float slope = dx/dy;
    float slope_in_red = dr/dy;
    float slope_in_green = dg/dy;
    float slope_in_blue = db/dy;
    if (steep == 1) {
        if (p0.pos[Y] > p1.pos[Y]) {
            SWAP(p0.pos[X], p1.pos[X]);
            SWAP(p0.pos[Y], p1.pos[Y]);
        }
        float x = p0.pos[X];
        float r = p0.color[R];
        float g = p0.color[G];
        float b = p0.color[B];
        for (int y = p0.pos[Y]; y < p1.pos[Y]; y++) {
            set_color(r, g, b, 1);
            draw_point(x, y);
            r += slope_in_red;
            g += slope_in_green;
            b += slope_in_blue;
            x += slope;
        }
    }
    else {
        float slope = dy/dx;
        float inv_slope = dx/dy;
        float slope_in_red = dr/dx;
        float slope_in_green = dg/dx;
        float slope_in_blue = db/dx;
        if (p0.pos[X] > p1.pos[X]) {
            SWAP(p0.pos[X], p1.pos[X]);
            SWAP(p0.pos[Y], p1.pos[Y]);
        }
        float y = p0.pos[Y];
        float r = p0.color[R];
        float g = p0.color[G];
        float b = p0.color[B];
        for (int x = p0.pos[X]; x < p1.pos[X]; x++) {
           set_color(r, g, b, 1);
            draw_point(x, y);
            r += slope_in_red;
            g += slope_in_green;
            b += slope_in_blue;
            y += slope;
            }    
        }
}

void random_interpolate() {
    POINT p0;
    POINT p1;
    for (int i = 0; i < 100; i++) {
        p0.pos[X] = (int)random_float(-400, 400);
        p0.pos[Y] = (int)random_float(-400, 400);
        p1.pos[X] = (int)random_float(-400, 400);
        p1.pos[Y] = (int)random_float(-400, 400);
        p0.color[R] = random_float(0, 1);
        p0.color[G] = random_float(0, 1);
        p0.color[B] = random_float(0, 1);
        p1.color[R] = random_float(0, 1);
        p1.color[G] = random_float(0, 1);
        p1.color[B] = random_float(0, 1);
        draw_line_color(p0, p1);
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
        draw_random();
    }
    else if (test == 1) {
        draw_fan();
    }
    else if (test == 2) {
        draw_coordinate();
    }
    else if (test == 3) {
        random_bresenham();
    }
    else if (test == 4) {
        random_wu();
    }
    else if (test == 5) {
        random_interpolate();
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
		case 't':       test++; if( test > 5 ) test = 0; draw_one_frame = 1;     glutPostRedisplay();    break;
        case 'y':       test--; if (test < 0 ) test = 5; draw_one_frame = 1;     glutPostRedisplay();    break;
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
