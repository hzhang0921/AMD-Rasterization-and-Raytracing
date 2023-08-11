
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
} POINT; // if there is no POINT , you would have to put struct point


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
float cur_heading = 90; //straight up

POINT tcur = {.pos[X] = 0, .pos[Y] = 0};


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
 * draw_point_struct()
 */
void draw_point_struct(POINT p) {
    glBegin(GL_POINTS);
    glVertex2f(p.pos[X], p.pos[Y]); // + or - to translate x and y by pixels
    glEnd();
}

/*
 * set_struct_color()
 */
void set_color( float red, float green, float blue, float a) {
    glColor4f(red, green, blue, a);
}

/*
 * draw_point_struct()
 */
void draw_whitepoint_struct(POINT p) {
    /*
    * draw points
    */
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
    glVertex2f(p.pos[X], p.pos[Y]); // + or - to translate x and y by pixels
    glEnd();
}

/* 
 * draw_horizontal()
 */
void draw_horizontal(POINT p0, POINT p1) {
    POINT p;
    if (p0.pos[Y] != p1.pos[Y]) {
        return;
    }
    if (p0.pos[X] > p1.pos[X]) {
        SWAPPOINT(p0, p1);
    }
    for (int x = p0.pos[X]; x < p1.pos[X]; x++) {
        p.pos[X] = x;
        p.pos[Y] = p0.pos[Y];
        draw_whitepoint_struct(p);
    }
}

/*
 * draw_vertical()
 */
void draw_vertical(POINT p0, POINT p1) {
    POINT p;
    if (p0.pos[X] != p1.pos[X]) {
        return;
    }
    if (p0.pos[Y] > p1.pos[Y]) {
        SWAPPOINT(p0, p1);
    }
    for (int y = p0.pos[Y]; y < p1.pos[Y]; y++) {
        p.pos[Y] = y;
        p.pos[X] = p0.pos[X];
        draw_whitepoint_struct(p);
    }
}


/*
 * draw_line_colorless
 */
void draw_line_colorless(POINT p0, POINT p1) {
    float dx = p1.pos[X] - p0.pos[X];
    float dy = p1.pos[Y] - p0.pos[Y];
    int steep = 0;
    if (fabs(dy) > fabs(dx)) {
        steep = 1;
    }
    float slope = dx/dy;
    if (steep == 1) {
        if (p0.pos[Y] > p1.pos[Y]) {
            SWAPPOINT(p0, p1);
        }
        for (int y = p0.pos[Y]; y < p1.pos[Y]; y++) {
            draw_point_struct(p0); // change this to draw_point_struct(p0) for gradient
            p0.pos[X] += slope;
            p0.pos[Y] = y+1;
        }
    }
    else {
        float slope = dy/dx;
        if (p0.pos[X] > p1.pos[X]) {
            SWAPPOINT(p0, p1);
        }
        for (int x = p0.pos[X]; x < p1.pos[X]; x++) {
            draw_point_struct(p0); // change this to draw_point_struct(p0) for gradient
            p0.pos[Y] += slope;
            p0.pos[X] = x+1;
        }    
    }
}

/*
 * draw_line_color_transition()
 */
void draw_line_color_transition(POINT p0, POINT p1) {
    float dx = p1.pos[X] - p0.pos[X];
    float dy = p1.pos[Y] - p0.pos[Y];
    float color_delta[4];
    float color_inc[4];
    subtract_vector(p1.color, p0.color, color_delta);
    int steep = 0;
    if (fabs(dy) > fabs(dx)) {
        steep = 1;
    }
    float slope = dx/dy;
    scalar_multiply_vector(1/slope, color_delta, color_inc);
    if (steep == 1) {
        if (p0.pos[Y] > p1.pos[Y]) {
            SWAPPOINT(p0, p1);
        }
        for (int y = p0.pos[Y]; y < p1.pos[Y]; y++) {
            draw_point_struct(p0); // change this to draw_point_struct(p0) for gradient
            add_vector(p0.color, color_inc, p0.color);
            p0.pos[X] += slope;
            p0.pos[Y] = y+1;
        }
    }
    else {
        float slope = dy/dx;
        subtract_vector(p1.color, p0.color, color_delta);
        scalar_multiply_vector(1/slope, color_delta, color_inc);
        if (p0.pos[X] > p1.pos[X]) {
            SWAPPOINT(p0, p1);
        }
        for (int x = p0.pos[X]; x < p1.pos[X]; x++) {
            draw_point_struct(p0); // change this to draw_point_struct(p0) for gradient
            add_vector(p0.color, color_inc, p0.color);
            p0.pos[Y] += slope;
            p0.pos[X] = x+1;
        }    
    }
}

/*
 * draw_line()
 */
void draw_line(POINT p0, POINT p1) {
    if (p0.pos[Y] == p1.pos[Y]) {
        draw_horizontal(p0, p1);
    }
    else if(p0.pos[X] == p1.pos[X]) {
        draw_vertical(p0, p1);
    }
    else {
        draw_line_colorless(p0, p1);
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
        draw_line(p0, p1);
    }
}

/*************************************************************************/
/* TURTLE Driver functions                                               */
/*************************************************************************/
/* 
 * void turtle_left()
 */
void turtle_left (float angle) {
    cur_heading += angle;
    while (cur_heading < 0) {
        cur_heading += 360;
    }
    while (cur_heading >= 360) {
        cur_heading -= 360;
    }
}

/*
 * void turtle_left()
 */
void turtle_right (float angle) {
    cur_heading -= angle;
    while (cur_heading < 0) {
        cur_heading += 360;
    }
    while (cur_heading >= 360) {
        cur_heading -= 360;
    }
}

/*
 * void turtle_forward() 
 */
void turtle_forward(float length) {
    POINT p1;
    p1.pos[X] = tcur.pos[X] + (length * cos(cur_heading / 360 * 2 * PI));
    p1.pos[Y] = tcur.pos[Y] + (length * sin(cur_heading / 360 * 2 * PI));
    draw_line(tcur, p1);
    tcur.pos[X] = p1.pos[X];
    tcur.pos[Y] = p1.pos[Y];
}

/*
 * void turtle_backward() 
 */
void turtle_backward(float length) {
    POINT p1;
    p1.pos[X] = tcur.pos[X] + (length * -cos(cur_heading / 360 * 2 * PI));
    p1.pos[Y] = tcur.pos[Y] + (length * -sin(cur_heading / 360 * 2 * PI));
    draw_line(tcur, p1);
    tcur.pos[X] = p1.pos[X];
    tcur.pos[Y] = p1.pos[Y];
}

/*
 * void turtle_home()
 */
void turtle_home() {
    cur_heading = 90;
    tcur.pos[X] = 0;
    tcur.pos[Y] = 0;
}


/* 
 * turtle_draw_square()
 */
void turtle_draw_square() {
    turtle_forward(100);
    turtle_right(90);
    turtle_forward(100);
    turtle_right(90);
    turtle_forward(100);
    turtle_right(90);
    turtle_forward(100);
    turtle_right(90);
    turtle_home();
}

/*
 * turtle_draw_circle()
 */
void turtle_draw_circle() {
    for (int i = 0; i < 360; i++) {
        turtle_forward(5);
        turtle_right(1);
    }
    turtle_home();
}

/*
 * turtle_draw_spiral()
 */
void turtle_draw_spiral() {
    for (int i = 0; i < 100; i++) {
        turtle_forward(5 + i);
        turtle_right(15);
    }
    turtle_home();
}

/*
 * turtle_draw_star()
 */
void turtle_draw_star() {
    turtle_right(27);
    for (int i = 0; i < 5; i++) {
        turtle_forward(100);
        turtle_right(144);
    }
    turtle_home();
}

/*
 * turtle_draw_V()
 */
void turtle_draw_v() {
    turtle_left(30);
    turtle_forward(100);
    turtle_backward(100);
    turtle_right(60);
    turtle_forward(100);
    turtle_home();
}

/*
 * turtle_draw_mountain()
 */
void turtle_draw_mountain() {
    turtle_right(90);
    turtle_forward(100);
    turtle_left(45);
    turtle_forward(100);
    turtle_right(90);
    turtle_forward(100);
    turtle_left(45);
    turtle_forward(100);
    turtle_home();
}

/*
 * turtle_draw_triangle()
 */
void turtle_draw_triangle() {
    turtle_right(90);
    for (int i = 0; i < 3; i++) {
        turtle_forward(100);
        turtle_right(240);
    }
    turtle_home();
}

/* 
 * turtle_draw_pentagon()
 */
void turtle_draw_pentagon() {
    turtle_right(90);
    for (int i = 0; i < 5; i++) {
    turtle_forward(100);
    turtle_right(288);
    }
    turtle_home();
}

/*
 * turtle_draw_polygon()
 */
void turtle_draw_polygon(int sides) {
    turtle_right(90);
    int total_degrees = (sides - 2) * 180;
    for (int i = 0; i < sides; i++) {
        turtle_forward(50);
        turtle_right(180 + (total_degrees/sides));
    }
    turtle_home();
}

/* 
 * random_turtle()
 */
void random_turtle() {
    POINT p0;
    POINT p1;
    for (int i = 0; i < 20; i++) {
        set_color(random_float(0, 1), random_float(0, 1), random_float(0, 1), 1);
        tcur.pos[X] = random_float(-400, 400);
        tcur.pos[Y] = random_float(-400, 400);
        int sides = random_float(3, 20);
        turtle_draw_polygon(sides);
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
        turtle_draw_square();
    }
    else if (test == 1) {
        turtle_draw_circle();
    }
    else if (test == 2) {
        turtle_draw_spiral();
    }
    else if (test == 3) {
        turtle_draw_star();
    }
    else if (test == 4) {
        turtle_draw_v();
    }
    else if (test == 5) {
        turtle_draw_mountain();
    }
    else if (test == 6) {
        turtle_draw_triangle();
    }
    else if (test == 7) {
        turtle_draw_pentagon();
    }
    else if (test == 8) {
        turtle_draw_polygon(10);
    }
    else if (test == 9) {
        random_turtle();
    }
    else if (test == 10) {
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
		case 't':       test++; if( test > 10 ) test = 0; draw_one_frame = 1;     glutPostRedisplay();    break;
        case 'y':       test--; if (test < 0 ) test = 10; draw_one_frame = 1;     glutPostRedisplay();    break;
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
