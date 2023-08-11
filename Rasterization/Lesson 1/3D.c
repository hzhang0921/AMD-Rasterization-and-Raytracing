/*
 *
 * point.c - simple GLUT app that draws one frame with a single point at origin
 *
 * To build:  gcc -framework OpenGL -framework GLUT 3D.c -o 3D
 *
 */
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

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int half_window_size = 400;
int Mojave_WorkAround = 1;
int draw_one_frame = 1;
int test = 0;

/*************************************************************************/
/* utility functions                                                     */
/*************************************************************************/
/*
 * random_float()
 */
float random_float( float low, float high )
{
    return( (float)(low + (rand()/(float)RAND_MAX)*(high - low)) );
}

/*
 * set_color()
 */
void set_color( float r, float g, float b, float a )
{
    glColor4f( r, g, b, a );
}

void Lorenz() 
{
    double x = 45;
    double y = 45;
    double z = 45;
    double a = 10.0;
    double b = 28.0;
    double c = 8.0 / 3.0;
    double t = 0.01;
    int lorenzIterationCount = 50000;
    //Iterate and update x,y and z locations
    //based upon the Lorenz equations
    glBegin(GL_POINTS);
    for (int i = 0; i < lorenzIterationCount; i++ ){
        double xt = x + t * a * (y - x);
        double yt = y + t * (x * (b - z) - y);
        double zt = z + t * (x * y - c * z);
        glVertex3d(x, y, z);
        x = xt;
        y = yt;
        z = zt;
    }
    glEnd();
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
	
	Lorenz();
	
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
        case 'a':       								draw_one_frame = 1;     glutPostRedisplay();    break;
		case 't':       test++; if( test > 1) test = 0; draw_one_frame = 1;     glutPostRedisplay();    break;
		case 'q':       exit( 0 );                             											break;
        case '\033':    exit( 0 );																		break;
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
    glutCreateWindow( "Lorenz Attractor" );

    /*
     * setup display(), Key() funtions
     */
    glutDisplayFunc( display );
    glutKeyboardFunc( Key );

    /*
     * setup OpenGL state
     */
    glClearColor( 0.0, 0.0, 0.0, 1.0 );
    /* changed from 2D to 3D */
    glOrtho( -half_window_size, half_window_size ,-half_window_size, half_window_size, -half_window_size, half_window_size);
    glPointSize( 3.0 );
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