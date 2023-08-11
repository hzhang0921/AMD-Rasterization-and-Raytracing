/*
 *
 * point.c - simple GLUT app that draws one frame with a single point at origin
 *
 * To build:  gcc -framework OpenGL -framework GLUT 2D.c -o 2D
 *
 */
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
int current[800][800];
int next[800][800];
int generations = 0;

/*************************************************************************/
/* utility functions                                                     */
/*************************************************************************/

/*
 * random_float()
 */
float random_float( float low, float high ) {
    return( (float)(low + (rand()/(float)RAND_MAX)*(high - low)) );
}

int random_int() {
    if (random_float(0, 1) > 0.5) {
        return 1;
    }
    else {
        return 0;
    }
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
 * set_color()
 */
void set_color( float r, float g, float b, float a ) {
    glColor4f( r, g, b, a );
}

/*
 * draw_one_point()
 */
void draw_one_point( void ) {
	set_color( 1, 1, 1, 1 );
	draw_point( 0, 0 );
}

void customize_current(void) {
    for (int y = 0; y < 11; y++) {
        for (int x = 0; x < 38; x++) {
            current[x][800 - y] = 0;
        }
    }
    current[1][795] = 1;
    current[1][794] = 1;
    current[2][795] = 1;
    current[2][794] = 1;
    current[11][795] = 1;
    current[11][794] = 1;
    current[11][793] = 1;
    current[12][796] = 1;
    current[12][792] = 1;
    current[13][797] = 1;
    current[13][791] = 1;
    current[14][797] = 1;
    current[800 - 9][14] = 1;
    current[800 - 6][15] = 1;
    current[800 - 4][16] = 1;
    current[800 - 8][16] = 1;
    current[800 - 5][17] = 1;
    current[800 - 6][17] = 1;
    current[800 - 7][17] = 1;
    current[800 - 6][18] = 1;
    current[800 - 3][21] = 1;
    current[800 - 4][21] = 1;
    current[800 - 5][21] = 1;
    current[800 - 3][22] = 1;
    current[800 - 4][22] = 1;
    current[800 - 5][22] = 1;
    current[800 - 2][23] = 1;
    current[800 - 6][23] = 1;
    current[800 - 2][25] = 1;
    current[800 - 3][25] = 1;
    current[800 - 6][25] = 1;
    current[800 - 7][25] = 1;
    current[800 - 3][35] = 1;
    current[800 - 4][35] = 1;
    current[800 - 3][36] = 1;
    current[800 - 4][36] = 1;

}

void init_current(void) {
    for (int y = 0; y < 800; y++) {
        for (int x = 0; x < 800; x++) {
            current[y][x] = random_int();
        }
    }
    customize_current();
}

void draw_current(void) {
    for (int y = 0; y < 800; y++) {
        for (int x = 0; x < 800; x++) {
            if (current[y][x] == 1) {
                draw_point(y, x);
            }
        }
    }
}


int count_neighbors(int y, int x) {
    int num_neighbors = current[y-1][x-1] + current[y][x-1] + current[y+1][x-1] + current[y-1][x] + current[y+1][x] + current[y-1][x+1] + current[y][x+1] + current[y+1][x+1];
    return num_neighbors;
}

void calc_next() {
    for (int y = 1; y < 799; y++) {
        for (int x = 1; x < 799; x++) {
            int num_neighbors = count_neighbors(y, x);
            if (current[y][x] == 1) {
                {
                    if (num_neighbors == 2 || num_neighbors == 3) {
                        next[y][x] = 1;
                    }
                    else {
                        next[y][x] = 0;
                    }
                }
            }
            else{
                if (num_neighbors == 3) {
                    next[y][x] = 1;
                }
                else {
                    next[y][x] = 0;
                }
            }
        }
    }
}

void draw_game() {
    draw_current();
    calc_next();
    for (int y = 1; y < 799; y++) {
        for (int x = 1; x < 799; x++) {
            current[y][x] = next[y][x];
        }
    }
    generations++;
    printf("%d", generations);
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
	
    draw_game();
	
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
        case 'a':       		                        draw_one_frame = 1;     glutPostRedisplay();    break;
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
    gluOrtho2D( 0, 2 * half_window_size , 0, 2 * half_window_size );
    glPointSize( 2.0 );
    glColor4f( 1.0, 1.0, 1.0, 1.0 );

	if( Mojave_WorkAround )
	{
		// Necessary for Mojave. Has to be different dimensions than in glutInitWindowSize();
		glutReshapeWindow( half_window_size * 2, half_window_size * 2 );
	}

    init_current();
    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}
