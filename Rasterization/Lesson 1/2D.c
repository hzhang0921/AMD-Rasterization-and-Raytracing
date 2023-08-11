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
int times_clicked = 0;
int mouse_x = 3;
int mouse_y = 3;
float colortable[100][4];

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
 * draw_point()
 */
void draw_point( float x, float y )
{
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
void set_color( float r, float g, float b, float a )
{
    glColor4f( r, g, b, a );
}

/*
 * draw_one_point()
 */
void draw_one_point( void )
{
	set_color( 1, 0, 0, 1 );
	draw_point( 0, 0 );
}

/*
 * draw_random_point()
 */
void draw_random_point(void)
{
    set_color(1, 0, 0, 1);
    draw_point((int)random_float(-400,400), (int)random_float(-400,400));
}

void draw_random_points_and_colors(void)
{
    for (int x = -400; x < 400; x++) {
        for (int y = -400; y < 400; y++) {
            if (random_float(0,1.0) > 0.5) {
                set_color(random_float(0,1.0), random_float(0,1.0), random_float(0,1.0), random_float(0,1.0));
                draw_point(x, y);
            }
        }
    }
}

void draw_tv_static(void) 
{
    for (int x = -400; x < 400; x++) {
        for (int y = -400; y < 400; y++) {
            float random_gray_value = random_float(0,1.0);
            set_color(random_gray_value, random_gray_value, random_gray_value, random_gray_value);
            draw_point(x, y);
        }
    }
}

/* 
 * fill in Radius (in pixels)
 */
void draw_in_circle(float radius)
{
    for (int x = -400; x < 400; x++) {
        for (int y = -400; y < 400; y++) {
            float distance = sqrt(x * x + y * y);
            if (distance < radius) {
                set_color(random_float(0,1.0), random_float(0,1.0), random_float(0,1.0), random_float(0,1.0));
                draw_point(x,y);
            }
        }
    }
}


/*
 * R, G, B, A = [0,1]
 */
void fill_color(float r, float g, float b, float a) 
{
    set_color(r, g, b, a);
    for (int x = -400; x < 400; x++) {
        for (int y = -400; y < 400; y++) {
            draw_point(x,y);
        }
    }
}

/*
 * Set RGB, gradients from Left to Right
 */
void fill_gradient_LR(float start_r, float start_g, float start_b, float end_r, float end_g, float end_b)
{
    float increment_r = (end_r - start_r)/800;
    float increment_g = (end_g - start_g)/800;
    float increment_b = (end_b - start_b)/800;
    float track_r = start_r;
    float track_g = start_g;
    float track_b = start_b;
    for (int x = -400; x < 400; x++) {
        track_r += increment_r;
        track_g += increment_g;
        track_b += increment_b;
        for (int y = -400; y < 400; y++) {
            set_color(track_r, track_g, track_b, 1);
            draw_point(x,y);
        }
    }
}

void mandelbrot(float acorner, float bcorner, float xside, float yside)
{
    float a, b;
    float x, y, next_x, next_y; float shade;
    float xinc = xside/800; float yinc = yside/800;
    int count, j, k;
    a = acorner;
    for( j = -400; j < 400; j++ )
    {
        b = bcorner;
        for( k = -400; k < 400; k++ )
            {
                x = 0.0;
                y = 0.0;
                count = 0;
                while(count < 100 && (x * x + y * y) < 4)
                    {
                    next_x = x * x - y * y + a;
                    next_y = 2 * x * y + b;
                    x = next_x;
                    y = next_y;
                    count = count + 1;
                    }
                // shade = count/100.0;
                // glColor4f(1, shade, shade, 1.0);
                glColor4f(colortable[count][0], colortable[count][1], colortable[count][2], colortable[count][3]);
                draw_point(j+0.5, k+0.5); 
                b -= yinc; 
            }
        a += xinc; 
    }
}

/*
 * Draws triangles based on arrays of length 2, [x, y]
 */
void draw_triangle(float* a, float* b, float* c) {
    glVertex2fv(a);
    glVertex2fv(b);
    glVertex2fv(c);
}

/*
 * Divides the triangles into 3 triangles, keeping track of the vertices
 * a is bottom left, b is bottom right, c is middle
 */
void divide_triangle(float* a, float* b, float* c, int m) {
    float v[3][2];
    if(m > 0) {
        for(int j = 0; j < 2; j++) // First loop updates all x values, second loop updates all y values
        {
            v[0][j] = (a[j] + b[j]) / 2; // point between bottom points
            v[1][j] = (a[j] + c[j]) / 2; // point between left and middle points
            v[2][j] = (b[j] + c[j]) / 2; // point between right and middle points
        }
        
        divide_triangle(a, v[0], v[1], m - 1); // bottom left, bottom middle, and left middle points, m - 1 to keep track of recursion
        divide_triangle(c, v[1], v[2], m - 1);
        divide_triangle(b, v[2], v[0], m - 1);
    } else {
        set_color(1, 1, 1, 1);
        draw_triangle(a, b, c);
    }
}

void Sierpinski(float length, int depth)
{
    glBegin(GL_TRIANGLES);
    float start_vertices[3][2] = {{-length/2, -length * sqrt(3)/4}, {length/2, -length * sqrt(3)/4}, {0, length * sqrt(3)/4}};
    divide_triangle(start_vertices[0], start_vertices[1], start_vertices[2], depth);
    glEnd();
}

void Bifurcation(void) 
{
    double x;
    double y = 0;
    double r = 2;
    double inc = 0.004;

    for (int i = 0; i < 800; i++) 
    {
        x = 0.5;
        r += inc;
        for (int j = 0; j < 200; j++)
        {
            x = r * x * (1 - x);
            set_color(1, 1, 1, 1);
            draw_point(x * 800 - 400, y - 400);
        }
        y = y + 1;
    }
}

void x_squared() 
{
    set_color(1, 1, 1, 1);
    for (int x = -400; x < 400; x++) 
    {
        draw_point(x, x*x);
    }
}

void two_to_x() 
{
    set_color(1, 1, 1, 1);
    for (int x = -400; x < 400; x++)
    {
        draw_point(x, pow(2, x));
    }
}

void log_x()
{
    set_color(1, 1, 1, 1);
    for (int x = -400; x < 400; x++)
    {
        draw_point(x, log(x));
    }
}

void sin_func()
{
    set_color(1, 1, 1, 1);
    for (int x = -400; x < 400; x++) 
    {
        set_color(random_float(0,1.0), random_float(0,1.0), random_float(0,1.0), random_float(0,1.0));
        float radians_x = x / 400.0 * 2 * PI;
        draw_point(x, 400 * sin(radians_x + times_clicked) ); //400 is multiplier for clarity
    }
}

void cos_func()
{
    set_color(1, 1, 1, 1);
    for (int x = -400; x < 400; x++)
    {
        float radians_x = x / 400 * 2 * PI;
        draw_point(x, 400 * cos(radians_x)); //400 is multiplier for clarity
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
	
    Sierpinski(400, 8);
	
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
 * Mouse routine
 */
void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        mouse_x = x;
        mouse_y = y;
        printf("%d ", mouse_x);
        printf("%d \n", mouse_y);
        times_clicked += 1;
        // Refresh the screen
        draw_one_frame = 1;

        // Redraw the scene
        glutPostRedisplay();
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
    glutMouseFunc( mouse );

    /*
     * setup OpenGL state
     */
    glClearColor( 0.0, 0.0, 0.0, 1.0 );
    gluOrtho2D( -half_window_size, half_window_size ,-half_window_size, half_window_size );
    glPointSize( 2.0 );
    glColor4f( 1.0, 0.0, 0.0, 1.0 );

	if( Mojave_WorkAround )
	{
		// Necessary for Mojave. Has to be different dimensions than in glutInitWindowSize();
		glutReshapeWindow( half_window_size * 2, half_window_size * 2 );
	}
	
    for (int i = 0; i < 100; i++)
    {
        colortable[i][0] = random_float(0, 1);
        colortable[i][1] = random_float(0, 1);
        colortable[i][2] = random_float(0, 1);
        colortable[i][3] = 1;
    }
    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}
