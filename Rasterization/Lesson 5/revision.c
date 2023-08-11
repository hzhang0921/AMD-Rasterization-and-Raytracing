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
#include <stdbool.h>
#include <strings.h> // library for strings
#include <sys/time.h>



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
#define MILLION             1000000
#define MESH_WIDTH 32
#define MESH_HEIGHT 32
#define X 0
#define Y 1
#define Z 2
#define W 3
#define R 0
#define G 1
#define B 2
#define A 3
#define U 0
#define V 1

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
    float normal_vector[4]; // normal vector for the triangle 
} TRIANGLE;

typedef struct timer {
    struct timeval start;
    struct timeval end;
} TIMER;

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int half_window_size = 400; // used for standarizing coordinate plane
int Mojave_WorkAround = 1; // used for a bug in Mojave
int draw_one_frame = 1;
int test = 0;
int textureIndex = 0;
int texOrRGBA = 0;
int current[800][800];
float depth[800][800];
float colorbuffer[800][800][4];
float light_vector[4] = {-1, -1, 1, 0}; // sets the infinite light vector to 1, 1, 1
IMAGE Textures[20]; // allocates memory for 20 images of 4096 x 4096 RGBA pixels
POINT vertex_list[1000000]; // vertices of the triangles
TRIANGLE triangle_list[1000000]; // list of triangles

float y_angle = 0;
int num_triangles = 0;
int num_verts = 0;

float model_center_x = 0;
float model_center_y = 0;
float model_center_z = 0;

TIMER vertex_timer;
TIMER pixel_timer;
TIMER frame_timer;
TIMER gl_timer;

float vertex_time;
float pixel_time;
float frame_time;
float gl_time;

int total_verts = 0;
int total_tris = 0;

int show_stats = 0;

/*************************************************************************/
/* Texture functions                                                     */
/*************************************************************************/
void texture_sample(IMAGE *tex, float s, float t, float tex_color[4])
{                           // passing address is good for memory management. 4 bytes vs 16 MB
    s = CLAMP(s, 0.0, 1.0);
    t = CLAMP(t, 0.0, 1.0);
    int v = s * (tex->height - 1);
    int u = t * (tex->width - 1);

    tex_color[R] = tex->data[v][u][R] / 255.0; 
    tex_color[G] = tex->data[v][u][G] / 255.0; 
    tex_color[B] = tex->data[v][u][B] / 255.0; 
    tex_color[A] = tex->data[v][u][A] / 255.0; 
}

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

void cross_vector(float a[4], float b[4], float result[4])
{
    result[0] = a[Y] * b[Z] - a[Z] * b[Y]; // x
    result[1] = a[Z] * b[X] - a[X] * b[Z]; // y
    result[2] = a[X] * b[Y] - a[Y] * b[X]; // z
    result[3] = 0;
}

float dotProduct(float a[4], float b[4])
{
    return (a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z]);
}

float magnitude(float a[4])
{
    return (sqrt(a[X] * a[X] + a[Y] * a[Y] + a[Z] * a[Z]));
}

void normalize(float a[4], float result[4])
{
    float mag = magnitude(a);
    if (mag == 0) // edge case for a zero vector
    {
        result[X] = 0;
        result[Y] = 0;
        result[Z] = 0;
    }
    result[X] = a[X] / mag;
    result[Y] = a[Y] / mag;
    result[Z] = a[Z] / mag;
}

void vector_cross(float a[4], float b[4], float result[4])
{
    result[0] = a[1] * b[2] - a[2] * b[1]; // x
    result[1] = a[2] * b[0] - a[0] * b[2]; // y
    result[2] = a[0] * b[1] - a[1] * b[0]; // z
    result[3] = 0;
}

/*************************************************************************/
/* buffer functions                                                      */
/*************************************************************************/
/*
 * clear_colorbuffer() resets all pixels to a default value of black (0, 0, 0, 0)
 */
void clear_colorbuffer() {
    for (int y = 0; y < 800; y++) {
        for (int x = 0; x < 800; x++) {
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
    for (int y = 0; y < 800; y++) {
        for (int x = 0; x < 800; x++) {
            depth[x][y] = 1000;
        }
    }
}

/*
 * draw_buffer()
 */
void draw_colorbuffer() {
    for (int y = 0; y < 800; y++) {
        for (int x = 0; x < 800; x++) {
            glBegin(GL_POINTS);
            glColor4f(colorbuffer[x][y][R], colorbuffer[x][y][G], colorbuffer[x][y][B], colorbuffer[x][y][A]);
            glVertex2f(x - 400, y - 400); // subtract 400 since color buffer is from 0 - 800
            glEnd();
        }
    }
}

/*************************************************************************/
/* texture functions                                                     */
/*************************************************************************/
void read_ppm(IMAGE *image, char *name)
{
    FILE *fp;
    int c;
    char buffer[16];
    int w;
    int h;
    int r, g, b;
    int max;
    unsigned char *data, *p;
    float scale;

    if ((fp = fopen(name, "rb")) == NULL)
    {
        image->width = 0;
        image->height = 0;
        return;
    }
    fscanf(fp, "%s\n", buffer);
    /*
     * skip lines that start with '#'
     */
    while (1)
    {
        c = fgetc(fp);
        if (c == '#')
        {
            do
            {
                c = fgetc(fp);
            } while (c != '\n');
        }
        else
        {
            ungetc(c, fp);
            break;
        }
    }
    fscanf(fp, "%d %d\n", &w, &h);
    fscanf(fp, "%d\n", &max);
    scale = 255.0 / max;

    image->width = w;
    image->height = h;

    if (buffer[0] == 'P' && buffer[1] == '3')
    {
        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++)
            {
                fscanf(fp, "%d %d %d", &r, &g, &b);

                r *= scale;
                g *= scale;
                b *= scale;

                image->data[j][i][R] = CLAMP(r, 0, 255);
                image->data[j][i][G] = CLAMP(g, 0, 255);
                image->data[j][i][B] = CLAMP(b, 0, 255);
                image->data[j][i][A] = 255;
            }
        }
    }
    else
    {
        data = (unsigned char *)malloc(w * h * 3);

        if (data == NULL)
        {
            image->width = 0;
            image->height = 0;

            fclose(fp);
            return;
        }

        p = data;

        fread((void *)data, 3, w * h, fp);

        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++)
            {
                r = *p++ * scale;
                g = *p++ * scale;
                b = *p++ * scale;

                image->data[h - j - 1][i][R] = CLAMP(r, 0, 255);
                image->data[h - j - 1][i][G] = CLAMP(g, 0, 255);
                image->data[h - j - 1][i][B] = CLAMP(b, 0, 255);
                image->data[h - j - 1][i][A] = 255;
            }
        }
        free(data);
    }

    fclose(fp);
}

/*
 * write_ppm()
 */
void write_ppm(IMAGE *image, char *name, int id)
{
    FILE *fp;
    int c;
    char filename[1024];
    int w;
    int h;
    int r, g, b;
    int max;
    unsigned char *data, *p;
    float scale;

    sprintf(filename, "%s_%d.ppm", name, id);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        printf("write_ppm(%s) failed to open\n", filename);
        return;
    }

    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", image->width, image->height);
    fprintf(fp, "255\n");

    for (int j = 0; j < image->height; j++)
    {
        for (int i = 0; i < image->width; i++)
        {
            fprintf(fp, "%d %d %d ", image->data[j][i][R], image->data[j][i][G], image->data[j][i][B]);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
}

void draw_checkerboard(IMAGE *image) {  // *image is where the IMAGE is saved
    image->width = 800;
    image->height = 800;
    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            int size = 50;
            if ((x / size) % 2 == 0) {
                if ((y / size) % 2 == 0) {
                    UCset_vector(image->data[x][y], 255, 255, 255, 255);
                }
                else {
                    UCset_vector(image->data[x][y], 0, 0, 0, 255);
                }
            }
            else {
                if ((y / size) % 2 == 0) {
                    UCset_vector(image->data[x][y], 0, 0, 0, 255);
                }
                else {
                    UCset_vector(image->data[x][y], 255, 255, 255, 255); 
                }
            }
        }
    }
    // write_ppm(&Textures[0], "../Images/checkerboard", 0); // saves the IMAGE to a given file location with file name
}


/*************************************************************************/
/* utility functions                                                     */
/*************************************************************************/
/*
 * start_timer()
 */
void start_timer( TIMER *timer )
{
    gettimeofday( &timer->start, NULL );
}

/*
 * stop_timer()
 */
void stop_timer( TIMER *timer )
{
    double usec;

    gettimeofday( &timer->end, NULL);
}

/*
 * elapsed_time()
 */
double elapsed_time( TIMER *timer )
{
    double usec;
    
    // Calculating total time taken by the program.
    usec = ((timer->end.tv_sec * MILLION) + timer->end.tv_usec) - ((timer->start.tv_sec * MILLION) + timer->start.tv_usec);
    
    return( usec/MILLION );
}


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
    int x, y;
    x = p.pos[X] + 400;
    y = p.pos[Y] + 400;
    float angle;
    if (x >= 0 && x < 800 && y >= 0 && y < 800) { // checks for out of bounds values
        if (p.world[Z] < depth[x][y]) { // checks if the point is closer than the previous point
            depth[x][y] = p.world[Z];
            if (texOrRGBA == 0) {
                set_vector(colorbuffer[x][y], p.color[R], p.color[G], p.color[B], p.color[A]); // sets the color buffer to the color of the point
            }
            else if (texOrRGBA == 1) {
                texture_sample(&Textures[textureIndex], p.tex_coords[U], p.tex_coords[V], p.color); // samples the texture
                set_vector(colorbuffer[x][y], p.color[R], p.color[G], p.color[B], p.color[A]); // sets the color buffer to the color of the point
            }
        }
    }
}

void draw_line(POINT p0, POINT p1, bool is_vertical, bool is_horizontal) { // general case to interpolate between two points
    float pos_delta[4];
    float color_delta[4];
    float tex_delta[4];
    float pos_inc[4];
    float color_inc[4];
    float tex_inc[4];
    POINT p;
    p = p0;

    subtract_vector( p1.pos, p0.pos, pos_delta);
    subtract_vector( p1.color, p0.color, color_delta);
    subtract_vector( p1.tex_coords, p0.tex_coords, tex_delta);

    float inv_slope;
    if (is_horizontal) {
        inv_slope = 1.0/pos_delta[X];
    } else if (is_vertical) {
        inv_slope = 1.0/pos_delta[Y];
    }
    
    scalar_multiply_vector( inv_slope, pos_delta, pos_inc);
    scalar_multiply_vector( inv_slope, color_delta, color_inc);
    scalar_multiply_vector( inv_slope, tex_delta, tex_inc);
    
    if(is_horizontal) {
        for( int x = p0.pos[X]; x <= p1.pos[X]; x++ ) {
            draw_RGBApoint_struct(p);
            add_vector(p.pos, pos_inc, p.pos);
            add_vector(p.color, color_inc, p.color);
            add_vector(p.tex_coords, tex_inc, p.tex_coords);
        }
    } else if (is_vertical) {
        for( int y = p0.pos[Y]; y <= p1.pos[Y]; y++ ) {
            draw_RGBApoint_struct(p);
            add_vector(p.pos, pos_inc, p.pos);
            add_vector(p.color, color_inc, p.color);
            add_vector(p.tex_coords, tex_inc, p.tex_coords);
        }
    }
}

void draw_line_color( POINT p0, POINT p1 ) { // selection of specific interpolation mode
    if( (int)(p1.pos[X] - p0.pos[X]) == 0 && (int)(p1.pos[Y] - p0.pos[Y])== 0 ) { // degenerate case
        draw_RGBApoint_struct( p0 );
    } else if( (int)p1.pos[Y] == (int)p0.pos[Y] ) { // horizontal line
        if( p0.pos[X] > p1.pos[X] ) {
            SWAPPOINT(p0,p1);
        }
        draw_line(p0, p1, false, true);
    } else if( (int)p1.pos[X] == (int)p0.pos[X] ) { // vertical line
        if( p0.pos[Y] > p1.pos[Y] ) {
            SWAPPOINT(p0,p1);
        }
        draw_line(p0, p1, true, false);
    } else if( fabs(p1.pos[X] - p0.pos[X]) > fabs(p1.pos[Y] - p0.pos[Y]) ) {
        if( p0.pos[X] > p1.pos[X] ) {
            SWAPPOINT(p0,p1);
        }
        draw_line(p0, p1, false, true);
    } else {
        if( p0.pos[Y] > p1.pos[Y] ) {
            SWAPPOINT(p0,p1);
        }
        draw_line(p0, p1, true, false);
    }
}

/*************************************************************************/
/* Triangle utilities                                                    */
/*************************************************************************/
TRIANGLE calc_normal(TRIANGLE t) {
    float e0[4], e1[4], normal[4];
    POINT p0, p1, p2;
    p0 = vertex_list[t.index_list[0]]; // access the first vertex of the triangle
    p1 = vertex_list[t.index_list[1]];
    p2 = vertex_list[t.index_list[2]];
    subtract_vector(p1.world, p0.world, e0); // Edge vector for edge 1
    subtract_vector(p2.world, p0.world, e1); // Edge vector for edge 2
    vector_cross(e0, e1, normal); // sets the normal vector in normal
    normalize(normal, normal); // normalizes the normal vector
    t.normal_vector[X] = normal[X]; // copies the normal vector over to the triangle
    t.normal_vector[Y] = normal[Y];
    t.normal_vector[Z] = normal[Z];
    return t;
}

void set_triangles(void) {
    num_triangles = 0;
    for (int col = 0; col < MESH_WIDTH - 1; col++) { // sets triangles
        for (int row = 0; row < MESH_HEIGHT - 1; row++) {
            triangle_list[num_triangles].index_list[0] = (row * MESH_WIDTH) + col; // flat top triangle
            triangle_list[num_triangles].index_list[1] = (row * MESH_WIDTH) + col + 1;
            triangle_list[num_triangles].index_list[2] = ((row + 1) * MESH_WIDTH) + col + 1;
            triangle_list[num_triangles] = calc_normal(triangle_list[num_triangles]); // sets the normal vector for a given triangle
            triangle_list[num_triangles+1].index_list[0] = (row * MESH_WIDTH) + col; // flat bot triangle
            triangle_list[num_triangles+1].index_list[1] = ((row + 1) * MESH_WIDTH) + col + 1;
            triangle_list[num_triangles+1].index_list[2] = ((row + 1) * MESH_WIDTH) + col;
            triangle_list[num_triangles+1] = calc_normal(triangle_list[num_triangles+1]); // sets the normal vector for the second triangle
            num_triangles += 2;
        }
    }
}

// float triangle_area(POINT p0, POINT p1, POINT p2) {
//      return(p2.pos[X] - p0.pos[X]) * (p1.pos[Y] - p0.pos[Y]) - (p2.pos[Y] - p2.pos[Y]) * (p1.pos[X] - p0.pos[X]);
// }


// void calculate_barycentric(POINT p0, POINT p1, POINT p2, POINT p, float *alpha, float *beta, float *gamma) {
//     float area = triangle_area(p0, p1, p2);
//     if (area <= 0) {
//         *alpha = 0;
//         *beta = 0;
//         *gamma = 0;
//         return;
//     }
//     *alpha = triangle_area(p, p1, p2) / area;
//     *beta = triangle_area(p0, p, p2) / area;
//     *gamma = triangle_area(p0, p1, p) / area;
// }



// void draw_barycentrictriangle(POINT p0, POINT p1, POINT p2) {
//     float alpha, beta, gamma;
//     POINT p;
//     int minX = MIN(MIN(p0.pos[X], p1.pos[X]), p2.pos[X]);
//     int minY = MIN(MIN(p0.pos[Y], p1.pos[Y]), p2.pos[Y]);
//     int maxX = MAX(MAX(p0.pos[X], p1.pos[X]), p2.pos[X]);
//     int maxY = MAX(MAX(p0.pos[Y], p1.pos[Y]), p2.pos[Y]);

//     for (int y = 0; y < maxY - minY; y++) {
//         for (int x = 0; x < maxX - minX; x++) {
//             p.pos[X] = x + minX;
//             p.pos[Y] = y + minY;
//             calculate_barycentric(p0, p1, p2, p, &alpha, &beta, &gamma);
//             if (alpha >= 0 && beta >= 0 && gamma >= 0) {
//                 p.pos[Z] = alpha * p0.pos[Z] + beta * p1.pos[Z] + gamma * p2.pos[Z];
//                 p.color[R] = alpha * p0.color[R] + beta * p1.color[R] + gamma * p2.color[R];
//                 p.color[G] = alpha * p0.color[G] + beta * p1.color[G] + gamma * p2.color[G];
//                 p.color[B] = alpha * p0.color[B] + beta * p1.color[B] + gamma * p2.color[B];
//                 p.color[A] = alpha * p0.color[A] + beta * p1.color[A] + gamma * p2.color[A];
//                 p.tex_coords[U] = alpha * p0.tex_coords[U] + beta * p1.tex_coords[U] + gamma * p2.tex_coords[U];
//                 p.tex_coords[V] = alpha * p0.tex_coords[V] + beta * p1.tex_coords[V] + gamma * p2.tex_coords[V];
//                 draw_RGBApoint_struct(p);
//             }
//             else {
//                 printf("alpha: %f, beta: %f, gamma: %f\n", alpha, beta, gamma);
//             }
//         }
//     }       
// }

/*
 * edgeFunction()
 */
float edgeFunction( float a[4], float b[4], float c[4] )
{
    return( c[X] - a[X]) * (b[Y] - a[Y]) - (c[Y] - a[Y]) * (b[X] - a[X] );
}

/*
 * draw_triangle_barycentric()
 */
void draw_barycentrictriangle( POINT p0, POINT p1, POINT p2 )
{
    POINT* v0;
    POINT* v1;
    POINT* v2;
    v0 = &p0;
    v1 = &p1;
    v2 = &p2;
    int minx = MIN(MIN(p0.pos[X], p1.pos[X]), p2.pos[X]);
    int miny = MIN(MIN(p0.pos[Y], p1.pos[Y]), p2.pos[Y]);
    int maxx = MAX(MAX(p0.pos[X], p1.pos[X]), p2.pos[X]);
    int maxy = MAX(MAX(p0.pos[Y], p1.pos[Y]), p2.pos[Y]);
    float   inv_area = 1.0/edgeFunction(v0->pos, v1->pos, v2->pos);
    float   w[4];
    POINT   p;
    
    /*
     * for all the points in triangle bounding box
     */
    for( int y = miny-1; y < maxy+2; y++ )
    {
        for( int x = minx-1; x < maxx+2; x++ )
        {
            set_vector( p.pos, x, y, 0, 1 );
       
            /*
             * compute barycentric weights
             */
            w[0] = edgeFunction(v1->pos, v2->pos, p.pos);
            w[1] = edgeFunction(v2->pos, v0->pos, p.pos);
            w[2] = edgeFunction(v0->pos, v1->pos, p.pos);
            w[3] = 0.0;
            
            if ((w[0] + w[1] + w[2]) < 0) {
                w[0] = -w[0];
                w[1] = -w[1];
                w[2] = -w[2];
                inv_area *= -1;
            }

            if( w[0] >= 0 && w[1] >= 0 && w[2] >= 0 )
            {
                /*
                 * if point is inside triangle, compute barycentric weighting of vertex attributes (e.g. z, color, s, t)
                 */
				scalar_multiply_vector( inv_area, w, w );
				
                p.pos[Z] = w[0] * v0->pos[Z] + w[1] * v1->pos[Z] + w[2] * v2->pos[Z];

				/*
				 * if per-vertex lighting, just copy vertex color to point
				 */
					p.color[R]  = w[0] * v0->color[R] + w[1] * v1->color[R] + w[2] * v2->color[R];
					p.color[G]  = w[0] * v0->color[G] + w[1] * v1->color[G] + w[2] * v2->color[G];
					p.color[B]  = w[0] * v0->color[B] + w[1] * v1->color[B] + w[2] * v2->color[B];
					p.color[A]  = w[0] * v0->color[A] + w[1] * v1->color[A] + w[2] * v2->color[A];
                    p.tex_coords[U] = w[0] * v0->tex_coords[U] + w[1] * v1->tex_coords[U] + w[2] * v2->tex_coords[U];
                    p.tex_coords[V] = w[0] * v0->tex_coords[V] + w[1] * v1->tex_coords[V] + w[2] * v2->tex_coords[V];
                    draw_RGBApoint_struct( p );
            }
        }
    }
}

void draw_with_GL(POINT p0, POINT p1, POINT p2) {
    glBegin(GL_TRIANGLES);
    glColor3f(1, 1, 1);
    glVertex2f(p0.pos[X] + 400, p0.pos[Y] + 400);
    glVertex2f(p1.pos[X] + 400, p1.pos[Y] + 400);
    glVertex2f(p2.pos[X] + 400, p2.pos[Y] + 400);
    glEnd();
}

void calc_vector_angle(float normal_vector[4], float light_vector[4], float *angle) {
    float dot = dotProduct(normal_vector, light_vector); // gets the dot product of the normal vector and the light vector
    float mag_normal = magnitude(normal_vector); 
    float mag_light = magnitude(light_vector);
    *angle = acos(dot / (mag_normal * mag_light));
}


/*************************************************************************/
/* 3D object movement functions                                          */
/*************************************************************************/
void rotate_model(float angle_x, float angle_y, float angle_z) {
    float theta_x = angle_x / 360 * 2 * PI;
    float theta_y = angle_y / 360 * 2 * PI;
    float theta_z = angle_z / 360 * 2 * PI;

    for (int i = 0; i < num_verts; i++) {
        float x = vertex_list[i].world[X] - model_center_x;
        float y = vertex_list[i].world[Y] - model_center_y;
        float z = vertex_list[i].world[Z] - model_center_z;

        // Rotate around X-axis:
        float y1 = y * cos(theta_x) - z * sin(theta_x);
        float z1 = y * sin(theta_x) + z * cos(theta_x);

        // Rotate around Y-axis:
        float x2 = x * cos(theta_y) + z1 * sin(theta_y);
        float z2 = -x * sin(theta_y) + z1 * cos(theta_y);

        // Rotate around Z-axis:
        float x3 = x2 * cos(theta_z) - y1 * sin(theta_z);
        float y3 = x2 * sin(theta_z) + y1 * cos(theta_z);

        vertex_list[i].world[X] = x3 + model_center_x;
        vertex_list[i].world[Y] = y3 + model_center_y;
        vertex_list[i].world[Z] = z2 + model_center_z;
    }
}

void translate_model(float x, float y, float z) {
    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].world[X] += x;
        vertex_list[i].world[Y] += y;
        vertex_list[i].world[Z] += z;
    }
        model_center_x += x;
        model_center_y += y;
        model_center_z += z;
}

void scale_model(int scale_x, int scale_y, int scale_z) {
    float original_center_x = model_center_x;
    float original_center_y = model_center_y;
    float original_center_z = model_center_z;

    translate_model(-original_center_x, -original_center_y, -original_center_z); // translates the model to the origin

    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].world[X] *= scale_x;
        vertex_list[i].world[Y] *= scale_y;
        vertex_list[i].world[Z] *= scale_z;
    }

    translate_model(-original_center_x, -original_center_y, -original_center_z); // translate le model back to its original position
}

void scale_viewport(int scale_x, int scale_y, int scale_z) {
    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].pos[X] *= scale_x;
        vertex_list[i].pos[Y] *= scale_y;
        vertex_list[i].pos[Z] *= scale_z;
    }
}


/*************************************************************************/
/* 3D -> 2D functions                                                    */
/*************************************************************************/

void project_model(void) {  
    for (int i = 0; i < num_verts; i++) {
        copy_vector(vertex_list[i].world, vertex_list[i].pos); 
    }
}

void perspectiveproject_model(float  nearplane, float farplane) { // nearplane is the hypothetical screen
    for (int i = 0; i < num_verts; i++) {
        float x = vertex_list[i].world[X];
        float y = vertex_list[i].world[Y];
        float z = vertex_list[i].world[Z];

        vertex_list[i].pos[X] = (nearplane * x) / z;
        vertex_list[i].pos[Y] = (nearplane * y) / z;
        vertex_list[i].pos[Z] = z / (farplane - nearplane);
    }
}

void wireframe_model(void) { // draws a wireframe model
    for (int t = 0; t < num_triangles; t++) {
        int i0 = triangle_list[t].index_list[0];
        int i1 = triangle_list[t].index_list[1];
        int i2 = triangle_list[t].index_list[2];
        POINT p0 = vertex_list[i0];
        POINT p1 = vertex_list[i1];
        POINT p2 = vertex_list[i2];
        set_vector(p0.color, 1, 1, 1, 1);
        set_vector(p1.color, 1, 1, 1, 1);
        set_vector(p2.color, 1, 1, 1, 1);
        draw_line_color(p0, p1);
        draw_line_color(p1, p2);
        draw_line_color(p2, p0);
    }
}

void draw_model(void) { // draws a filled model
    float angle;
    for (int t = 0; t < num_triangles; t++) {
        int i0 = triangle_list[t].index_list[0];
        int i1 = triangle_list[t].index_list[1];
        int i2 = triangle_list[t].index_list[2];
        POINT p0 = vertex_list[i0];
        POINT p1 = vertex_list[i1];
        POINT p2 = vertex_list[i2];
        // set_vector(p0.color, 1, 0, 0, 1);
        // set_vector(p1.color, 0, 1, 0, 1);
        // set_vector(p2.color, 0, 0, 1, 1);
        calc_vector_angle(triangle_list[t].normal_vector, light_vector, &angle); // stores in angle as radian
        // closer to 0 the more aligned the light is with the normal vector
        if (angle > (PI / 2.0) && angle < (PI * 3.0 / 2.0)) {
            set_vector(p0.color, 0, 0, 0, 1);
            set_vector(p1.color, 0, 0, 0, 1);
            set_vector(p2.color, 0, 0, 0, 1);
        }
        else {
            float ratio = ((PI / 2.0) - angle) / (PI / 2.0); 
            set_vector(p0.color, ratio * p0.color[R], ratio * p0.color[G], ratio * p0.color[B], 1);
            set_vector(p1.color, ratio * p1.color[R], ratio * p1.color[G], ratio * p1.color[B], 1);
            set_vector(p2.color, ratio * p2.color[R], ratio * p2.color[G], ratio * p2.color[B], 1);
        }
        draw_barycentrictriangle(p0, p1, p2);
    }
}

/*************************************************************************/
/* Shape functions                                                       */
/*************************************************************************/
void init_sphere(int scale, float center_x, float center_y, float center_z) { // creates a sphere with radius 1
    num_verts = 0;
    for (int col = 0; col < MESH_WIDTH; col++) { // sets points
        for (int row = 0; row < MESH_HEIGHT; row++) {
            float u, v, x, y, z;
            u = (float)col / MESH_HEIGHT - 1;
            v = (float)row / MESH_WIDTH - 1;
            x = cos(u * PI) + center_x; // (0, 1) to (pi, -1) this is treated as the Z
            y = cos(v * 2 * PI) * sin(u * PI) + center_y; // this is treated as the x
            z = sin(v * 2 * PI) * sin(u * PI) + center_z; // this is treated as the y
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 1, 1, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    model_center_x = center_x;
    model_center_y = center_y;
    model_center_z = center_z;
    set_triangles();
}

void init_quad() {
    POINT p0, p1, p2, p3;
    set_vector(p0.world, -1, -1, 0, 1);
    set_vector(p0.tex_coords, 0, 0, 0, 0);
    set_vector(p1.world, 1, -1, 0, 1);
    set_vector(p1.tex_coords, 1, 0, 0, 0);
    set_vector(p2.world, 1, 1, 0, 1);
    set_vector(p2.tex_coords, 1, 1, 0, 0);
    set_vector(p3.world, -1, 1, 0, 1);
    set_vector(p3.tex_coords, 0, 1, 0, 0);
    vertex_list[0] = p0;
    vertex_list[1] = p1;
    vertex_list[2] = p2;
    vertex_list[3] = p3;
    num_verts = 4;
    triangle_list[0].index_list[0] = 0;
    triangle_list[0].index_list[1] = 3;
    triangle_list[0].index_list[2] = 1;
    triangle_list[1].index_list[0] = 3;
    triangle_list[1].index_list[1] = 2;
    triangle_list[1].index_list[2] = 1;
    num_triangles = 2;
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
    start_timer( &frame_timer );
    glClear(GL_COLOR_BUFFER_BIT);
    clear_colorbuffer();
    clear_depthbuffer();
    // draw_checkerboard(&Textures[0]);
    read_ppm(&Textures[0], "../Images/earth.ascii.ppm");

    if (test == 0) 
    { 
        start_timer( &vertex_timer );
        init_sphere(1, 0, 0, 1);
        scale_model(150, 150, 150);
        rotate_model(0, y_angle, 0);
        translate_model(0, 0, 0);
        project_model();
        scale_viewport(1, 1, 1);
        stop_timer( &vertex_timer );
        start_timer( &pixel_timer );
        draw_model();
        stop_timer( &pixel_timer );
    }
    else if (test == 1) 
    {
        init_quad();
        scale_model(1, 1, 1);
        rotate_model(0, y_angle, 0);
        translate_model(0, 0, 1);
        perspectiveproject_model(4, 50);
        scale_viewport(50, 50, 50);
        draw_model();
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
    stop_timer( &frame_timer );
    frame_time = elapsed_time(&frame_timer);
    vertex_time = elapsed_time(&vertex_timer);
    pixel_time = elapsed_time(&pixel_timer);
    printf( "vertex = %5.4f, pixel = %5.4f, frame = %5.4f, fps = %5.4f\n", vertex_time, pixel_time, frame_time, 1.0/frame_time );
    printf( "num_vertices = %5.2fk, num_triangles = %5.2fk, vertices/sec = %5.2fk, triangles/sec = %5.2fk\n", num_verts/1000.0, num_triangles/1000.0, num_verts/vertex_time/1000.0, num_triangles/pixel_time/1000.0 );

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
        y_angle += 1;
        if (y_angle > 360) y_angle = 0;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 'Y':
        y_angle -= 1;
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
 