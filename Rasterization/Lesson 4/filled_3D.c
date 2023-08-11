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
    float normal_vector[4]; // normal vector for the triangle 
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
float light_vector[4] = {1, 1, 1, 0}; // sets the infinite light vector to 1, 1, 1
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
 * init_normal() // calculates and stores the normal vector for a given triangle into the normal coordinates.
 */
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

/*
 * calculate angle between normal and light vector
 */
void calc_vector_angle(float normal_vector[4], float light_vector[4], float *angle) {
    float dot = dotProduct(normal_vector, light_vector); // gets the dot product of the normal vector and the light vector
    float mag_normal = magnitude(normal_vector); 
    float mag_light = magnitude(light_vector);
    *angle = acos(dot / (mag_normal * mag_light));
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
    int x, y;
    x = p.pos[X] + 400;
    y = p.pos[Y] + 400;
    float angle;
    if (x >= 0 && x < 800 && y >= 0 && y < 800) { // checks for out of bounds values
        if (p.world[Z] < depth[y][x]) { // checks if the point is closer than the previous point
            depth[y][x] = p.world[Z];
            set_vector(colorbuffer[y][x], p.color[R], p.color[G], p.color[B], p.color[A]); // sets the color buffer to the color of the point
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

void draw_triangles() {
    num_triangles = 0;
    for (int col = 0; col < 31; col++) { // sets triangles
        for (int row = 0; row < 31; row++) {
            triangle_list[num_triangles].index_list[0] = (row * 32) + col; // flat top triangle
            triangle_list[num_triangles].index_list[1] = (row * 32) + col + 1;
            triangle_list[num_triangles].index_list[2] = ((row + 1) * 32) + col + 1;
            triangle_list[num_triangles] = calc_normal(triangle_list[num_triangles]); // sets the normal vector for a given triangle
            triangle_list[num_triangles+1].index_list[0] = (row * 32) + col; // flat bot triangle
            triangle_list[num_triangles+1].index_list[1] = ((row + 1) * 32) + col + 1;
            triangle_list[num_triangles+1].index_list[2] = ((row + 1) * 32) + col;
            triangle_list[num_triangles+1] = calc_normal(triangle_list[num_triangles+1]); // sets the normal vector for the second triangle
            num_triangles += 2;
        }
    }
}

// initializes the points and triangles of a sphere
void init_sphere(int scale, float center_x, float center_y, float center_z, float radius) {
    num_verts = 0;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            x = radius * cos(u * PI) + center_x; // (0, 1) to (pi, -1) this is treated as the Z
            y = radius * cos(v * 2 * PI) * sin(u * PI) + center_y; // this is treated as the x
            z = radius * sin(v * 2 * PI) * sin(u * PI) + center_z; // this is treated as the y
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

// Reads and initalizes a 3D object from an .obj file
void read_obj_file(char *filename) {
    FILE *fp; // file pointer but doesn't need to be called *fp
    num_verts = 0;
    num_triangles = 0;
    float x, y, z;
    int i0, i1, i2;
    fp = fopen(filename, "r"); // opens the file
    if (fp == NULL) return; // checks if the file is opened successfully
    while ( fscanf(fp, "v %f %f %f\n", &x, &y, &z) == 3) { // looks for stuff in that format, and saves it into the variables x, y, and z
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, 1, 1, 0, 0);
            num_verts++;
    }
    fscanf(fp, " \n");
    while ( fscanf(fp, "f %d %d %d\n", &i0, &i1, &i2) == 3) {
        triangle_list[num_triangles].index_list[0] = i0 - 1; // index starts @ 1
        triangle_list[num_triangles].index_list[1] = i1 - 1;
        triangle_list[num_triangles].index_list[2] = i2 - 1;
        calc_normal(triangle_list[num_triangles]); // sets the normal vector for a given triangle
        num_triangles++;
    }
    fclose(fp);
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
                p.world[Z] = alpha * p0.world[Z] + beta * p1.world[Z] + gamma * p2.world[Z]; // sets the depth at that point to the weighted average of the depths of the vertices 
                draw_RGBApoint_struct(p);
            }
        }
    }
}  

// Rotates the model
void rotate_model(float angle_x, float angle_y, float angle_z) {
    for (int i = 0; i < num_verts; i++) {
        float x = vertex_list[i].world[X];
        float y = vertex_list[i].world[Y];
        float z = vertex_list[i].world[Z];
        float theta = angle_y / 360 * 2 * PI;
        vertex_list[i].world[X] = x * cos(theta) - z * sin(theta);
        vertex_list[i].world[Y] = y;
        vertex_list[i].world[Z] = x * sin(theta) + z * cos(theta);
    }
}

// Scales the model by 100
void scale_model(void) {
    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].world[X] *= 100;
        vertex_list[i].world[Y] *= 100;
        vertex_list[i].world[Z] *= 100;
    }
}

// projects the 3D model onto the 2D screen
void project_model(void) {  
    for (int i = 0; i < num_verts; i++) {
        copy_vector(vertex_list[i].world, vertex_list[i].pos); 
    }
}

// Draws the 3D model
void draw_model(void) {
    for (int t = 0; t < num_triangles; t++) {
        int i0 = triangle_list[t].index_list[0];
        int i1 = triangle_list[t].index_list[1];
        int i2 = triangle_list[t].index_list[2];
        POINT p0 = vertex_list[i0];
        POINT p1 = vertex_list[i1];
        POINT p2 = vertex_list[i2];
        float greyscale = 1 / t;
        set_vector(p0.color, 1, 1, 1, 1);
        set_vector(p1.color, greyscale, greyscale, greyscale, 1);
        set_vector(p2.color, greyscale, greyscale, greyscale, 1);
        barycentric_triangle(p0, p1, p2);
    }
}

// sets the perspective of the model. Before this, everything is using world coordinates. After this, everything is using perspective coordinates
void perspective_model(void) {
    int plane_distance = 50; // hard-coded for now. This is the distance from the eyeball to the projection plane
    int far_plane = 500;
    float eyeball[4], d[4];
    float dx, dy;
    set_vector(eyeball, 0, 0, -450, 1);
        for (int i = 0; i < num_verts; i++) {
            // subtract_vector(vertex_list[i].world, eyeball, d);
            copy_vector(vertex_list[i].world, d);
            dx = d[X] / d[Z]; // include a check state for 0
            dy = d[Y] / d[Z];
            vertex_list[i].pos[X] += dx * plane_distance;
            vertex_list[i].pos[Y] += dy * plane_distance;
            vertex_list[i].pos[Z] = d[Z];
            // vertex_list[i].pos[Z] = (d[Z] - plane_distance) / (far_plane - plane_distance); 
        }

    }


void translate_model(float tx, float ty, float tz) {
    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].world[X] += tx;
        vertex_list[i].world[Y] += ty;
        vertex_list[i].world[Z] += tz;
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

    if (test == 0) // confirms that draw_line is working
    {
        POINT p0, p1;
        set_vector(p0.pos, -100, -100, 0, 1);
        set_vector(p0.color, 1, 0, 0, 1);
        set_vector(p1.pos, 100, 100, 0, 1);
        set_vector(p1.color, 0, 1, 0, 1);
        draw_line_color(p0, p1);
    }
    else if (test == 1) // confirms that barycentric triangle drawing is working
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
    else if (test == 2) // overlapping test
    {
        overlapping_triangleTest();
    }
    else if (test == 3) // tests sphere drawing is working
    {
        init_sphere(1, 0, 0, 0, 100);
        rotate_model(0, y_angle, 0);
        translate_model(0, 0, 5);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 4)
    {
        read_obj_file("teapot.obj");
        scale_model();
        rotate_model(0, y_angle, 0);
        translate_model(0, 0, 5);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 5)
    {
        POINT p;
        set_vector(p.pos, 0, 0, 0, 1);
        set_vector(p.color, 1, 1, 1, 1);
        draw_RGBApoint_struct(p);
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