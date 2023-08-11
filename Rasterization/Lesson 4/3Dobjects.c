/*
 *
 * lines.c - simple GLUT app that draws one frame with a single point at origin
 *
 * To build:  gcc -framework OpenGL -framework GLUT lines.c -o lines
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
    unsigned char data[4096][4096][4]; // just for making code easy and flexible. char is 0-255
} IMAGE;

typedef struct triangle
{
    int index_list[3];
} TRIANGLE;


/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int half_window_size = 400;
int Mojave_WorkAround = 1;
int draw_one_frame = 1;
int test = 0;
int textureIndex = 0;
int current[800][800];
IMAGE Textures[20];
int num_verts = 0;
POINT vertex_list[1000000]; // vertices
TRIANGLE triangle_list[1000000];
int num_triangles;
float y_angle = 0;
float z_angle = 0;
float z_translate = 0;

/*************************************************************************/
/* Texture functions                                                     */
/*************************************************************************/
void texture_sample(IMAGE *tex, float s, float t, float tex_color[4])
{                           // passing address is good for memory management. 4 bytes vs 16 MB
    int u = s * tex->width; // -> = (*tex).
    int v = t * tex->height;
    u = CLAMP(u, 0, tex->width - 1);
    v = CLAMP(v, 0, tex->height - 1);
    tex_color[R] = tex->data[v][u][R] / 255.0; // Since R is a unigned char
    tex_color[G] = tex->data[v][u][G] / 255.0; // Since R is a unigned char
    tex_color[B] = tex->data[v][u][B] / 255.0; // Since R is a unigned char
    tex_color[A] = tex->data[v][u][A] / 255.0; // Since R is a unigned char
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

/*************************************************************************/
/* utility functions                                                     */
/*************************************************************************/
/*
 * read_ppm()
 */
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

/*
 * random_float()
 */
float random_float(float low, float high)
{
    return ((float)(low + (rand() / (float)RAND_MAX) * (high - low)));
}
/*
 * draw_point_struct()
 */
void draw_point_struct(POINT p)
{
    /*
     * draw points
     */
    glColor4f(1, 1, 1, 1);
    // glColor4f(p.color[R], p.color[G], p.color[B], 1.0); // this gets overridden by the texture_sample
    // float tex_color[4];
    // texture_sample(&Textures[textureIndex], p.tex_coords[0],p.tex_coords[1], tex_color);
    // glColor4f(tex_color[R], tex_color[G], tex_color[B], 1.0); // effectively only pulls colours from a given IMAGE
    glBegin(GL_POINTS);
    glVertex2f(p.pos[X], p.pos[Y]);
    // printf("p.pos[X]: %f\n", p.pos[X]);
    // printf("p.pos[Y]: %f\n", p.pos[Y]);
    glEnd();
}


/*
 * draw_line_color()
 */
void draw_line_color( POINT p0, POINT p1 )
{
    float     pos_delta[4];
    float     color_delta[4];
    float     tex_delta[4];
    float     pos_inc[4];
    float     color_inc[4];
    float     tex_inc[4];
    POINT    p;
    
    if( (int)(p1.pos[X] - p0.pos[X]) == 0 && (int)(p1.pos[Y] - p0.pos[Y])== 0 )
    {
        /*
         * degenerate case - just draw point
         */
        draw_point_struct( p0 );
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
        subtract_vector( p1.tex_coords,   p0.tex_coords,      tex_delta   );

        float inv_slope = 1.0/pos_delta[X];

        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );
        scalar_multiply_vector( inv_slope, tex_delta,   tex_inc   );

        for( int x = p0.pos[X]; x <= p1.pos[X]; x++ )
        {
            draw_point_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
            add_vector( p.tex_coords,   tex_inc,   p.tex_coords   );
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
        subtract_vector( p1.tex_coords,   p0.tex_coords,      tex_delta   );

        float inv_slope = 1.0/pos_delta[Y];
        
        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );
        scalar_multiply_vector( inv_slope, tex_delta,   tex_inc   );

        for( int y = p0.pos[Y]; y <= p1.pos[Y]; y++ )
        {
            draw_point_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
            add_vector( p.tex_coords,   tex_inc,   p.tex_coords   );
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
        subtract_vector( p1.tex_coords,   p0.tex_coords,      tex_delta   );

        float inv_slope = 1.0/pos_delta[X];

        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );
        scalar_multiply_vector( inv_slope, tex_delta,   tex_inc   );

        for( int x = p0.pos[X]; x <= p1.pos[X]; x++ )
        {
            draw_point_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
            add_vector( p.tex_coords,   tex_inc,   p.tex_coords   );
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
        subtract_vector( p1.tex_coords,   p0.tex_coords,      tex_delta   );

        float inv_slope = 1.0/pos_delta[Y];

        scalar_multiply_vector( inv_slope, pos_delta,   pos_inc   );
        scalar_multiply_vector( inv_slope, color_delta, color_inc );
        scalar_multiply_vector( inv_slope, tex_delta,   tex_inc   );

        for( int y = p0.pos[Y]; y <= p1.pos[Y]; y++ )
        {
            draw_point_struct( p );

            add_vector( p.pos,   pos_inc,   p.pos   );
            add_vector( p.color, color_inc, p.color );
            add_vector( p.tex_coords,   tex_inc,   p.tex_coords   );
        }
    }
}

void init_cube(void) {
    num_triangles = 12;
    num_verts = 8;
    set_vector(vertex_list[0].world, -0.5, 0.5, -0.5, 1);
    set_vector(vertex_list[1].world, 0.5, 0.5, -0.5, 1);
    set_vector(vertex_list[2].world, 0.5, -0.5, -0.5, 1);
    set_vector(vertex_list[3].world, -0.5, -0.5, -0.5, 1);
    set_vector(vertex_list[4].world, -0.5, 0.5, 0.5, 1);
    set_vector(vertex_list[5].world, 0.5, 0.5, 0.5, 1);
    set_vector(vertex_list[6].world, 0.5, -0.5, 0.5, 1);
    set_vector(vertex_list[7].world, -0.5, -0.5, 0.5, 1);
    triangle_list[0].index_list[0] = 0; // front face top triangle
    triangle_list[0].index_list[1] = 1;
    triangle_list[0].index_list[2] = 2;
    triangle_list[1].index_list[0] = 0; // front face bottom triangle
    triangle_list[1].index_list[1] = 2;
    triangle_list[1].index_list[2] = 3;
    triangle_list[2].index_list[0] = 4; // left face top triangle
    triangle_list[2].index_list[1] = 0;
    triangle_list[2].index_list[2] = 3;
    triangle_list[3].index_list[0] = 4; // left face bottom triangle
    triangle_list[3].index_list[1] = 3;
    triangle_list[3].index_list[2] = 7;
    triangle_list[4].index_list[0] = 1; // right face top triangle
    triangle_list[4].index_list[1] = 5;
    triangle_list[4].index_list[2] = 6;
    triangle_list[5].index_list[0] = 1; // right face bottom triangle
    triangle_list[5].index_list[1] = 6;
    triangle_list[5].index_list[2] = 2;
    triangle_list[6].index_list[0] = 5; // back face top triangle 
    triangle_list[6].index_list[1] = 4;
    triangle_list[6].index_list[2] = 7;
    triangle_list[7].index_list[0] = 5; // back face bottom triangle
    triangle_list[7].index_list[1] = 7;
    triangle_list[7].index_list[2] = 6;
    triangle_list[8].index_list[0] = 4; // top face top triangle
    triangle_list[8].index_list[1] = 5;
    triangle_list[8].index_list[2] = 1;
    triangle_list[9].index_list[0] = 4; // top face bottom triangle
    triangle_list[9].index_list[1] = 1;
    triangle_list[9].index_list[2] = 0;
    triangle_list[10].index_list[0] = 3; // bottom face top triangle
    triangle_list[10].index_list[1] = 6;
    triangle_list[10].index_list[2] = 2;
    triangle_list[11].index_list[0] = 7; // bottom face bottom triangle
    triangle_list[11].index_list[1] = 2;
    triangle_list[11].index_list[2] = 3;
}

void draw_triangles() {
    num_triangles = 0;
    for (int col = 0; col < 31; col++) { // sets triangles
        for (int row = 0; row < 31; row++) {
            triangle_list[num_triangles].index_list[0] = (row * 32) + col; // flat top triangle
            triangle_list[num_triangles].index_list[1] = (row * 32) + col + 1;
            triangle_list[num_triangles].index_list[2] = ((row + 1) * 32) + col + 1;
            triangle_list[num_triangles].index_list[0] = (row * 32) + col; // flat bot triangle
            triangle_list[num_triangles].index_list[1] = ((row + 1) * 32) + col + 1;
            triangle_list[num_triangles].index_list[2] = ((row + 1) * 32) + col;
            num_triangles += 2;
        }
    }
}


void init_mesh(int scale, float center_x, float center_y, float center_z) { // flat 2d plane with points and triangles
    num_verts = 0;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            x = center_x + (u - 0.5) * scale;
            y = center_y + (v - 0.5) * scale;
            z = center_z;
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

void init_cylinder(int scale, float center_x, float center_y, float center_z, float radius) { // wrap the 2d plane into the shape of a cylinder, ignoring the top and bottom
    num_verts = 0;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            x = radius * cos(u * 2 * PI) + center_x; // controls the x value based on the row
            y = radius * sin(u * 2 * PI) + center_y; // controls the y value based on the column
            z = scale * v + center_z; // controls the height of the cylinder since it goes into the z axis
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

void init_cone(int scale, float center_x, float center_y, float center_z, float radius) {
    num_verts = 0;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            float r = radius * (1.0 - v); // as row gets larger, the radius gets smaller and smaller. This assumes that we draw the cone with the point in the back
            x = r * cos(u * 2 * PI) + center_x;
            y = r * sin(u * 2 * PI) + center_y;
            z = scale * v + center_z;
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

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

void init_torus(int scale, float center_x, float center_y, float center_z, float tube_radius, float hole_radius) {
    num_verts = 0;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            x = (tube_radius * cos(v * 2 * PI)) + (hole_radius * cos(u * 2 * PI)) + center_x; // play with radiuses
            y = (tube_radius * sin(v * 2 * PI)) + (hole_radius * sin(u * 2 * PI)) + center_y;
            z = (tube_radius * sin(v * 2 * PI) + center_z);
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

void init_helix(int scale, float center_x, float center_y, float center_z, float tube_radius, float hole_radius) {
    num_verts = 0;
    int num_of_turns = 3; // changes how many rotations the helix makes
    for (int col = 0; col < 32; col++) { 
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31; 
            x = (tube_radius * cos(v * 2 * PI * num_of_turns)) + (hole_radius * cos(u * 2 * PI * num_of_turns)) + center_x;
            y = (tube_radius * sin(v * 2 * PI * num_of_turns)) + (hole_radius * sin(u * 2 * PI * num_of_turns)) + center_y;
            z = tube_radius * sin(v * 2 * PI) + center_z + (4 * u); // Makes the helix height longer
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

void init_horn(int scale, float center_x, float center_y, float center_z, float radius) {
    num_verts = 0;
    int curve = 2;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            float r = radius * pow((1.0 - v), curve); // as row gets larger, the radius gets smaller and smaller. This assumes that we draw the cone with the point in the back
            x = r * cos(u * 2 * PI) + center_x;
            y = r * sin(u * 2 * PI) + center_y;
            z = scale * v + center_z;
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}

void init_parabola(int scale, float center_x, float center_y, float center_z, float tube_radius) {
    num_verts = 0;
    for (int col = 0; col < 32; col++) { // sets points
        for (int row = 0; row < 32; row++) {
            float u, v, x, y, z;
            u = (float)col / 31;
            v = (float)row / 31;
            x = (tube_radius * cos(v * 2 * PI)) + (u) + center_x;
            y = (tube_radius * sin(v * 2 * PI)) + (u * u) + center_y;
            z = (tube_radius * sin(v * 2 * PI) + center_z);
            set_vector(vertex_list[num_verts].world, x, y, z, 1);
            set_vector(vertex_list[num_verts].color, 1, 0, 0, 1);
            set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
            num_verts++;
        }
    }
    draw_triangles();
}



void translate_model(float dx, float dy, float dz) { // bugged since because you are translating in the z axis only, sometimes when the z axis becomes the y axis, you aren't simply moving in and out
    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].world[X] += dx;
        vertex_list[i].world[Y] += dy;
        vertex_list[i].world[Z] += dz;
    }
}

void draw_model(void) {
    for (int t = 0; t < num_triangles; t++) {
        int i0 = triangle_list[t].index_list[0];
        int i1 = triangle_list[t].index_list[1];
        int i2 = triangle_list[t].index_list[2];
        POINT p0 = vertex_list[i0];
        POINT p1 = vertex_list[i1];
        POINT p2 = vertex_list[i2];
        glColor4f(1.0, 1.0, 1.0, 1.0);
        draw_line_color(p0, p1);
        draw_line_color(p1, p2);
        draw_line_color(p2, p0);
    }
}

void scale_model(void) {
    for (int i = 0; i < num_verts; i++) {
        vertex_list[i].world[X] *= 100;
        vertex_list[i].world[Y] *= 100;
        vertex_list[i].world[Z] *= 100;
    }
}

void project_model(void) {
    for (int i = 0; i < num_verts; i++) {
        copy_vector(vertex_list[i].world, vertex_list[i].pos); 
    }
}

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

void perspective_model(void) {
    int plane_distance = 50; // hard-coded for now. This is the distance from the eyeball to the projection plane
    float eyeball[4], d[4];
    float dx, dy;
    set_vector(eyeball, 0, 0, -450, 1);
        for (int i = 0; i < num_verts; i++) {
            subtract_vector(vertex_list[i].world, eyeball, d);
            dx = d[X] / d[Z]; // include a check state for 0
            dy = d[Y] / d[Z];
            vertex_list[i].world[X] += dx * plane_distance;
            vertex_list[i].world[Y] += dy * plane_distance;
        }

    }

void set_vertex(float x, float y, float z, float r, float g, float b, float a, float u, float v) {
    set_vector(vertex_list[num_verts].world, x, y, z, 1);
    set_vector(vertex_list[num_verts].color, r, g, b, a);
    set_vector(vertex_list[num_verts].tex_coords, u, v, 0, 0);
    num_verts++;
}

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
        num_triangles++;
    }
    fclose(fp);
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

    if (test == 0)
    {
        init_cube();
        scale_model();
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 1)
    {
        init_mesh(5, 0, 0, 0);
        scale_model();
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 2)
    {
        init_cylinder(1, 0, 0, 0, 1);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 3)
    {
        init_cone(1, 0, 0, 0, 1);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 4)
    {
        init_sphere(1, 0, 0, 0, 1);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 5)
    {
        init_torus(1, 0, 0, 0, 0.25, 1.5);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 6)
    {
        init_helix(1, 0, 0, 0, 0.25, 1.5);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 7)
    {
        init_horn(1, 0, 0, 0, 1);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 8)
    {
        init_parabola(1, 0, 0, 0, 0.25);
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
    }
    else if (test == 9)
    {
        read_obj_file("teapot.obj");
        scale_model();
        translate_model(0, 0, z_translate);
        rotate_model(0, y_angle, 0);
        perspective_model();
        project_model();
        draw_model();
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
    case 'z':
        z_angle++;
        if (z_angle > 360) z_angle = 0;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 'Z':
        z_angle--;
        if (z_angle < 0) z_angle = 360;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case '+':
        z_translate++;;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case '-':
        z_translate--;
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
    gluOrtho2D(-half_window_size, half_window_size, -half_window_size, half_window_size);
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
