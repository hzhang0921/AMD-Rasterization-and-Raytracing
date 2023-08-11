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
} POINT;

typedef struct image
{
    int width;
    int height;
    unsigned char data[4096][4096][4]; // just for making code easy and flexible. char is 0-255
} IMAGE;

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
    // glColor4f(p.color[R], p.color[G], p.color[B], 1.0); // this gets overridden by the texture_sample
    float tex_color[4];
    texture_sample(&Textures[textureIndex], p.tex_coords[0],p.tex_coords[1], tex_color);
    glColor4f(tex_color[R], tex_color[G], tex_color[B], 1.0); // effectively only pulls colours from a given IMAGE
    glBegin(GL_POINTS);
    glVertex2f(p.pos[X], p.pos[Y]); 
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

/*
 * draw textured quad()
 */
void draw_quad(POINT tl, POINT tr, POINT bl, POINT br)
{
    POINT p0 = bl;
    POINT p1 = br;
    float left_delta[4], left_inc[4], right_delta[4], right_inc[4];
    float left_texD[4], left_texInc[4], right_texD[4], right_texInc[4];

    subtract_vector(tl.color, bl.color, left_delta);
    subtract_vector(tr.color, br.color, right_delta);
    subtract_vector(tl.tex_coords, bl.tex_coords, left_texD);
    subtract_vector(tr.tex_coords, br.tex_coords, right_texD);

    int y0 = bl.pos[Y];
    int y1 = tl.pos[Y];
    float dy = y1 - y0;

    scalar_multiply_vector(1.0 / dy, left_delta, left_inc);
    scalar_multiply_vector(1.0 / dy, right_delta, right_inc);
    scalar_multiply_vector(1.0 / dy, left_texD, left_texInc);
    scalar_multiply_vector(1.0 / dy, right_texD, right_texInc);

    for (int y = y0; y <= y1; y++)
    {
        draw_line_color(p0, p1);
        p0.pos[Y] += 1;
        p1.pos[Y] += 1; 

        add_vector(p0.color, left_inc, p0.color);
        add_vector(p1.color, right_inc, p1.color);
        add_vector(p0.tex_coords, left_texInc, p0.tex_coords);
        add_vector(p1.tex_coords, right_texInc, p1.tex_coords);
    }
}


/*
 * draws and saves a checkerbaord image at the given IMAGE pointer
 */
void draw_checkerboard(IMAGE *image) {  // *image is where the IMAGE is saved
    image->width = 800;
    image->height = 800;
    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            int size = 8;
            if ((x / size) % 2 == 0) {
                if ((y / size) % 2 == 0) {
                    UCset_vector(image->data[y][x], 255, 255, 255, 255);
                }
                else {
                    UCset_vector(image->data[y][x], 0, 0, 0, 255);
                }
            }
            else {
                if ((y / size) % 2 == 0) {
                    UCset_vector(image->data[y][x], 0, 0, 0, 255);
                }
                else {
                    UCset_vector(image->data[y][x], 255, 255, 255, 255); 
                }
            }
        }
    }
    // write_ppm(&Textures[1], "../Images/checkerboard", 0); // saves the IMAGE to a given file location with file name
}


/*
 * draw_gradient() Draws a gradient from left to right. If the lefet and right side are the same color, it will draw a solid color
 */
void draw_gradient(IMAGE *image, float r1, float g1, float b1, float r2, float g2, float b2) {
    image->width = 800;
    image->height = 800;
    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            unsigned char r = (int)((r2 - r1) * (x / (float)image->width) * 255);
            unsigned char g = (int)((g2 - g1) * (x / (float)image->width) * 255);
            unsigned char b = (int)((b2 - b1) * (x / (float)image->width) * 255);
            image->data[y][x][R] = r ;
            image->data[y][x][G] = g ;
            image->data[y][x][B] = b ;
            image->data[y][x][A] = 255;
        }
    }
    // write_ppm(&Textures[2], "../Images/gradient", 0); // saves the IMAGE to a given file location with file name
}

/*
 * copy Image() // Copies from image1 to image2
 */
void copyImage(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 0; y < image2->height; y++) {
        for (int x = 0; x < image2->width; x++) {
            image2->data[y][x][R] = image1->data[y][x][R];
            image2->data[y][x][G] = image1->data[y][x][G];
            image2->data[y][x][B] = image1->data[y][x][B];
            image2->data[y][x][A] = image1->data[y][x][A];
        }
    }
    // write_ppm(&Textures[3], "../Images/copyImage", 0); // saves the IMAGE to a given file location with file name
}

/*
 * negative_image() 
 */
void negativeImage(IMAGE *image1, IMAGE *image2) { // copies a negative image from image1 to image2
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 0; y < image2->height; y++) {
        for (int x = 0; x < image2->width; x++) {
            image2->data[y][x][R] = 255 - image1->data[y][x][R];
            image2->data[y][x][G] = 255 - image1->data[y][x][G];
            image2->data[y][x][B] = 255 - image1->data[y][x][B];
            image2->data[y][x][A] = 255 - image1->data[y][x][A];
        }
    }
    // write_ppm(&Textures[3], "../Images/NegativeImage", 0); // saves the IMAGE to a given file location with file name
}

/*
 * rotate()
 */
void rotate(IMAGE *image1, IMAGE *image2, int angle) {
    image2->height = image1->height;
    image2->width = image1->width;
    double centerX = image2->width / 2.0;
    double centerY = image2->height / 2.0;

    for (int y = 0; y < image2->height; y++) {
        for (int x = 0; x < image2->width; x++) {
            // Translate the point to the origin before rotation
            double x1 = x - centerX;
            double y1 = y - centerY;
            // Calculates the polar coordinates
            double r = sqrt((x1 * x1) + (y1 * y1));
            double theta = atan2(y1, x1); // atan2 returns the angle in radians
            // printf("theta: %f %f, %f \n", theta, x1, y1);
            // Apply rotation
            double rotated_x = r * cos(theta + ((angle) * PI / 180.0)); // Converts angle to radians and adds it onto the theta, then converts it back to cartesian. 
            double rotated_y = r * sin(theta + ((angle) * PI / 180.0));
            // Trying Rotation Matrix
            // int matrix_x = x1 * cos(theta) - y1 * sin(theta);
            // int matrix_y = x1 * sin(theta) + y1 * cos(theta);
            // Translate the point back after rotation
            int converted_x = (int)(rotated_x + centerX);
            int converted_y = (int)(rotated_y + centerY);
            // Apply clamping to ensure indices are within bounds
            if ( converted_x < 0 || converted_x >= image1->width || converted_y < 0 || converted_y >= image1->height ) {
                image2->data[y][x][R] = 0;
                image2->data[y][x][G] = 0;
                image2->data[y][x][B] = 0;
                image2->data[y][x][A] = 0;
            }
            else {
                image2->data[y][x][R] = image1->data[CLAMP(converted_y, 0, image2->height-1)][CLAMP(converted_x, 0, image2->width-1)][R];
                image2->data[y][x][G] = image1->data[CLAMP(converted_y, 0, image2->height-1)][CLAMP(converted_x, 0, image2->width-1)][G];
                image2->data[y][x][B] = image1->data[CLAMP(converted_y, 0, image2->height-1)][CLAMP(converted_x, 0, image2->width-1)][B];
                image2->data[y][x][A] = image1->data[CLAMP(converted_y, 0, image2->height-1)][CLAMP(converted_x, 0, image2->width-1)][A];
            }
        }
    }
}

/*
 * flipVertical()
 */
void flipvertical(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 0; y < image2->height; y++) {
        for (int x = 0; x < image2->width; x++) {
            image2->data[y][x][R] = image1->data[image2->height - y - 1][x][R];
            image2->data[y][x][G] = image1->data[image2->height - y - 1][x][G];
            image2->data[y][x][B] = image1->data[image2->height - y - 1][x][B];
            image2->data[y][x][A] = image1->data[image2->height - y - 1][x][A];
        }
    }
}

/*
 * luminosity() // calculates the average of red + green + blue aka grayscale
 */
void luminosity(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 0; y < image2->height; y++) {
        for (int x = 0; x < image2->width; x++) {
            int avg = (image1->data[y][x][R] + image1->data[y][x][G] + image1->data[y][x][B]) / 3;
            image2->data[y][x][R] = avg;
            image2->data[y][x][G] = avg;
            image2->data[y][x][B] = avg;
            image2->data[y][x][A] = 255;
        }
    }
}

/*
 * sepia tone() // adds a sepia filter to an image
 */
void sepia(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 0; y < image2->height; y++) {
        for (int x = 0; x < image2->width; x++) {
            image2->data[y][x][R] = CLAMP(((image1->data[y][x][R] * .393) + (image1->data[y][x][G] * .769) + (image1->data[y][x][B] * .189)), 0, 255);
            image2->data[y][x][G] = CLAMP(((image1->data[y][x][R] * .349) + (image1->data[y][x][G] * .686) + (image1->data[y][x][B] * .168)), 0, 255);
            image2->data[y][x][B] = CLAMP(((image1->data[y][x][R] * .272) + (image1->data[y][x][G] * .534) + (image1->data[y][x][B] * .131)), 0, 255);
            image2->data[y][x][A] = 255;
        }
    }
}

/*
 * blur() // blurs an image
 */
void blur(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 1; y < (image2->height)-1; y++) {
        for (int x = 1; x < (image2->width)-1; x++) {
            int avgred = ((image1->data[y][x][R] + image1->data[y][x-1][R] + image1->data[y][x+1][R] + image1->data[y-1][x][R] + image1->data[y+1][x][R] + image1->data[y-1][x-1][R] + image1->data[y-1][x+1][R] + image1->data[y+1][x-1][R] + image1->data[y+1][x+1][R]) / 9.0);
            int avggre = ((image1->data[y][x][G] + image1->data[y][x-1][G] + image1->data[y][x+1][G] + image1->data[y-1][x][G] + image1->data[y+1][x][G] + image1->data[y-1][x-1][G] + image1->data[y-1][x+1][G] + image1->data[y+1][x-1][G] + image1->data[y+1][x+1][G]) / 9.0);
            int avgblu = ((image1->data[y][x][B] + image1->data[y][x-1][B] + image1->data[y][x+1][B] + image1->data[y-1][x][B] + image1->data[y+1][x][B] + image1->data[y-1][x-1][B] + image1->data[y-1][x+1][B] + image1->data[y+1][x-1][B] + image1->data[y+1][x+1][B]) / 9.0);
            image2->data[y][x][R] = avgred;
            image2->data[y][x][G] = avggre;
            image2->data[y][x][B] = avgblu;
        }
    }
}

/*
 * min() // finds the minimum value of 9 given values
 */
unsigned char min_char(unsigned char a, unsigned char b, unsigned char c, unsigned char d, unsigned char e, unsigned char f, unsigned char g, unsigned char h, unsigned char i) {
    unsigned char min = a;
    if (b < min) {
        min = b;
    }
    if (c < min) {
        min = c;
    }
    if (d < min) {
        min = d;
    }
    if (e < min) {
        min = e;
    }
    if (f < min) {
        min = f;
    }
    if (g < min) {
        min = g;
    }
    if (h < min) {
        min = h;
    }
    if (i < min) {
        min = i;
    }
    return min;
}

/*
 * min_neighbor() // finds the minimum value of the 8 surrounding pixels
 */
 void min_neighbor(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 1; y < (image2->height)-1; y++) {
        for (int x = 1; x < (image2->width)-1; x++) {
            int minred = min_char(image1->data[y][x][R], image1->data[y][x-1][R], image1->data[y][x+1][R], image1->data[y-1][x][R], image1->data[y+1][x][R], image1->data[y-1][x-1][R], image1->data[y-1][x+1][R], image1->data[y+1][x-1][R], image1->data[y+1][x+1][R]);
            int mingre = min_char(image1->data[y][x][G], image1->data[y][x-1][G], image1->data[y][x+1][G], image1->data[y-1][x][G], image1->data[y+1][x][G], image1->data[y-1][x-1][G], image1->data[y-1][x+1][G], image1->data[y+1][x-1][G], image1->data[y+1][x+1][G]);
            int minblu = min_char(image1->data[y][x][B], image1->data[y][x-1][B], image1->data[y][x+1][B], image1->data[y-1][x][B], image1->data[y+1][x][B], image1->data[y-1][x-1][B], image1->data[y-1][x+1][B], image1->data[y+1][x-1][B], image1->data[y+1][x+1][B]);
            image2->data[y][x][R] = minred;
            image2->data[y][x][G] = mingre;
            image2->data[y][x][B] = minblu;
            // printf("%d %d %d\n", minred, mingre, minblu);
        }
    }
 }

/*
 * max_char() // finds the maximum value of 9 given values
 */
unsigned char max_char(unsigned char a, unsigned char b, unsigned char c, unsigned char d, unsigned char e, unsigned char f, unsigned char g, unsigned char h, unsigned char i) {
    unsigned char max = a;
    if (b > max) {
        max = b;
    }
    if (c > max) {
        max = c;
    }
    if (d > max) {
        max = d;
    }
    if (e > max) {
        max = e;
    }
    if (f > max) {
        max = f;
    }
    if (g > max) {
        max = g;
    }
    if (h > max) {
        max = h;
    }
    if (i > max) {
        max = i;
    }
    return max;
}

/*
 * max_neighbor() 
 */
void max_neighbor(IMAGE *image1, IMAGE *image2) {
    image2->width = image1->width;
    image2->height = image1->height;
    for (int y = 1; y < (image2->height)-1; y++) {
        for (int x = 1; x < (image2->width)-1; x++) {
            int maxred = max_char(image1->data[y][x][R], image1->data[y][x-1][R], image1->data[y][x+1][R], image1->data[y-1][x][R], image1->data[y+1][x][R], image1->data[y-1][x-1][R], image1->data[y-1][x+1][R], image1->data[y+1][x-1][R], image1->data[y+1][x+1][R]);
            int maxgre = max_char(image1->data[y][x][G], image1->data[y][x-1][G], image1->data[y][x+1][G], image1->data[y-1][x][G], image1->data[y+1][x][G], image1->data[y-1][x-1][G], image1->data[y-1][x+1][G], image1->data[y+1][x-1][G], image1->data[y+1][x+1][G]);
            int maxblu = max_char(image1->data[y][x][B], image1->data[y][x-1][B], image1->data[y][x+1][B], image1->data[y-1][x][B], image1->data[y+1][x][B], image1->data[y-1][x-1][B], image1->data[y-1][x+1][B], image1->data[y+1][x-1][B], image1->data[y+1][x+1][B]);
            image2->data[y][x][R] = maxred;
            image2->data[y][x][G] = maxgre;
            image2->data[y][x][B] = maxblu;
            // printf("%d %d %d\n", maxred, maxgre, maxblu);
        }
    }
}



/*
 * Draw function for outputting to OpenGL window
 */
void draw_GLquad() // draws the actual quad on GL. To control the sampling texture, change TextureIndex
{
    POINT p0, p1, p2, p3;
    set_vector(p0.pos, -300, 300, 0, 0);
    set_vector(p1.pos, 300, 300, 0, 0);
    set_vector(p2.pos, -300, -300, 0, 0);
    set_vector(p3.pos, 300, -300, 0, 0);
    // set_vector(p0.color, 1, 0, 0, 1);
    // set_vector(p1.color, 0, 0, 1, 1);
    // set_vector(p2.color, 0, 1, 0, 1);
    // set_vector(p3.color, 1, 0, 1, 1);
    set_vector(p0.tex_coords, 0, 0, 0, 0);
    set_vector(p1.tex_coords, 1, 0, 0, 0);
    set_vector(p2.tex_coords, 0, 1, 0, 0);
    set_vector(p3.tex_coords, 1, 1, 0, 0);
    draw_quad(p0, p1, p2, p3);
}

void fisheye(IMAGE *input, IMAGE *output, float angle) { // distorts the image to a fisheye view, and rotates it by the angle
    output->height = input->height;
    output->width = input->width;
    double centerX = output->width / 2.0;
    double centerY = input->height / 2.0;
    for (int y = 0; y < output->height; y++) {
        for (int x = 0; x < input->width; x++) {
            // Translate the point to the origin before rotation
            double x1 = x - centerX;
            double y1 = y - centerY;
            // Calculates the polar coordinates
            double r = sqrt((x1 * x1) + (y1 * y1));
            r = r * (r / (sqrt((centerX/2 * centerX/2) + (centerY/2 * centerY/2))));
            double theta = atan2(y1, x1); // atan2 returns the angle in radians
            // printf("theta: %f %f, %f \n", theta, x1, y1);
            // Apply rotation
            double rotated_x = r * cos(theta + ((angle) * PI / 180.0)); // Converts angle to radians and adds it onto the theta, then converts it back to cartesian. 
            double rotated_y = r * sin(theta + ((angle) * PI / 180.0));
            // Translate the point back after rotation
            int converted_x = (int)(rotated_x + centerX);
            int converted_y = (int)(rotated_y + centerY);
            // Apply clamping to ensure indices are within bounds
            if ( converted_x < 0 || converted_x >= input->width || converted_y < 0 || converted_y >= input->height ) {
                output->data[y][x][R] = 0;
                output->data[y][x][G] = 0;
                output->data[y][x][B] = 0;
                output->data[y][x][A] = 0;
            }
            else {
                output->data[y][x][R] = input->data[CLAMP(converted_y, 0, output->height-1)][CLAMP(converted_x, 0, output->width-1)][R];
                output->data[y][x][G] = input->data[CLAMP(converted_y, 0, output->height-1)][CLAMP(converted_x, 0, output->width-1)][G];
                output->data[y][x][B] = input->data[CLAMP(converted_y, 0, output->height-1)][CLAMP(converted_x, 0, output->width-1)][B];
                output->data[y][x][A] = input->data[CLAMP(converted_y, 0, output->height-1)][CLAMP(converted_x, 0, output->width-1)][A];
            }
        }
    }
}

void swirl(IMAGE *input, IMAGE *output) {
    output->height = input->height;
    output->width = input->width;
    double centerX = output->width / 2.0;
    double centerY = input->height / 2.0;
    for (int y = 0; y < output->height; y++) {
        for (int x = 0; x < input->width; x++) {
            double x1 = x - centerX;
            double y1 = y - centerY;
            double r = sqrt((x1 * x1) + (y1 * y1));
            double theta = atan2(y1, x1); // atan2 returns the angle in radians
            double max = sqrt((centerX * centerX / 4) + (centerY * centerY / 4));
            theta = theta + .8 * PI * (r / max);
            int rotated_x = (int)(r * cos(theta) + centerX);
            int rotated_y = (int)(r * sin(theta) + centerY);
            if (rotated_x < 0 || rotated_x >= input->width || rotated_y < 0 || rotated_y >= input->height) {
                output->data[y][x][R] = 0;
                output->data[y][x][G] = 0;
                output->data[y][x][B] = 0;
                output->data[y][x][A] = 0;
            }
            else {
                output->data[y][x][R] = input->data[CLAMP(rotated_y, 0, output->height-1)][CLAMP(rotated_x, 0, output->width-1)][R];
                output->data[y][x][G] = input->data[CLAMP(rotated_y, 0, output->height-1)][CLAMP(rotated_x, 0, output->width-1)][G];
                output->data[y][x][B] = input->data[CLAMP(rotated_y, 0, output->height-1)][CLAMP(rotated_x, 0, output->width-1)][B];
                output->data[y][x][A] = input->data[CLAMP(rotated_y, 0, output->height-1)][CLAMP(rotated_x, 0, output->width-1)][A];
            }
            
        }
    }
}

void solarizing(IMAGE *input, IMAGE *output) {
    output->height = input->height;
    output->width = input->width;
    for (int y = 0; y < input->height; y++) {
        for (int x = 0; x < input->width; x++) {
            output->data[y][x][R] = (input->data[y][x][R] > 128) ? (2 * (input->data[y][x][R] - 128)) : (2 * (128 - input->data[y][x][R]));
            output->data[y][x][G] = (input->data[y][x][G] > 128) ? (2 * (input->data[y][x][G] - 128)) : (2 * (128 - input->data[y][x][G]));
            output->data[y][x][B] = (input->data[y][x][B] > 128) ? (2 * (input->data[y][x][B] - 128)) : (2 * (128 - input->data[y][x][B]));
        }
    }
}

void saturation(IMAGE *input, IMAGE *output) {
    output->height = input->height;
    output->width = input->width;
    for (int y = 0; y < input->height; y++) {
        for (int x = 0; x < input->width; x++) {
            float v = (input->data[y][x][R] + input->data[y][x][G] + input->data[y][x][B]) / 3.0;
            output->data[y][x][R] = input->data[y][x][R] + (0.4 * (input->data[y][x][R] - v));
            output->data[y][x][G] = input->data[y][x][G] + (0.4 * (input->data[y][x][G] - v));
            output->data[y][x][B] = input->data[y][x][B] + (0.4 * (input->data[y][x][B] - v));
        }
    }
}

void subtract(IMAGE *input1, IMAGE *input2, IMAGE *output) {
    output->height = input1->height;
    output->width = input1->width;
    for (int y = 0; y < input1->height; y++) {
        for (int x = 0; x < input1->width; x++) {
            output->data[y][x][R] = input1->data[y][x][R] - input2->data[y][x][R];
            output->data[y][x][G] = input1->data[y][x][G] - input2->data[y][x][G];
            output->data[y][x][B] = input1->data[y][x][B] - input2->data[y][x][B];
        }
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

    if (test == 0)
    {
        read_ppm(&Textures[0], "../Images/earth.ascii.ppm");
        textureIndex = 0;
        draw_GLquad();
    }
    else if (test == 1)
    {
        draw_checkerboard( &Textures[1] );
        textureIndex = 1;
        draw_GLquad();
    }
    else if (test == 2)
    {
        draw_gradient( &Textures[2], 0, 0, 0, 0, 0, 1 );
        textureIndex = 2;
        draw_GLquad();
    }
    else if (test == 3)
    {
        copyImage(&Textures[0], &Textures[3]);
        textureIndex = 3;
        draw_GLquad();
    }
    else if (test == 4)
    {
        negativeImage(&Textures[0], &Textures[4]);
        textureIndex = 4;
        draw_GLquad();
    }
    else if (test == 5)
    {
        rotate(&Textures[0], &Textures[5], 45);
        textureIndex = 5;
        draw_GLquad();
    }
    else if (test == 6)
    {
        flipvertical(&Textures[0], &Textures[6]);
        textureIndex = 6;
        draw_GLquad();
    }
    else if (test == 7)
    {
        luminosity(&Textures[0], &Textures[7]);
        textureIndex = 7;
        draw_GLquad();
    }
    else if (test == 8)
    {
        sepia(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 9)
    {
        blur(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 10)
    {
        min_neighbor(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 11)
    {
        max_neighbor(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 12)
    {
        fisheye(&Textures[0], &Textures[test], 0);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 13)
    {
        swirl(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 14)
    {
        solarizing(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    }
    else if (test == 15)
    {
        saturation(&Textures[0], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
    } 
    else if (test == 16)
    {
        subtract(&Textures[0], &Textures[1], &Textures[test]);
        textureIndex = test;
        draw_GLquad();
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
    case 'y':
        test--;
        if (test < 0)
            test = 16;
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
