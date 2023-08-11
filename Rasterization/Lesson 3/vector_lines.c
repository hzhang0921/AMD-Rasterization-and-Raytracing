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
int current[800][800];
IMAGE Textures[10];

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
 * draw_point()
 */
void draw_point(float x, float y)
{
    /*
     * draw points
     */
    glBegin(GL_POINTS);
    glVertex2f(x, y); // + or - to translate x and y by pixels
    glEnd();
}

/*
 * draw_point_struct()
 */
void draw_point_struct(POINT p)
{
    /*
     * draw points
     */
    glColor4f(p.color[R], p.color[G], p.color[B], 1.0);
    // float tex_color[4];
    // texture_sample(&Textures[0], p.tex_coords[0], p.tex_coords[1], tex_color);
    // glColor4f(tex_color[R], tex_color[G], tex_color[B], 1.0);
    glBegin(GL_POINTS);
    glVertex2f(p.pos[X], p.pos[Y]); // + or - to translate x and y by pixels
    glEnd();
}

/*
 * set_color()
 */
void set_color(float r, float g, float b, float a)
{
    glColor4f(r, g, b, a);
}

/*
 * draw_horizontal()
 */
void draw_horizontal(int x0, int y0, int x1, int y1)
{
    if (y0 != y1)
    {
        return;
    }
    if (x0 > x1)
    {
        SWAP(x0, x1);
    }
    for (int x = x0; x < x1; x++)
    {
        draw_point(x, y0);
    }
}

/*
 * draw_vertical()
 */
void draw_vertical(int x0, int y0, int x1, int y1)
{
    if (x0 != x1)
    {
        return;
    }
    if (y0 > y1)
    {
        SWAP(y0, y1);
    }
    for (int y = y0; y < y1; y++)
    {
        draw_point(x0, y);
    }
}

/*
 * draw_diagonal()
 */
void draw_diagonal(int x0, int y0, int x1, int y1)
{
    float dx = x1 - x0;
    float dy = y1 - y0;
    float slope = dy / dx;
    float inv_slope = dx / dy;
    float y = y0;
    float x = x0;
    if (fabs(dy) > fabs(dx))
    {
        if (y0 > y1)
        {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        x = x0;
        for (int y = y0; y < y1; y++)
        {
            draw_point(x, y);
            x += inv_slope;
        }
    }
    else
    {
        if (x0 > x1)
        {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        y = y0;
        for (int x = x0; x < x1; x++)
        {
            draw_point(x, y);
            y += slope;
        }
    }
}

/*
 * draw_line()
 */
void draw_line(int x0, int y0, int x1, int y1)
{
    if (y0 == y1)
    {
        draw_horizontal(x0, y0, x1, y1);
    }
    else if (x0 == x1)
    {
        draw_vertical(x0, y0, x1, y1);
    }
    else
    {
        draw_diagonal(x0, y0, x1, y1);
    }
}

/*
 * draw_random
 */
void draw_random()
{
    for (int i = 0; i < 100; i++)
    {
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
void draw_fan()
{
    for (int angle = 0; angle < 360; angle += 10)
    {
        float y1 = sin(angle / 360.0 * 2 * PI);
        float x1 = cos(angle / 360.0 * 2 * PI);
        int radius = 200;
        draw_line(0, 0, radius * x1, radius * y1);
    }
}

/*
 * draw_coordinate()
 */
void draw_coordinate()
{
    int width = 10;
    for (int x = -400; x < 400; x += width)
    {
        draw_line(x, -400, x, 400);
    }
    for (int y = -400; y < 400; y += width)
    {
        draw_line(-400, y, 400, y);
    }
}

/*
 * draw_Bresenham()
 */
void draw_bresenham(int x0, int y0, int x1, int y1)
{
    int y = y0;
    int dx = x1 - x0;
    int dy = y1 - y0;
    int error = 2 * dy - dx;
    for (int x = x0; x < x1; x++)
    {
        draw_point(x, y);
        error += 2 * dy;
        if (error > 0)
        {
            y++;
            error -= 2 * dx;
        }
    }
}

/*
 * random_Bresenham()
 */
void random_bresenham()
{
    for (int i = 0; i < 100; i++)
    {
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
void draw_wu(int x0, int y0, int x1, int y1)
{
    float dx = x1 - x0;
    float dy = y1 - y0;
    float slope = dy / dx;
    float inv_slope = dx / dy;
    int neg_slope = 1;
    if (slope < 0)
    {
        neg_slope = -1;
    }
    if (fabs(dy) > fabs(dx))
    { // tests for steepness
        if (y0 > y1)
        {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        float x = x0;
        for (int y = y0; y < y1; y++)
        {
            set_color(fpart(x), fpart(x), fpart(x), 1);
            draw_point(ipart(x), y);
            set_color(1 - fpart(x), 1 - fpart(x), 1 - fpart(x), 1);
            draw_point(ipart(x) + neg_slope, y);
            x += inv_slope;
        }
    }
    else
    {
        if (x0 > x1)
        {
            SWAP(x0, x1);
            SWAP(y0, y1);
        }
        float y = y0;
        for (int x = x0; x < x1; x++)
        {
            set_color(fpart(y), fpart(y), fpart(y), 1);
            draw_point(x, ipart(y));
            set_color(1 - fpart(y), 1 - fpart(y), 1 - fpart(y), 1);
            draw_point(x, ipart(y) + neg_slope);
            y += slope;
        }
    }
}

/*
 * draw_antialiasing ()
 */
void draw_antialising(int x0, int y0, int x1, int y1)
{
    if (y0 == y1)
    {
        draw_horizontal(x0, y0, x1, y1);
    }
    else if (x0 == x1)
    {
        draw_vertical(x0, y0, x1, y1);
    }
    else
    {
        draw_wu(x0, y0, x1, y1);
    }
}

/*
 * random_wu()
 */
void random_wu()
{
    for (int i = 0; i < 100; i++)
    {
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
 * draw_line_color()
 */
void draw_line_color( POINT p0, POINT p1 )
{
    float 	pos_delta[4];
    float 	color_delta[4];
    float 	tex_delta[4];
    float 	pos_inc[4];
    float 	color_inc[4];
    float 	tex_inc[4];
	POINT	p;
	
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
		// subtract_vector( p1.tex_coords,   p0.tex_coords, 	 tex_delta   );

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
		subtract_vector( p1.tex_coords,   p0.tex_coords, 	 tex_delta   );

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
		subtract_vector( p1.tex_coords,   p0.tex_coords, 	 tex_delta   );

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
		subtract_vector( p1.tex_coords,   p0.tex_coords, 	 tex_delta   );

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


void random_interpolate()
{
    POINT p0;
    POINT p1;
    for (int i = 0; i < 100; i++)
    {
        p0.pos[X] = (int)random_float(-400, 400);
        p0.pos[Y] = (int)random_float(-400, 400);
        p1.pos[X] = (int)random_float(-400, 400);
        p1.pos[Y] = (int)random_float(-400, 400);
        p0.color[R] = random_float(0, 1);
        p0.color[G] = random_float(0, 1);
        p0.color[B] = random_float(0, 1);
        p0.color[A] = 1;
        p1.color[R] = 1 - p0.color[R];
        p1.color[G] = 1 - p0.color[R];
        p1.color[B] = 1 - p0.color[R];
        p1.color[A] = 1;
        draw_line_color(p0, p1);
    }
}

/*
 * draw quad()
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

    for (int y = y1; y >= y0; y--)
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

void draw_random_quad()
{
    POINT p0, p1, p2, p3;
    p0.pos[X] = -200;
    p0.pos[Y] = 200;
    p1.pos[X] = 200;
    p1.pos[Y] = 200;
    p2.pos[X] = -200;
    p2.pos[Y] = -200;
    p3.pos[X] = 200;
    p3.pos[Y] = -200;
    set_vector(p0.color, 1, 0, 0, 1);
    set_vector(p1.color, 0, 0, 1, 1);
    set_vector(p2.color, 0, 1, 0, 1);
    set_vector(p3.color, 1, 0, 1, 1);
    set_vector(p0.tex_coords, 0, 0, 0, 0);
    set_vector(p1.tex_coords, 0, 1, 0, 0);
    set_vector(p2.tex_coords, 0, 1, 0, 0);
    set_vector(p3.tex_coords, 1, 1, 0, 0);
    draw_quad(p0, p1, p2, p3);
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
        random_interpolate();
    }
    else if (test == 1)
    {
        draw_fan();
    }
    else if (test == 2)
    {
        draw_coordinate();
    }
    else if (test == 3)
    {
        random_bresenham();
    }
    else if (test == 4)
    {
        random_wu();
    }
    else if (test == 5)
    {
        draw_random();
    }
    else if (test == 6)
    {
        draw_random_quad();
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
        if (test > 6)
            test = 0;
        draw_one_frame = 1;
        glutPostRedisplay();
        break;
    case 'y':
        test--;
        if (test < 0)
            test = 6;
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

    read_ppm(&Textures[0], "../Images/earth.ascii.ppm");

    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}
