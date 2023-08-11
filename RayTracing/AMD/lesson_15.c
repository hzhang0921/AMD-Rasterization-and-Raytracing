/*
 *
 * point.c - simple GLUT app that draws one frame with a single point at origin
 *
 * To build:  gcc -framework OpenGL -framework GLUT lesson_8.c -o lesson_8
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
#include <strings.h>

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/

#define WIN_WIDTH			800
#define WIN_HEIGHT			800

#define PI					3.14159265358979323846264
#define EPSILON				0.000000001
#define MILLION             1000000

#define SWAP(a,b)			{ float tmp = (a); (a) = (b); (b) = tmp; }
#define ABS(x)				(((x) < 0) ? -(x) : (x))
#define SQR(a)		 		((a)*(a))
#define MIN(a, b)           (((a) < (b)) ? (a) : (b))
#define MAX(a, b)           (((a) > (b)) ? (a) : (b))
#define CLAMP(a, low, high) (((a) < (low)) ? (low) : (((a) > (high)) ? (high) : (a) ))
#define MIN3(a,b,c)         (((a) < (b)) ? (((a) < (c)) ? (a) : (c)) : (((b) < (c)) ? (b) : (c)))
#define MAX3(a,b,c)         (((a) > (b)) ? (((a) > (c)) ? (a) : (c)) : (((b) > (c)) ? (b) : (c)))

#define X					0
#define Y					1
#define Z					2
#define W					3

#define R					0
#define G					1
#define B					2
#define A					3

#define S					0
#define T					1

#define NO_TEXTURE			0
#define TEX_ENABLE			1
#define TEX_MODULATE		2

#define MAX_TEXTURE			4

#define BIAS				0.0001

#define OPAQUE					0
#define REFLECT					1
#define REFLECT_AND_REFRACT		2
#define EMISSIVE                3

#define NO_SHADOWS				0
#define HARD_SHADOWS			1
#define SOFT_SHADOWS			2

/*************************************************************************/
/* structs                                                               */
/*************************************************************************/
typedef struct image {
    int   width;
    int   height;
	unsigned char data[4096][4096][4];
} IMAGE;

typedef struct ray {
	float start[4]; // starting coordinates
    float direction[4]; // direction vector
} RAY;

typedef struct sphere_RT {
	int		material;
	int		textured;
	int		tex_index;
	float	index_of_refraction;
	float 	center[4];
	float 	radius;
	float 	color[4];
	float 	fuzz;	
} SPHERE_RT;

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
int half_window_size = (WIN_WIDTH/2);
int Mojave_WorkAround = 1;
int draw_one_frame = 1;

float color_buffer[WIN_HEIGHT][WIN_WIDTH][4];

float near = 1.0;
float far = 100.0;

float center[4] 	= { 0, 0,  3, 0 };
float L[4] 			= { 1, 1,  0, 0 };
float V[4] 			= { 0, 0, -1, 0 };
float H[4] 			= { 0, 0,  0, 0 };
float eye_vect[4] 	= { 0, 0, -1, 0 };
float eye_pos[4] 	= { 0, 0,  0, 1 };
float light_pos[4] 	= { 8, 5,-10, 1 };
float zero_vect[4] 	= { 0, 0,  0, 0 };

IMAGE textures[MAX_TEXTURE];
IMAGE *cur_texture = &textures[0];
int texturing = 0;
int phong_lighting = 1;
int shadows = 0;

float white[4] 	= { 1,   1,   1, 1 };
float blue[4] 	= { 0.5, 0.7, 1, 1 };

SPHERE_RT sphere_list[100];
int num_spheres = 0;

int max_ray_depth 		= 10;
int samples_per_pixel 	= 64;

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

/*************************************************************************/
/* vector functions                                                      */
/*************************************************************************/
void vect_print( char *name, float v[4] )
{
    printf( "%s = [%f,%f,%f,%f]\n", name, v[0], v[1], v[2], v[3] );
}

void vect_set( float result[4], float x, float y, float z, float w )
{
    result[0] = x;
    result[1] = y;
    result[2] = z;
    result[3] = w;
}

void vect_copy( float v[4], float result[4] )
{
    result[0] = v[0];
    result[1] = v[1];
    result[2] = v[2];
    result[3] = v[3];
}

void vect_add( float v0[4], float v1[4], float result[4] )
{
    result[0] = v0[0] + v1[0];
    result[1] = v0[1] + v1[1];
    result[2] = v0[2] + v1[2];
    result[3] = v0[3] + v1[3];
}

void vect_add_three( float v0[4], float v1[4], float v2[4], float result[4] )
{
    result[0] = v0[0] + v1[0] + v2[0];
    result[1] = v0[1] + v1[1] + v2[1];
    result[2] = v0[2] + v1[2] + v2[2];
    result[3] = v0[3] + v1[3] + v2[3];
}

void vect_subtract( float v0[4], float v1[4], float result[4] )
{
    result[0] = v0[0] - v1[0];
    result[1] = v0[1] - v1[1];
    result[2] = v0[2] - v1[2];
    result[3] = v0[3] - v1[3];
}

void vect_add_scalar( float s, float v1[4], float result[4] )
{
    result[0] = v1[0] + s;
    result[1] = v1[1] + s;
    result[2] = v1[2] + s;
    result[3] = v1[3] + s;
}

void vect_sub_scalar( float s, float v1[4], float result[4] )
{
    result[0] = v1[0] - s;
    result[1] = v1[1] - s;
    result[2] = v1[2] - s;
    result[3] = v1[3] - s;
}

void vect_multiply_scalar( float s, float v1[4], float result[4] )
{
    result[0] = s * v1[0];
    result[1] = s * v1[1];
    result[2] = s * v1[2];
    result[3] = s * v1[3];
}

void vect_multiply( float v0[4], float v1[4], float result[4] )
{
    result[0] = v0[0] * v1[0];
    result[1] = v0[1] * v1[1];
    result[2] = v0[2] * v1[2];
    result[3] = v0[3] * v1[3];
}

void vect_div_scalar( float a, float v[4], float result[4] )
{
    float tmp;
    
    if( ABS(a) < EPSILON )
    {
        tmp = 1.0/((a > 0.0) ? EPSILON : -EPSILON);
    }
    else
    {
        tmp = 1.0/a;
    }
    
    result[X] = tmp * v[X];
    result[Y] = tmp * v[Y];
    result[Z] = tmp * v[Z];
    result[W] = tmp * v[W];
}

float vect_length( float v[4] )
{
    return( sqrt( SQR(v[X]) + SQR(v[Y]) + SQR(v[Z]) ) );
}

void vect_normalize( float v[4], float result[4] )
{
    float length = vect_length( v );
    
    vect_div_scalar( length, v, result );
    result[W] = 0.0;
}

float vect_dot_product( float a[4], float b[4] )
{
    return( (a[X] * b[X]) + (a[Y] * b[Y]) + (a[Z] * b[Z]) );
}

/*
 * vect_cross_product()
 */
void vect_cross_product( float a[4], float b[4], float r[4] )
{
    /*
     *   Set r to the vector cross product of v1 and v2
     */
    r[X] = (a[Y] * b[Z]) - (a[Z] * b[Y]);
    r[Y] = (a[Z] * b[X]) - (a[X] * b[Z]);
    r[Z] = (a[X] * b[Y]) - (a[Y] * b[X]);
    r[W] = 0.0;
}

/*
 * vect_interpolate()
 */
void vect_interpolate( float factor, float a[4], float b[4], float r[4] )
{
    r[X] = (a[X] * factor) + (b[X] * (1-factor));
    r[Y] = (a[Y] * factor) + (b[Y] * (1-factor));
    r[Z] = (a[Z] * factor) + (b[Z] * (1-factor));
    r[W] = (a[W] * factor) + (b[W] * (1-factor));
}

/*
 * vect_random()
 */
void vect_random( float result[4] )
{
    result[0] = random_float(0,1);
    result[1] = random_float(0,1);
    result[2] = random_float(0,1);
    result[3] = 0;
}

/*
 * vect_random_min_max()
 */
void vect_random_min_max( float result[4], float min, float max )
{
    result[0] = random_float( min, max );
    result[1] = random_float( min, max );
    result[2] = random_float( min, max );
    result[3] = 0;
}

/*
 * vect_sphere_map()
 */
void vect_sphere_map( float v[4], float tex[4] )
{
    float r     = sqrt( (v[X]*v[X]) + (v[Z]*v[Z]) );
    float theta = asin( v[Z]/r );
    float phi   = asin( v[Y] );
    
    tex[S] = (theta + PI/2.0)/(2.0 * PI) + ((v[X] > 0.0) ? 0.0 : 0.5);
    tex[T] =   (phi + PI/2.0)/(PI);
    tex[2] = 0;
    tex[3] = 0;
}

/*
 * vect_reflect()
 */
void vect_reflect( float I[4], float N[4], float r[4] )
{
    float factor = 2.0 * vect_dot_product( I, N );
    float tmp[4];
 
    /*
     * reflected_vector = incident_vector - 2 * incident_vector . normal * normal
     */
	vect_multiply_scalar( factor, N, tmp );
    vect_subtract( I, tmp, r );
}

/*
 * vect_refract()
 */ 			
int vect_refract( float I[4], float N[4], float index_of_refraction, float r[4] )
{
    float cosi = vect_dot_product(I, N);
    float etat = index_of_refraction;
    float etai = 1.0;
    float n[4];
    float tmp[4];

   
    vect_copy( N, n );

    if( cosi < 0 )
    {
    	cosi = -cosi; 
    }
    else
    {
    	SWAP( etai, etat );
    	vect_multiply_scalar( -1, n, n ); 
    }
    
    float eta = etai / etat;
    float k = 1 - SQR(eta) * (1 - SQR(cosi));
    
    if( k < 0 )
    {
    	return( 0 );
    }

    vect_multiply_scalar( eta, I, r );
    vect_multiply_scalar( eta * cosi - sqrtf(k), n, tmp );
    vect_add( r, tmp, r );
    
    return( 1 );
}

/*************************************************************************/
/* buffer functions                                                      */
/*************************************************************************/
/*
 * clear_color_buffer()
 */
void clear_color_buffer( float r, float g, float b, float a )
{
	for( int y = 0 ; y < WIN_HEIGHT; y++ )
	{ 
		for( int x = 0 ; x < WIN_WIDTH; x++ )
		{
			vect_set( color_buffer[y][x], r, g, b, a );
		}
	}
}

/*
 * draw_color_buffer()
 */
void draw_color_buffer( void )
{
	float x_shift = WIN_WIDTH/2;
	float y_shift = WIN_HEIGHT/2;

	/*
	 * send all points to OpenGL
	 */
	glBegin(GL_POINTS);

	for( int y = 0 ; y < WIN_HEIGHT; y++ )
	{ 
		for( int x = 0 ; x < WIN_WIDTH; x++ )
		{
			glColor4fv( color_buffer[y][x] );
			glVertex2f( (x - x_shift)+0.5, (y - y_shift)+0.5 );
		}
	}
	
	glEnd();
}

/*************************************************************************/
/* image functions                                                       */
/*************************************************************************/
/*
 * image_checkerboard()
 */
void image_checkerboard( IMAGE *image, int width, int height )
{
	unsigned char c;
	
	image->width  = width;
	image->height = height;
	
	for( int j = 0; j < height; j++ )
	{
		for( int i = 0; i < width; i++ )
		{
			if( j/16 & 1 )
				c = (i/16 & 1) ? 0 : 255;
			else
				c = (i/16 & 1) ? 255 : 0;
				
			image->data[j][i][R] = c;
			image->data[j][i][G] = c; 
			image->data[j][i][B] = c; 
			image->data[j][i][A] = 255;
		}
	}
}

/*
 * image_checkerboard_color()
 */
void image_checkerboard_color( IMAGE *image, int width, int height, char color1[4], char color2[4] )
{
	unsigned char c;
	
	image->width  = width;
	image->height = height;
	
	for( int j = 0; j < height; j++ )
	{
		for( int i = 0; i < width; i++ )
		{
			if( j/16 & 0x1 )
			{	
				if( i/16 & 0x1 )
				{
					image->data[j][i][R] = color1[R];
					image->data[j][i][G] = color1[G];
					image->data[j][i][B] = color1[B]; 
					image->data[j][i][A] = color1[A];
				}
				else
				{
					image->data[j][i][R] = color2[R];
					image->data[j][i][G] = color2[G];
					image->data[j][i][B] = color2[B]; 
					image->data[j][i][A] = color2[A];
				}
			}
			else
			{
				if( i/16 & 0x1 )
				{
					image->data[j][i][R] = color2[R];
					image->data[j][i][G] = color2[G];
					image->data[j][i][B] = color2[B]; 
					image->data[j][i][A] = color2[A];
				}
				else
				{
					image->data[j][i][R] = color1[R];
					image->data[j][i][G] = color1[G];
					image->data[j][i][B] = color1[B]; 
					image->data[j][i][A] = color1[A];
				}
			}
		}
	}
}

/*
 * image_copy()
 */
void image_copy( IMAGE *input, IMAGE *output )
{
	output->width  = input->width;
	output->height = input->height;
	
	for( int j = 0; j < input->height; j++ )
	{
		for( int i = 0; i < input->width; i++ )
		{
			output->data[j][i][R] = input->data[j][i][R];
			output->data[j][i][G] = input->data[j][i][G];
			output->data[j][i][B] = input->data[j][i][B];
			output->data[j][i][A] = input->data[j][i][A];
		}
	}
}

/*
 * read_ppm()
 */
void read_ppm( IMAGE *image, char *name )
{
    FILE            *fp;
    int				c;
    char            buffer[16];
    int				w;
    int				h;
    int             r, g, b;
    int             max;
    unsigned char   *data, *p;
    float 			scale;
    
    if( (fp = fopen(name, "rb")) == NULL )
    {
        image->width    = 0;
        image->height   = 0;
        return;
    }

    fscanf( fp, "%s\n", buffer );

	/*
	 * skip lines that start with '#'
	 */
    while( 1 )
    {
		c = fgetc( fp );
		if( c == '#' )
		{
			do{ c = fgetc( fp ); } while( c != '\n' );
		}
		else
		{
			ungetc( c, fp );
			break;
		}
 	}
    
    fscanf( fp, "%d %d\n", &w, &h );
    fscanf( fp, "%d\n", &max );

	scale = 255.0/max;
	
	image->width    = w;
	image->height   = h;
     	   
    if( buffer[0] == 'P' && buffer[1] == '3' )
    {
        for( int j = 0; j < h; j++ )
        {
            for( int i = 0; i < w; i++ )
            {
                fscanf( fp, "%d %d %d", &r, &g, &b );
                
                r *= scale;
				g *= scale;
				b *= scale;
				                
                image->data[j][i][R] = CLAMP( r, 0, 255 );
                image->data[j][i][G] = CLAMP( g, 0, 255 );
                image->data[j][i][B] = CLAMP( b, 0, 255 );
                image->data[j][i][A] = 255;
            }
        }
    }
    else
    {
        data = (unsigned char *)malloc( w * h * 3 );
        
        if( data == NULL )
		{
			image->width    = 0;
			image->height   = 0;
			
			fclose( fp );
			return;
		}
        
        p = data;
        
        fread( (void *)data, 3, w * h, fp );
 
        for( int j = 0; j < h; j++ )
        {
            for( int i = 0; i < w; i++ )
            {
                r = *p++ * scale;
                g = *p++ * scale;
                b = *p++ * scale;
 
                image->data[h-j-1][i][R] = CLAMP( r, 0, 255 );
                image->data[h-j-1][i][G] = CLAMP( g, 0, 255 );
                image->data[h-j-1][i][B] = CLAMP( b, 0, 255 );
                image->data[h-j-1][i][A] = 255;
            }
        }
        
        free( data );
    }
    
    fclose( fp );
}

/*
 * sample_texture()
 */
void sample_texture( IMAGE *tex, float s, float t, float tex_color[4] )
{
	/*
	 * nearest texture sampling
	 */
	int u = s * (tex->width-1);
	int v = t * (tex->height-1);
		
	u = CLAMP( u, 0, tex->width-1  );
	v = CLAMP( v, 0, tex->height-1 );

	tex_color[R] = tex->data[v][u][R]/255.0;
	tex_color[G] = tex->data[v][u][G]/255.0;
	tex_color[B] = tex->data[v][u][B]/255.0;
	tex_color[A] = tex->data[v][u][A]/255.0;
}

/*************************************************************************/
/* Graphics Lesson ray tracing functions                                 */
/*************************************************************************/
/*
 * ray_calc_hit_point()
 */
void ray_calc_hit_point( float t, RAY *r, float p[4] )
{
	float tmp[4];
	
	vect_multiply_scalar( t, r->direction, tmp );
	vect_add( r->start, tmp, p );
}

/*
 * solveQuadratic_RT()
 */ 			
int solveQuadratic_RT( float a, float b, float c, float *t )
{
    float discr = b * b - 4 * a * c;
    float t0, t1;
    
    if( discr < 0 )
    {
    	return( 0 );
    }
    
    if( discr == 0 )
    {
    	t0 = - 0.5 * b / a;
    	t1 = t0;
    }
    else
    {
        float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
        t0 = q / a;
        t1 = c / q;
    }
    
    *t = MIN( t0, t1 );
        
    return( 1 );
}

/*
 * ray_hits_sphere()
 */ 			
int ray_hits_sphere( RAY *ray, SPHERE_RT *sphere, float *t )
{
	float t0, t1; // solutions for t if the ray intersects
	float L[4];

	vect_subtract( ray->start, sphere->center, L );
	float a = vect_dot_product( ray->direction, ray->direction );
	float b = 2 * vect_dot_product( ray->direction, L );
	float c = vect_dot_product( L, L ) - (sphere->radius * sphere->radius);
	
	return( solveQuadratic_RT(a, b, c, t) );
}

/*
 * init_spheres()
 */
void init_spheres( void )
{
	SPHERE_RT 	*sphere;
	int			i = 0;

    sphere = &sphere_list[i++];
    sphere->material	= OPAQUE;
	sphere->textured	= texturing;
	sphere->tex_index	= 0;
	sphere->radius 		= 999.5;
	sphere->index_of_refraction = 1.0;
	sphere->fuzz		= 0.0;
	vect_set( sphere->center, 0, -1000, 5, 1 );
	vect_set( sphere->color, 0, 1, 0, 1 );


    sphere = &sphere_list[i++];	
	sphere->material	= REFLECT_AND_REFRACT;
	sphere->textured	= texturing;
	sphere->tex_index	= 0;
	sphere->radius 		= 1;
	sphere->index_of_refraction = 2.5;
	sphere->fuzz		= 0.0;
	vect_set( sphere->center, 0, 0.5, 5, 1 );
	vect_set( sphere->color, 1, 0, 0, 1 );

    sphere = &sphere_list[i++];
    sphere->material	= REFLECT;
	sphere->textured	= texturing;
	sphere->tex_index	= 0;
	sphere->radius 		= 1;
	sphere->index_of_refraction = 1.0;
	sphere->fuzz		= 0.5;
	vect_set( sphere->center, 2.5, 0.5, 5, 1 );
	vect_set( sphere->color, 0, 1, 0, 1 );

    sphere = &sphere_list[i++];
    sphere->material	= OPAQUE;
	sphere->textured	= texturing;
	sphere->tex_index	= 0;
	sphere->radius 		= 1;
	sphere->index_of_refraction = 1.0;
	sphere->fuzz		= 0.0;
	vect_set( sphere->center, 1.5, 0.5, 10, 1 );
	vect_set( sphere->color, 0, 0, 1, 1 );


	num_spheres = i;
}

/*
 * vect_random_in_unit_cube()
 */ 			
void vect_random_in_unit_cube(float v[4]) 
{
	v[X] = random_float(-1, 1);
	v[Y] = random_float(-1, 1);
	v[Z] = random_float(-1, 1);
	v[W] = 0.0;
}

/*
 * vect_random_in_unit_sphere()
 */ 			
void vect_random_in_unit_sphere(float v[4]) 
{
	while (1)
	{
		vect_random_in_unit_cube(v);
		float length_squared = vect_dot_product(v, v);
		if (length_squared <= 1.0)
		{
			return;
		}
	}
}

/*
 * vect_random_in_unit_hemisphere()
 */ 			
void vect_random_in_unit_hemisphere(float n[4], float v[4])
{
	vect_random_in_unit_sphere(v);
	if (vect_dot_product(v, n) > 0.0) 
	{
		vect_multiply_scalar(-1, v, v);
	}
}

/*
 * find_color_of_closest_hit()
 */ 			
void find_color_of_closest_hit( RAY *ray, float color[4], int depth )
{
	float	closest_so_far = MAXFLOAT;
	float	t;
	int		hit = -1;
	
	if( depth == 0 )
	{
		vect_set( color, 0, 0, 0, 1 );
		return;
	}
	
	/*
	 * find closest intersection (if any)
	 */
    for( int i = 0; i < num_spheres; i++ )
    {
    	SPHERE_RT *sphere = &sphere_list[i];
    	
		if( ray_hits_sphere( ray, sphere, &t ) && t > BIAS && t < closest_so_far )
		{
			hit 			= i;
			closest_so_far	= t;
		}
	}

	/*
	 * if no intersection, return sky blue gradient
	 */
	if( hit == -1 )
	{
		vect_interpolate( 0.5*(ray->direction[Y] + 1.0), blue, white, color );
		return;
	}

	/*
	 * handle the intersection of ray with sphere
	 */
	SPHERE_RT 	*sphere = &sphere_list[hit];
	float 		hit_point[4];
	float		N[4];
    float       reflectance = 1.0;
    float       I = 1.0;

	/*
	 * calculate hit point and normal
	 */
	ray_calc_hit_point( closest_so_far, ray, hit_point );
	vect_subtract( hit_point, sphere->center, N );
	vect_normalize( N, N );

	/*
	 * calculate color of surface we hit
	 */
	vect_copy( sphere->color, color );
	
	if( sphere->material == EMISSIVE )
	{
		return;
	}

	if( phong_lighting )
	{
		/*
		 * calculate ambient, diffuse and specular lighting
		 */
		float Ka = 0.1;
		float Kd = vect_dot_product( N, L );
		float Ks = vect_dot_product( N, H );
	
		Kd = MAX( Kd, 0 );
		Ks = MAX( Ks, 0 );
		Ks = pow( Ks, 128 );
	
		I = CLAMP( Ka + Kd + Ks, 0, 1 );
	}
	else
	{
		RAY		scattered;
		float	scattered_color[4];
		float	d[4];
	
		vect_copy( hit_point, scattered.start );
		vect_random_in_unit_hemisphere( N, d );
		vect_normalize( d, scattered.direction );
	
		find_color_of_closest_hit( &scattered, scattered_color, depth-1 );

		vect_multiply( color, scattered_color, color );
	}
	
	if( texturing && sphere->textured )
	{
		IMAGE 	*image = &textures[sphere->tex_index];
		float	tex[4];
		float	tex_color[4];
		
		vect_sphere_map( N, tex );
		sample_texture( image, tex[S], tex[T], tex_color );
		
		if( texturing == TEX_MODULATE )
		{
			vect_multiply( color, tex_color, color );
		}
		else
		{
			vect_copy( tex_color, color );
		}
	}

	/*
	 * calculate reflections
	 */
	if( sphere->material == REFLECT )
	{
		RAY		reflected;
		float	reflect_color[4];

		vect_copy( hit_point, reflected.start );
		vect_reflect( ray->direction, N, reflected.direction );
		if (sphere-> fuzz > 0.0) 
		{
			float rand[4];
			vect_random_in_unit_sphere(rand);
			vect_multiply_scalar(sphere->fuzz, rand, rand);
			vect_add(reflected.direction, rand, rand);
			vect_normalize(rand, reflected.direction);
		}
	
		find_color_of_closest_hit( &reflected, reflect_color, depth-1 ); // recursively track those reflections
	
		vect_interpolate( 0.1, color, reflect_color, color );
	}

	/*
	 * calculate refractions
	 */
	if( sphere->material == REFLECT_AND_REFRACT )
	{
		RAY		refracted;
		float	refract_color[4];
		float	kr;
		
		vect_copy( hit_point, refracted.start );
		
		if( vect_refract( ray->direction, N, sphere->index_of_refraction, refracted.direction ) )
		{
			find_color_of_closest_hit( &refracted, refract_color, depth-1 );

			vect_interpolate( 0.1, color, refract_color, color );
		}
	}

    if (shadows == HARD_SHADOWS)
    {
        /*
        * calculate shadows
        */
        RAY		shadow_ray;
        float	distance_to_light;
        int		is_shadowed = 0;

        vect_copy( hit_point, shadow_ray.start );
        vect_subtract( light_pos, hit_point, shadow_ray.direction );
        distance_to_light = vect_length( shadow_ray.direction );
        vect_normalize( shadow_ray.direction, shadow_ray.direction );
        
        for( int i = 0; i < num_spheres; i++ )
        {
            SPHERE_RT *sphere = &sphere_list[i];
            
            if( ray_hits_sphere( &shadow_ray, sphere, &t ) && t > BIAS && t < distance_to_light )
            {
                is_shadowed = 1;
                break;
            }
        }

        if( is_shadowed )
        {
            I *= 0.6;
	    }
    }

    if (shadows == SOFT_SHADOWS)
    {
        /*
        * calculate soft shadows()
        */
        RAY shadow_ray;
        float distance_to_light;
        float jitter_light[4];
        int shadow_count = 0;
        int num_shadow_rays = 10;

        for (int k = 0; k < num_shadow_rays; k++) 
        {
            vect_random_in_unit_sphere(jitter_light);
            vect_add(light_pos, jitter_light, jitter_light);
            vect_copy(hit_point, shadow_ray.start);
            vect_subtract(jitter_light, hit_point, shadow_ray.direction);
            distance_to_light = vect_length(shadow_ray.direction);
            vect_normalize(shadow_ray.direction, shadow_ray.direction);
            for (int i = 0; i < num_spheres; i++) {
                SPHERE_RT *sphere = &sphere_list[i];
                if (ray_hits_sphere(&shadow_ray, sphere, &t) && t > BIAS && t < distance_to_light) {
                    shadow_count++;
                    break;
                }
            }
        }
        I *= 1.0 - shadow_count / (float)num_shadow_rays;			
    }

	vect_multiply_scalar( I, color, color );
}

/*
 * ray_trace_graphics_lessons()
 */
void ray_trace_graphics_lessons( void )
{
	RAY 	ray;	
	float 	sum[4];
	float 	color[4];
	float	x_jitter = 0;
	float	y_jitter = 0;
	
	/*
	 * set up list of sphere in "world"
	 */
	init_spheres();
	
	vect_set( ray.start, 0, 0, 0, 1 );
	
	/*
	 * for every pixel, generate a ray and see what objects it intersects
	 */ 
	for( int j = 0; j < WIN_HEIGHT; j++ )
	{
		for( int i = 0; i < WIN_WIDTH; i++ )
		{	
			vect_set( sum, 0, 0, 0, 0 );
			
			for( int s = 0; s < samples_per_pixel; s++ )
			{
				if( samples_per_pixel > 1 )
				{
					x_jitter = random_float(0,1);
					y_jitter = random_float(0,1);
				}
				
				float u = (i+x_jitter)/(float)(WIN_WIDTH-1);		// normalize to range 0.0 -> 1.0
				float v = (j+y_jitter)/(float)(WIN_HEIGHT-1);		// normalize to range 0.0 -> 1.0
			
				vect_set( ray.direction, (2*u)-1, (2*v)-1, near, 0 );
				vect_normalize( ray.direction, ray.direction );

				find_color_of_closest_hit( &ray, color, max_ray_depth );
				
				vect_add( sum, color, sum );
			}
			
			vect_multiply_scalar( 1.0/samples_per_pixel, sum, color_buffer[j][i] );
		}
	
		printf( "%6.2f%%\n", 100.0 * j/(float)(WIN_HEIGHT-1) );
	}
}

/*************************************************************************/
/* Opengl functions                                                      */
/*************************************************************************/

/*
 * init_gl_state
 */
void init_gl_state( void )
{
    float clear_color[4] = { 0, 0, 0, 1 };
    
    /*
     * GL draw state
     */
    glClearColor( clear_color[R], clear_color[G], clear_color[B], clear_color[A] );
    glPointSize( 2.0 );
    glColor4f( 1.0,0.0,0.0,1.0 );
    
    /*
     * GL depth state
     */
    glClearDepth( 1.0 );
    glDisable( GL_DEPTH_TEST );
    glDepthFunc( GL_LEQUAL );
    
    /*
     * GL view state
     */
    glViewport( 0, 0, WIN_WIDTH, WIN_HEIGHT );

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( -half_window_size, half_window_size, -half_window_size, half_window_size, -40, 40 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

    /*
     * GL texture state
     */ 
	glBindTexture( GL_TEXTURE_2D, 0 );
    glDisable( GL_TEXTURE_2D );
  	
    /*
     * culling state
     */
    glDisable( GL_CULL_FACE );
    
    /*
     * GL fog state
     */
    glFogfv( GL_FOG_COLOR, clear_color );
    glDisable( GL_FOG );
    
    /*
     * GL blend state
     */
    glBlendEquation( GL_FUNC_ADD );
    glBlendColor( 0.5, 0.5, 0.5, 1.0 );
    glBlendFunc( GL_CONSTANT_COLOR, GL_CONSTANT_COLOR );
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
 
	// clear pixels
 	clear_color_buffer( 0, 0, 0, 0 );
 	
	// do ray tracing
 	ray_trace_graphics_lessons();
 	
	// copy pixels to GL
 	draw_color_buffer();
 	
 	// tell GLUT to display pixels
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
        case 't':       texturing = (texturing+1)%3;		draw_one_frame = 1;		glutPostRedisplay();    break;
        case 'p':		phong_lighting = 1-phong_lighting;	draw_one_frame = 1;		glutPostRedisplay();    break;
        case 's':		shadows = (shadows+1)%3;			draw_one_frame = 1;		glutPostRedisplay();    break;
    	case ' ':											draw_one_frame = 1;     glutPostRedisplay();    break;
		case 'q':       exit( 0 );																			break;
        case '\033':    exit( 0 );																			break;
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

	
	if( Mojave_WorkAround )
	{
		// Necessary for Mojave. Has to be different dimensions than in glutInitWindowSize();
		glutReshapeWindow( half_window_size * 2, half_window_size * 2 );
	}

	/*
	 * set up texture for future use
	 */	
	image_checkerboard( &textures[0], 256, 256 );
	
	/*
	 * set up light vectors
	 */	
	vect_normalize( light_pos, L );
	vect_add( L, V, H );
	vect_normalize( H, H );
		
    /*
     * setup OpenGL state
     */
	init_gl_state();
	
    /*
     * start loop that calls display() and Key() routines
     */
    glutMainLoop();

    return 0;
}
