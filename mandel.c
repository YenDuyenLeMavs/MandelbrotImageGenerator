// Yen Duyen Amy Le
// 1001827177

#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h> // added this library to use threading
#include <sys/time.h> // added this library to measure how long program takes

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
// added parameter numThreads to pass in new command line argument to do multi-threading
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int numThreads ); 
// transfered code of original compute_image() into separate function compute_image_helper() for the threads to invoke
void * compute_image_helper(void * arg);

// added struct parameters because can only pass one argument to compute_image_helper()
struct parameters 
{
	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
	// for loop iterator initialized to different value (istart)
	// for each thread because each thread works on different set of pixels
	int istart; 
	// for loop halting condition different value for each thread (widthHalt)
	// because each thread works on a different set of pixels
	int widthHalt;
};

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	// add description of new command line argument -n to specify number of threads
	printf("-n <threads> Number of threads to compute the image. (default=1)\n"); 
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main( int argc, char *argv[] )
{
	// get time program begin to calculate
	// how long program takes to run later
	struct timeval begin_time;
	struct timeval end_time;
	
	gettimeofday(&begin_time, NULL);
	
	// actual mandel program starts below
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	// added new configuration value for new command line argument to do multi-threading
	// if -n is not given, assume a default of one thread
	int numThreads = 1; 

	// For each command line argument given,
	// override the appropriate configuration value.

	// added n to parse the new command line argument to do multi-threading
	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:n:h"))!=-1) { 
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			// add new case for new command line argument -n
			// atoi b/c numThreads is int not char
			case 'n': 
				numThreads = atoi(optarg); 
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	// added numThreads ???????
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s numThreads=%d\n",xcenter,ycenter,scale,max,outfile,numThreads);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	// added argument numThreads to pass in new command line argument to do multi-threading
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,numThreads); 

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	// get end time to calculate and print out
	// how long program takes to run
	gettimeofday(&end_time, NULL);
	
	long time_to_execute = (end_time.tv_sec*1000000 + end_time.tv_usec) - (begin_time.tv_sec*1000000 + begin_time.tv_usec);
	printf("This code took %ld microseconds to execute\n", time_to_execute);
	
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

// transfered code of original compute_image() into separate function compute_image_helper() for the threads to invoke
void * compute_image_helper(void * arg)
{
	struct parameters * params = (struct parameters *) arg;
	// dereferencing members of parameters struct
	struct bitmap *bm = params -> bm;
	double xmin = params -> xmin;
	double xmax = params -> xmax;
	double ymin = params -> ymin;
	double ymax = params -> ymax;
	int max = params -> max;
	int istart = params -> istart;
	int widthHalt = params -> widthHalt;
	
	// original code in compute_image() below
	// For every pixel in the image...
	int width = bitmap_width(bm);
	int height = bitmap_height(bm);
	
	int i, j;
	for(j=0;j<height;j++) {
		// i = istart instead of i = 0 and i < widthHalt instead of i < width
		// because each thread invokes compute_image_helper()
		// such that each thread can work on different vertical lines of pixels
		// to speed up generating the image (each thread's parameters has a 
		// different value for istart and a different value for widthHalt
		for(i=istart;i<widthHalt;i++) { 

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
	return NULL;
}

// added parameter numThreads to pass in new command line argument to do multi-threading
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int numThreads ) 
{
	int width = bitmap_width(bm);
	
	// put rest of original code in compute_image()
	// (aka the nested for loops) into compute_image_helper()
	// following code is to generate multiple threads
	// written based on multi-threaded.c from github
	
	// added array to hold thread IDs
	// added array to hold the parameters for each thread execution
	pthread_t tid[numThreads]; 
	struct parameters params[numThreads]; 
	
	// give each thread the arguments passed into compute_image()
	// that will be used in compute_image_helper() 
	int k;
	for(k = 0; k < numThreads; k++)
	{
		params[k].bm = bm;
		params[k].xmin = xmin;
		params[k].xmax = xmax;
		params[k].ymin = ymin;
		params[k].ymax = ymax;
		params[k].max = max;
	}
	
	// each thread processes a different set of vertical lines (columns) of pixels
	// the first line in each set is istart, and each set goes up to (not including) widthHalt
	// widthIncrement is how many vertical lines of pixels each thread gets to process
	// set istart = 0 for first thread
	// set widthHalt = widthIncrement for first thread
	int widthIncrement = width / numThreads;
	int istart = 0;
	int widthHalt = widthIncrement;
	
	for(k = 0; k < numThreads; k++)
	{
		params[k].istart = istart;
		params[k].widthHalt = widthHalt;
		// set istart for next thread created
		istart += widthIncrement;
		// set widthHalt for next thread created
		widthHalt += widthIncrement;
	}
	
	// make the widthHalt for the last thread created one higher than the pixel width
	// so the last column of pixels gets processed by the for loops in 
	// compute_image_helper()
	params[numThreads - 1].widthHalt = width + 1;

	// creating numThreads threads, each with their own set of parameters
	for(k = 0; k < numThreads; k++)
	{
		pthread_create(&tid[k], NULL, compute_image_helper, (void*) &params[k]);	
	}
	
	// waits for each thread to terminate before ending the entire program
	for(k = 0; k < numThreads; k++)
	{
		pthread_join(tid[k], NULL);
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int R = 255*i/max; // was 255
	int G = 175*i/max;
	int B = 255*i/max;
	int A = 0.8*i/max; // was 0.8
	return MAKE_RGBA(R,G,B,A);
}




