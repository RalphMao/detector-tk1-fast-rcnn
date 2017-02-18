#ifndef __EXTERN__
#define __EXTERN__

//#include "SSE2NEON.h"
#include <cstdlib>
#include <cstring>
#include "sse.hpp"
#include "wrappers.hpp"
#include <cmath>
#include <typeinfo>
#include <iostream>
#include "model.h"

#define PI 3.14159265359f

void imPad( float *A, float *B, int h, int w, int d, int pt, int pb,
  int pl, int pr, int flag, float val );

double* rgb2luv_setup( double z, double *mr, double *mg, double *mb, double &minu, double &minv, double &un, double &vn );
void rgb2luv( unsigned char *I, double *J, int n, double nrm );
void rgb2hsv( unsigned char *I, double *J, int n, double nrm ) ;
void rgb2gray( unsigned char *I, double *J, int n, double nrm );
void normalize( unsigned char *I, double *J, int n, double nrm );
double* rgbConvert( unsigned char *I, int n, int d, int flag, double nrm );

//float* rgb2luv_setup( float z, float *mr, float *mg, float *mb, float &minu, float &minv, float &un, float &vn );
void rgb2luv( unsigned char *I, float *J, int n, float nrm );
void rgb2hsv( unsigned char *I, float *J, int n, float nrm );
void rgb2gray( unsigned char *I, float *J, int n, float nrm );
void rgb2gray( double *I, float *J, int n, float nrm );
void normalize( unsigned char *I, float *J, int n, float nrm );
float* rgbConvert( unsigned char *I, int n, int d, int flag, float nrm );

 float* rgb2luv_setup( float z, float *mr, float *mg, float *mb, float &minu, float &minv, float &un, float &vn );
 void rgb2luv( float *I, float *J, int n, float nrm ) ;
 void rgb2luv_sse( float *I, float *J, int n, float nrm );
 void rgb2hsv( float *I, float *J, int n, float nrm );
 void rgb2gray( float *I, float *J, int n, float nrm );
 void rgb2gray( double *I, float *J, int n, float nrm );
 void normalize( float *I, float *J, int n, float nrm );
 float* rgbConvert( float *I, int n, int d, int flag, float nrm );
 
 void resampleCoef( int ha, int hb, int &n, int *&yas,int *&ybs, double *&wts, int bd[2], int pad);
 void imResample( float *A, float *B, int ha, int hb, int wa, int wb, int d, float r );
 void grad2( float *I, float *Gx, float *Gy, int h, int w, int d );
 void gradMag( float *I, float *M, float *O, int h, int w, int d, bool full );
 void gradMagNorm( float *M, float *S, int h, int w, float norm );
 void gradQuantize( float *O, float *M, int *O0, int *O1, float *M0, float *M1, 
	 int nb, int n, float norm, int nOrients, bool full, bool interpolate );
 void gradHist( float *M, float *O, float *H, int h, int w, int bin, int nOrients, int softBin, bool full );

 void convTri( float *I, float *O, int h, int w, int d, int r, int s );

 uint32* buildLookup( int *dims, int w );
 void buildLookupSs( uint32 *&cids1, uint32 *&cids2, int *dims, int w, int m );
 float interp( float *I, int h, int w, float x, float y );
 
#endif
