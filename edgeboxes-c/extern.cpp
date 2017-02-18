#include "extern.h"
extern int debug_count;
void imPad( float *A, float *B, int h, int w, int d, int pt, int pb,
  int pl, int pr, int flag, float val ){
  int h1=h+pt, hb=h1+pb, w1=w+pl, wb=w1+pr, x, y, z, mPad;
  int ct=0, cb=0, cl=0, cr=0;
  if(pt<0) { ct=-pt; pt=0; } if(pb<0) { h1+=pb; cb=-pb; pb=0; }
  if(pl<0) { cl=-pl; pl=0; } if(pr<0) { w1+=pr; cr=-pr; pr=0; }
  int *xs, *ys; x=pr>pl?pr:pl; y=pt>pb?pt:pb; mPad=x>y?x:y;
  bool useLookup = ((flag==2 || flag==3) && (mPad>h || mPad>w))
    || (flag==3 && (ct || cb || cl || cr ));
  // helper macro for padding
  #define PAD(XL,XM,XR,YT,YM,YB) \
  for(x=0;  x<pl; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XL+cl)*h+YT+ct]; \
  for(x=0;  x<pl; x++) for(y=pt; y<h1; y++) B[x*hb+y]=A[(XL+cl)*h+YM+ct]; \
  for(x=0;  x<pl; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XL+cl)*h+YB-cb]; \
  for(x=pl; x<w1; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XM+cl)*h+YT+ct]; \
  for(x=pl; x<w1; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XM+cl)*h+YB-cb]; \
  for(x=w1; x<wb; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XR-cr)*h+YT+ct]; \
  for(x=w1; x<wb; x++) for(y=pt; y<h1; y++) B[x*hb+y]=A[(XR-cr)*h+YM+ct]; \
  for(x=w1; x<wb; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XR-cr)*h+YB-cb];
  // build lookup table for xs and ys if necessary
  if( useLookup ) {
 	    xs = (int*) malloc(wb*sizeof(int)); int h2=(pt+1)*2*h;
    ys = (int*) malloc(hb*sizeof(int)); int w2=(pl+1)*2*w;
    //xs = (int*) wrMalloc(wb*sizeof(int)); int h2=(pt+1)*2*h;
    //ys = (int*) wrMalloc(hb*sizeof(int)); int w2=(pl+1)*2*w;
    if( flag==2 ) {
      for(x=0; x<wb; x++) { z=(x-pl+w2)%(w*2); xs[x]=z<w ? z : w*2-z-1; }
      for(y=0; y<hb; y++) { z=(y-pt+h2)%(h*2); ys[y]=z<h ? z : h*2-z-1; }
    } else if( flag==3 ) {
      for(x=0; x<wb; x++) xs[x]=(x-pl+w2)%w;
      for(y=0; y<hb; y++) ys[y]=(y-pt+h2)%h;
    }
  }
  // pad by appropriate value
  for( z=0; z<d; z++ ) {
    // copy over A to relevant region in B
    for( x=0; x<w-cr-cl; x++ )
      memcpy(B+(x+pl)*hb+pt,A+(x+cl)*h+ct,sizeof(float)*(h-ct-cb));
    // set boundaries of B to appropriate values
    if( flag==0 && val!=0 ) { // "constant"
      for(x=0;  x<pl; x++) for(y=0;  y<hb; y++) B[x*hb+y]=val;
      for(x=pl; x<w1; x++) for(y=0;  y<pt; y++) B[x*hb+y]=val;
      for(x=pl; x<w1; x++) for(y=h1; y<hb; y++) B[x*hb+y]=val;
      for(x=w1; x<wb; x++) for(y=0;  y<hb; y++) B[x*hb+y]=val;
    } else if( useLookup ) { // "lookup"
      PAD( xs[x], xs[x], xs[x], ys[y], ys[y], ys[y] );
    } else if( flag==1 ) {  // "replicate"
      PAD( 0, x-pl, w-1, 0, y-pt, h-1 );
    } else if( flag==2 ) { // "symmetric"
      PAD( pl-x-1, x-pl, w+w1-1-x, pt-y-1, y-pt, h+h1-1-y );
    } else if( flag==3 ) { // "circular"
      PAD( x-pl+w, x-pl, x-pl-w, y-pt+h, y-pt, y-pt-h );
    }
    A += h*w;  B += hb*wb;
  }
  //if( useLookup ) { wrFree(xs); wrFree(ys); }
  if( useLookup ) { free(xs); free(ys); }
  #undef PAD
}

// Constants for rgb2luv conversion and lookup table for y-> l conversion
 double* rgb2luv_setup( double z, double *mr, double *mg, double *mb,
  double &minu, double &minv, double &un, double &vn )
{
  // set constants for conversion
	const double y0=(double) ((6.0/29)*(6.0/29)*(6.0/29));
  const double a= (double) ((29.0/3)*(29.0/3)*(29.0/3));
  un=(double) 0.197833; vn=(double) 0.468331;
  mr[0]=(double) 0.430574*z; mr[1]=(double) 0.222015*z; mr[2]=(double) 0.020183*z;
  mg[0]=(double) 0.341550*z; mg[1]=(double) 0.706655*z; mg[2]=(double) 0.129553*z;
  mb[0]=(double) 0.178325*z; mb[1]=(double) 0.071330*z; mb[2]=(double) 0.939180*z;
  double maxi=(double) 1.0/270; minu=-88*maxi; minv=-134*maxi;
  // build (padded) lookup table for y->l conversion assuming y in [0,1]
  static double lTable[1064]; static bool lInit=false;
  if( lInit ) return lTable; double y, l;
  for(int i=0; i<1025; i++) {
    y = (double) (i/1024.0);
    l = y>y0 ? 116*(double)pow((double)y,1.0/3.0)-16 : y*a;
    lTable[i] = l*maxi;
  }
  for(int i=1025; i<1064; i++) lTable[i]=lTable[i-1];
  lInit = true; return lTable;
}

// Convert from rgb to luv
void rgb2luv( unsigned char *I, double *J, int n, double nrm ) {
  double minu, minv, un, vn, mr[3], mg[3], mb[3];
  double *lTable = rgb2luv_setup(nrm,mr,mg,mb,minu,minv,un,vn);
  double *L=J, *U=L+n, *V=U+n; unsigned char *R=I, *G=R+n, *B=G+n;
  for( int i=0; i<n; i++ ) {
    double r, g, b, x, y, z, l;
    r=(double)*R++; g=(double)*G++; b=(double)*B++;
    x = mr[0]*r + mg[0]*g + mb[0]*b;
//	if(i<10)
//		printf("x[%d]=%lf, ",i, (double)x);
    y = mr[1]*r + mg[1]*g + mb[1]*b;
    z = mr[2]*r + mg[2]*g + mb[2]*b;
    l = lTable[(int)(y*1024)];
    *(L++) = l; z = 1/(x + 15*y + 3*z + (double)1e-35);
    *(U++) = l * (13*4*x*z - 13*un) - minu;
    *(V++) = l * (13*9*y*z - 13*vn) - minv;
  }
}

// Convert from rgb to hsv
void rgb2hsv( unsigned char *I, double *J, int n, double nrm ) {
  double *H=J, *S=H+n, *V=S+n;
  unsigned char *R=I, *G=R+n, *B=G+n;
  for(int i=0; i<n; i++) {
    const double r=(double)*(R++), g=(double)*(G++), b=(double)*(B++);
    double h, s, v, minv, maxv;
    if( r==g && g==b ) {
      *(H++) = 0; *(S++) = 0; *(V++) = r*nrm; continue;
    } else if( r>=g && r>=b ) {
      maxv = r; minv = g<b ? g : b;
      h = (g-b)/(maxv-minv)+6; if(h>=6) h-=6;
    } else if( g>=r && g>=b ) {
      maxv = g; minv = r<b ? r : b;
      h = (b-r)/(maxv-minv)+2;
    } else {
      maxv = b; minv = r<g ? r : g;
      h = (r-g)/(maxv-minv)+4;
    }
    h*=(double) (1/6.0); s=1-minv/maxv; v=maxv*nrm;
    *(H++) = h; *(S++) = s; *(V++) = v;
  }
}

// Convert from rgb to gray
void rgb2gray( unsigned char *I, double *J, int n, double nrm ) {
  double *GR=J; unsigned char *R=I, *G=R+n, *B=G+n; int i;
  double mr=(double).2989360213*nrm, mg=(double).5870430745*nrm, mb=(double).1140209043*nrm;
  for(i=0; i<n; i++) *(GR++)=(double)*(R++)*mr + (double)*(G++)*mg + (double)*(B++)*mb;
}



// Copy and normalize only
void normalize( unsigned char *I, double *J, int n, double nrm ) {
  for(int i=0; i<n; i++) *(J++)=(unsigned char)*(I++)*nrm;
}

// Convert rgb to various colorspaces
double* rgbConvert( unsigned char *I, int n, int d, int flag, double nrm ) {
  double *J = (double*) malloc(n*(flag==0 ? (d==1?1:d/3) : d)*sizeof(double));
  int i, n1=d*(n<1000?n/10:100); double thr = double(1.001);
  if(flag>1 && nrm==1) for(i=0; i<n1; i++) if(I[i]>thr)
    std::cerr<<"For floats all values in I must be smaller than 1."<<std::endl;
  bool useSse = 0;//n%4==0 && typeid(double)==typeid(float);
  //if( flag==2 && useSse )
  //  for(i=0; i<d/3; i++) rgb2luv_sse(I+i*n*3,(float*)(J+i*n*3),n,(float)nrm);
  if( (flag==0 && d==1) || flag==1 ) normalize(I,J,n*d,nrm);
  else if( flag==0 ) for(i=0; i<d/3; i++) rgb2gray(I+i*n*3,J+i*n*1,n,nrm);
  else if( flag==2 ) for(i=0; i<d/3; i++) rgb2luv(I+i*n*3,J+i*n*3,n,nrm);
  else if( flag==3 ) for(i=0; i<d/3; i++) rgb2hsv(I+i*n*3,J+i*n*3,n,nrm);
  else std::cerr<<"Unknown flag."<<std::endl;
  return J;
}

//float* rgb2luv_setup( float z, float *mr, float *mg, float *mb,
//  float &minu, float &minv, float &un, float &vn )
//{
//  // set constants for conversion
//	const float y0=(float) ((6.0/29)*(6.0/29)*(6.0/29));
//  const float a= (float) ((29.0/3)*(29.0/3)*(29.0/3));
//  un=(float) 0.197833; vn=(float) 0.468331;
//  mr[0]=(float) 0.430574*z; mr[1]=(float) 0.222015*z; mr[2]=(float) 0.020183*z;
//  mg[0]=(float) 0.341550*z; mg[1]=(float) 0.706655*z; mg[2]=(float) 0.129553*z;
//  mb[0]=(float) 0.178325*z; mb[1]=(float) 0.071330*z; mb[2]=(float) 0.939180*z;
//  float maxi=(float) 1.0/270; minu=-88*maxi; minv=-134*maxi;
//  // build (padded) lookup table for y->l conversion assuming y in [0,1]
//  static float lTable[1064]; static bool lInit=false;
//  if( lInit ) return lTable; float y, l;
//  for(int i=0; i<1025; i++) {
//    y = (float) (i/1024.0);
//    l = y>y0 ? 116*(float)pow((double)y,1.0/3.0)-16 : y*a;
//    lTable[i] = l*maxi;
//  }
//  for(int i=1025; i<1064; i++) lTable[i]=lTable[i-1];
//  lInit = true; return lTable;
//}

// Convert from rgb to luv
void rgb2luv( unsigned char *I, float *J, int n, float nrm ) {
  float minu, minv, un, vn, mr[3], mg[3], mb[3];
  float *lTable = rgb2luv_setup(nrm,mr,mg,mb,minu,minv,un,vn);
  float *L=J, *U=L+n, *V=U+n; unsigned char *R=I, *G=R+n, *B=G+n;
  for( int i=0; i<n; i++ ) {
    float r, g, b, x, y, z, l;
    r=(float)*R++; g=(float)*G++; b=(float)*B++;
    x = mr[0]*r + mg[0]*g + mb[0]*b;
    y = mr[1]*r + mg[1]*g + mb[1]*b;
    z = mr[2]*r + mg[2]*g + mb[2]*b;
    l = lTable[(int)(y*1024)];
    *(L++) = l; z = 1/(x + 15*y + 3*z + (float)1e-35);
    *(U++) = l * (13*4*x*z - 13*un) - minu;
    *(V++) = l * (13*9*y*z - 13*vn) - minv;
  }
}


// Convert from rgb to hsv
void rgb2hsv( unsigned char *I, float *J, int n, float nrm ) {
  float *H=J, *S=H+n, *V=S+n;
  unsigned char *R=I, *G=R+n, *B=G+n;
  for(int i=0; i<n; i++) {
    const float r=(float)*(R++), g=(float)*(G++), b=(float)*(B++);
    float h, s, v, minv, maxv;
    if( r==g && g==b ) {
      *(H++) = 0; *(S++) = 0; *(V++) = r*nrm; continue;
    } else if( r>=g && r>=b ) {
      maxv = r; minv = g<b ? g : b;
      h = (g-b)/(maxv-minv)+6; if(h>=6) h-=6;
    } else if( g>=r && g>=b ) {
      maxv = g; minv = r<b ? r : b;
      h = (b-r)/(maxv-minv)+2;
    } else {
      maxv = b; minv = r<g ? r : g;
      h = (r-g)/(maxv-minv)+4;
    }
    h*=(float) (1/6.0); s=1-minv/maxv; v=maxv*nrm;
    *(H++) = h; *(S++) = s; *(V++) = v;
  }
}

// Convert from rgb to gray
void rgb2gray( unsigned char *I, float *J, int n, float nrm ) {
  float *GR=J; unsigned char *R=I, *G=R+n, *B=G+n; int i;
  float mr=(float).2989360213*nrm, mg=(float).5870430745*nrm, mb=(float).1140209043*nrm;
  for(i=0; i<n; i++) *(GR++)=(float)*(R++)*mr + (float)*(G++)*mg + (float)*(B++)*mb;
}

//// Convert from rgb (double) to gray (float)
//void rgb2gray( double *I, float *J, int n, float nrm ) {
//  float *GR=J; double *R=I, *G=R+n, *B=G+n; int i;
//  double mr=.2989360213*nrm, mg=.5870430745*nrm, mb=.1140209043*nrm;
//  for(i=0; i<n; i++) *(GR++) = (float) (*(R++)*mr + *(G++)*mg + *(B++)*mb);
//}

// Copy and normalize only
void normalize( unsigned char *I, float *J, int n, float nrm ) {
  for(int i=0; i<n; i++) *(J++)=(float)*(I++)*nrm;
}

// Convert rgb to various colorspaces
float* rgbConvert( unsigned char *I, int n, int d, int flag, float nrm ) {
  float *J = (float*) malloc(n*(flag==0 ? (d==1?1:d/3) : d)*sizeof(float));
  int i, n1=d*(n<1000?n/10:100); float thr = float(1.001);
  if(flag>1 && nrm==1) for(i=0; i<n1; i++) if(I[i]>thr)
    std::cerr<<"For floats all values in I must be smaller than 1.";
  bool useSse = 0;//n%4==0 && typeid(float)==typeid(float);
  //if( flag==2 && useSse )
  //  for(i=0; i<d/3; i++) rgb2luv_sse(I+i*n*3,(float*)(J+i*n*3),n,(float)nrm);
  if( (flag==0 && d==1) || flag==1 ) normalize(I,J,n*d,nrm);
  else if( flag==0 ) for(i=0; i<d/3; i++) rgb2gray(I+i*n*3,J+i*n*1,n,nrm);
  else if( flag==2 ) for(i=0; i<d/3; i++) rgb2luv(I+i*n*3,J+i*n*3,n,nrm);
  else if( flag==3 ) for(i=0; i<d/3; i++) rgb2hsv(I+i*n*3,J+i*n*3,n,nrm);
  else std::cerr<<"Unknown flag.";
  return J;
}




float* rgb2luv_setup( float z, float *mr, float *mg, float *mb,
  float &minu, float &minv, float &un, float &vn )
{
  // set constants for conversion
  const float y0=(float) ((6.0/29)*(6.0/29)*(6.0/29));
  const float a= (float) ((29.0/3)*(29.0/3)*(29.0/3));
  un=(float) 0.197833; vn=(float) 0.468331;
  mr[0]=(float) 0.430574*z; mr[1]=(float) 0.222015*z; mr[2]=(float) 0.020183*z;
  mg[0]=(float) 0.341550*z; mg[1]=(float) 0.706655*z; mg[2]=(float) 0.129553*z;
  mb[0]=(float) 0.178325*z; mb[1]=(float) 0.071330*z; mb[2]=(float) 0.939180*z;
  float maxi=(float) 1.0/270; minu=-88*maxi; minv=-134*maxi;
  // build (padded) lookup table for y->l conversion assuming y in [0,1]
  static float lTable[1064]; static bool lInfloat=false;
  if( lInfloat ) return lTable; float y, l;
  for(int i=0; i<1025; i++) {
    y = (float) (i/1024.0);
    l = y>y0 ? 116*(float)pow((double)y,1.0/3.0)-16 : y*a;
    lTable[i] = l*maxi;
  }
  for(int i=1025; i<1064; i++) lTable[i]=lTable[i-1];
  lInfloat = true; return lTable;
}

// Convert from rgb to luv
void rgb2luv( float *I, float *J, int n, float nrm ) {
  float minu, minv, un, vn, mr[3], mg[3], mb[3];
  float *lTable = rgb2luv_setup(nrm,mr,mg,mb,minu,minv,un,vn);
  float *L=J, *U=L+n, *V=U+n; float *R=I, *G=R+n, *B=G+n;
  for( int i=0; i<n; i++ ) {
    float r, g, b, x, y, z, l;
    r=(float)*R++; g=(float)*G++; b=(float)*B++;
    x = mr[0]*r + mg[0]*g + mb[0]*b;
    y = mr[1]*r + mg[1]*g + mb[1]*b;
    z = mr[2]*r + mg[2]*g + mb[2]*b;
    l = lTable[(int)(y*1024)];
    *(L++) = l; z = 1/(x + 15*y + 3*z + (float)1e-35);
    *(U++) = l * (13*4*x*z - 13*un) - minu;
    *(V++) = l * (13*9*y*z - 13*vn) - minv;
  }
}

// Convert from rgb to luv using sse
void rgb2luv_sse( float *I, float *J, int n, float nrm ) {
  const int k=256; float R[k], G[k], B[k];
  if( (size_t(R)&15||size_t(G)&15||size_t(B)&15||size_t(I)&15||size_t(J)&15)
    || n%4>0 ) { rgb2luv(I,J,n,nrm); return; }
  int i=0, i1, n1; float minu, minv, un, vn, mr[3], mg[3], mb[3];
  float *lTable = rgb2luv_setup(nrm,mr,mg,mb,minu,minv,un,vn);
  while( i<n ) {
    n1 = i+k; if(n1>n) n1=n; float *J1=J+i; float *R1, *G1, *B1;
    // convert to floats (and load input into cache)
    if( typeid(float) != typeid(float) ) {
      R1=R; G1=G; B1=B; float *Ri=I+i, *Gi=Ri+n, *Bi=Gi+n;
      for( i1=0; i1<(n1-i); i1++ ) {
        R1[i1] = (float) *Ri++; G1[i1] = (float) *Gi++; B1[i1] = (float) *Bi++;
      }
    } else { R1=((float*)I)+i; G1=R1+n; B1=G1+n; }
    // compute RGB -> XYZ
    for( int j=0; j<3; j++ ) {
      __m128 _mr, _mg, _mb, *_J=(__m128*) (J1+j*n);
      __m128 *_R=(__m128*) R1, *_G=(__m128*) G1, *_B=(__m128*) B1;
      _mr=SET(mr[j]); _mg=SET(mg[j]); _mb=SET(mb[j]);
      for( i1=i; i1<n1; i1+=4 ) *(_J++) = ADD( ADD(MUL(*(_R++),_mr),
        MUL(*(_G++),_mg)),MUL(*(_B++),_mb));
    }
    { // compute XZY -> LUV (wfloathout doing L lookup/normalization)
      __m128 _c15, _c3, _cEps, _c52, _c117, _c1024, _cun, _cvn;
      _c15=SET(15.0f); _c3=SET(3.0f); _cEps=SET(1e-35f);
      _c52=SET(52.0f); _c117=SET(117.0f), _c1024=SET(1024.0f);
      _cun=SET(13*un); _cvn=SET(13*vn);
      __m128 *_X, *_Y, *_Z, _x, _y, _z;
      _X=(__m128*) J1; _Y=(__m128*) (J1+n); _Z=(__m128*) (J1+2*n);
      for( i1=i; i1<n1; i1+=4 ) {
        _x = *_X; _y=*_Y; _z=*_Z;
        _z = RCP(ADD(_x,ADD(_cEps,ADD(MUL(_c15,_y),MUL(_c3,_z)))));
        *(_X++) = MUL(_c1024,_y);
        *(_Y++) = SUB(MUL(MUL(_c52,_x),_z),_cun);
        *(_Z++) = SUB(MUL(MUL(_c117,_y),_z),_cvn);
      }
    }
    { // perform lookup for L and finalize computation of U and V
      for( i1=i; i1<n1; i1++ ) J[i1] = lTable[(int)J[i1]];
      __m128 *_L, *_U, *_V, _l, _cminu, _cminv;
      _L=(__m128*) J1; _U=(__m128*) (J1+n); _V=(__m128*) (J1+2*n);
      _cminu=SET(minu); _cminv=SET(minv);
      for( i1=i; i1<n1; i1+=4 ) {
        _l = *(_L++);
        *_U = SUB(MUL(_l,*_U),_cminu); _U++;
        *_V = SUB(MUL(_l,*_V),_cminv); _V++;
      }
    }
    i = n1;
  }
}

// Convert from rgb to hsv
void rgb2hsv( float *I, float *J, int n, float nrm ) {
  float *H=J, *S=H+n, *V=S+n;
  float *R=I, *G=R+n, *B=G+n;
  for(int i=0; i<n; i++) {
    const float r=(float)*(R++), g=(float)*(G++), b=(float)*(B++);
    float h, s, v, minv, maxv;
    if( r==g && g==b ) {
      *(H++) = 0; *(S++) = 0; *(V++) = r*nrm; continue;
    } else if( r>=g && r>=b ) {
      maxv = r; minv = g<b ? g : b;
      h = (g-b)/(maxv-minv)+6; if(h>=6) h-=6;
    } else if( g>=r && g>=b ) {
      maxv = g; minv = r<b ? r : b;
      h = (b-r)/(maxv-minv)+2;
    } else {
      maxv = b; minv = r<g ? r : g;
      h = (r-g)/(maxv-minv)+4;
    }
    h*=(float) (1/6.0); s=1-minv/maxv; v=maxv*nrm;
    *(H++) = h; *(S++) = s; *(V++) = v;
  }
}

// Convert from rgb to gray
void rgb2gray( float *I, float *J, int n, float nrm ) {
  float *GR=J; float *R=I, *G=R+n, *B=G+n; int i;
  float mr=(float).2989360213*nrm, mg=(float).5870430745*nrm, mb=(float).1140209043*nrm;
  for(i=0; i<n; i++) *(GR++)=(float)*(R++)*mr + (float)*(G++)*mg + (float)*(B++)*mb;
}

// Convert from rgb (double) to gray (float)
void rgb2gray( double *I, float *J, int n, float nrm ) {
  float *GR=J; double *R=I, *G=R+n, *B=G+n; int i;
  double mr=.2989360213*nrm, mg=.5870430745*nrm, mb=.1140209043*nrm;
  for(i=0; i<n; i++) *(GR++) = (float) (*(R++)*mr + *(G++)*mg + *(B++)*mb);
}

// Copy and normalize only
void normalize( float *I, float *J, int n, float nrm ) {
  for(int i=0; i<n; i++) *(J++)=(float)*(I++)*nrm;
}

// Convert rgb to various colorspaces
float* rgbConvert( float *I, int n, int d, int flag, float nrm ) {
  float *J = (float*) wrMalloc(n*(flag==0 ? (d==1?1:d/3) : d)*sizeof(float));
 // std::cout<<"n="<<n<<",d="<<d<<std::endl;
  int i, n1=d*(n<1000?n/10:100); float thr = float(1.001);
  if(flag>1 && nrm==1) for(i=0; i<n1; i++) if(I[i]>thr)
    wrError("For floats all values in I must be smaller than 1.");
  bool useSse = n%4==0 && typeid(float)==typeid(float);
  useSse = false;
  if( flag==2 && useSse )
    for(i=0; i<d/3; i++) rgb2luv_sse(I+i*n*3,(float*)(J+i*n*3),n,(float)nrm);
  else if( (flag==0 && d==1) || flag==1 ) normalize(I,J,n*d,nrm);
  else if( flag==0 ) for(i=0; i<d/3; i++) rgb2gray(I+i*n*3,J+i*n*1,n,nrm);
  else if( flag==2 ) for(i=0; i<d/3; i++) rgb2luv(I+i*n*3,J+i*n*3,n,nrm); 
  else if( flag==3 ) for(i=0; i<d/3; i++) rgb2hsv(I+i*n*3,J+i*n*3,n,nrm);
  else wrError("Unknown flag.");
  return J;
}

// compute interpolation values for single column for resapling
void resampleCoef( int ha, int hb, int &n, int *&yas,
  int *&ybs, float *&wts, int bd[2], int pad=0 )
{
  const float s = float(hb)/float(ha), sInv = 1/s; float wt, wt0=float(1e-3)*s;
  bool ds=ha>hb; int nMax; bd[0]=bd[1]=0;
  if(ds) { n=0; nMax=ha+(pad>2 ? pad : 2)*hb; } else { n=nMax=hb; }
  // initialize memory
  wts = (float*)alMalloc(nMax*sizeof(float),16);
  yas = (int*)alMalloc(nMax*sizeof(int),16);
  ybs = (int*)alMalloc(nMax*sizeof(int),16);
  if( ds ) for( int yb=0; yb<hb; yb++ ) {
    // create coefficients for downsampling
    float ya0f=yb*sInv, ya1f=ya0f+sInv, W=0;
    int ya0=int(ceil(ya0f)), ya1=int(ya1f), n1=0;
    for( int ya=ya0-1; ya<ya1+1; ya++ ) {
      wt=s; if(ya==ya0-1) wt=(ya0-ya0f)*s; else if(ya==ya1) wt=(ya1f-ya1)*s;
      if(wt>wt0 && ya>=0) { ybs[n]=yb; yas[n]=ya; wts[n]=wt; n++; n1++; W+=wt; }
    }
    if(W>1) for( int i=0; i<n1; i++ ) wts[n-n1+i]/=W;
    if(n1>bd[0]) bd[0]=n1;
    while( n1<pad ) { ybs[n]=yb; yas[n]=yas[n-1]; wts[n]=0; n++; n1++; }
  } else for( int yb=0; yb<hb; yb++ ) {
    // create coefficients for upsampling
    float yaf = (float(.5)+yb)*sInv-float(.5); int ya=(int) floor(yaf);
    wt=1; if(ya>=0 && ya<ha-1) wt=1-(yaf-ya);
    if(ya<0) { ya=0; bd[0]++; } if(ya>=ha-1) { ya=ha-1; bd[1]++; }
    ybs[yb]=yb; yas[yb]=ya; wts[yb]=wt;
  }
}

// resample A using bilinear interpolation and and store result in B
void imResample( float *A, float *B, int ha, int hb, int wa, int wb, int d, float r ) {
	//	printf("ha=%d, hb=%d, wa=%d, wb=%d\n",ha,hb,wa,wb);
  int hn, wn, x, x1, y, z, xa, xb, ya; float *A0, *A1, *A2, *A3, *B0, wt, wt1;
  float *C = (float*) alMalloc((ha+4)*sizeof(float),16); for(y=ha; y<ha+4; y++) C[y]=0;
  bool sse = (typeid(float)==typeid(float)) && !(size_t(A)&15) && !(size_t(B)&15);
  // get coefficients for resampling along w and h
  int *xas, *xbs, *yas, *ybs; float *xwts, *ywts; int xbd[2], ybd[2];
  resampleCoef( wa, wb, wn, xas, xbs, xwts, xbd, 0 );
  resampleCoef( ha, hb, hn, yas, ybs, ywts, ybd, 4 );
 // printf("resampling coefficients: wa=%d, wb=%d, wn=%d, xas=%d, xbs=%d, xwts=%f, xbd[0]=%d, xbd[1]=%d\n", wa,wb,wn,*xas,*xbs,*xwts,xbd[0],xbd[1]);
 // printf("resampling coefficients: ha=%d, hb=%d, hn=%d, yas=%d, ybs=%d, ywts=%f, ybd[0]=%d, ybd[1]=%d\n", ha,hb,hn,*yas,*ybs,*ywts,ybd[0],ybd[1]);
  if( wa==2*wb ) r/=2; if( wa==3*wb ) r/=3; if( wa==4*wb ) r/=4;
  r/=float(1+1e-6); for( y=0; y<hn; y++ ) ywts[y] *= r;
  // resample each channel in turn
  for( z=0; z<d; z++ ) for( x=0; x<wb; x++ ) {
    if(x==0) x1=0; xa=xas[x1]; xb=xbs[x1]; wt=xwts[x1]; wt1=1-wt; y=0;
    A0=A+z*ha*wa+xa*ha; A1=A0+ha, A2=A1+ha, A3=A2+ha; B0=B+z*hb*wb+xb*hb;
    // variables for SSE (simple casts to float)
    float *Af0, *Af1, *Af2, *Af3, *Bf0, *Cf, *ywtsf, wtf, wt1f;
    Af0=(float*) A0; Af1=(float*) A1; Af2=(float*) A2; Af3=(float*) A3;
    Bf0=(float*) B0; Cf=(float*) C;
    ywtsf=(float*) ywts; wtf=(float) wt; wt1f=(float) wt1;
    // resample along x direction (A -> C)
  //  #define FORs(X) if(sse) for(; y<ha-4; y+=4) STR(Cf[y],X);
    #define FORr(X) for(; y<ha; y++) C[y] = X;
    if( wa==2*wb ) {
    //  FORs( ADD(LDu(Af0[y]),LDu(Af1[y])) );
      FORr( A0[y]+A1[y] ); x1+=2;
    } else if( wa==3*wb ) {
     // FORs( ADD(LDu(Af0[y]),LDu(Af1[y]),LDu(Af2[y])) );
      FORr( A0[y]+A1[y]+A2[y] ); x1+=3;
    } else if( wa==4*wb ) {
   //   FORs( ADD(LDu(Af0[y]),LDu(Af1[y]),LDu(Af2[y]),LDu(Af3[y])) );
      FORr( A0[y]+A1[y]+A2[y]+A3[y] ); x1+=4;
    } else if( wa>wb ) {
      int m=1; while( x1+m<wn && xb==xbs[x1+m] ) m++; float wtsf[4];
      for( int x0=0; x0<(m<4?m:4); x0++ ) wtsf[x0]=float(xwts[x1+x0]);
      #define U(x) MUL( LDu(*(Af ## x + y)), SET(wtsf[x]) )
      #define V(x) *(A ## x + y) * xwts[x1+x]
      if(m==1) { 
		  //FORs(U(0));                    
		  FORr(V(0)); }
      if(m==2) { 
		  //FORs(ADD(U(0),U(1)));        
		  FORr(V(0)+V(1)); }
      if(m==3) { 
		  //FORs(ADD(U(0),U(1),U(2)));
		  FORr(V(0)+V(1)+V(2)); }
      if(m>=4) { 
		  //FORs(ADD(U(0),U(1),U(2),U(3))); 
		  FORr(V(0)+V(1)+V(2)+V(3)); }
      #undef U
      #undef V
      for( int x0=4; x0<m; x0++ ) {
        A1=A0+x0*ha; wt1=xwts[x1+x0]; Af1=(float*) A1; wt1f=float(wt1); y=0;
        //FORs(ADD(LD(Cf[y]),MUL(LDu(Af1[y]),SET(wt1f)))); 
		FORr(C[y]+A1[y]*wt1);
      }
      x1+=m;
    } else {
      bool xBd = x<xbd[0] || x>=wb-xbd[1]; x1++;
      if(xBd) memcpy(C,A0,ha*sizeof(float));
    //  if(!xBd) FORs(ADD(MUL(LDu(Af0[y]),SET(wtf)),MUL(LDu(Af1[y]),SET(wt1f))));
      if(!xBd) FORr( A0[y]*wt + A1[y]*wt1 );
    }
    #undef FORs
    #undef FORr
    // resample along y direction (B -> C)
    if( ha==hb*2 ) {
      float r2 = r/2; int k=((~((size_t) B0) + 1) & 15)/4; y=0;
      for( ; y<k; y++ )  B0[y]=(C[2*y]+C[2*y+1])*r2;
   //   if(sse) for(; y<hb-4; y+=4) STR(Bf0[y],MUL((float)r2,_mm_shuffle_ps(ADD(
   //     LDu(Cf[2*y]),LDu(Cf[2*y+1])),ADD(LDu(Cf[2*y+4]),LDu(Cf[2*y+5])),136)));
      for( ; y<hb; y++ ) B0[y]=(C[2*y]+C[2*y+1])*r2;
    } else if( ha==hb*3 ) {
      for(y=0; y<hb; y++) B0[y]=(C[3*y]+C[3*y+1]+C[3*y+2])*(r/3);
    } else if( ha==hb*4 ) {
      for(y=0; y<hb; y++) B0[y]=(C[4*y]+C[4*y+1]+C[4*y+2]+C[4*y+3])*(r/4);
    } else if( ha>hb ) {
      y=0;
      //if( sse && ybd[0]<=4 ) for(; y<hb; y++) // Requires SSE4
      //  STR1(Bf0[y],_mm_dp_ps(LDu(Cf[yas[y*4]]),LDu(ywtsf[y*4]),0xF1));
      #define U(o) C[ya+o]*ywts[y*4+o]
      if(ybd[0]==2) for(; y<hb; y++) { ya=yas[y*4]; B0[y]=U(0)+U(1); }
      if(ybd[0]==3) for(; y<hb; y++) { ya=yas[y*4]; B0[y]=U(0)+U(1)+U(2); }
      if(ybd[0]==4) for(; y<hb; y++) { ya=yas[y*4]; B0[y]=U(0)+U(1)+U(2)+U(3); }
      if(ybd[0]>4)  for(; y<hn; y++) { B0[ybs[y]] += C[yas[y]] * ywts[y]; }
      #undef U
    } else {
      for(y=0; y<ybd[0]; y++) B0[y] = C[yas[y]]*ywts[y];
      for(; y<hb-ybd[1]; y++) B0[y] = C[yas[y]]*ywts[y]+C[yas[y]+1]*(r-ywts[y]);
      for(; y<hb; y++)        B0[y] = C[yas[y]]*ywts[y];
    }
  }
  alFree(xas); alFree(xbs); alFree(xwts); alFree(C);
  alFree(yas); alFree(ybs); alFree(ywts);

}

// compute x and y gradients for just one column (uses sse)
void grad1( float *I, float *Gx, float *Gy, int h, int w, int x ) {
  int y, y1; float *Ip, *In, r; __m128 *_Ip, *_In, *_G, _r;
  // compute column of Gx
  Ip=I-h; In=I+h; r=.5f;
  if(x==0) { r=1; Ip+=h; } else if(x==w-1) { r=1; In-=h; }
  if( h<4 || h%4>0 || (size_t(I)&15) || (size_t(Gx)&15) ) {
    for( y=0; y<h; y++ ) 
		*Gx++=(*In++-*Ip++)*r;
  } else {
    _G=(__m128*) Gx; _Ip=(__m128*) Ip; _In=(__m128*) In; _r = SET(r);
    for(y=0; y<h; y+=4) *_G++=MUL(SUB(*_In++,*_Ip++),_r);
  }
  // compute column of Gy
  #define GRADY(r) *Gy++=(*In++-*Ip++)*r;
  Ip=I; In=Ip+1;
  // GRADY(1); Ip--; for(y=1; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
  y1=((~((size_t) Gy) + 1) & 15)/4; 
  if(y1==0) y1=4; if(y1>h-1) y1=h-1;
  GRADY(1); Ip--; for(y=1; y<y1; y++) GRADY(.5f);
  _r = SET(.5f); _G=(__m128*) Gy;
  for(; y+4<h-1; y+=4, Ip+=4, In+=4, Gy+=4)
    *_G++=MUL(SUB(LDu(*In),LDu(*Ip)),_r);
  for(; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
  #undef GRADY
}

// compute x and y gradients at each location (uses sse)
void grad2( float *I, float *Gx, float *Gy, int h, int w, int d ) {
  int o, x, c, a=w*h; for(c=0; c<d; c++) for(x=0; x<w; x++) {
    o=c*a+x*h; grad1( I+o, Gx+o, Gy+o, h, w, x );
  }
}

float* acosTable() {
  const int n=10000, b=10; int i;
  static float a[n*2+b*2]; static bool init=false;
  float *a1=a+n+b; if( init ) return a1;
  for( i=-n-b; i<-n; i++ )   a1[i]=PI;
  for( i=-n; i<n; i++ )      a1[i]=float(acos(i/float(n)));
  for( i=n; i<n+b; i++ )     a1[i]=0;
  for( i=-n-b; i<n/10; i++ ) if( a1[i] > PI-1e-6f ) a1[i]=PI-1e-6f;
  init=true; return a1;
}

void gradMag( float *I, float *M, float *O, int h, int w, int d, bool full ) {
  int x, y, y1, c, h4, s; float *Gx, *Gy, *M2; __m128 *_Gx, *_Gy, *_M2, _m;
  float *acost = acosTable(), acMult=10000.0f;
  // allocate memory for storing one column of output (padded so h4%4==0)
  h4=(h%4==0) ? h : h-(h%4)+4; s=d*h4*sizeof(float);
  M2=(float*) alMalloc(s,16); _M2=(__m128*) M2;
  Gx=(float*) alMalloc(s,16); _Gx=(__m128*) Gx;
  Gy=(float*) alMalloc(s,16); _Gy=(__m128*) Gy;
  // compute gradient magnitude and orientation for each column
  for( x=0; x<w; x++ ) {
    // compute gradients (Gx, Gy) with maximum squared magnitude (M2)
    for(c=0; c<d; c++) {
      grad1( I+x*h+c*w*h, Gx+c*h4, Gy+c*h4, h, w, x );
     // if(debug_count ==0){
	//	  std::cout << "I am Gx" <<std::endl;
	//	  std::cout << *(Gx+c*h4) <<std::endl;
	//	  debug_count++;
	  //}
	  for( y=0; y<h4/4; y++ ) {
        y1=h4/4*c+y;
        _M2[y1]=ADD(MUL(_Gx[y1],_Gx[y1]),MUL(_Gy[y1],_Gy[y1]));
        if( c==0 ) continue; _m = CMPGT( _M2[y1], _M2[y] );
        _M2[y] = OR( AND(_m,_M2[y1]), ANDNOT(_m,_M2[y]) );
        _Gx[y] = OR( AND(_m,_Gx[y1]), ANDNOT(_m,_Gx[y]) );
        _Gy[y] = OR( AND(_m,_Gy[y1]), ANDNOT(_m,_Gy[y]) );
      }
    }
//	if(debug_count ==0){
//		for(int i=0;i<10;i++)
//			std::cout << Gx[i] << ",";
//		std::cout << std::endl;
//		debug_count ++;
//	}
    // compute gradient mangitude (M) and normalize Gx
    for( y=0; y<h4/4; y++ ) {
      _m = __MIN( RCPSQRT(_M2[y]), SET(1e10f) );
//		if(debug_count == 0){
//			std::cout << M2[0] <<std::endl;
//		}
//	  _M2[y] = _m);
	  _M2[y] = RCP(_m);
	//	if(debug_count == 0){
	//		std::cout << M2[0] <<std::endl;
	//	}
		
      if(O) _Gx[y] = MUL( MUL(_Gx[y],_m), SET(acMult) );
/*	if(debug_count == 0){
		for(int i=0;i<10;i++)
			std::cout << Gx[i] <<",";
		std::cout << std::endl;
		for(int i=0;i<10;i++)
			std::cout << Gy[i] <<",";
		std::cout << std::endl;
		debug_count ++;
	}
  */    if(O) _Gx[y] = XOR( _Gx[y], AND(_Gy[y], SET(-0.f)) );
	}
/*	if(debug_count == 1){
		for(int i=0;i<10;i++)
			std::cout << Gx[i] <<",";
		std::cout << std::endl;
		for(int i=0;i<10;i++)
			std::cout << Gy[i] <<",";
		std::cout << std::endl;
		debug_count ++;
	}*/
    memcpy( M+x*h, M2, h*sizeof(float) );
    // compute and store gradient orientation (O) via table lookup
    if( O!=0 ) for( y=0; y<h; y++ ) O[x*h+y] = acost[(int)Gx[y]];
    if( O!=0 && full ) {
      y1=((~size_t(O+x*h)+1)&15)/4; y=0;
      for( ; y<y1; y++ ) O[y+x*h]+=(Gy[y]<0)*PI;
      for( ; y<h-4; y+=4 ) STRu( O[y+x*h],
        ADD( LDu(O[y+x*h]), AND(CMPLT(LDu(Gy[y]),SET(0.f)),SET(PI)) ) );
      for( ; y<h; y++ ) O[y+x*h]+=(Gy[y]<0)*PI;
    }
  }
  alFree(Gx); alFree(Gy); alFree(M2);
}

// convolve one column of I by a 2rx1 triangle filter
void convTriY( float *I, float *O, int h, int r, int s ) {
  r++; float t, u; int j, r0=r-1, r1=r+1, r2=2*h-r, h0=r+1, h1=h-r+1, h2=h;
  u=t=I[0]; for( j=1; j<r; j++ ) u+=t+=I[j]; u=2*u-t; t=0;
  if( s==1 ) {
    O[0]=u; j=1;
    for(; j<h0; j++) O[j] = u += t += I[r-j]  + I[r0+j] - 2*I[j-1];
    for(; j<h1; j++) O[j] = u += t += I[j-r1] + I[r0+j] - 2*I[j-1];
    for(; j<h2; j++) O[j] = u += t += I[j-r1] + I[r2-j] - 2*I[j-1];
  } else {
    int k=(s-1)/2; h2=(h/s)*s; if(h0>h2) h0=h2; if(h1>h2) h1=h2;
    if(++k==s) { k=0; *O++=u; } j=1;
    for(;j<h0;j++) { u+=t+=I[r-j] +I[r0+j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
    for(;j<h1;j++) { u+=t+=I[j-r1]+I[r0+j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
    for(;j<h2;j++) { u+=t+=I[j-r1]+I[r2-j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
  }
}

// convolve I by a 2rx1 triangle filter (uses SSE)
void convTri( float *I, float *O, int h, int w, int d, int r, int s ) {
  r++; float nrm = 1.0f/(r*r*r*r); int i, j, k=(s-1)/2, h0, h1, w0;
  if(h%4==0) h0=h1=h; else { h0=h-(h%4); h1=h0+4; } w0=(w/s)*s;
  float *T=(float*) alMalloc(2*h1*sizeof(float),16), *U=T+h1;
  while(d-- > 0) {
    // initialize T and U
    for(j=0; j<h0; j+=4) STR(U[j], STR(T[j], LDu(I[j])));
    for(i=1; i<r; i++) for(j=0; j<h0; j+=4) INC(U[j],INC(T[j],LDu(I[j+i*h])));
    for(j=0; j<h0; j+=4) STR(U[j],MUL(nrm,(SUB(MUL(2,LD(U[j])),LD(T[j])))));
    for(j=0; j<h0; j+=4) STR(T[j],0);
    for(j=h0; j<h; j++ ) U[j]=T[j]=I[j];
    for(i=1; i<r; i++) for(j=h0; j<h; j++ ) U[j]+=T[j]+=I[j+i*h];
    for(j=h0; j<h; j++ ) { U[j] = nrm * (2*U[j]-T[j]); T[j]=0; }
    // prepare and convolve each column in turn
    k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    for( i=1; i<w0; i++ ) {
      float *Il=I+(i-1-r)*h; if(i<=r) Il=I+(r-i)*h; float *Im=I+(i-1)*h;
      float *Ir=I+(i-1+r)*h; if(i>w-r) Ir=I+(2*w-r-i)*h;
      for( j=0; j<h0; j+=4 ) {
        INC(T[j],ADD(LDu(Il[j]),LDu(Ir[j]),MUL(-2,LDu(Im[j]))));
        INC(U[j],MUL(nrm,LD(T[j])));
      }
      for( j=h0; j<h; j++ ) U[j]+=nrm*(T[j]+=Il[j]+Ir[j]-2*Im[j]);
      k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    }
    I+=w*h;
  }
  alFree(T);
}

// normalize gradient magnitude at each location (uses sse)
void gradMagNorm( float *M, float *S, int h, int w, float norm ) {
  __m128 *_M, *_S, _norm; int i=0, n=h*w, n4=n/4;
  _S = (__m128*) S; _M = (__m128*) M; _norm = SET(norm);
  bool sse = !(size_t(M)&15) && !(size_t(S)&15);
  sse=false;	
  // I am not sure whether this change is right. by ttq
  if(sse) for(; i<n4; i++) { *_M=MUL(*_M,RCP(ADD(*_S++,_norm))); _M++; }
  if(!sse) i*=4; for(; i<n; i++) M[i] /= (S[i] + norm);
}

// helper for gradHist, quantize O and M into O0, O1 and M0, M1 (uses sse)
void gradQuantize( float *O, float *M, int *O0, int *O1, float *M0, float *M1,
  int nb, int n, float norm, int nOrients, bool full, bool interpolate )
{
  // assumes all *OUTPUT* matrices are 4-byte aligned
  int i, o0, o1; float o, od, m;
  __m128i _o0, _o1, *_O0, *_O1; __m128 _o, _od, _m, *_M0, *_M1;
  // define useful constants
  const float oMult=(float)nOrients/(full?2*PI:PI); const int oMax=nOrients*nb;
  const __m128 _norm=SET(norm), _oMult=SET(oMult), _nbf=SET((float)nb);
  const __m128i _oMax=SET(oMax), _nb=SET(nb);
  // perform the majority of the work with sse
  _O0=(__m128i*) O0; _O1=(__m128i*) O1; _M0=(__m128*) M0; _M1=(__m128*) M1;
  if( interpolate ) for( i=0; i<=n-4; i+=4 ) {
    _o=MUL(LDu(O[i]),_oMult); _o0=CVT(_o); _od=SUB(_o,CVT(_o0));
    _o0=CVT(MUL(CVT(_o0),_nbf)); _o0=AND(CMPGT(_oMax,_o0),_o0); *_O0++=_o0;
    _o1=ADD(_o0,_nb); _o1=AND(CMPGT(_oMax,_o1),_o1); *_O1++=_o1;
    _m=MUL(LDu(M[i]),_norm); *_M1=MUL(_od,_m); *_M0++=SUB(_m,*_M1); _M1++;
  } else for( i=0; i<=n-4; i+=4 ) {
    _o=MUL(LDu(O[i]),_oMult); _o0=CVT(ADD(_o,SET(.5f)));
    _o0=CVT(MUL(CVT(_o0),_nbf)); _o0=AND(CMPGT(_oMax,_o0),_o0); *_O0++=_o0;
    *_M0++=MUL(LDu(M[i]),_norm); *_M1++=SET(0.f); *_O1++=SET(0);
  }
  // compute trailing locations without sse
  if( interpolate ) for(; i<n; i++ ) {
    o=O[i]*oMult; o0=(int) o; od=o-o0;
    o0*=nb; if(o0>=oMax) o0=0; O0[i]=o0;
    o1=o0+nb; if(o1==oMax) o1=0; O1[i]=o1;
    m=M[i]*norm; M1[i]=od*m; M0[i]=m-M1[i];
  } else for(; i<n; i++ ) {
    o=O[i]*oMult; o0=(int) (o+.5f);
    o0*=nb; if(o0>=oMax) o0=0; O0[i]=o0;
    M0[i]=M[i]*norm; M1[i]=0; O1[i]=0;
  }
}


// compute nOrients gradient histograms per bin x bin block of pixels
void gradHist( float *M, float *O, float *H, int h, int w,
  int bin, int nOrients, int softBin, bool full )
{
//	std::cout << "h=" <<h<<", w=" <<w<<", bin="<<bin<<", nOrients="<<nOrients<<", softBin="<<softBin<<", full="<<full<<std::endl;
  const int hb=h/bin, wb=w/bin, h0=hb*bin, w0=wb*bin, nb=wb*hb;
  const float s=(float)bin, sInv=1/s, sInv2=1/s/s;
  float *H0, *H1, *M0, *M1; int x, y; int *O0, *O1; float xb, init;
  O0=(int*)alMalloc(h*sizeof(int),16); M0=(float*) alMalloc(h*sizeof(float),16);
  O1=(int*)alMalloc(h*sizeof(int),16); M1=(float*) alMalloc(h*sizeof(float),16);
  // main loop
  for( x=0; x<w0; x++ ) {
    // compute target orientation bins for entire column - very fast
    gradQuantize(O+x*h,M+x*h,O0,O1,M0,M1,nb,h0,sInv2,nOrients,full,softBin>=0);

    if( softBin<0 && softBin%2==0 ) {
      // no interpolation w.r.t. either orienation or spatial bin
      H1=H+(x/bin)*hb;
      #define GH H1[O0[y]]+=M0[y]; y++;
      if( bin==1 )      for(y=0; y<h0;) { GH; H1++; }
      else if( bin==2 ) 
		  for(y=0; y<h0;) 
		  { GH; GH; H1++; }
      else if( bin==3 ) for(y=0; y<h0;) { GH; GH; GH; H1++; }
      else if( bin==4 ) for(y=0; y<h0;) { GH; GH; GH; GH; H1++; }
      else for( y=0; y<h0;) { for( int y1=0; y1<bin; y1++ ) { GH; } H1++; }
      #undef GH

    } else if( softBin%2==0 || bin==1 ) {
      // interpolate w.r.t. orientation only, not spatial bin
      H1=H+(x/bin)*hb;
      #define GH H1[O0[y]]+=M0[y]; H1[O1[y]]+=M1[y]; y++;
      if( bin==1 )      for(y=0; y<h0;) { GH; H1++; }
      else if( bin==2 ) 
		  for(y=0; y<h0;) 
		  { GH; GH; H1++; }
      else if( bin==3 ) for(y=0; y<h0;) { GH; GH; GH; H1++; }
      else if( bin==4 ) for(y=0; y<h0;) { GH; GH; GH; GH; H1++; }
      else for( y=0; y<h0;) { for( int y1=0; y1<bin; y1++ ) { GH; } H1++; }
      #undef GH

    } else {
      // interpolate using trilinear interpolation
      float ms[4], xyd, yb, xd, yd; __m128 _m, _m0, _m1;
      bool hasLf, hasRt; int xb0, yb0;
      if( x==0 ) { init=(0+.5f)*sInv-0.5f; xb=init; }
      hasLf = xb>=0; xb0 = hasLf?(int)xb:-1; hasRt = xb0 < wb-1;
      xd=xb-xb0; xb+=sInv; yb=init; y=0;
      // macros for code conciseness
      #define GHinit yd=yb-yb0; yb+=sInv; H0=H+xb0*hb+yb0; xyd=xd*yd; \
        ms[0]=1-xd-yd+xyd; ms[1]=yd-xyd; ms[2]=xd-xyd; ms[3]=xyd;
      #define GH(H,ma,mb) H1=H; STRu(*H1,ADD(LDu(*H1),MUL(ma,mb)));
      // leading rows, no top bin
      for( ; y<bin/2; y++ ) {
        yb0=-1; GHinit;
        if(hasLf) { H0[O0[y]+1]+=ms[1]*M0[y]; H0[O1[y]+1]+=ms[1]*M1[y]; }
        if(hasRt) { H0[O0[y]+hb+1]+=ms[3]*M0[y]; H0[O1[y]+hb+1]+=ms[3]*M1[y]; }
      }
      // main rows, has top and bottom bins, use SSE for minor speedup
      if( softBin<0 ) for( ; ; y++ ) {
        yb0 = (int) yb; if(yb0>=hb-1) break; GHinit; _m0=SET(M0[y]);
        if(hasLf) { _m=SET(0,0,ms[1],ms[0]); GH(H0+O0[y],_m,_m0); }
        if(hasRt) { _m=SET(0,0,ms[3],ms[2]); GH(H0+O0[y]+hb,_m,_m0); }
      } else for( ; ; y++ ) {
        yb0 = (int) yb; if(yb0>=hb-1) break; GHinit;
        _m0=SET(M0[y]); _m1=SET(M1[y]);
        if(hasLf) { _m=SET(0,0,ms[1],ms[0]);
          GH(H0+O0[y],_m,_m0); GH(H0+O1[y],_m,_m1); }
        if(hasRt) { _m=SET(0,0,ms[3],ms[2]);
          GH(H0+O0[y]+hb,_m,_m0); GH(H0+O1[y]+hb,_m,_m1); }
      }
      // final rows, no bottom bin
      for( ; y<h0; y++ ) {
        yb0 = (int) yb; GHinit;
        if(hasLf) { H0[O0[y]]+=ms[0]*M0[y]; H0[O1[y]]+=ms[0]*M1[y]; }
        if(hasRt) { H0[O0[y]+hb]+=ms[2]*M0[y]; H0[O1[y]+hb]+=ms[2]*M1[y]; }
      }
      #undef GHinit
      #undef GH
    }
  }
  alFree(O0); alFree(O1); alFree(M0); alFree(M1);
  // normalize boundary bins which only get 7/8 of weight of interior bins
  if( softBin%2!=0 ) for( int o=0; o<nOrients; o++ ) {
    x=0; for( y=0; y<hb; y++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    y=0; for( x=0; x<wb; x++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    x=wb-1; for( y=0; y<hb; y++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    y=hb-1; for( x=0; x<wb; x++ ) H[o*nb+x*hb+y]*=8.f/7.f;
  }
}

// construct lookup array for mapping fids to channel indices
uint32* buildLookup( int *dims, int w ) {
  int c, r, z, n=w*w*dims[2]; uint32 *cids=new uint32[n]; n=0;
  for(z=0; z<dims[2]; z++) for(c=0; c<w; c++) for(r=0; r<w; r++)
    cids[n++] = z*dims[0]*dims[1] + c*dims[0] + r;
  return cids;
}

// construct lookup arrays for mapping fids for self-similarity channel
void buildLookupSs( uint32 *&cids1, uint32 *&cids2, int *dims, int w, int m ) {
  int i, j, z, z1, c, r; int locs[1024];
  int m2=m*m, n=m2*(m2-1)/2*dims[2], s=int(w/m/2.0+.5);
  cids1 = new uint32[n]; cids2 = new uint32[n]; n=0;
  for(i=0; i<m; i++) locs[i]=uint32((i+1)*(w+2*s-1)/(m+1.0)-s+.5);
  for(z=0; z<dims[2]; z++) for(i=0; i<m2; i++) for(j=i+1; j<m2; j++) {
    z1=z*dims[0]*dims[1]; n++;
    r=i%m; c=(i-r)/m; cids1[n-1]= z1 + locs[c]*dims[0] + locs[r];
    r=j%m; c=(j-r)/m; cids2[n-1]= z1 + locs[c]*dims[0] + locs[r];
  }
}

// return I[x,y] via bilinear interpolation
float interp( float *I, int h, int w, float x, float y ) {
  x = x<0 ? 0 : (x>w-1.001 ? w-1.001 : x);
  y = y<0 ? 0 : (y>h-1.001 ? h-1.001 : y);
  int x0=int(x), y0=int(y), x1=x0+1, y1=y0+1;
  float dx0=x-x0, dy0=y-y0, dx1=1-dx0, dy1=1-dy0;
  return I[x0*h+y0]*dx1*dy1 + I[x1*h+y0]*dx0*dy1 +
    I[x0*h+y1]*dx1*dy0 + I[x1*h+y1]*dx0*dy0;
}
