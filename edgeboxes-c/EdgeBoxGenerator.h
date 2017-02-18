# ifndef __EDGEBOXGENERATOR__

#include "extern.h"
#include "Image.h"
#include <vector>
#include <algorithm>

//int clamp( int v, int a, int b ) { return v<a?a:v>b?b:v; }
int clamp( int v, int a, int b );// { return v<a?a:v>b?b:v; }
// trivial array class encapsulating pointer arrays
template <class T> class Array
{
public:
  Array() { _h=_w=0; _x=0; _free=0; }
  virtual ~Array() { clear(); }
  void clear() { if(_free) delete [] _x; _h=_w=0; _x=0; _free=0; }
  void init(int h, int w) { clear(); _h=h; _w=w; _x=new T[h*w](); _free=1; }
  T& val(size_t c, size_t r) { return _x[c*_h+r]; }
  int _h, _w; T *_x; bool _free;
};

// convenient typedefs
typedef std::vector<float> vectorf;
typedef std::vector<int> vectori;
typedef Array<float> arrayf;
typedef Array<int> arrayi;

typedef struct { int c, r, w, h; float s; } Box;
typedef std::vector<Box> Boxes;
//bool boxesCompare( const Box &a, const Box &b ) { return a.s<b.s; }
bool boxesCompare( const Box &a, const Box &b ); //{ return a.s<b.s; }

float boxesOverlap( Box &a, Box &b );
void boxesNms( Boxes &boxes, float thr, float eta, int maxBoxes );
void edgeBoxesWrapper(Image &E_img, Image &O_img, Boxes & boxes, float alpha, float beta, float eta,
	float minScore,	float maxBoxes, float edgeMinMag,
	float edgeMergeThr, float clusterMinMag, float maxAspectRatio,
	float minBoxArea, float gamma, float kappa);

// main class for generating edge boxes
class EdgeBoxGenerator
{
public:
  // method parameters (must be manually set)
  float _alpha, _beta, _eta, _minScore; int _maxBoxes;
  float _edgeMinMag, _edgeMergeThr, _clusterMinMag;
  float _maxAspectRatio, _minBoxArea, _gamma, _kappa;

  // main external routine (set parameters first)
  void generate( Boxes &boxes, arrayf &E, arrayf &O, arrayf &V );

private:
  // edge segment information (see clusterEdges)
  int h, w;                         // image dimensions
  int _segCnt;                      // total segment count
  arrayi _segIds;                   // segment ids (-1/0 means no segment)
  vectorf _segMag;                  // segment edge magnitude sums
  vectori _segR, _segC;             // segment lower-right pixel
  std::vector<vectorf> _segAff;          // segment affinities
  std::vector<vectori> _segAffIdx;       // segment neighbors

  // data structures for efficiency (see prepDataStructs)
  arrayf _segIImg, _magIImg; arrayi _hIdxImg, _vIdxImg;
  std::vector<vectori> _hIdxs, _vIdxs; vectorf _scaleNorm;
  float _scStep, _arStep, _rcStepRatio;

  // data structures for efficiency (see scoreBox)
  arrayf _sWts; arrayi _sDone, _sMap, _sIds; int _sId;

  // helper routines
  void clusterEdges( arrayf &E, arrayf &O, arrayf &V );
  void prepDataStructs( arrayf &E );
  void scoreAllBoxes( Boxes &boxes );
  void scoreBox( Box &box );
  void refineBox( Box &box );
  void drawBox( Box &box, arrayf &E, arrayf &V );
};

# define __EDGEBOXGENERATOR__
# endif
