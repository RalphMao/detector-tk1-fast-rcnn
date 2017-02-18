# include "opencv2/opencv.hpp"
# include "extern.h"
# include "model.h"


using namespace cv;
class Option;
class Model;

# ifndef _IMAGE_
# define _IMAGE_
class Image{
public:
	Image(Mat& cv_image);
	Image(){
		this->numChannels = 0;
		this->numCols = 0;
		this->numRows = 0;
	}
	~Image(){	
		if(this->numChannels!=0 && this->numCols!=0 && this->numRows!=0)
			delete []imgData;
	}

	void Set(int numRows, int numCols, int numChannels);

	void myNormalize(int foo);
	
	void imPadWrapper(Image & img_warpped, int  pt, int pb, int pl, int pr, char*type);

	void imChnsWrapper(Image &img,  Option &opts);

	void rgbConvertWrapper(Image &img, string colorspace, bool useSingle);

	void imgResampleWrapper(Image &img_shrink, float alpha, std::string method, float norm);
	
	void gradWrapper(Image &H, Image &M, bool full, int normRad, float normConst, int shrink, int nOrients, int softBin);

	void log(int channel);

//	void edgeDectionExternal(Model &model, Option& opts, float* chns, float*chnsSs);
	void edgesNmsWrapper(Image &E0, Image &E, int r, int s, int m, int nThreads);


	float *imgData;
	int numCols;
	int numRows;
	int numChannels;
};
#endif
