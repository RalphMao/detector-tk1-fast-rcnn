#include <string>
//#include <iostream>
#include "Image.h"
#include "model.h"
class Image;
class Model;

//typedef unsigned int uint32;
//typedef unsigned short uint16;
//typedef unsigned char uint8;



#ifndef _OPTION_
#define _OPTION_
class Option{
public:
	Option():imWidth(32), gtWidth(16),nPos(500000),nNeg(500000),nImgs(0),
		nTrees(8),fracFtrs(0.25),minCount(1),minChild(8),
		maxDepth(64),discretize("pca"), nSample(256),nClasses(2),
		split("gini"),nOrients(4),grdSmooth(0),chnSmooth(2),
		simSmooth(8),normRad(4),shrink(2),nCells(5),rgbd(0),
		stride(2),multiscale(0),sharpen(2),nTreesEval(4),
		nThreads(4),nms(0),seed(1),useParfor(0),modelDir("./models/"),
		modelFnm("model"),bsdsDir("./BSR/BSDS500/data/"),
		alpha(0.65),beta(0.75),eta(1),minScore(0.01),maxBoxes(10000),
		edgeMinMag(0.1), edgeMergeThr(0.5), clusterMinMag(0.5),
		maxAspectRatio(3),minBoxArea(1000),gamma(2),kappa(1.5)
	{
		
		int nChnsGrad = (this->nOrients+1)*2;
		int nChnsColor=3;
		if(this->rgbd == 1)
			nChnsColor =1;
		else if(this->rgbd==2){
			nChnsGrad=nChnsGrad*2; 
			nChnsColor=nChnsColor+1;
		}
		this->nChns = nChnsGrad + nChnsColor;
		this->nChnFtrs = this->imWidth*this->imWidth*this->nChns/this->shrink/this->shrink;
		this->nSimFtrs = (this->nCells*this->nCells)*(this->nCells*this->nCells-1)/2*this->nChns;
		this->nTotFtrs = this->nChnFtrs + this->nSimFtrs;

	}
	~Option() {}
	void edgeDection( Image& img, Model & model, Image &E_out, Image & G);
	void edgeDectionExternal(Model &model, float* img, float* chns, float*chnsSs, int h, int w, int Z,
		float *E, uint32*ind);
	void edgeChns(Image &img_before, Image &img_after);
	

//private:
	int imWidth;
	int gtWidth;
	int nPos;	int nNeg;	//	decrease to speedup training
	int nImgs;
	int nTrees;
	double fracFtrs;
	int minCount;
	int minChild;
	int maxDepth;
	std::string discretize;
	int nSample;
	int nClasses;
	std::string split;
	int nOrients;
	int grdSmooth;
	int chnSmooth;
	int simSmooth;
	int normRad;
	int shrink;
	int nCells;
	int rgbd;
	int stride;
	int multiscale;
	int sharpen;
	int nTreesEval;
	int nThreads;
	int nms;
	int seed;
	int useParfor;	//	parallelize if sufficient memory
	int nChns;
	int nChnFtrs;
	int nSimFtrs;
	int nTotFtrs;

	float alpha;
	float beta;
	float eta;
	float minScore;
	float maxBoxes;
	float edgeMinMag;
	float edgeMergeThr;
	float clusterMinMag;
	float maxAspectRatio;
	float minBoxArea;
	float gamma;
	float kappa;

	std::string modelDir;
	std::string modelFnm;
	std::string bsdsDir;

};
#endif