//#include "opencv2\opencv.hpp"
#include "Option.h"
//#include "model.h"
#include "time.h"
#include "Image.h"
#include "EdgeBoxGenerator.h"
#include "ttq.h"
//#include <iostream>
#include <cstdio>
//using namespace cv;
int debug_count = 0;
int main(int argc, char* argv[])
{
	
	
	Mat image;
	image = imread(argv[1]);
    int num = 10000;
    if (argc > 2)
    {
        num = atoi(argv[2]);
    }



	//upload the image
	Image myimg(image);
	Image chnsReg, chnsSim;
	Option myopts;
	

	time_t start = clock();	
	//for(int i=0;i<100;i++){
	char* filename = "./models/forest/modelBsds.dat";
	Model mymodel;
	
//	ttqDebug("%s\n", "Here we go");	
	mymodel.load(filename);
	
//	ttqDebug("%s\n", "Here we go");	

	//int pt=myopts.imWidth/2+(4-(myimg.numRows+myopts.imWidth)%4)%4;
	//int pb=myopts.imWidth/2+(4-(myimg.numRows+myopts.imWidth)%4)%4;
	//int pl=myopts.imWidth/2+(4-(myimg.numCols+myopts.imWidth)%4)%4;
	//int pr=myopts.imWidth/2+(4-(myimg.numCols+myopts.imWidth)%4)%4;
	//const int h1 = (int) ceil(double(myimg.numRows+pt+pb-myopts.imWidth)/myopts.stride);
	//const int w1 = (int) ceil(double(myimg.numCols+pl+pr-myopts.imWidth)/myopts.stride);
	//const int h2 = h1*myopts.stride+myopts.gtWidth;
	//const int w2 = w1*myopts.stride+myopts.gtWidth;	
	//printf("in main:h1=%d, w1=%d, h2=%d, w2=%d\n",h1,w1,h2,w2);
	//float *G =  new float [(h2-myopts.gtWidth)*(w2-myopts.gtWidth)];
	Image E, G, E_new;
	

	myopts.edgeDection(myimg, mymodel, E, G);
	G.edgesNmsWrapper(E, E_new, 2,0,1, myopts.nThreads);
	//printf("%f\n", E.imgData[182271-1]);
	//printf("%f\n", G.imgData[182271-1]);
	//printf("%f\n", E_new.imgData[182271-1]);
	//std::cout << E.numRows << "," << E.numCols <<", "<<E.numChannels<<std::endl;
	//std::cout << G.numRows << "," << G.numCols <<", "<<G.numChannels<<std::endl;

	//for(int i=0;i<myopts.nChns;i++){
	//	std::cout << i+1 << ":" <<std::endl;
	//	chnsSim.log(i);
	//	std::cout<<std::endl;
	//}
	//std::istream infile;
	Boxes boxes;
	edgeBoxesWrapper(E_new, G, boxes, myopts.alpha, myopts.beta, myopts.eta,
	myopts.minScore, myopts.maxBoxes, myopts.edgeMinMag,
	myopts.edgeMergeThr, myopts.clusterMinMag, myopts.maxAspectRatio,
	myopts.minBoxArea, myopts.gamma, myopts.kappa);
	
	//FILE* fid_result=fopen(strcat(argv[1], ".txt"), "w");
	for(int i=0;i<boxes.size() && i < num;i++)
    {
		printf("%d %d %d %d\n", boxes[i].c, boxes[i].r, boxes[i].c + boxes[i].w, boxes[i].r + boxes[i].h);
        }
	//fclose(fid_result);
//	time_t end = clock();
//	std::cout <<"time:" <<  double(end -start)/CLOCKS_PER_SEC << "s\n" ;

	
	//imshow("Hello, OpenCV",image);
	//waitKey(-1);
	return 0;
}
