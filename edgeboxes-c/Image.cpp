# include "Image.h"

Image::Image(Mat& cv_image){
	this->imgData = new float [cv_image.rows*cv_image.cols*3];
	this->numChannels = 3;
	this->numCols = cv_image.cols;
	this->numRows = cv_image.rows;
	for(int i=0; i<cv_image.cols; i++)
		for(int j=0; j<cv_image.rows; j++){
			imgData[i*cv_image.rows+j] = cv_image.at<Vec3b>(j,i)[2];	// 'b'
			imgData[cv_image.rows*cv_image.cols+i*cv_image.rows+j] = cv_image.at<Vec3b>(j,i)[1]; // 'g'
			imgData[2*cv_image.rows*cv_image.cols+i*cv_image.rows+j] = cv_image.at<Vec3b>(j,i)[0]; // 'r'
		}
}

void Image::Set(int numRows, int numCols, int numChannels){
	this->numChannels = numChannels;
	this->numCols = numCols;
	this->numRows = numRows;
	this->imgData = new float [this->numChannels*this->numCols*this->numRows];
}

void Image::myNormalize(int foo){
	for(int i=0;i< this->numRows*this->numCols*this->numChannels;i++)
		this->imgData[i] = this->imgData[i]/255;
}

void Image::imPadWrapper(Image & img_warpped, int  pt, int pb, int pl, int pr, char*type){
	int flag=0;
	float val=0;
	if(!strcmp(type,"replicate")) flag=1;
	else if(!strcmp(type,"symmetric")) flag=2;
	else if(!strcmp(type,"circular")) flag=3;
	img_warpped.numChannels = this->numChannels;
	img_warpped.numRows = this->numRows+pt+pb;
	img_warpped.numCols = this->numCols+pl+pr;
	img_warpped.imgData = new float[img_warpped.numChannels*img_warpped.numRows*img_warpped.numCols];
	imPad( this->imgData, img_warpped.imgData,this->numRows, this->numCols,  this->numChannels,  
		pt, pb, pl,  pr, flag, val );
	img_warpped.numCols = this->numCols + pl + pr;
	img_warpped.numRows = this->numRows + pt + pb;
	return;
}

void Image::rgbConvertWrapper(Image &J, string colorspace="luv", bool useSingle=true){
	int flag = 0;
	string outClass;
	if(colorspace == "rgb"){
		flag = 1;
	}else if(colorspace == "luv"){
		flag = 2;
	}else{
		std::cout << "wrong colour space"<<std::endl;
	}
	if(useSingle==false){
		outClass = "double";
	}else{
		outClass = "single";
	}



	if(useSingle==true )
		J.imgData = rgbConvert(this->imgData, this->numCols *this->numRows, this->numChannels, flag, 1.0 );
	J.numCols = this->numCols;
	J.numRows = this->numRows;
	J.numChannels = this->numChannels;
}


void Image::imChnsWrapper(Image &img,  Option &opts){
}

void Image::imgResampleWrapper(Image &img_shrink, float alpha, std::string method="", float norm=1){
	bool bilinear = true;
	if(method == ""){
		bilinear = true;
	}
	if(bilinear){
		if(int(this->numRows*alpha)!=this->numRows*alpha && int(this->numCols*alpha)!=this->numCols*alpha)
			std::cerr<<"imgResample Error"<<std::endl;
		img_shrink.Set(this->numRows*alpha, this->numCols*alpha, this->numChannels);
		imResample( this->imgData, img_shrink.imgData, this->numRows, int(this->numRows*alpha), 
			 this->numCols, int(this->numCols*alpha),
			this->numChannels, norm );
	}
}

void Image::gradWrapper(Image &H, Image &M, bool full, int normRad, float normConst =0.005, int shrink=1, int nOrients=9, int softBin=1){
	int binSize = max(shrink, 1);
	H.Set(this->numRows/binSize, this->numCols/binSize, nOrients);
	H.imgData = new float[this->numCols/binSize*this->numRows/binSize*nOrients];
	memset(H.imgData,0,this->numCols/binSize*this->numRows/binSize*nOrients*sizeof(float));   
	M.Set(this->numRows, this->numCols,1);
	M.imgData = new float[this->numCols*this->numRows];
	float *O_imgData = new float[this->numCols*this->numRows];
	float* M_convTri = new float[this->numCols*this->numRows];

	gradMag( this->imgData, M.imgData, O_imgData,  this->numRows, this->numCols, this->numChannels, full );
	
/*	std::cout<<"M:"<<std::endl;
	for(int i=0;i<10;i++){
		for (int j=0;j<10;j++)
			std::cout<<M.imgData[this->numCols*this->numRows*0+j*this->numRows+i]<<' ';
		std::cout<<std::endl;
	}
*/
	//std::cout<<"O:"<<std::endl;	
	//for(int i=0;i<10;i++){
	//	for (int j=0;j<10;j++)
	//		std::cout<<O_imgData[this->numCols*this->numRows*0+j*this->numRows+i]<<' ';
	//	std::cout<<std::endl;
	//}
	//
	convTri( M.imgData, M_convTri, this->numRows, this->numCols,  1, normRad, 1 );
	//std::cout<<"M_convTri:"<<std::endl;
	//for(int i=0;i<10;i++){
	//	for (int j=0;j<10;j++)
	//		std::cout<<M_convTri[this->numCols*this->numRows*0+j*this->numRows+i]<<' ';
	//	std::cout<<std::endl;
	//}
	gradMagNorm( M.imgData, M_convTri, this->numRows, this->numCols, normConst );
	//	std::cout<<"M:"<<std::endl;
	//for(int i=0;i<10;i++){
	//	for (int j=0;j<10;j++)
	//		std::cout<<M.imgData[this->numCols*this->numRows*0+j*this->numRows+i]<<' ';
	//	std::cout<<std::endl;
	//}
	gradHist( M.imgData, O_imgData, H.imgData, this->numRows, this->numCols, binSize, nOrients, softBin, full );
	/*
	std::cout<<"H in the function:"<<std::endl;
	for(int i=0;i<10;i++){
		for (int j=0;j<10;j++)
			std::cout<<H.imgData[H.numCols*H.numRows*0+j*H.numRows+i]<<' ';
		std::cout<<std::endl;
	}*/
	delete []M_convTri;
	delete []O_imgData;
	



}






void Image::log(int channel=0){
	for(int i=0;i<10;i++){
		for (int j=0;j<10;j++)
			std::cout<<this->imgData[this->numCols*this->numRows*channel+j*this->numRows+i]<<' ';
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void Image::edgesNmsWrapper(Image &E0, Image &E, int r, int s, int m, int nThreads){
	int h=E0.numRows;
	int w=E0.numCols;
	E.Set(h,w,E0.numChannels);
	E.imgData = new float[h*w*E0.numChannels];

		// suppress edges where edge is stronger in orthogonal direction
	#ifdef USEOMP
	nThreads = nThreads<omp_get_max_threads() ? nThreads : omp_get_max_threads();
	#pragma omp parallel for num_threads(nThreads)
	#endif

	for( int x=0; x<w; x++ ) for( int y=0; y<h; y++ ) {
		float e=E.imgData[x*h+y]=E0.imgData[x*h+y]; if(!e) continue; e*=m;
		float coso=cos(this->imgData[x*h+y]), sino=sin(this->imgData[x*h+y]);
		for( int d=-r; d<=r; d++ ) if( d ) {
			float e0 = interp(E0.imgData,h,w,x+d*coso,y+d*sino);
			if(e < e0) { E.imgData[x*h+y]=0; break; }
		}
	}
	// suppress noisy edge estimates near boundaries
	s=s>w/2?w/2:s; s=s>h/2? h/2:s;
	for( int x=0; x<s; x++ ) for( int y=0; y<h; y++ ) {
		E.imgData[x*h+y]*=x/float(s); E.imgData[(w-1-x)*h+y]*=x/float(s); }
	for( int x=0; x<w; x++ ) for( int y=0; y<s; y++ ) {
		E.imgData[x*h+y]*=y/float(s); E.imgData[x*h+(h-1-y)]*=y/float(s); }



}
