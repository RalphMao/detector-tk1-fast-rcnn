# include "Option.h"
#include <cstdio>
void Option::edgeDection(Image &img, Model & model,Image &E_out, Image & G){
	nTreesEval = (nTreesEval<nTrees)?nTreesEval:nTrees;
	stride= (stride>shrink)?stride:shrink;
	if(multiscale){
	}
	else{
		int pt=imWidth/2+(4-(img.numRows+imWidth)%4)%4;
		int pb=imWidth/2+(4-(img.numRows+imWidth)%4)%4;
		int pl=imWidth/2+(4-(img.numCols+imWidth)%4)%4;
		int pr=imWidth/2+(4-(img.numCols+imWidth)%4)%4;
		char* type = "symmetric";
		
		Image img_pad;
		img.imPadWrapper(img_pad,pt, pb, pl, pr, type);
		img_pad.myNormalize(1);
		//std::cout<<"imPad:"<<std::endl;
		//img_pad.log(0);
		//img_pad.log(1);
		//img_pad.log(2);

		Image img_cs;
		img_pad.rgbConvertWrapper(img_cs, "luv", 1);
		//std::cout<<"rgbConvert"<<std::endl;
		//img_cs.log(0);
		//img_cs.log(1);
		//img_cs.log(2);
		
		Image img_shrink;
		img_cs.imgResampleWrapper(img_shrink, float(1.0)/float(this->shrink), "", 1);
		//std::cout << "opts.shrink=" << this->shrink <<std::endl;		// the defition of shrink is int
		//img_shrink.log(0);
		//img_shrink.log(1);
		//img_shrink.log(2);
	//	for (int s=1;s<3;s++){
		Image* temp_img;
		Image *chns = new Image[4];
//		std::cout << "nChns=" <<this->nChns <<std::endl;
		
		int s=1;
		if(s==this->shrink)
			temp_img = &img_shrink;
		else
			temp_img = &img_cs;

		Image H1, M1;
		(*temp_img).gradWrapper(H1, M1, 0 , this->normRad, float(0.01), float(this->shrink)/float(s), this->nOrients, 0);
	//	std::cout << "H="  <<std::endl;		
	//	H1.log(0);


		M1.imgResampleWrapper(chns[0], float(s)/float(this->shrink), "", 1);
		//std::cout << "chns[1]="  <<std::endl;		
		//chns[0].log(0);

		H1.imgResampleWrapper(chns[1], max(float(1.0), float(s)/float(this->shrink)), "", 1);
		//std::cout << "chns[2]="  <<std::endl;		
		//chns[1].log(3);
		
		s=2;
		if(s==this->shrink)
			temp_img = &img_shrink;
		else
			temp_img = &img_cs;

		Image H2, M2;
		(*temp_img).gradWrapper(H2, M2, 0 , this->normRad, float(0.01), float(this->shrink)/float(s), this->nOrients, 0);
		//std::cout << "H2="  <<std::endl;		
		//H2.log(0);
		//std::cout << "M2="  <<std::endl;		
		//M2.log(0);

		M2.imgResampleWrapper(chns[2], float(s)/float(this->shrink), "", 1);
		//std::cout << "chns[3]="  <<std::endl;		
		//chns[2].log(0);

		H2.imgResampleWrapper(chns[3], max(float(1.0), float(s)/float(this->shrink)), "", 1);
		//std::cout << "chns[4]="  <<std::endl;		
		//chns[3].log(3);

		Image chns_temp;
		int stride = 0, step = 0;
		chns_temp.Set(img_shrink.numRows, img_shrink.numCols, this->nChns);
		chns_temp.imgData = new float [img_shrink.numRows*img_shrink.numCols*this->nChns];
		step = img_shrink.numRows*img_shrink.numCols*img_shrink.numChannels;
		memcpy(chns_temp.imgData, img_shrink.imgData, step*sizeof(float));
		
		for(int i=0;i<4;i++){
			stride = stride + step;
			step = chns[i].numCols*chns[i].numRows*chns[i].numChannels;
			memcpy(chns_temp.imgData+stride, chns[i].imgData, step*sizeof(float));
		}
		Image chnsReg;
		Image chnsSim;

		int chnSm=this->chnSmooth/this->shrink; 
		int simSm=this->simSmooth/this->shrink;
		chnsReg.Set(chns_temp.numRows, chns_temp.numCols, chns_temp.numChannels);
		chnsSim.Set(chns_temp.numRows, chns_temp.numCols, chns_temp.numChannels);
		chnsReg.imgData = new float [chns_temp.numRows* chns_temp.numCols* chns_temp.numChannels];
		chnsSim.imgData = new float [chns_temp.numRows* chns_temp.numCols* chns_temp.numChannels];
		convTri(chns_temp.imgData,chnsReg.imgData, chns_temp.numRows, chns_temp.numCols, chns_temp.numChannels, chnSm, 1); 
		convTri(chns_temp.imgData,chnsSim.imgData, chns_temp.numRows, chns_temp.numCols, chns_temp.numChannels, simSm, 1);

		s=this->sharpen;
		Image img_rgb, img_conv;
		
		if(s){
			img_pad.rgbConvertWrapper(img_rgb, "rgb", 1);
			img_conv.Set(img_rgb.numRows, img_rgb.numCols, img_rgb.numChannels);
			img_conv.imgData = new float[img_rgb.numRows* img_rgb.numCols* img_rgb.numChannels];
			convTri(img_rgb.imgData, img_conv.imgData, img_rgb.numRows, img_rgb.numCols, img_rgb.numChannels, 1, 1);
		}
		//img_rgb.log(0);
		//img_conv.log(0);
		// create outputs
		const int h1 = (int) ceil(double(img_rgb.numRows-this->imWidth)/this->stride);
		const int w1 = (int) ceil(double(img_rgb.numCols-this->imWidth)/this->stride);
		const int h2 = h1*this->stride+this->gtWidth;
		const int w2 = w1*this->stride+this->gtWidth;
		float *E = new float [h2*w2*1];	
		memset(E,0,h2*w2*1*sizeof(float));
		uint32 *ind = new uint32 [img_rgb.numRows*img_rgb.numCols*this->nTreesEval];
		memset(ind,0,img_rgb.numRows*img_rgb.numCols*this->nTreesEval*sizeof(uint32));
		//int *segsOut = new int [this->gtWidth*this->gtWidth*h1*w1*this->nTreesEval];
		//memset(segsOut,0,this->gtWidth*this->gtWidth*h1*w1*this->nTreesEval*sizeof(int));

		edgeDectionExternal(model, img_conv.imgData, chnsReg.imgData, chnsSim.imgData, 
			img_rgb.numRows, img_rgb.numCols, img_rgb.numChannels, 
			E, ind);

		float t= float(this->stride*this->stride)/float(this->gtWidth*this->gtWidth)/float(this->nTreesEval);
		int r= this->gtWidth/2;
		if(this->sharpen == 0)
			t = t*2;
		else if(this->sharpen == 1)
			t = t*1.8f;
		else
			t = t*1.66f;

		float *E_resize = new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		float *E_convconv = new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		E_out.imgData =  new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		E_out.numRows = h2-this->gtWidth;
		E_out.numCols = w2-this->gtWidth;
		E_out.numChannels = 1;
		float *G_x =  new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		float *G_y =  new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		float *G_xx =  new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		float *G_xy =  new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		float *G_yy =  new float [(h2-this->gtWidth)*(w2-this->gtWidth)];
		G.numRows = h2-this->gtWidth;
		G.numCols = w2-this->gtWidth;
		G.numChannels = 1;
		G.imgData = new float [(h2-this->gtWidth)*(w2-this->gtWidth)];

		for(int i=0;i<w2-this->gtWidth;i++){
			memcpy(E_resize+i*(h2-this->gtWidth), E+r+(i+r)*h2, (h2-this->gtWidth)*sizeof(float));
			for(int j=0;j<h2-this->gtWidth;j++)
				E_resize[i*(h2-this->gtWidth)+j] = E_resize[i*(h2-this->gtWidth)+j]*t;
		}
		convTri(E_resize, E_out.imgData, (h2-this->gtWidth),(w2-this->gtWidth), 1, 1, 1);
		convTri(E_out.imgData, E_convconv, (h2-this->gtWidth),(w2-this->gtWidth), 1, 4, 1);

		
		grad2(E_convconv, G_x, G_y, (h2-this->gtWidth),(w2-this->gtWidth), 1);
		grad2(G_x, G_xx, G_xy, (h2-this->gtWidth),(w2-this->gtWidth), 1);
		grad2(G_y, G_xy, G_yy, (h2-this->gtWidth),(w2-this->gtWidth), 1);
		
		for(int i=0;i<G.numRows;i++)
			for(int j=0;j<G.numCols;j++){
				float temp = atan(G_yy[i*G.numCols+j]*
				((G_xy[i*G.numCols+j]>0)?-1:(G_xy[i*G.numCols+j]==0)?0:1)/(G_xx[i*G.numCols+j]+1/1e5f));
		//		if(i*G.numCols+j == 182271-1)
		//			printf("temp = %f\n", temp);
				G.imgData[i*G.numCols+j] = temp - int(temp/PI)*PI;
				if(G.imgData[i*G.numCols+j]<0)
					G.imgData[i*G.numCols+j] = G.imgData[i*G.numCols+j] + PI;
			}

		
		delete []G_yy;
		delete []G_xy;
		delete []G_xx;
		delete []G_y;
		delete []G_x;
		delete []E_convconv;
		delete []E_resize;
		delete []ind;
		delete []E;
		delete []chns;
		//O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);


//		printf("%d, %d\n", (h2-this->gtWidth),(w2-this->gtWidth));
//		printf("%f, ",E_resize[182271-1]);
//		printf("%f, %f, %f, %f, %f, %f, %f\n",E_out.imgData[182271-1], G_x[182271-1], G_y[182271-1], G_xx[182271-1], G_xy[182271-1], G_yy[182271-1], G.imgData[182271-1]);
//		printf("%f, ",G.imgData[182271-2]);

	}
}



void Option::edgeChns(Image &img_before, Image &img_after){
	img_before.rgbConvertWrapper(img_after,"luv", 1);
	
}

void Option::edgeDectionExternal(Model &model, float* img, float* chns, float*chnsSs, int h, int w, int Z, 
								 float *E, uint32*ind){
	// get model
	float* thrs = model.thresh;
	uint32 *fids = model.fids;
	uint32 *child = model.child;
	uint8 *segs = model.segs;
	uint8 *nSegs = model.nSegs;
	uint16 *eBins = model.eBins;
	uint32 *eBnds = model.eBnds;

	const int nBnds = (model.eBnds_ele_num-1)/model.thresh_ele_num;
	
	// get opts
	const int shrink = this->shrink;
	const int imWidth = this->imWidth; 
	const int gtWidth = this->gtWidth;
	const int nChns = this->nChns; 
	const int nCells = this->nCells; 
	const uint32 nChnFtrs = this->nChnFtrs;	//std::cout <<"nChnFtrs="<<nChnFtrs<<std::endl;
	const int stride = this->stride;
	const int nTreesEval = this->nTreesEval;
	int sharpen = this->sharpen; 
	int nThreads = this->nThreads;
	const char *msgSharpen="Model supports sharpening of at most %i pixels!\n";
	if( sharpen>nBnds-1 ) {
		sharpen=nBnds-1; 
		printf(msgSharpen, sharpen);
	}

	// get dimensions and constants
	const int *fidsSize = model.fids_size;
	const int nTreeNodes = (int) fidsSize[0];
	const int nTrees = (int) fidsSize[1];
	const int h1 = (int) ceil(double(h-imWidth)/stride);
	const int w1 = (int) ceil(double(w-imWidth)/stride);
	const int h2 = h1*stride+gtWidth;
	const int w2 = w1*stride+gtWidth;
	const int imgDims[3] = {h,w,Z};
	const int chnDims[3] = {h/shrink,w/shrink,nChns};
	const int indDims[3] = {h1,w1,nTreesEval};
	const int outDims[3] = {h2,w2,1};
	const int segDims[5] = {gtWidth,gtWidth,h1,w1,nTreesEval};
			    //printf("size=%d,%d,%d,%d,%d,%d,%d,%d,%d\n",h,w,Z,h1,w1,h2,w2,stride,nTreesEval);

	// construct lookup tables
	uint32 *iids, *eids, *cids, *cids1, *cids2;
	iids = buildLookup( (int*)imgDims, gtWidth );
	eids = buildLookup( (int*)outDims, gtWidth );
	cids = buildLookup( (int*)chnDims, imWidth/shrink );
	buildLookupSs( cids1, cids2, (int*)chnDims, imWidth/shrink, nCells );



	// apply forest to all patches and store leaf inds
	#ifdef USEOMP
	nThreads = min(nThreads,omp_get_max_threads());
	#pragma omp parallel for num_threads(nThreads)
	#endif
	for( int c=0; c<w1; c++ ) for( int t=0; t<nTreesEval; t++ ) {
		for( int r0=0; r0<2; r0++ ) for( int r=r0; r<h1; r+=2 ) {
			int o = (r*stride/shrink) + (c*stride/shrink)*h/shrink;
			// select tree to evaluate
			int t1 = ((r+c)%2*nTreesEval+t)%nTrees; uint32 k = t1*nTreeNodes;
			while( child[k] ) {
			// compute feature (either channel or self-similarity feature)
			uint32 f = fids[k]; float ftr;
			if( f<nChnFtrs ) ftr = chns[cids[f]+o]; else
				ftr = chnsSs[cids1[f-nChnFtrs]+o]-chnsSs[cids2[f-nChnFtrs]+o];
			// compare ftr to threshold and move left or right accordingly
			if( ftr < thrs[k] ) k = child[k]-1; else k = child[k];
			k += t1*nTreeNodes;
			}
			// store leaf index and update edge maps
			ind[ r + c*h1 + t*h1*w1 ] = k;
		}
	}
	//  for(int i=0;i<10;i++)
	//printf("%d, ",ind[i]);
	//printf("\n");
	 //printf("stride=%d, imWidth=%d, gtWidth=%d, h=%d, w1=%d, h1=%d\n", stride, imWidth, gtWidth,h, w1, h1);


		// compute edge maps (avoiding collisions from parallel executions)
	if( !sharpen ) for( int c0=0; c0<gtWidth/stride; c0++ ) {
		#ifdef USEOMP
		#pragma omp parallel for num_threads(nThreads)
		#endif
		for( int c=c0; c<w1; c+=gtWidth/stride ) {
			for( int r=0; r<h1; r++ ) for( int t=0; t<nTreesEval; t++ ) {
				uint32 k = ind[ r + c*h1 + t*h1*w1 ];
				float *E1 = E + (r*stride) + (c*stride)*h2;
				int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+1]; if(b0==b1) continue;
				for( int b=b0; b<b1; b++ ) E1[eids[eBins[b]]]++;
			//	memcpy(segsOut+(r+c*h1+t*h1*w1)*gtWidth*gtWidth,
			//		segs+k*gtWidth*gtWidth,gtWidth*gtWidth*sizeof(int));
			}
		}
	}

		// computed sharpened edge maps, snapping to local color values
	if( sharpen ) {
		// compute neighbors array
		const int g=gtWidth; uint16 N[4096*4];
		for( int c=0; c<g; c++ ) for( int r=0; r<g; r++ ) {
			int i=c*g+r; uint16 *N1=N+i*4;
			N1[0] = c>0 ? i-g : i; N1[1] = c<g-1 ? i+g : i;
			N1[2] = r>0 ? i-1 : i; N1[3] = r<g-1 ? i+1 : i;
		}
		#ifdef USEOMP
		#pragma omp parallel for num_threads(nThreads)
		#endif
		int count=0;
		for( int c=0; c<w1; c++ ) for( int r=0; r<h1; r++ ) {
			for( int t=0; t<nTreesEval; t++ ) {
			// get current segment and copy into S
			
				uint32 k = ind[ r + c*h1 + t*h1*w1 ];
				int m = nSegs[k]; if( m==1 ) continue;
				//	if(count<10){
				////	printf("%d, %d\n", b0, b1);
				//	count = count+1;
				//}
				//	if(count == 5)
				//		printf("I am at count ==5\n");
				uint8 S0[4096], *S= S0;//*S= segsOut+(r+c*h1+t*h1*w1)*g*g;
				memcpy(S,segs+k*g*g, g*g*sizeof(uint8));
				/*if(count==5){
					printf("I am S\n");
					for(int i=0;i<16;i++){
						for(int j=0;j<16;j++)
							printf("%d, ", S[i*16+j]);
						printf("\n");
					}
				}
*/
				// compute color model for each segment using every other pixel
				int ci, ri, s, z; float ns[100], mus[1000];
				const float *I1 = img+(c*stride+(imWidth-g)/2)*h+r*stride+(imWidth-g)/2;
				for( s=0; s<m; s++ ) { ns[s]=0; for( z=0; z<Z; z++ ) mus[s*Z+z]=0; }
				for( ci=0; ci<g; ci+=2 ) for( ri=0; ri<g; ri+=2 ) {
					s = S[ci*g+ri]; ns[s]++;
					for( z=0; z<Z; z++ ) mus[s*Z+z]+=I1[z*h*w+ci*h+ri];
				}
				for(s=0; s<m; s++) for( z=0; z<Z; z++ ) mus[s*Z+z]/=ns[s];
				// update segment S according to local color values
				int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+sharpen];
				//if(count<10){
				//	printf("%d, %d\n", b0, b1);
				////	count = count+1;
				//}
				for( int b=b0; b<b1; b++ ) {
					float vs[10], d, e, eBest=1e10f; int i, sBest=-1, ss[4];
					for( i=0; i<4; i++ ) ss[i]=S[N[eBins[b]*4+i]];
					for( z=0; z<Z; z++ ) {
						vs[z]=I1[iids[eBins[b]]+z*h*w];
						/*if(b==b0 && count == 5)
							printf("vs[z]=%f, ",vs[z]);*/
					}
					//if(b==b0 && count == 5)
						//printf("\n");
					for( i=0; i<4; i++ ) {
						s=ss[i]; if(s==sBest) continue;
						e=0; for( z=0; z<Z; z++ ) { d=mus[s*Z+z]-vs[z]; e+=d*d; }
						if( e<eBest ) { eBest=e; sBest=s; }
					}
					S[eBins[b]]=sBest;
					/*if(count == 5){
						printf("b=%d, eBins[b]=%d, S[eBins[b]]=%d, iids[eBins[b]]=%d, sBest=%d, d=%f\n",
							b, eBins[b], S[eBins[b]], iids[eBins[b]], sBest, d);
					}*/
				}
		/*		if(count<10){
					for(int i=0;i<10;i++)
						printf("%d, ",(int)S0[i]);
					printf("I am here!\n");
				}*/
				// convert mask to edge maps (examining expanded set of pixels)
				float *E1 = E + c*stride*h2 + r*stride; b1=eBnds[k*nBnds+sharpen+1];
				for( int b=b0; b<b1; b++ ) {
					int i=eBins[b]; uint8 s=S[i]; uint16 *N1=N+i*4;
					if( s!=S[N1[0]] || s!=S[N1[1]] || s!=S[N1[2]] || s!=S[N1[3]] )
					E1[eids[i]]++;
				}
			}
		}
	}
	

	// free memory
	delete [] iids; delete [] eids;
	delete [] cids; delete [] cids1; delete [] cids2;
}
