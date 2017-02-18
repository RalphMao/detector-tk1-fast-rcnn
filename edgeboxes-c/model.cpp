# include "model.h"
// #include <fcntl.h>
// # include <unistd.h>
//# include <iostream>
#include "ttq.h"
#include <cstdio>
#include <omp.h>

#define _TTQ_READ(_dst, _size, _numb, _fd)         \
		(read((_fd), (_dst), (_numb)*(_size)))


void Model::load(char* filename){
	FILE*  fid = fopen(filename, "rb");

	
	this->thresh_dim= new int(0);
	
	fread((void*)this->thresh_dim, sizeof(int), 1, fid);	
	this->thresh_size = new int [*this->thresh_dim];
	
	fread((void*)this->thresh_size, sizeof(int), *this->thresh_dim, fid);
	//std::cout << "dim=" << *this->thresh_dim <<", size="<< this->thresh_size[0] <<"," <<this->thresh_size[1]<<std::endl;
	this->thresh_ele_num = 1;

	for(int i=0; i<*this->thresh_dim; i++){
		this->thresh_ele_num *= this->thresh_size[i];
		// ttqDebug("size %d: %d, %d\n", i, this->thresh_ele_num, this->thresh_size[i]);	
	}
	
	
	this->thresh = new float[this->thresh_ele_num];
	//std::cout<<this->thresh_ele_num << std::endl;
	
	fread((void*)this->thresh, sizeof(float), this->thresh_ele_num, fid);
	//std::cout << n <<std::endl;
	/*for(int i=this->thresh_ele_num-8;i<this->thresh_ele_num;i++)
		std::cout << (float)this->thresh[i]<<", ";
	std::cout << std::endl;*/
	
	
	//ttqDebug("%s\n", "Here we go");	

	this->fids_dim= new int(0);
	fread((void*)this->fids_dim, sizeof(int),1, fid);	
	this->fids_size = new int [*this->fids_dim];
	fread((void*)this->fids_size, sizeof(int),*this->fids_dim, fid);
	//std::cout << "dim=" << *this->fids_dim <<", size="<< this->fids_size[0] <<"," <<this->fids_size[1]<<std::endl;
	this->fids_ele_num = 1;
	for(int i=0; i<*this->fids_dim; i++)
		this->fids_ele_num = this->fids_ele_num*this->fids_size[i];
	this->fids = new uint32 [this->fids_ele_num];
	fread((void*)this->fids, sizeof( uint32),this->fids_ele_num, fid);
	//for(int i=this->fids_ele_num-8;i<this->fids_ele_num;i++){
	//	std::cout << (int)this->fids[i]<<", ";

	//}
	//std::cout << std::endl;

		this->child_dim= new int(0);
	fread((void*)this->child_dim, sizeof(int),1, fid);	
	this->child_size = new int [*this->child_dim];
	fread((void*)this->child_size, sizeof(int),*this->child_dim, fid);
	//std::cout << "dim=" << *this->child_dim <<", size="<< this->child_size[0] <<"," <<this->child_size[1]<<std::endl;
	this->child_ele_num = 1;
	for(int i=0; i<*this->child_dim; i++)
		this->child_ele_num = this->child_ele_num*this->child_size[i];
	this->child = new uint32 [this->child_ele_num];
	fread((void*)this->child, sizeof( uint32),this->child_ele_num, fid);
	//for(int i=this->child_ele_num-8;i<this->child_ele_num;i++)
	//	std::cout << (int)this->child[i]<<", ";
	//std::cout << std::endl;

	this->count_dim= new int(0);
	fread((void*)this->count_dim, sizeof(int),1, fid);	
	this->count_size = new int [*this->count_dim];
	fread((void*)this->count_size, sizeof(int),*this->count_dim, fid);
	//std::cout << "dim=" << *this->count_dim <<", size="<< this->count_size[0] <<"," <<this->count_size[1]<<std::endl;
	this->count_ele_num = 1;
	for(int i=0; i<*this->count_dim; i++)
		this->count_ele_num = this->count_ele_num*this->count_size[i];
	this->count = new uint32 [this->count_ele_num];
	fread((void*)this->count, sizeof( uint32),this->count_ele_num, fid);
	//for(int i=0;i<8;i++)
	//	std::cout << this->count[i]<<", ";
	//std::cout << std::endl;

		this->depth_dim= new int(0);
	fread((void*)this->depth_dim, sizeof(int),1, fid);	
	this->depth_size = new int [*this->depth_dim];
	fread((void*)this->depth_size, sizeof(int),*this->depth_dim, fid);
	//std::cout << "dim=" << *this->depth_dim <<", size="<< this->depth_size[0] <<"," <<this->depth_size[1]<<std::endl;
	this->depth_ele_num = 1;
	for(int i=0; i<*this->depth_dim; i++)
		this->depth_ele_num = this->depth_ele_num*this->depth_size[i];
	this->depth = new uint32 [this->depth_ele_num];
	fread((void*)this->depth, sizeof( uint32),this->depth_ele_num, fid);
	//for(int i=0;i<8;i++)
	//	std::cout << this->depth[i]<<", ";
	//std::cout << std::endl;

	this->segs_dim= new int(0);
	fread((void*)this->segs_dim, sizeof(int),1, fid);	
	this->segs_size = new int [*this->segs_dim];
	fread((void*)this->segs_size, sizeof(int),*this->segs_dim, fid);
	//std::cout << "dim=" << *this->segs_dim <<", size="<< this->segs_size[0] <<"," <<this->segs_size[1]<<std::endl;
	this->segs_ele_num = 1;
	for(int i=0; i<*this->segs_dim; i++)
		this->segs_ele_num = this->segs_ele_num*this->segs_size[i];
	this->segs = new uint8 [this->segs_ele_num];
	fread((void*)this->segs, sizeof( uint8),this->segs_ele_num, fid);
	//for(int i=0;i<8;i++)
	//	std::cout << this->segs[i]<<", ";
	//std::cout << std::endl;

			this->nSegs_dim= new int(0);
	fread((void*)this->nSegs_dim, sizeof(int),1, fid);	
	this->nSegs_size = new int [*this->nSegs_dim];
	fread((void*)this->nSegs_size, sizeof(int),*this->nSegs_dim, fid);
	//std::cout << "dim=" << *this->nSegs_dim <<", size="<< this->nSegs_size[0] <<"," <<this->nSegs_size[1]<<std::endl;
	this->nSegs_ele_num = 1;
	for(int i=0; i<*this->nSegs_dim; i++)
		this->nSegs_ele_num = this->nSegs_ele_num*this->nSegs_size[i];
	this->nSegs = new uint8 [this->nSegs_ele_num];
	fread((void*)this->nSegs, sizeof( uint8),this->nSegs_ele_num, fid);
	//for(int i=0;i<8;i++)
	//	std::cout << this->nSegs[i]<<", ";
	//std::cout << std::endl;

			this->eBins_dim= new int(0);
	fread((void*)this->eBins_dim, sizeof(int),1, fid);	
	this->eBins_size = new int [*this->eBins_dim];
	fread((void*)this->eBins_size, sizeof(int),*this->eBins_dim, fid);
	//std::cout << "dim=" << *this->eBins_dim <<", size="<< this->eBins_size[0] <<"," <<this->eBins_size[1]<<std::endl;
	this->eBins_ele_num = 1;
	for(int i=0; i<*this->eBins_dim; i++)
		this->eBins_ele_num = this->eBins_ele_num*this->eBins_size[i];
	this->eBins = new uint16 [this->eBins_ele_num];
	fread((void*)this->eBins, sizeof( uint16),this->eBins_ele_num, fid);
	//for(int i=0;i<8;i++)
	//	std::cout << this->eBins[i]<<", ";
	//std::cout << std::endl;

		this->eBnds_dim= new int(0);
	fread((void*)this->eBnds_dim, sizeof(int),1, fid);	
	this->eBnds_size = new int [*this->eBnds_dim];
	fread((void*)this->eBnds_size, sizeof(int),*this->eBnds_dim, fid);
	//std::cout << "dim=" << *this->eBnds_dim <<", size="<< this->eBnds_size[0] <<"," <<this->eBnds_size[1]<<std::endl;
	this->eBnds_ele_num = 1;
	for(int i=0; i<*this->eBnds_dim; i++)
		this->eBnds_ele_num = this->eBnds_ele_num*this->eBnds_size[i];
	this->eBnds = new uint32 [this->eBnds_ele_num];
	fread((void*)this->eBnds, sizeof( uint32),this->eBnds_ele_num, fid);
	//for(int i=0;i<8;i++)
	//	std::cout << this->eBnds[i]<<", ";
	//std::cout << std::endl;

	fclose(fid);
}

Model:: ~Model(){
	delete []this->thresh; delete [] this->thresh_size; delete this->thresh_dim; 
	delete []this->fids; delete []this->fids_size; delete this->fids_dim; 
	delete []this->child; delete []this->child_size; delete this-> child_dim; 
	delete []this->count; delete []this->count_size; delete this->count_dim; 
	delete []this->depth; delete []this->depth_size; delete this->depth_dim; 
	delete []this->segs; delete []this->segs_size; delete this->segs_dim; 
	delete []this->nSegs; delete []this->nSegs_size; delete this-> nSegs_dim; 
	delete []this->eBins; delete []this->eBins_size; delete this->eBins_dim; 
	delete []this->eBnds; delete []this-> eBnds_size; delete this->eBnds_dim; 

}
