# ifndef __MODEL__
typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;
class Model{
public:
	Model() {}
	~Model();
	void load(char* filename);

	float* thresh; int* thresh_size; int* thresh_dim; int thresh_ele_num;
	uint32* fids; int*fids_size; int* fids_dim; int fids_ele_num;
	uint32* child; int * child_size; int* child_dim; int child_ele_num;
	uint32* count; int * count_size; int* count_dim; int count_ele_num;
	uint32* depth; int * depth_size; int* depth_dim; int depth_ele_num;
	uint8* segs; int * segs_size; int* segs_dim; int segs_ele_num;
	uint8* nSegs; int * nSegs_size; int* nSegs_dim; int nSegs_ele_num;
	uint16* eBins; int *eBins_size; int* eBins_dim; int eBins_ele_num;
	uint32* eBnds; int *eBnds_size; int* eBnds_dim; int eBnds_ele_num;



};

# define __MODEL__
# endif