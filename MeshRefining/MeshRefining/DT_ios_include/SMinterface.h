
/*
	@设置smooth参数，分配内存，在smooth操作前调用
*/
void SetSmoothParameter(void **pvSmoothData, int itech='f', int igoal=5, double dThresh = 30.0);

/*
	@与ResetSmooth成对调用，每执行完一遍smooth调用一次
*/
void InitSmooth(void *pvSmoothData);

/*
	@与InitSmooth成对调用，每执行完一遍smooth调用一次
*/
void ResetSmooth(void *pvSmoothData);

void OpenMessageFile(char *filename);


/*
	@(in)num_incident_vtx: number of neighboring vertex 邻接点个数 
	@(in)num_tet：number of neighboring tetrahedron 邻接四面体单元个数
	@(in)free_vtx：coordinate of current vertex 当前操作顶点坐标
	@(in)vtx_list：index list of neighboring vertex 邻接顶点坐标
	@(in)vtx_connectivity：list of connectivity 邻接关系
	@(out)qual: minimal quality of the local mesh 局部网格的最小质量
	@(out)cellqual: the minimal quality of every neighboring cells 每个邻接单元的质量最小值
	@(out)new_vtx: new coordinates of the vertex 顶点新位置坐标

	tips; the output values are valid only when the function returns successfully(return 1).
	return 0: failed, 1: succeed 
*/
int VertexSmooth(int num_incident_vtx, int num_tet, double *free_vtx, 
	double **vtx_list, int **vtx_connectivity, double *qual, double *cellqual, double *new_vtx, void *pvSmoothData);