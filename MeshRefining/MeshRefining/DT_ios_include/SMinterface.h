
/*
	@����smooth�����������ڴ棬��smooth����ǰ����
*/
void SetSmoothParameter(void **pvSmoothData, int itech='f', int igoal=5, double dThresh = 30.0);

/*
	@��ResetSmooth�ɶԵ��ã�ÿִ����һ��smooth����һ��
*/
void InitSmooth(void *pvSmoothData);

/*
	@��InitSmooth�ɶԵ��ã�ÿִ����һ��smooth����һ��
*/
void ResetSmooth(void *pvSmoothData);

void OpenMessageFile(char *filename);


/*
	@(in)num_incident_vtx: number of neighboring vertex �ڽӵ���� 
	@(in)num_tet��number of neighboring tetrahedron �ڽ������嵥Ԫ����
	@(in)free_vtx��coordinate of current vertex ��ǰ������������
	@(in)vtx_list��index list of neighboring vertex �ڽӶ�������
	@(in)vtx_connectivity��list of connectivity �ڽӹ�ϵ
	@(out)qual: minimal quality of the local mesh �ֲ��������С����
	@(out)cellqual: the minimal quality of every neighboring cells ÿ���ڽӵ�Ԫ��������Сֵ
	@(out)new_vtx: new coordinates of the vertex ������λ������

	tips; the output values are valid only when the function returns successfully(return 1).
	return 0: failed, 1: succeed 
*/
int VertexSmooth(int num_incident_vtx, int num_tet, double *free_vtx, 
	double **vtx_list, int **vtx_connectivity, double *qual, double *cellqual, double *new_vtx, void *pvSmoothData);