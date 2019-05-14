/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#ifndef __iso3d_polytri_h__
#define __iso3d_polytri_h__

typedef int BOOL;
#define TRUE 1
#define FALSE 0

typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;
struct tVertexStructure
{
	double x, y;
	int vnum;
	tVertex prev, next;
	BOOL ear;
	int out_tri;
	int prt;
};

typedef struct Tri2D
{
	int forms[3];
	int neigs[3];
} Tri2D;

BOOL Left(double x1, double y1,
		  double x2, double y2,
		  double x3, double y3);

BOOL LeftOn(double x1, double y1,
		    double x2, double y2,
		    double x3, double y3);

BOOL Collinear(double x1, double y1,
		       double x2, double y2,
		       double x3, double y3);


BOOL InCone(tVertex a, tVertex b);

BOOL Xor(BOOL x, BOOL y);


BOOL IntersectProp(double x1, double y1, double x2, double y2,
				   double x3, double y3, double x4, double y4);

BOOL Between(double x1, double y1, 
			 double x2, double y2, 
			 double x3, double y3);

BOOL Intersect(double x1, double y1, double x2, double y2,
			   double x3, double y3, double x4, double y4);

BOOL Diagonalie(tVertex a, tVertex b);

BOOL Diagonal(tVertex a, tVertex b);

void EarInit(void);

int EarRemoval(void);

void Unitize();

void Init(double x[], double y[], int nn, int nf);

void Swap();

int SwapBase(int it1, int it2);

void calcTriParams(
				double x1, double y1, 
				double x2, double y2, 
				double x3, double y3,
				double *doubArea, 
				double *cenx, double *ceny,
				double *sqRad
				);

void Output(int prt[], int forms[], int neigs[], int nf);

void FreeMemory();

#  ifdef __cplusplus
extern "C" {
#  endif /* __cplusplus */

#ifdef __PolyTri_Test__

int readFr2(const char *fname, double x[], double y[], int *nn);
int writeFr2(const char *fname, double x[], double y[], int nn);
int writePL2(const char *fname, double x[], double y[], int prt[], int nn,
			 int forms[], int nf);

#endif /* __PolyTri_Test__ */

int triangulate(double x[], double y[], int nn,
	int prt[], int forms[], int neigs[], int nf);

#  ifdef __cplusplus
}
#  endif /* __cplusplus */
#endif /* __iso3d_polytri_h__ */