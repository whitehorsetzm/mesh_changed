/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2008年06月24日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, June. 24, 2008
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953166
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：ferguson_surface.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对曲线的通用接口
 *     Define a set of interface for Ferguson surfaces
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/


#ifndef __eemas_ferguson_surface_h__
#define __eemas_ferguson_surface_h__

#include "vector.h"
#include "surface.h"
#include "GeomT.h"


namespace EEMAS {

	typedef double Mat1by1d[1];
	typedef double Mat1by2d[2];
	typedef double Mat2by1d[2];
	typedef double Mat2by2d[4];

class FergusonSurface: public Surface
{
public:
  int flag;
	/* --------------------------------------------------------------------
	 * 构造/析构函数
	 * -------------------------------------------------------------------*/
	FergusonSurface();
	~FergusonSurface();

	/* --------------------------------------------------------------------
	 * 释放空间
	 * -------------------------------------------------------------------*/
	int destroy();

	/* --------------------------------------------------------------------
	 * 初始化一条Furgeson曲面
	 * 输入参数
	 *     data_points, nu, nv			型值点数组
	 * 描述:
	 *     支持自由端点条件
	 * ------------------------------------------------------------------*/
	int initialize(
	double data_points[], int nu, int nv		  /* 型值点数组 */);

	/* --------------------------------------------------------------------
	 * 初始化一条Furgeson曲面
	 * 输入参数
	 *     data_points, nu, nv			型值点数组
	 * 描述:
	 *     支持自由端点条件
	 * ------------------------------------------------------------------*/
	int initialize(Vector data_points[], int nu, int nv		  /* 型值点数组 */);
     void get_box(double x[2],double y[2],double z[2]);
	/* --------------------------------------------------------------------
	 * 通过点坐标获取相应的参数值
	 * 输入参数：
	 *   coord		点的坐标
	 *   u, v		参数坐标(返回值)
	 *   tol		点离曲面距离的最大允许值，小于此值，认为点在曲面上
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	virtual int coord_to_param(Vector coord,
				double *u, double *v, 
				double tol = 1.0e-4);

	/* --------------------------------------------------------------------
	 * 通过参数值获取相应点坐标的
	 * 输入参数：
	 *   u, v		参数坐标
	 *   coord		点的坐标(返回值)
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	int	param_to_coord(double u, double v, Vector *coord);

	int param_to_plane_coord(double u, double v, double *x, double *y);

	/* --------------------------------------------------------------------
	 * 获取某点在曲面上的投影点
	 * 输入参数：
	 *   coord		点的坐标
	 *   u, v		最近点(垂足)的参数坐标(返回值)
	 *   tol		点离曲面距离的最大允许值，小于此值，认为点在曲面上
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	virtual int project(Vector coord, double *u, double *v, double tol = 1.0e-4);
	
	virtual int projectSTNB(Vector coord, double *u, double *v);

	virtual int projectSTNB(int iN, Vector *pCoord, double *&pOutUV);
	
	/* --------------------------------------------------------------------
	 * 获取一阶偏导数
	 * 输入参数：
	 *   u, v		参数坐标
	 *   d1u, d1v	u向偏导/v向偏导
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	int	param_to_uv_d1(double u, double v, Vector *d1u, Vector *d1v);

	/* --------------------------------------------------------------------
	 * 获取二阶偏导数
	 * 输入参数：
	 *   u, v		参数坐标
	 *   d2uv		二阶偏导
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * ------------------------------------------------------------------*/
	int	param_to_uv_d2(double u, double v, Vector *d2uv);

	/* --------------------------------------------------------------------
	 * 获取二阶偏导数
	 * 输入参数：
	 *   u, v		参数坐标
	 *   d2u2
	 *   d2v2
	 *   d2uv		二阶偏导
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * ------------------------------------------------------------------*/
	int	param_to_uv_d2(double u, double v, Vector *d2u2, Vector *d2v2, Vector *d2uv); // added [2/7/2009 leon]

	int param_to_principal_curvature(double u, double v, double *kmax, double *kmin); // added [2/7/2009 leon]

	int param_to_d_EFG(double u, double v, 
					   double *du_E, double *du_F, double *du_G,
					   double *dv_E, double *dv_F, double *dv_G); // added [2/10/2009 leon]
	/* --------------------------------------------------------------------
	 * 排序
	 * 输入参数：
	 *   a    待排序数组 
	 *   num    数组大小
	 * 
	 * ------------------------------------------------------------------*/
	void sort(double a[], int num);


	/* --------------------------------------------------------------------
	 * 获取最近patch中心点
	 * 输入参数：
	 * coord 基点坐标
	 * ------------------------------------------------------------------*/
	void fguess(Vector coord,double uv[]);


	/* --------------------------------------------------------------------
	 * 在_duv方向上一维搜索，获取最小值
	 * 输入参数：
	 * _coord 基点坐标
	 * _uv 初始参数坐标
	 * _duv 搜索方向
	 * _ax，_bx 搜索范围
	 * func 待求最小值的函数
	 * ------------------------------------------------------------------*/
	double line_search(Vector _coord,double _uv[],double _duv[],
		               double _ax, double _bx,double func(double p[],int n),
					   double _tol, double &_xmin);
	
	/* --------------------------------------------------------------------
	 * 确定搜索方向的有效性
	 * 输入参数：
	 * _uv 初始参数坐标
	 * _duv 搜索方向
	 * _box 参数平面范围
	 * ------------------------------------------------------------------*/
	bool feasib(double _uv[], double _duv[],
		        double _box[][2], double _tol);


	/* ---------------------------------------------------------------------
	 * 获取新的搜索方向
	 * 输入参数：
	 */

	bool newduv(double _uv[], double _gk[], double _duv[], 
				double _box[][2], double _tol);

	/* --------------------------------------------------------------------
	 * 直接获取最近点
	 * 输入参数：
	 * _coord 基点坐标
	 * _uv 初始参数坐标
	 */
	void locu2(Vector _coord, double _uv[],Vector &_proj, 
			   double &_sq_dist, int &_flag);

	// added [1/4/2009 ly]
	// 计算以参数坐标给出的3个点的外心位置参数坐标
	void calc_riemannian_circumcenter(double u[3], double v[3], double *cu, double *cv);
	// added [2/15/2009 ly]
	void calc_riemannian_circumcenter(double u[3], double v[3], double *cu, double *cv, Mat2by2d metric);

	bool calc_circumcenter_on_plane(Vec2d pnt[3], Vec2d &cen);	// added [1/22/2009 leon]	

	int is_delaunay_broken(Vec2d verts[3], Vec2d cen, Vec2d pnt, double rad); // added [2/8/2009 leon]
	// 1: broken; -1: kept; 0: four points on same circle

	void get_surface_size(double *su, double *sv) { // added [1/28/2009 leon]
		*su = size_u;
		*sv = size_v;
	}

	void param_to_riemannian_metric(double u, double v, Mat2by2d mt); // added [2/10/2009 leon]
	
	double calc_riemannian_distance(Vec2d from, Vec2d to, int n); // added [2/10/2009 leon]
	double calc_riemannian_distance(double uf, double vf, double ut, double vt, int n); //  [2/13/2009 leon]

private:
	static double dist_to_proj_point_func(double uv[], int n);
//	static void dist_to_proj_point_dfunc(double uv[], double duv[], int n);

	// added [1/7/2009 ly]
	void calc_riemannian_distance(Mat1by1d val, Mat2by1d du, Mat2by2d mt);

protected:
	int num_u, num_v;
	Vector *d0_vector;		/* 型值点数组 */
	Vector *d1_u_vector;	/* 型值点u向切向量数组 */
	Vector *d1_v_vector;	/* 型值点u向切向量数组 */
	Vector *d2_uv_vector;	/* 型值点uv向二阶偏导数组 */

	// added [1/28/2009 leon]
	double size_u;			// maximum length of u curves, the length is approximated by polylines
	double size_v;			// maximum length of v curves, the length is approximated by polylines
	/* ------------------------------------------------------------------
	 * 根据朱心雄的《自由曲线曲面造型技术》，对于每个曲面片，其表达式为
	 * P(u,w) = UMBM'W'  (M'和W'表示矩阵M和W的转置)
	 * 其中
	 *     _                _
	 *    |  2   -2   1    1 |
	 * M= | -3    3  -2   -1 |
	 *    |  0    0   1    0 |
	 *    |  1    0   0    0 |
	 *     -                -
	 *    |   P(0,0)   P(0,1)   Pw(0,0)    Pw(0,1) |
	 * B= |   P(1,0)   P(1,1)   Pw(1,0)    Pw(1,1) |
	 *    |  Pu(0,0)  Pu(0,1)  Puw(0,0)   Puw(0,1) |
	 *    |  Pu(1,0)  Pu(1,1)  Puw(1,0)   Puw(1,1) |
	 *
	 * 记C = MBM', 则P(u,w) = UCW'，展开得：
	 * P(u,w) = c11*u3*v3 + c12*u3*v2 + c13*u3*v + c14*u3 +
	            c21*u2*v3 + c22*u2*v2 + c23*u2*v + c24*u2 +
				c31*u*v3  + c32*u*v2  + c33*u*v  + c34*u  +
				c41*v3    + c42*v2    + c43*v    + c44    
	 * patch_coeffs成员变量记录每个曲面片的上述参数，其size为：
	 * 16 * num_u * num_v
	 *-------------------------------------------------------------------*/
	Vector *patch_coeffs;
#ifdef _DEBUG
	Vector *bound_info;
#endif

	static FergusonSurface *myself_object;	  
	static Vector proj_point;				  
};

}/* namespace EEMAS */

#endif /* __eemas_ferguson_surface_h__ */
