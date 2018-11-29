/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2008��06��24��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, June. 24, 2008
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953166
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��ferguson_surface.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ��������ߵ�ͨ�ýӿ�
 *     Define a set of interface for Ferguson surfaces
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
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
	 * ����/��������
	 * -------------------------------------------------------------------*/
	FergusonSurface();
	~FergusonSurface();

	/* --------------------------------------------------------------------
	 * �ͷſռ�
	 * -------------------------------------------------------------------*/
	int destroy();

	/* --------------------------------------------------------------------
	 * ��ʼ��һ��Furgeson����
	 * �������
	 *     data_points, nu, nv			��ֵ������
	 * ����:
	 *     ֧�����ɶ˵�����
	 * ------------------------------------------------------------------*/
	int initialize(
	double data_points[], int nu, int nv		  /* ��ֵ������ */);

	/* --------------------------------------------------------------------
	 * ��ʼ��һ��Furgeson����
	 * �������
	 *     data_points, nu, nv			��ֵ������
	 * ����:
	 *     ֧�����ɶ˵�����
	 * ------------------------------------------------------------------*/
	int initialize(Vector data_points[], int nu, int nv		  /* ��ֵ������ */);
     void get_box(double x[2],double y[2],double z[2]);
	/* --------------------------------------------------------------------
	 * ͨ���������ȡ��Ӧ�Ĳ���ֵ
	 * ���������
	 *   coord		�������
	 *   u, v		��������(����ֵ)
	 *   tol		�������������������ֵ��С�ڴ�ֵ����Ϊ����������
	 * ����ֵ��
	 *   �������
	 * ������
	 * -------------------------------------------------------------------*/
	virtual int coord_to_param(Vector coord,
				double *u, double *v, 
				double tol = 1.0e-4);

	/* --------------------------------------------------------------------
	 * ͨ������ֵ��ȡ��Ӧ�������
	 * ���������
	 *   u, v		��������
	 *   coord		�������(����ֵ)
	 * ����ֵ��
	 *   �������
	 * ������
	 * -------------------------------------------------------------------*/
	int	param_to_coord(double u, double v, Vector *coord);

	int param_to_plane_coord(double u, double v, double *x, double *y);

	/* --------------------------------------------------------------------
	 * ��ȡĳ���������ϵ�ͶӰ��
	 * ���������
	 *   coord		�������
	 *   u, v		�����(����)�Ĳ�������(����ֵ)
	 *   tol		�������������������ֵ��С�ڴ�ֵ����Ϊ����������
	 * ����ֵ��
	 *   �������
	 * ������
	 * -------------------------------------------------------------------*/
	virtual int project(Vector coord, double *u, double *v, double tol = 1.0e-4);
	
	virtual int projectSTNB(Vector coord, double *u, double *v);

	virtual int projectSTNB(int iN, Vector *pCoord, double *&pOutUV);
	
	/* --------------------------------------------------------------------
	 * ��ȡһ��ƫ����
	 * ���������
	 *   u, v		��������
	 *   d1u, d1v	u��ƫ��/v��ƫ��
	 * ����ֵ��
	 *   �������
	 * ������
	 * -------------------------------------------------------------------*/
	int	param_to_uv_d1(double u, double v, Vector *d1u, Vector *d1v);

	/* --------------------------------------------------------------------
	 * ��ȡ����ƫ����
	 * ���������
	 *   u, v		��������
	 *   d2uv		����ƫ��
	 * ����ֵ��
	 *   �������
	 * ������
	 * ------------------------------------------------------------------*/
	int	param_to_uv_d2(double u, double v, Vector *d2uv);

	/* --------------------------------------------------------------------
	 * ��ȡ����ƫ����
	 * ���������
	 *   u, v		��������
	 *   d2u2
	 *   d2v2
	 *   d2uv		����ƫ��
	 * ����ֵ��
	 *   �������
	 * ������
	 * ------------------------------------------------------------------*/
	int	param_to_uv_d2(double u, double v, Vector *d2u2, Vector *d2v2, Vector *d2uv); // added [2/7/2009 leon]

	int param_to_principal_curvature(double u, double v, double *kmax, double *kmin); // added [2/7/2009 leon]

	int param_to_d_EFG(double u, double v, 
					   double *du_E, double *du_F, double *du_G,
					   double *dv_E, double *dv_F, double *dv_G); // added [2/10/2009 leon]
	/* --------------------------------------------------------------------
	 * ����
	 * ���������
	 *   a    ���������� 
	 *   num    �����С
	 * 
	 * ------------------------------------------------------------------*/
	void sort(double a[], int num);


	/* --------------------------------------------------------------------
	 * ��ȡ���patch���ĵ�
	 * ���������
	 * coord ��������
	 * ------------------------------------------------------------------*/
	void fguess(Vector coord,double uv[]);


	/* --------------------------------------------------------------------
	 * ��_duv������һά��������ȡ��Сֵ
	 * ���������
	 * _coord ��������
	 * _uv ��ʼ��������
	 * _duv ��������
	 * _ax��_bx ������Χ
	 * func ������Сֵ�ĺ���
	 * ------------------------------------------------------------------*/
	double line_search(Vector _coord,double _uv[],double _duv[],
		               double _ax, double _bx,double func(double p[],int n),
					   double _tol, double &_xmin);
	
	/* --------------------------------------------------------------------
	 * ȷ�������������Ч��
	 * ���������
	 * _uv ��ʼ��������
	 * _duv ��������
	 * _box ����ƽ�淶Χ
	 * ------------------------------------------------------------------*/
	bool feasib(double _uv[], double _duv[],
		        double _box[][2], double _tol);


	/* ---------------------------------------------------------------------
	 * ��ȡ�µ���������
	 * ���������
	 */

	bool newduv(double _uv[], double _gk[], double _duv[], 
				double _box[][2], double _tol);

	/* --------------------------------------------------------------------
	 * ֱ�ӻ�ȡ�����
	 * ���������
	 * _coord ��������
	 * _uv ��ʼ��������
	 */
	void locu2(Vector _coord, double _uv[],Vector &_proj, 
			   double &_sq_dist, int &_flag);

	// added [1/4/2009 ly]
	// �����Բ������������3���������λ�ò�������
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
	Vector *d0_vector;		/* ��ֵ������ */
	Vector *d1_u_vector;	/* ��ֵ��u������������ */
	Vector *d1_v_vector;	/* ��ֵ��u������������ */
	Vector *d2_uv_vector;	/* ��ֵ��uv�����ƫ������ */

	// added [1/28/2009 leon]
	double size_u;			// maximum length of u curves, the length is approximated by polylines
	double size_v;			// maximum length of v curves, the length is approximated by polylines
	/* ------------------------------------------------------------------
	 * ���������۵ġ����������������ͼ�����������ÿ������Ƭ������ʽΪ
	 * P(u,w) = UMBM'W'  (M'��W'��ʾ����M��W��ת��)
	 * ����
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
	 * ��C = MBM', ��P(u,w) = UCW'��չ���ã�
	 * P(u,w) = c11*u3*v3 + c12*u3*v2 + c13*u3*v + c14*u3 +
	            c21*u2*v3 + c22*u2*v2 + c23*u2*v + c24*u2 +
				c31*u*v3  + c32*u*v2  + c33*u*v  + c34*u  +
				c41*v3    + c42*v2    + c43*v    + c44    
	 * patch_coeffs��Ա������¼ÿ������Ƭ��������������sizeΪ��
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
