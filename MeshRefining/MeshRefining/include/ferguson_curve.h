/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2008��05��19��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, May. 18, 2007
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953166
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��ferguson_curve.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ�����Ferguson���ߵ�ͨ�ýӿ�
 *     Define a set of interface for Ferguson curve
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __eemas_ferguson_curve_h__
#define __eemas_ferguson_curve_h__

//#include "vector.h"
#include "curve.h"
#include <assert.h>
#include <memory.h>
#include <math.h>

namespace EEMAS {

class FergusonSurface;

/* �������ߵı߽����� */
enum EndCondition 
{
	SpecifiedD2 = 0,  // free form end
	SpecifiedD1,      // clipping form end
	Parabolic,        // parabolic form end
};

class FergusonCurve : public Curve
{
public:
  int flag;
	/* -------------------------------------------------------------------
	 * ���캯������������
	 * -----------------------------------------------------------------*/
	FergusonCurve();
	virtual ~FergusonCurve();
	
	/* --------------------------------------------------------------------
	 * ��ʼ��һ��Furgeson����
	 * �������
	 *     data_points, npnts			��ֵ������
	 *     e_cond						�˵�����
	 *     v0, vn	                    �˵���������
	 * ����:
	 *     ֧��3��˵�����
	 *     1.   ָ����ĩ��Ķ��׵���(�洢��v0/vn)
	 *     2.   ָ����ĩ���һ�׵���(�洢��v0/vn)
	 *     3.   �����߶˵�����: ��ĩ���߶�Ϊ������
	 *
	 *     ȷʡ״̬�£��˵�����Ϊ����1�������˶��׵���Ϊ0 
	 * ------------------------------------------------------------------*/
	int initialize(
		const Vector data_points[], int npnts,		  /* ��ֵ������ */
		EndCondition e_cond = SpecifiedD2,			  /* ָ���˵����� */
		Vector v0 = Vector(0.0, 0.0, 0.0),		      /* �׵���� */
		Vector vn = Vector(0.0, 0.0, 0.0)             /* ĩ����� */);

	/* -------------------------------------------------------------------
	 * ��ȡ�˵�����
	 * -----------------------------------------------------------------*/
	Vector start_coord();
	Vector end_coord();
    void get_box(double x[2],double y[2],double z[2]);

	/* --------------------------------------------------------------------
	 * ͨ������ֵ��ȡ��Ӧ�������
	 * ���������
	 *   u			��������
	 *   coord		�������(����ֵ)
	 * ����ֵ��
	 *   �������
	 * ������
	 * -------------------------------------------------------------------*/
	virtual int	param_to_coord(double u, Vector *coord);


    /*get the dx,dy,dz by the u
     *
     *
     * */
    virtual int get_d_coord(Vector coord,Vector &dcoord);

	/* --------------------------------------------------------------------
	 * ��ȡ�ܻ���
	 * -------------------------------------------------------------------*/
	virtual double get_total_arc_length();

	/* --------------------------------------------------------------------
	 * ͨ���������������һ������
	 * -------------------------------------------------------------------*/
	virtual int	 arc_length_to_coord(
								double al,              /* ���� */
								Vector *coord,		    /* ������� */
								double tol = 1.0e-4,    /* ������� */
								int max_iter = 20	    /* ���������� */
								);

	/* -------------------------------------------------------------------------------
	 * ��ȡĳ���������ϵ�ͶӰ�㣨����㣩
	 * ���������
	 *     coord    ����
	 *     cls_pnt  ͶӰ��(����ֵ)
	 *     cls_u	ͶӰ��Ĳ���ֵ(����ֵ)
	 * ����ֵ��
	 *     �������
	 * ������
	 *         
	 * -----------------------------------------------------------------------------*/
	virtual int project(
		Vector coord,			/* ���� */
		Vector *cls_pnt,		/* ����� */
		double *cls_u			/* �����Ĳ���ֵ */
		);

	/*------------------------------------------------------------------------------
	 *��ȡ������ֵ������������
	 *���������
	 *          d1_vector_temp  ��ֵ������������
	 *          npts       ��ֲ�����
	 *����ֲ��
	 *          
	 *������
	 *-------------------------------------------------------------------------------*/
	void get_d1_vector(Vector* d1_vector_temp, int npts)
	{
		int i;

		for(i=0; i < data_point_num && i < npts; i++)
		{
			d1_vector_temp[i] = d1_vector[i];
		}
	}
	/*------------------------------------------------------------------------------
	 *��ȡ���߻�������
	 *���������
	 *          arc_length ��������
	 *          npts       ��ֲ�����
	 *����ֲ��
	 *          
	 *������
	 *-------------------------------------------------------------------------------*/
	void get_arc_length(double* arc_length_temp, int npts)
	{
		int i;

		for(i=0; i < data_point_num && i < npts; i++)
		{
			arc_length_temp[i] = arc_length[i];
		}
	}
	
protected:

	/* ---------------------------------------------------------------------
	 * �ͷ�������Դ�����ָ�����ʼ״̬
	 * -------------------------------------------------------------------*/
	virtual void destroy();

	/* ---------------------------------------------------------------------
	 * ����׷�Ϸ����Ferguson���߸���ֵ�㴦��������ֵ
	 * ��֪����:
	 *     end_condition		    �˵�����
	 *     _v0, _vn		           �˵���������
	 * ����:
	 *     ֧��3��˵�����
	 *     1.   ָ����ĩ��Ķ��׵���(�洢��_v0/_vn)
	 *     2.   ָ����ĩ���һ�׵���(�洢��_v0/_vn)
	 *     3.   �����߶˵�����: ��ĩ���߶�Ϊ������
	 *
	 *     ȷʡ״̬�£��˵�����Ϊ����1�������˶��׵���Ϊ0
	 * --------------------------------------------------------------------*/
	void calc_d1_vector();

	/* ---------------------------------------------------------------------
	 * ����ĳ�������߶��ϵĲ���ֵ��ȡ���Ӧ�Ļ���
	 * �������:
	 *		iseg,			�ֶ�������� 
	 *      u			    ���صĲ�������ֵ 
	 *		sal				���ڷֶ������ϵĻ���(����ֵ)
	 *      tol				�������
	 *      max_iter		����������
	 * ����ֵ
	 *      �������
	 * ����:
	 * 
	 * --------------------------------------------------------------------*/
	 int seg_param_to_arc_length(
					int iseg,			/* �ֶ�������� */
					double u,			/* ���صĲ�������ֵ */
					double *sal,		/* ���ڷֶ������ϵĻ��� */
					double tol = 1.0e-4,/* ������� */
					int max_iter = 20	/* ���������� */
					);

	/* ---------------------------------------------------------------------
	 * ����ĳ�������߶��ϵĻ���������Ӧ�Ĳ�������ֵ
	 * �������:
	 *		iseg,			�ֶ�������� 
	 *		sal				���ڷֶ������ϵĻ���
	 *      u			    ���صĲ�������ֵ 
	 *      tol				�������
	 *      max_iter		����������
	 * ����ֵ
	 *      �������
	 * ����:
	 * 
	 * --------------------------------------------------------------------*/
	 int seg_arc_length_to_param(
					int iseg,			/* �ֶ�������� */
					double sal,			/* ���ڷֶ������ϵĻ��� */
					double *u,			/* ���صĲ�������ֵ */
					double tol = 1.0e-4,/* ������� */
					int max_iter = 20	/* ���������� */
					);

	 /* ------------------------------------------------------------------
	  * ���ݵ��ڲ������߶��ϵĲ�������ֵȷ������ά����ֵ
	  * ���������
	  *		iseg		�ֶ�������� 
	  *		u			��������ֵ 
	  *		coord		���ص���ά����ֵ
	  * ����ֵ��
	  *     �������
	  * ������
	  *
	  * ----------------------------------------------------------------*/
	int seg_param_to_coord(
					int iseg,			/* �ֶ�������� */
					double u,			/* ��������ֵ */
					Vector *coord		/* ���ص���ά����ֵ */
					);

    /*get the seg_param_to_d1coord
     *
     *
     * */
    int seg_param_to_d1coord(int iseg, double u, Vector *d1coord);

protected:

#if 0
	/* ---------------------------------------------------------------------
	 * ���õ�ǰ���߶α��
	 * ���������
	 *     iseg    ��ǰ���߶ε����(����0)
	 * ����ֵ��
	 *     �������
	 * ������
	 *        
	 * -------------------------------------------------------------------*/
	int set_current_segment(int iseg)
	{
		if (iseg < 0 || iseg >= data_point_num)
		{
			return 2; /* error arguments */
		}

		current_segment = iseg;
	}

	/* --------------------------------------------------------------------------
	 * ����ָ������ֵ(specified_seg_arc_length)
	 * ���������
	 *     sal       ָ������ֵ
	 * ����ֵ��
	 *     �������
	 * ������
	 *     ָ������������specified_seg_arc_length��
	 * --------------------------------------------------------------------------*/
	static int specify_seg_arc_length(double sal)
	{
		if (sal < 0.0)
			return 2; /* error argument */

		specified_seg_arc_length = sal;
		return 0;
	}

	/* --------------------------------------------------------------------------
	 * ��Ӧ��ǰ���߶�, �������ֵΪu�ĵ�Ļ�����ָ������(specified_seg_arc_length)
	 * �Ĳ�ֵ
	 * ���������
	 *     u       ����ֵ
	 * ����ֵ��
	 *     ������ֵ
	 * ������
	 *     ָ������������specified_seg_arc_length��
	 * --------------------------------------------------------------------------*/
	static double seg_arc_length_subtract(double u)
	{
		Vector cvpnt;
		assert(u >= 0.0 && u <= data_point_num - 1);

		param_to_coord(u, &cvpnt);

		return (cvpnt.x - base_point.x)*(cvpnt.x - base_point.x)
			   (cvpnt.y - base_point.y)*(cvpnt.y - base_point.y)
			   (cvpnt.z - base_point.z)*(cvpnt.z - base_point.z);
	}

#endif 

	/* ---------------------------------------------------------------------
	 * ���õ�ǰ���߶�������ģֵ��ϵ������
	 * ���������
	 *     iseg    ��ǰ���߶ε����(����0)
	 * ����ֵ��
	 *     �������
	 * ������
	 *     �ڵ���curve_seg_d1_module֮ǰ����Ҫ���øú������ö�Ӧ��ϵ��
	 * -------------------------------------------------------------------*/
    static int set_cur_d1_modue_coeff(FergusonCurve * curve, int iseg)
	{
		assert(curve && curve->cv_seg_d1_module_coeff);
		memcpy(cur_seg_d1_module_coeff, &(curve->cv_seg_d1_module_coeff[5*iseg]), 
			5 * sizeof(cur_seg_d1_module_coeff[0]));

		return 1;
	}

	/* --------------------------------------------------------------------
	 * ���㵱ǰ���߶β���ֵΪu�ĵ��������ģֵ
	 * ���������
	 *     u       ����ֵ
	 * ����ֵ��
	 *     ������ģֵ
	 * ������
	 *     ���߶ε�����������ģϵ����������ľ�̬����cur_seg_d1_coeff��
	 * --------------------------------------------------------------------*/
	static double curve_seg_d1_module(double u)
	{
		/* ---------------------------------------------------------
		 * f = sqrt(a4 * t**4 + a3 * t**3 + a2 * t**2 + a1 * t + a0)
		 *                     ________________________
		 * cur_seg_d1_coeff = |a0 | a1 | a2 | a3 | a4 |
		 *                     ------------------------
		 * --------------------------------------------------------*/
		double u_2 = u * u, u_3 = u_2 * u, u_4 = u_3 * u;
		double square = 0.0;

		square = 
			cur_seg_d1_module_coeff[0] + cur_seg_d1_module_coeff[1] * u + 
			cur_seg_d1_module_coeff[2] * u_2 + cur_seg_d1_module_coeff[3] * u_3 +
			cur_seg_d1_module_coeff[4] * u_4;
		
		if (square < 0.0)
			square = 0.0;
		//assert(square >= 0.0);

		return sqrt(square);
	}
	
	/* --------------------------------------------------------------------
	 * �������Ͳ���Ϊu�����ߵ�֮��ľ���ƽ��
	 * ���������
	 *     uc		���ߵ����ֵ
	 * ����ֵ��
	 *     ����ֵ   
	 * ������
	 *     ������Ϣ������base_point�У����ߵ���Ϣ
	 * --------------------------------------------------------------------*/
	static double calc_pnt_curve_sq_dist(const double uc)
	{
		Vector cpnt;

		assert(base_curve);
		base_curve->param_to_coord(uc, &cpnt);

		return (cpnt.x - base_point.x)*(cpnt.x - base_point.x) +
			   (cpnt.y - base_point.y)*(cpnt.y - base_point.y) +
			   (cpnt.z - base_point.z)*(cpnt.z - base_point.z);
	}

protected:
	Vector *d0_vector;		/* ��ֵ������ */
	Vector *d1_vector;		/* ��ֵ������������ */
	double *arc_length;		/* ��������, ǰdata_point_num-1����ԱΪ�����ֶ����ߵĻ�����
	                           ���һ����ԱΪ�ܻ��� */
	int data_point_num;		/* ��ֵ����Ŀ */
	
	EndCondition end_condition;	/* �˵����� */
	Vector _v0, _vn;			/* ��ĩ��Ķ˵���������ֵ */

	Vector *cv_seg_d0_coeff;	/* �������߶ε�λ������ϵ�� size = data_point_num * 4 */
	Vector *cv_seg_d1_coeff;    /* �������߶ε�������ϵ�� size = data_point_num * 3 */
	Vector *cv_seg_d2_coeff;    /* �������߶εĶ�������ϵ�� size = data_point_num * 2 */

	double *cv_seg_d1_module_coeff;		       /* �������߶ε�������ģ������ϵ������ size = data_point_num * 5 */
	static double cur_seg_d1_module_coeff[5];  /* ��ǰ���߶ε�������ģ������ϵ������ */
	
	static FergusonCurve *base_curve;		  /* �����ߣ����ھ�̬��������calc_pnt_curve_sq_dist */
	static Vector base_point;				  /* �����ߣ����ھ�̬��������calc_pnt_curve_sq_dist */
#if 0
	static double specified_seg_arc_length;	   /* �;�̬����seg_arc_length_func��Ӧ */
#endif

	friend class FergusonSurface;
};


								

} /* namespace EEMAS */

#endif /* __eemas_ferguson_curve_h__ */
