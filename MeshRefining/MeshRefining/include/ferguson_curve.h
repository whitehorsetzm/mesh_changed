/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2008年05月19日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, May. 18, 2007
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953166
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：ferguson_curve.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对Ferguson曲线的通用接口
 *     Define a set of interface for Ferguson curve
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
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

/* 样条曲线的边界条件 */
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
	 * 构造函数和析构函数
	 * -----------------------------------------------------------------*/
	FergusonCurve();
	virtual ~FergusonCurve();
	
	/* --------------------------------------------------------------------
	 * 初始化一条Furgeson曲线
	 * 输入参数
	 *     data_points, npnts			型值点数组
	 *     e_cond						端点条件
	 *     v0, vn	                    端点条件参数
	 * 描述:
	 *     支持3类端点条件
	 *     1.   指定首末点的二阶导数(存储在v0/vn)
	 *     2.   指定首末点的一阶导数(存储在v0/vn)
	 *     3.   抛物线端点条件: 首末曲线段为抛物线
	 *
	 *     确省状态下，端点条件为类型1，且两端二阶导数为0 
	 * ------------------------------------------------------------------*/
	int initialize(
		const Vector data_points[], int npnts,		  /* 型值点数组 */
		EndCondition e_cond = SpecifiedD2,			  /* 指明端点条件 */
		Vector v0 = Vector(0.0, 0.0, 0.0),		      /* 首点参数 */
		Vector vn = Vector(0.0, 0.0, 0.0)             /* 末点参数 */);

	/* -------------------------------------------------------------------
	 * 获取端点坐标
	 * -----------------------------------------------------------------*/
	Vector start_coord();
	Vector end_coord();
    void get_box(double x[2],double y[2],double z[2]);

	/* --------------------------------------------------------------------
	 * 通过参数值获取相应点坐标的
	 * 输入参数：
	 *   u			参数坐标
	 *   coord		点的坐标(返回值)
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	virtual int	param_to_coord(double u, Vector *coord);


    /*get the dx,dy,dz by the u
     *
     *
     * */
    virtual int get_d_coord(Vector coord,Vector &dcoord);

	/* --------------------------------------------------------------------
	 * 获取总弧长
	 * -------------------------------------------------------------------*/
	virtual double get_total_arc_length();

	/* --------------------------------------------------------------------
	 * 通过弧长求解曲线上一点坐标
	 * -------------------------------------------------------------------*/
	virtual int	 arc_length_to_coord(
								double al,              /* 弧长 */
								Vector *coord,		    /* 点的坐标 */
								double tol = 1.0e-4,    /* 容许误差 */
								int max_iter = 20	    /* 最大迭代次数 */
								);

	/* -------------------------------------------------------------------------------
	 * 获取某点在曲线上的投影点（最近点）
	 * 输入参数：
	 *     coord    基点
	 *     cls_pnt  投影点(返回值)
	 *     cls_u	投影点的参数值(返回值)
	 * 返回值：
	 *     错误代码
	 * 描述：
	 *         
	 * -----------------------------------------------------------------------------*/
	virtual int project(
		Vector coord,			/* 基点 */
		Vector *cls_pnt,		/* 最近点 */
		double *cls_u			/* 最近点的参数值 */
		);

	/*------------------------------------------------------------------------------
	 *获取曲线型值点切向量数组
	 *输入参数：
	 *          d1_vector_temp  型值点切向量数组
	 *          npts       型植点个数
	 *返回植：
	 *          
	 *描述：
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
	 *获取曲线弧长数组
	 *输入参数：
	 *          arc_length 弧长数组
	 *          npts       型植点个数
	 *返回植：
	 *          
	 *描述：
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
	 * 释放所有资源，并恢复到初始状态
	 * -------------------------------------------------------------------*/
	virtual void destroy();

	/* ---------------------------------------------------------------------
	 * 利用追赶法求解Ferguson曲线各型值点处的切向量值
	 * 已知条件:
	 *     end_condition		    端点条件
	 *     _v0, _vn		           端点条件参数
	 * 描述:
	 *     支持3类端点条件
	 *     1.   指定首末点的二阶导数(存储在_v0/_vn)
	 *     2.   指定首末点的一阶导数(存储在_v0/_vn)
	 *     3.   抛物线端点条件: 首末曲线段为抛物线
	 *
	 *     确省状态下，端点条件为类型1，且两端二阶导数为0
	 * --------------------------------------------------------------------*/
	void calc_d1_vector();

	/* ---------------------------------------------------------------------
	 * 利用某点在曲线段上的参数值获取其对应的弧长
	 * 输入参数:
	 *		iseg,			分段曲线序号 
	 *      u			    返回的参数坐标值 
	 *		sal				点在分段曲线上的弧长(返回值)
	 *      tol				容许误差
	 *      max_iter		最大迭代次数
	 * 返回值
	 *      错误代码
	 * 描述:
	 * 
	 * --------------------------------------------------------------------*/
	 int seg_param_to_arc_length(
					int iseg,			/* 分段曲线序号 */
					double u,			/* 返回的参数坐标值 */
					double *sal,		/* 点在分段曲线上的弧长 */
					double tol = 1.0e-4,/* 容许误差 */
					int max_iter = 20	/* 最大迭代次数 */
					);

	/* ---------------------------------------------------------------------
	 * 利用某点在曲线段上的弧长求解其对应的参数坐标值
	 * 输入参数:
	 *		iseg,			分段曲线序号 
	 *		sal				点在分段曲线上的弧长
	 *      u			    返回的参数坐标值 
	 *      tol				容许误差
	 *      max_iter		最大迭代次数
	 * 返回值
	 *      错误代码
	 * 描述:
	 * 
	 * --------------------------------------------------------------------*/
	 int seg_arc_length_to_param(
					int iseg,			/* 分段曲线序号 */
					double sal,			/* 点在分段曲线上的弧长 */
					double *u,			/* 返回的参数坐标值 */
					double tol = 1.0e-4,/* 容许误差 */
					int max_iter = 20	/* 最大迭代次数 */
					);

	 /* ------------------------------------------------------------------
	  * 根据点在参数曲线段上的参数坐标值确定其三维坐标值
	  * 输入参数：
	  *		iseg		分段曲线序号 
	  *		u			参数坐标值 
	  *		coord		返回的三维坐标值
	  * 返回值：
	  *     错误代码
	  * 描述：
	  *
	  * ----------------------------------------------------------------*/
	int seg_param_to_coord(
					int iseg,			/* 分段曲线序号 */
					double u,			/* 参数坐标值 */
					Vector *coord		/* 返回的三维坐标值 */
					);

    /*get the seg_param_to_d1coord
     *
     *
     * */
    int seg_param_to_d1coord(int iseg, double u, Vector *d1coord);

protected:

#if 0
	/* ---------------------------------------------------------------------
	 * 设置当前曲线段编号
	 * 输入参数：
	 *     iseg    当前曲线段的序号(基于0)
	 * 返回值：
	 *     错误代码
	 * 描述：
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
	 * 设置指定弧长值(specified_seg_arc_length)
	 * 输入参数：
	 *     sal       指定弧长值
	 * 返回值：
	 *     错误代码
	 * 描述：
	 *     指定弧长保存在specified_seg_arc_length中
	 * --------------------------------------------------------------------------*/
	static int specify_seg_arc_length(double sal)
	{
		if (sal < 0.0)
			return 2; /* error argument */

		specified_seg_arc_length = sal;
		return 0;
	}

	/* --------------------------------------------------------------------------
	 * 对应当前曲线段, 计算参数值为u的点的弧长和指定弧长(specified_seg_arc_length)
	 * 的差值
	 * 输入参数：
	 *     u       参数值
	 * 返回值：
	 *     弧长差值
	 * 描述：
	 *     指定弧长保存在specified_seg_arc_length中
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
	 * 设置当前曲线段切向量模值的系数数组
	 * 输入参数：
	 *     iseg    当前曲线段的序号(基于0)
	 * 返回值：
	 *     错误代码
	 * 描述：
	 *     在调用curve_seg_d1_module之前，需要调用该函数设置对应的系数
	 * -------------------------------------------------------------------*/
    static int set_cur_d1_modue_coeff(FergusonCurve * curve, int iseg)
	{
		assert(curve && curve->cv_seg_d1_module_coeff);
		memcpy(cur_seg_d1_module_coeff, &(curve->cv_seg_d1_module_coeff[5*iseg]), 
			5 * sizeof(cur_seg_d1_module_coeff[0]));

		return 1;
	}

	/* --------------------------------------------------------------------
	 * 计算当前曲线段参数值为u的点的切向量模值
	 * 输入参数：
	 *     u       参数值
	 * 返回值：
	 *     切向量模值
	 * 描述：
	 *     曲线段的切向量函数模系数保存在类的静态变量cur_seg_d1_coeff中
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
	 * 计算基点和参数为u的曲线点之间的距离平方
	 * 输入参数：
	 *     uc		曲线点参数值
	 * 返回值：
	 *     距离值   
	 * 描述：
	 *     基点信息保存在base_point中，曲线点信息
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
	Vector *d0_vector;		/* 型值点数组 */
	Vector *d1_vector;		/* 型值点切向量数组 */
	double *arc_length;		/* 弧长数组, 前data_point_num-1个成员为各个分段曲线的弧长，
	                           最后一个成员为总弧长 */
	int data_point_num;		/* 型值点数目 */
	
	EndCondition end_condition;	/* 端点条件 */
	Vector _v0, _vn;			/* 首末点的端点条件参数值 */

	Vector *cv_seg_d0_coeff;	/* 所有曲线段的位置向量系数 size = data_point_num * 4 */
	Vector *cv_seg_d1_coeff;    /* 所有曲线段的切向量系数 size = data_point_num * 3 */
	Vector *cv_seg_d2_coeff;    /* 所有曲线段的二阶向量系数 size = data_point_num * 2 */

	double *cv_seg_d1_module_coeff;		       /* 所有曲线段的切向量模函数的系数数组 size = data_point_num * 5 */
	static double cur_seg_d1_module_coeff[5];  /* 当前曲线段的切向量模函数的系数数组 */
	
	static FergusonCurve *base_curve;		  /* 基曲线，用于静态函数，如calc_pnt_curve_sq_dist */
	static Vector base_point;				  /* 基曲线，用于静态函数，如calc_pnt_curve_sq_dist */
#if 0
	static double specified_seg_arc_length;	   /* 和静态函数seg_arc_length_func对应 */
#endif

	friend class FergusonSurface;
};


								

} /* namespace EEMAS */

#endif /* __eemas_ferguson_curve_h__ */
