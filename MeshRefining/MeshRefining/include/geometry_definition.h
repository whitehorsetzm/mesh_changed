/************************************************************************/
/* define some geometry classes                                         */
/* first creator: 孙力胜__2007_11_30                                    */
/************************************************************************/
#ifndef GEOMETRY_DEFINITION_H
#define GEOMETRY_DEFINITION_H

#include "vector.h"
#include "spacing.h"
#include <iostream>
#include "global.h"
#include <math.h>

//#pragma comment (lib, "dfor")

typedef struct FergusonCoef
{
	double coef[5];
}FergusonCoef;

/* 样条曲线的边界条件 */
enum EndCondition 
{
	ZeroSecondDerivative = 0,  // free form end
	Parabolic,                 // parabolic form end
	SpecifiedTangent           // clipping form end
};

// extern "C" {
// 	extern void __stdcall CFUNC(double v[3]);
// }


enum DerivativeNumber
{
    Zero = 0,
	One   ,
	Two
};

//父类
class CCurve
{

};

class FergusonCurve : public CCurve
{
public:
	FergusonCurve(int n,Vector *vec)
	{
		int *iTmp = NULL;
		ncp = n;
		Vector vec0(0.0,0.0,0.0);
		control_point = new Vector[ncp];
		tangent_vector = new Vector[ncp];
//		printf( "size:%d\n", sizeof(FergusonCoef) );
//		assert( ncp > 1 );
		coef = new FergusonCoef[ncp-1];
//		iTmp = new int[ncp*10];
//		delete []iTmp;
		arc_length = new double[ncp];
		for(int i = 0; i<ncp; i++)
		{
			control_point[i] = vec[i];
			tangent_vector[i] =vec0;
		}
	}

	~FergusonCurve()
	{
		delete[] control_point;
		delete[] tangent_vector;
		delete[] coef;
		delete[] arc_length;
	}

	void calculate_ferguson_curve_tangent_vector(	/* 计算ferguson曲线控制点切向量的函数 */
		EndCondition endcondition,					/* 指明端点条件 */
		Vector _v0 = Vector(0.0,0.0,0.0),			/*  */
		Vector _vn = Vector(0.0,0.0,0.0)
		);

	void calculate_segment_coefficient();			/* 计算ferguson曲线各段的系数 */

	void calculate_segment_length(double tol);		/* 计算各控制点到起始点的弧长 */

	void calculate_partial_segment_length(			/* 计算ferguson曲线某段中指定区间的弧长 */
		const double *_coef,			/* 此段的系数 */
		double _u0,						/* 起始处参数 */
        double _u1,						/* 末尾处参数 */
		double& _len,					/* 计算得到的区间的长度 */
		double _tol						/* 允许的误差 */
		);
	
	int discretization_with_psue(					/* 根据psue程序中的方法离散此ferguson曲线 */
		Vector **discretized_point,		/* 保存最终得到的离散点的三维坐标 */
		int **disc_points_hints,
		GlobalSpacing *globalspacing	/* 输入的背景网格信息 */
		);

	// added [12/25/2008 ly]
	int sample_with_psue(					/* 根据psue程序中的方法进行采样 */
		Vector **discretized_point,		/* 保存采样点的三维坐标 */
		int **disc_points_hints,
		double **size,					/* 保存采样点的理想曲率尺寸信息 */	
		GlobalSpacing *globalspacing	/* 输入的背景网格信息 */
		);

	int discretization_with_proportional_2d(
		double** gcoord,				/* 存放离散点坐标的数组 */
		int n,							/* 指定的分段数 */
		double q						/* 等比数列的比值 */
		);

	int discretization_with_proportional_3d(
		double** gcoord,				/* 存放离散点坐标的数组 */
		int n,							/* 指定的分段数 */
		double q						/* 等比数列的比值 */
		);

	int get_control_point_number(){return ncp;}		/* 返回控制顶点数 */

	int calculate_discretized_point_vector(			/* 计算曲线离散点向量 */
		double *sample_arc_length,		/* 各采样点到起始点的弧长 */
		int *sample_points_hints,
		double *sample_density,			/* 各采样点处密度 */
		double single_arc,				/* 单段采样弧长 */
		int &spn,						/* 返回采样点的数目 */
		Vector **discretized_points,		/* 最终的离散结果 */
		int **disc_points_hints
		);

	void calculate_sample_arc_length(				/* 计算各采样点到起始点的弧长 */
		double **sample_arc_length,		/* 返回各采样点处的弧长 */
		double &single_arc,				/* 返回采样时每段的弧长 */
		int &spn,						/* 返回采样点数目 */
		double spmin					/* 输入的全局最小尺寸 */
		);

	void calculate_sample_points_density(			/* 计算采样点的空间密度 */
		int **sample_points_hints,
		double **sample_points_density,	/* 返回的采样点密度数组 */
		double *sample_arc_length,		/* 输入采样点弧长数组 */
		int &spn,						/* 输入采样点数目 */
		GlobalSpacing *globalspacing	/* 背景网格信息 */
		);

	// added by ly [12/24/2008 ly]
	void calculate_sample_points_coordinate_size(
		Vector **sample_points_coordinate, /* 返回的采样点三维坐标数组 */
		double **size,					/* 返回的采样点的理想曲率尺寸信息 */
		int **disc_points_hints,
		double *sample_arc_length,		/* 输入采样点弧长数组 */
		int &spn,						/* 输入采样点数目 */
		GlobalSpacing *globalspacing	/* 背景网格信息 */
		);

	void calculate_disnum_before_sam_points(		/* 计算每个采样点之前应该有多少离散点 */
		double *dis_num_before_sam_points,	/* 输出每个采样点之前应该有多少个离散点 */
		int &dpn,						/* 输出最终离散点的个数 */
		double &average_size,			/* 输出平均每两个离散点之间的弧上有多少个离散点，*/
										 /* 一个很接近1的数字 */
		double single_arc,				/* 输入采样时每段的弧长 */
		int &spn,						/* 输入采样点数目 */
		double *sample_points_density	/* 输入每个采样点处的密度 */
		);

	void calculate_3Dcoord_according_arc_length(	/* 根据弧长计算三维坐标 */
		double arc_length,				/* 弧长 */
		Vector &vec,					/* 存放结果 */
		double tol						/* 允许的误差 */
		);

	void calculate_parameter_according_arc_length(	/* 根据曲线某段上的一段弧求相应参数值 */
		double *_coef,					/* 某段的系数 */
		double seg_length,				/* 某段的总弧长 */
		double arc_length,				/* 输入的弧长 */
		double _u0,						/* 输入的弧长的起始端参数 */
		double &_u1,					/* 求得的参数值 */
		double tol=1.0e-6				/* 允许的误差 */
		);

	void get_3D_coord(								/* 根据参数求取三维坐标 */
		Vector &vec,					/* 返回点的三维坐标 */
		int i,							/* 指明第几段 */
		double u1						/* 指明参数值 */
		);

	// added [1/2/2009 ly]
	void get_D1_vector(							/* 根据参数求取一阶导矢 */
		Vector &vec,					/* 返回点的一阶导矢 */
		int i,							/* 指明第几段 */
		double u1						/* 指明参数值 */
		); 

	// added [1/2/2009 ly]
	void get_D2_vector(							/* 根据参数求取二阶导矢 */
		Vector &vec,					/* 返回点的二阶导矢 */
		int i,							/* 指明第几段 */
		double u1						/* 指明参数值 */
		); 

	// added [1/2/2009 ly]
	double get_Curvature(							/* 根据参数求取曲率 */
		int i,							/* 指明第几段 */
		double u1						/* 指明参数值 */
		); 	

protected:
	Vector *control_point;	//3-D coordinates of control points
	Vector *tangent_vector;	/* 控制点上的切向量 */
	int ncp;				//number of control points
	FergusonCoef *coef;		/* Ferguson系数的数组 */
	double *arc_length;		/* 曲线各段的弧长 */
};

#endif