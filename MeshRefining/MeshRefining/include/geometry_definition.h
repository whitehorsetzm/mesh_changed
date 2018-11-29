/************************************************************************/
/* define some geometry classes                                         */
/* first creator: ����ʤ__2007_11_30                                    */
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

/* �������ߵı߽����� */
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

//����
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

	void calculate_ferguson_curve_tangent_vector(	/* ����ferguson���߿��Ƶ��������ĺ��� */
		EndCondition endcondition,					/* ָ���˵����� */
		Vector _v0 = Vector(0.0,0.0,0.0),			/*  */
		Vector _vn = Vector(0.0,0.0,0.0)
		);

	void calculate_segment_coefficient();			/* ����ferguson���߸��ε�ϵ�� */

	void calculate_segment_length(double tol);		/* ��������Ƶ㵽��ʼ��Ļ��� */

	void calculate_partial_segment_length(			/* ����ferguson����ĳ����ָ������Ļ��� */
		const double *_coef,			/* �˶ε�ϵ�� */
		double _u0,						/* ��ʼ������ */
        double _u1,						/* ĩβ������ */
		double& _len,					/* ����õ�������ĳ��� */
		double _tol						/* �������� */
		);
	
	int discretization_with_psue(					/* ����psue�����еķ�����ɢ��ferguson���� */
		Vector **discretized_point,		/* �������յõ�����ɢ�����ά���� */
		int **disc_points_hints,
		GlobalSpacing *globalspacing	/* ����ı���������Ϣ */
		);

	// added [12/25/2008 ly]
	int sample_with_psue(					/* ����psue�����еķ������в��� */
		Vector **discretized_point,		/* ������������ά���� */
		int **disc_points_hints,
		double **size,					/* �����������������ʳߴ���Ϣ */	
		GlobalSpacing *globalspacing	/* ����ı���������Ϣ */
		);

	int discretization_with_proportional_2d(
		double** gcoord,				/* �����ɢ����������� */
		int n,							/* ָ���ķֶ��� */
		double q						/* �ȱ����еı�ֵ */
		);

	int discretization_with_proportional_3d(
		double** gcoord,				/* �����ɢ����������� */
		int n,							/* ָ���ķֶ��� */
		double q						/* �ȱ����еı�ֵ */
		);

	int get_control_point_number(){return ncp;}		/* ���ؿ��ƶ����� */

	int calculate_discretized_point_vector(			/* ����������ɢ������ */
		double *sample_arc_length,		/* �������㵽��ʼ��Ļ��� */
		int *sample_points_hints,
		double *sample_density,			/* �������㴦�ܶ� */
		double single_arc,				/* ���β������� */
		int &spn,						/* ���ز��������Ŀ */
		Vector **discretized_points,		/* ���յ���ɢ��� */
		int **disc_points_hints
		);

	void calculate_sample_arc_length(				/* ����������㵽��ʼ��Ļ��� */
		double **sample_arc_length,		/* ���ظ������㴦�Ļ��� */
		double &single_arc,				/* ���ز���ʱÿ�εĻ��� */
		int &spn,						/* ���ز�������Ŀ */
		double spmin					/* �����ȫ����С�ߴ� */
		);

	void calculate_sample_points_density(			/* ���������Ŀռ��ܶ� */
		int **sample_points_hints,
		double **sample_points_density,	/* ���صĲ������ܶ����� */
		double *sample_arc_length,		/* ��������㻡������ */
		int &spn,						/* �����������Ŀ */
		GlobalSpacing *globalspacing	/* ����������Ϣ */
		);

	// added by ly [12/24/2008 ly]
	void calculate_sample_points_coordinate_size(
		Vector **sample_points_coordinate, /* ���صĲ�������ά�������� */
		double **size,					/* ���صĲ�������������ʳߴ���Ϣ */
		int **disc_points_hints,
		double *sample_arc_length,		/* ��������㻡������ */
		int &spn,						/* �����������Ŀ */
		GlobalSpacing *globalspacing	/* ����������Ϣ */
		);

	void calculate_disnum_before_sam_points(		/* ����ÿ��������֮ǰӦ���ж�����ɢ�� */
		double *dis_num_before_sam_points,	/* ���ÿ��������֮ǰӦ���ж��ٸ���ɢ�� */
		int &dpn,						/* ���������ɢ��ĸ��� */
		double &average_size,			/* ���ƽ��ÿ������ɢ��֮��Ļ����ж��ٸ���ɢ�㣬*/
										 /* һ���ܽӽ�1������ */
		double single_arc,				/* �������ʱÿ�εĻ��� */
		int &spn,						/* �����������Ŀ */
		double *sample_points_density	/* ����ÿ�������㴦���ܶ� */
		);

	void calculate_3Dcoord_according_arc_length(	/* ���ݻ���������ά���� */
		double arc_length,				/* ���� */
		Vector &vec,					/* ��Ž�� */
		double tol						/* �������� */
		);

	void calculate_parameter_according_arc_length(	/* ��������ĳ���ϵ�һ�λ�����Ӧ����ֵ */
		double *_coef,					/* ĳ�ε�ϵ�� */
		double seg_length,				/* ĳ�ε��ܻ��� */
		double arc_length,				/* ����Ļ��� */
		double _u0,						/* ����Ļ�������ʼ�˲��� */
		double &_u1,					/* ��õĲ���ֵ */
		double tol=1.0e-6				/* �������� */
		);

	void get_3D_coord(								/* ���ݲ�����ȡ��ά���� */
		Vector &vec,					/* ���ص����ά���� */
		int i,							/* ָ���ڼ��� */
		double u1						/* ָ������ֵ */
		);

	// added [1/2/2009 ly]
	void get_D1_vector(							/* ���ݲ�����ȡһ�׵�ʸ */
		Vector &vec,					/* ���ص��һ�׵�ʸ */
		int i,							/* ָ���ڼ��� */
		double u1						/* ָ������ֵ */
		); 

	// added [1/2/2009 ly]
	void get_D2_vector(							/* ���ݲ�����ȡ���׵�ʸ */
		Vector &vec,					/* ���ص�Ķ��׵�ʸ */
		int i,							/* ָ���ڼ��� */
		double u1						/* ָ������ֵ */
		); 

	// added [1/2/2009 ly]
	double get_Curvature(							/* ���ݲ�����ȡ���� */
		int i,							/* ָ���ڼ��� */
		double u1						/* ָ������ֵ */
		); 	

protected:
	Vector *control_point;	//3-D coordinates of control points
	Vector *tangent_vector;	/* ���Ƶ��ϵ������� */
	int ncp;				//number of control points
	FergusonCoef *coef;		/* Fergusonϵ�������� */
	double *arc_length;		/* ���߸��εĻ��� */
};

#endif