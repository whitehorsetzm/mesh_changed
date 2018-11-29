/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ά����ͬ��Delaunay���������� (�汾�ţ�0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2005��9��15��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * ��ϵ��ʽ
 *   �绰��+86-571-87953165
 *   ���棺+86-571-87953167
 *   ���䣺zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __iso3d_matrix_h__
#define __iso3d_matrix_h__
#include "vector.h"

namespace iso3d
{


/* 
 * ������ m = a
 * copy a matrix m = a
 */
void matrixCopy(double *m, const double *a, const int row=4, const int col=4); 

/* 
 * ������� m = a + b
 * add two matrices m = a + b
 */
void matrixAdd(double *m, const double *a, const double *b,
		   const int  row=4, const int  col=4); 

/* 
 * ������� m = a - b
 * subtract two matrices m = a - b
 */
void matrixSub(double *m, const double *a, const double *b,
		   const int  row=4, const int  col=4); 

/* 
 * ������� m = a * b ����col1 = row2ʱ��Ч
 * multiply two matrices m = a * b ( valid only while col1 = row2)
 */
bool matrixMultiply(double *m, const double *a, const double *b, 
			    const int  row1=4, const int  col1=4, 
			    const int  row2=4, const int  col2=4);

/* 
 * ����ת�� (m = transpose(a))
 * transpose a matrix (m = transpose(a))
 */
void matrixTranspose(double *m, const double *a, const int row=4, const int col=4); 

/* 
 * ����ת�� (m = transpose(m))
 * transpose a matrix (m = transpose(m))
 */
void matrixTranspose(double *m, const int row=4, const int col=4); 

/*
 * �������� m = invert(a)
 * true ��ʾ�ɹ���false ��ʾ������󲻴���
 * invert a matrix m = invert(a)
 * true: successful; false: no inverse matrix
 */
bool matrixInverse(double *m, const double *a, const int n=4);

/*
 * �������� m = invert(m)
 * true ��ʾ�ɹ���false ��ʾ������󲻴���
 * invert a matrix m = invert(m)
 * true: successful; false: no inverse matrix
 */
bool matrixInverse(double *m, const int n=4);

/*
 * �����һ�� 
 * identify a matrix
 */
void matrixIdentity(double *m, const int n=4);


/* v ���� m �任��Ϊ v', �������� v' */
Vector vectorTransform(const Vector& v, const double m[16]);
}
#endif /* __iso3d_matrix_h__ */
