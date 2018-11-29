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
#ifndef __iso3d_matrix_h__
#define __iso3d_matrix_h__
#include "vector.h"

namespace iso3d
{


/* 
 * 矩阵复制 m = a
 * copy a matrix m = a
 */
void matrixCopy(double *m, const double *a, const int row=4, const int col=4); 

/* 
 * 矩阵相加 m = a + b
 * add two matrices m = a + b
 */
void matrixAdd(double *m, const double *a, const double *b,
		   const int  row=4, const int  col=4); 

/* 
 * 矩阵相减 m = a - b
 * subtract two matrices m = a - b
 */
void matrixSub(double *m, const double *a, const double *b,
		   const int  row=4, const int  col=4); 

/* 
 * 矩阵相乘 m = a * b 仅当col1 = row2时有效
 * multiply two matrices m = a * b ( valid only while col1 = row2)
 */
bool matrixMultiply(double *m, const double *a, const double *b, 
			    const int  row1=4, const int  col1=4, 
			    const int  row2=4, const int  col2=4);

/* 
 * 矩阵转置 (m = transpose(a))
 * transpose a matrix (m = transpose(a))
 */
void matrixTranspose(double *m, const double *a, const int row=4, const int col=4); 

/* 
 * 矩阵转置 (m = transpose(m))
 * transpose a matrix (m = transpose(m))
 */
void matrixTranspose(double *m, const int row=4, const int col=4); 

/*
 * 矩阵求逆 m = invert(a)
 * true 表示成功，false 表示该逆矩阵不存在
 * invert a matrix m = invert(a)
 * true: successful; false: no inverse matrix
 */
bool matrixInverse(double *m, const double *a, const int n=4);

/*
 * 矩阵求逆 m = invert(m)
 * true 表示成功，false 表示该逆矩阵不存在
 * invert a matrix m = invert(m)
 * true: successful; false: no inverse matrix
 */
bool matrixInverse(double *m, const int n=4);

/*
 * 矩阵归一化 
 * identify a matrix
 */
void matrixIdentity(double *m, const int n=4);


/* v 经过 m 变换后为 v', 函数返回 v' */
Vector vectorTransform(const Vector& v, const double m[16]);
}
#endif /* __iso3d_matrix_h__ */
