/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2008年06月30日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, June. 30, 2008
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953166
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：math_util.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对曲线的通用接口
 *     Define a set of interface for math utilities
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/


#ifndef __eemas_math_util_h__
#define __eemas_math_util_h__

#include <math.h>

namespace EEMAS {
/* -------------------------------------------------------------------------------
 * 矩阵运算
 * -----------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------------
 * Given a matrix a[0..n-1][0..n-1], this routine replace it by the LU decompositon
 * of a rowwise permutation of itself. a is input. On output, it is arranged as in
 * equation (2.3.14) (see <Numerical Recipes in C++> Second Edition, authored by 
 * Willam H. Press, etc.); indx[0..n-1] is an output vector that records the row 
 * permutation effected by the partial pivoting; d is output as +/-1 depending on 
 * whether the number of row interchanges was even or odd, respectively. This routine
 * is used in combination with lubksb to solv linear equations or invert a matrix.
 * -----------------------------------------------------------------------------*/
int ludcmp(int n, double a[], int indx[], double *d);

/* -------------------------------------------------------------------------------
 * Solves the set of n linear equations AX = B. Here a[0..n-1][0..n-1] is input, not 
 * as the matrix A but rather as its LU decompositon, determined by the outine ludcmp.
 * index[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is
 * input as the right-hand side vector B, and returns with the solution vector X. a and
 * indx are not modified by this routine and can be left in place for successive calls
 * with different right-hand sides b. This routine takes into account the possibility 
 * that b will begin with many zero elements, so its efficient for use in matrix inversion
 * -----------------------------------------------------------------------------*/
int lubksb(int n, double a[], int indx[], double b[]);

/* -------------------------------------------------------------------------------
 * 方程根求解函数
 * -----------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------------
 * Given a intial guess x[0..n-1], for a root in n dimensions, take ntrial Newton
 * -Raphson steps to improve the root. Stop if the root converges in either summed 
 * abosolute variable increments tolx or summed absolute function values tolf
 * -----------------------------------------------------------------------------*/
int mnewt(int ntrial, int n, double x[], 
		  int func(double x[], int n, double fv[]), 
	  	  int jac_func(double x[], int n, double jfv[]),
		  double tolx, double tolf);

/* -------------------------------------------------------------------------------
 * Given a function func, and an initial guessed range x1 to x2, the routine expands
 * the range geometrically until a root is bracket by the returned values x1 and x2 
 * (in which case zbrac returns true) or until the range becomes unacceptably large
 * (in which case zbrac returns false).
 * -----------------------------------------------------------------------------*/
bool zbrac(double func(const double), double &x1, double &x2);

/* -------------------------------------------------------------------------------
 * Given a function func defined on the interval from x1 to x2 subdivided into n 
 * equally spaced segments, and search for zero crossings of the function. The array
 * xb1[0 ... nb - 1] and xb2[0 ... nb - 1] will be filled sequentially with any 
 * bracketing pairs that are found, and must be provided with a size nb that is 
 * sufficient to hold the maximum number of roots sought. nroot will be set to the 
 * number of bracketing pairs actually found.
 * -----------------------------------------------------------------------------*/
void zbrak(double func(const double), const double x1, const double x2, const int n,
		   double *xb1, double *xb2, int nb, int &nroot);

/* -------------------------------------------------------------------------------
 * Using bisection, find the root of a function func known to lie between x1 and x2.
 * The root, returned as rtbis, will be refined until its accuracy is +/-xacc
 * -----------------------------------------------------------------------------*/
double rtbis(double func(const double), const double x1, const double x2, 
			 const double xacc);


#if 0
/* -------------------------------------------------------------------------------
 * Given a intial guess x[0..n-1] for a root in n dimensions, take ntrial Newton-Raphson
 * steps to improve the root. Stop if the root converge in either summed absolute 
 * variable increments tolx or summed absolute function values tolf
 * -----------------------------------------------------------------------------*/
void mnewt(int (*func)(int, double*, double*),
		   int (*jaco_func)(int, double*, double*),
		   int n, double x[],
		   int ntrial,  double tolx,
		   double tolf);
#endif

/* -------------------------------------------------------------------------------
 * Given a intial guess (x[0] x[1]) for a root in 2 dimensions, take ntrial 
 * Newton-Raphson steps to improve the root. Stop if the root converge in either 
 * summed absolute variable increments tolx or summed absolute function values tolf
 * -----------------------------------------------------------------------------*/
int iterative_newton_raphson_2D(
		int (*func)(double *, double*),
		int (*jaco_func)(double*,double*),
		double x[], int ntrial, double tolx, double tolf);

/* -------------------------------------------------------------------------------
 * 一维搜索函数
 * -----------------------------------------------------------------------------*/

/* -------------------------------------------------------------------------------
 * 辅助函数
 * -----------------------------------------------------------------------------*/
inline double DSIGN(double x1, double x2)
{
	return x2 == 0.0 ? 0.0 : (x2 > 0.0 ? fabs(x1) : -fabs(x1)); 
}

inline double DMAX(double x1, double x2)
{
	return x1 > x2 ? x1 : x2;
}

inline void shft2(double &a, double &b, const double &c)
{
	a = b;
	b = c;
}

inline void shft3(double &a, double &b, double &c, const double &d)
{
	a = b;
	b = c;
	c = d;
}

/* -------------------------------------------------------------------------------
 * Given a function func, and given distinct initial points ax and bx, this routine 
 * searches in the downhill direction (defined by the function as evaluated at the
 * initial points) and returns new point ax, bx, cx that bracket a minimum of the
 * function. Also returned are the function values at the three points, fa, fb, and fc
 * -----------------------------------------------------------------------------*/
void mnbrak(double &ax, double &bx, double &cx, 
			double &fa, double &fb, double &fc,
			double func(const double));

/* -------------------------------------------------------------------------------
 * Given a function func f, and given a bracketing triplet of abscissas ax, bx, cx
 * such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx),
 * this routine performs a golden search for the minimum, isolating it to a fractional 
 * precision about tol. The abscissa of the minimum is returned as xmin, and the 
 * minimum function value is returned as the returned function value 
 * -----------------------------------------------------------------------------*/
double golden_search(const double ax, const double bx, const double cx, 
					 double func(const double), const double tol,
					 double &xmin);


/* -------------------------------------------------------------------------------
 * Given a function func f, and given a bracketing triplet of abscissas ax, bx, cx
 * such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx),
 * this routine isolates the minimum to a fractional precision of about tol using
 * Brent's method. The abscissa of the minimum is returned  as xmin, and the minimum 
 * function value is returned as the returned function value 
 * -----------------------------------------------------------------------------*/
double brent_search(const double ax, const double bx, const double cx, 
				    double func(const double), const double tol,
					double &xmin);

/* --------------------------------------------------------------------------------
 * 一维迭代复合Simpson积分
 * 输入参数:
 *     a, b               积分域
 *     (*func)(double)    积分函数
 *     tol                容许的积分误差
 *     max_iter           最大迭代次数
 * 返回值:
 *	   函数(*func)(double)在[a b]上的积分值
 * 描述:
 *
 *     当迭代次数超过max_iter或前后两次积分值的差值和前次积分值之比小于tol时，
 *     积分结束。
 * -------------------------------------------------------------------------------*/
double iterative_composite_simpson_1D(double a, double b,  /* 积分域 */
									  double (*func)(double), /* 一维积分函数 */
									  double tol = 1.0e-6, /* 容许的积分误差 */
									  int max_iter = 30 /* 最大迭代次数 */
									  );

/* --------------------------------------------------------------------------------
 * 一维迭代Newton-Raphson方程根求解算法
 * 输入参数:
 *     x_init					初值
 *     (*func)(double)			原函数
 *	   (*deriv_func)(double)    导函数 
 *     tol						容许的积分误差
 *     max_iter					最大迭代次数
 * 返回值:
 *	   (*func)(double) = 0的根
 * 描述:
 *
 *     当迭代次数超过max_iter或前后两次方程根的差值小于tol时，计算结束
 * -------------------------------------------------------------------------------*/
double iterative_newton_raphson_1D(
						double x_init,							/* 初值 */
						double (*func)(double),					/* 原函数 */
						double (*deriv_func)(double),			/* 导函数 */
						double tol = 1.0e-4,                    /* 容许误差 */
						int max_iter = 20					/* 最大迭代次数 */
						);

}

#endif /* __eemas_math_util_h__ */