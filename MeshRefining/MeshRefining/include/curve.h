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
 * 文件名称 (File Name)：curve.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对曲线的通用接口
 *     Define a set of interface for curve
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/


#ifndef __eemas_curve_h__
#define __eemas_curve_h__

#include "vector.h"
namespace EEMAS {

class Curve
{
public:
	/* -------------------------------------------------------------------
	 * 获取端点坐标
	 * -----------------------------------------------------------------*/
	virtual Vector start_coord() = 0;
	virtual Vector end_coord() = 0; 

	/* --------------------------------------------------------------------
	 * 获取总弧长
	 * -------------------------------------------------------------------*/
	virtual double get_total_arc_length() = 0;

	/* --------------------------------------------------------------------
	 * 通过弧长求解曲线上一点坐标
	 * -------------------------------------------------------------------*/
	virtual int	 arc_length_to_coord(
								double al,              /* 弧长 */
								Vector *coord,		    /* 点的坐标 */
								double tol = 1.0e-4,    /* 容许误差 */
								int max_iter = 20	    /* 最大迭代次数 */
								) = 0;
protected:
};

}/* namespace EEMAS */

#endif /* __eemas_curve_h__ */
