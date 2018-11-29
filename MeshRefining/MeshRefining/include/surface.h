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
 * Copyright reserved, May. 24, 2008
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
 *     Define a set of interface for surfaces
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/


#ifndef __eemas_surface_h__
#define __eemas_surface_h__

namespace EEMAS {

class Surface
{
public:
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
				double tol = 1.0e-6) = 0;

	/* --------------------------------------------------------------------
	 * 通过点坐标获取相应的参数值
	 * 输入参数：
	 *   u, v		参数坐标
	 *   coord		点的坐标(返回值)
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	virtual int	param_to_coord(double u, double v, Vector *coord) = 0;

	/* --------------------------------------------------------------------
	 * 获取曲面上离某点最近的点的坐标
	 * 输入参数：
	 *   coord		点的坐标
	 *   u, v		最近点(垂足)的参数坐标(返回值)
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	virtual int project(Vector coord, double *u, double *v, double tol = 1.0e-6) = 0;

	/* --------------------------------------------------------------------
	 * 获取一阶偏导数
	 * 输入参数：
	 *   u, v		参数坐标
	 *   d1u, d1v	u向偏导/v向偏导
	 * 返回值：
	 *   错误代码
	 * 描述：
	 * -------------------------------------------------------------------*/
	virtual int	param_to_uv_d1(double u, double v, Vector *d1u, Vector *d1v) = 0;

protected:
};

}/* namespace EEMAS */

#endif /* __eemas_surface_h__ */
