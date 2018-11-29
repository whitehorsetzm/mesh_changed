/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2008年06月11日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, June. 11, 2007
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953166
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：geometryquery.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对边界网格的排序算法
 *     Define a set of geometry query functions
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#ifndef __eemas_geometry_query_h__
#define __eemas_geometry_query_h__

#include <vector>

using namespace std;

namespace EEMAS
{

/* ---------------------------------------------------------
 * Support surfaces 
 * ---------------------------------------------------------*/
extern int query_data_points(GBFace *face, int *u, int *v, double **uv);
extern double calculate_curve_length(GBCurve *curve);
extern int query_curve_connected(GBCurve *curve,vector<GBCurve*> &curve_connected_start,vector<GBCurve*> &curve_connected_end);//查询与curve相连的curve,在删除重复点后才可用此函数
extern int query_loop(GBCurve *curve,vector<GBLoop*> &loop_vec);	//查询curve所在的环
};

#endif /* __eemas_geometry_query_h__ */