/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2008��06��11��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, June. 11, 2007
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953166
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��geometryquery.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ����Ա߽�����������㷨
 *     Define a set of geometry query functions
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
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
extern int query_curve_connected(GBCurve *curve,vector<GBCurve*> &curve_connected_start,vector<GBCurve*> &curve_connected_end);//��ѯ��curve������curve,��ɾ���ظ����ſ��ô˺���
extern int query_loop(GBCurve *curve,vector<GBLoop*> &loop_vec);	//��ѯcurve���ڵĻ�
};

#endif /* __eemas_geometry_query_h__ */