/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2008��05��19��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, May. 18, 2007
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953166
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��curve.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ��������ߵ�ͨ�ýӿ�
 *     Define a set of interface for curve
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
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
	 * ��ȡ�˵�����
	 * -----------------------------------------------------------------*/
	virtual Vector start_coord() = 0;
	virtual Vector end_coord() = 0; 

	/* --------------------------------------------------------------------
	 * ��ȡ�ܻ���
	 * -------------------------------------------------------------------*/
	virtual double get_total_arc_length() = 0;

	/* --------------------------------------------------------------------
	 * ͨ���������������һ������
	 * -------------------------------------------------------------------*/
	virtual int	 arc_length_to_coord(
								double al,              /* ���� */
								Vector *coord,		    /* ������� */
								double tol = 1.0e-4,    /* ������� */
								int max_iter = 20	    /* ���������� */
								) = 0;
protected:
};

}/* namespace EEMAS */

#endif /* __eemas_curve_h__ */
