/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2009��10��30��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 30, 2009
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953166
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��gridedgehash.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ����߱���GridEdgeHash (��ϣ��)�����ڼ�¼����֮����ڽӹ�ϵ
 *     Define class GridEdgeHash, to record the neighboring info. between grid nodes
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#ifndef __eemas_gridedgehash_h__
#define __eemas_gridedgehash_h__

#include <vector>
#include <map>

/* ��ϣ����ÿ���ڵ������ */
typedef struct GEdgeHashNode 
{
	int node2;						/* node1 Ϊ�ؼ��֣����ظ���¼ */
	std::vector<int>	vecElems;	/* ����������ߵĵ�Ԫ�б� optional, ��m_bRecordElemsΪtrue����Ч */
} GEdgeHashNode;

typedef std::multimap<int, GEdgeHashNode*> HashNodeMultiMap;

/* -----------------------------------------------------------------------------
 * ����߱���GridEdgeHash (��ϣ��)�����ڼ�¼����֮����ڽӹ�ϵ
 * ----------------------------------------------------------------------------*/
class GridEdgeHash: public HashNodeMultiMap
{
public:
	
	/* ---------------------------------------------------------------------------
	 * ˽����������
	 * ---------------------------------------------------------------------------*/
	/* ���ñߵ������ڵ�������Сֵ����ֵ���������ڵ�������Ϊ��ϣ��Ĺؼ��� */
	enum HashType {MIN_KEY, MAX_KEY, DOUB_KEY};	


	/* ---------------------------------------------------------------------------
	 * �������������
	 * ---------------------------------------------------------------------------*/
	GridEdgeHash()	{
		m_eHashType = MIN_KEY;
		m_bElemRecorded = false;
	}

	virtual ~GridEdgeHash();

	/* ---------------------------------------------------------------------------
	 * �������úͲ�ѯ����
	 * ---------------------------------------------------------------------------*/
	void set_hash_type(HashType eHashType) {
		m_eHashType = eHashType;
	}
	HashType get_hash_type() {
		return m_eHashType;
	}
	void enable_elem_record(bool bEnable) {
		m_bElemRecorded = bEnable;
	}
	bool is_elem_recorded() {
		return m_bElemRecorded;
	}


	/* ---------------------------------------------------------------------------
	 * ��������������Ϣ��ʼ����ϣ��
	 * ---------------------------------------------------------------------------*/
	int initialize(
		std::vector<int> vecETopus,	/* ��Ԫ���� */
		std::vector<int> vecENodes	/* ��Ԫ���� */
		);

	/* ---------------------------------------------------------------------------
	 * ���ݹ�ϣ���ѯ�ڵ�������ڽӽڵ�
	 * ---------------------------------------------------------------------------*/
	int find_neighbors(int nodeKey, std::vector<int>& vecNeigNodes);
	
	/* ---------------------------------------------------------------------------
	 * ���ݹ�ϣ���ѯ�����ڵ�����е�Ԫ
	 * ---------------------------------------------------------------------------*/
	int find_node_elems(int nodeKey, std::vector<int>& vecNeigElems);

	/* ---------------------------------------------------------------------------
	 * ���ݹ�ϣ���ѯ�����ߵ����е�Ԫ
	 * ---------------------------------------------------------------------------*/
	int find_edge_elems(int nodeKey, int node2, std::vector<int>& vecNeigElems);

protected:
	/* ---------------------------------------------------------------------------
	 * �ܱ�������
	 * ---------------------------------------------------------------------------*/
	void release();	/* �ͷ���Դ */
	
	int add_edge(
		int nodeKey,		/* �ؼ��� */
		int node2,			/* �ڶ��ڵ� */
		int elem,			/* ��Ԫ */
		int *addCount		/* ���ӵ����ͣ�
							    0. ��һ�����Ӹùؼ��֣�
							    1. ��һ�����Ӹñ�
							   >1. �ٴ����Ӹñ� */
		);

protected:
	
	HashType m_eHashType;	/* ��ϣ�������, ȱʡΪMIN_KEY */
	bool m_bElemRecorded;	/* �Ƿ��¼������Ԫ�ıߵ���Ϣ, ȱʡΪfalse */
};

#endif /* __eemas_gridedgehash_h__ */