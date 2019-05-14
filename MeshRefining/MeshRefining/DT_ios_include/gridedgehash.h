/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2009年10月30日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 30, 2009
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953166
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：gridedgehash.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义边表类GridEdgeHash (哈希表)，用于记录顶点之间的邻接关系
 *     Define class GridEdgeHash, to record the neighboring info. between grid nodes
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#ifndef __eemas_gridedgehash_h__
#define __eemas_gridedgehash_h__

#include <vector>
#include <map>

/* 哈希表中每个节点的类型 */
typedef struct GEdgeHashNode 
{
	int node2;						/* node1 为关键字，不重复记录 */
	std::vector<int>	vecElems;	/* 包含这条表边的单元列表 optional, 当m_bRecordElems为true是有效 */
} GEdgeHashNode;

typedef std::multimap<int, GEdgeHashNode*> HashNodeMultiMap;

/* -----------------------------------------------------------------------------
 * 定义边表类GridEdgeHash (哈希表)，用于记录顶点之间的邻接关系
 * ----------------------------------------------------------------------------*/
class GridEdgeHash: public HashNodeMultiMap
{
public:
	
	/* ---------------------------------------------------------------------------
	 * 私有数据类型
	 * ---------------------------------------------------------------------------*/
	/* 利用边的两个节点索引的小值、大值，或两个节点索引作为哈希表的关键字 */
	enum HashType {MIN_KEY, MAX_KEY, DOUB_KEY};	


	/* ---------------------------------------------------------------------------
	 * 构造和析构函数
	 * ---------------------------------------------------------------------------*/
	GridEdgeHash()	{
		m_eHashType = MIN_KEY;
		m_bElemRecorded = false;
	}

	virtual ~GridEdgeHash();

	/* ---------------------------------------------------------------------------
	 * 属性设置和查询函数
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
	 * 根据输入网格信息初始化哈希表
	 * ---------------------------------------------------------------------------*/
	int initialize(
		std::vector<int> vecETopus,	/* 单元拓扑 */
		std::vector<int> vecENodes	/* 单元索引 */
		);

	/* ---------------------------------------------------------------------------
	 * 根据哈希表查询节点的所有邻接节点
	 * ---------------------------------------------------------------------------*/
	int find_neighbors(int nodeKey, std::vector<int>& vecNeigNodes);
	
	/* ---------------------------------------------------------------------------
	 * 根据哈希表查询包含节点的所有单元
	 * ---------------------------------------------------------------------------*/
	int find_node_elems(int nodeKey, std::vector<int>& vecNeigElems);

	/* ---------------------------------------------------------------------------
	 * 根据哈希表查询包含边的所有单元
	 * ---------------------------------------------------------------------------*/
	int find_edge_elems(int nodeKey, int node2, std::vector<int>& vecNeigElems);

protected:
	/* ---------------------------------------------------------------------------
	 * 受保护函数
	 * ---------------------------------------------------------------------------*/
	void release();	/* 释放资源 */
	
	int add_edge(
		int nodeKey,		/* 关键点 */
		int node2,			/* 第二节点 */
		int elem,			/* 单元 */
		int *addCount		/* 增加的类型：
							    0. 第一次增加该关键字，
							    1. 第一次增加该边
							   >1. 再次增加该边 */
		);

protected:
	
	HashType m_eHashType;	/* 哈希表的类型, 缺省为MIN_KEY */
	bool m_bElemRecorded;	/* 是否记录包含单元的边的信息, 缺省为false */
};

#endif /* __eemas_gridedgehash_h__ */