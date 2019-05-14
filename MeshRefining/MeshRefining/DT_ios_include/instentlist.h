#ifndef __inst_ent_list_h
#define __inst_ent_list_h

#ifdef _USING_STD_LIB
#include <vector>
#else
#include "myvector.h"
#endif

#define MAX_INST_ENT_SIZE 1024*10//1024

#define IET_INACTIVE_FLAG			1
#define IET_INVALID_FLAG			(1<<1)
#define IET_IRRECOVERABLE_FLAG		(1<<2)
#define IET_NEWBORN_FLAG			(1<<3)
#define IET_NOEDIT_FLAG				(1<<4)

enum InstEntType {IET_NUL=0, IET_NOD, IET_EDG, IET_FAC};

typedef struct InstEntity
{
	InstEntType type;
	int i1, i2, i3; 
	int prt1, prt2;	/* 面片的2个父亲单元 */
	int idx;
	int flag;
	/* --------------------------------------------------------
	 * 有时需要采用counter而非flag来标志一个实体的引用情况，例如：
	 * counter = 0：表示资源不被引用；
	 * counter > 0：表示资源被占用
	 * -------------------------------------------------------*/
	int counter;	
}InstEntity;

typedef struct InstEntityLnk
{
	InstEntity instEnt;
	int next;
} InstEntityLnk;

/* 将其组织成哈希表 */
class InstEntList
{
public:
	/* 构造、析构函数 */
	InstEntList();
	~InstEntList();

	/* 初始化 */
	int initialise(int maxIdx);
	int resetDirtNodes();

	/* 查询、增加和减少 */
	int findEdge(int i1, int i2);
	int findFace(int i1, int i2, int i3);

	int insertEdge(int i1, int i2);
	int insertFace(int i1, int i2, int i3);

	int removeEdge(int i1, int i2);
	int removeFace(int i1, int i2, int i3);
	int removeEnt(int iEnt);

	/* 遍历所有InstEntities */
	int getFirstEnt();
	int getNextEnt(int iEnt);
	bool isLocked();
	void lock();
	void unlock();

	bool isEmpty();
	int getSize();

	/* 状态标志 */
	void setRecoverable(int iEnt, bool flag);
	void setValid(int iEnt, bool flag);
	void setActive(int iEnt, bool flag);
	void setNewBorn(int iEnt, bool flag);
	void setNoEdit(int iEnt, bool flag);
	bool isRecoverable(int iEnt);
	bool isValid(int iEnt);
	bool isActive(int iEnt);
	bool isNewBorn(int iEnt);
	bool isNoEdit(int iEnt);
	int counterIncre(int iEnt);
	int counterDecre(int iEnt);
	int entCounter(int iEnt);
	void setEntCounter(int iEnt, int entCnt);
	InstEntType entType(int iEnt);
	int nodeIndices(int iEnt, int indices[]);

protected:
	int *m_pNodeHash;  /* 哈希表 */
	int m_nHashMaxIdx; /* 最大索引值 */
	InstEntityLnk m_arrInstEntLnks[MAX_INST_ENT_SIZE]; /* 所有InstEntities */
	int m_nInstEntArrSize; /* 数组的有效大小 */
	int m_nValidInstEnts;  /* 真正有效的相交几何元素的个数 */		
	/* 在遍历数组时，可能还会增加元素，为保证遍历是对原数组的遍历，可对其加锁 */
	bool m_bLocked;
	int m_nLockedEntArrSize;

#ifdef _USING_STD_LIB
	/*a vector of tested elements*/
	std::vector<INTEGER> m_vecDirtNodes;
#else
	MyVector<int> m_vecDirtNodes; /* 在哈希表中拥有关键字的节点数组 */
#endif
};

typedef InstEntList InflatationEntList;

#endif