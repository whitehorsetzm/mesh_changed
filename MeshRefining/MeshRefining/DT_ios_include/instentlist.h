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
	int prt1, prt2;	/* ��Ƭ��2�����׵�Ԫ */
	int idx;
	int flag;
	/* --------------------------------------------------------
	 * ��ʱ��Ҫ����counter����flag����־һ��ʵ���������������磺
	 * counter = 0����ʾ��Դ�������ã�
	 * counter > 0����ʾ��Դ��ռ��
	 * -------------------------------------------------------*/
	int counter;	
}InstEntity;

typedef struct InstEntityLnk
{
	InstEntity instEnt;
	int next;
} InstEntityLnk;

/* ������֯�ɹ�ϣ�� */
class InstEntList
{
public:
	/* ���졢�������� */
	InstEntList();
	~InstEntList();

	/* ��ʼ�� */
	int initialise(int maxIdx);
	int resetDirtNodes();

	/* ��ѯ�����Ӻͼ��� */
	int findEdge(int i1, int i2);
	int findFace(int i1, int i2, int i3);

	int insertEdge(int i1, int i2);
	int insertFace(int i1, int i2, int i3);

	int removeEdge(int i1, int i2);
	int removeFace(int i1, int i2, int i3);
	int removeEnt(int iEnt);

	/* ��������InstEntities */
	int getFirstEnt();
	int getNextEnt(int iEnt);
	bool isLocked();
	void lock();
	void unlock();

	bool isEmpty();
	int getSize();

	/* ״̬��־ */
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
	int *m_pNodeHash;  /* ��ϣ�� */
	int m_nHashMaxIdx; /* �������ֵ */
	InstEntityLnk m_arrInstEntLnks[MAX_INST_ENT_SIZE]; /* ����InstEntities */
	int m_nInstEntArrSize; /* �������Ч��С */
	int m_nValidInstEnts;  /* ������Ч���ཻ����Ԫ�صĸ��� */		
	/* �ڱ�������ʱ�����ܻ�������Ԫ�أ�Ϊ��֤�����Ƕ�ԭ����ı������ɶ������ */
	bool m_bLocked;
	int m_nLockedEntArrSize;

#ifdef _USING_STD_LIB
	/*a vector of tested elements*/
	std::vector<INTEGER> m_vecDirtNodes;
#else
	MyVector<int> m_vecDirtNodes; /* �ڹ�ϣ����ӵ�йؼ��ֵĽڵ����� */
#endif
};

typedef InstEntList InflatationEntList;

#endif