#ifndef __dihe_ang_h__
#define __dihe_ang_h__

#include <stdio.h>
#include <stdlib.h>

//class DTIso3D; /* �����iso3d.h�������ã���������Ķ��� */

/* --------------------------------------------------------------------------------------
 * Ϊ�˺ͱ�����DTIso3D�е��������ݽ����νӣ���Ҫ����һ��ȫ��ָ��
 * ----------------------------------------------------------------------------------*/
//extern DTIso3D *g_DTIso3D = NULL;

/* ---------------------------------------------------------------------------------------
 * ���ȣ�Ҫ�����е�Ԫ�б������ж���ǣ��Ǵ洢˳��ͱ���һ�£����ڱ任ʱע�⼰ʱ����
 * -------------------------------------------------------------------------------------*/
typedef struct
{
	float sinValue;
	int heapIndex;	//����λ�洢���ȼ������ȼ�0��ߣ�6���
}ElemAngle;

typedef struct
{
	float sinValue;

	//���һλ�洢�Ƿ�Ϊ��ǣ��ӵ�������λ��ʼ����λ�洢���ȼ������ȼ�0��ߣ�6��ͣ�ǰ��ʣ��λ�洢��heap�е�index
	int heapIndex;	
}ElemQual;

typedef ElemQual TetraElemQual;
//typedef ElemAngle TetraAngles[6];
//typedef float TetraAngleSins[6];

#define Type TetraElemQual
#define BlockedArray_Type BlockedArray_TetQuals
#include "blockedarray.h"

//extern BlockedArray_TetraAngles g_TetraAngles;	//�漰ɾ������ӣ����Ƿ���п�λ����

/* ---------------------------------------------------------------------------------------
 * ÿ��������Ķ�������¼��Ԫ�š��Ƕ�Ӧ�ı��룬��ˣ�ֻ��һ���������ɣ������ĺ�3λ��¼
 * ���룬ǰ������λ��¼��Ԫ��
 * 1. ����Ӧ�ĵ�Ԫ��ɾ�����ö������ȻҲ��ɾ��
 * --------------------------------------------------------------------------------------*/
typedef struct
{
	int Elem;
	ElemQual *elmQual;
}MeshElment;

typedef struct
{
	int angleCode;
	ElemAngle *angleVal;
}DiheAngle;

/* ---------------------------------------------------------------------------------------
 * ���ڵ������ǣ���Ԫ���ܱ�ɾ�����ֱ�����ռλ����ˣ��ڵ�Ԫ��ɾ����ʱ�����Ӧ�Ķ����
 * Ҳ��Ҫ��ɾ�����ʶ�����Ҫ��ͨ����Ԫ������ѯ����Ӧ������ڶ��������
 * 
 * ��ÿ����Ԫ��6���Ƕ���Ҫ�����������������ڴ���Ϊ6 x ��Ԫ�� x sizeof(int)�ֽ�
 * 32λϵͳ�£�1����Ԫ���ڴ�ԼΪ24MB
 *
 * ��Ӧ�Ż����⣬����ֻ��10%�ĽǶ�С��30�Ȼ����150�ȣ���Ϊ��Ծ����ǣ��������Ҫ�Ż���
 * ������Ҫ���Ǹ���ʡ�ڴ�����ݽṹ
 * ActAnglePos���飺��¼ÿ����Ԫ�ж��ٸ��������Ҫ����(��3λ)�����ڴ洢��ActAngleIDs�е�
 *                  ��ʼλ�ã�ǰ������λ��
 * ActAngleIDs���飺�洢��Ԫ������ڶ��������
 *
 * ����10%�Ľ�Ϊ��Ծ�����Ϊ�����������ڴ���Ϊ 10% x 6 x ��Ԫ�� x sizeof(int)�ֽ� + 
 * ��Ԫ�� x sizeof(int); 32λϵͳ�£�1����Ԫ���ڴ�ԼΪ6.4MB
 * --------------------------------------------------------------------------------------*/
//extern BlockedArray_Int g_ActAnglePos;  
//extern BlockedArray_Int g_ActAngleIDs;	//ÿ�ζ��ں�����ӣ�ǰ����Чλ����δ���

/* ---------------------------------------------------------------------------------------
 * ����g_DTIso3D�����������Ի�ȡ�����Ķ����SINֵ
 * -------------------------------------------------------------------------------------*/
float getAngleSin(DiheAngle ang);
/* ---------------------------------------------------------------------------------------
 * ����getAngleSin���Ƚ�ang1��ang2
 * -------------------------------------------------------------------------------------*/
int compareDiheAng(DiheAngle ang1, DiheAngle ang2);

int setIndex(int index, int priority);
int getIdx(int index);
int getPriority(int index);
int getElem(DiheAngle ang);
int getAngleCode(DiheAngle ang);
int getEdge(DiheAngle ang, int ep[2]);
int getEdgeFace(DiheAngle ang, int fc[2]);

/* ------------------------------------------------------------------------------------
 * ���е�ԪӦ�ð����ȼ����гɶ�
 * ����Ҫ�߱�������
 * 1. ����Ԫ�ذ��Ƕȴ�С��С��������
 * 2. ������LOG(N)��ʱ��������ɾ��һ��Ԫ��
 * -----------------------------------------------------------------------------------*/
#define Type MeshElment
#define BlockedArray_Type BlockedArray_ActElem
#include "blockedarray.h"
#define LOG2_SKIRTPOLY_PER_BLOCK 12

class ActElemHeap
{
public:
	ActElemHeap():nHeapSize(0), nMaxHeapSize(LOG2_SKIRTPOLY_PER_BLOCK)
	{ 
		nChanged = 0;
		changedAngle = NULL;
		m_actElmHeap.init(nMaxHeapSize);
		m_actElmHeap.allocMemory(nMaxHeapSize); 
	}
	~ActElemHeap(){};

public:
	MeshElment getMinElem();
	void rmMinElem();
	int insertElem(MeshElment ang);
	void removeElem(int idx);
	void minHeapifyDown(int beg);
	void minHeapifyUp(int beg);
	void addPriority(int idx); //����idx����Ԫ�����ȼ���ֵ

	bool isEmpty(){ return nHeapSize == 0; }
	bool isFull(){ return nHeapSize == nMaxHeapSize; };
	int heapSize(){ return nHeapSize; }
	void allocHeap()
	{
		nMaxHeapSize = 1.25*nMaxHeapSize;
		m_actElmHeap.allocMemory(nMaxHeapSize);
	}

	void printfinfo();

private:
	int nHeapSize;
	int nMaxHeapSize;
	int nChanged;
	int *changedAngle;
	BlockedArray_ActElem m_actElmHeap;
};


#define Type DiheAngle
#define BlockedArray_Type BlockedArray_ActAngle
#include "blockedarray.h"
#define LOG2_SKIRTPOLY_PER_BLOCK 12

class ActAngleHeap
{
public:
	ActAngleHeap():nHeapSize(0), nMaxHeapSize(LOG2_SKIRTPOLY_PER_BLOCK)
	{ 
		nChanged = 0;
		changedAngle = NULL;
		m_actAngleHeap.init(nMaxHeapSize);
		m_actAngleHeap.allocMemory(nMaxHeapSize); 
	}
	~ActAngleHeap(){};

public:
	DiheAngle getMinAngle();
	void rmMinAngle();
	int insertAngle(DiheAngle ang);
	void removeAngle(int idx);		//���Ϊ��Ч
	//DiheAngle rootAngle();
	void minHeapifyDown(int beg);
	void minHeapifyUp(int beg);
	void addPriority(int idx); //����idx����Ԫ�����ȼ���ֵ

	bool isEmpty(){ return nHeapSize == 0; }
	bool isFull(){ return nHeapSize == nMaxHeapSize; };
	int heapSize(){ return nHeapSize; }
	void allocHeap()
	{
		nMaxHeapSize = 1.25*nMaxHeapSize;
		m_actAngleHeap.allocMemory(nMaxHeapSize);
	}

	void printfinfo();

private:
	int nHeapSize;
	int nMaxHeapSize;
	int nChanged;
	int *changedAngle;
	BlockedArray_ActAngle m_actAngleHeap;
};
#endif 