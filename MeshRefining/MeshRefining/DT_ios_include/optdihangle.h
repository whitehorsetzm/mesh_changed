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
	//���һλ�洢�Ƿ�Ϊ��ǣ��ӵ�������λ��ʼ����λ�洢���ȼ���
	//���ȼ�0��ߣ�6��ͣ��ӵ�������λ��ʼ����λ�洢�ǵı�ţ�
	//ǰ��ʣ��λ�洢��heap�е�index
	int iReserved;	
}ElemQual;

typedef ElemQual TetraElemQual;

#define Type TetraElemQual
#define BlockedArray_Type BlockedArray_TetQuals
#include "blockedarray.h"


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
float getAngleSin(MeshElment me);
/* ---------------------------------------------------------------------------------------
 * ����getAngleSin���Ƚ�ang1��ang2
 * -------------------------------------------------------------------------------------*/
int compareDiheAng(MeshElment me1, MeshElment me2);

//setter function
void setAcute(int &val, int acute);
void setPriority(int &val, int priority);
void setAngleCode(int &val, int angCode);
void setHeapIndex(int &val, int idx);

//getter function
int getAcute(int val);
int getPriority(int val);
int getAngleCode(int val);
int getHeapIndex(int val);

int getEdge(MeshElment ang, int ep[2]);
int getEdgeFace(MeshElment ang, int fc[2]);

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
	int insertElem(MeshElment elm);
	void removeElem(int idx);
	void minHeapifyDown(int beg);
	void minHeapifyUp(int beg);
	void addPriority(int idx); //����idx����Ԫ�����ȼ���ֵ

	bool isEmpty(){ return nHeapSize == 0; }
	bool isFull(){ return nHeapSize == nMaxHeapSize; };
	int heapSize();
	void allocHeap()
	{
		nMaxHeapSize = (int)(1.25*nMaxHeapSize);
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

#endif 