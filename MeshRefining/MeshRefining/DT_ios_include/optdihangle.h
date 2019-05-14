#ifndef __dihe_ang_h__
#define __dihe_ang_h__

#include <stdio.h>
#include <stdlib.h>

//class DTIso3D; /* 避免和iso3d.h交叉引用，故申明类的定义 */

/* --------------------------------------------------------------------------------------
 * 为了和保存在DTIso3D中的网格数据进行衔接，需要定义一个全局指针
 * ----------------------------------------------------------------------------------*/
//extern DTIso3D *g_DTIso3D = NULL;

/* ---------------------------------------------------------------------------------------
 * 首先，要在所有单元中保存所有二面角，角存储顺序和边码一致，并在变换时注意及时更新
 * -------------------------------------------------------------------------------------*/
typedef struct
{
	float sinValue;
	//最后一位存储是否为锐角，从倒数第四位开始的三位存储优先级，
	//优先级0最高，6最低，从倒数第七位开始的三位存储角的编号，
	//前面剩余位存储在heap中的index
	int iReserved;	
}ElemQual;

typedef ElemQual TetraElemQual;

#define Type TetraElemQual
#define BlockedArray_Type BlockedArray_TetQuals
#include "blockedarray.h"


/* ---------------------------------------------------------------------------------------
 * 每个待处理的二面角需记录单元号、角对应的编码，因此，只需一个整数即可，整数的后3位记录
 * 编码，前面所有位记录单元号
 * 1. 若对应的单元被删除，该二面角自然也被删除
 * --------------------------------------------------------------------------------------*/
typedef struct
{
	int Elem;
	ElemQual *elmQual;
}MeshElment;

/* ---------------------------------------------------------------------------------------
 * 现在的问题是：单元可能被删除，又被重新占位，因此，在单元被删除的时候，其对应的二面角
 * 也需要被删除，故而，需要能通过单元索引查询到对应二面角在堆里的索引
 * 
 * 如每个单元的6个角都需要保留索引，则所需内存数为6 x 单元数 x sizeof(int)字节
 * 32位系统下，1百万单元需内存约为24MB
 *
 * 对应优化问题，假设只有10%的角度小于30度或大于150度（称为活跃二面角），因而需要优化，
 * 我们需要考虑更节省内存的数据结构
 * ActAnglePos数组：记录每个单元有多少个二面角需要保留(后3位)及其在存储在ActAngleIDs中的
 *                  起始位置（前面其它位）
 * ActAngleIDs数组：存储单元二面角在堆里的索引
 *
 * 仍以10%的角为活跃二面角为例，则所需内存数为 10% x 6 x 单元数 x sizeof(int)字节 + 
 * 单元数 x sizeof(int); 32位系统下，1百万单元需内存约为6.4MB
 * --------------------------------------------------------------------------------------*/
//extern BlockedArray_Int g_ActAnglePos;  
//extern BlockedArray_Int g_ActAngleIDs;	//每次都在后面添加，前面无效位置如何处理

/* ---------------------------------------------------------------------------------------
 * 访问g_DTIso3D的网格数据以获取真正的二面角SIN值
 * -------------------------------------------------------------------------------------*/
float getAngleSin(MeshElment me);
/* ---------------------------------------------------------------------------------------
 * 调用getAngleSin，比较ang1和ang2
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
 * 所有单元应该按优先级排列成堆
 * 堆需要具备的特征
 * 1. 所有元素按角度大小从小到大排列
 * 2. 可以在LOG(N)的时间里插入和删除一个元素
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
	void addPriority(int idx); //增加idx个单元的优先级数值

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