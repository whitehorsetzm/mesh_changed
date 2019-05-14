#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#define MAX_SKIRT_POLY 128 
#define _OUTPUT_LEVEL_III

enum QUALTYPE
{
	DIHEDRAL_ANGLE,
	DIHEDRAL_ANGLE_SINE,
	SOLID_ANGLE,
	ASPECT_RATIO
};

class HashPolyNode;
class HashPolyList;

//链表节点，记录每条边的壳信息
class SkirtPolyNode
{
	friend HashPolyNode;

public:
	//constructor
	SkirtPolyNode() : next(NULL){ }
	SkirtPolyNode(int edgeNode[2], int eleIdx, int skirtPolyNode[], int npt = 0,
				float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0)
	{
		this->edgeNode[0] = edgeNode[0];
		this->edgeNode[1] = edgeNode[1];

		this->eleIdx = eleIdx;
#ifdef _DEBUG
		assert(npt < MAX_SKIRT_POLY && npt >= 0);
#endif
		this->nSkirtNodes = npt;
		for (int i=0; i<npt; i++)
			this->skirtPoly[i] = skirtPolyNode[i];

		this->minQual = minQual;
		this->maxQual = maxQual;
		this->aveQual = aveQual;

		this->next = NULL;
		bDeleted = false;
	}

	SkirtPolyNode(int edgeNode[2])
	{
		SkirtPolyNode(edgeNode, -1, NULL);
	}

	~SkirtPolyNode(){}

	//setter/getter functions
	void setQual(float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0);
	void getQual(float *minQual, float* maxQual, float* aveQual);
	float getMinQual(){ return minQual; };

	void setElemQuals(float quals[], int nElem);
	float getMinElemsQual();

	void setSkirtPoly(int skirtPoly[], int npt);
	void getSkirtPoly(int skirtPoly[MAX_SKIRT_POLY], int* npt);	// npt 为实际长度

	void setEdge(int edgeNode[2], int eleIdx);
	void getEdge(int edgeNode[2], int* eleIdx);

	void setDeleted(bool bdel){ bDeleted = bdel; }
	bool isDeleted() { return bDeleted; };

public:
	SkirtPolyNode* getSkirtPolyNode(int edgeNode[2], int eleIdx, int skirtPolyNode[], int npt = 0,
		float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0, SkirtPolyNode* next = NULL);

	void insertAfter(SkirtPolyNode* polynode);
	SkirtPolyNode* removeAfter();

	
	/* print info */
	void printinfo();

private:
	int edgeNode[2];
	int eleIdx;
	float minQual;
	float maxQual;
	float aveQual;
	int nSkirtNodes;
	int skirtPoly[MAX_SKIRT_POLY];
	float polyEleQual[MAX_SKIRT_POLY];
	bool bDeleted;
	SkirtPolyNode* next;
};

//哈希节点，每个节点是一个链表，记录具有相同关键字的边壳信息
class HashPolyNode
{
	friend class HashPolyList;

public:
	// constructor
	HashPolyNode():next(NULL){}
	HashPolyNode(int ptIdx)
	{
		this->ptIdx = ptIdx;
		nPolyCnt = 0;
		pSkirtPolyList = NULL;
		next = NULL;
	}
	~HashPolyNode(){}

	// setter/getter functions
	int getPtidx(){ return ptIdx; }
	int getPolyCnt(){ return nPolyCnt; }
	SkirtPolyNode* skirtPolylist(){ return pSkirtPolyList; }

	void insertAfter(HashPolyNode* node);
	HashPolyNode* removeAfter();

public:
	SkirtPolyNode * addSkirtPoly(int edgeNode[2], int eleIdx, int skirtPolyNode[], float elemQuals[],
			int npt = 0, float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0);

	void deleteSkirtPoly(int edgeNode[2]);

	bool updateSkirtPoly(int edgeNode[2], int eleIdx, int skirtPolyNode[], float elemQuals[],
		int npt = 0, float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0);

	/*按照质量值从小到大的顺序调整polynode在链表中的位置*/
	void resortSkirtPoly(SkirtPolyNode* polynode);		//有问题！！！

	bool isPolyExist(int edgeNode[2]);

	SkirtPolyNode* getSkirtPoly(int edgeNode[2]);
	
	/* print info */
	void printinfo();

private:
	void addSkirtPoly(SkirtPolyNode* skirtPoly); //按质量顺序添加

private:
	int ptIdx;		//key
	int nPolyCnt;
	SkirtPolyNode* pSkirtPolyList;
	HashPolyNode* next;
};

//哈希链表，链表的每个节点具有一个唯一的关键字
class HashPolyList
{
public:
	HashPolyList():first(NULL){};
	~HashPolyList(){};

public:
	/* add */
	SkirtPolyNode * addSkirtPoly(int edgeNode[2], int eleIdx, int skirtPolyNode[], float elemQuals[], int npt = 0,
				float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0);

	void deleteSkirtPoly(int edgeNode[2]);

	void updateSkirtPoly(int edgeNode[2], int eleIdx, int skirtPolyNode[], float elemQuals[],
		int npt = 0, float minQual = 0.0, float maxQual = 0.0, float aveQual = 0.0);

	bool isPolyExist(int edgeNode[2]);

	//edgeNode为输入，其它为输出
	void getSkirtPoly(int edgeNode[2], int *eleIdx, int skirtPolyNode[], int *npt,
		float *minQual, float *maxQual, float *aveQual);

	float getMinTetQual(int beg, int end); //等待实现

	SkirtPolyNode* getSkirtPoly(int edgeNode[2]);

	void removeZeroSkirtPoly();

	/* print info */
	void printinfo();

private:
	HashPolyNode *first;
};

#define Type SkirtPolyNode*
#define BlockedArray_Type BlockedArray_SkirtPoly
#include "blockedarray.h"
#define LOG2_SKIRTPOLY_PER_BLOCK 12

class SkirtPolyHeap
{
public:
	SkirtPolyHeap():nHeapSize(0), nMaxHeapSize(LOG2_SKIRTPOLY_PER_BLOCK)
	{ 
		m_pSkirtPolyHeap.init(nMaxHeapSize);
		m_pSkirtPolyHeap.allocMemory(nMaxHeapSize); 
	}
	~SkirtPolyHeap(){};

public:
	int insertSkirtPoly(SkirtPolyNode *pPolyNode);
	SkirtPolyNode *getMinSkirtPoly();
	void rmMinSkirtPoly();
	void minHeapifyDown(int beg);
	void minHeapifyUp(int beg);

	bool isEmpty(){ return nHeapSize == 0; }
	bool isFull(){ return nHeapSize == nMaxHeapSize; };
	void allocHeap()
	{
		nMaxHeapSize = 2*nMaxHeapSize;
		m_pSkirtPolyHeap.allocMemory(nMaxHeapSize);
	}

	void printfinfo();

private:
	int nHeapSize;
	int nMaxHeapSize;
	BlockedArray_SkirtPoly m_pSkirtPolyHeap;
};