#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#define MAX_SKIRT_POLY 128 
#define _OUTPUT_LEVEL_III

enum QUALTYPE
{
	SMALL_DIHEDRAL_ANGLE,
	LARGE_DIHEDRAL_ANGLE,
	ALL_DIHEDRAL_ANGLES,
	DIHEDRAL_ANGLE_SINE,
	SOLID_ANGLE,
	ASPECT_RATIO
};

#if 1
class HashPolyNode;
class HashPolyList;

//����ڵ㣬��¼ÿ���ߵĿ���Ϣ
class SkirtPolyNode
{
	friend HashPolyNode;

public:
	//constructor
	SkirtPolyNode() : next(NULL){ }
	SkirtPolyNode(int edgeNode[2], int sideFace[2], float minQual = 0.0)
	{
		this->edgeNode[0] = edgeNode[0];
		this->edgeNode[1] = edgeNode[1];
		this->sideFace[0] = sideFace[0];
		this->sideFace[1] = sideFace[1];
		this->minQual = minQual;
		this->next = NULL;
		bDeleted = false;
	}

	~SkirtPolyNode(){}

	//setter/getter functions
	void setQual(float minQual = 0.0);
	void getQual(float *minQual);
	float getMinQual(){ return minQual; };

	void setEdge(int edgeNode[2]);
	void getEdge(int edgeNode[2]);
	void setFace(int sideFace[2]);
	void getFace(int sideFace[2]);

	void setDeleted(bool bdel){ bDeleted = bdel; }
	bool isDeleted() { return bDeleted; };

public:
	SkirtPolyNode* getSkirtPolyNode(int edgeNode[2], int sideFace[2], float minQual = 0.0, SkirtPolyNode* next = NULL);

	void insertAfter(SkirtPolyNode* polynode);
	SkirtPolyNode* removeAfter();

	
	/* print info */
	void printinfo();

private:
	int edgeNode[2];
	int sideFace[2];
	//int eleIdx;
	float minQual;
	//float maxQual;
	//float aveQual;
	//int nSkirtNodes;
	//int skirtPoly[MAX_SKIRT_POLY];
	//float polyEleQual[MAX_SKIRT_POLY];
	bool bDeleted;
	SkirtPolyNode* next;
};

//��ϣ�ڵ㣬ÿ���ڵ���һ��������¼������ͬ�ؼ��ֵı߿���Ϣ
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
	SkirtPolyNode * addSkirtPoly(int edgeNode[2], int sideFace[2], float minQual = 0.0);

	void deleteSkirtPoly(int edgeNode[2], int sideFace[2]);

	bool updateSkirtPoly(int edgeNode[2], int sideFace[2], float minQual = 0.0);

	/*��������ֵ��С�����˳�����polynode�������е�λ��*/
	void resortSkirtPoly(SkirtPolyNode* polynode);		//�����⣡����

	bool isPolyExist(int edgeNode[2], int sideFace[2], SkirtPolyNode ** pNode);

	SkirtPolyNode* getSkirtPoly(int edgeNode[2], int sideFace[2]);
	
	/* print info */
	void printinfo();

private:
	void addSkirtPoly(SkirtPolyNode* skirtPoly); //������˳�����

private:
	int ptIdx;		//key
	int nPolyCnt;
	SkirtPolyNode* pSkirtPolyList;
	HashPolyNode* next;
};

//��ϣ���������ÿ���ڵ����һ��Ψһ�Ĺؼ���
class HashPolyList
{
public:
	HashPolyList():first(NULL){};
	~HashPolyList(){};

public:
	/* add */
	SkirtPolyNode * addSkirtPoly(int edgeNode[2], int sideFace[2], float minQual = 0.0);

	void deleteSkirtPoly(int edgeNode[2], int sideFace[2]);

	void updateSkirtPoly(int edgeNode[2], int sideFace[2], float minQual = 0.0);

	bool isPolyExist(int edgeNode[2], int sideFace[2], SkirtPolyNode ** pNode);

	//edgeNodeΪ���룬����Ϊ���
	void getSkirtPoly(int edgeNode[2], int sideFace[2], float *minQual);

	SkirtPolyNode* getSkirtPoly(int edgeNode[2], int sideFace[2]);

	void removeZeroSkirtPoly();

	/* print info */
	void printinfo();

private:
	HashPolyNode *first;
};
#else
class HashPolyNode;
class HashPolyList;

//����ڵ㣬��¼ÿ���ߵĿ���Ϣ
class SkirtPolyNode
{
	friend HashPolyNode;

public:
	//constructor
	SkirtPolyNode() : next(NULL){ }
	SkirtPolyNode(int edgeNode[2], float minQual = 0.0)
	{
		this->edgeNode[0] = edgeNode[0];
		this->edgeNode[1] = edgeNode[1];
		this->minQual = minQual;
		this->next = NULL;
		bDeleted = false;
	}

	~SkirtPolyNode(){}

	//setter/getter functions
	void setQual(float minQual = 0.0);
	void getQual(float *minQual);
	float getMinQual(){ return minQual; };

	void setEdge(int edgeNode[2]);
	void getEdge(int edgeNode[2]);

	void setDeleted(bool bdel){ bDeleted = bdel; }
	bool isDeleted() { return bDeleted; };

public:
	SkirtPolyNode* getSkirtPolyNode(int edgeNode[2],float minQual = 0.0, SkirtPolyNode* next = NULL);

	void insertAfter(SkirtPolyNode* polynode);
	SkirtPolyNode* removeAfter();


	/* print info */
	void printinfo();

private:
	int edgeNode[2];
	//int eleIdx;
	float minQual;
	//float maxQual;
	//float aveQual;
	//int nSkirtNodes;
	//int skirtPoly[MAX_SKIRT_POLY];
	//float polyEleQual[MAX_SKIRT_POLY];
	bool bDeleted;
	SkirtPolyNode* next;
};

//��ϣ�ڵ㣬ÿ���ڵ���һ��������¼������ͬ�ؼ��ֵı߿���Ϣ
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
	SkirtPolyNode * addSkirtPoly(int edgeNode[2], float minQual = 0.0);

	void deleteSkirtPoly(int edgeNode[2]);

	bool updateSkirtPoly(int edgeNode[2], float minQual = 0.0);

	/*��������ֵ��С�����˳�����polynode�������е�λ��*/
	void resortSkirtPoly(SkirtPolyNode* polynode);		//�����⣡����

	bool isPolyExist(int edgeNode[2]);

	SkirtPolyNode* getSkirtPoly(int edgeNode[2]);

	/* print info */
	void printinfo();

private:
	void addSkirtPoly(SkirtPolyNode* skirtPoly); //������˳�����

private:
	int ptIdx;		//key
	int nPolyCnt;
	SkirtPolyNode* pSkirtPolyList;
	HashPolyNode* next;
};

//��ϣ���������ÿ���ڵ����һ��Ψһ�Ĺؼ���
class HashPolyList
{
public:
	HashPolyList():first(NULL){};
	~HashPolyList(){};

public:
	/* add */
	SkirtPolyNode * addSkirtPoly(int edgeNode[2], float minQual = 0.0);

	void deleteSkirtPoly(int edgeNode[2]);

	void updateSkirtPoly(int edgeNode[2], float minQual = 0.0);

	bool isPolyExist(int edgeNode[2]);

	//edgeNodeΪ���룬����Ϊ���
	void getSkirtPoly(int edgeNode[2], float *minQual);

	SkirtPolyNode* getSkirtPoly(int edgeNode[2]);

	void removeZeroSkirtPoly();

	/* print info */
	void printinfo();

private:
	HashPolyNode *first;
};
#endif

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