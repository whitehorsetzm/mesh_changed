/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ά����ͬ��Delaunay���������� (�汾�ţ�0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2005��9��15��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * ��ϵ��ʽ
 *   �绰��+86-571-87953165
 *   ���棺+86-571-87953167
 *   ���䣺zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
 
#ifndef __iso3d_h__
#define __iso3d_h__

#pragma warning(disable:4786)

#include <vector>
#include <list>
#include <map>
#include <stack>
#include "myvector.h"
#include "IntStack.h"
#include "IntIntMap.h"
#include "iso3d_define.h"
#include "space.h"
#include "spr.h"
#include "SMinterface.h"
#include "instentlist.h"
#include "blockedarray_instance.h"
#include "skirtpolylist.h"
#include "optdihangle.h"

//#pragma comment (lib, "dfor")

using namespace std;

//extern "C"  void _stdcall  OUTPUTPL3_FLOAT(const char *, int, int *, int *, int *, int *, float *, int *, int *);
//extern "C"  void _stdcall  OUTPUTPL3_DOUBLE(const char *, int, int *, int *, int *, int *, double *, int *, int *);

extern INTEGER g_nInnTotalCaviElem;
extern INTEGER g_nInnInstPoint;
extern INTEGER g_nInnTryPoints;
extern INTEGER g_nInnInitElems;

extern int g_nSPRCallings;
extern int g_nSPRCallSucc;
extern int g_nSPRPolySize;

extern int g_nGlobalSmoothingCallings;
extern int g_nGlobalSmoothingCallSucc;
extern int g_nEdContSmoothingCallings;
extern int g_nEdContSmoothingCallSucc;
extern int g_nEdSpltSmoothingCallings;
extern int g_nEdSpltSmoothingCallSucc;

/* ����ɹ�edge contraction(EC)�����ڵ��Ż���������ǰ��EC����ִ��ǰ�����������ı��� */
#define MAX_EC_QUALITY_RATIO_DIVIDE 20
extern double g_dECQualityRatio_MinRatio;
extern double g_dECQualityRatio_MaxRatio;
extern double g_dECQualityRatio_Divide[MAX_EC_QUALITY_RATIO_DIVIDE+1];
extern int g_nECQualityRatio_SuccStatistics[MAX_EC_QUALITY_RATIO_DIVIDE];
extern int g_nECQualityRatio_TotalStatistics[MAX_EC_QUALITY_RATIO_DIVIDE];

extern int g_nEdContCallings;
extern int g_nEdContCallSucc;
extern int g_nEdContCallNewP;

/* target type of shell transformation for quality improvement */
enum QI_ST_TARGET
{
	TARGET_CUR_SHELL = 0, /* attempt to improve the quality of the covering mesh of the current shell */
	TARGET_ONE_TETRA      /* attempt to remove one bad element, ps: this elem. may not be included in the current shell */
};

/*
 * declaration of class DTIso3D
 */
class DTIso3D
{
public:
	/*
	 * ���캯��
	 * constructors & destructor
	 */
	DTIso3D();
	virtual ~DTIso3D();

	INTEGER steinerPntBegin;
	INTEGER steinerPntEnd;
	map<int, int> sPntToEdge; 
	/* End of Add */

	
	/* �ͷ��ڴ� free memory */
	void freeMemory();

	/*
	 * ��ȡpls�ļ�
	 * read data from *.pls file
	 */
	int readPLS(const char* fname);
	
	/*
	 * ��дba3�ļ�
	 * read/write data from/to a *.pls file
	 */
	int readBA3(const char* fname);

	/*
	 * ��ȡpl3�ļ�
	 * read data from *.pls file
	 */
	int readPL3(const char* fname);
	
	/* �ؽ���Ԫ���ڹ�ϵ */
	int setupCellNeig();

	/* ---------------------------------------
     * ���м���(��������)���뵽.ele�ļ���
     * --------------------------------------*/
	int writeTetgenELE(const char* fname, const char *comments = NULL);
	
	/* -------------------------------------
	 * ����tetgen�ĵ�Delaunay����� 
	 * ----------------------------------*/
	int readTetgenELE(const char* fname, bool withOuterBox = false);

	/*
	 * �ɲ����ж���ԭ����BA3�ļ���PLS�ļ��е���Ϣ
	 */
	int initialInputInfFromParameters(
	/* ---------------------------�߽綨����Ϣ---------------------------------*/
	int nbp,				/* number of boundary points */
	double *lcoord,         /* coord. of boundary points , size = 3 * nbp */
	int nbf,								/* number of boundary surface */
	int *ibnd,              /* topu. of boundary entities, size = 4 * nbf (i1, i2, i3, igeom) */

	/* ---------------------------�ܶȿ�����Ϣ---------------------------------*/
	  /* background meshes */
	int       nbmn,			/* number of background mesh nodes */		 
	double    *ndinfo,      /* 4 * nbmn, (x, y, z, space) */
	int       nbme,			/* number of background mesh elements */
	int       *bm_topu,      /* 8 * nbme (i1, i2, i3, i4, iNeig1, iNeig2, iNeig3, iNeig4) */
	  /* grid sources */
	int       nsrc[3],      /* number of (point, line, triangle) sources */
	double    *src_coord    /* 6 *(nsrc[0] + 2 * nsrc[1] + 3 * nsrc[3]), 
			                   (x, y, z, space, innner radius, outer raidus)*/
	);

	int writeBA3(const char* fname);

	/* -----------------------------------------------------------
	 * �������Ƿ����ཻ 
	 * ---------------------------------------------------------*/
	bool isSelfIntersect();

	/*
	 * calculate outer box
	 */
	bool calcBox(MYPOINT minW, MYPOINT maxW, int startNode = INIT_NOD_NUM);
	
	/*
	 *	scale the geometry
	 */
	bool scaGeom(int startNode = INIT_NOD_NUM);
	
	/*
	 * scale the background mesh & source
	 */
	bool scaleSpace();

	/*
	 * �������������������任����������
	 * transform the background mesh from the WCS to NCS
	 */
	int wToN_BGMesh(BGMesh *pBGMesh);

	/*
	 * ����������ӹ�������任����������
	 * transform the background mesh from the NCS to WCS
	 */
	int nToW_BGMesh(BGMesh *pBGMesh);

	/*
	 * ��Դ����������任����������
	 * transform the source from the WCS to NCS
	 */
	int wToN_Source(Source *pSource);

	/*
	 * ��Դ�ӹ�������任����������
	 * transform the source from the NCS to WCS
	 */
	int nToW_Source(Source *pSource);

	/*
	 * calculate density values 
	 */
	bool calcDens();
	
	/*
	 * setup initial triangulation
	 */
	int setupInitTri();
	
	bool isOuterBoxEdge(INTEGER a, INTEGER b)
	{
		INTEGER t;
		if (a > b)
		{
			t = a;
			a = b;
			b = t;
		}

		return ((a < 8 && b < 8 && 
		        ((a == 0 && (b == 1 || b == 2 || b == 3 || b == 4 || b == 5 || b == 7)) ||
				 (a == 1 && (b == 2 || b == 5)) ||
				 (a == 2 && (b == 3 || b == 5 || b == 6 || b == 7)) ||
				 (a == 3 &&  b == 7) ||
				 (a == 4 && (b == 5 || b == 7)) ||
				 (a == 5 && (b == 6 || b == 7)) ||
				 (a == 6 &&  b == 7))));
		
	}
	/*
	 * find first elements including the point
	 */
	int findFirstEle(MYPOINT pnt, INTEGER *ele);
	
	/* -----------------------------------------------------------------------------------------------
	 * Ѱ�Ұ����ڵ������Base Elements������bases�Ĵ�С��������ǰȷ���ģ����ṩ������С��maxSize��,
	 * �Ա�֤������ȷ�ԡ�maxSize��ȱʡ����ֵΪMAX_SHELL_SIZE; maxSize = -1��ʾ��������Խ����
	 * -----------------------------------------------------------------------------------------------*/
	int findBaseElements(MYPOINT pnt, INTEGER bases[], int *baseSize, int maxSize = -1);

	/* -----------------------------------------------------------------------------------------------
	 * ���ϰ汾�ı߽������㷨�У����ڸ���������㹲����������ݵ�Ӱ�죬���������������Ŷ����������
	 * �߽��Ĳ��룬���ڴ󲿷������������õģ���ĳЩ�����£������������
	 * ------------------------------------------------------------------------------------------------
	 */
	int bndPntInst();
	
	
	/* -----------------------------------------------------------------------------------------------
	 * ���ǳ�����George�������Robust B-W Kernel����ɱ߽����룬�����ԭ�еĻ��ڵ��Ŷ��Ĳ������
	 * -----------------------------------------------------------------------------------------------*/
	int addInnerNode_Robust(INTEGER iNod);

	/* -----------------------------------------------------------------------------------------------
	 * ���ǳ�����George�������Robust B-W Kernel����ɱ߽����룬�����ԭ�еĻ��ڵ��Ŷ��Ĳ������
	 * -----------------------------------------------------------------------------------------------*/
	int addBoundNode_Robust(INTEGER iNod);

	/* 
	 * add boundary point
	 */
	int addBndPnt(MYPOINT pt, INTEGER iNod);
	 



	/* ----------------------------------------
	 |			  							 |					
	 | functions for tree searching          |
	 |										 |
	 * ----------------------------------------*/
	/* Clear the stack of tree searching
	 * Parameters:
	 *   void
	 * Return:
	 *   the size of the stack before clearing it
	 */
	void clrTreeSearEles();
	/*
	 * Get the current searching element
	 * Parameters:
	 *    void
	 * Result:
	 *    the top of the stack  the stack is not empty
     *    NULL_ElEM				otherwise
	 */
	INTEGER pickTreeSearch();
	/*
	 * Add an element for searching
	 * Paramerters:
	 *     iEle    the element added 
	 * Result:
	 *     the size of the stack after adding iEle
	 */
	void addTreeSearch(INTEGER iEle);

	/*
	 * Check if the stack is empty
	 * Parameters:
	 *     void
	 * Result:
	 *     true		empty stack
	 *     false    otherwise
	 */
	bool isTreeSearchEmpty();

	/* ----------------------------------------
	 |			  							 |					
	 | functions for deleted elements        |
	 |										 |
	 * ----------------------------------------*/
	/*
	 * Flag an element as deleted but not to add it to the array of deleted elements
	 * Parameters:
	 *      iEle	deleted element
	 * Result:
	 *      1									iEle is valid
	 *      0									otherwise
	 *
	 * Notes:
	 *   1. All deleted elements is flagged by reversing its radiuses
	 */
	int setDeleted(INTEGER iEle, int flag = 1);
	
	/*
	 * Flag an element as deleted & add it to the array of deleted elements
	 * Parameters:
	 *      iEle	deleted element
	 * Result:
	 *      the size of the array of deleted	iEle is valid
	 *      elements after deleting iEle
	 *      -1									otherwise
	 *
	 * Notes:
	 *   1. All deleted elements is flagged by reversing its radiuses
	 */
	int addDeleted(INTEGER iEle);
	/*
	 * Check if an element is deleted
	 * Parameters:
	 *       iEle	checked element
	 * Result:
	 *       true	deleted
	 *       false  otherwise
	 */
//	bool isDelEle(INTEGER iEle);
	bool isDelEle(Elem *pElem)
	{
		return pElem->iReserved2 & 0x1;
	}
	
	/*
	 * Recover all deleted elements to their initial states 
	 * Parameters:
	 *        void
	 * Result
	 *        the size of deleted elements
	 * Notes:
	 *   1. the function is called if an error happens while inserting ia node.
	 *   2. all deleted elements are recovered by re-reversing their radiuses.
	 */
	int recoverDelEles();

	/* ----------------------------------------
	 |			  							 |					
	 | functions for empty locations         |
	 |										 |
	 * ----------------------------------------*/
	/*
	 * Get ia location for the new element
	 * Parameters:
	 *    void
	 * Result:
	 *    ia location in the array of elements (m_pElems)
	 * Notes:
	 *    1. If the empty locations inner the array of elements 
	 *       (whose indices < m_nElems) are not consumed, get them 
	 *       first from the tail;
	 *       otherwise, get the location of the tail of the array of 
	 *       elements
	 */
	INTEGER getNewEleLoc();

	int erasTailLocs();
	/*
	 *	Update the array of empty locations 
	 *  Parameters:
	 *     void
	 *  Result:
	 *     the size of the array of empty locations after the updation
	 *  Notes:
	 *     1. The function is called while inserting ia node correctly
	 *     2. All the positions of deleted elements arised in the insertion 
	 *        of ia node could be used as the preferred candidate locations 
	 *        for new elements produced in next insertions
	 */
	int updateEmpLocs(bool eraseEmpSpace = true);
	/*
	 *  Recover the array of empty locations to the initial state before the
	 *  insertion of ia node
	 *  Parameters:
	 *      void
	 *  Result:
	 *      the size of the array of empty locations after the recovery
	 *  Notes:
	 *      1. The function is called if an error happens while inserting ia node
	 *      2. The free indicator m_nLocInd is moved to the end of the array of 
	 *         empty locations to perform the recovery
	 */
	int recoverEmpLocs();

	/*
	 * add an element directly into the array of empty locations
	 * CAUTION: the elements must be deleted
	 */
	int addEmpLoc(INTEGER iEle);
	
	

	/* ----------------------------------------
	 |			  							 |					
	 | functions for new elements            |
	 |										 |
	 * ----------------------------------------*/
	/*
	 * Create an element with given info.
	 * Parameters:
	 *      iLoc		location for the new element
	 *		form[]		forming points 
	 *      cen			center of the circumsphere
	 *		rad			square of the radius of the circumsphere
	 * Result:
	 *      the size of the array of elements
	 */
	INTEGER crtEle(INTEGER iLoc, INTEGER form[], MYPOINT cen, REAL rad);
	/*
	 * Recover the new elements to their initial states 
	 * Parameters:
	 *     void
	 * Result:
	 *     the size of the array of new elements 
	 * Notes:
	 *     1. The function is called if an error happens while inserting ia node
	 *     2. All new elements is flagged as deleted in the recovery
	 *     3. As m_nAftTail records the number of elements whose positions are
	 *        taken from the tail of the array of elements rather than the 
	 *        the empty locations formed by previous deleted elements, the m_nElems
	 *        should subtract m_nAftTail in the recoverry if m_nAftTail > 0
	 */
	int recoverNewEles();
	/*
	 * Clear the array of new elements & reset m_nAftTail as 0
	 * Parameters:
	 *    void
	 * Result
	 *    the size of the array of new elements before clearing it
	 */
	int clrNewEles();

	int addNewCreEles();

	/* ----------------------------------------
	 |			  							 |					
	 | functions for tested elements         |
	 |										 |
	 * ----------------------------------------*/
	/* Check if an element is tested
	 * Parameters:
	 *     iEle		checked element
	 * Result
	 *     true		tested
	 *     false    otherwise
	 * Notes:
	 *     1. All accessed elements are flagged as tested while 
	 *        traversing all elements by neighboring searching, which
	 *        prevents repetitive accesses & infinite searching
	 *     2. All tested element are flagged by setting their iReserved members 
	 *        as 1s (default values are 0s)
	 */
//	bool isTstEle(INTEGER iEle);
	bool isTstEle(Elem *pElem);

	/*
	 * Set the tested flag of an element
	 * Parameters:
	 *     iEle		the element
	 *     flag     tested flag. true tested; false otherwise
	 * Result
	 *     true		tested
	 *     false    otherwise
	 */
//	void enableEleTstFlag(INTEGER iEle);
	void enableEleTstFlag(Elem *pElem);
//	void disableEleTstFlag(INTEGER iEle);
	void disableEleTstFlag(Elem *pElem);

	/*
	 * Flag an elment as tested & add it to the array of tested elements
	 * Parameters:
	 *     iEle		tested element
	 * Result
	 *     the size of the array of tested 
	 *     elements after adding iEle              iEle is valid
	 *      -1								       otherwise
	 */
	void addTstEle(INTEGER iEle, Elem *pEle);
	/*
	 * Clear flags of all tested elements & clear the array of tested elements
	 * Parameters:
	 *    void
	 * Result:
	 *    the size of the array of tested elements before clearing it
	 * Notes:
	 *    1. A traversing operator to the array of elements could be implemented to
	 *       clear all flags of tested elements. However, we attemp to implement ia local
	 *       rather than global operator for the insertion of ia node. The array of tested
	 *       elements helps fulfill our goal.
	 */
	void clrTstEles();

	/* ----------------------------------------
	 |			  							 |					
	 | functions for cavity sides            |
	 |										 |
	 * ----------------------------------------*/
	/*
	 * Add ia cavity side with given info.
	 * Parameters:
	 *     i1		the start node of the side
	 *     i2		the end node of the side
	 *     iEle		the element including the side
	 *     iNei		the inicator of the neigbor who shares 
	                the face of iEle represented by the side 
	 *     iLast	the tail position of the list of cavity sides who share the
	 *              common start node i1
	 * Result:
	 *     the size of the array of cavity sides after adding the side
	 */
	int addCavSide(INTEGER i1, INTEGER i2, INTEGER iEle, int iNei, int iLast);
	/*
	 * Get info. about ia cavity side
	 * Parameters:
	 *     i1		the start node of the side
	 *     i2		the end node of the side
	 *    *iLast    returned variable, the tail position of the list of cavity sides 
	 *              who share the common start node i1
	 *	Resut:
	 *		> 0     the position of the side i1->i2 while it exits
	 *     0        no sides with the start point i1 exit
	 *    -1		the side i1->i2 doesnot exist, but sides with i1 as the start point exist
	 *  Notes:
	 *     1. All cavity sides with the common start point are organized as ia list, which
	 *        is connected by the iNei member of struct CavSide.
	 *     2. The head of the list is stored in the last 30 bits of the iReserved member of
	 *        the start node
	 *     3. If -1 is returned, iLast brings ia value represents the tail of the list
	 *     4. The valid positions of the array of cavity sides donot include 0, as 0 represents
	 *        ia case where no sides with i1 as the start point exist. So take care to the first
	 *        element of the array of cavity sides
	 */
	int getCavSide(INTEGER i1, INTEGER i2, int *iLast);
	/*
	 * Set the head of the list of cavity sides with i1 as the start point,
	 * which is stored in the last 30 bits of iReserved of node i1, as cs_ind
	 * Parameters:
	 *    i1		the start point of the side
	 *    cs_ind	the head of the list of cavity sides with i1 as the start point
	 * Return:
	 *    1			successful
	 *    0         otherwise
	 * Notes:
	 *    1. Take care to the sign bit
	 */
	int setCavSideInd(INTEGER i1, int cs_ind);
	/*
	 * Get the head of the list of cavity sides with i1 as the start point,
	 * which is stored in the last 30 bits of iReserved of node i1
	 * Parameters:
	 *    i1		the start point of the side
	 * Return:
	 *    1			successful
	 *    0         otherwise
	 * Notes:
	 *    1. Take care to the sign bit
	 */
	int getCavSideInd(INTEGER i1);


	/* ----------------------------------------
	 |			  							 |					
	 | functions for cavity nodes            |
	 |										 |
	 * ----------------------------------------*/
	/*
	 * Add ia cavity node
     * Parameters:
     *    i		the cavity node
     * Return:
     *    the size of the array 
	 *    of cavity nodes after adding i      if i is valid
     *    -1								  otherwise
     */
     int addCavNode(INTEGER i);
	/*
	 * Get the flag value of ia node, which is stored in the 2 bits ahead of the iReserved member
	 * Parameters:
	 *		i      the node
	 * Result:
	 *      flag value of the node    if i is valid
	 *      -1						  otherwise
	 * Notes:
	 *      1. The flag ranges in 0~3
	 */
	int getNodeFlag(INTEGER i);
	/*
	 * Set the flag value of ia node, which is stored in the 2 bits ahead of the iReserved member
	 * Parameters:
	 *		i      the node
	 *      flag   flag value
	 * Result
	 *      flag value of the node/*
	 * Get the flag value of ia node, which is stored in the 2 bits ahead of the iReserved member
	 * Parameters:
	 *		i      the node
	 * Result:
	 *      1     successful
	 *      0     otherwise
	 * Notes:
	 *      1. the flag ranges in 0~3
	 */
	int setNodeFlag(INTEGER i, int flag);
	/* Reset the iReserved of all nodes to 0
	 * Parameters:
	 *    void
	 * Result:
	 *     the size of all nodes flagged
	 * Notes:
	 *     1. All cavity nodes encountered while inserting ia node are flagged, which 
	 *        is ia superset of all nodes, and the cavity sides with them as start poings
	 *        exist. Therefore, to clear iReserved of the array of flagged nodes is 
	 *        able to clear the iReserved of all nodes whose iReserved members indicate
	 *        lists of cavity sides in their last 30 bits.
	 */        
	int clrNodeReserved();

	/*
	 *update neighboring info. of elements near the cavity
	 */
	bool updateNeigInfo();


	/*
	 * calcuate element parameters
	 * tv: triple volume
	 * cen: center of circumcenter
	 * rad: radius of circumcenter
	 */
	bool calcElePar(MYPOINT pnt[], REAL *tv, MYPOINT *cen, REAL *rad);

	/* ---------------------------------
	 * for disturbance point info.     |
	 * ---------------------------------*/
	bool isDisturbed(INTEGER iNod);
	bool addDistInfo(INTEGER iNod, MYPOINT pnt);
	bool setBackDistNode();

	/* ---------------------------------
	 * for normalizing the coordinates |
	 * ---------------------------------*/
	/*
	 * scale the coordinates range [minW, maxW] to ia range of [minN, maxN],
	 * and return according center & scale value
	 */
	bool scaleFactor(const MYPOINT minW, const MYPOINT maxW, 
					 const MYPOINT minN, const MYPOINT maxN);

	bool simpleScale(const MYPOINT minW, const MYPOINT maxW);
						
	/* 
	 * world cordinates to normalization coordinates 
	 */
	bool wToN(MYPOINT wp, MYPOINT *np); 
	
	/* 
	 * normalization cordinates to world coordinates 
	 */
	bool nToW(MYPOINT np, MYPOINT *wp);	

	bool wToN(double wf, double *nf);
	bool nToW(double nf, double *wf);
	

	/* ***************************************************
	 * functions for edge recovery
	 * **************************************************/
	bool isNodInc(INTEGER iNod, Elem *pElem);
	int getSphereSize(INTEGER iNod, INTEGER *iHint = NULL);
	int findSphere(INTEGER iNod, Sphere sph, int *nSph, INTEGER *iHint = NULL);  /* iHint�ǰ�����ĵ�һ����Ԫ��ΪNULLʱ���ѯ��ϣ���� */
	int findShell(INTEGER s, INTEGER e, Shell she, int *nShe, INTEGER *iHint = NULL); /* iHint�ǰ����ߵĵ�һ����Ԫ��ΪNULLʱ���ѯ��ϣ���� */
	bool isMeshEdge(INTEGER i1, INTEGER i2);
	bool isMeshFace(INTEGER i1, INTEGER i2, INTEGER i3);

	int findShellFaces(INTEGER s, INTEGER e, Shell she, int nShe, INTEGER thirdNodes[], INTEGER faceParents[], int *faceNum);
	int findShellFaces_InOrder(INTEGER s, INTEGER e, INTEGER iInitElem, INTEGER iInitNode, INTEGER thirdNodes[], INTEGER faceParents[], int *faceNum);
	INTEGER addNode(MYPOINT pnt, REAL space);
	INTEGER delNode(INTEGER iNod);
	int calcEdgTetsInt(INTEGER iEdg, PipeSrch* pSrch, PipeSrch* pAftSrch,
		Pipe pip, int *nPipSize);
	bool isValidEFPipel(int fc, int ec);
	int setPipel(INTEGER iEle, int nPip);
	int getPipel(INTEGER iEle);
	int findPipe(INTEGER iEdg, Sphere sph, int nSph, Pipe pip, int *size);
	int findPipe_Head(INTEGER iS, INTEGER iE, Sphere sph, int nSph, INTEGER elems[], int *size, INTEGER commNds[3]);
	int locDec(Pipe pip, int size, INTEGER stnBeg, INTEGER stnEnd);

	/* ����ɾ��Steiner�� */
	int removeEdgeSteinerPnt_CBR(int s, int e, int nStnBeg, Pipe pip, int size);
	int removeFaceSteinerPnt_CBR(int faceIdx, int iface[], int nStnBeg, int nStnEnd, Cluster clus, int nClus);

	/* -----------------------------------------------------------------
	 * ������չ��ǻ��ɾ��Steiner��
	 * ----------------------------------------------------------------*/
	inline int isSubTet(INTEGER iSubTet) 
	{
#ifdef _ONE_LONG_ARRAY
		assert(iSubTet >= 0 && iSubTet < m_nElems);
#else
		assert(iSubTet >= 0 && iSubTet < m_pElems.getArraySize());
#endif
		return m_pElems[iSubTet].iReserved2 & (1<<15);//0x20;
	}
	inline void flagSubTet(INTEGER iSubTet, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iSubTet >= 0 && iSubTet < m_nElems);
#else
		assert(iSubTet >= 0 && iSubTet < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iSubTet].iReserved2 &=~(1<<15)) :
					(m_pElems[iSubTet].iReserved2 |= (1<<15)); //0x20
	}
	inline void addSubTet(INTEGER iSubTet)
	{
		flagSubTet(iSubTet,true); 
		m_vecSubTets.push_back(iSubTet);
	} 
	int initialiseSubTets(Pipe pip, int nPip);
	int initialiseSubTets(Cluster cluss, int nClus);
	int initialiseSubTets(INTEGER elems[], int nElems);
	int cleanSubTets();
	int calcSubTetsNum();

	/* -------------------------------------------------------------
	 * ����Steiner��ָ��߽��ʱ�����ӷֹ��̿��ܻ�����������Ԫ
	 * ��ˣ��������¼��Щ�������Ԫ��������ͨ���Ŷ�������֤��Ԫ
	 * �������
	 * -------------------------------------------------------------*/
	/* -----------------------------------------------------------------
	 * ������չ��ǻ��ɾ��Steiner��
	 * ----------------------------------------------------------------*/
	inline bool isNegativeVElem(INTEGER iElem) 
	{
#ifdef _ONE_LONG_ARRAY
		assert(iElem >= 0 && iElem < m_nElems);
#else
		assert(iElem >= 0 && iElem < m_pElems.getArraySize());
#endif
		return m_pElems[iElem].iReserved2 & ~(1<<16);//0x40;
	}
	inline void flagNegativeVElem(INTEGER iElem, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iElem >= 0 && iElem < m_nElems);
#else
		assert(iElem >= 0 && iElem < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iElem].iReserved2 &=~(1<<16))  :
					(m_pElems[iElem].iReserved2 |= (1<<16)); //0x40
	}
	inline double calcElemVol(INTEGER iElem)
	{
		INTEGER *pForm = NULL;

#ifdef _ONE_LONG_ARRAY
		assert(iElem >= 0 && iElem < m_nElems);
#else
		assert(iElem >= 0 && iElem < m_pElems.getArraySize());
#endif

		pForm = m_pElems[iElem].form;

		return GEOM_FUNC::orient3d(m_pNodes[pForm[0]].pt, m_pNodes[pForm[1]].pt, 
			m_pNodes[pForm[3]].pt, m_pNodes[pForm[2]].pt);
	}
	int removeNegativeVElem_Edge(INTEGER iEdg);
	int cleanNegativeVElems()
	{
		for (int i = 0; i < m_vecNegativeVElems.size(); i++)
			flagNegativeVElem(m_vecNegativeVElems[i], false);
		m_vecNegativeVElems.clear();

		return 1;
	}

	int removeUndesirablePipels_SPR_V1(Sphere sphPipe, int nSphPipe, INTEGER delPips[], int nDelPips, INTEGER badPips[][4], int nBadPips);
		
	int updateNeig(INTEGER iEle, INTEGER iNeiBef, INTEGER iNeiAft);
	int updateNeig(INTEGER iEle, INTEGER f1, INTEGER f2, INTEGER f3, 
		INTEGER iNeiAft);
	int specCase1(Pipe pip, int size);
	int specCase2(Pipe pip, int size);
	int edgCase(Pipel* pPipel, Pipe pip, int size);
	int douEdgCase(Pipel* pPipel, Pipe pip, int size);
	int nodFacCase(Pipel* pPipel, Pipe pip, int size);
	int edgFacCase(Pipel* pPipel, Pipe pip, int size);
	int douFacCase(Pipel* pPipel, Pipe pip, int size);
	INTEGER findPipNeig(INTEGER iEle, int iNeig, INTEGER f1, INTEGER f2, INTEGER f3, 
		Pipe pip, int size);
	int findPipeSN(INTEGER iABCD, int ic, int *sn, Pipe pip, int size);
	int setPipelSN(Pipel *pPipel, int sn);
	int getPipelSN(Pipel* pPipel, int *sn);
	INTEGER abstEdges();
	INTEGER getEdge(INTEGER i1, INTEGER i2, INTEGER *iLast);
	INTEGER addEdge(INTEGER i1, INTEGER i2, INTEGER iLast);
	int resetSurEdgHash(INTEGER maxSize);
	INTEGER abstNodFirEle();
	int checkEdgeSP(INTEGER *nLost = NULL, INTEGER *nSPs = NULL, INTEGER *nMaxSPs = NULL);
	int checkFaceSP(INTEGER *nLost = NULL, INTEGER *nCut = NULL, INTEGER *nSPs = NULL, INTEGER *nMaxSPs = NULL);
	int recvEdges();
	int recvEdge(INTEGER iEdg, Sphere sph, int nSph);
	int recvEdges_NoStnPnt();
	int recvEdge_NoStnPnt(INTEGER iEdg, Sphere sph, int nSph, int *nInit = NULL);
	void setBndProtectionFlag(bool bFlag)
	{
		m_bBndProtection = bFlag;
	}
	int updateInstEntsStatus();
	int classifyBndEdges(INTEGER *nLost, INTEGER *nSPs = NULL);
	int classifyBndFaces(INTEGER *nLost, INTEGER *nSPs = NULL);
	int recvBndEdge_NoStnPnt();
	int recvBndEdge_NoStnPnt_V2();
	int recvBndFace_NoStnPnt_V2();
	int recvBndEdge_NoStnPnt_V2_OneLoop();
	int recvBndFace_NoStnPnt_V2_OneLoop();
	int findAllInstEntitiesOfBE(INTEGER iEdg, InstEntList& listInstEnts, bool bRecordInstEnts = false);
	int findAllInstEntitiesOfBF(INTEGER iFac, InstEntList& listInstEnts, bool bAllEdgeRecved = false, bool bRecordInstEnts = false);
	int removeInstEntOfBE(INTEGER iEdg, InstEntList& listInstEnts, int iEnt);
	int removeInstEntOfBE_V2(INTEGER iEdg, InstEntList& listInstEnts, int iEnt);
	int removeInstEntOfBF_V2(INTEGER iFac, InstEntList& listInstEnts, int iEnt);
	int swap23(INTEGER indices[3], INTEGER iElems[3], bool tvCheck = true);
	int swap32(INTEGER indices[2], INTEGER iElems[3], bool tvCheck = true);
	int swap32_LabelNodes(INTEGER indices[2], INTEGER iElems[3], INTEGER skirtNodes[3]);
	int removeAnEdge_DFS(INTEGER a, INTEGER b, INTEGER lastFacNode, Shell she, int *pnShe, INTEGER cstNodes[], int cstType, 
		int iEnt, InstEntList& listInstEnts, InstEntList& listDependentEnts, int level, int iInitNode = -1);
	int removeAnEdge_BFS(INTEGER a, INTEGER b, Shell she, int *pnShe, INTEGER cstNodes[], int cstType, 
		int iEnt, InstEntList& listInstEnts, InstEntList& listDependentEnts, int level, int iInitNode = -1);
	int remeshAShell_SPR(INTEGER a, INTEGER b, Shell she, int *pnShe, INTEGER iNodB, INTEGER iNodE, 
		int iEnt, InstEntList& listInstEnts, InstEntList& listInflatedEnts, InstEntList& listDependentEnts);
	bool isEar_SkirtDigonal(INTEGER a, INTEGER b, int faceNum, INTEGER faceNodes[], INTEGER faceNodesDep[], InstEntList& listDependentEnts,
		INTEGER iiP, INTEGER iiC, INTEGER iiN, bool *bSwap44Zero);
	int prevUnClippedNode(int iiC, int faceNum, int faceNodesClp[]);
	int nextUnClippedNode(int iiC, int faceNum, int faceNodesClp[]);
	int includeDependentEdgeCnt(INTEGER iElem, InstEntList& listDependentEnts);
	int remeshAShell_EarRemoval(INTEGER a, INTEGER b, Shell she, int *pnShe, INTEGER cstNodes[], int cstType, 
		int iEnt, InstEntList& listInstEnts, InstEntList& listDependentEnts, INTEGER iInitNode = -1);
	int remeshAShell_DynProg(INTEGER a, INTEGER b, INTEGER c, Shell she, int *pnShe, INTEGER cstNodes[], int cstType, 
		int iEnt, InstEntList& listInstEnts, InstEntList& listDependentEnts, INTEGER iInitNode = -1);
	int swap23Recursive(INTEGER a, INTEGER b, INTEGER matrix_k[][MAX_SKIRT_POLY_SIZE], int ii, int jj, 
		int faceNum, INTEGER faceNodes[], INTEGER faceParents[], InstEntList& listDependentEnts);
	int calcDepNodeRef_Recursive(INTEGER matrix_k[][MAX_SKIRT_POLY_SIZE], int i, int j, int l, int faceNum, int *nodeRef);
	int qualArrayRecursive(int faceNum, INTEGER matrix_k[][MAX_SKIRT_POLY_SIZE], float *matrix_qualFaces, 
		int ii, int jj, float qualityArray[], int *size);

	/* ��¼1����Ԫ�����ü�����Ŀǰֻ�õ�iReserved2���2λ��18/19λ */
	int elemCounterIncre(int iElem)
	{
		int cnt;
		Elem *pElem = &m_pElems[iElem];
		assert(!isDelEle(pElem));
		cnt = elemCounter(iElem) + 1;
		assert(cnt <= 6);
		pElem->iReserved2 &= ~(0x7 << 17); /* ��18/19/20λ���Ϊ 0 */
		pElem->iReserved2 |= (cnt << 17);	/* ����ֵ����18/19/20λ�� */
		return cnt;
	}
	int elemCounterDecre(int iElem)
	{
		int cnt;
		Elem *pElem = &m_pElems[iElem];
		assert(!isDelEle(pElem));
		cnt = elemCounter(iElem) - 1;
		assert(cnt >= 0);
		pElem->iReserved2 &= ~(0x7 << 17); /* ��18/19/20λ���Ϊ 0 */
		pElem->iReserved2 |= (cnt << 17);	/* ����ֵ����18/19/20λ�� */
		return cnt;
	}
	int elemCounter(int iElem)
	{
		return (m_pElems[iElem].iReserved2 >> 17) & 0x7;
	}
	/* �����ཻ�߻���Ƭ�б�(ֻ������) */
	int updateInstEntities(InstEntList& listInstEnts);
	int setEdgeConstraintsofSPR(const polyhedron& poly, 
		INTEGER i1, INTEGER i2, InstEntList& listInstEnts, int iActEnt, 
		std::map<INTEGER,int>& mapG2L, INTEGER l2g[], bool bActEntAllowed = false);

	int recvEdges_LessStnPnt();
	int recvEdge_OneByOne(INTEGER iEdg, Sphere sph, int nSph, int *nInit, int *nActl);
	int recvEdgeByTopuImpr(INTEGER elems[], int elemSize, INTEGER s, INTEGER e, int level, 
		INTEGER intEdges[][3], int nIntE, INTEGER intFaces[][4], int nIntF);
	int recvEdge_LessStnPnt_DivideAndConquer(INTEGER iEdg, Sphere sph, int nSph, int *nInit, int *nActl, int *nUndo);
	int recvEdge_LessStnPnt(INTEGER iEdg, Sphere sph, int nSph, int *nInit = NULL, int *nActl = NULL, int *nUndo = NULL);

	void setBERecovered(INTEGER iEdg, bool flag)
	{
		flag ? 	m_pSurEdgs[iEdg].subHead |= 0x1 :
			    m_pSurEdgs[iEdg].subHead &= ~(0x1);
	}
	void setBESubHead(INTEGER iEdg, INTEGER subHead)
	{
		m_pSurEdgs[iEdg].subHead &= 0x1;	/* ����ǰ1λ���������λ��0 */
	    m_pSurEdgs[iEdg].subHead |= (subHead << 1); /* ��subHead��䵽��������λ */
	}
//	void setBEInstEntCnt(INTEGER iEdg, int cnt)
//	{
//		m_pSurEdgs[iEdg].subHead = cnt;
//	}
	bool isBERecovered(INTEGER iEdg)
	{
		return m_pSurEdgs[iEdg].subHead & 0x1;
	}
	INTEGER getBESubHead(INTEGER iEdg)
	{
		return m_pSurEdgs[iEdg].subHead >> 1;
	}
//	int getBEInstEntCnt(INTEGER iEdg)
//	{
//		return m_pSurEdgs[iEdg].subHead;
//	}
	INTEGER begNdOfBE(INTEGER iEdg)
	{
		return m_pSurEdgs[iEdg].iStart;
	}
	INTEGER endNdOfBE(INTEGER iEdg)
	{
		return m_pSurEdgs[iEdg].iEnd;
	}

	void setBFFullRecovered(INTEGER iFac, bool flag)
	{
		flag ? m_pSurTris[iFac].subHead |= 0x1 :
			   m_pSurTris[iFac].subHead &= ~(0x1);
	}
	void setBFInteriorlyRecovered(INTEGER iFac, bool flag)
	{
		flag ? m_pSurTris[iFac].subHead |= 0x2 :
			   m_pSurTris[iFac].subHead &= ~(0x2);
	}
	void setBFSubHead(INTEGER iFac, INTEGER subHead)
	{
		m_pSurTris[iFac].subHead &= 0x3;	/* ����ǰ2λ���������λ��0 */
	    m_pSurTris[iFac].subHead |= (subHead << 2); /* ��subHead��䵽��������λ */
	}

	bool isBFFullRecovered(INTEGER iFac)
	{
		return m_pSurTris[iFac].subHead & 0x1;
	}
	bool isBFInteriorlyRecovered(INTEGER iFac)
	{
		return m_pSurTris[iFac].subHead & 0x2;
	}
	INTEGER getBFSubHead(INTEGER iFac)
	{
		return m_pSurTris[iFac].subHead >> 2;
	}

	/* ***************************************************
	 * functions for face recovery
	 * **************************************************/
	int recvFaces();
	int recvFace(INTEGER iFac, Sphere sph, int nSph);
	int recvFaces_NoStnPnt();
	int recvFace_NoStnPnt(INTEGER iFac, Sphere sph, int nSph);
	int ouIn(INTEGER iEle, INTEGER iTri, int *iLost);
	int setOuIn(INTEGER iEle, int flag);
	int getOuIn(INTEGER iEle);
	int pntLn(MYPOINT pnt, MYPOINT ln1, MYPOINT ln2);
#ifdef _NULL_DEFINED
	int findFirCluEdg(INTEGER iFac, Sphere sph, int nSph, 
				CluEdgArr edgs, int *nEdg, INTEGER *iEle);
#endif
	int creClus(Clusterel* pClus, INTEGER iFac,
		CluEdgArr cedgs, int *nEdg);
	int checNewEdg(Clusterel* pClus, int ec, INTEGER iFac, 
		INTEGER mi, INTEGER ma,
		CluEdgArr cedgs, int *nEdg);
	bool isCluEdgChked(INTEGER es, INTEGER ee, 
		INTEGER *iNod, CluEdgArr cedgs, int nEdg);
	int addCluEdg(INTEGER es, INTEGER ee, 
		INTEGER iNod, CluEdgArr cedgs, int *nEdg);
	int addClusCands(Clusterel* pClus, ClusCand cands, int *nCand);
	bool isFRChkEle(INTEGER iEle);
	bool setFRChkEleFlag(INTEGER iEle, bool flag);
	int addFRChkEle(INTEGER iEle);
	int clrFRChkEles();
	bool isCPClusNeig(INTEGER iCand, Cluster cluss, int nClus);
	int coPlanCase(Clusterel *pClus, INTEGER iFac, VECTOR normf, DesFacet *pDesFacet);
	int oneEdgCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, DesFacet *pDesFacet);
	int twoEdgCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, DesFacet *pDesFacet);
	int thrEdgCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, DesFacet *pDesFacet);
	int fouEdgCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, DesFacet *pDesFacet);
	INTEGER findClusNeig(INTEGER iEle, int iNeig, 
		INTEGER f1, INTEGER f2, INTEGER f3, 
		Cluster cluss, int nClus);
	int findClusSN(INTEGER iABCD, int ic, int *sn, Cluster cluss, int nClus);
	int setClusSN(Clusterel *pClus, int nFac, int sn);
	int getClusSN(Clusterel* pClus, int nFac, int *sn);
	int zeroSCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus,
		int id, int ibot[3], int sns[3], DesFacet *pDesFacet);
	int oneSCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus,
		int id, int ibot[3], int sns[3], DesFacet *pDesFacet);
	int doubSCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus,
		int id, int ibot[3], int sns[3], DesFacet *pDesFacet);
	int triSCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus,
		int id, int ibot[3], int sns[3], DesFacet *pDesFacet);
	int zeroNCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, int iabcd[4], int sns[3], DesFacet *pDesFacet);
	int oneNCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, int iabcd[4], int sns[3], DesFacet *pDesFacet);
	int doubNeigNCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, int iabcd[4], int sns[3], DesFacet *pDesFacet);
	int doubOppoNCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, int iabcd[4], int sns[3], DesFacet *pDesFacet);
	int triNCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, int iabcd[4], int sns[3], DesFacet *pDesFacet);
	int fouNCase(Clusterel *pClus, INTEGER iFac, VECTOR normf,
		Cluster cluss, int nClus, int iabcd[4], int sns[3], DesFacet *pDesFacet);

	
	/* 
	 * �ָ�һ����Ƭʱ����Ƭ�ı߿����Ѿ�����ϣ�����Щ���ϵĴ�ϵ��
	 * �������ڼ���ʱ������ֵ������Ƭ�ָ�ʱ��������Щ�Ѿ�������
	 * �������ֵ����ֵ����ʱ����Ҫ�����ߺ�����ཻ�㷨�����ܳ���
	 * �����۵����󣬴Ӷ�Ӱ���жϣ���ˣ����ǿ��Ը��ݱߵĽڵ���
	 * ��Ĺ�ϵ׼ȷ�ƶϳ��ߺ����һЩ��ϵ���������㣩���Ӷ���֤�㷨
	 * �Ľ�׳��
	 * ���º���ͨ���жϸ����Ľڵ㣨��ţ��Ƿ�����Ƭ���γɽڵ�����
	 * ��Ͻڵ����жϽڵ����Ĺ�ϵ
	 * Edges of a facet might be splitted before recovering it. Floating-point
	 * roundoff erros is appended to the coordinates of splitting points. Therefore,
	 * this kind of errors cound accumulate while using these coordinate values in
	 * new floating-point computations such as the check of the relation between
	 * an edge & a facet. The accumulated floating-point errors could finally reverse
	 * our judgement. 
	 * The following function check if a node is in a facet by compare the ID of the node
	 * & all node IDs of the facet, including the spliting nodes in edges belonging to it.
	 * The pure integer operations will avoid the floating-point roundoff errors described
	 * above.
	 */
	bool isNodeOfFacet(INTEGER iNod, INTEGER iFac);

	/* 
	 * ��������Ƭ��Ϣ
	 */
	int updateSubFacets(INTEGER iFacet, INTEGER iDesFacet);

	/* 
	 * �߽���Ƭ���ֽ�ɶ������Ƭ��ÿ������Ƭ����2�������嵥Ԫ��
	 * �ж���2��Ԫ����������ʱ��Ҫ����ȷ������Ƭ�����򣨳���/���⣩��
	 * ����ԭʼ��Ƭ����������֪�ģ���˹ؼ���ȷ������Ƭ��ԭʼ��Ƭ
	 * �����һ���ԡ�
	 *
	 * ��ͳ�ķ��������ڸ������㣬���ֱ��������Ƭ�ķ����ԭʼ��Ƭ�ķ���
	 * ������߷�����ͬ������Ϊ����Ƭ��ԭʼ��Ƭ����һ�¡�������Ƭ�ӽ���
	 * һ��ֱ��ʱ��������������
	 * 
	 * ���õİ취�Ǳ��⸡�����㣬���߼����£�
	 * 1. ȷ��ĳ������Ƭ�����򣨻����������㣩
	 * 2. ��������Ƭ�����ڹ�ϵȷ��������������Ƭ������
	 */
	int typeAttachedEles(INTEGER iFacet, INTEGER iDesFacet);

	/* ----------------------------------------------------
	 * Ѱ�Һ�һ�����ָ���Ƭ�ཻ�����е�Ԫ
	 * iFac: ��Ƭ���
	 * sph/nSph:  ��Ƭ��1���˵����
	 * cluss: ���к���Ƭ�ཻ�ĵ�Ԫ��������һ���ߺ���Ƭ�ཻ���������Σ���
	 * --------------------------------------------------*/
	int findCluster(INTEGER iFac, Sphere sph, int nSph,
		Cluster cluss, int *size);

	/* --------------------------------------------------------------------
	 * ����Ƭ�ı�δ���ָ�ʱ��Ѱ�����а������ָ���Ƭ�ĵ�Ԫ���㷨
	 * �����б��Ѿ��ָ������µ�FindCluster�㷨�ǲ�һ����
	 * -----------------------------------------------------------------*/
	int findCluster_NoEdgeRecved(INTEGER iFac, Sphere sph, int nSph,
		Cluster cluss, int *size);

	/*
	 * decompose all clusterels & 
	 * set types of elements sharing boundary faces
	 */
	int decClus(Cluster cluss, int nClus, INTEGER iFac, DesFacet *pDesFacet);

	/* flag all elements as outer or inner */
	int typeEles();

	/* flag all outer elements as deleted */
	int typeDelEles();

	/*
	 * remove all redundant elements
	 */
	int rmvOuterEles();

	/* ***************************************************
	 * functions for insertion of inner points
	 * **************************************************/
#ifdef _ERROR_CHK
	INTEGER getEleCntWithInnerNods();
	INTEGER getUndiscNodCnt(INTEGER iNodS, INTEGER iNodE);
#endif /* _ERROR_CHK */	
	int beforeInnerPntInst();
	int innerPntInst();
	INTEGER autCreNodes();
	
	/* ����һ�����������洢�´����Ľڵ� ��m_pCreateNodes)
	 * bEraseOldMem		�Ƿ�ɾ���ɵ��ڴ棨������С�������ڵ�����
	 */
	INTEGER autCreNodes_NewArr(bool bEraseOldMem = true);
	INTEGER getPrtEle(INTEGER iNod);
	int setPrtEle(INTEGER iNod, INTEGER iEle);
	int setCldNod(INTEGER iEle, INTEGER iNod);
	int setDiscard(INTEGER iNod, bool bDiscard);
	bool isDiscard(INTEGER iNod);
	INTEGER getCldNod(INTEGER iEle);
	int addInnerNode(INTEGER iNod);
	int bakEraseNod(INTEGER iNodS, INTEGER *iNodE);
	int distNode(INTEGER iNod);
	INTEGER findSwapNode(INTEGER iNod, INTEGER *iNodE);
	int swapNode(INTEGER iNod, INTEGER iSwap);
	int assiNode(INTEGER iNod, INTEGER iAssi);
	INTEGER swapNodSpacs(INTEGER iNod, INTEGER iDegS, INTEGER iDegE);
	int updateDiscNods();

	int clrNodAndEleRsvr();

	/*
	 * remove deleted elements & redundant nodes
	 */
	int rmvNodsAndEles();

	/*
	 * fill the parent domains of all surface
	 */
	int fillParents(bool bConformal = false);
	INTEGER findParent(INTEGER i1, INTEGER i2, INTEGER i3);
	INTEGER findParent_BruteForce(INTEGER i1, INTEGER i2, INTEGER i3);

	int evalSphereSize();
	
	/*
	 * Edit surface grids
	 */
	int editSurfs();

	/* ***************************************************
	 * Ϊһ�±߽�ָ����㷨��ĺ���
	 * functions for constrained boundary recovery
	 * **************************************************/
	/* 
	 * ��ӡԼ���߽�ָ���Ϣ
	 */
	int printConformalRecvInfo();

	/* 
	 * Ϊ��ȫ�ָ����ƻ��߽�������׼�� 
	 * prepare to recover destroyed boundaries (edges/facets)
	 */
	int prepareRecvDesBnds();

	/* �ָ����б��ƻ��ߣ����ϵĵ㡡recovery all extra-added node in the destroyed edges & facets */
	int recvDesEdges();
	int recvDesFacets();

	/* 
	 * Node��iReserved��ĵ�1λΪ�߽�ڵ��־������λ��ʾ�����ڵ��ĳ����Ԫ��
	 * �߽�ڵ����
	 * (1) ��ʼ��������ڵ�;
	 * (2) ִ��һ�±߽�ָ��㷨ʱ�ڱ���ߺ���Ƭ������ĸ����㡣
	 * The iReserved member of a Node object is divided into two parts, 
	 * the first bit flags if the node is a boundary node, and the other
	 * bits index an element(tetrahedron) including the node.
	 * The boundary nodes is composed of both (1) nodes of the surface mesh,
	 * and (2) extra added in the stage of conformal boundary recovery.
	 */
	bool isBndNode(INTEGER iNod);
	int setBndNode(INTEGER iNod, bool bBnd);
	INTEGER firstSph(INTEGER iNod);
	int setFirstSph(INTEGER iNod, INTEGER iEle);

	/* 
	 * ��������������ӵĸ�����(Extra-Added Node in the surface Edge (EANE))
	 * ������
	 *     iNod			  ����������
	 *     iFac1 & iFac2  ��������ߵ�����������Ƭ
	 *     iPre & iNxt	  ��ϱ���ߵ�һϵ�е����һ�����У�iPre��iNxtΪiNod
				 		  ����������е�ǰ��ڵ㣬���������ǻ��ǻ���ΪԼ����
	 * remove an Extra-Added Node in the surface Edge (EANE)
	 * Parameters
	 *     iNod				index of the extra-added node
	 *     iFac1 & iFac2	two surface facets including the surface edge
	 *     iPre & iNxt		A sequence of nodes is constructed while splitting an surface edge for conformal boundary recovery. iPre & iNxt is the previous & next node of iNod in the sequence, respectively. 
							 The link with them as end points is a constraint of the shell loop
	 */
	int removeEANE(INTEGER iNod, INTEGER iDesEdg, INTEGER iPre, INTEGER iNxt);

	/* 
	* ��������Ƭ�����ӵĸ�����(Extra-Added Node in the surface Facet (EANF))
	* remove an Extra-Added Node in the surface Facet (EANF)
	*/
	int removeEANF(INTEGER iNod, INTEGER iDesFac); 

	/*
	 * ���㹲���ƻ��ߵ���������ķ�����
	 * calc. the normals of two facets sharing the destroyed edge
	 */
	int calcFacetNormalOfDesEdge(Sphere halfOu, int nHalfOu, 
		INTEGER iNod1, INTEGER iNod2, VECTOR vf1, VECTOR vf2);

	/* ��һ����ֳ��������� divided a sphere into two half spheres */
	int divideSphere(Sphere sph, int nSph, Sphere halfOu, int *nHalfOu,
	Sphere halfIn, int *nHalfIn);

	/* 
 	 * ����ֽ����λ��
	 * ������
	 *     sph_cen	    ��������
	 *     norm			�ֽⷽ������
	 *     p_init		���룺��ʼѡ��λ�ã���������շ���Լ�������ķֽ��
	 *     deta_init	���룺��ʼ��norm����չ�������: ������չ��
	 *     half			����
	 *     nHalf		�����е�Ԫ��Ŀ
 	 * �豣֤�ֽ����λ�ö԰����ǿ��ӵ�
	 * calculate one point resulted from splitting the extra-added node sph_cen
	 * Parameters��
	 *     sph_cen	    index of the sphere center
	 *     norm			normal of the splitting direction
	 *     p_init		Input:  initially given position; 
						Output: finial qualified position
	 *     deta_init	Input:  initial stretch length along the norm
						Output: final stretch length
	 *     half			a half sphere
	 *     nHalf		number of eles.(tetrahedron) of the half sphere
	 * the finial position must be visible for the half ball
	 */
	int calcSplit(INTEGER sph_cen, VECTOR norm, MYPOINT p_init, 
	REAL *deta_init, Sphere half, int nHalf);

	/* 
	 * �õ��ǻ��ı߽������ 
	 * obtain an array of boundary edges of the Shell Loop
	 */
	int obtainShellEDs(INTEGER sph_cen, Sphere halfOu, int nHalfOu,
	ShellED shellEDs[], int *nShellEDs);

	/* 
	 * �Կǻ��ı߽�������������
	 * sort the array of boundary edges of the Shell Loop
	 */
	int sortShellEDs(ShellED shellEDs[], int nShellEDs);

	/*
	 * ���ǻ�һ��Ϊ��
	 * divide the shell loop (described by its boundary edges) into
	 * two ones
	 */
	int divideShellEDs(INTEGER iNod1, INTEGER iNod2, 
	ShellED shellEDs[], int nShellEDs, 
	ShellED lftSEDs[], int *nLftSEDs,
	ShellED rgtSEDs[], int *nRgtSEDs,
	int *i1, int *i2);

	/* ��ȡ�����Ͻڵ����� */
	int obtainShellNDs(ShellED shellEDs[], int nShellEDs,
		INTEGER shellNDs[], int *nShellNDs);

	/* ��ȡ���ߣ����ߴ��ָ���iEdg�������ߣ���ı�� */
	int obtainShellNDs_ColinIdx(INTEGER iEdg, INTEGER iPre, INTEGER iNxt,
		INTEGER lftSNDs[], int nLftSNDs, 
		INTEGER rgtSNDs[], int nRgtSNDs,
		int lftSNDs_ColinIdx[2][MAX_SHELL_ED], int rgtSNDs_ColinIdx[2][MAX_SHELL_ED],
		int *nLftColIdx, int *nRgtColIdx);

	int obtainShellNDs_ColinIdxInEdg(INTEGER iEdg, INTEGER shellSNDs[], int nSNDs, int colIdx[], int *colNDSize);


	/* ���ǻ��ǻ� triangulate a shell loop */
	int triShellLoop(INTEGER shellNDs[], int nShellNDs, VECTOR norm,
	int prt[], int forms[], int neigs[], int nf, int shellNDs_ColinIdx[][MAX_SHELL_ED], int colNDSize); /* ������3���� */


	/* �����µ�Ԫ add new elements */
	int addNewEles(INTEGER iFirst, INTEGER iSecond, 
		ShellED shellEDs[], int nShellEDs,
		int prt[], int forms[], int neigs[],
		INTEGER new_out_tet[], INTEGER new_in_tet[], int nf);

	/* 
	 * ���¹���߽�������ĸ���Ԫ���ڽӹ�ϵ
	 * update neighboring info. of four new elements formed near the shared
	 * boundary edges
	 */
	int updateNeigInfoOfSharedSED(INTEGER lftOu, INTEGER lftIn, 
	INTEGER rgtOu, INTEGER rgtIn);

	/* 
	 * ����������д��ΪPLS�ļ� (��������) 
	 */
	int writePLS(const char *fname);

	/* for test only */
	int savePLS_IntersectCheck(MYPOINT ln[2], MYPOINT fac[3], MYPOINT pnt, MYPOINT pntCpy);
	/*
	 * write .pl3 file
	 */
	int writePL3(const char *fname);
	/*
	 * write volume mesh information to parameters
	 */
    int writeOutputParameters(
	/* ---------------------------�����Ϣ---------------------------------*/
	int *ngp,               /* number of grid points */
	double  **g_coord,      /* gird point coords., 3 * ngp */
	int *nel,               /* number of elements */
	int **gtopu,             /* (j1, j2, j3, i4), 4 * nel */
	int *prt_arr          /* parent elements, size = nbp) */
    );
	/*
	 * write .ngb file
	*/
	int writeNGB(const char *fname);
	/*
	* write .ngb file
	*/
	int readNGB(const char *fname);

	/* 
	 * print the meshing result in the screen
	 */
	int printResult();

	/*
	 * read .pl3 file
	 */
	int readParPL3(const char *fname,
		INTEGER *nElem, INTEGER *nNode, INTEGER *nSurf);
	int readGriPL3(const char *fname,
		float *coords, INTEGER *elem, INTEGER *surf);
	void setGridData(const char *fname);

	/* 
	 * write .obj file
	 */
	int writeOBJ(const char *fname);
	
	/* space control */
	REAL spacFrmPnt(PntSource *src, MYPOINT pt);
	REAL spacFrmLne(PntSource *src1, PntSource *src2, MYPOINT pt);
	REAL spacFrmTri(PntSource *src1, PntSource *src2, PntSource *src3, MYPOINT pt);
	REAL fromSrc(MYPOINT pnt);


	int checkEmptyEles();
	/*
	 * Check the validity of the data after ia point insertion
	 * Parameters:
	 *	 void
	 * Result:
	 *   and error code
	 */
	int checkGlobalData(bool bInvHole = true);
	
	int checkNeigs(bool bFinal = false);

	int checkNodeData();

	/* --------------------------------
	 * ���ڼ�鵥Ԫ���������������Ƿ�
	 * ��ȷ����recvFace֮�󱻵���
	 * bInitFacet = true    ���α߽�ָ�����ã���������Ƭ
	 * bInitFacet = false   Լ���߽�ָ�����ã�����������Ƭ
	 * -------------------------------*/
	int checkOuIn(bool bInitFacet);

	/* ----------------------------------
	 * ����Ƿ�߽������ȫ�ָ�(���λָ�)
	 * --------------------------------*/
	int checkConfEdges();

	/* ----------------------------------
	 * ����Ƿ�߽������ȫ�ָ�(���λָ�)
	 * --------------------------------*/
	int checkConfFaces();

	/* Ѱ�Ұ���ĳ����Ƭ��������Ԫ */
	int findParents(INTEGER i1, INTEGER i2, INTEGER i3, INTEGER prts[], int *num);

	/* --------------------------------
	 * ���ڼ�鵥Ԫ���������������Ƿ�
	 * ��ȷ����recvFace֮�󱻵���
	 * -------------------------------*/
//	int checkInOut();

#ifdef _ERROR_CHK
	int writeOBJ(const char *fname,
			int cav_fc[], int nfc,
			MYPOINT cav_vt[], int npc);

	int chkEleResv();
	int chkNodResv(INTEGER end, INTEGER iNodE);

	int chkAllFacs(INTEGER iNodE);

#endif 

	/* 
	 * write a .pls file for the 3D shell orbit 
	 */
	int writePLS_ShellOrbit(char *fname, INTEGER shellNDs[], INTEGER nShellNDs);
	
	/*
	 * write ia .obj file after inserting boundary points
	 */
	int writeOBJAftInstBP(const char *fname);

#ifdef _TEST
	/*
	 * remove elements without recovery
	 */
	int rmvEleWoutRecv();

	/*
	 * write an .obj file of the surface mesh (without recovery)
	 */
	int writeOBJOfSurfWoutRecv(const char *fname);
	
	/*
	 * write an .obj file of the surface mesh after reading the PLS file
	 */
	int writeOBJAftReaPLS(const char *fname);

	/*
	 * write an .pip file while find a pipe when recover a edge
	 */
	int writePipe(Pipe pip, int nPip, const char *fname);

	/*
	 * write an object file after remove all outer elements
	 */
	int writeOBJAftRmvOuterEles(const char *fname);
	
	int writeOBJOfInner(const char *fname);
	int writeOBJOfOuter(const char *fname);
	int writeOBJOfUndef(const char *fname);

	int writePL3_Subset_Global(const char *fname, int elems[], int nElems);
#endif //_TEST

protected:
	/* *********************************************************************************
	 * �޸Ŀ�ǻ��ʹ����Ч(���͡����ӡ����ա�)
	 * �㷨��Դ��INRIA����Ľ�׳�ĵ�����ں�
	 *
	 * ��֪����: 
	 *     ����Ԫ��       �����������ĵ�Ԫ(��������һ������ʱ�������ñߵ�"��"�����е�Ԫ
	 *                    ��Ϊ����Ԫ
	 *     ��ǻ��         �����ƻ�Delaunay׼��ĵ�Ԫ
	 *     ��ǻ�߽磺	  ��ǻ�ı߽��棬��һ����Ԫ������ʾ(iElem, iCode)��iElem��ʾ���ڿ�
	 *                    ǻ�Ұ�����ǻ�߽�ĵ�Ԫ����Ʊ߽�Ԫ��iCode��ʾ��ǻ�߽��ڱ߽�Ԫ��
	 *                    �ı��(0/1/2/3)
	 *     
	 * ********************************************************************************/
	typedef struct CaviBnd
	{
		INTEGER	iElem;	/* ���ڿ�ǻ�Ұ�����ǻ�߽�ĵ�Ԫ����Ʊ߽�Ԫ */
		INTEGER iCode;	/* ��ǻ�߽��ڱ߽�Ԫ�еı��(0/1/2/3) */
		REAL   volume, cen[DIM], rad;	/* ���Ӵ������͸ñ߽�Ԫ�γɵ��µ�Ԫ����� */
		bool disable_param_calc;		/* �Ƿ���Ҫ���в������� */
	} CaviBnd;

	typedef struct CaviEdg
	{
		INTEGER iEnd;
		INTEGER hashNxt;
	} CaviEdg;

	void delCaviBond(INTEGER iCaviE, INTEGER iCode);
	void addCaviBond(INTEGER iCaviE, INTEGER iCode);

	int getCaviEdge(INTEGER iNod1, INTEGER iNod2);
	int addCaviEdge(INTEGER iNod1, INTEGER iNod2);

	inline bool isCaviBond(INTEGER iCaviE, INTEGER iCode)
	{
		assert(iCode >= 0 && iCode <= DIM);
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		return m_pElems[iCaviE].iReserved2 & (1 << (5 + iCode*2));
	}

	inline bool isCaviBondOverset(INTEGER iCaviE, INTEGER iCode)
	{
		assert(iCode >= 0 && iCode <= DIM);
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif

		return m_pElems[iCaviE].iReserved2 & (1 << (6 + iCode*2));
	}

	inline void flagCaviBond(INTEGER iCaviE, INTEGER iCode, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iCaviE].iReserved2 &= ~(1 << (5 + iCode*2))) :
					(m_pElems[iCaviE].iReserved2 |= (1 << (5 + iCode*2)));
	}

	inline void flagCaviBondOverset(INTEGER iCaviE, INTEGER iCode, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iCaviE].iReserved2 &= ~(1 << (6 + iCode*2))) :
					(m_pElems[iCaviE].iReserved2 |= (1 << (6 + iCode*2)));
	}

	inline int getCaviBondNodeRef(INTEGER iNod)
	{
	#if 0
		std::map<INTEGER, INTEGER>::iterator mapit = 
			m_mapCBNodeRef.find(iNod);

		return mapit == m_mapCBNodeRef.end() ? 0 : mapit->second;
	#endif
#ifdef _ONE_LONG_ARRAY
		assert(iNod >= 0 && iNod < m_nNodes);
#else
		assert(iNod >= 0 && iNod < m_pNodes.getArraySize());
#endif
		return m_arrCBNodeRef[iNod];
	}

	inline int minusCaviBondNodeRef(INTEGER iNod)
	{
	#if 0
		std::map<INTEGER, INTEGER>::iterator mapit = 
			m_mapCBNodeRef.find(iNod);
	
		if (mapit != m_mapCBNodeRef.end())
		{
			if (mapit->second > 0)
				mapit->second--;
			return mapit->second;
		}
	#endif
#ifdef _ONE_LONG_ARRAY
		assert(iNod >= 0 && iNod < m_nNodes);
#else
		assert(iNod >= 0 && iNod < m_pNodes.getArraySize());
#endif
		assert(m_arrCBNodeRef[iNod] > 0);
		return --m_arrCBNodeRef[iNod];
	}

	inline int plusCaviBondNodeRef(INTEGER iNod)
	{
	#if 0
		std::map<INTEGER, INTEGER>::iterator mapit = 
			m_mapCBNodeRef.find(iNod);
	
		if (mapit != m_mapCBNodeRef.end())
		{
			assert(mapit->second >= 0);
			mapit->second++;
			return mapit->second;
		}
		else
		{/* add one cavity boundary node */
			m_mapCBNodeRef.insert(make_pair(iNod,1));
			return 1;
		}
	#endif

#ifdef _ONE_LONG_ARRAY
		assert(iNod >= 0 && iNod < m_nNodes);
#else
		assert(iNod >= 0 && iNod < m_pNodes.getArraySize());
#endif
		return ++m_arrCBNodeRef[iNod];
	}

	void delCaviElem(INTEGER iCaviE);
	void addCaviElem(INTEGER iCaviE);
	void addBaseElem(INTEGER iBaseE);

	inline bool isCaviElem(INTEGER iCaviE) 
	{
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		return m_pElems[iCaviE].iReserved2 & 0x8;
	}
	inline bool isCaviElemOverset(INTEGER iCaviE)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		return m_pElems[iCaviE].iReserved2 & 0x10;
	}
	inline void flagCaviElem(INTEGER iCaviE, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iCaviE].iReserved2 &=~(0x8)) :
					(m_pElems[iCaviE].iReserved2 |= (0x8));
	}
	inline void flagCaviElemOverset(INTEGER iCaviE, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iCaviE >= 0 && iCaviE < m_nElems);
#else
		assert(iCaviE >= 0 && iCaviE < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iCaviE].iReserved2 &=~(0x10)) :
					(m_pElems[iCaviE].iReserved2 |= (0x10));
	} 
	inline bool isBaseElem(INTEGER iBaseE)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iBaseE >= 0 && iBaseE < m_nElems);
#else
		assert(iBaseE >= 0 && iBaseE < m_pElems.getArraySize());
#endif
		return m_pElems[iBaseE].iReserved2 & 0x2;
	}

	inline bool isBaseElemOverset(INTEGER iBaseE)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iBaseE >= 0 && iBaseE < m_nElems);
#else
		assert(iBaseE >= 0 && iBaseE < m_pElems.getArraySize());
#endif
		return m_pElems[iBaseE].iReserved2 & 0x4;
	}
	inline void flagBaseElem(INTEGER iBaseE, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iBaseE >= 0 && iBaseE < m_nElems);
#else
		assert(iBaseE >= 0 && iBaseE < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iBaseE].iReserved2 &=~(0x2)) :
					(m_pElems[iBaseE].iReserved2 |= (0x2));
	}
	inline void flagBaseElemOverset(INTEGER iBaseE, INTEGER flag)
	{
#ifdef _ONE_LONG_ARRAY
		assert(iBaseE >= 0 && iBaseE < m_nElems);
#else
		assert(iBaseE >= 0 && iBaseE < m_pElems.getArraySize());
#endif
		flag == 0 ? (m_pElems[iBaseE].iReserved2 &=~(0x4)) :
					(m_pElems[iBaseE].iReserved2 |= (0x4));
	}
	inline int getCaviElemNodeRef(INTEGER iNod)
	{
	#if 0
		std::map<INTEGER, INTEGER>::iterator mapit = 
			m_mapCENodeRef.find(iNod);

		return mapit == m_mapCENodeRef.end() ? 0 : mapit->second;
	#endif
#ifdef _ONE_LONG_ARRAY
		assert(iNod >= 0 && iNod < m_nNodes);
#else
		assert(iNod >= 0 && iNod < m_pNodes.getArraySize());
#endif
		return m_arrCENodeRef[iNod];
	}

	inline int minusCaviElemNodeRef(INTEGER iNod)
	{
	#if 0
		std::map<INTEGER, INTEGER>::iterator mapit = 
			m_mapCENodeRef.find(iNod);
	
		if (mapit != m_mapCENodeRef.end())
		{
			if (mapit->second > 0)
				mapit->second--;
			return mapit->second;
		}
	#endif

#ifdef _ONE_LONG_ARRAY
		assert(iNod >= 0 && iNod < m_nNodes);
#else
		assert(iNod >= 0 && iNod < m_pNodes.getArraySize());
#endif
		assert(m_arrCENodeRef[iNod] > 0);
		return --m_arrCENodeRef[iNod];
	}

	inline int plusCaviElemNodeRef(INTEGER iNod)
	{
	#if 0
		std::map<INTEGER, INTEGER>::iterator mapit = 
			m_mapCENodeRef.find(iNod);
	
		if (mapit != m_mapCENodeRef.end())
		{
			assert(mapit->second >= 0);
			mapit->second++;
			return mapit->second;
		}
		else
		{/* add one cavity boundary node */
			m_mapCENodeRef.insert(make_pair(iNod,1));
			return 1;
		}
	#endif
#ifdef _ONE_LONG_ARRAY
		assert(iNod >= 0 && iNod < m_nNodes);
#else
		assert(iNod >= 0 && iNod < m_pNodes.getArraySize());
#endif
		return ++m_arrCENodeRef[iNod];
	}

	int formInitCaviData(INTEGER baseE[], int baseSize, INTEGER iNod, bool instInner = true);
	int modifyCavityToValid(INTEGER iInstNode);
	int modifyCavityToNoIsolatedPnt();
	bool rejectFieldPoint(INTEGER iNod);
	int updateTriangulation(INTEGER iNod);
	int addInnerNode_Robust(INTEGER iNod, INTEGER iBaseE);
	int resetCavityNodeRef(INTEGER maxSize);
	int resetCavityEdgHash(INTEGER maxSize);
	int updateCavityBound();
	int updateCavityNodeRef();
	int updateCavityEdge();
	void clearCavityData();
	void updateDeletElems();

	int cleanCaviElems();
	int cleanCaviBonds();
	int checkCaviData();

	/* ���ڽ�һЩ��Ƭɾ���ĺ��� */
	int deleteUnusedSurfTris(double xmin, double xmax, 
		double ymin, double ymax, 
		double zmin, double zmax);

	int deleteUnusedSurfNodes();

	int sphereToPolyhedron(Sphere sph, int nSph, polyhedron *poly,
						   INTEGER *numv, APOINT_ARRAY vtxArray, 
						   INTEGER l2g[], std::map<int, int>& mapG2L);

	int sphereToPolyhedron(Sphere sph, int nSph, polyhedron *poly,
						   INTEGER *numv, APOINT_ARRAY vtxArray, 
						   INTEGER l2g[], IntIntMap& mapG2L);

	/* ------------------------------------------------------------------------------
	 * �������߽���ȡ����
	 * ----------------------------------------------------------------------------*/
	int sphereToPolyhedron(INTEGER iNod, Sphere sph, int nSph, polyhedron *poly,
						INTEGER *numv, APOINT_ARRAY vtxArray, 
						INTEGER l2g[], std::map<int, int>& mapG2L);
	/* �������� */
	int calc_sphere_quality(Sphere sph, int nSph, 
		std::map<int, int>& mapG2L, double *qual, INTEGER l2g_nod[], INTEGER &numv);
	int calc_sphere_quality(Sphere sph, int nSph, 
		IntIntMap &mapG2L, double *qual, INTEGER l2g_nod[], INTEGER &numv, bool bCalQuality);
	/* ���ֽ������뵽��ǰ�����廯�� */
	int addDivideResult2(divide *divi, INTEGER l2g_nod[], INTEGER numv, 
		Sphere sph, INTEGER nSph, INTEGER iNod);
	int addDivideResult(divide *divi, INTEGER l2g_nod[], INTEGER numv,
		Sphere sph, INTEGER nSph, int inOut = INNER, bool isSubT = false, INTEGER *newElems = NULL);
	/* ------------------------------------------------------------------------------
 	 * ���ĳ����ǻ��֤�����޷������廯������ǡ�������͸ÿ�ǻ
	 * ----------------------------------------------------------------------------*/
	 int flatSphereByFaceNeig(Sphere sph, int& nSph, int inOu = INNER, bool noSutTets = false);
	 int flatSphereByFaceNeig_SubTetOnly(Sphere sph, int& nSph);
	 int flatSphereByFaceNeig_NoIsolatedPnt(Sphere sph, int& nSph, bool noSutTets = false); /* �������ǻ���ڲ��� */
	 int flatSphereByEdgeNeig_NoIsolatedPnt(Sphere sph, int& nSph); /* �������ǻ���ڲ��� */
	 /* -------------------------------------------------------------------------------
	  * ����һ������չ��Ƭ�б�Կ�ǻ������չ 
	  * ----------------------------------------------------------------------------*/
	 int flatSphereByCandidateFaceList_NoIsolatedPnt(Sphere sph, int &nSph, InflatationEntList &listInflationEnts, 
		 INTEGER consEdges[][2], int consEdgeSize);
	 bool isFaceNeedInflated(INTEGER indices[], INTEGER facePrt1, INTEGER facePrt2, int infDir[], 
		 INTEGER flatEdges[], int *flatEdgeSize, INTEGER consEdges[][2], int consEdgeSize);

	 int initCavityData(Sphere sph, int nSph, int nBase); /* ��ʼ��һ����ǻ�����޸�����ʹ���������µ� */

	 /* ------------------------------------------------------------------------------
	  * ��������SPR(Small Polyhedron Reconnection)����ɾ��һ���ڲ��㣻
	  * ----------------------------------------------------------------------------*/
	int removeInnnerNode(INTEGER iNod, VECTOR norm, Sphere sph, int nSph);

	/* ------------------------------------------------------------------------------
	 * ������������Steiner���λ��
	 * ----------------------------------------------------------------------------*/
	int repositionInnnerNode(INTEGER iNod, VECTOR norm, Sphere sph, int nSph);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ��һ���ڲ��㣻
	 * �������·���1����:
	 * �ȳ���ɾ��/�ٳ������� 
	 * ----------------------------------------------------------------------------*/
	int removeThenRepositionInnnerNode(INTEGER iNod, VECTOR norm, Sphere sph, int nSph, bool isNeedReposition = false);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ��һ���ڲ��㣻
	 * �������·���1����
	 * ----------------------------------------------------------------------------*/
	int removeAndRepositionInnnerNode(INTEGER iNod, VECTOR norm, Sphere sph, int nSph, bool isNeedReposition = false);

	/* 
 	 * ������һ���������з���һ���㣬ʹ��Ͷ���������б���Ƭ���������
	 */
	int calcSphereCenter(INTEGER sph_cen, VECTOR norm, MYPOINT p_init, 
					REAL *deta_init, polyhedron *poly);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ��һ���ڲ���
	 * ----------------------------------------------------------------------------*/
	int removeInnnerNode_SPR(INTEGER iNod);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ��һ���ڲ���
	 * ----------------------------------------------------------------------------*/
	int removeInnnerNode_SPR_QualImprv(INTEGER iNod, bool bSmooth = false, INTEGER iSmoothedNd = -1);

	/* ------------------------------------------------------------------------------
	 * ��������Flip����ɾ��һ���ڲ���
	 * ----------------------------------------------------------------------------*/
	int removeInnnerNode_Flip(INTEGER iNod);

	int getDesEdges(DesFacet *pDF, DesEdge *pDSArray[]);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ��һ���ڲ���
	 * ----------------------------------------------------------------------------*/
	int removeInnnerNodeArray_SPR(INTEGER nodeArray[], INTEGER nodeSize);

public:
	int findSegEndInCaviB(INTEGER caviBNs[], int nCaviBNs, INTEGER steinerArr[], int nSteiner, INTEGER segEnds[], int *nSegEnds);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ�����ϵ�һ��Steiner��
	 * segBeg<-->segEnd: Ҫ�����ı�
	 * ----------------------------------------------------------------------------*/
	int removeEdgeSteinerPoint_SPR(INTEGER iNod, INTEGER segBeg, INTEGER segEnd, INTEGER steinerArr[], int nSteiner);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ�����ϵ�һ��Steiner��
	 * iFace[]: Ҫ��������
	 * ----------------------------------------------------------------------------*/
	int removeFaceSteinerPoint_SPR(INTEGER iNod, INTEGER fIdx, INTEGER iface[], 
		INTEGER nStnBeg, INTEGER nStnEnd);

	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ�������ڲ�Steiner��
	 * ----------------------------------------------------------------------------*/
	int removeInnerSteinerPoint_SPR();

	/* ------------------------------------------------------------------------------
	 * ��������flip����ɾ�������ڲ�Steiner��
	 * ----------------------------------------------------------------------------*/
	int removeInnerSteinerPoint_Flip();
	
	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ�����б��ϵ��������ڲ�Steiner��
	 * (һ����һ���ߵķ�ʽ)
	 * ----------------------------------------------------------------------------*/
	int removeEdgeSteinerPoint_SPR(int iDesEdg);
	int removeEdgeSteinerPoint_SPR(bool bNotTryOne = true); /* ������ֻ��һ���ڵ�ʱ�������� */
	/* ------------------------------------------------------------------------------
	 * ��������SPR(Small Polyhedron Reconnection)����ɾ���������ϵ��������ڲ�Steiner��
	 * (һ����һ���ߵķ�ʽ)
	 * ----------------------------------------------------------------------------*/
	int removeFacSteinerPoint_SPR(bool bTryInrOrTwoEdg = true); /* ������2��������Steiner�㣬�����2��������Steiner��ʱ */
	
	int removeAllLeftNodes(); 

	int printSteinerPntInfo();

	/* --------------------------------------------------------------------------------------
	 * �ڲ����ڲ����ɾ��Steiner�㣬��Ϊ�ڲ����ڲ���Ĺ����У��ڵ�͵�Ԫ��״̬��Ҫ���±궨
	 * ���iReserved�򣺲������õĵ㣬��ֵΪ-1����ʾ��ɾ�����������õĵ㣬��ֵָ��������ĳ����Ԫ
	 * ���iReserved�����±궨��Ԫ��in/ou���ԣ�����δ��ɾ���ĵ�Ԫ�궨ΪIN
	 * ----------------------------------------------------------------------------------------*/
	int prepareRmvSteinerPntAftInstInnNode();

protected:
	int cleanSteinerPointVector();

	std::vector<INTEGER> m_vecSteinerPnt; /* Steiner�� */

	//�����Ż�
public:
	void getactiveset(Sphere sph,
                  int nSph,
				  REAL quals[],
				  REAL qualgrads[][3],
                  REAL activegrads[][3],
                  int *numactive,
                  REAL worstqual,
				  int qualmeasure);
	REAL getinitialalpha(INTEGER iNod,
                     Sphere sph,
					 int nSph,
					 REAL quals[],
				     REAL qualgrads[][3],
                     REAL d[3],
                     REAL r,
                     REAL worstqual,
					 int qualmeasure);
	void nonsmoothlinesearch(INTEGER iNod,
	                     Sphere sph,
						 int nSph,
                         REAL d[],
                         REAL inworstqual,
						 REAL *ouworstqual,
                         REAL *alpha,
                         REAL r, 
						 int qualmeasure);
	void getoptinfo(INTEGER iNod,
                 INTEGER iElem,
                 REAL *qual,
                 REAL qualgrad[][3],
                 REAL *volume,
                 REAL volumegrad[3],
                 int qualmeasure);
	bool combinedSmoothing(INTEGER iNod,
	           Sphere sph,
			   int nSph,
			   bool bSmartLap,
			   REAL inworstqual,
			   REAL qualthreshold,
			   REAL *outworstqual,
			   REAL *outquals,
			   int qualmeasure);
	void initoptimize(int qualmeasure);
	void meshopt();
	bool meshcomp(int &icnt, float &fmin, float &faver, bool bfirst = false);
	void evalMesh(int *numOfBadElems, REAL *minQual, REAL *aveQualOfBadElems, int qualmeasure, REAL qualLimit);

	// smooth����
	void setsmoothpar(char *filename);
	int smooth(int imaxpasses, bool noQualLimit, REAL qualLimit, int qualmeasure);
	int neighborhood(INTEGER node, Sphere sph, int nSph, int *num_incident_vtx, int **neigVerts, int *num_tet,
		double** free_vtx, double*** vtx_list, int*** vtx_connectivity, int*** cell_vtx, int** cells);
	bool isMeshBetter(int ncell, int* cells, double *newQuals, double angDegree);
	void updateMeshQual(int ncell, int* cells, double *newQuals);
	int lookupidx(int *pNodes, int npt, int idx);
	void setNodeIllShaped(int ptIdx);
	void setNodeWellShaped(int ptIdx);
	bool isNodeWellShaped(int ptIdx);

	void enableNodeSmoothed(INTEGER iNod);
	void disableNodeSmoothed(INTEGER iNod);
	bool isNodeSmoothed(INTEGER iNod);
	
	void enableSuppressionTried(INTEGER iNod);
	void disableSuppressionTried(INTEGER iNod);
	bool isSuppressionTried(INTEGER iNod);

#if 1
	int countEC(INTEGER iNod);
	int countECIncre(INTEGER iNod);
	int countECDecre(INTEGER iNod);
	int analyseECStatistics();
#endif
	
	int multifaceRemovalPass(int nMaxPass, REAL qualLimit = 0.5, int qualmeasure = QUALWARPEDMINSINE);
	int multifaceRemoval(INTEGER facep[3], INTEGER a, INTEGER b, INTEGER iElem, INTEGER iNeig, 
		REAL qualLimit = 0.5, int qualmeasure = QUALWARPEDMINSINE);

	//���˱任���Ǳ任��
	void prepareElemInfo(int qualmeasure);
	void updateElemInfo();
	int shelltransform(int nMaxPass, int maxResurLevel, REAL qualLimit = 0.5, int qualmeasure = QUALWARPEDMINSINE);	//������ֵ����threshold������Ԫ���пǱ任
	int tranverseST(int iPass = 1);	//����������Ԫ���пǱ任

	int shelltransform(INTEGER a, INTEGER b, INTEGER c, Shell she, int *pnShe, 
		int* cstNodes, InstEntList& listDependentEnts, int level,  REAL qualLimit = 0.5, int qualmeasure = QUALWARPEDMINSINE);
	int shelltransform_multiface_removal(INTEGER a, INTEGER b, INTEGER skirtNodes[], INTEGER faceElems[][2], int skirtSize, REAL qold, REAL qualLimit = 0.5, int qualmeasure = QUALWARPEDMINSINE);
	/* �滻������Ϊ������ */
	int replaceLocalMesh(INTEGER oldElemsArray[], int oldElemSize, INTEGER newElemsConns[], int newElemSize, INTEGER *newElemsArray);

	int recursiveST(INTEGER a, INTEGER b, INTEGER c, int iElm, Shell she, int *pnShe, 
		int* cstNodes, InstEntList& listDependentEnts, int level, 
		REAL qualLimit = 0.5, int qualmeasure = QUALWARPEDMINSINE);
	int removeedge(int iElm, int a, int b);
	int swap23Recursive(
		INTEGER a, INTEGER b, INTEGER matrix_k[][MAX_SKIRT_POLY_SIZE], 
		float matrix_b[][MAX_SKIRT_POLY_SIZE], float *matrix_q_A, float *matrix_q_B, 
		float qualLimit, int qualmeasure,
		int ii, int jj, int faceNum, INTEGER faceNodes[], INTEGER faceParents[], InstEntList& listDependentEnts);

	int removeInnnerNode_opt(INTEGER iNod);
	int removeInnerNode_OptSPR(int iNod, bool bSmooth, int iSmoothNd);
	int removeInnerNode_OptSPR(int iNod, REAL qualLimit, int qualmeasure);
	int splitEdge_Opt(INTEGER s, INTEGER e, INTEGER iElm, REAL qualLimit, int qualmeasure, int smoothOpt = 1);

	int collapseedge(int s, int e, int iele);
	bool splitqoption(float oqual, float nqual);
	int splitedge(int s, int e, int iele);
	void setShellPipe(int s, int e, Shell she, int nshe, int p, Pipe pip, int &nPip);
	void nodeedge(int iNod, int edgeNodes[], int &nedge);

	void flagdepfacenode(int a, int b, int *nodeFlag, int *faceNodes, int nface, InstEntList& listDependentEnts);
	//����������������������
	float tetqual(INTEGER a, INTEGER b, INTEGER c, INTEGER d, float* quals, bool &bAcute, QUALTYPE qtype = QUALTYPE::DIHEDRAL_ANGLE_SINE);
	//boutputΪtureʱ������������������diahedral��
	float tetdihedral(INTEGER a, INTEGER b, INTEGER c, INTEGER d, bool &bAcute, bool boutput = false, float* dihedral = NULL, bool bangle = false);
	float dihedral(int npt, float coord[][3], bool &bAcute, bool bangle = false);
	float dihedral(int a, int b, int c, int d, bool &bAcute);
	void updatequal(INTEGER a, INTEGER b, INTEGER c, INTEGER d, INTEGER e, int tElems[2], 
		float threshold = MESH_QUAL_THRESHOLD);
	void setoptangle(int iElm);
#ifdef _DEBUG
	void setoptangle(int iElm, REAL qual, REAL qualLimit, int qualmeasure, bool bCheck = true);
#else
	void setoptangle(int iElm, REAL qual, REAL qualLimit, int qualmeasure);
#endif
	void setoptangle(int a, int b, int c, int d, REAL qualLimit, int qualmeasure);
	int getfaceidx(int a, int b, int *faceNodes, int faceNum, int *nxtFace);
	int getfaceidx(int a, int b, Shell she, int nShe, int *faceNodes, int faceNum, int *nxtFace);
	int sortShellFaces(int a, int b, Shell she, int nShe, int *faceNodes, int *faceParents, int faceNum, int *sortedNodes);
	int getfaceidx(int fn[2], int faceNum, int *nxtFace);

	void printshell(int a, int b, int *faceNodes, int faceNum);
	void printtet(int ielm);
	void printBadCells(const char *chFileName, int qualmeasure, float qualLimit);

	void getshellqual(int a, int b, Shell she, int nShe, float *qual);
	void getshellqual(Shell she, int nShe, float qualArray[]);
	float qualoption(float qold, int level);
	float qualrefer(int level, float qold);
	void cyclingShellQual(int a, int b, int faceNum, int* faceNodes, int (*matrix_k)[128], 
		int* ijCycling, int ijCyclingCnt, std::set<float>& cylingQuals, int &cyclingTetNum);
	void qualRecursive(int a, int b, int (*matrix_k)[128], int i, int j, int faceNum, int *faceNodes, 
		std::set<float>& cylingQuals, int &cyclingTetNum);
	void cyclingShellQual(int a, int b, int faceNum, int* faceNodes, int (*matrix_k)[128], 
		int* ijCycling, int ijCyclingCnt, REAL& minQual, int &cyclingTetNum, int qualmeasure);
	void qualRecursive(int a, int b, int (*matrix_k)[128], int i, int j, int faceNum, int *faceNodes, 
		REAL& minQual, int &cyclingTetNum, int qualmeasure);
	void fillInShellElems(int a, int b, int skirtSize, int *skirtNodes, int (*matrix_k)[128], 
		int *ijCycled, int ijCycledCnt, int newElemsConns[], int &newElemsSize);
	void fillInShellElemsRecursive(int a, int b, int (*matrix_k)[128], int i, int j, int skirtSize, int *skirtNodes, 
		int *ijCycled, int ijCycledCnt, int newElemsConns[], int &newElemsSize);

	bool qualCompare(int cnt, float *quals, std::set<float>cylingQuals);
	void copyQuals(std::set<float>cylingQuals, int &cnt, float *quals);

	void attackbadcells_V1();
	void contractEdges(REAL qualLimit, int qualmeasure);
	void splitlong_V1();
	void splitlong(REAL qualLimit, int qualmeasure);
	void checkbadcells();
	void printElemQuality(const char *chFileName, int qualmeasure, float qualLimit, bool bStored = false);
	void evaluateQuality();
	bool checkElemQuality(int qualmeasure);

	//new datastructure
	int addElmAngle(TetraElemQual tetEleQual, int elmIdx);
	int updateElmAngle(TetraElemQual tetEleQual, int elmIdx);
	int codeAngle(int elmIdx, int angIdx);
	//ɾ����ԪiElm�ڶ��еĽ�
	void removeElmAng(int iElm);
	void setAnglePolicy(float threshold);
	float getAngleUpper(float val);

	//��任
	void faceswap(int iPass);
	bool isSwaped(int iEle, int iFace);
	void setSwaped(int iEle, int iFace, bool bSwaped);

	void updatebdnpnts();


public:
	/* ----------------------------------------------------- 
	 * functions for library
	 * ----------------------------------------------------*/
	/* -------------------------------------------
	 * ����Source
	 * ------------------------------------------*/
	int crtSources(
		double		pdCX[],		/* x coord. of center points of sources */
		double		pdCY[],		/* y coord. of center points of sources */
		double		pdCZ[],		/* z coord. of center points of sources */
		double		pdIn[],		/* inner radius of sources */
		double		pdOu[],		/* out radius of sources */
		double		pdSp[],		/* space values of sources */
		int			nPS,		/* number of point source */
		int			nLS,		/* number of line source */
		int			nTS 		/* number of triangle source */);

	/* -------------------------------------------
	 * ������������
	 * ------------------------------------------*/
	int crtBkgMsh(	
		double		pdBMNX[],	/* x coord. of background mesh */
		double		pdBMNY[],	/* y coord. of background mesh */
		double		pdBMNZ[],	/* z coord. of background mesh */
		double		pdBMNSpc[],	/* space values of background mesh nodes */
		int			nBMN,		/* number of background mesh nodes */
		int			pnBMEFm[],	/* forming points of background mesh elements */
		int			pnBMENg[],	/* neighboring eles. of background mesh elements */
		int			nBME		/* number of background mesh elements */);

	/* -------------------------------------------
	 * ��������
	 * ------------------------------------------*/
	int crtSurface(
		double		pdSNX[],	/* x coord. of boundary nodes */
		double		pdSNY[],	/* y coord. of boundary nodes */
		double		pdSNZ[],	/* z coord. of boundary nodes */
		int			nSN,		/* number of boundary nodes */
		int			pnSFFm[],	/* forming points of surface facets */
		int			pnSFPt[],	/* patches of surface facets */
		int         nSF 		/* number of facets */);

	/* -------------------------------------------
	 * ���������
	 * ------------------------------------------*/
	int outVolMsh(
		double	  **ppdMNX,		/* x coord. of 3D mesh nodes */
		double    **ppdMNY,		/* y coord. of 3D mesh nodes */
		double    **ppdMNZ,		/* z coord. of 3D mesh nodes */
		int        *pnMN,		/* number of 3D mesh nodes */
		int       **ppnMEFm,	/* forming points of 3D mesh elements */
		int       **ppnMENg,	/* neighboring eles. of 3D mesh elements */
		int        *pnME,		/* number of 3D mesh elements */
		int       **ppnPrt		/* parents of surface facets */);

protected:
	std::vector<INTEGER> m_vecBaseElems;			/* ����Ԫ���� */
	std::vector<INTEGER> m_vecCaviElems;			/* ��ǻ��Ԫ���� */	
	std::vector<CaviBnd> m_vecCaviBonds;			/* ��ǻ�߽����� */
	std::vector<CaviEdg> m_vecCaviEdges;			/* ��ǻ�߽����Ϣ(�ڻָ��߽����ʹ��) */
	//std::map<INTEGER, INTEGER> m_mapCENodeRef;	/* ��ǻ��Ԫ�Խڵ������ */
	//std::map<INTEGER, INTEGER> m_mapCBNodeRef;	/* ��ǻ�߽�Խڵ������ */
	short int *m_arrCENodeRef;						/* ��ǻ��Ԫ�Խڵ������ */
	short int *m_arrCBNodeRef;						/* ��ǻ�߽�Խڵ������ */
	int *m_pCaviEdgeHash;							/* ��ǻ�߽�Ĺ�ϣ�� */
	std::vector<INTEGER> m_vecNodeRef;				/* ���б����õĽڵ� */
	std::vector<INTEGER> m_vecCaviENodeRef;		    /* ��ǻ�߽��ϣ���Ӧ��NodeRef��ֵ */ 
	int m_arrCENodeRefSize;	
	int m_arrCBNodeRefSize;
	int m_nCaviEdgeHashSize;					
	INTEGER m_nDistinctC, m_nDistinctB;			   /* ��ǻ��Ԫ�е����Ŀ/��ǻ�߽��е����Ŀ */

protected:
	/* *********************************************************************************
	 * ʵ�ֽ�׳�ĵ�����ں�����ĸ������ݽṹ
	 * ��֪����: 
	 *     ����Ԫ��       �����������ĵ�Ԫ(��������һ������ʱ�������ñߵ�"��"�����е�Ԫ
	 *                    ��Ϊ����Ԫ
	 *     ��ǻ��         �����ƻ�Delaunay׼��ĵ�Ԫ
	 *     ��ǻ�߽磺	  ��ǻ�ı߽��棬��һ����Ԫ������ʾ(iElem, iCode)��iElem��ʾ���ڿ�
	 *                    ǻ�Ұ�����ǻ�߽�ĵ�Ԫ����Ʊ߽�Ԫ��iCode��ʾ��ǻ�߽��ڱ߽�Ԫ��
	 *                    �ı��(0/1/2/3)
	 *     
	 * ****************************************************************************** */
	/*basic data structures*/
#ifdef _ONE_LONG_ARRAY
	Elem *m_pElems;			/* elements array */
	INTEGER m_nElems;       /* the number of elements */
	INTEGER m_nAllocElems;		/*allocated data size for elements*/

	Node *m_pNodes;         /* nodes array */
	INTEGER m_nNodes;       /* the number of nodes */
	INTEGER m_nAllocNodes;		/*allocated data size for nodes*/
#else
	BlockedArray_Elem m_pElems;
	BlockedArray_Node m_pNodes;
#endif

	SurTri *m_pSurTris;		/* surface triangles*/
	SurEdg *m_pSurEdgs;		/* surface edges */
	
	INTEGER m_nSurNodes;	/* the number of surface nodes */
	INTEGER m_nSurTris;		/* the number of surface triangles */
	INTEGER m_nSurEdgs;		/* the number of surface edges */
	
	/* for dynamical memory allocation */
	INTEGER m_nAllocSurTris;	/*allocated data size for surface triangles*/
	INTEGER m_nAllocSurEdgs;	/*allocated data size for surface edges*/
	
	INTEGER m_nPeakElems;		/* the maximal size of the occupied space of the element array */
	INTEGER m_nPeakNodes;		/* the maximal size of the occupied space of the node array */

	/*
	 * ΪԼ���߽�ָ��㷨��������ݳ�Ա
     * data membersfor Constrained Boundary Recovery 
     */
	DesEdge *m_pDesEdges;
	INTEGER m_nDesEdges;
	INTEGER m_nAllocDesEdges;

	DesFacet *m_pDesFacets;
	INTEGER m_nDesFacets;
	INTEGER m_nAllocDesFacets;

	/* Ϊԭʼ�ı߽�ߴ�����ϣ�� */
	INTEGER *m_pSurEdgHash;
	INTEGER  m_nSurEdgeHash_MaxIdx;
	/* scale factor & scale center*/
	REAL m_scale;						/*scale factor*/
	REAL m_scale_recip;
	MYPOINT m_cenW;	          /*world coordinates of center points*/
	MYPOINT m_cenN; 					/*normalization coordinates of center points*/
	MYPOINT m_minW, m_maxW;   /*world corrdinates range*/
	MYPOINT m_minN, m_maxN;   /*normalization coordinates range*/
	/* 
	 * ���嵥Ԫ�ܶȿ��Ƶ����ַ�ʽ���������� & Դ 
	 * two kinds of space controlling is defined: 
	 * background mesh & sources
	 */
	BGMesh m_BGMesh;
	Source m_Source;

	Node *m_pCreateNodes;	/* �Զ��������ڲ��㣬��ʱ���ڴ棬������1��ʼ */
	INTEGER m_nCreateNodes;	/* �Զ��������ڲ�����Ŀ */

#ifdef _USING_STD_LIB
	/*a stack of elements for tree search*/
	std::stack<INTEGER> m_stkTreeSear;
#else
	IntStack m_stkTreeSear;
#endif
			
	/*a vector of deleted elements*/
	std::vector<INTEGER> m_vecDelEles;
		
	/*a vector of empty locations*/
	std::vector<INTEGER> m_vecEmpLocs;
	int m_nLocInd; /*location indicator*/
				
	/*a vector of added elements*/
	std::vector<INTEGER> m_vecNewEles;
	int m_nAftTail; //new elements after tail

#ifdef _USING_STD_LIB
	/*a vector of tested elements*/
	std::vector<INTEGER> m_vecTstEles;
#else
	MyVector<INTEGER> m_vecTstEles; 
#endif

	/*a vector of checked elements (used in the stage of face recovery*/
	std::vector<INTEGER> m_vecFREles;

	/*ia vector of cavity nodes*/
	std::vector<INTEGER> m_vecCavNodes;

	/*ia vector of cavity edges found*/
	std::vector<CavSide> m_vecCavSides;

	/*ia list of boundary node list (sorted by the inserting order)*/
	std::list<INTEGER> m_lstInstBndNods;

	/* ia list of disturbance information */
	std::list<DistInfo*> m_lstDistInfo;

	/* Boundary steiner points begin index & end index */
	INTEGER m_nBndSteinerBeg, m_nBndSteinerEnd;

	REAL m_fAlpha;	/* ȫ������ */

	std::vector<INTEGER> m_vecNewCreEles; /* ֻ����Щ�´����ĵ�Ԫ�ϼӵ� */

	bool m_bFacRigInner;

	/* ��ʱ��¼subTet */
	std::vector<INTEGER> m_vecSubTets; /* ��¼��/��ָ�ʱ�ֲ���ǻ�ڲ��ĵ�Ԫ */

	/* ��ʱ��¼�������Ԫ */
	std::vector<INTEGER> m_vecNegativeVElems;

protected:
	/* --------------------------------------------------------
	 * ���º������ڹ����ǻ�к���Ƭ���ཻ��ʵ�����Ŀ 
	 * -------------------------------------------------------*/
	int initPolyIntEntities(INTEGER s, INTEGER e, Pipe pips, int nPip);
	int updatePolyIntEntities(INTEGER s, INTEGER e, INTEGER flagElems[], int elemSize);
	
	/* ----------------------------------------------------
	 * �������ݳ�Ա��¼һ����ǻ�к���Ƭ���ཻ��ʵ�����Ŀ 
	 * -------------------------------------------------*/
	IntEHash m_arrPolyIntEdges[MAX_POLY_INT_EDGE];
	IntFHash m_arrPolyIntFaces[MAX_POLY_INT_FACE];
	INTEGER *m_arrPolyIntEHead;
	INTEGER *m_arrPolyIntFHead;
	INTEGER m_nPolyIntEdges, m_nPolyIntFaces;
	INTEGER m_nPolyIntEHead, m_nPolyIntFHead;
	bool m_bBndProtection;	/* �ڻָ��߽�ʱ�Ƿ�����ѻָ��߽磬��ʼֵΪfalse */
	int m_nMaxInflationLevel;
	int m_nLostBndEdges;
	int m_nLostBndFaces;
	int m_nCutBndFaces; /* �����м�̴�����Ƭ�� */
	InstEntList m_listBndIntEnts;
	bool m_bRecordInstEnts;
	bool m_bHillClimbingSPs;
	bool m_bLessInstOneBE;

	/* -----------------------------------------------------------------------
	 * ����2���������topology improvement for boundary recovery��Ч��
	 * ---------------------------------------------------------------------*/
	int m_nShellTransCallingTimes;
	int m_nShellTransDidSomeTimes;
	int m_nShellTransSuccessTimes;
	int m_nMultiRemovComplexTimes;
	int m_nMultiRemovShTransTimes;

	bool m_bRecordNewEles; /* ��¼�´����ĵ�Ԫ */
	bool m_bRecurFlip_FaceLock;	/* �ڵݹ鷭תʱ��������Ƭ(����Ч������)����������Ԫ���죬Ч���Բ */
	bool m_bBndESheTransEnabled;

	bool m_bWithOutBoxNodes;


	//smooth
	void *pvSmoothData;

	//shell transformation
	BlockedArray_TetQuals m_tetraElmQuals;
	ActElemHeap* m_elmHeap;
	HashPolyList* m_hashpoly;
	SkirtPolyHeap* m_heappoly;
	SkirtPolyHeap* m_heappoly2;
	bool m_bopt;
	bool m_blower;
	bool m_btranverse;
	bool m_brepair;
	int iloop;
	int m_n32swapTimes;
	int m_n23swapTimes;
	float threshold;
	float curqual;
	int curSwapElm;
	int m_qelmcnt[5]; //С��30�������ǵĸ���
	bool m_biter;
	QI_ST_TARGET m_eSTTarget;

	//�Ƕ�����
	int m_angInteVal;
	int m_nAngInter;
	float* m_angValues;
	int *m_eleCnt;	//����Ӧ����ǵĸ���

	/* �����Ż�ʱ����Ľڵ���ʼ��� */
	INTEGER m_nEdgeSplitNodeStartID;
	INTEGER m_nInitMeshNodes;
	INTEGER m_nCurMeshNodes;
	REAL m_dNodeDiffThreshold;

	/* �����Ż��������� */
	int m_eQualityMeasure;	/* �Ż����� */
	REAL m_dQualThreshold;	/* ��ֵ */
	int m_nImprovementCycle;

	bool m_bEnableEdgeRemoval;
	bool m_bEnableFaceRemoval;
	bool m_bEnableEdgeContract;
	bool m_bEnableEdgeSplitting;
};
/* ע��qualArray1��qualArray2�е�ֵ�����������е� */
extern int compareQualArray(float qualArray1[], float qualArray2[], int size1, int size2, float epsilon);
#endif /* __iso3d_h__ */
