/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#ifndef __iso3d_define_h__
#define __iso3d_define_h__

//前三个点的右手规则指向第四个点,中间两个点为二面角边
const int tetra_pnt[6][4] = 
{ {0,1,2,3}, {0,2,3,1}, {0,3,1,2}, {1,2,0,3}, {1,0,3,2}, {2,0,1,3} };

const int dihed_edge[6] =
{ 5, 3, 4, 0, 2, 1 };

//前两个点表示一条边，前三个点的右手规则指向第四个点
const int tetra_edge[6][3] = 
{ {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 2, 0}, {1, 3, 2}, {2, 3, 0} };


//前三个点的右手规则指向第四个点
const int tetra_face[4][4] = 
{ {0,1,2,3}, {0,2,3,1}, {0,3,1,2}, {1,3,2,0}};

enum BADCELLTYPE
{
	ROUND = 1,		//Good
	NEEDLE = 2,
	WEDGE = 3,
	SPINDLE = 4,
	SLIVER = 5,
	CAP = 6
};

//#define _ANGLE_INTERVAL_COMP

/*
 * 基于内存和精度的不同需求，可编译为支持单精度或
 * 双精度两个不同版本
 * use double or float type to define float variables according to
 * memory & precision considerations
 */
#ifdef _DOUBLE
typedef double REAL;
#else 
typedef float REAL;
#endif /* _DOUBLE */

/*
 * 基于问题规模，可编译为支持普通整型和长整型两个不同版本
 * while more than billions elements are required, 
 * "long" type  is required
 */
#ifdef _LONG
typedef long INTEGER;
#else
typedef int INTEGER;
#endif /* _LONG */

/* *********************************************************************************
 * 定义一些常量和变量控制全局内存的动态分配
 * 由于输入数据中，表面三角片的数量是确定的，
 * 因此，根据作为基变量，确定其信息的内存分配大小
 * define some constants controlling memory allocation
 * As the number of facets is known at the beginning, it is
 * selected as the base variable to determined the memory sizes of other variables
 * **********************************************************************************/
/*
 * 将PROB_SIZE_UNIT个表面三角片的问题大小定义为1
 * the size of problem including PROB_SIZE_UNIT facets is defined as 1
 */
#define PROB_SIZE_UNIT ((INTEGER)(2000))

/*
 * 基于表面三角片确定其它信息内存分配大小
 * determing the momory size of other variables based on the number of facets
 */

/* 问题大小 problem size */
#define PROB_SIZE_BASED_ON_FAC(x) ((REAL) (x / PROB_SIZE_UNIT))

/* 
 * 一致边界恢复时，表面三角片的数量可能会增加，由此需要确定一个放大系数
 * 分配表面三角片的存储空间
 * The number of facets maybe increases for conforming boundary recovery procedure,
 * so an amplified factor is required to determine the memory size of facets
 */
#define FAC_AMPLIFY_FACTOR (1.2)
#define INIT_ALLOC_FAC_NUM(x) ((INTEGER) (FAC_AMPLIFY_FACTOR * x))

/* 表面约束边的数量: 根据欧式定理，边的数量＝3/2 * 面的数量
 * memory size of the surface edges 
 * it equals 3/2 multiplying the number of facets 
 */
#define SUR_EDG_RATIO (1.5)
#define INIT_ALLOC_SUR_EDG_NUM(x) ((INTEGER) (SUR_EDG_RATIO * x + 1))

/* 网格节点的数量 memory size of the mesh nodes */
#define NOD_RATIO (5.0)
#define INIT_ALLOC_NOD_NUM(x) ((INTEGER) (NOD_RATIO * x))

/* 单元的数量 memory size of the elements */
#define ELE_RATIO (20.0)
#define INIT_ALLOC_ELE_NUM(x) ((INTEGER) (ELE_RATIO * x))

/* 
 * 当内存不够时，重新分配空间和原空间的比例
 * ratio between newly-added memory size to the orginal one when memory bound is touched
 */
#define NUM_ADD_RATIO (0.25)

#define CHECK_FREE(x) if(x){ delete []x; x = NULL; }


/* *********************************************************************************
 * 小值定义(控制浮点误差)
 * definition of Epsilon values
 * **********************************************************************************/
#ifdef _DOUBLE
#define EPS_ABS_ZERO (1.0e-18)
#define EPS_ZERO    (1.e-6)
#define EPS_ZERO_SQ (1.e-12)
#else
#define EPS_ABS_ZERO (1.0e-12)
#define EPS_ZERO    (1.e-5)
#define EPS_ZERO_SQ (1.e-6)
#endif /* _DOUBLE */

/* 问题维度 problem dimensions */
#define DIM 3                  

/* ************************************************************************************
 * 基本数据结构定义 
 * Definitions for basic data structures 
 **************************************************************************************/

/* 坐标维度定义 coord. dimensions */
#define X 0
#define Y 1
#define Z 2

//优化相关
#define MESH_QUAL_THRESHOLD sin(ANGLE2RADIO(15.0))

#define NULL_NEIG -1 /* 无效邻居索引 invalid neighbor index */
#define NULL_ELEM -1 /* 无效单元索引 invalid elements index */

typedef REAL MYPOINT[DIM]; /* 坐标点 point */
typedef MYPOINT VECTOR;    /* 向量 vector */

typedef INTEGER NEIG_ARR[DIM+1]; /* 相邻单元索引 array of indices of neighboring elements */
typedef NEIG_ARR FORM_ARR;		 /* 形成节点索引 array of indices of forming points */
typedef MYPOINT FORM_PNTS[DIM+1];  /* 形成节点 array of forming points */

/* 网格节点 NODE */
typedef struct Node
{
	MYPOINT pt;			/* 坐标 coords. */
	REAL space;			/* 尺寸控制 space control */
	INTEGER iReserved;	/* 保留域，辅助信息 reserved domain, help info. */
	INTEGER iReserved2; /* 保留域，辅助信息 第1位标记点是否为wellshaped，第2位标记点是否为边界点 */
} Node;

/*
 * 网格单元(二维：三角形；三维：四面体)
 * element (triangle in 2D & tetrahedral in 3D) 
 */
typedef struct Elem
{
	MYPOINT cen;			/* 外接球球心 center of circumsphere */
	REAL rad;			/* 外接球半径 radius of circumsphere */
	NEIG_ARR neig;		/* 相邻单元 neighboring elements */
	FORM_ARR form;		/* 形成节点 forming points */  
	INTEGER iReserved;	/* 保留域，辅助信息 reserved domain, help info. */
	INTEGER iReserved2; /* 保留域，辅助信息 reserved domain, help info. */
	INTEGER iReserved3; /* 保留域，第1\2\3\4位标示面i是否尝试过做2-3变换 */
} Elem;

/* 表面三角形 surface triangle */
typedef struct SurTri
{
	INTEGER form[3];	/* 形成节点 forming points */
	INTEGER edgs[3];    /* 表面边 surface edges */
	INTEGER parent;		/* 母单元 */
	INTEGER iPatch;     /* 曲面片编号 */
	INTEGER subHead;	/* 第0位存储是否为已恢复边的信息；其它位保留对应打断后的子面信息 */
} SurTri;

/* 
 * 表面边
 * iNxt域的使用：
 * (1) 提取边数据时：以iStart为首点形成的边链表中的下一条边；
 * (2) 边界恢复时：边界边对应的破坏边的索引 
 * surface edge
 * How to use iNxt ?
 * (1) while obtaining the data of surface edge
 *       next edge in the list of edges with iStart as the first node
 * (2) while recovering boundaries (edges & faces)
 *       index for the destroyed edge
 */
typedef struct SurEdg
{
	INTEGER iStart, iEnd; /* 起点 & 终点 start & end point */
	INTEGER iLftF, iRgtF; /* 左面片/右面片 */
	INTEGER iNxt;		  /* 辅助信息 help info，哈希表指针，指向下一条边 */
	INTEGER subHead;	  /* 第0位存储是否已完全恢复的信息；
						   * 第1为存储是否没有边刺穿其内部
						   * 其它位保留对应打断后的子边信息 */
} SurEdg;

/* ************************************************************************************
 * 初始三角化数据 
 * Data for Initial Triangulation (IT)
 * ***********************************************************************************/

#define INIT_ELE_NUM 5    /* 初始单元数目 num. of eles. of the IT */
#define INIT_NOD_NUM 8    /* 初始节点数目 num. of nodes of the IT */
#define BOX_SIZE (4.0)	  /* 包围盒尺寸 scale of the box of the IT */	
	
/* 坐标数据 coords. */
const MYPOINT g_cors[] = 
{
	{-BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, {BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, 
	{BOX_SIZE, BOX_SIZE, -BOX_SIZE},   {-BOX_SIZE, BOX_SIZE, -BOX_SIZE},
	{-BOX_SIZE, -BOX_SIZE, BOX_SIZE},  {BOX_SIZE, -BOX_SIZE, BOX_SIZE},
	{BOX_SIZE, BOX_SIZE, BOX_SIZE},    {-BOX_SIZE, BOX_SIZE, BOX_SIZE}
};

/* 形成节点 forming points */
const FORM_ARR g_form[] =
{ 
	{0, 1, 2, 5}, 
	{2, 5, 6, 7},
	{0, 4, 5, 7},
	{0, 2, 3, 7},
	{0, 5, 2, 7}
};

/* 相邻关系　neighbor elements */
const NEIG_ARR g_neig[] =
{
	{NULL_NEIG, 4, NULL_NEIG, NULL_NEIG},
	{NULL_NEIG, NULL_NEIG, 4, NULL_NEIG},
	{NULL_NEIG, 4, NULL_NEIG, NULL_NEIG},
	{NULL_NEIG, NULL_NEIG, 4, NULL_NEIG},
	{1, 3, 2, 0}
};

/* ************************************************************
 * 归一化空间 
 * unitized space 
 * ************************************************************/
const MYPOINT g_minN = {-0.5, -0.5, -0.5};/* {0.0, 0.0, 0.0}; */
const MYPOINT g_maxN = {0.5, 0.5, 0.5}; /* {1.0, 1.0, 1.0} */

/* *************************************************************
 * 插入边界点 & 内部点：控制量定义 & 数据结构
 * insert boundary & inner points: definiton & data structures
 * ************************************************************/

 /* 
 * 局部数据结构: 空腔边
 * 更新新生成单元之间的相邻关系时使用
 * Local data structures: side in the cavity
 * Used to help update neigboring info. of a newly created element
 */
typedef struct CavSide
{
	INTEGER iEnd;	/* 终点 end point */
	int iNxt;		/* 下一条边 next side */
	INTEGER iEle;	/* 包含该边的第一个空腔单元 first cavity element including the side */
	int iNei;		/* 相邻位置　neighbor position */
} CavSide;

/* 扰动点信息 info. of distrubed point */
typedef struct DistInfo
{
	INTEGER iNod;  /* disturbing point index */
	MYPOINT old_pt;  /* old position of disturbing point, 
				      utilized for restoring the disturbance */
} DistInfo;

/* 
 * 节点坐标扰动比例(相对节点space值) 
 * ratio of disturbance of a node's coordinate submember to its space value
 */
#define DISTURB_RATIO (0.01)

/* 
 * 边界节点插入时的扰动控制: 
 * (1) 当一个循环中增加的节点和剩余节点数的比值低于MAX_ADD_RATIO时，直接扰动而不是推迟插入
 * (2) 当插入循环次数超过MAX_ADD_CIRCLE时，直接扰动而不是推迟插入
 * disturbance control of boundary point insertion
 * disturb it directly instead of postpone the operation if
 * (1) num. of successfully added points to that of remaining points less than MAX_ADD_RATIO
 * (2) num. of circles of boundary point insertion exceeds MAX_ADD_CIRCLE
 */
#define MAX_ADD_RATIO  0.1
#define MAX_ADD_CIRCLE 5

/* 
 * 当成功插入多少节点时在屏幕上显示信息
 * print info. while how many more points are inserted successfully inserted
 */
#define PRINT_PER_NOD 2000

/* *************************************************************
 * 边界边恢复：控制量定义 & 数据结构
 * recovery boundary edges: definiton & data structures
 * ************************************************************/

/* 管道大小最大值 max. size of a pipe */
#define MAX_PIPE_SIZE 512
/* 管道元分解后形成单元最大值 max. size of elements decomposed from a pipel */
#define MAX_DEC_SIZE  16
/* 球大小最大值 max. size of a sphere */
#define MAX_SPHERE_SIZE 10182//5096//
/* 壳大小最大值 max. size of a shell */
#define MAX_SHELL_SIZE 1024
/* 寻找管道时，候选单元数量的最大值 max. size of candidate eles. while searching a pipe */
#define MAX_SRCH_SIZE 10182//5096
/* 利用拓扑变换进行单元替换时，最大的单元大小容许值 */
#define MAX_LOCAL_MESH_SIZE 1024

#if (MAX_SRCH_SIZE < MAX_SPHERE_SIZE)
#define MAX_SRCH_SIZE MAX_SPHEE_SIZE
#endif

/* 球 sphere */
typedef INTEGER Sphere[MAX_SPHERE_SIZE];

/* 壳 shell */
typedef INTEGER Shell[MAX_SHELL_SIZE];

/* 管道元类型 type of a pipel */
#define DEG 0	
#define NOD 1	
#define EDG 2   
#define FAC 3	
#define NOD_NOD (((NOD) << 2) | NOD)	
#define EDG_NOD (((NOD) << 2) | EDG)
#define FAC_NOD (((NOD) << 2) | FAC)
#define NOD_EDG (((EDG) << 2) | NOD)
#define EDG_EDG (((EDG) << 2) | EDG)
#define FAC_EDG (((EDG) << 2) | FAC)
#define NOD_FAC (((FAC) << 2) | NOD)
#define EDG_FAC (((FAC) << 2) | EDG)
#define FAC_FAC (((FAC) << 2) | FAC)
#define NOD_DEG (((DEG) << 2) | NOD)
#define EDG_DEG (((DEG) << 2) | EDG)

/* 
 * 构建管道的辅助数据结构
 * help data structures for searching(constructing) a pipe 
 */
typedef struct PipeSrch
{
	INTEGER eles[MAX_SRCH_SIZE];	/* 候选单元 candidate elements */
	int codes[MAX_SRCH_SIZE];		/* 特征实体编码 codes for characteristic entities(node/edge/face) */
	int pres[MAX_SRCH_SIZE];		/* 搜索前后的比对信息 comparison info. bef. & aft. searching */
	int nEles;						/* 单元数目 num. of candidate elements */
	int type;						/* 类型 type */
	INTEGER iNod;					/* 遗失边上的当前打断点 separate point */
	int ec;							/* 由于浮点运算误差，和边相交可能会被误判为和面相交，保留
	                                   该边编码，以利回退处理 edge code for possible setback operation */
} PipeSrch;

/*  管道元分解后形成单元集合 a set of elements after decom. a pipel */
typedef INTEGER DecTets[MAX_DEC_SIZE];

/* 管道元 pipel */
typedef struct Pipel
{
	INTEGER iEle;					 /* 单元索引 element indices */
	int type;						 /* 类型 type */
	INTEGER iNod1, iNod2;			 /* 相关节点索引 indices of two node involved */
	int cod1, cod2;					 /* 相关实体编码 codes of two entities involved */ 
	int ec;							 /* 由于浮点运算误差，和边相交可能会被误判为和面相交，保留
	                                    该边编码，以利回退处理 edge code for possible setback operation */
	DecTets dectets;				 /* 分解后的单元索引 element indices after decomposition */
	int nTets;						 /* 分解后单元数量 number of eles. obtained from decom. */
	int flag;						 /* 单元有效标志（缺省值为0；如果其分解产生负体积单元，则flag = 1） element validity flag */
} Pipel;

/* 管道 pipe */
typedef Pipel Pipe[MAX_PIPE_SIZE];

/*
 * 边－边类型管道元的子类型：两条边是对边/邻边
 * subcases of a pipel typed EDG_EDG
 * SC_NEIG: two edges are neighboring
 * SC_CONT: two edges are opposite
 */
enum {SC_NEIG = 0, SC_OPPO};

/*
 * 分解一个两条边被打断的面时,存在两种分解方式：S型/N型
 * types of decomposing a face with two edges seperated: S or N
 */
enum {TYPE_S = 0, TYPE_N};

/* *************************************************************
 * 边界面恢复：控制量定义 & 数据结构
 * recovery boundary faces: definiton & data structures
 * ************************************************************/

/* 簇元类型 type of a clusterel */
#define CO_PLAN 0	/* 共面 co-planar */
#define ONE_EDG 1	/* 1条打断边 1 separated edge */
#define TWO_EDG 2	/* 2条打断边 2 separated edges */
#define THR_EDG 3	/* 3条打断边 3 separated edges */
#define FOU_EDG 4	/* 4条打断边 4 separated edges */

/* 簇大小的最大值 max. size of a cluster */
#define MAX_CLUSTER_SIZE 1024
/* 候选簇元大小最大值 max. size of candidate clusterels */
#define MAX_CLUS_CAND_SIZE 2048

/* 簇元 clusterel */
typedef struct Clusterel
{
	INTEGER iEle;				/* 单元索引 element index */
	int type;					/* 类型 type */
	int codes[4];				/* 相交边编码 codes for intersection edges */
	INTEGER nodes[4]; 			/* 节点编号 node indices */
	int ntypes[4];				/* 节点类型 node types*/
	int nTets;					/* 分解后的单元索引 element indices after decomposition */
	DecTets dectets;			/* 分解后单元数量 number of eles. obtained from decom. */
	int flag;						 /* 单元有效标志（缺省值为0；如果其分解产生负体积单元，则flag = 1） element validity flag */
} Clusterel;

/* 簇 cluster */
typedef Clusterel Cluster[MAX_CLUSTER_SIZE];

/* 候选簇元 candidate clusterels */
typedef INTEGER ClusCand[MAX_CLUS_CAND_SIZE];

/* 簇元边 Cluster Edge */
typedef struct CluEdg
{
	INTEGER iStart, iEnd;	/* 簇元边起点&终点 start & end point */
	INTEGER iNod;			/* 簇元边和遗失面的交点 */
} CluEdg;

/* 
 * 和一个遗失面相交的簇元边的最大数目 
 * max. size of cluster edges intersecting with a missing face 
 */
#define MAX_CLU_EDG_SIZE 512

/* 簇元边数组 array of cluster edges */
typedef CluEdg CluEdgArr[MAX_CLU_EDG_SIZE];

/* 边和遗失面的相交类型 type of an edge intersecting with a missing face */
#define NOD_NUL 0
#define NOD_EXT 1
#define NOD_BEG 2
#define NOD_END 3
#define NOD_MID 4

/* ***********************************************************************
 * 辅助点剔除：控制量定义 & 数据结构
 * supress extra-added points in the boundary: definiton & data structures
 * ************************************************************************/

/* 壳边数目最大值 max. size of shell edges */
#define MAX_SHELL_ED 128

/* 
 * 假设网格一条边AB，其对应的壳记为Shell(AB)，壳中除A、B以外的所有点构成一
 * 个三维环，每个ShellED对象表征环上一条边
 * Given an edge AB, its shell is denoted as Shell(AB), all mesh nodes * of which excluding A & B constructs a loop. Each object of ShellED * represents an edge of the loop. 
 */
typedef struct ShellED
{
	INTEGER iNod1, iNod2;			/* 两个端点 two nodes of the edge */
	INTEGER old_ou_tet, old_in_tet; /* 共享该边的两个四面体 
	                                   two tetrahedra (outer & inner) of the shell sharing the edge */
} ShellED;  

/* 
 * 破坏边和破坏面最多的加点个数 
 * max. number of extra-added points in a missing entity (edge/face) 
 */
#define MAX_ADD_NOD 512

/* 
 * 破坏边/破坏面和表面边和表面面片数目的比值, 用以
 * 指导破坏边和破坏面片的内存分配
 * (predicted) ratio between the num. of destroyed edges/faces to
 * that of total edges/faces, which is used to guide the memory
 * allocation of destroyed edges & faces
 */
#define DES_EDG_RATIO ((REAL)(0.01))
#define DES_FAC_RATIO ((REAL)(0.01))
#define MIN_DES_EDG ((INTEGER)20)
#define MIN_DES_FAC ((INTEGER)20)
#define DEF_DES_EDG(x) ((INTEGER) (DES_EDG_RATIO * x))
#define DEF_DES_FAC(x) ((INTEGER) (DES_FAC_RATIO * x))
/* x为表面边界总数 x represents the num. of edges */
#define INIT_DES_EDG(x) (DEF_DES_EDG(x) < MIN_DES_EDG ? MIN_DES_EDG : DEF_DES_EDG(x)) 
/* x为表面面片总数 x represents the num. of faces */
#define INIT_DES_FAC(x) (DEF_DES_FAC(x) < MIN_DES_FAC ? MIN_DES_FAC : DEF_DES_FAC(x))

/* 
 * 被破坏的表面边 destroyed surface edge 
 * 注意：边的两个端点也存在iAddNodes中
 * Note: end points of the edge are also included in iAddNodes
 */
typedef struct DesEdge
{
	INTEGER iEdg;					/* 表面边索引 index of the surface edge */
	INTEGER iAddNodes[MAX_ADD_NOD];	/* 边上节点序列 a sequence of nodes in the edge */
	int nNodes;
} DesEdge;

/* 被破坏的表面面片 destroyed surface facets */
typedef struct DesFacet
{	
	INTEGER iFac;					/* 表面面片索引 index of the surface edge */					
	INTEGER iAddNodes[MAX_ADD_NOD]; /* 面内节点 all extra-added nodes inner the missing facet */
	int nNodes;

	INTEGER subFts[MAX_CLUSTER_SIZE][3];
	INTEGER subNgs[MAX_CLUSTER_SIZE][3];
	INTEGER subEls[MAX_CLUSTER_SIZE][2];
	int nSubFts;
} DesFacet;

/* ***********************************************************************
 * 多余单元删除：控制量定义 & 数据结构
 * remove outer elements: definiton & data structures
 * ************************************************************************/
/*
 * 单元iReserved域的2～3位标志单元的几何属性：
 * flagging if an element is inner or outer the domain concerned
 * (two bits(No. 2~3) after the test bit in the iReserved of an Elem object)
 */
#define UNDEF	 (0)	/* 未定义 undefined */
#define	OUTER	 (1)	/* 区域外 outer the domain */
#define INNER	 (2)	/* 区域内 inner the domain */

/* ----------------------------------------------------
 * 以下数据结构记录一个空腔中和面片边相交的实体的数目 
 * -------------------------------------------------*/
#define MAX_POLY_INT_EDGE 1024
#define MAX_POLY_INT_FACE 1024
typedef struct IntEHash
{
	INTEGER iNod1, iNod2;
	INTEGER hashNxt;
} IntEHash;

typedef struct IntFHash
{
	INTEGER iNod1, iNod2, iNod3;
	INTEGER hashNxt;
} IntFHash;

/* 在采用动态编程法做edge removal操作时，我们使用了静态数组保持一些信息
 * MAX_SKIRT_POLY_SIZE限定了数组的最大大小
 */
#define MAX_SKIRT_POLY_SIZE 128 
#define MAX_ANGLE_INTERVAL_SIZE 45

/* types of quality measures that may be used */
enum QualityMetrics
{
    QUALMINSINE = 0,
    QUALRADIUSRATIO,
    QUALVLRMS3RATIO,
    QUALMEANSINE,
    QUALMINSINEANDEDGERATIO,
    QUALWARPEDMINSINE,
    QUALMINANGLE,/* merits below for quality evaluation only */
    QUALMAXANGLE,
	QUALALLANGLES,
	QUALSINEANGLE
};
#define WRAPPED_SINE_RATIO (1.0)

#endif /* __iso3d_define_h__ */
