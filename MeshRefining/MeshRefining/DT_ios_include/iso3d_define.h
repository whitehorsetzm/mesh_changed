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

#ifndef __iso3d_define_h__
#define __iso3d_define_h__

//ǰ����������ֹ���ָ����ĸ���,�м�������Ϊ����Ǳ�
const int tetra_pnt[6][4] = 
{ {0,1,2,3}, {0,2,3,1}, {0,3,1,2}, {1,2,0,3}, {1,0,3,2}, {2,0,1,3} };

const int dihed_edge[6] =
{ 5, 3, 4, 0, 2, 1 };

//ǰ�������ʾһ���ߣ�ǰ����������ֹ���ָ����ĸ���
const int tetra_edge[6][3] = 
{ {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 2, 0}, {1, 3, 2}, {2, 3, 0} };


//ǰ����������ֹ���ָ����ĸ���
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
 * �����ڴ�;��ȵĲ�ͬ���󣬿ɱ���Ϊ֧�ֵ����Ȼ�
 * ˫����������ͬ�汾
 * use double or float type to define float variables according to
 * memory & precision considerations
 */
#ifdef _DOUBLE
typedef double REAL;
#else 
typedef float REAL;
#endif /* _DOUBLE */

/*
 * ���������ģ���ɱ���Ϊ֧����ͨ���ͺͳ�����������ͬ�汾
 * while more than billions elements are required, 
 * "long" type  is required
 */
#ifdef _LONG
typedef long INTEGER;
#else
typedef int INTEGER;
#endif /* _LONG */

/* *********************************************************************************
 * ����һЩ�����ͱ�������ȫ���ڴ�Ķ�̬����
 * �������������У���������Ƭ��������ȷ���ģ�
 * ��ˣ�������Ϊ��������ȷ������Ϣ���ڴ�����С
 * define some constants controlling memory allocation
 * As the number of facets is known at the beginning, it is
 * selected as the base variable to determined the memory sizes of other variables
 * **********************************************************************************/
/*
 * ��PROB_SIZE_UNIT����������Ƭ�������С����Ϊ1
 * the size of problem including PROB_SIZE_UNIT facets is defined as 1
 */
#define PROB_SIZE_UNIT ((INTEGER)(2000))

/*
 * ���ڱ�������Ƭȷ��������Ϣ�ڴ�����С
 * determing the momory size of other variables based on the number of facets
 */

/* �����С problem size */
#define PROB_SIZE_BASED_ON_FAC(x) ((REAL) (x / PROB_SIZE_UNIT))

/* 
 * һ�±߽�ָ�ʱ����������Ƭ���������ܻ����ӣ��ɴ���Ҫȷ��һ���Ŵ�ϵ��
 * �����������Ƭ�Ĵ洢�ռ�
 * The number of facets maybe increases for conforming boundary recovery procedure,
 * so an amplified factor is required to determine the memory size of facets
 */
#define FAC_AMPLIFY_FACTOR (1.2)
#define INIT_ALLOC_FAC_NUM(x) ((INTEGER) (FAC_AMPLIFY_FACTOR * x))

/* ����Լ���ߵ�����: ����ŷʽ�����ߵ�������3/2 * �������
 * memory size of the surface edges 
 * it equals 3/2 multiplying the number of facets 
 */
#define SUR_EDG_RATIO (1.5)
#define INIT_ALLOC_SUR_EDG_NUM(x) ((INTEGER) (SUR_EDG_RATIO * x + 1))

/* ����ڵ������ memory size of the mesh nodes */
#define NOD_RATIO (5.0)
#define INIT_ALLOC_NOD_NUM(x) ((INTEGER) (NOD_RATIO * x))

/* ��Ԫ������ memory size of the elements */
#define ELE_RATIO (20.0)
#define INIT_ALLOC_ELE_NUM(x) ((INTEGER) (ELE_RATIO * x))

/* 
 * ���ڴ治��ʱ�����·���ռ��ԭ�ռ�ı���
 * ratio between newly-added memory size to the orginal one when memory bound is touched
 */
#define NUM_ADD_RATIO (0.25)

#define CHECK_FREE(x) if(x){ delete []x; x = NULL; }


/* *********************************************************************************
 * Сֵ����(���Ƹ������)
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

/* ����ά�� problem dimensions */
#define DIM 3                  

/* ************************************************************************************
 * �������ݽṹ���� 
 * Definitions for basic data structures 
 **************************************************************************************/

/* ����ά�ȶ��� coord. dimensions */
#define X 0
#define Y 1
#define Z 2

//�Ż����
#define MESH_QUAL_THRESHOLD sin(ANGLE2RADIO(15.0))

#define NULL_NEIG -1 /* ��Ч�ھ����� invalid neighbor index */
#define NULL_ELEM -1 /* ��Ч��Ԫ���� invalid elements index */

typedef REAL MYPOINT[DIM]; /* ����� point */
typedef MYPOINT VECTOR;    /* ���� vector */

typedef INTEGER NEIG_ARR[DIM+1]; /* ���ڵ�Ԫ���� array of indices of neighboring elements */
typedef NEIG_ARR FORM_ARR;		 /* �γɽڵ����� array of indices of forming points */
typedef MYPOINT FORM_PNTS[DIM+1];  /* �γɽڵ� array of forming points */

/* ����ڵ� NODE */
typedef struct Node
{
	MYPOINT pt;			/* ���� coords. */
	REAL space;			/* �ߴ���� space control */
	INTEGER iReserved;	/* �����򣬸�����Ϣ reserved domain, help info. */
	INTEGER iReserved2; /* �����򣬸�����Ϣ ��1λ��ǵ��Ƿ�Ϊwellshaped����2λ��ǵ��Ƿ�Ϊ�߽�� */
} Node;

/*
 * ����Ԫ(��ά�������Σ���ά��������)
 * element (triangle in 2D & tetrahedral in 3D) 
 */
typedef struct Elem
{
	MYPOINT cen;			/* ��������� center of circumsphere */
	REAL rad;			/* �����뾶 radius of circumsphere */
	NEIG_ARR neig;		/* ���ڵ�Ԫ neighboring elements */
	FORM_ARR form;		/* �γɽڵ� forming points */  
	INTEGER iReserved;	/* �����򣬸�����Ϣ reserved domain, help info. */
	INTEGER iReserved2; /* �����򣬸�����Ϣ reserved domain, help info. */
	INTEGER iReserved3; /* �����򣬵�1\2\3\4λ��ʾ��i�Ƿ��Թ���2-3�任 */
} Elem;

/* ���������� surface triangle */
typedef struct SurTri
{
	INTEGER form[3];	/* �γɽڵ� forming points */
	INTEGER edgs[3];    /* ����� surface edges */
	INTEGER parent;		/* ĸ��Ԫ */
	INTEGER iPatch;     /* ����Ƭ��� */
	INTEGER subHead;	/* ��0λ�洢�Ƿ�Ϊ�ѻָ��ߵ���Ϣ������λ������Ӧ��Ϻ��������Ϣ */
} SurTri;

/* 
 * �����
 * iNxt���ʹ�ã�
 * (1) ��ȡ������ʱ����iStartΪ�׵��γɵı������е���һ���ߣ�
 * (2) �߽�ָ�ʱ���߽�߶�Ӧ���ƻ��ߵ����� 
 * surface edge
 * How to use iNxt ?
 * (1) while obtaining the data of surface edge
 *       next edge in the list of edges with iStart as the first node
 * (2) while recovering boundaries (edges & faces)
 *       index for the destroyed edge
 */
typedef struct SurEdg
{
	INTEGER iStart, iEnd; /* ��� & �յ� start & end point */
	INTEGER iLftF, iRgtF; /* ����Ƭ/����Ƭ */
	INTEGER iNxt;		  /* ������Ϣ help info����ϣ��ָ�룬ָ����һ���� */
	INTEGER subHead;	  /* ��0λ�洢�Ƿ�����ȫ�ָ�����Ϣ��
						   * ��1Ϊ�洢�Ƿ�û�бߴ̴����ڲ�
						   * ����λ������Ӧ��Ϻ���ӱ���Ϣ */
} SurEdg;

/* ************************************************************************************
 * ��ʼ���ǻ����� 
 * Data for Initial Triangulation (IT)
 * ***********************************************************************************/

#define INIT_ELE_NUM 5    /* ��ʼ��Ԫ��Ŀ num. of eles. of the IT */
#define INIT_NOD_NUM 8    /* ��ʼ�ڵ���Ŀ num. of nodes of the IT */
#define BOX_SIZE (4.0)	  /* ��Χ�гߴ� scale of the box of the IT */	
	
/* �������� coords. */
const MYPOINT g_cors[] = 
{
	{-BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, {BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, 
	{BOX_SIZE, BOX_SIZE, -BOX_SIZE},   {-BOX_SIZE, BOX_SIZE, -BOX_SIZE},
	{-BOX_SIZE, -BOX_SIZE, BOX_SIZE},  {BOX_SIZE, -BOX_SIZE, BOX_SIZE},
	{BOX_SIZE, BOX_SIZE, BOX_SIZE},    {-BOX_SIZE, BOX_SIZE, BOX_SIZE}
};

/* �γɽڵ� forming points */
const FORM_ARR g_form[] =
{ 
	{0, 1, 2, 5}, 
	{2, 5, 6, 7},
	{0, 4, 5, 7},
	{0, 2, 3, 7},
	{0, 5, 2, 7}
};

/* ���ڹ�ϵ��neighbor elements */
const NEIG_ARR g_neig[] =
{
	{NULL_NEIG, 4, NULL_NEIG, NULL_NEIG},
	{NULL_NEIG, NULL_NEIG, 4, NULL_NEIG},
	{NULL_NEIG, 4, NULL_NEIG, NULL_NEIG},
	{NULL_NEIG, NULL_NEIG, 4, NULL_NEIG},
	{1, 3, 2, 0}
};

/* ************************************************************
 * ��һ���ռ� 
 * unitized space 
 * ************************************************************/
const MYPOINT g_minN = {-0.5, -0.5, -0.5};/* {0.0, 0.0, 0.0}; */
const MYPOINT g_maxN = {0.5, 0.5, 0.5}; /* {1.0, 1.0, 1.0} */

/* *************************************************************
 * ����߽�� & �ڲ��㣺���������� & ���ݽṹ
 * insert boundary & inner points: definiton & data structures
 * ************************************************************/

 /* 
 * �ֲ����ݽṹ: ��ǻ��
 * ���������ɵ�Ԫ֮������ڹ�ϵʱʹ��
 * Local data structures: side in the cavity
 * Used to help update neigboring info. of a newly created element
 */
typedef struct CavSide
{
	INTEGER iEnd;	/* �յ� end point */
	int iNxt;		/* ��һ���� next side */
	INTEGER iEle;	/* �����ñߵĵ�һ����ǻ��Ԫ first cavity element including the side */
	int iNei;		/* ����λ�á�neighbor position */
} CavSide;

/* �Ŷ�����Ϣ info. of distrubed point */
typedef struct DistInfo
{
	INTEGER iNod;  /* disturbing point index */
	MYPOINT old_pt;  /* old position of disturbing point, 
				      utilized for restoring the disturbance */
} DistInfo;

/* 
 * �ڵ������Ŷ�����(��Խڵ�spaceֵ) 
 * ratio of disturbance of a node's coordinate submember to its space value
 */
#define DISTURB_RATIO (0.01)

/* 
 * �߽�ڵ����ʱ���Ŷ�����: 
 * (1) ��һ��ѭ�������ӵĽڵ��ʣ��ڵ����ı�ֵ����MAX_ADD_RATIOʱ��ֱ���Ŷ��������Ƴٲ���
 * (2) ������ѭ����������MAX_ADD_CIRCLEʱ��ֱ���Ŷ��������Ƴٲ���
 * disturbance control of boundary point insertion
 * disturb it directly instead of postpone the operation if
 * (1) num. of successfully added points to that of remaining points less than MAX_ADD_RATIO
 * (2) num. of circles of boundary point insertion exceeds MAX_ADD_CIRCLE
 */
#define MAX_ADD_RATIO  0.1
#define MAX_ADD_CIRCLE 5

/* 
 * ���ɹ�������ٽڵ�ʱ����Ļ����ʾ��Ϣ
 * print info. while how many more points are inserted successfully inserted
 */
#define PRINT_PER_NOD 2000

/* *************************************************************
 * �߽�߻ָ������������� & ���ݽṹ
 * recovery boundary edges: definiton & data structures
 * ************************************************************/

/* �ܵ���С���ֵ max. size of a pipe */
#define MAX_PIPE_SIZE 512
/* �ܵ�Ԫ�ֽ���γɵ�Ԫ���ֵ max. size of elements decomposed from a pipel */
#define MAX_DEC_SIZE  16
/* ���С���ֵ max. size of a sphere */
#define MAX_SPHERE_SIZE 10182//5096//
/* �Ǵ�С���ֵ max. size of a shell */
#define MAX_SHELL_SIZE 1024
/* Ѱ�ҹܵ�ʱ����ѡ��Ԫ���������ֵ max. size of candidate eles. while searching a pipe */
#define MAX_SRCH_SIZE 10182//5096
/* �������˱任���е�Ԫ�滻ʱ�����ĵ�Ԫ��С����ֵ */
#define MAX_LOCAL_MESH_SIZE 1024

#if (MAX_SRCH_SIZE < MAX_SPHERE_SIZE)
#define MAX_SRCH_SIZE MAX_SPHEE_SIZE
#endif

/* �� sphere */
typedef INTEGER Sphere[MAX_SPHERE_SIZE];

/* �� shell */
typedef INTEGER Shell[MAX_SHELL_SIZE];

/* �ܵ�Ԫ���� type of a pipel */
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
 * �����ܵ��ĸ������ݽṹ
 * help data structures for searching(constructing) a pipe 
 */
typedef struct PipeSrch
{
	INTEGER eles[MAX_SRCH_SIZE];	/* ��ѡ��Ԫ candidate elements */
	int codes[MAX_SRCH_SIZE];		/* ����ʵ����� codes for characteristic entities(node/edge/face) */
	int pres[MAX_SRCH_SIZE];		/* ����ǰ��ıȶ���Ϣ comparison info. bef. & aft. searching */
	int nEles;						/* ��Ԫ��Ŀ num. of candidate elements */
	int type;						/* ���� type */
	INTEGER iNod;					/* ��ʧ���ϵĵ�ǰ��ϵ� separate point */
	int ec;							/* ���ڸ����������ͱ��ཻ���ܻᱻ����Ϊ�����ཻ������
	                                   �ñ߱��룬�������˴��� edge code for possible setback operation */
} PipeSrch;

/*  �ܵ�Ԫ�ֽ���γɵ�Ԫ���� a set of elements after decom. a pipel */
typedef INTEGER DecTets[MAX_DEC_SIZE];

/* �ܵ�Ԫ pipel */
typedef struct Pipel
{
	INTEGER iEle;					 /* ��Ԫ���� element indices */
	int type;						 /* ���� type */
	INTEGER iNod1, iNod2;			 /* ��ؽڵ����� indices of two node involved */
	int cod1, cod2;					 /* ���ʵ����� codes of two entities involved */ 
	int ec;							 /* ���ڸ����������ͱ��ཻ���ܻᱻ����Ϊ�����ཻ������
	                                    �ñ߱��룬�������˴��� edge code for possible setback operation */
	DecTets dectets;				 /* �ֽ��ĵ�Ԫ���� element indices after decomposition */
	int nTets;						 /* �ֽ��Ԫ���� number of eles. obtained from decom. */
	int flag;						 /* ��Ԫ��Ч��־��ȱʡֵΪ0�������ֽ�����������Ԫ����flag = 1�� element validity flag */
} Pipel;

/* �ܵ� pipe */
typedef Pipel Pipe[MAX_PIPE_SIZE];

/*
 * �ߣ������͹ܵ�Ԫ�������ͣ��������ǶԱ�/�ڱ�
 * subcases of a pipel typed EDG_EDG
 * SC_NEIG: two edges are neighboring
 * SC_CONT: two edges are opposite
 */
enum {SC_NEIG = 0, SC_OPPO};

/*
 * �ֽ�һ�������߱���ϵ���ʱ,�������ַֽⷽʽ��S��/N��
 * types of decomposing a face with two edges seperated: S or N
 */
enum {TYPE_S = 0, TYPE_N};

/* *************************************************************
 * �߽���ָ������������� & ���ݽṹ
 * recovery boundary faces: definiton & data structures
 * ************************************************************/

/* ��Ԫ���� type of a clusterel */
#define CO_PLAN 0	/* ���� co-planar */
#define ONE_EDG 1	/* 1����ϱ� 1 separated edge */
#define TWO_EDG 2	/* 2����ϱ� 2 separated edges */
#define THR_EDG 3	/* 3����ϱ� 3 separated edges */
#define FOU_EDG 4	/* 4����ϱ� 4 separated edges */

/* �ش�С�����ֵ max. size of a cluster */
#define MAX_CLUSTER_SIZE 1024
/* ��ѡ��Ԫ��С���ֵ max. size of candidate clusterels */
#define MAX_CLUS_CAND_SIZE 2048

/* ��Ԫ clusterel */
typedef struct Clusterel
{
	INTEGER iEle;				/* ��Ԫ���� element index */
	int type;					/* ���� type */
	int codes[4];				/* �ཻ�߱��� codes for intersection edges */
	INTEGER nodes[4]; 			/* �ڵ��� node indices */
	int ntypes[4];				/* �ڵ����� node types*/
	int nTets;					/* �ֽ��ĵ�Ԫ���� element indices after decomposition */
	DecTets dectets;			/* �ֽ��Ԫ���� number of eles. obtained from decom. */
	int flag;						 /* ��Ԫ��Ч��־��ȱʡֵΪ0�������ֽ�����������Ԫ����flag = 1�� element validity flag */
} Clusterel;

/* �� cluster */
typedef Clusterel Cluster[MAX_CLUSTER_SIZE];

/* ��ѡ��Ԫ candidate clusterels */
typedef INTEGER ClusCand[MAX_CLUS_CAND_SIZE];

/* ��Ԫ�� Cluster Edge */
typedef struct CluEdg
{
	INTEGER iStart, iEnd;	/* ��Ԫ�����&�յ� start & end point */
	INTEGER iNod;			/* ��Ԫ�ߺ���ʧ��Ľ��� */
} CluEdg;

/* 
 * ��һ����ʧ���ཻ�Ĵ�Ԫ�ߵ������Ŀ 
 * max. size of cluster edges intersecting with a missing face 
 */
#define MAX_CLU_EDG_SIZE 512

/* ��Ԫ������ array of cluster edges */
typedef CluEdg CluEdgArr[MAX_CLU_EDG_SIZE];

/* �ߺ���ʧ����ཻ���� type of an edge intersecting with a missing face */
#define NOD_NUL 0
#define NOD_EXT 1
#define NOD_BEG 2
#define NOD_END 3
#define NOD_MID 4

/* ***********************************************************************
 * �������޳������������� & ���ݽṹ
 * supress extra-added points in the boundary: definiton & data structures
 * ************************************************************************/

/* �Ǳ���Ŀ���ֵ max. size of shell edges */
#define MAX_SHELL_ED 128

/* 
 * ��������һ����AB�����Ӧ�ĿǼ�ΪShell(AB)�����г�A��B��������е㹹��һ
 * ����ά����ÿ��ShellED�����������һ����
 * Given an edge AB, its shell is denoted as Shell(AB), all mesh nodes * of which excluding A & B constructs a loop. Each object of ShellED * represents an edge of the loop. 
 */
typedef struct ShellED
{
	INTEGER iNod1, iNod2;			/* �����˵� two nodes of the edge */
	INTEGER old_ou_tet, old_in_tet; /* ����ñߵ����������� 
	                                   two tetrahedra (outer & inner) of the shell sharing the edge */
} ShellED;  

/* 
 * �ƻ��ߺ��ƻ������ļӵ���� 
 * max. number of extra-added points in a missing entity (edge/face) 
 */
#define MAX_ADD_NOD 512

/* 
 * �ƻ���/�ƻ���ͱ���ߺͱ�����Ƭ��Ŀ�ı�ֵ, ����
 * ָ���ƻ��ߺ��ƻ���Ƭ���ڴ����
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
/* xΪ����߽����� x represents the num. of edges */
#define INIT_DES_EDG(x) (DEF_DES_EDG(x) < MIN_DES_EDG ? MIN_DES_EDG : DEF_DES_EDG(x)) 
/* xΪ������Ƭ���� x represents the num. of faces */
#define INIT_DES_FAC(x) (DEF_DES_FAC(x) < MIN_DES_FAC ? MIN_DES_FAC : DEF_DES_FAC(x))

/* 
 * ���ƻ��ı���� destroyed surface edge 
 * ע�⣺�ߵ������˵�Ҳ����iAddNodes��
 * Note: end points of the edge are also included in iAddNodes
 */
typedef struct DesEdge
{
	INTEGER iEdg;					/* ��������� index of the surface edge */
	INTEGER iAddNodes[MAX_ADD_NOD];	/* ���Ͻڵ����� a sequence of nodes in the edge */
	int nNodes;
} DesEdge;

/* ���ƻ��ı�����Ƭ destroyed surface facets */
typedef struct DesFacet
{	
	INTEGER iFac;					/* ������Ƭ���� index of the surface edge */					
	INTEGER iAddNodes[MAX_ADD_NOD]; /* ���ڽڵ� all extra-added nodes inner the missing facet */
	int nNodes;

	INTEGER subFts[MAX_CLUSTER_SIZE][3];
	INTEGER subNgs[MAX_CLUSTER_SIZE][3];
	INTEGER subEls[MAX_CLUSTER_SIZE][2];
	int nSubFts;
} DesFacet;

/* ***********************************************************************
 * ���൥Ԫɾ�������������� & ���ݽṹ
 * remove outer elements: definiton & data structures
 * ************************************************************************/
/*
 * ��ԪiReserved���2��3λ��־��Ԫ�ļ������ԣ�
 * flagging if an element is inner or outer the domain concerned
 * (two bits(No. 2~3) after the test bit in the iReserved of an Elem object)
 */
#define UNDEF	 (0)	/* δ���� undefined */
#define	OUTER	 (1)	/* ������ outer the domain */
#define INNER	 (2)	/* ������ inner the domain */

/* ----------------------------------------------------
 * �������ݽṹ��¼һ����ǻ�к���Ƭ���ཻ��ʵ�����Ŀ 
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

/* �ڲ��ö�̬��̷���edge removal����ʱ������ʹ���˾�̬���鱣��һЩ��Ϣ
 * MAX_SKIRT_POLY_SIZE�޶������������С
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
