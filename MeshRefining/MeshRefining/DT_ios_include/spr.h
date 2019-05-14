#ifndef __spr_h__
#define __spr_h__

#pragma warning(disable:4786)
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <map>
#include <stack>
#include <queue>
#include <set>
#include <assert.h>

#include "celltopology.h"
#include "gridedgehash.h"
#include "geom_func.h"
#include "instentlist.h"

#if defined(WIN32) && !defined(__cplusplus)
#define inline __inline
#endif

typedef double REAL;
typedef double APOINT[3];	/* 定义坐标点POINT，这里的POINT即表征了一个double的三维数组*/
typedef int TRI[3];			//定义一个三角形，增加的那一个元素是保存这个三角形是属于哪个整体面的三角形
typedef int NEG[3];			/* 定义相邻关系 */ 
typedef int TETRA[4];		//定义一个四面体

#define NOD 1	
#define EDG 2   
#define FAC 3
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define ANGLE2RADIO(x) (x*PI)/180
#define RADIO2ANGLE(x) (x*180.0)/PI
#define SPR_ZERO_QUALITY 1.0e-16			/* default: 1.0e-7小于这个值的单元质量被截断为该值(正面效果：提高效率；负面效果：降低质量) */
#define SPR_UNALLOWED_QUALITY 0.0//1.0e-18		/* 取绝对值0会引起一些问题，除非有特别的拓扑限制 0.0			/* 小于这个值的单元被视为0体积单元 */
#define SPR_EPS_ZERO 1e-8
#define SPR_EPS_ZERO_SQ 1.e-8
#define MAX_SPR_POLY_SIZE 400//default is 300
#define MAX_CONSTRAINT_SIZE 300
#define MAX_COPLANE_NODES 512 // = MAX_ADD_NODE (defined in iso3d_define.h)
#define MAX_ALLOWED_POLY_SIZE 41
#define MAX_SIMP_POLY_NUM	32 /* 最多的子多边形数目 */	
#define MAX_SUB_PROBLEM_NUM 128
#define MAX_VOLUME 1.e16
#define MAX_SPR_INT_EDGE 128
#define MAX_SPR_INT_FACE 128

// #define MAX_S 20
#define MAX_VALUE(x, y) ((x) > (y) ? (x) : (y))
#define MIN_VALUE(x, y) ((x) < (y) ? (x) : (y))
#define SWAP_VALUE(x, y, t)	\
	{						\
		t = x;				\
		x = y;				\
		y = t;				\
	}

extern int count;		//四面体有效性的判断次数
extern int real;		//实际的递归次数
extern int box_n;		//用box消掉的面面相交的次数
extern int face_n;		//总共的面面相交的判断次数
extern int line_n;		//总共的线面相交的判断次数
extern int box_l_n;	//用box消掉的线面相交的次数
extern int float_count;	//总共的用浮点计算判断的面面相交的次数
extern long orient3d_n;		//计算体积的次数
extern long Call_V;	//试图查找次数
extern long Call_V_R;	//从已知数据中调用的次数
extern long Call[10];

extern double start, stop;	//记录整体时间
extern double duration;

extern double starttmp[10], stoptmp[10];	//记录细节时间
extern double dura[10];
extern double quality_calc_time;

#if 0
double Visible[MAX_S][MAX_S][MAX_S][MAX_S];
#endif

enum SPR_OPTIMAL_TYPE {
	MAX_MIN_CELL_QUALITY = 0,
	MIN_INT_NODE_NUMBER,
	EDGE_REMOVAL
};

typedef  APOINT APOINT_ARRAY[MAX_SPR_POLY_SIZE];
extern APOINT_ARRAY vertices;
extern int num_vertices;
extern bool spr_opt_solution;
extern double spr_time_limit;		/* 单次操作的时间限制(指单个多面体，一个问题可能包含多个多面体，time_limit = 0.0表示不限制) */
extern int spr_optimal_type;
extern int spr_max_int_node;
extern bool spr_verbose_output;
extern bool spr_no_more_recur;

/* for test only */
extern int spr_l2g[MAX_SPR_POLY_SIZE];

typedef struct polyhedron
{
        TRI t[MAX_SPR_POLY_SIZE];			/* 面片的形成节点 */
	NEG n[MAX_SPR_POLY_SIZE];			/* 面片的相邻关系 */
	int used_v[MAX_SPR_POLY_SIZE];		/* 包含节点的相邻面片数目，为0时，表示对应节点不在当前多面体中 */
	int vprt[MAX_SPR_POLY_SIZE];		/* 前2位为父亲节点，指示节点在单元中的编码，包含节点的一个单元 */
	int surface[MAX_SPR_POLY_SIZE];		/* 面片所在的曲面片编号，实际没有用到 */
	int num_t;							/* 面片数目 */
	int num_vable;						/* used_v值大于0的顶点个数 */
	int sort_f[MAX_SPR_POLY_SIZE];		/* 对面片进行排序的结果(用于加速SPR操作) */ 
	int sort_n[MAX_SPR_POLY_SIZE];		/* 对顶点进行排序的结果(用于加速SPR操作)，实际没有用到 */
	int newf_fg[MAX_SPR_POLY_SIZE];		/* 新增加的面的flag 0：非新增加的面；1：新增加的面 */
	int newf_id[3];						/* 新增加的面，不超过3个 */
	int newf_ct;						/* 新增加面的数目 */
	int allowedFaceNum;					/* 允许包含边的面片最大值，当spr_optimal_type = EDG_REMOVAL时启用 */ 
} polyhedron;

typedef struct divide
{
	TETRA te[MAX_SPR_POLY_SIZE];
	int num;
	double q;
	/* 为统计实际的交点个数，需要设计哈希表 */
	int intEdges[MAX_SPR_INT_EDGE][2], intFaces[MAX_SPR_INT_FACE][3];
	int intEdgeNum, intFaceNum;
	int intNum;				/* 实际的交点个数，当spr_optimal_type = MIN_INT_NODE_NUMBER时启用 */
	int facNum;				/* 包含边的面片个数，当spr_optimal_type = EDG_REMOVAL时启用 */ 
} divide;

typedef struct divide_record
{
	int num;
	int div_rd[MAX_SPR_POLY_SIZE];
} divide_record;

typedef struct SubDivideProblem
{
	polyhedron *poly;	/* 父亲多面体 */
	int tetraID;		/* 四面体编号, 由其四个节点组成 */
	int fIdx;			/* 面片编号  */
	int nIdx;			/* 节点编号  */
	double initQ;		/* 质量阀值  */
	double tetrQ;		/* 四面体单元的质量 */
	double diviQ;		/* 分解的质量，用于排序 */
#ifdef _SPR_USING_STD_VECTOR
	std::vector<polyhedron*> vec_simp_poly;	/* 描述子问题的多面体 */
#else
	polyhedron *simpPolyArray[MAX_SIMP_POLY_NUM];
	int simpPolyNum;
#endif
	divide_record	divR;					/* 用于记录以前分解结果的信息 */
	int intNum;			/* 交点的个数(当spr_optimal_type = MIN_INT_NODE_NUMBER时有效) */
	int facNum;			/* 包含边的面片个数，当spr_optimal_type = EDG_REMOVAL时启用 */ 
} SubDivideProblem;

/* 将SPR算法用于Steiner点删除时，需要增加一些约束 */
typedef struct Constraints
{
	int lineCsts[MAX_CONSTRAINT_SIZE][2];		/* 线约束的首末端点 */
	int lineCst_size;							/* 线约束的个数 */
	int faceCsts[MAX_CONSTRAINT_SIZE][3];		/* 面约束的3个端点 */
	int faceCst_size;							/* 面约束的个数 */
	int colineNodes[MAX_COPLANE_NODES];			/* 共线节点列表 */
	int colineN_size;							/* 共线节点个数 */
	int cofaceNodes[MAX_COPLANE_NODES];			/* 和面约束不冲突的共面节点列表（面约束所在的面包含的其它Steiner点列表）*/
	int cofaceN_size;							/* 共面节点的个数 */
	int badVolmCsts[MAX_CONSTRAINT_SIZE][4];	/* 不希望存在于面中的四面体单元 */
	int badVolmCst_size;						/* 体约束的个数 */
} Constraints; 

extern Constraints spr_csts;

bool operator == (const divide_record & t1, const divide_record& t2);

bool operator < (const divide_record & t1, const divide_record& t2);

bool operator > (const divide_record & t1, const divide_record& t2);

/* -----------------------------------------------------------------------------
 * 实现set类的运行时比较功能
 * ----------------------------------------------------------------------------*/
template <class T>
class RuntimeCmp 
{
public:
	enum cmp_mode {normal, reverse};
	
private:
	cmp_mode mode;
	
public:
	// constructor for sorting criterion
	// - default criterion uses value normal
	
	RuntimeCmp(cmp_mode m = normal): mode(m) {}
	
	// comparison of elements
	bool operator() (const T& t1, const T&t2) const {
		return mode == normal ? t1 < t2 : t2 < t1;
	}
	
	//comparison of sorting criteria
	bool operator == (const RuntimeCmp& rc) {
		return mode == rc.mode;
	}
};

typedef struct InsertedTetra
{
        TRI t[4];			/* 4个面片 */
	NEG n[4];			/* 4个面片之间的相邻关系 */
	int idx[4];			/* 4个面片的编号，小于0意味着会被增加, 大于0会被删除 */
	int num_del;		/* 4个面片中有多少个面片会被删除 */
	int flg[MAX_SPR_POLY_SIZE];		/* 确定哪些面是否被删除的标志 -1: 不被删除, 0~2:被删除面的局部索引 */
	bool complex;		/* 是否产生复边 */
} InsertedTetra;
 
typedef struct IDX_QUALITY
{
	int idx;
	double quality;
} IDX_QUALITY; 

/* ----------------------------------------------------------------------------
 * 为增强边界恢复算法，我们需要增强SPR算法
 * --------------------------------------------------------------------------*/
enum EDG_CONS_TYPE {MUST_EXIST = 0, MUST_BE_EXCLUDED, DEPEND_ON_LIST};
typedef struct EdgeConstraint
{
	int i1, i2;						   /* 边所在的2个端点 */
	int j1, j2, j3;					   /* 必须消除的边或面，如果i3 = -1表示是一条边 */
	EDG_CONS_TYPE type;				   /* 类型 */
	InstEntList *pListAllowed;		   /* 允许相交的实体 */
	InstEntList *pListExcluded;		   /* 不允许相交的实体 */
} EdgeConstraint;
extern EdgeConstraint edgeConsArray[MAX_CONSTRAINT_SIZE];
extern int edgeConsArraySize;
extern int spr_shell_edge_nodeS, spr_shell_edge_nodeE;

/* spr_optimal_type = EDGE_REMOVAL时，初始壳包含边的数目 */
extern int spr_shell_size;
extern int spr_allowed_face_num;	/* 记录子问题允许包含的面片数最大值，临时变量 */
extern int spr_subprob_face_num;	/* 记录子问题实际包含的面片数，临时变量 */
extern int spr_dependent_edges[MAX_CONSTRAINT_SIZE][3]; /* 第3个成员为类型。=0：没有顶点为壳边顶点；=1：有一个顶点为壳边端点；*/
extern int spr_dependent_edges_cnt;

#ifdef _SPR_MAP_TETRA_VOLUME
extern std::map<int, double> mapTetraVolume;
#endif
//std::map<int, int>	  mapFaceInst;

void init(double *q, struct polyhedron *P, struct divide *T);
int optimal(double init_q, struct polyhedron *P, struct divide *T, divide_record *divRd);
int optimal(double init_q, std::vector<struct polyhedron*> &vec_simp_poly, struct divide *T, divide_record *divRd);
double distance_square(APOINT p1, APOINT p2);
int istetra(int ELE[], struct polyhedron *P, double q, double *quality_ELE,
				   struct polyhedron *Q, int f, Constraints *csts = NULL);	

#ifdef _SPR_USING_STD_VECTOR
int istetra(int ELE[], struct polyhedron *P, double q, double *quality_ELE,
				   struct polyhedron *Q,
				   std::vector<struct polyhedron *>& vec_simp_poly, int f, bool *freeSimpPoly,
				   Constraints *csts = NULL, int maxIntNum = 0, int *intNum = NULL);
#else
int istetra(int ELE[], struct polyhedron *P, double q, double *quality_ELE,
				   struct polyhedron *Q,
				   struct polyhedron *simpPolyArray[], int *simpPolyNum,
				   int f, bool *freeSimpPoly,
				   Constraints *csts = NULL, int maxIntNum = 0, int *intNum = NULL);
#endif /* _SPR_USING_STD_VECTOR */

void sphere(struct polyhedron *P);
double areas(double xa, double ya, double za, 
		   double xb, double yb, double zb,
		   double xc, double yc, double zc);
double isVisible_SPR(APOINT p1, APOINT p2, APOINT p3, APOINT p4);

int delELE(int ELE[4], struct polyhedron *P, struct polyhedron *Q, int f, int *newCnt);

#ifdef _SPR_USING_STD_VECTOR
int delELE(int ELE[4], struct polyhedron *P, 
			struct polyhedron *Q,
			std::vector<struct polyhedron *>& vec_simp_poly, int f, bool *freeSimpPoly,
				   Constraints *csts = NULL, int maxIntNum = 0, int *intNum = NULL);
#else
int delELE(int ELE[4], struct polyhedron *P, 
			struct polyhedron *Q,
			struct polyhedron* simpPolyArray[], int *simpPolyNum,
			int f, bool *freeSimpPoly,
			Constraints *csts = NULL, int maxIntNum = 0, int *intNum = NULL);
#endif

int isNewFacetsValid_BruteForce(polyhedron *poly, int newCnt,
								Constraints *csts = NULL, int maxIntNum = 0, int *intNum = NULL,
								InsertedTetra *instTetra = NULL, polyhedron *prtpoly = NULL);

int isNewFacetsValid_BFS(polyhedron *poly, int newCnt);
void init_poly(struct polyhedron *P);
void mergeT(struct divide *T, int ELE[]);
void mergeT(struct divide *sT, struct divide *tT);
void octahedron(struct polyhedron *P);
void copyT(struct divide *Tc, struct divide *T);
void init_tetra(struct divide *T);
/*
void init_nodef(struct node_f *header);
*/
int saveConstraints(const char *fname, Constraints *csts, double iniQ);
int readConstraints(const char *fname, Constraints *csts, double *iniQ);
int openpls(char *address, struct polyhedron *P, bool reverse = true);
int savepl3(char *address, struct polyhedron *P, struct divide *T);
void findparent(int parent[], struct polyhedron *P, struct divide *T);
int isintersect(int iface[], int iline[], APOINT facep[], APOINT linep[], int colineNodes[], int colineNSize);
int isintersect(int iface[], APOINT facep[], Constraints *cst);
int isIntersectEdgeConstraints(int iface[], APOINT facep[]);
int isintersect(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], struct polyhedron *P);
int isintersect_old(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], struct polyhedron *P);
void rearrage(TRI t[], int n);
int savepls(char *address, struct polyhedron *P, bool bAllNode = false);
int saveneg(char *address, struct polyhedron *P);
int isborder(TRI face1, TRI face2, int m, int n, int *l, int *k);
int islinecross(int index1, int index2, int index3, int index4, APOINT A, APOINT B, APOINT C, APOINT D, struct polyhedron *P);

int islinecross_old(APOINT A, APOINT B, APOINT C, APOINT D);
int islinecross_new(APOINT A, APOINT B, APOINT C, APOINT D);
int islinecross_newer(int index1, int index2, int index3, int index4, struct polyhedron *P);

int isboxinter(APOINT facep1[], APOINT facep2[]);
int isboxinter_l(APOINT A, APOINT B, APOINT facep2[]);
int isintercross(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], int m, int n, int *res, struct polyhedron *P);
int islinecrossface(TRI face1, TRI face2, APOINT facep1[], APOINT facep2[], int m, int n, int res, struct polyhedron *P);
/*这里面的ln1和ln2是线的两个端点，fac1，fac2，fac3是那个面的三个坐标点，d是这条线的长度，areaa是这个面的面积，
通过参量可以返回的值是这个pnt，返回的是如果线和面相交，那么交点的情况。val是相交的具体情况
w1，w2，w3具体是指什么？
在函数中，当返回一个NOD（1）的时候，表明这条线和这个面相交的情况是这条线穿过了面的端点
当返回一个EDG（2）的时候，表明这条线和这个面相交的情况是这条线穿过了面的边界线
当返回一个FAC（3）的时候，表明这条线和这个面直接相交
当返回-1的时候，表明这条线和这个面不相交
当返回-2的时候，表明程序出现严重错误，但是具体还是不大懂*/
int lnFacInt(double ln1[], double ln2[], 
				double fac1[], double fac2[], double fac3[], 
				double d, double areaa,
				double pnt[], int *val,
				double *w1, double *w2, double *w3);
int lnFacInt_own(int line1, int line2, TRI face, APOINT ln1, APOINT ln2, APOINT fac1, APOINT fac2, APOINT fac3, struct polyhedron *P);
double isVisible_Search(int index1, int index2, int index3, int index4, double *volume = NULL);
double isVisible_Search(APOINT p1, APOINT p2, APOINT p3, APOINT p4);
int isNeigborAllowed(int fidx, int nidx, polyhedron *poly);

double SPRlogTime(void);
void SPRlogTime_Init();

inline int sort_quad_pair(int *j1, int *j2, int *j3, int *j4);
inline int make_quad_pair(int j1, int j2, int j3, int j4);

/* ----------------------------------------------------------------------
 * 建立复杂多面体的相邻关系，如果GridEdgeHash存在复边，假定复边对应的
 * 单元是有序的
 * -------------------------------------------------------------------*/
int build_poly_neig(struct polyhedron *poly, GridEdgeHash *eHash);

int build_manifold_poly_neig(struct polyhedron *poly);

/* ----------------------------------------------------------------------
 * 分解复杂多面体为多个简单多面体
 * 基于相邻关系的着色算法
 * -------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
int divide_poly(struct polyhedron *cmpx_poly, 
				std::vector<struct polyhedron*>& vec_simp_poly);
#else
int divide_poly(struct polyhedron *cmpx_poly, 
				struct polyhedron *simpPolyArray[], int *simpPolyNum);
#endif /* _SPR_USING_STD_VECTOR */

/* ----------------------------------------------------------------------
 * 获得首个包含某节点的父亲单元
 * -------------------------------------------------------------------*/
inline void get_node_first_parent(struct polyhedron *poly,
					int nidx, int *tidx, int *tcod);

inline int node_face_code(struct polyhedron *poly, int nidx, int fidx);

inline int edge_face_code(struct polyhedron *poly, 
						  int idx1, int idx2, int fidx,
						  int *ecod, int *eort);
inline int add_face_edge_pair(struct polyhedron *poly, 
						int idx1, int idx2, int fidx,
						std::vector<int>& vecp);
/* ----------------------------------------------------------------------
 * 获得所有包含某节点的父亲单元
 * -------------------------------------------------------------------*/
inline void get_node_parents(struct polyhedron *poly, int nidx, 
				std::vector<int>& vecp);

inline void get_node_parents(struct polyhedron *poly, int nidx, int *vecp, int *num);

/* ----------------------------------------------------------------------
 * 获得所有包含某条边的所有父亲单元
 * -------------------------------------------------------------------*/
inline void get_edge_parents(struct polyhedron *poly, int idx1, int idx2,
				std::vector<int>& vecp);

inline int get_edge_parents(int i1, int i2, std::vector<int>& vecp, 
							 polyhedron *poly, int *face1, int *face2, 
							 int *code1, int *code2, int *nume);

inline int get_edge_parents(int i1, int i2, int nodeParents[], int numParents, 
							 polyhedron *poly, int *face1, int *face2, 
							 int *code1, int *code2, int *nume);

inline int get_edge_parents(int i1, int i2, std::vector<int>& vecp, 
							 polyhedron *poly, int *face1, int *face2, 
							 int *code1, int *code2);

/* ----------------------------------------------------------------------
 * 已知：存在一个复杂多面体，每个面都是有方向的，且法向指向形体内部。其某
 * 些边由4个面共享，这些面形成两对面，每对面各属于一个简单多面体。简单多面
 * 体中的每条边由且仅由两个面片共享。
 * 算法：决定面片对
 * -------------------------------------------------------------------*/
int judge_pair_from_four_facet(
	int edgeS,  int edgeE,	 /* 公共边的首末点 */
	int othe[], bool orit[], /* othe: 四个面片中除首末点外的另外一个点，
							  * orit: 边在面中的绕向 */
	int retn[]
	);

/* ----------------------------------------------------------------------
 * 根据单元信息填充一个InsertedTetra对象
 * -------------------------------------------------------------------*/
inline int fill_inserted_tetra(
			int ELE[], int f,			/* 节点和基面 */
			InsertedTetra *instTetra,
			polyhedron *poly,
			polyhedron *sub_poly);

/* 确定4个面的信息, 保证索引最小的点是第1个点，且保证每个面的法向指向区域内部 */
inline void fill_inst_tetra_facets(int ELE[], InsertedTetra *instTetra);

/* ----------------------------------------------------------------------
 * 对复边进行排序
 * -------------------------------------------------------------------*/
int sort_comp_edges(struct polyhedron*poly, GridEdgeHash *eHash);

int sort_quad_pair_test();
int four_facet_pair_test(struct polyhedron* poly);
int assert_poly_neig(struct polyhedron *poly);
int assert_simp_poly_edge(struct polyhedron *simp_poly);

/* ----------------------------------------------------------------------
 * 对中间结果进行记录
 * -------------------------------------------------------------------*/
inline int is_nonoverlap_tet(divide_record *oldRd, int key, int *pos);
inline int create_div_record(divide_record *oldRd, int key, divide_record *newRd, int *found);
inline void binarySearch(int * loc, int * found, int Item, int Array[], int right , int left);
inline int is_tested_div_record(divide_record *cur_record);
inline void add_div_record(divide_record *cur_record);

//inline int sort_TRIp_pair(int *j1, int *j2, int *j3);
//inline int make_TRIp_pair(int *j1, int *j2, int *j3);

extern std::set<divide_record, RuntimeCmp<divide_record> > setDivideRecord;	/* 记录节点在环内的编号 */

/* ----------------------------------------------------------------------
 * 对面片和节点进行排序
 * -------------------------------------------------------------------*/
inline void sort_facets(polyhedron *poly);
inline int get_sorted_facet(polyhedron *poly, int fIter);
inline void sort_nodes(polyhedron *poly, int fDig);
inline int get_sorted_node(polyhedron *poly, int nIter);

inline int calc_worst_facet(polyhedron *poly);

/* ----------------------------------------------------------------------
 * 为防止面片和复边但非对面上的点形成对应关系，排除这些点
 * -------------------------------------------------------------------*/
inline int is_mutex_facet_node(polyhedron *poly, int fIdx, int nIdx);
inline int get_mutex_facet_node(polyhedron *poly, int fIdx, std::set<int>& setMutexNodes);
inline int get_around_nodes(polyhedron *poly, int fIdx, std::set<int>& setAroundNodes);
inline int get_neighb_nodes(polyhedron *poly, int fIdx, std::set<int>& setAroundNodes);


/* ---------------------------------------------------------------------------------------
 * 分两步做递归
 * 1. 记录当前子问题的所有可能分解方案（增加1个四面体单元情形下）
 * 2. 启用某种寻优策略来确定优先的分解方案
 * 由于现有横向的比较，而不是盲目地递归，有可能会获得更好的求解效率
 * --------------------------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
inline int optimal_twoStep_firstCall(double iniQ, std::vector<polyhedron *>& vec_simp_poly, divide *dive, divide_record *divR,
					  Constraints *csts, int maxIntNum = 0);
inline int optimal_twoStep(double iniQ, std::vector<polyhedron *>& vec_simp_poly, divide *dive, divide_record *divR, 
					  Constraints *csts = NULL, int maxIntNum = 0);
#else
inline int optimal_twoStep_firstCall(double iniQ, polyhedron *simpPolyArray[], int simpPolyNum, divide *dive, divide_record *divR,
					  Constraints *csts, int maxIntNum = 0);
inline int optimal_twoStep(double iniQ,  polyhedron *simpPolyArray[], int simpPolyNum, divide *dive, divide_record *divR, 
					  Constraints *csts = NULL, int maxIntNum = 0);
#endif /* _SPR_USING_STD_VECTOR */

inline int optimal_twoStep(double iniQ, polyhedron *poly, divide *dive, divide_record *divR, 
					  Constraints *csts = NULL, int maxIntNum = 0);
inline int optimal_fourFace(double iniQ, polyhedron *poly, divide *dive, divide_record *divR, Constraints *csts = NULL, int maxIntNum = 0);

#ifdef _SPR_USING_STD_VECTOR
inline int sort_subproblems(std::vector<SubDivideProblem*>& vec_sub_prob);
#else
inline int sort_subproblems(SubDivideProblem *subProbArray[], int subProbNum);
#endif /* _SPR_USING_STD_VECTOR */

/* ---------------------------------------------------------------------------------------
 * 分两步做递归
 * 1. 记录当前子问题的所有可能分解方案（增加1个四面体单元情形下）
 * --------------------------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
inline int record_subproblems(double iniQ, struct polyhedron *poly, 
							  divide_record *divR,
							  std::vector<SubDivideProblem*>& vec_sub_prob, 
							  Constraints *csts = NULL, int maxIntNum = 0);
inline int record_subproblems(double iniQ, struct polyhedron *poly, int f, int maxSize,
							  divide_record *divR,
							  std::vector<SubDivideProblem*>& vec_sub_prob, 
							  Constraints *csts = NULL, int maxIntNum = 0);
#else
inline int record_subproblems(double iniQ, struct polyhedron *poly, 
							  divide_record *divR,
							  SubDivideProblem *subProbArray[], int *subProbNum, 
							  Constraints *csts = NULL, int maxIntNum = 0);
inline int record_subproblems(double iniQ, struct polyhedron *poly, int f, int maxSize,
							  divide_record *divR,
							 SubDivideProblem *subProbArray[], int *subProbNum,
							  Constraints *csts = NULL, int maxIntNum = 0);
#endif /* _SPR_USING_STD_VECTOR */
							 
/* 计算相对值 */
inline double calc_facet_quality_rel(polyhedron *poly, int fIdx, polyhedron *pld_poly, double mindst); //zhaodawei

/* 计算相对值 */
inline double calc_facet_quality_rel(polyhedron *poly, int fIdx);
/* 计算绝对值 */
inline double calc_facet_quality(polyhedron *poly, int fIdx);
inline int compare_idx_quality(const void *arg1, const void *arg2 );

/* 初始化 */
int init_spr();

int small_poly_reconn(polyhedron *poly, divide *divi, double iniQ, 
					  Constraints *csts = NULL, int *smallPolySize = NULL, int maxIntNum = 0);

int small_ball_reconn(polyhedron *poly, divide *divi, double iniQ, int *iAdjNode);
int small_ball_reconn(polyhedron *poly, divide *divi, 
					  double minQ_forAll, double minQ_forBnd, 
					  int *nodeFlags, int *iAdjNode);
bool is_valid_poly(struct polyhedron *poly);
int remove_overlap_faces(struct polyhedron *poly);

#endif /* #ifndef __spr_h__ */
