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
typedef double APOINT[3];	/* ���������POINT�������POINT��������һ��double����ά����*/
typedef int TRI[3];			//����һ�������Σ����ӵ���һ��Ԫ���Ǳ�������������������ĸ��������������
typedef int NEG[3];			/* �������ڹ�ϵ */ 
typedef int TETRA[4];		//����һ��������

#define NOD 1	
#define EDG 2   
#define FAC 3
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define ANGLE2RADIO(x) (x*PI)/180
#define RADIO2ANGLE(x) (x*180.0)/PI
#define SPR_ZERO_QUALITY 1.0e-16			/* default: 1.0e-7С�����ֵ�ĵ�Ԫ�������ض�Ϊ��ֵ(����Ч�������Ч�ʣ�����Ч������������) */
#define SPR_UNALLOWED_QUALITY 0.0//1.0e-18		/* ȡ����ֵ0������һЩ���⣬�������ر���������� 0.0			/* С�����ֵ�ĵ�Ԫ����Ϊ0�����Ԫ */
#define SPR_EPS_ZERO 1e-8
#define SPR_EPS_ZERO_SQ 1.e-8
#define MAX_SPR_POLY_SIZE 400//default is 300
#define MAX_CONSTRAINT_SIZE 300
#define MAX_COPLANE_NODES 512 // = MAX_ADD_NODE (defined in iso3d_define.h)
#define MAX_ALLOWED_POLY_SIZE 41
#define MAX_SIMP_POLY_NUM	32 /* �����Ӷ������Ŀ */	
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

extern int count;		//��������Ч�Ե��жϴ���
extern int real;		//ʵ�ʵĵݹ����
extern int box_n;		//��box�����������ཻ�Ĵ���
extern int face_n;		//�ܹ��������ཻ���жϴ���
extern int line_n;		//�ܹ��������ཻ���жϴ���
extern int box_l_n;	//��box�����������ཻ�Ĵ���
extern int float_count;	//�ܹ����ø�������жϵ������ཻ�Ĵ���
extern long orient3d_n;		//��������Ĵ���
extern long Call_V;	//��ͼ���Ҵ���
extern long Call_V_R;	//����֪�����е��õĴ���
extern long Call[10];

extern double start, stop;	//��¼����ʱ��
extern double duration;

extern double starttmp[10], stoptmp[10];	//��¼ϸ��ʱ��
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
extern double spr_time_limit;		/* ���β�����ʱ������(ָ���������壬һ��������ܰ�����������壬time_limit = 0.0��ʾ������) */
extern int spr_optimal_type;
extern int spr_max_int_node;
extern bool spr_verbose_output;
extern bool spr_no_more_recur;

/* for test only */
extern int spr_l2g[MAX_SPR_POLY_SIZE];

typedef struct polyhedron
{
        TRI t[MAX_SPR_POLY_SIZE];			/* ��Ƭ���γɽڵ� */
	NEG n[MAX_SPR_POLY_SIZE];			/* ��Ƭ�����ڹ�ϵ */
	int used_v[MAX_SPR_POLY_SIZE];		/* �����ڵ��������Ƭ��Ŀ��Ϊ0ʱ����ʾ��Ӧ�ڵ㲻�ڵ�ǰ�������� */
	int vprt[MAX_SPR_POLY_SIZE];		/* ǰ2λΪ���׽ڵ㣬ָʾ�ڵ��ڵ�Ԫ�еı��룬�����ڵ��һ����Ԫ */
	int surface[MAX_SPR_POLY_SIZE];		/* ��Ƭ���ڵ�����Ƭ��ţ�ʵ��û���õ� */
	int num_t;							/* ��Ƭ��Ŀ */
	int num_vable;						/* used_vֵ����0�Ķ������ */
	int sort_f[MAX_SPR_POLY_SIZE];		/* ����Ƭ��������Ľ��(���ڼ���SPR����) */ 
	int sort_n[MAX_SPR_POLY_SIZE];		/* �Զ����������Ľ��(���ڼ���SPR����)��ʵ��û���õ� */
	int newf_fg[MAX_SPR_POLY_SIZE];		/* �����ӵ����flag 0���������ӵ��棻1�������ӵ��� */
	int newf_id[3];						/* �����ӵ��棬������3�� */
	int newf_ct;						/* �����������Ŀ */
	int allowedFaceNum;					/* ��������ߵ���Ƭ���ֵ����spr_optimal_type = EDG_REMOVALʱ���� */ 
} polyhedron;

typedef struct divide
{
	TETRA te[MAX_SPR_POLY_SIZE];
	int num;
	double q;
	/* Ϊͳ��ʵ�ʵĽ����������Ҫ��ƹ�ϣ�� */
	int intEdges[MAX_SPR_INT_EDGE][2], intFaces[MAX_SPR_INT_FACE][3];
	int intEdgeNum, intFaceNum;
	int intNum;				/* ʵ�ʵĽ����������spr_optimal_type = MIN_INT_NODE_NUMBERʱ���� */
	int facNum;				/* �����ߵ���Ƭ��������spr_optimal_type = EDG_REMOVALʱ���� */ 
} divide;

typedef struct divide_record
{
	int num;
	int div_rd[MAX_SPR_POLY_SIZE];
} divide_record;

typedef struct SubDivideProblem
{
	polyhedron *poly;	/* ���׶����� */
	int tetraID;		/* ��������, �����ĸ��ڵ���� */
	int fIdx;			/* ��Ƭ���  */
	int nIdx;			/* �ڵ���  */
	double initQ;		/* ������ֵ  */
	double tetrQ;		/* �����嵥Ԫ������ */
	double diviQ;		/* �ֽ���������������� */
#ifdef _SPR_USING_STD_VECTOR
	std::vector<polyhedron*> vec_simp_poly;	/* ����������Ķ����� */
#else
	polyhedron *simpPolyArray[MAX_SIMP_POLY_NUM];
	int simpPolyNum;
#endif
	divide_record	divR;					/* ���ڼ�¼��ǰ�ֽ�������Ϣ */
	int intNum;			/* ����ĸ���(��spr_optimal_type = MIN_INT_NODE_NUMBERʱ��Ч) */
	int facNum;			/* �����ߵ���Ƭ��������spr_optimal_type = EDG_REMOVALʱ���� */ 
} SubDivideProblem;

/* ��SPR�㷨����Steiner��ɾ��ʱ����Ҫ����һЩԼ�� */
typedef struct Constraints
{
	int lineCsts[MAX_CONSTRAINT_SIZE][2];		/* ��Լ������ĩ�˵� */
	int lineCst_size;							/* ��Լ���ĸ��� */
	int faceCsts[MAX_CONSTRAINT_SIZE][3];		/* ��Լ����3���˵� */
	int faceCst_size;							/* ��Լ���ĸ��� */
	int colineNodes[MAX_COPLANE_NODES];			/* ���߽ڵ��б� */
	int colineN_size;							/* ���߽ڵ���� */
	int cofaceNodes[MAX_COPLANE_NODES];			/* ����Լ������ͻ�Ĺ���ڵ��б���Լ�����ڵ������������Steiner���б�*/
	int cofaceN_size;							/* ����ڵ�ĸ��� */
	int badVolmCsts[MAX_CONSTRAINT_SIZE][4];	/* ��ϣ�����������е������嵥Ԫ */
	int badVolmCst_size;						/* ��Լ���ĸ��� */
} Constraints; 

extern Constraints spr_csts;

bool operator == (const divide_record & t1, const divide_record& t2);

bool operator < (const divide_record & t1, const divide_record& t2);

bool operator > (const divide_record & t1, const divide_record& t2);

/* -----------------------------------------------------------------------------
 * ʵ��set�������ʱ�ȽϹ���
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
        TRI t[4];			/* 4����Ƭ */
	NEG n[4];			/* 4����Ƭ֮������ڹ�ϵ */
	int idx[4];			/* 4����Ƭ�ı�ţ�С��0��ζ�Żᱻ����, ����0�ᱻɾ�� */
	int num_del;		/* 4����Ƭ���ж��ٸ���Ƭ�ᱻɾ�� */
	int flg[MAX_SPR_POLY_SIZE];		/* ȷ����Щ���Ƿ�ɾ���ı�־ -1: ����ɾ��, 0~2:��ɾ����ľֲ����� */
	bool complex;		/* �Ƿ�������� */
} InsertedTetra;
 
typedef struct IDX_QUALITY
{
	int idx;
	double quality;
} IDX_QUALITY; 

/* ----------------------------------------------------------------------------
 * Ϊ��ǿ�߽�ָ��㷨��������Ҫ��ǿSPR�㷨
 * --------------------------------------------------------------------------*/
enum EDG_CONS_TYPE {MUST_EXIST = 0, MUST_BE_EXCLUDED, DEPEND_ON_LIST};
typedef struct EdgeConstraint
{
	int i1, i2;						   /* �����ڵ�2���˵� */
	int j1, j2, j3;					   /* ���������ı߻��棬���i3 = -1��ʾ��һ���� */
	EDG_CONS_TYPE type;				   /* ���� */
	InstEntList *pListAllowed;		   /* �����ཻ��ʵ�� */
	InstEntList *pListExcluded;		   /* �������ཻ��ʵ�� */
} EdgeConstraint;
extern EdgeConstraint edgeConsArray[MAX_CONSTRAINT_SIZE];
extern int edgeConsArraySize;
extern int spr_shell_edge_nodeS, spr_shell_edge_nodeE;

/* spr_optimal_type = EDGE_REMOVALʱ����ʼ�ǰ����ߵ���Ŀ */
extern int spr_shell_size;
extern int spr_allowed_face_num;	/* ��¼�����������������Ƭ�����ֵ����ʱ���� */
extern int spr_subprob_face_num;	/* ��¼������ʵ�ʰ�������Ƭ������ʱ���� */
extern int spr_dependent_edges[MAX_CONSTRAINT_SIZE][3]; /* ��3����ԱΪ���͡�=0��û�ж���Ϊ�Ǳ߶��㣻=1����һ������Ϊ�Ǳ߶˵㣻*/
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
/*�������ln1��ln2���ߵ������˵㣬fac1��fac2��fac3���Ǹ������������㣬d�������ߵĳ��ȣ�areaa�������������
ͨ���������Է��ص�ֵ�����pnt�����ص�������ߺ����ཻ����ô����������val���ཻ�ľ������
w1��w2��w3������ָʲô��
�ں����У�������һ��NOD��1����ʱ�򣬱��������ߺ�������ཻ������������ߴ�������Ķ˵�
������һ��EDG��2����ʱ�򣬱��������ߺ�������ཻ������������ߴ�������ı߽���
������һ��FAC��3����ʱ�򣬱��������ߺ������ֱ���ཻ
������-1��ʱ�򣬱��������ߺ�����治�ཻ
������-2��ʱ�򣬱�������������ش��󣬵��Ǿ��廹�ǲ���*/
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
 * �������Ӷ���������ڹ�ϵ�����GridEdgeHash���ڸ��ߣ��ٶ����߶�Ӧ��
 * ��Ԫ�������
 * -------------------------------------------------------------------*/
int build_poly_neig(struct polyhedron *poly, GridEdgeHash *eHash);

int build_manifold_poly_neig(struct polyhedron *poly);

/* ----------------------------------------------------------------------
 * �ֽ⸴�Ӷ�����Ϊ����򵥶�����
 * �������ڹ�ϵ����ɫ�㷨
 * -------------------------------------------------------------------*/
#ifdef _SPR_USING_STD_VECTOR
int divide_poly(struct polyhedron *cmpx_poly, 
				std::vector<struct polyhedron*>& vec_simp_poly);
#else
int divide_poly(struct polyhedron *cmpx_poly, 
				struct polyhedron *simpPolyArray[], int *simpPolyNum);
#endif /* _SPR_USING_STD_VECTOR */

/* ----------------------------------------------------------------------
 * ����׸�����ĳ�ڵ�ĸ��׵�Ԫ
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
 * ������а���ĳ�ڵ�ĸ��׵�Ԫ
 * -------------------------------------------------------------------*/
inline void get_node_parents(struct polyhedron *poly, int nidx, 
				std::vector<int>& vecp);

inline void get_node_parents(struct polyhedron *poly, int nidx, int *vecp, int *num);

/* ----------------------------------------------------------------------
 * ������а���ĳ���ߵ����и��׵�Ԫ
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
 * ��֪������һ�����Ӷ����壬ÿ���涼���з���ģ��ҷ���ָ�������ڲ�����ĳ
 * Щ����4���湲����Щ���γ������棬ÿ���������һ���򵥶����塣�򵥶���
 * ���е�ÿ�������ҽ���������Ƭ����
 * �㷨��������Ƭ��
 * -------------------------------------------------------------------*/
int judge_pair_from_four_facet(
	int edgeS,  int edgeE,	 /* �����ߵ���ĩ�� */
	int othe[], bool orit[], /* othe: �ĸ���Ƭ�г���ĩ���������һ���㣬
							  * orit: �������е����� */
	int retn[]
	);

/* ----------------------------------------------------------------------
 * ���ݵ�Ԫ��Ϣ���һ��InsertedTetra����
 * -------------------------------------------------------------------*/
inline int fill_inserted_tetra(
			int ELE[], int f,			/* �ڵ�ͻ��� */
			InsertedTetra *instTetra,
			polyhedron *poly,
			polyhedron *sub_poly);

/* ȷ��4�������Ϣ, ��֤������С�ĵ��ǵ�1���㣬�ұ�֤ÿ����ķ���ָ�������ڲ� */
inline void fill_inst_tetra_facets(int ELE[], InsertedTetra *instTetra);

/* ----------------------------------------------------------------------
 * �Ը��߽�������
 * -------------------------------------------------------------------*/
int sort_comp_edges(struct polyhedron*poly, GridEdgeHash *eHash);

int sort_quad_pair_test();
int four_facet_pair_test(struct polyhedron* poly);
int assert_poly_neig(struct polyhedron *poly);
int assert_simp_poly_edge(struct polyhedron *simp_poly);

/* ----------------------------------------------------------------------
 * ���м������м�¼
 * -------------------------------------------------------------------*/
inline int is_nonoverlap_tet(divide_record *oldRd, int key, int *pos);
inline int create_div_record(divide_record *oldRd, int key, divide_record *newRd, int *found);
inline void binarySearch(int * loc, int * found, int Item, int Array[], int right , int left);
inline int is_tested_div_record(divide_record *cur_record);
inline void add_div_record(divide_record *cur_record);

//inline int sort_TRIp_pair(int *j1, int *j2, int *j3);
//inline int make_TRIp_pair(int *j1, int *j2, int *j3);

extern std::set<divide_record, RuntimeCmp<divide_record> > setDivideRecord;	/* ��¼�ڵ��ڻ��ڵı�� */

/* ----------------------------------------------------------------------
 * ����Ƭ�ͽڵ��������
 * -------------------------------------------------------------------*/
inline void sort_facets(polyhedron *poly);
inline int get_sorted_facet(polyhedron *poly, int fIter);
inline void sort_nodes(polyhedron *poly, int fDig);
inline int get_sorted_node(polyhedron *poly, int nIter);

inline int calc_worst_facet(polyhedron *poly);

/* ----------------------------------------------------------------------
 * Ϊ��ֹ��Ƭ�͸��ߵ��Ƕ����ϵĵ��γɶ�Ӧ��ϵ���ų���Щ��
 * -------------------------------------------------------------------*/
inline int is_mutex_facet_node(polyhedron *poly, int fIdx, int nIdx);
inline int get_mutex_facet_node(polyhedron *poly, int fIdx, std::set<int>& setMutexNodes);
inline int get_around_nodes(polyhedron *poly, int fIdx, std::set<int>& setAroundNodes);
inline int get_neighb_nodes(polyhedron *poly, int fIdx, std::set<int>& setAroundNodes);


/* ---------------------------------------------------------------------------------------
 * ���������ݹ�
 * 1. ��¼��ǰ����������п��ֽܷⷽ��������1�������嵥Ԫ�����£�
 * 2. ����ĳ��Ѱ�Ų�����ȷ�����ȵķֽⷽ��
 * �������к���ıȽϣ�������äĿ�صݹ飬�п��ܻ��ø��õ����Ч��
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
 * ���������ݹ�
 * 1. ��¼��ǰ����������п��ֽܷⷽ��������1�������嵥Ԫ�����£�
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
							 
/* �������ֵ */
inline double calc_facet_quality_rel(polyhedron *poly, int fIdx, polyhedron *pld_poly, double mindst); //zhaodawei

/* �������ֵ */
inline double calc_facet_quality_rel(polyhedron *poly, int fIdx);
/* �������ֵ */
inline double calc_facet_quality(polyhedron *poly, int fIdx);
inline int compare_idx_quality(const void *arg1, const void *arg2 );

/* ��ʼ�� */
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
