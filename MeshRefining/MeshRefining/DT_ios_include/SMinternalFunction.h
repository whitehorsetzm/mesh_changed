#ifndef SM_INTERNAL_FUNC_H
#define SM_INTERNAL_FUNC_H 1

#include "SM_config.h"
/* Function Declarations */

/* INTERNAL FUNCTIONS */
/* routines for mallocing and initializing space */
SMlocal_mesh    *SMmallocLocalMesh(int num_pts, SMparam *smooth_param);
SMlap_info      *SMmallocLap(int num_values);
SMoptimal       *SMmallocOpt(int num_tri);
SMlp            *SMmallocLP(int num_active, int num_tri);
SMactive        *SMmallocActive(int num_values);
SMquality_table *SMmallocQualityTable(void);
void             SMinitLocalMesh(int num_pts, int num_tet, 
                                 double *free_vtx, double **vtx_list,
                                 int **vtx_connectivity, SMlocal_mesh *local_mesh, 
                                 SMparam *smooth_param);
void             SMinitSmoothParam(int technique, int FunctionID,
                                   double AcceptFunction, void *ext_smooth_data);  
SMstats         *SMinitStats(SMstats *smooth_stats);
SMprocinfo      *SMinitProcinfo(void);
void             SMinitLap(int num_values, SMlap_info *lap_info);
void             SMinitOpt(int num_values, int maxit, SMoptimal *opt_info);
void             SMinitLP(SMlocal_mesh *local_mesh);

/* function-gradient routines -- 2D */
void             SMcomputeFunction(SMlocal_mesh *local_mesh, 
                         SMparam *smooth_param, double *function);
void             SMcomputeGradient(SMlocal_mesh *local_mesh, 
                         SMparam *smooth_param, double **gradient);

void             SMcomputeTriCosines(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeCosGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2],int *num_values);

void             SMcomputeInteriorTriCosines(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeInteriorCosGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2],int *num_values);

void             SMcomputeNegTriCosines(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeNegCosGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeTriAngles(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeAngGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2],int *num_values);

void             SMcomputeInteriorTriAngles(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeInteriorAngGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeNegTriAngles(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeNegAngGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2],int *num_values);

void             SMcomputeTriSines(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeSineGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2],int *num_values);

void             SMcomputeInteriorTriSines(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeInteriorSineGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeTriJacobians(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeJacobianGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeScaledTriJacobians(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeScaledJacobianGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeInteriorScaledTriJacobians(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeInteriorScaledJacobianGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeAreaLengthRatio(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeAreaLengthRatioGradients(double *vtx1, double *vtx2, 
                         double *vtx3, double (*gradient)[2], int *num_values);

void             SMcomputeLengthAreaRatio(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);
void             SMcomputeTriArea(double *vtx1, double *vtx2, 
                         double *vtx3, double *function, int *num_values);

void SMNormJacSquared2D(double *vtx1, double *vtx2, double *vtx3, 
		       double *function, int *num_values);
void SMcomputeNormJacSquaredGradients2D(double *vtx1, double *vtx2, 
		double *vtx3, double (*gradient)[2], int *num_values);

void             SMcopyFunction(int num_values, double *function1, double *function2);
void             SMcopyActive(SMactive *active1, SMactive *active2);

/* function - gradient routines 3D */
#include "SMderiv.h"
#include "SMintrinsic.h"
#include "SMdihed_func.h"

/* New initial point routines */
void             SMlaplaceSmooth(SMlocal_mesh *local_mesh,
                                 SMparam *smooth_param, SMprocinfo *procinfo);
void             SMcalcNewInitPt(SMlocal_mesh *local_mesh, 
                                 SMparam *smooth_param, SMprocinfo *procinfo);
void             SMcentroidSmoothMesh(int num_incident,double **incident_vtx,
				      double *free_vtx, int dimension,
				      const int is_surface_vertex);

/* optimization routines */
void             SMfindActiveSet(int num_values, double *function, 
                            double active_eps, SMactive *active_info);
double           **SMgetActiveDirections(int num_active, double **gradient, 
                                  int *active_ind, int dimension);
void             SMsearchDirection(SMlocal_mesh *local_mesh);
void             SMsearchEdgesFaces(int num_active, double **G, double **dir, 
                                    SMoptimal *opt_info);
void             SMprintActiveSet(SMactive *active, double *function);
void             SMstepAcceptance(SMlocal_mesh *local_mesh, 
                                  SMparam *smooth_param);
void             SMstepToCusp(SMlocal_mesh *local_mesh,SMparam *smooth_param,
                              SMprocinfo *procinfo);
void             SMcomputeVerticalStep(SMlocal_mesh *local_mesh,
                              SMparam *smooth_param);
double           *SMformVerticalRHS(int n,double *function, int ind[MAX_DIM], 
                                   double max_value);

void             SMgetGradientProjections(SMoptimal *opt_info);
double           SMgetMinEstimate(SMoptimal *opt_info);
int              SMcheckEquilibrium(SMoptimal *opt_info);
void             SMcomputeAlpha(SMoptimal *opt_info);

int              SMconvexHullTest(double **vec, int num_vec);
void             SMfindPlaneNormal(double pt1[3], double pt2[3], 
                                    double pt3[3],double *cross);
int              SMcheckVectorDots(double **vec,int num_vec,double *normal);
void             SMfindPlanePoints(int dir1, int dir2, double **vec, 
                       int num_vec, double *pt1, 
                       double *pt2, double *pt3, int *status);

/* untangling routines */
int    SMuntangle_mesh(SMlocal_mesh *local_mesh);
int    SMdegenerate(int num_active, int num_constraints, double **A, double *b);
void SMcomputeConstraintMatrix(SMlocal_mesh *local_mesh, int num_constriants,
                                                        double **Amat, double *b);
int    SMphaseOneLP(int num_constraints, int num_active, double **A, double *x, SMlp *lp_info);
int    SMsolveLP(int num_constraints, int num_active, double **A, double *x, double *c, double *pi,
                           SMlp *lp_info);
void SMgetActiveMatrix(double **AAmat,int num_active,int *active_ind, double *Bmat);
void SMgetActiveRHS(double *c,int num_active,int *active_ind, double *c_active);
int SMremoveIdenticalVtx(int dimension, int *num_incident_vtx,int *num_tri, 
                         double ***vtx_list, int ***vtx_connectivity);

/* stats routines */
void             SMregisterEvents(void);
void             SMaccumulateStats(SMlocal_mesh *local_mesh,
                               SMparam *smooth_param, 
		               SMstats *smooth_stats);
void             SMprintStats(SMsmooth_data *);
void             SMwriteStatsFile(SMstats *smooth_stats, int smooth_count);

/* matrix routines */
double      *SMsolve2x2(double a11, double a12, double a21, double a22, 
		             double b1, double b2);
void         SMsolve2x2_t(double a11, double a12, double a21, double a22, 
		             double b1, double b2, double *x);
double      *SMsolveSPD3x3(double **A, double *B);
void         SMsolve3x3(double **A, double *B, double *x);
int          SMformGrammian(int num_vecs, double **vecs, double **G, int dimension);
void         SMformPDGrammian(SMoptimal *opt_info);
double     **SMformReducedMatrix(int num_active, double **G);
double     **SMformVerticalMatrix(int num_active, double **PDG);
int          SMnonSingularTest(int n, double **A);
void         SMtransposeMatrix(double *mat, int n, int m, double *mat_T);
void         SMtransposeMatrix2(double **mat, int n, int m, double **mat_T);
double       SMcondition3x3(double **A);
double       SMdeterminant3x3(double a1[3], double a2[3], double a3[3]);
void         SMmultiply3x3(double a1[3], double a2[3], double a3[3],
                   double b1[3], double b2[3], double b3[3],
		   double r1[3], double r2[3], double r3[3]);
void SMtranspose3x3(double a1[3], double a2[3], double a3[3],
		    double b1[3], double b2[3], double b3[3]);
void SMadjoint3x3(double a1[3], double a2[3], double a3[3],
                  double b1[3], double b2[3], double b3[3]);
double SMfrobenius_norm_squared3x3(double a1[3], double a2[3], double a3[3] );
double SMfrobenius_norm_squared_adj3x3(double a1[3], double a2[3], double a3[3] );
double SMfrobenius_norm_squared2x2(double a1[2], double a2[2]);


/* miscellaneous routines */
void             SMrecordIterValues(SMoptimal *opt_info, double *pt);
void             SMinsertQualityInfo(SMquality_table *quality_table, int measure_id, 
                                      double *function, int num_values);
int              SMvalidMesh(SMlocal_mesh *local_mesh);
int SMvalidityCheck(SMlocal_mesh *local_mesh);
int SMorient2D(double *vtx1, double *vtx2, double *vtx3);
int SMorient3D(double *vtx1, double *vtx2, double *vtx3, double *free_vtx); 

/* write to matlab file routines */
void             SMwriteLocalMesh(FILE *fp, SMlocal_mesh *local_mesh);
void             SMwriteLocalAxes(FILE *fp, SMlocal_mesh *mesh);
void             SMwriteLocalTriangleList(FILE *fp, SMlocal_mesh *local_mesh);
void             SMwriteActiveSet(FILE *fp,SMlocal_mesh *local_mesh);
void             SMwritePoint(FILE *fp, double x, double y);
void             SMwriteSearch(FILE *fp, SMlocal_mesh *local_mesh);

/* free routines */
void             SMfreeOptimal(SMlocal_mesh *local_mesh);
void             SMfreeLP(SMlocal_mesh *local_mesh, int num_active,
			  int num_constraints);
void             SMfreeActive(SMactive *active);
void             SMfreeParam(SMparam *smooth_param);
void             SMfreeLocalMesh(SMlocal_mesh *local_mesh);
void             SMfreeProcinfo(SMprocinfo *procinfo);
void             SMfreeQualityTable(SMquality_table *quality_table);
void             SMfreeActiveDirections(const int num_active, double **dir);

/* error routines */
int SMwrite_ordered_points(SMlocal_mesh *local_mesh);

/* stuff for surface smoothing */
void             SMunitNormal(const double location[MAX_DIM],
			      double normal[MAX_DIM]);  


/* lapack stuff */
#ifdef rs6000
#define DGESV dgesv
#else
#define DGESV dgesv_
#endif

void dgesv_(int *n, int *nrhs, double *A, int *lda, int *IPIV, double *B,
	    int *LDB, int *INFO);
#endif
