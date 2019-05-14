#ifndef SM_EXTERNAL_FUNC_H
#define SM_EXTERNAL_FUNC_H 1

#include "SMuserDefs.h"
#include "SMerror.h"

/* EXTERNAL FUNCTIONS */
/* functions needed in the user code */

#ifdef __cplusplus
extern "C" {
#endif
  /* Initialization Routines */
void             SMinitSmoothing(int dimension,
                     int technique, int FunctionID, double AcceptFunction, 
                     void **smooth_data);
void             SMsetProblemDimension(void *smooth_data, int dimension);
void             SMsetSmoothTechnique(void *smooth_data, int technique);
void             SMinitGlobalMinValue(void *smooth_data);
void             SMsetSmoothThreshold(void *smooth_data, double accept);
void             SMsetSmoothFunction(void *smooth_data, int FunctionID);
void             SMinitLogging(int argc, char** argv);
  
  /* Mesh Improvement Routines */
int              SMsmooth(int num_pts, int num_tet, double *free_vtx, 
                          double **vtx_list, int **vtx_connectivity,
                          void *smooth_data, int is_surface_vertex);

int              SMuntangle(int num_pts, int num_tet, double *free_vtx, 
                          double **vtx_list, int **vtx_connectivity,
                          void *smooth_data);

/* void             SMsetUserQualityFunction2D(void *ext_smooth_data, 
                                          int values_per_tri,
                                          SMfunction_ptr2D userQualityFunc,
                                          SMgradfunc_ptr2D userQualityGrad); 
*/

  /* Statistics Routines */
void             SMinitSmoothStats(void *smooth_data);
void             SMprintSmoothStats(void *smooth_data);
void             SMconvertToDegrees(int,double *);

  /* Quality Routines */
void             SMinitQualityTable(void *smooth_data);
void             SMaccumulateQualityInformation(void *smooth_data, double **vtx); 
void             SMprintQualityInformation(void *smooth_data);
int                SMinvalidMesh(void *smooth_data);

/* Error Routines - Accessed through CHKERR(ierr) */
int              SMerror(int line,char *func,char* file,char *dir,int n,int p,char *mess);

  /* Clean-up routines */
void             SMfinalizeSmoothing(void *smooth_data);
void             SMfinalizeLogging(void *smooth_data);
  
  /* Routines for surface smoothing */
void             SMsetNormal(const double new_norm[]);

  /* User customization */
void SMsetMeshValidity(int mesh_validity, void *ext_smooth_data);
  
#ifdef __cplusplus
}
#endif

#endif 
