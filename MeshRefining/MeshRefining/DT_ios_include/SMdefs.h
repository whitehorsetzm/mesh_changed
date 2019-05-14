#ifndef SM_DEFS_H
#define SM_DEFS_H 1

#include "SM_config.h"

//#define MAX_NUM_PTS 150
#define MAX_NUM_PTS 500
#define MAX_DIM 3

#define XDIR 0
#define YDIR 1
#define ZDIR 2

#define TRUE 1
#define FALSE 0

/* various constants to make the code more readable */
#define CCW               1
#define CW                  0
#define A_BIG_POS_NUMBER    1E300
#define A_BIG_NEG_NUMBER -1E300
#define VALID_MESH                 1
#define INVALID_MESH             0
#define MAX_SM_ITER               20
#define MAX_SM_INTS               50
#define END_PT_TOL                 1E-10
#define SINGULAR_SYSTEM_TEST 1E-14
#define MACHINE_EPS               1E-15
#define LESS_THAN_MACHINE_EPS(x)   ( ((fabs(x)+1.0) > 1.0) ? 0 : 1 )
#define INIT_MIN_VALUE         1E300
#define MAX_G_NUM                2500		//150

#ifndef WIN_VC32
#include <unistd.h>
#else
#include <io.h>
#endif
#include <fcntl.h>
#include <stdlib.h>

#if defined(PARCH_sun4) && !defined(__cplusplus) && defined(_Gnu_)
    extern int  open(const char *, int, ...);
    extern int  creat(const char *, unsigned short);
    extern int  write(int, const void *, unsigned int); 
    extern int  close(int);
    extern int  fprintf(FILE*,const char*,...);
#endif
#if defined(solaris) && !defined(__cplusplus) && defined(_Gnu_)
    extern int  open(const char *, int, ...);
    extern int  write(int, const void *, unsigned int); 
    extern int  close(int);
/*    extern int  fprintf(FILE*,const char*,...);*/
    extern void *malloc(long unsigned int);
    extern void free(void *);
#endif

/* the debugging system
     Level 0 provides no information and the debug macros are empty
     Level 1 provides user function information only
             e.g. threshold set, function used, etc.
     Level 2 provides basic algorithmic information for each local
             submesh
     Level 3 provides more information, data structures, and details
             than most users would want to know about 

     The default is Level 0
*/

#if SM_DEBUG_LEVEL != 0
#define SM_DEBUG_PRINT(level, statement)\
{\
   if ((level <= SM_DEBUG_LEVEL))\
     {\
     fprintf(stdout,statement);\
     fflush(stdout);\
     }\
}
#define SM_DEBUG_ACTION(level, statement)\
{\
   if ((level <= SM_DEBUG_LEVEL))\
     {\
     statement\
     fflush(stdout);\
     }\
}
#else
#define SM_DEBUG_PRINT(level, statement)\
{\
}
#define SM_DEBUG_ACTION(level, statement)\
{\
}
#endif

#define SM_ERROR(msg,action)\
{ \
   fprintf(stderr,msg);\
   action\
}

/* This should be used with debugging level 3.  It writes out a
series of matlab files that illustrate the initial local submesh,
the search direction at each step, and the final local submesh */

#if defined(LOCALTEST)
#define MATLAB_ON(a) a
#else
#define MATLAB_ON(a)
#endif
#define MATLAB_OFF

/* The enables statistics gathering;  the user must initialize the
statistics datastructures in the main code and probably wants to print
them out after every smoothing pass 
   The information provided with stats include the number of grid points
   smoothed, how many required optimizaiton, the avg number of optimization
   steps, the a breakdown of the optimization termination criteria
*/

#ifdef SM_STATS
#define STATS_ON(a) a
#else
#define STATS_ON(a)
#endif
#define STATS_OFF

/* if maximizing the min value */
#ifdef SM_MAXIMIZING
#define SM_MAXIMIZE(a) a
#else
#define SM_MAXIMIZE(a)
#endif

#include "SMuserDefs.h"

/* new initial point constants */
#define NONE            -1 
#define CENTROID         1

/* optimization step constants */
#define STEP_DONE        101
#define STEP_NOT_DONE    102

/* optimization termination constants */
#define STEP_ACCEPTED    100
#define IMP_TOO_SMALL    101
#define FLAT_NO_IMP      102
#define STEP_TOO_SMALL   103
#define EQUILIBRIUM       104
#define ZERO_SEARCH       105
#define MAX_ITER_EXCEEDED 106
#define LAP_ENOUGH        107

#define MIN_IMP          .001

#define SM_PSISROOT(procinfo) ((procinfo->myid == 0) ? 1 : 0)

#include "SMdata_structs.h"

#define      MEM_ERROR -25

/* memory allocation macros */
#ifndef MY_MALLOC
#define MY_MALLOC(a,b,c,d) { \
if (c == 0) { \
    printf("Size zero in MY_MALLOC!\n"); \
    a = NULL; \
} else { \
    a = b malloc(c); \
} \
}
#endif
#ifndef MY_MALLOCN
#define MY_MALLOCN(a,b,c,d) { \
if (c == 0) { \
    a = NULL; \
} else { \
    a = b malloc(c); \
} \
}
#endif
#ifndef MY_FREE
#define MY_FREE(a) \
{ \
    if (a != NULL) free(a); \
}
#endif
#ifndef MY_FREEN
#define MY_FREEN(a) \
{  \
    if (a != NULL) free(a); \
}
#endif
#ifndef SAFE_ACOS
#ifdef IRIX6
#include <values.h>
#endif
#define SAFE_ACOS(a) \
(((a) > 1.0) ? 0 : (((a) < -1.0) ? M_PI : acos(a)))
#endif

#ifndef SAFE_ASIN
#define SAFE_ASIN(a) \
(((a) > 1.0) ? M_PI/2 : (((a) < -1.0) ? -M_PI/2 : asin(a)))
#endif

#ifndef __cplusplus
/* Don't define these for C++, where you can use std::min and std::max
   instead. */ 
#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif
#ifndef MIN
#define MIN(a,b) (a < b ? a : b)
#endif
#endif

/* list macros */
#define INIT_LIST(a) \
{ \
  (a)->ptr = NULL; \
}

#define GET_NEXT_LIST_ITEM(list,data) \
{ \
  ADV_LIST(list); \
  data = (list)->ptr; \
}	

/* printing macros used in debugging */
#define PRINT_ORDERED_PTS(local_mesh) \
{ \
  int i99,j99; \
  for (i99=0;i99<local_mesh->dimension;i99++) \
      printf(" free_vtx[%d] = %f; ",i99,local_mesh->free_vtx[i99]); \
  printf("\n"); \
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) \
          printf(" vtx_list[%d][%d]= %f;",i99,j99,local_mesh->incident_vtx[i99][j99]); \
      printf("\n"); \
  } \
}

#define WRITE_ORDERED_PTS(fp,local_mesh) \
{ \
  int i99,j99; \
  fprintf(fp,"%d  %d\n",local_mesh->num_incident_vtx, local_mesh->num_tri);\
  for (i99=0;i99<local_mesh->dimension;i99++) \
      fprintf(fp,"%f  ",local_mesh->free_vtx[i99]); \
  fprintf(fp,"\n"); \
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) \
          fprintf(fp,"%f  ",local_mesh->incident_vtx[i99][j99]); \
      fprintf(fp,"\n"); \
  } \
  for (i99=0;i99<local_mesh->num_tri;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) \
          fprintf(fp,"%d  ",local_mesh->vtx_connectivity[i99][j99]); \
      fprintf(fp,"\n"); \
  } \
}
#define WRITE_BINARY_ORDERED_PTS(local_mesh) \
{ \
  int i99,j99; \
  int fd99; \
  char filename99[128]; \
  double temp99; \
  sprintf(filename99,"test.data"); \
  if ((fd99 = creat(filename99, 0666)) == -1) { \
     printf("cannot create filename for writing\n"); \
     exit(0); \
  } \
  temp99 = (double) local_mesh->num_incident_vtx; \
  write(fd99,&temp99,sizeof(double)); \
  temp99 = (double) local_mesh->num_tri;\
  write(fd99,&temp99,sizeof(double)); \
  for (i99=0;i99<local_mesh->dimension;i99++) {\
      temp99 = local_mesh->original_pt[i99]; \
      write(fd99,&temp99,sizeof(double)); \
  } \
  for (i99=0;i99<local_mesh->num_incident_vtx;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++){ \
          temp99 = local_mesh->incident_vtx[i99][j99]; \
          write(fd99,&temp99,sizeof(double)); \
      }\
  } \
  if (local_mesh->dimension ==3 ) {\
    for (i99=0;i99<local_mesh->num_tri;i99++) { \
      for (j99=0;j99<local_mesh->dimension;j99++) {\
          temp99 = (double) local_mesh->vtx_connectivity[i99][j99]; \
          write(fd99,&temp99,sizeof(double)); \
      } \
    } \
  } \
  close(fd99); \
}

#define LOCAL_MIN_VOLUME(local_mesh) \
{ \
  int i99,j99; \
  int num_tet99; \
  int num_values99; \
  double vtx99[4][3]; \
  double function99[6]; \
  double min_function99; \
  min_function99=1E300;\
  num_tet99 = local_mesh->num_tri;\
  vtx99[0][0] = local_mesh->free_vtx[0]; \
  vtx99[0][1] = local_mesh->free_vtx[1]; \
  vtx99[0][2] = local_mesh->free_vtx[2]; \
  for (i99=0;i99<num_tet99;i99++) { \
      for (j99=0;j99<3;j99++) {\
         vtx99[j99+1][0]=local_mesh->incident_vtx[local_mesh->vtx_connectivity[i99][j99]][0];\
         vtx99[j99+1][1]=local_mesh->incident_vtx[local_mesh->vtx_connectivity[i99][j99]][1];\
         vtx99[j99+1][2]=local_mesh->incident_vtx[local_mesh->vtx_connectivity[i99][j99]][2];\
      }\
      vComputeTetVolume(vtx99[0],vtx99[1],vtx99[2],vtx99[3],function99,&num_values99); \
      if (function99[0]<min_function99) min_function99=function99[0];\
  }\
  printf("Minimum volume in the local submesh is %f\n",min_function99);\
}

#define PRINT_FUNCTION_VALUES(opt_info) \
{ \
  int i99; \
  for (i99=0;i99<opt_info->num_values;i99++) { \
      printf("Index %d Function Value %f \n",i99,opt_info->function[i99]); \
  } \
}

#define SM_RECORD_ITER_VALUE(opt_info) \
{ \
    opt_info->prev_active_values[opt_info->iter_count] = opt_info->active->true_active_value; \
}

#define SM_DOT(c,a,b,n) {\
  if (n==2) c = a[0]*b[0] + a[1]*b[1]; \
  else if (n==3) c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];\
  else { \
    int i99; \
    c = 0; \
    for (i99=0;i99<n;i99++) c += a[i99]*b[i99]; \
  } \
}

#define SM_NORMALIZE(v,n) {\
    int i99; \
    double mag99; \
    if (n==2){ \
       mag99 = sqrt(v[0]*v[0] + v[1]*v[1]) ; \
       if (mag99 != 0) { \
          v[0] = v[0]/mag99; \
          v[1] = v[1]/mag99; \
       } \
    } else if (n==3) {\
     mag99 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ; \
     if (mag99 != 0) { \
         v[0] = v[0]/mag99; \
         v[1] = v[1]/mag99; \
         v[2] = v[2]/mag99; \
     } \
   } else { \
     mag99 = 0; \
     for (i99=0;i99<n;i99++) mag99+=v[i99]+v[i99]; \
     if (mag99 != 0) { \
       for (i99=0;i99<n;i99++) v[i99] = v[i99]/mag99;\
     } \
   }\
}

#define SM_COPY_VECTOR(a,b,n) { \
  int i99; \
  if (n==2) { \
     a[0] = b[0];  a[1] = b[1];  \
  } else if (n==3) {\
     a[0] = b[0];  a[1] = b[1];  a[2] = b[2]; \
  } else { \
     for (i99=0;i99<n;i99++) a[i99] = b[i99]; \
  } \
}

#define PRINT_MATRIX(nn,mm,AA) { \
  int i99, j99; \
  for (i99=0;i99<nn;i99++) { \
    for (j99=0;j99<mm;j99++) { \
      printf("A[%d][%d]=%f",i99,j99,AA[i99][j99]);\
    } \
    printf("\n");\
  }\
}

#endif
