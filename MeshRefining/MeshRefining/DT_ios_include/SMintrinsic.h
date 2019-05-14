#ifndef SM_INTRINSIC_H
#define SM_INTRINSIC_H 1

#ifndef ad_intrinsic_h
#define ad_intrinsic_h 1
#include "SM_config.h"

#ifdef IN_AD_INTRINSIC_C
#define EXTERN 
#else
#define EXTERN extern
#endif


#ifdef __STDC__
#define USE_PROTOTYPES
#define USE_FN_ARGS
#endif

/* Maybe some C++ does not define __STDC__ ? */

#ifdef __cplusplus 
#define USE_PROTOTYPES
#define USE_FN_ARGS
#endif

#ifdef USE_PROTOTYPES
#define Proto(x) x
#else
#define Proto(x) ()
#endif

#define ONE_ARG_FN(x) EXTERN void x Proto((DERIV_TYPE *,DERIV_TYPE))
#define TWO_ARG_FN(x) EXTERN void x Proto((DERIV_TYPE *,DERIV_TYPE, DERIV_TYPE))

#ifdef __cplusplus
extern "C" {
#endif

ONE_ARG_FN (g_ad_log);
ONE_ARG_FN (g_ad_log10);
ONE_ARG_FN (g_ad_sqrt);
ONE_ARG_FN (g_ad_cos);
ONE_ARG_FN (g_ad_sin);
ONE_ARG_FN (g_ad_tan);
ONE_ARG_FN (g_ad_cosh);
ONE_ARG_FN (g_ad_sinh);
ONE_ARG_FN (g_ad_tanh);
ONE_ARG_FN (g_ad_exp);
TWO_ARG_FN (g_ad_pow);
TWO_ARG_FN (g_ad_fmin);
TWO_ARG_FN (g_ad_fmax);
TWO_ARG_FN (g_ad_fmod);
ONE_ARG_FN (g_ad_fabs);
ONE_ARG_FN (g_ad_acos);
ONE_ARG_FN (g_ad_asin);
ONE_ARG_FN (g_ad_atan);
TWO_ARG_FN (g_ad_atan2);

#ifdef __cplusplus
}  /* end extern "C" */
#endif

#undef EXTERN
#endif
#endif
