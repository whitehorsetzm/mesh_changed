#ifndef GR_AdaptPred
#define GR_AdaptPred 1

#include "GR_config.h"

#ifdef __cplusplus
extern "C" {
#endif
extern bool qAdaptPred;

void exactinit(void);
double orient2d_shew(const double * const pa, const double * const pb,
		     const double * const pc);
double orient3d_shew(const double * const pa, const double * const pb,
		     const double * const pc, const double * const pd);
double incircle_shew(const double * const pa, const double * const pb,
		     const double * const pc, const double * const pd);
double insphere_shew(const double * const pa, const double * const pb,
		     const double * const pc, const double * const pd,
		     const double * const pe);

#ifdef __cplusplus
}
#endif
  
#endif
