#ifndef GR_Vec
#define GR_Vec 1

#include "GR_config.h"
#include <math.h>

#define dMAG2D(a) sqrt(a[0]*a[0] + a[1]*a[1])
#define dMAG3D(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])

#define dDIST2D(a, b) sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
#define dDIST_SQ_2D(a, b) ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))

#define dDIST3D(a, b) sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) +  (a[2]-b[2])*(a[2]-b[2]))
#define dDIST_SQ_3D(a, b) ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) +  (a[2]-b[2])*(a[2]-b[2]))

#define adDIFF3D(a, b) {a[0]-b[0], a[1]-b[1], a[2]-b[2]}
#define adDIFF2D(a, b) {a[0]-b[0], a[1]-b[1]}
#define vSWAP3D(a, b) do {register double adTemp__; \
adTemp__ = b[0]; b[0] = a[0]; a[0] = adTemp__; \
adTemp__ = b[1]; b[1] = a[1]; a[1] = adTemp__; \
adTemp__ = b[2]; b[2] = a[2]; a[2] = adTemp__; \
} while(0)

#define vSCALE2D(a,d) {a[0]*=d; a[1]*=d;}
#define vSCALE3D(a,d) {a[0]*=d; a[1]*=d; a[2]*=d;}

#define vNORMALIZE2D(a) do {double dInvMag = 1./dMAG2D(a); \
a[0]*=dInvMag; a[1]*=dInvMag;} while(0)
#define vNORMALIZE3D(a) do {double dInvMag = 1./dMAG3D(a); \
a[0]*=dInvMag; a[1]*=dInvMag; a[2]*=dInvMag;} while(0)

#define dDOT2D(a,b) (a[0]*b[0] + a[1]*b[1])
#define dDOT3D(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define dCROSS2D(a,b) (a[0]*b[1] - a[1]*b[0])
#define vCROSS3D(a,b,res) do { \
res[0] = a[1]*b[2] - a[2]*b[1]; \
res[1] = a[2]*b[0] - a[0]*b[2]; \
res[2] = a[0]*b[1] - a[1]*b[0]; \
} while(0)

#endif
