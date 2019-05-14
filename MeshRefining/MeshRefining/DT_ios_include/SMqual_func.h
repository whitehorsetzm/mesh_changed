#ifndef SM_QUAL_FUNC_H
#define SM_QUAL_FUNC_H 1

typedef void (*SMfunction_ptr2D)(double *, double *, double *, double *, int *);
typedef void (*SMgradfunc_ptr2D)(double *, double *, double *, double(*)[2], int *);

typedef void (*SMfunction_ptr3D)(const double adCoord0[3], 
                        const double adCoord1[3], const double adCoord2[3], 
                        const double adCoord3[3], double *adResult, 
                        int* const piNGradient);

typedef void (*SMgradfunc_ptr3D)(const double adCoord0[3], 
                        const double adCoord1[3],const double adCoord2[3], 
                        const double adCoord3[3],double (*adGradient)[3], 
                        int* const piNGradient);
#ifdef __cplusplus
extern "C" {
#endif
void SMset2DUserQualityFunction(void *ext_smooth_data, 
				int values_per_cell,
				SMfunction_ptr2D userQualityFunc,
				SMgradfunc_ptr2D userQualityGrad);
void SMset3DUserQualityFunction(void *ext_smooth_data, 
				int values_per_cell,
				SMfunction_ptr3D userQualityFunc,
				SMgradfunc_ptr3D userQualityGrad);
#ifdef __cplusplus
}
#endif


#endif

