#ifndef SM_DIHED_FUNC_H
#define SM_DIHED_FUNC_H 1

#include "SM_config.h"

double dComputeTetJacobian(const double adCoord0[3], const double adCoord1[3],
                 	             const double adCoord2[3], const double adCoord3[3]);

void vSolids(const double adCoord0[3], const double adCoord1[3],
	     const double adCoord2[3], const double adCoord3[3],
	     double adResult[4], int* const piNResult);
void vDihedrals(const double adCoord0[3], const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult);
void vNegateDihedrals(const double adCoord0[3], 
                const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult);
void vSineDihedrals(const double adCoord0[3], 
                    const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
		    double adResult[6], int* const piNResult);
void vCosineDihedrals(const double adCoord0[3], 
                      const double adCoord1[3],
		      const double adCoord2[3], const double adCoord3[3],
		      double adResult[6], int* const piNResult);
void vNegateCosineDihedrals(const double adCoord0[3], 
                      const double adCoord1[3],
		      const double adCoord2[3], const double adCoord3[3],
		      double adResult[6], int* const piNResult);
void vGradSineDihedrals(const double adCoord0[3], 
                        const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double (*adGradient)[3], int* const piNGradient);
void vGradCosineDihedrals(const double adCoord0[3], 
                        const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double (*adGradient)[3], int* const piNGradient);
void vNegateGradCosineDihedrals(const double adCoord0[3], 
                        const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double (*adGradient)[3], int* const piNGradient);
void vGradSolids(const double adCoord0[3], const double adCoord1[3],
		 const double adCoord2[3], const double adCoord3[3],
		 double (*adGradient)[3], int* const piNGradient);
void vGradDihedrals(const double adCoord0[3],  const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
		    double (*adGradient)[3], int* const piNGradient);
void vNegateGradDihedrals(const double adCoord0[3], 
                        const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double (*adGradient)[3], int* const piNGradient);
void vScaledJacobian(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult);
void vGradScaledJacobian(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double (*adGradient)[3], int* const piNGradient);

void vSMRSVolumeRatio(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult);
void vGradSMRSVolumeRatio(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double (*adGradient)[3], int* const piNGradient);
void vComputeTetVolume(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult);
void vNormJacSquared(const double adCoord0[3], 
                const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		     double adResult[6], int* const piNResult);
void vGradNormJacSquared(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
			 double (*adGradient)[3], int* const piNGradient);
void vCondition(const double adCoord0[3], const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult);
void vGradCondition(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
    		    double (*adGradient)[3], int* const piNGradient);


void  g_ad_vCross(DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3], DERIV_TYPE adResult[3]);
void  g_ad_dMagnitude(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVec[3]);
void  g_ad_dDot(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3]);
void  g_ad_dNegDot(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3]);
void  g_ad_vUnitNormal(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                       DERIV_TYPE adCoord2[3], DERIV_TYPE adResult[3]);
void g_ad_vSineDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                         DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                         DERIV_TYPE adResult[6], int  *piNResult);
void  g_ad_vCosineDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                                   DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                                   DERIV_TYPE adResult[6], int  *piNResult);
void  g_ad_vDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                             DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                             DERIV_TYPE adResult[6], int  *piNResult);

#endif

