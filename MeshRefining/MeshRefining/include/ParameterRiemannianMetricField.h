// ParameterRiemannianMetricField.h: interface for the ParameterRiemannianMetricField class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARAMETERRIEMANNIANMETRICFIELD_H__2DB7E3E6_4F1A_4B76_A7D7_29B83782B585__INCLUDED_)
#define AFX_PARAMETERRIEMANNIANMETRICFIELD_H__2DB7E3E6_4F1A_4B76_A7D7_29B83782B585__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "MetricField.h"
#include "ferguson_surface.h"
#include "spacing.h"

class ParameterRiemannianMetricField : public MetricField
{
private:
	EEMAS::FergusonSurface * surface;
	SurfBndBKGMesh *geometrySizeField;
	GlobalSpacing * globalSpacing;
public:
	ParameterRiemannianMetricField(EEMAS::FergusonSurface *fs, SurfBndBKGMesh * field, GlobalSpacing *gs) 
		: surface(fs), geometrySizeField(field), globalSpacing(gs) {}
	virtual ~ParameterRiemannianMetricField() {}
	virtual void getMetric(double u, double v, Mat2by2d mt);

};

#endif // !defined(AFX_PARAMETERRIEMANNIANMETRICFIELD_H__2DB7E3E6_4F1A_4B76_A7D7_29B83782B585__INCLUDED_)
