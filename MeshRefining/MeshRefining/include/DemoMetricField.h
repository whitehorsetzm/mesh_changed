// DemoMetricField.h: interface for the DemoMetricField class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DEMOMETRICFIELD_H__63C7A639_7E98_4031_9AD6_8A08C7D7E478__INCLUDED_)
#define AFX_DEMOMETRICFIELD_H__63C7A639_7E98_4031_9AD6_8A08C7D7E478__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "MetricField.h"

class DemoMetricField : public MetricField  
{
public:
	DemoMetricField() {}
	virtual ~DemoMetricField() {}

	virtual void getMetric(double x, double y, Mat2by2d mt);
};

#endif // !defined(AFX_DEMOMETRICFIELD_H__63C7A639_7E98_4031_9AD6_8A08C7D7E478__INCLUDED_)
