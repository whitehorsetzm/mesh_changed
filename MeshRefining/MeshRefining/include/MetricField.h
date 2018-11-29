// MetricField.h: interface for the MetricField class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_METRICFIELD_H__8DEE7F9B_3817_45DC_8C46_9763DE42BA69__INCLUDED_)
#define AFX_METRICFIELD_H__8DEE7F9B_3817_45DC_8C46_9763DE42BA69__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

typedef double Mat1by1d[1];
typedef double Mat1by2d[2];
typedef double Mat2by1d[2];
typedef double Mat2by2d[4];

class MetricField
{
public:
	
	MetricField(void) {}
	
	virtual ~MetricField(void) {}
	
	virtual void getMetric(double x, double y, Mat2by2d mt) = 0;
	
	virtual double getDirectionalSpacing(double x, double y, double dx, double dy);

	virtual double calcDistance(double x1, double y1, double x2, double y2);

	virtual int isDelaunayBroken(double vx[3], double vy[3], double x, double y);

	virtual void calcCircumCenter(Mat2by2d mt, double vx[3], double vy[3], double *cx, double *cy);

private:
	void multiply(Mat2by1d du, Mat2by2d mt, Mat1by1d val);
};

#endif // !defined(AFX_METRICFIELD_H__8DEE7F9B_3817_45DC_8C46_9763DE42BA69__INCLUDED_)
