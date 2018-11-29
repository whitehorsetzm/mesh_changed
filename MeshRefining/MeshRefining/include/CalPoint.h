#include<iostream>
#include"math.h"
using namespace std;
///modify by zhvliu
class CalPoint
{

public:
	CalPoint()
	{
		x=0;
		y=0;
		z=0;
	}
	CalPoint(double xx,double yy,double zz)
	{
		x=xx;
		y=yy;
		z=zz;

	}


public:
	double x;
	double y;
	double z;



};

double Dot(CalPoint p,CalPoint a,CalPoint b)
{
	return(a.x-p.x)*(b.x-p.x)+(a.y-p.y)*(b.y-p.y);
}//pa与pb的点积 (2维)
double Cross(CalPoint p,CalPoint a,CalPoint b){
	return (a.y-p.y)*(b.x-p.x)-(a.x-p.x)*(b.y-p.y);
}//pa与pb的叉积(2维)

double Cross(CalPoint a,CalPoint b){
	return a.y*b.x-a.x*b.y;
}//向量a与b的叉积(2维)

double Dot(CalPoint a,CalPoint b){
	return a.x*b.x+a.y*b.y;
}//向量a和b的点积(2维)



double ThreeDot(CalPoint p,CalPoint a,CalPoint b){
	return(a.x-p.x)*(b.x-p.x)+(a.y-p.y)*(b.y-p.y)+(a.z-p.z)*(b.z-p.z);
}//pa与pb的点积 (3维)    


CalPoint ThreeCross(CalPoint p,CalPoint a,CalPoint b){
	CalPoint C;
	C.x=(a.y-p.y)*(b.z-p.z)-(a.z-p.z)*(b.y-p.y);
	C.y=(a.z-p.z)*(b.x-p.x)-(a.x-p.x)*(b.z-p.z);
	C.z=(a.x-p.x)*(b.y-p.y)-(a.y-p.y)*(b.x-p.x);
	return C;
}//pa与pb的叉积(3维)

double ThreeDot(CalPoint a,CalPoint b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}//向量a与b的点积 (3维)

CalPoint ThreeCross(CalPoint a,CalPoint b){           
	CalPoint C;
	C.x=a.y*b.z-a.z*b.y;
	C.y=a.z*b.x-a.x*b.z;
	C.z=a.x*b.y-a.y*b.x;
	return C;
}//向量a与b的叉积(3维)

double TetrahedronArea(CalPoint a,CalPoint b,CalPoint c,CalPoint d){
	CalPoint temp(d.x-a.x,d.y-a.y,d.z-a.z);
	double ans=ThreeDot(ThreeCross(a,b,c),temp)/6;
	return ans;
}