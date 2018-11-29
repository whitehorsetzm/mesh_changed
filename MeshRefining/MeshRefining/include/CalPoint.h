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
}//pa��pb�ĵ�� (2ά)
double Cross(CalPoint p,CalPoint a,CalPoint b){
	return (a.y-p.y)*(b.x-p.x)-(a.x-p.x)*(b.y-p.y);
}//pa��pb�Ĳ��(2ά)

double Cross(CalPoint a,CalPoint b){
	return a.y*b.x-a.x*b.y;
}//����a��b�Ĳ��(2ά)

double Dot(CalPoint a,CalPoint b){
	return a.x*b.x+a.y*b.y;
}//����a��b�ĵ��(2ά)



double ThreeDot(CalPoint p,CalPoint a,CalPoint b){
	return(a.x-p.x)*(b.x-p.x)+(a.y-p.y)*(b.y-p.y)+(a.z-p.z)*(b.z-p.z);
}//pa��pb�ĵ�� (3ά)    


CalPoint ThreeCross(CalPoint p,CalPoint a,CalPoint b){
	CalPoint C;
	C.x=(a.y-p.y)*(b.z-p.z)-(a.z-p.z)*(b.y-p.y);
	C.y=(a.z-p.z)*(b.x-p.x)-(a.x-p.x)*(b.z-p.z);
	C.z=(a.x-p.x)*(b.y-p.y)-(a.y-p.y)*(b.x-p.x);
	return C;
}//pa��pb�Ĳ��(3ά)

double ThreeDot(CalPoint a,CalPoint b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}//����a��b�ĵ�� (3ά)

CalPoint ThreeCross(CalPoint a,CalPoint b){           
	CalPoint C;
	C.x=a.y*b.z-a.z*b.y;
	C.y=a.z*b.x-a.x*b.z;
	C.z=a.x*b.y-a.y*b.x;
	return C;
}//����a��b�Ĳ��(3ά)

double TetrahedronArea(CalPoint a,CalPoint b,CalPoint c,CalPoint d){
	CalPoint temp(d.x-a.x,d.y-a.y,d.z-a.z);
	double ans=ThreeDot(ThreeCross(a,b,c),temp)/6;
	return ans;
}