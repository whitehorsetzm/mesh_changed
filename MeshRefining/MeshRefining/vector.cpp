/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2007年10月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 15, 2007
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953165
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：vector.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套密度控制机制
 *     Define a set of element spacing controlling scheme.
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include "vector.h"
#include <math.h>


/* *********************************************************************************
 * 定义三维矢量
 * define 3-dimensional vector
 * **********************************************************************************/

/* 矢量求模 calc. magnitude */
double Vector::magnitude() const
{
	return sqrt(x * x + y * y + z * z);
}

/* 归一 normalization */
void Vector::normalize()
{
	double m = magnitude();
	x = x / m;
	y = y / m;
	z = z / m;
}

double Vector::getDistance(const Vector &a)
{
    double temp=(a.x-x)*(a.x-x)+(a.y-y)*(a.y-y)+(a.z-z)*(a.z-z);
    double dis =sqrt(temp);
    return dis;
}

/* **********************************
 * 重载操作符 overloaded operator 
 * **********************************/

/* 赋值 assignment */
const Vector& Vector::operator = (const Vector& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}

const Vector& Vector::operator = (double v)
{
	x = y = z = v;
	return *this;
}

/* 四则运算 simple mathamatic operator */
const Vector Vector::operator + (const Vector& v) const
{
	Vector vr;
	vr.x = x + v.x;
	vr.y = y + v.y;
	vr.z = z + v.z;
	return vr;
}

const Vector Vector::operator - (const Vector& v) const
{
	Vector vr;
	vr.x = x - v.x;
	vr.y = y - v.y;
	vr.z = z - v.z;
	return vr;
}

const Vector Vector::operator * (double s) const
{
	Vector vr;
	vr.x = x * s;
	vr.y = y * s;
	vr.z = z * s;
	return vr;
}
 
const Vector Vector::operator / (double s) const
{
	Vector vr;
	vr.x = x / s;
	vr.y = y / s;
	vr.z = z / s;
	return vr;
}


const Vector& Vector::operator += (const Vector& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

const Vector& Vector::operator -= (const Vector& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

const Vector& Vector::operator += (double delta)
{
	x += delta;
	y += delta;
	z += delta;
	return *this;
}

const Vector& Vector::operator -= (double delta)
{
	x -= delta;
	y -= delta;
	z -= delta;
	return *this;
}

bool Vector::operator==(const Vector &v) const
{
	return x == v.x && y == v.y && z == v.z;
}

bool Vector::operator!=(const Vector &v) const
{
	return !(*this == v);
}


/* 点积 & 叉积 dot & cross */
double Vector::operator * (const Vector& v) const
{
	return x * v.x + y * v.y + z * v.z;
}

const Vector Vector::operator ^ (const Vector& v) const
{
	Vector vr;
    vr.x = y * v.z - v.y * z;
	vr.y = v.x * z - x * v.z;
	vr.z = x * v.y - y * v.x;
	return vr;
}


const Vector Vector::operator - ()
{
	Vector vr;
	vr.x = -x;
	vr.y = -y;
	vr.z = -z;
	return vr;
}

const Vector operator * (double scale, const Vector& v)
{
	 return v * scale;
}

