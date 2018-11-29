#pragma once

#ifndef _GEOMT_H
#define _GEOMT_H

#include <cassert>
#include <math.h>

//NOTE: This is really a bad usage, but there's no good class template alias for now.
#define PointT VectorT
#define DimensionT VectorT

template<int dim, typename T = double>
class VectorT
{
public:
	VectorT() { 
		for (int i = 0; i < dim; ++i) {
			v[i] = 0;
		}
	}
	VectorT(double val) {
		for (int i = 0; i < dim; ++i)
		{
			v[i] = val;
		}
	}
	VectorT(double v1, double v2) {
		v[0] = v1;
		v[1] = v2;
		if (dim == 3) {
			v[2] = 0;
		}
	}
	VectorT(double v1, double v2, double v3) {
		assert(dim == 3);
		v[0] = v1;
		v[1] = v2;
		v[2] = v3;
	}
	VectorT(const VectorT& vec) {
		for (int i = 0; i < dim; ++i) {
			v[i] = vec.v[i];
		}
	}
	
	/* 矢量求模 calc. magnitude */
	double magnitude() const;
	
	/* 归一 normalization */
	void normalize();

	/* **********************************
	 * 重载操作符 overloaded operator 
	 * **********************************/

	/* 赋值 assignment */
	inline const VectorT& operator = (const VectorT& v);
	inline const VectorT& operator = (double v);

	/* 四则运算 simple mathematic operator */
    inline const VectorT operator + (const VectorT& v) const;
    inline const VectorT operator - (const VectorT& v) const;
    inline const VectorT operator * (double s) const;  
    inline const VectorT operator / (double s) const;
    
    inline const VectorT& operator += (const VectorT& v);
    inline const VectorT& operator -= (const VectorT& v);
    inline const VectorT& operator += (double delta);
    inline const VectorT& operator -= (double delta);

	/* 取值 */
	inline T& operator[](int i);

	/* 点积 dot */
	inline double operator * (const VectorT& v) const;
	/* 叉积 cross */
	inline const VectorT operator ^ (const VectorT& V) const;

	friend const VectorT operator - (const VectorT& p);
	friend const VectorT operator * (double scale, const VectorT& v);

	/* 视为2进制正整数进行转换 */
	static int vec2int(const VectorT&);
	static const VectorT int2vec(int i);

public:
	T v[dim];
};

/* 矢量求模 calc. magnitude */
template<int dim, typename T>
double VectorT<dim, T>::magnitude() const
{
	double s = 0; 
	for (int i = 0; i < dim; ++i)
	{
		s += v[i] * v[i];
	}
	return sqrt(s);
}

/* 归一 normalization */
template<int dim, typename T>
void VectorT<dim, T>::normalize()
{
	double m = magnitude();
	for (int i = 0; i < dim; ++i)
	{
		v[i] /= m;
	}
}

/* **********************************
* 重载操作符 overloaded operator 
* **********************************/

/* 赋值 assignment */
template<int dim, typename T>
const VectorT<dim, T>& VectorT<dim, T>::operator = (const VectorT<dim, T>& vec)
{
	for (int i = 0; i < dim; ++i) {
		v[i] = vec.v[i];
	}
	return *this;
}

template<int dim, typename T>
const VectorT<dim, T>& VectorT<dim, T>::operator = (double val)
{
	for (int i = 0; i < dim; ++i) {
		v[i] = val;
	}
	return *this;
}

/* 四则运算 simple mathematic operator */
template<int dim, typename T>
const VectorT<dim, T> VectorT<dim, T>::operator + (const VectorT<dim, T>& vec) const
{
	VectorT<dim> vr;
	for (int i = 0; i < dim; ++i) {
		vr.v[i] = v[i] + vec.v[i];
	}
	return vr;
}

template<int dim, typename T>
const VectorT<dim, T> VectorT<dim, T>::operator - (const VectorT<dim, T>& vec) const
{
	VectorT<dim> vr;
	for (int i = 0; i < dim; ++i) {
		vr.v[i] = v[i] - vec.v[i];
	}
	return vr;
}

template<int dim, typename T>
const VectorT<dim, T> VectorT<dim, T>::operator * (double s) const
{
	VectorT<dim> vec;
	for (int i = 0; i < dim; ++i) {
		vec[i] = v[i] * s;
	}
	return vec;
}

template<int dim, typename T>
const VectorT<dim, T> VectorT<dim, T>::operator / (double s) const
{
	VectorT<dim> vec;
	for (int i = 0; i < dim; ++i) {
		vec[i] = v[i] / s;
	}
	return vec;
}

template<int dim, typename T>
const VectorT<dim, T>& VectorT<dim, T>::operator += (const VectorT<dim, T>& vec)
{
	for (int i = 0; i < dim; ++i)
	{
		v[i] += vec.v[i];
	}
	return *this;
}

template<int dim, typename T>
const VectorT<dim, T>& VectorT<dim, T>::operator -= (const VectorT<dim, T>& vec)
{
	for (int i = 0; i < dim; ++i)
	{
		v[i] -= vec.v[i];
	}
	return *this;
}

template<int dim, typename T>
const VectorT<dim, T>& VectorT<dim, T>::operator += (double delta)
{
	for (int i = 0; i < dim; ++i)
	{
		v[i] += delta;
	}
	return *this;
}

template<int dim, typename T>
const VectorT<dim, T>& VectorT<dim, T>::operator -= (double delta)
{
	for (int i = 0; i < dim; ++i)
	{
		v[i] -= delta;
	}
	return *this;
}

/* 点积 dot */
template<int dim, typename T>
double VectorT<dim, T>::operator * (const VectorT<dim, T>& vec) const
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
	{
		s += v[i] * vec.v[i];
	}
	return s;
}

/* 叉积 cross */
template<int dim, typename T>
const VectorT<dim, T> VectorT<dim, T>::operator ^ (const VectorT<dim, T>& vec) const
{
	assert(dim == 3);
	VectorT<dim> vr;
	vr.v[0] = v[1] * vec.v[2] - vec.v[1] * v[2];
	vr.v[1] = vec.v[0] * v[2] - v[0] * vec.v[2];
	vr.v[2] = v[0] * vec.v[1] - v[1] * vec.v[0];
	return vr;
}

template<int dim, typename T>
T& VectorT<dim, T>::operator [] (int i) 
{
	return v[i];
}


template<int dim, typename T>
const VectorT<dim, T> operator - (const VectorT<dim, T>& vec)
{
	VectorT<dim> vr;
	for (int i = 0; i < dim; ++i)
	{
		vr.v[i] = - vec.v[i];
	}
	return vr;
}

template<int dim, typename T>
const VectorT<dim, T> operator * (double scale, const VectorT<dim, T>& vec)
{
	VectorT<dim> vr;
	for (int i = 0; i < dim; ++i)
	{
		vr.v[i] = vec.v[i] * scale;
	}
	return vr;
}

template<int dim, typename T>
int VectorT<dim, T>::vec2int(const VectorT<dim, T> & vec)
{
	int n = 0;
	for (int i = 0; i < dim; ++i)
	{
		if (vec.v[i]) {
			n += 1 << (dim - 1 - i);
		}
	}
	return n;
}

template<int dim, typename T>
const VectorT<dim, T> VectorT<dim, T>::int2vec(int n)
{
	VectorT<dim, T> vec;
	for (int j = dim - 1; j >= 0; --j)
	{
		int v = n & 1;
		vec[j] = v;
		n >>= 1;
	}
	return vec;
}

typedef VectorT<2> Vec2d;
typedef VectorT<3> Vec3d;
typedef VectorT<2> Point2d;
typedef VectorT<3> Point3d;
typedef VectorT<2> Dimension2d;
typedef VectorT<3> Dimension3d;

template<int dim>
struct BoxT {
	PointT<dim> location;
	DimensionT<dim> size;
	BoxT() : location(0), size(0) {}
	BoxT(PointT<dim> loc, DimensionT<dim> sz) : location(loc), size(sz) {}
	BoxT(double x1, double y1, double x2, double y2) {
		location[0] = x1;
		location[1] = y1;
		size[0] = x2 - x1;
		size[1] = y2 - y1;
	}
	BoxT(double x1, double y1, double z1, double x2, double y2, double z2) {
		location[0] = x1;
		location[1] = y1;
		location[2] = z1;
		size[0] = x2 - x1;
		size[1] = y2 - y1;
		size[2] = z2 - z1;
	}
	PointT<dim> center() const {
		return location + size / 2;
	}
	PointT<dim> vertex(int index) {
		assert(index < (2 << dim));
		PointT<dim> coord = location;
		VectorT<dim, int> vec = VectorT<dim, int>::int2vec(index);
		for (int i = 0; i < dim; ++i)
		{
			coord[i] += vec[i] * size[i];
		}
		return coord;
	}
};

typedef BoxT<3> Box3d;
typedef BoxT<2> Box2d;

#endif
