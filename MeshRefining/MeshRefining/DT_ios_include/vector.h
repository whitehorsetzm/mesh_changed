/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#ifndef __iso3d_vector_h__
#define	__iso3d_vector_h__

/* *********************************************************************************
 * 定义三维矢量
 * define 3-dimensional vector
 * **********************************************************************************/
class Vector
{
public:
	Vector() { x = y = z = 0.; }
	Vector(double v1, double v2, double v3) {
		x = v1; y = v2; z = v3;
	}
	Vector(const Vector& v) {
		x = v.x;  y = v.y; z = v.z;
	}
	
	void set(double x1, double y1, double z1)
	{
		x = x1; y = y1; z = z1;
	}

	/* 矢量求模 calc. magnitude */
	double magnitude();
	
	/* 归一 normalization */
	void normalize();

	double length_squared()
	{
		return x*x + y*y + z*z;
	}

	/* **********************************
	 * 重载操作符 overloaded operator 
	 * **********************************/

	/* 赋值 assignment */
	const Vector& operator = (const Vector& v);
	const Vector& operator = (double v);

	/* 四则运算 simple mathamatic operator */
    const Vector operator + (const Vector& v) const;
    const Vector operator - (const Vector& v) const;
    const Vector operator * (double s) const;  
    const Vector operator / (double s) const;
    
    const Vector& operator += (const Vector& v);
    const Vector& operator -= (const Vector& v);
    const Vector& operator += (double delta);
    const Vector& operator -= (double delta);

	/* 点积 & 叉积 dot & cross */
	double operator * (const Vector& v) const;
	const Vector operator ^ (const Vector& V) const;

//	friend const Vector operator - (const Vector& p);
	friend const Vector operator * (double scale, const Vector& v);

	//- vector cross product, non-commutative
	//friend Vector operator*(const Vector &v1, const Vector &v2);
	
public:
	double x, y, z;
};

#endif /* __iso3d_vector_h__*/