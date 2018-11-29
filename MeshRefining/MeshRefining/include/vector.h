/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ά����ͬ��Delaunay���������� (�汾�ţ�0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2005��9��15��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * ��ϵ��ʽ
 *   �绰��+86-571-87953165
 *   ���棺+86-571-87953167
 *   ���䣺zdchenjj@yahoo.com.cn
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
 * ������άʸ��
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
	
	/* ʸ����ģ calc. magnitude */
	double magnitude() const;
	
	/* ��һ normalization */
	void normalize();

    //calculate the distance
    double getDistance(const Vector &a);

	/* **********************************
	 * ���ز����� overloaded operator 
	 * **********************************/

	/* ��ֵ assignment */
	const Vector& operator = (const Vector& v);
	const Vector& operator = (double v);

	/* �������� simple mathamatic operator */
    const Vector operator + (const Vector& v) const;
    const Vector operator - (const Vector& v) const;
    const Vector operator * (double s) const;  
    const Vector operator / (double s) const;
    
    const Vector& operator += (const Vector& v);
    const Vector& operator -= (const Vector& v);
    const Vector& operator += (double delta);
    const Vector& operator -= (double delta);

	/* ��� & ��� dot & cross */
	double operator * (const Vector& v) const;
	const Vector operator ^ (const Vector& V) const;
	const Vector operator - ();

//	friend const Vector operator - (const Vector& p);
	friend const Vector operator * (double scale, const Vector& v);

public:
	double x, y, z;
};

#endif /* __iso3d_vector_h__*/
