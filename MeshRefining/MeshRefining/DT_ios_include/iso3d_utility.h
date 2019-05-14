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

#ifndef __iso3d_utility_h__
#define __iso3d_utility_h__

#include "iso3d_define.h"
#include "vector.h"
/* **********************************************************
 * ��Ԫʵ��(�ڵ�/��/��)�ı���/����
 * coding & decoding of elemental entities (nodes/edges/faces)
 * *********************************************************/

/* 
 * ������� DNC(Decoding Node Codes) 
 *    ��֪��A�����nc
 *    δ֪����Ԫ������3���ڵ����ib ic id
 *    Լ��:	abcd��һ���������ַ���������� 
 * DNC:	Decoding Node Codes
 *   Known:	nc, code of node A
 * Unknown:	ib, ic, id, codes of three other nodes
 *    s.t.:	abcd is a positive tet. 
 */
#define DNC(nc, ia, ib, ic, id) \
	ia = (nc);					\
	switch ((nc))				\
	{							\
	case 0:						\
		ib = 1; ic = 2; id = 3;	\
		break;					\
	case 1:						\
		ib = 3; ic = 2; id = 0;	\
		break;					\
	case 2:						\
		ib = 0; ic = 1; id = 3;	\
		break;					\
	case 3:						\
		ib = 2; ic = 1; id = 0;	\
		break;					\
	}	

/* 
 * ˫������� DDNC(Decoding Double Node Codes) 
 *    ��֪��C D�����ic id
 *    δ֪����Ԫ������2���ڵ����ia ib 
 *    Լ��:	abcd��һ���������ַ���������� 
 * DDNC: Decoding Double Node Codes
 *   Known:	ic/id, code of node C/D
 * Unknown:	ia/ib, codes of two other nodes
 *    s.t.:	abcd is a positive tet. 
 */
#define DDNC(ia, ib, ic, id)					\
	switch ((ic))								\
	{											\
	case 0:										\
		switch((id))							\
		{										\
		case 1:									\
			ia = 2; ib = 3;	break;				\
		case 2:									\
			ia = 3; ib = 1;	break;				\
		case 3:									\
			ia = 1; ib = 2;	break;				\
		}										\
		break;									\
	case 1:										\
		switch((id))							\
		{										\
		case 0:									\
			ia = 3; ib = 2;	break;				\
		case 2:									\
			ia = 0; ib = 3;	break;				\
		case 3:									\
			ia = 2; ib = 0;	break;				\
		}										\
		break;									\
	case 2:										\
		switch((id))							\
		{										\
		case 0:									\
			ia = 1; ib = 3;	break;				\
		case 1:									\
			ia = 3; ib = 0;	break;				\
		case 3:									\
			ia = 0; ib = 1;	break;				\
		}										\
		break;									\
	case 3:										\
		switch((id))							\
		{										\
		case 0:									\
			ia = 2; ib = 1;	break;				\
		case 1:									\
			ia = 0; ib = 2;	break;				\
		case 2:									\
			ia = 1; ib = 0;	break;				\
		}										\
		break;									\
	}		

/* 
 * ������� DEC(Decoding Edge Codes) 
 *    ��֪����BD��ֵib id
 *    δ֪����Ԫ������2���ڵ����ia ic 
 *    Լ��:	abcd��һ���������ַ���������� 
 * DEC: Decoding Edge Codes
 *   Known:	ib/id, code of an edge BD
 * Unknown:	ia/ic, codes of two other nodes
 *    s.t.:	abcd is a positive tet. 
 */
#define DEC(ec, ia, ib, ic, id)			\
	switch ((ec))						\
	{									\
	case 0:								\
		ia = 3; ib = 0; ic = 2; id = 1;	\
		break;							\
	case 1:								\
		ia = 1; ib = 0; ic = 3; id = 2;	\
		break;							\
	case 2:								\
		ia = 2; ib = 0; ic = 1; id = 3;	\
		break;							\
	case 3:								\
		ia = 3; ib = 1; ic = 0; id = 2;	\
		break;							\
	case 4:								\
		ia = 0; ib = 1; ic = 2; id = 3;	\
		break;							\
	case 5:								\
		ia = 1; ib = 2; ic = 0; id = 3;	\
		break;							\
	}

/* 
 * ˫������� DDEC(Decoding Double Edge Codes) 
 *    ��֪������Сֵ&�����ֵ miec & maec
 *    δ֪����Ԫ��4���ڵ���ֵia ib ic id
            ��������������ֵ ec1 ec2
			�����ߵĹ�ϵ ����/���
 *    Լ��:	abcd��һ���������ַ����������;
            ������ڱ� ec1 = AD ec2 = BD
			����ǶԱ� ec1 = AD ec2 = BC
 * DDEC: Decoding Double Edge Codes
 *   Known:	miec & maec, min. & max edge codes 
 * Unknown:	ia, ib, ic, id, codes of four nodes
            ec1, ec2, codes of two intersection edges 
			sc, relation of two edges, neighboring? opposite?
 *    s.t.:	abcd is a positive tet. 
            if two edges neighboring, ec1 = AD ec2 = BD
			if two edges opposite, ec1 = AD ec2 = BC
 */
#define DDEC(miec, maec, ia, ib, ic, id, ec1, ec2, sc)								\
	switch ((miec))																	\
	{																				\
	case 0:																			\
		switch ((maec))																\
		{																			\
		case 1:																		\
			ia = 2; ib = 1; ic = 3, id = 0; ec1 = 1; ec2 = 0; sc = SC_NEIG; break;	\
		case 2:																		\
			ia = 1; ib = 3; ic = 2, id = 0; ec1 = 0; ec2 = 2; sc = SC_NEIG; break;	\
		case 3:																		\
			ia = 0; ib = 2; ic = 3, id = 1; ec1 = 0; ec2 = 3; sc = SC_NEIG; break;	\
		case 4:																		\
			ia = 3; ib = 0; ic = 2, id = 1; ec1 = 4; ec2 = 0; sc = SC_NEIG; break;	\
		case 5:																		\
			ia = 0; ib = 2; ic = 3, id = 1; ec1 = 0; ec2 = 5; sc = SC_OPPO; break;	\
		}																			\
		break;																		\
	case 1:																			\
		switch ((maec))																\
		{																			\
		case 2:																		\
			ia = 3; ib = 2; ic = 1, id = 0; ec1 = 2; ec2 = 1; sc = SC_NEIG; break;	\
		case 3:																		\
			ia = 1; ib = 0; ic = 3, id = 2; ec1 = 3; ec2 = 1; sc = SC_NEIG; break;	\
		case 4:																		\
			ia = 0; ib = 3; ic = 1, id = 2; ec1 = 1; ec2 = 4; sc = SC_OPPO; break;	\
		case 5:																		\
			ia = 0; ib = 3; ic = 1, id = 2; ec1 = 1; ec2 = 5; sc = SC_NEIG; break;	\
		}																			\
		break;																		\
	case 2:																			\
		switch ((maec))																\
		{																			\
		case 3:																		\
			ia = 0; ib = 1; ic = 2, id = 3; ec1 = 2; ec2 = 3; sc = SC_OPPO; break;	\
		case 4:																		\
			ia = 0; ib = 1; ic = 2, id = 3; ec1 = 2; ec2 = 4; sc = SC_NEIG; break;	\
		case 5:																		\
			ia = 2; ib = 0; ic = 1, id = 3; ec1 = 5; ec2 = 2; sc = SC_NEIG; break;	\
		}																			\
		break;																		\
	case 3:																			\
		switch ((maec))																\
		{																			\
		case 4:																		\
			ia = 2; ib = 3; ic = 0, id = 1; ec1 = 3; ec2 = 4; sc = SC_NEIG; break;	\
		case 5:																		\
			ia = 3; ib = 1; ic = 0, id = 2; ec1 = 5; ec2 = 3; sc = SC_NEIG; break;	\
		}																			\
		break;																		\
	case 4:																			\
		switch ((maec))																\
		{																			\
		case 5:																		\
			ia = 1; ib = 2; ic = 0, id = 3; ec1 = 4; ec2 = 5; sc = SC_NEIG; break;	\
		}																			\
		break;																		\
	}

/* 
 * ������� DFC(Decoding Face Codes) 
 *    ��֪	��BCD��ֵfc
 *    δ֪����Ԫ4���ڵ����ia ib ic id
 *    Լ��:	abcd��һ���������ַ���������� 
 * DFC: Decoding Face Codes
 *   Known:	fc, code of face BCD
 * Unknown:	ia, ib, ic, id, codes of 4 nodes
 *    s.t.:	abcd is a positive tet. 
 */
#define DFC(fc, ia, ib, ic, id)			\
	ia = (fc);							\
	switch ((fc))						\
	{									\
	case 0:								\
		ib = 1; ic = 2; id = 3; break;	\
	case 1:								\
		ib = 2; ic = 0; id = 3; break;	\
	case 2:								\
		ib = 0; ic = 1; id = 3; break;	\
	case 3:								\
		ib = 0; ic = 2; id = 1; break;	\
	}

/* 
 * ������������ DEFC(Decoding Edge & Face Codes) 
 *    ��֪	ec ��AD��ֵ
            fc ��BCD��ֵ
 *    δ֪����Ԫ4���ڵ����ia ib ic id
 *    Լ��:	abcd��һ���������ַ���������� 
 * DEFC: Decoding Edge & Face Codes
 *   Known:	ec, code of edge AD
            fc, code of face BCD
 * Unknown:	ia, ib, ic, id, codes of 4 nodes
 *    s.t.:	abcd is a positive tet. 
 */
#define DEFC(ec, fc, ia, ib, ic, id)		\
	ia = (fc);								\
	switch ((fc))							\
	{										\
	case 0:									\
		switch ((ec))						\
		{									\
		case 0:								\
			ib = 2; ic = 3; id = 1; break;	\
		case 1:								\
			ib = 3; ic = 1; id = 2; break;	\
		case 2:								\
			ib = 1; ic = 2; id = 3; break;	\
		}									\
		break;								\
	case 1:									\
		switch ((ec))						\
		{									\
		case 0:								\
			ib = 3; ic = 2; id = 0; break;	\
		case 3:								\
			ib = 0; ic = 3; id = 2; break;	\
		case 4:								\
			ib = 2; ic = 0; id = 3; break;	\
		}									\
		break;								\
	case 2:									\
		switch ((ec))						\
		{									\
		case 1:								\
			ib = 1; ic = 3; id = 0; break;	\
		case 3:								\
			ib = 3; ic = 0; id = 1; break;	\
		case 5:								\
			ib = 0; ic = 1; id = 3; break;	\
		}									\
		break;								\
	case 3:									\
		switch ((ec))						\
		{									\
		case 2:								\
			ib = 2; ic = 1; id = 0; break;	\
		case 4:								\
			ib = 0; ic = 2; id = 1; break;	\
		case 5:								\
			ib = 1; ic = 0; id = 2; break;	\
		}									\
		break;								\
	}

/* 
 * ˫���������� DDFC(Decoding Double Face Codes) 
 *    ��֪	mifc ��ACD��ֵ
            mafc ��BCD��ֵ
 *    δ֪����Ԫ4���ڵ����ia ib ic id
 *    Լ��:	abcd��һ���������ַ���������� 
            mifc < mafc
 * DDFC: Decoding Double Face Codes
 *   Known:	mifc, code of face ACD
            mafc, code of face BCD
 * Unknown:	ia, ib, ic, id, codes of 4 nodes
 *    s.t.:	abcd is a positive tet. 
            mifc < mafc
 */
#define DDFC(mifc, mafc, ia, ib, ic, id) \
	ia = (mafc);						 \
	ic = (mifc);						 \
	switch ((mifc))						 \
	{								\
	case 0:							\
		switch ((mafc))				\
		{							\
		case 1:						\
			ib = 2; id = 3; break;	\
		case 2:						\
			ib = 3; id = 1; break;	\
		case 3:						\
			ib = 1; id = 2; break;	\
		}							\
		break;						\
	case 1:							\
		switch ((mafc))				\
		{							\
		case 2:						\
			ib = 0; id = 3; break;	\
		case 3:						\
			ib = 2; id = 0; break;	\
		}							\
		break;						\
	case 2:							\
		switch ((mafc))				\
		{							\
		case 3:						\
			ib = 0; id = 1; break;	\
		}							\
		break;						\
	}

/*
 * �������˵���빹������ֵ
 * CEC: Construct Edge Codes from codes of two end nodes
 */
#define CEC(mie, mae, ec)			\
	if ((mie) == 0)					\
		ec = (mie) + (mae) - 1;	\
	else							\
		ec = (mie) + (mae);		

/*
 * ������룬��ȡ�����˵�ĵ���n1 n2
 * DECON: Decoding Edge Codes to Obtain Nodes
 */
#define DECON(ec, n1, n2)		\
	switch ((ec))				\
	{							\
	case 0:						\
		n1 = 0; n2 = 1; break;	\
	case 1:						\
		n1 = 0; n2 = 2; break;	\
	case 2:						\
		n1 = 0; n2 = 3; break;	\
	case 3:						\
		n1 = 1; n2 = 2; break;	\
	case 4:						\
		n1 = 1; n2 = 3; break;	\
	case 5:						\
		n1 = 2; n2 = 3; break;	\
	}

/* 
 * ����i1, i2, i3, i4��һ�������嵥Ԫ���ĸ��ڵ����ֵ��
 * i2->i3->i4��һ��������Ƭ, i1->i2->i3->i4��һ�������������,
 * ��i1, i2, i3, i4��ֵ�ƶ��������������⻹�������ڲ�
 * �����_FAC_RIG_INNER�������Ƭ�����ַ���ָ�������ⲿ������ָ�������ڲ�
 * i1, i2, i3 and i4 are codes of four nodes of a tet.,
 * i2->i3->i4 is a facet, i1->i2->i3->i4 is a tet. with positive volume,
 * If the values of i1, i2, i3, i4 are given, flag the element as outer or inner the domain
 * Notice:
 * if macro _FAC_RIG_INNER is defined, normal based on the right-hand rule of facet i2->i3->i4
 * directs inside; otherwise, it directs outside
 */
#ifdef _FAC_RIG_INNER
#define FLAGE(i1, i2, i3, i4, flag)						\
	flag = OUTER;										\
	switch ((i1))										\
	{													\
	case 0:												\
		if (((i2) == 1 && (i3) == 3 && (i4) == 2) ||	\
		    ((i2) == 2 && (i3) == 1 && (i4) == 3) ||	\
		    ((i2) == 3 && (i3) == 2 && (i4) == 1))		\
			flag = INNER;								\
		break;											\
	case 1:												\
		if (((i2) == 2 && (i3) == 3 && (i4) == 0) ||	\
			((i2) == 0 && (i3) == 2 && (i4) == 3) ||	\
			((i2) == 3 && (i3) == 0 && (i4) == 2))		\
			flag = INNER;								\
		break;											\
	case 2:												\
		if (((i2) == 0 && (i3) == 3 && (i4) == 1) ||	\
			((i2) == 1 && (i3) == 0 && (i4) == 3) ||	\
			((i2) == 3 && (i3) == 1 && (i4) == 0))		\
			flag = INNER;								\
		break;											\
	case 3:												\
		if (((i2) == 2 && (i3) == 0 && (i4) == 1) ||	\
			((i2) == 1 && (i3) == 2 && (i4) == 0) ||	\
			((i2) == 0 && (i3) == 1 && (i4) == 2))		\
			flag = INNER;								\
		break;											\
	}	
	
#else
#define FLAGE(i1, i2, i3, i4, flag)						\
	flag = INNER;										\
	switch ((i1))										\
	{													\
	case 0:												\
		if (((i2) == 1 && (i3) == 3 && (i4) == 2) ||	\
		    ((i2) == 2 && (i3) == 1 && (i4) == 3) ||	\
		    ((i2) == 3 && (i3) == 2 && (i4) == 1))		\
			flag = OUTER;								\
		break;											\
	case 1:												\
		if (((i2) == 2 && (i3) == 3 && (i4) == 0) ||	\
			((i2) == 0 && (i3) == 2 && (i4) == 3) ||	\
			((i2) == 3 && (i3) == 0 && (i4) == 2))		\
			flag = OUTER;								\
		break;											\
	case 2:												\
		if (((i2) == 0 && (i3) == 3 && (i4) == 1) ||	\
			((i2) == 1 && (i3) == 0 && (i4) == 3) ||	\
			((i2) == 3 && (i3) == 1 && (i4) == 0))		\
			flag = OUTER;								\
		break;											\
	case 3:												\
		if (((i2) == 2 && (i3) == 0 && (i4) == 1) ||	\
			((i2) == 1 && (i3) == 2 && (i4) == 0) ||	\
			((i2) == 0 && (i3) == 1 && (i4) == 2))		\
			flag = OUTER;								\
		break;											\
	}		
#endif /*_FAC_RIG_INNER*/

/* 
 * ����i1, i2, i3, i4��һ�������嵥Ԫ���ĸ��ڵ����ֵ��
 * i2->i3->i4��һ��������Ƭ, i1->i2->i3->i4��һ�������������,
 * ��i1, i2, i3, i4��ֵ�ƶ��������������⻹�������ڲ�
 * �����_FAC_RIG_INNER�������Ƭ�����ַ���ָ�������ⲿ������ָ�������ڲ�
 * i1, i2, i3 and i4 are codes of four nodes of a tet.,
 * i2->i3->i4 is a facet, i1->i2->i3->i4 is a tet. with positive volume,
 * If the values of i1, i2, i3, i4 are given, flag the element as outer or inner the domain
 * Notice:
 * if macro _FAC_RIG_INNER is defined, normal based on the right-hand rule of facet i2->i3->i4
 * directs inside; otherwise, it directs outside
 */
#define ELEM_ORIENTATION(i1, i2, i3, i4, ort)						\
	ort = 1;										\
	switch ((i1))										\
	{													\
	case 0:												\
		if (((i2) == 1 && (i3) == 3 && (i4) == 2) ||	\
		    ((i2) == 2 && (i3) == 1 && (i4) == 3) ||	\
		    ((i2) == 3 && (i3) == 2 && (i4) == 1))		\
			ort = 0;								\
		break;											\
	case 1:												\
		if (((i2) == 2 && (i3) == 3 && (i4) == 0) ||	\
			((i2) == 0 && (i3) == 2 && (i4) == 3) ||	\
			((i2) == 3 && (i3) == 0 && (i4) == 2))		\
			ort = 0;								\
		break;											\
	case 2:												\
		if (((i2) == 0 && (i3) == 3 && (i4) == 1) ||	\
			((i2) == 1 && (i3) == 0 && (i4) == 3) ||	\
			((i2) == 3 && (i3) == 1 && (i4) == 0))		\
			ort = 0;								\
		break;											\
	case 3:												\
		if (((i2) == 2 && (i3) == 0 && (i4) == 1) ||	\
			((i2) == 1 && (i3) == 2 && (i4) == 0) ||	\
			((i2) == 0 && (i3) == 1 && (i4) == 2))		\
			ort = 0;								\
		break;											\
	}	
	

/*
 * �������ƽ��
 * square of distance between two points
 */
REAL squaDist(MYPOINT p1, MYPOINT p2);

/*
 * �����������ֵ
 * absolute value of area of a triangle 
 */
REAL areas(REAL xa, REAL ya, REAL za, 
		   REAL xb, REAL yb, REAL zb,
		   REAL xc, REAL yc, REAL zc);

/* 
 * �������ַ�����
 * normal of a face based on the right-hand rule
 */
int facNorm(MYPOINT p1, MYPOINT p2, MYPOINT p3, VECTOR norm);

/*
 * ��λ������
 * normalize a vector 
 */
int normVect(VECTOR v);

/*
 * �������
 * cross operation of two vectors
 */
int crossVect(VECTOR v1, VECTOR v2, VECTOR *v); /* v = v1^v2 */

/*
 * �Ƚ����������ķ��� compare directions of two vectors
 * -1 δ֪ undefined
 * 0  ���� of same directions
 * 1  ͬ�� of opposite directions
 */
int compVectDir(VECTOR v1, VECTOR v2);

/* 
 * ��ȡ�����������任ΪZ�ᣨ0��0��1������ı任����(4��4)
 * get the transform matrix to convert vector v to (0, 0, 1) (Z axis)
 */
int vectorToZ(Vector n, double m[16]);

/* 
 * p4���������p1p2p3�Ƿ���ӣ���p1p2p3p4������Ƿ�Ϊ��
 * check if p4 is visible to p1p2p3, i.e. the volume of p1p2p3p4 is positive 
 */
bool isVisible(MYPOINT p1, MYPOINT p2, MYPOINT p3, MYPOINT p4);

/* 
 * ia ib ic id �������嵥Ԫ���ĸ���ĵ���ֵ���ж������嵥Ԫ
 * A->B->C->D�Ƿ������(volume-positive)
 * ia, ib, ic and id are four node codes of four forming points of Tetrahedra 
 * ABCD, judge if the volume of ABCD is positive
 */
bool isVolumePositive(int ia, int ib, int ic, int id);

int generateRandArray(int *arr, int maxValue);

#endif /* __iso3d_utility_h__ */