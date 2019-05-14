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

#ifndef __iso3d_utility_h__
#define __iso3d_utility_h__

#include "iso3d_define.h"
#include "vector.h"
/* **********************************************************
 * 单元实体(节点/边/面)的编码/解码
 * coding & decoding of elemental entities (nodes/edges/faces)
 * *********************************************************/

/* 
 * 点码解码 DNC(Decoding Node Codes) 
 *    已知：A点编码nc
 *    未知：单元的另外3个节点编码ib ic id
 *    约束:	abcd是一个满足右手法则的四面体 
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
 * 双点码解码 DDNC(Decoding Double Node Codes) 
 *    已知：C D点编码ic id
 *    未知：单元的另外2个节点编码ia ib 
 *    约束:	abcd是一个满足右手法则的四面体 
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
 * 边码解码 DEC(Decoding Edge Codes) 
 *    已知：边BD码值ib id
 *    未知：单元的另外2个节点编码ia ic 
 *    约束:	abcd是一个满足右手法则的四面体 
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
 * 双边码解码 DDEC(Decoding Double Edge Codes) 
 *    已知：边码小值&编码大值 miec & maec
 *    未知：单元的4个节点码值ia ib ic id
            两条边排序后的码值 ec1 ec2
			两条边的关系 相邻/相对
 *    约束:	abcd是一个满足右手法则的四面体;
            如果是邻边 ec1 = AD ec2 = BD
			如果是对边 ec1 = AD ec2 = BC
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
 * 面码解码 DFC(Decoding Face Codes) 
 *    已知	面BCD码值fc
 *    未知：单元4个节点编码ia ib ic id
 *    约束:	abcd是一个满足右手法则的四面体 
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
 * 边面组合码解码 DEFC(Decoding Edge & Face Codes) 
 *    已知	ec 边AD码值
            fc 面BCD码值
 *    未知：单元4个节点编码ia ib ic id
 *    约束:	abcd是一个满足右手法则的四面体 
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
 * 双面组合码解码 DDFC(Decoding Double Face Codes) 
 *    已知	mifc 面ACD码值
            mafc 面BCD码值
 *    未知：单元4个节点编码ia ib ic id
 *    约束:	abcd是一个满足右手法则的四面体 
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
 * 从两个端点点码构建编码值
 * CEC: Construct Edge Codes from codes of two end nodes
 */
#define CEC(mie, mae, ec)			\
	if ((mie) == 0)					\
		ec = (mie) + (mae) - 1;	\
	else							\
		ec = (mie) + (mae);		

/*
 * 解码编码，获取两个端点的点码n1 n2
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
 * 假设i1, i2, i3, i4是一个四面体单元的四个节点的码值，
 * i2->i3->i4是一个表面面片, i1->i2->i3->i4是一个正体积四面体,
 * 由i1, i2, i3, i4的值推断四面体在区域外还是区域内部
 * 如果宏_FAC_RIG_INNER则表明面片的右手法向指向区域外部，否则，指向区域内部
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
 * 假设i1, i2, i3, i4是一个四面体单元的四个节点的码值，
 * i2->i3->i4是一个表面面片, i1->i2->i3->i4是一个正体积四面体,
 * 由i1, i2, i3, i4的值推断四面体在区域外还是区域内部
 * 如果宏_FAC_RIG_INNER则表明面片的右手法向指向区域外部，否则，指向区域内部
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
 * 两点距离平方
 * square of distance between two points
 */
REAL squaDist(MYPOINT p1, MYPOINT p2);

/*
 * 三点面积绝对值
 * absolute value of area of a triangle 
 */
REAL areas(REAL xa, REAL ya, REAL za, 
		   REAL xb, REAL yb, REAL zb,
		   REAL xc, REAL yc, REAL zc);

/* 
 * 表面右手法向量
 * normal of a face based on the right-hand rule
 */
int facNorm(MYPOINT p1, MYPOINT p2, MYPOINT p3, VECTOR norm);

/*
 * 单位化向量
 * normalize a vector 
 */
int normVect(VECTOR v);

/*
 * 向量叉乘
 * cross operation of two vectors
 */
int crossVect(VECTOR v1, VECTOR v2, VECTOR *v); /* v = v1^v2 */

/*
 * 比较两个向量的方向 compare directions of two vectors
 * -1 未知 undefined
 * 0  反向 of same directions
 * 1  同向 of opposite directions
 */
int compVectDir(VECTOR v1, VECTOR v2);

/* 
 * 获取将给定向量变换为Z轴（0，0，1）所需的变换矩阵(4×4)
 * get the transform matrix to convert vector v to (0, 0, 1) (Z axis)
 */
int vectorToZ(Vector n, double m[16]);

/* 
 * p4点对三角面p1p2p3是否可视，即p1p2p3p4的体积是否为正
 * check if p4 is visible to p1p2p3, i.e. the volume of p1p2p3p4 is positive 
 */
bool isVisible(MYPOINT p1, MYPOINT p2, MYPOINT p3, MYPOINT p4);

/* 
 * ia ib ic id 是四面体单元的四个点的点码值，判断四面体单元
 * A->B->C->D是否有体积(volume-positive)
 * ia, ib, ic and id are four node codes of four forming points of Tetrahedra 
 * ABCD, judge if the volume of ABCD is positive
 */
bool isVolumePositive(int ia, int ib, int ic, int id);

int generateRandArray(int *arr, int maxValue);

#endif /* __iso3d_utility_h__ */