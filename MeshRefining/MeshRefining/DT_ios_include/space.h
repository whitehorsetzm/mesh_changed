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

#ifndef __pdmg_space_h__
#define __pdmg_space_h__

#include "iso3d_define.h"

/* ��Դ point source*/
typedef struct PntSource
{
	MYPOINT cen;				/* ����λ�� center position */
	REAL ou_rad, in_rad;	/* �⾶ & �ھ� outer & inner radius */
	REAL space;             /* �ߴ�ֵ space value */
	MYPOINT bmin, bmax;		/* Ӱ������ box of the impacted region */
} PntSource;

/* ��Դ line source*/
typedef struct LinSource
{
	PntSource src[2];		/* ������Դ two point sources */
	MYPOINT bmin, bmax;		/* Ӱ������ box of the impacted region */
} LinSource;

/* ��Դ triangle source*/
typedef struct TriSource
{
	PntSource src[3];		/* ������Դ two point sources */
	MYPOINT bmin, bmax;		/* Ӱ������ box of the impacted region */

} TriSource;

typedef struct Source
{
	PntSource *pPntS;	/* һ���Դ a group of point sources */
	LinSource *pLinS;	/* һ����Դ a group of line sources */
	TriSource *pTriS;	/* һ����Դ a group of triangle sources */
	int nPntS;			/* ��Դ���� number of point sources */			
	int nLinS;			/* ��Դ���� number of line sources */	
	int nTriS;			/* ��Դ���� number of triangle sources */	
} Source;

/* �������� */
typedef struct BGMesh
{
	Elem *pElems;	 /* ��������Ԫ���� background element array*/
	Node *pNodes;	 /* ��������ڵ����� backgound nodes array */
	INTEGER nElems;	 /* ��������Ԫ��Ŀ number of background elements */
	INTEGER nNodes;  /* ��������ڵ���Ŀ number of background nodes */
} BGMesh;

/* ��ձ���������� free the object of background mesh */
int freeBGMesh(BGMesh *pBGMesh);

/* ���Դ���� free the object of source */
int freeSource(Source *pSource);
/*
 * ��ȡBA3�ļ�
 * read a .ba3 file
 */
int readBA3(const char *fname, BGMesh *pBGMesh, Source *pSource);

/*
 * дBA3�ļ�
 * write a .ba3 file
 */
int writeBA3(const char *fname, BGMesh *pBGMesh, Source *pSource);

#endif /*  __pdmg_space_h__ */
