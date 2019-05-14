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

#ifndef __pdmg_space_h__
#define __pdmg_space_h__

#include "iso3d_define.h"

/* 点源 point source*/
typedef struct PntSource
{
	MYPOINT cen;				/* 中心位置 center position */
	REAL ou_rad, in_rad;	/* 外径 & 内径 outer & inner radius */
	REAL space;             /* 尺寸值 space value */
	MYPOINT bmin, bmax;		/* 影响区域 box of the impacted region */
} PntSource;

/* 线源 line source*/
typedef struct LinSource
{
	PntSource src[2];		/* 两个点源 two point sources */
	MYPOINT bmin, bmax;		/* 影响区域 box of the impacted region */
} LinSource;

/* 面源 triangle source*/
typedef struct TriSource
{
	PntSource src[3];		/* 三个点源 two point sources */
	MYPOINT bmin, bmax;		/* 影响区域 box of the impacted region */

} TriSource;

typedef struct Source
{
	PntSource *pPntS;	/* 一组点源 a group of point sources */
	LinSource *pLinS;	/* 一组线源 a group of line sources */
	TriSource *pTriS;	/* 一组面源 a group of triangle sources */
	int nPntS;			/* 点源个数 number of point sources */			
	int nLinS;			/* 线源个数 number of line sources */	
	int nTriS;			/* 面源个数 number of triangle sources */	
} Source;

/* 背景网格 */
typedef struct BGMesh
{
	Elem *pElems;	 /* 背景网格单元数组 background element array*/
	Node *pNodes;	 /* 背景网格节点数组 backgound nodes array */
	INTEGER nElems;	 /* 背景网格单元数目 number of background elements */
	INTEGER nNodes;  /* 背景网格节点数目 number of background nodes */
} BGMesh;

/* 清空背景网格对象 free the object of background mesh */
int freeBGMesh(BGMesh *pBGMesh);

/* 清空源对象 free the object of source */
int freeSource(Source *pSource);
/*
 * 读取BA3文件
 * read a .ba3 file
 */
int readBA3(const char *fname, BGMesh *pBGMesh, Source *pSource);

/*
 * 写BA3文件
 * write a .ba3 file
 */
int writeBA3(const char *fname, BGMesh *pBGMesh, Source *pSource);

#endif /*  __pdmg_space_h__ */
