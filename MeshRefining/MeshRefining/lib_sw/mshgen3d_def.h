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

 * Copyright reserved, 2006, 03, 08

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

#ifndef __dtiso3d_meshgen3d_def_h__

#define __dtiso3d_meshgen3d_def_h__



// #  ifdef __cplusplus

// extern "C" {

// #  endif /* __cplusplus */

int meshImproved(int nNodes,double *x,double *y,double *z, int nTetras,int *tetras, int nFacets, int *facets,
                 int &outNodes,double*outx,double*outy,double*outz,int &outNtetras,int *ouTetras);

int meshGen3D_memo(

/* ------------------------- input arguments --------------------------*/

	double		*pdSNX,	/* x coord. of boundary nodes */

	double		*pdSNY,	/* y coord. of boundary nodes */

	double		*pdSNZ,	/* z coord. of boundary nodes */

	int			nSN,		/* number of boundary nodes */

	int			*pnSFFm,	/* forming points of surface facets */

	int			*pnSFPt,	/* patches of surface facets */

	int         nSF,		/* number of facets */

	double		*pdBMNX,	/* x coord. of background mesh */

	double		*pdBMNY,	/* y coord. of background mesh */

	double		*pdBMNZ,	/* z coord. of background mesh */

	double		*pdBMNSpc,	/* space values of background mesh nodes */

	int			nBMN,		/* number of background mesh nodes */

	int			*pnBMEFm,	/* forming points of background mesh elements */

	int			*pnBMENg,	/* neighboring eles. of background mesh elements */

	int			nBME,		/* number of background mesh elements */

	double		*pdCX,		/* x coord. of center points of sources */

	double		*pdCY,		/* y coord. of center points of sources */

	double		*pdCZ,		/* z coord. of center points of sources */

	double		*pdIn,		/* inner radius of sources */

	double		*pdOu,		/* out radius of sources */

	double		*pdSp,		/* space values of sources */

	int			nPS,		/* number of point source */

	int			nLS,		/* number of line source */

	int			nTS,		/* number of triangle source */

/* ------------------------- output arguments -------------------------*/

	double	  **ppdMNX,		/* x coord. of 3D mesh nodes */

	double    **ppdMNY,		/* y coord. of 3D mesh nodes */

	double    **ppdMNZ,		/* z coord. of 3D mesh nodes */

	int        *pnMN,		/* number of 3D mesh nodes */

	int       **ppnMEFm,	/* forming points of 3D mesh elements */

	int       **ppnMENg,	/* neighboring eles. of 3D mesh elements */

	int        *pnME,		/* number of 3D mesh elements */

        int        **ppnPrt	/* parents of surface facets */
	);



// #  ifdef __cplusplus

// }

// #  endif /* __cplusplus */



void printRmtTime();

void printRcvTime();



#endif /* __dtiso3d_meshgen3d_def_h__ */
