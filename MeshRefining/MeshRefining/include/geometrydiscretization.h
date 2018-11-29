/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 多学科应用模拟的赋能环境
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2007年10月11日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 11, 2007
 * 
 * 联系方式 (For further information, please conctact)
 *   电话 (Tel)：+86-571-87953165
 *   传真 (Fax)：+86-571-87953167
 *   邮箱 (Mail)：chenjj@zju.edu.cn
 *
 * 文件名称 (File Name)：geometrydiscretization.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对0维~3维几何体的离散方法和离散数据结构
 *     Define a set of data structures & discretization 
 *     algorithms for zero ~ three dimensions geometries.
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __eemas_geometrydiscretization_h__
#define __eemas_geometrydiscretization_h__

# include <stdio.h>
# include "global.h"
# include "entityiterator.h"
# include "geomhandle.h"
# include <math.h>
# include "spacing.h"

// #include "SpatialTree.h"
// 
// QuadTree * quadTree;

#include <map>

typedef EntityHandle ElemHandle;
typedef EntityHandle GPntHandle;

typedef std::vector<ElemHandle> ElemHandleArray;
typedef std::vector<GPntHandle> GPntHandleArray;

typedef ArrayIterator ElemArrayIterator;
typedef ArrayIterator GPntArrayIterator;

int get_grid_element_index(ElemHandle eh);
int get_grid_point_index(GPntHandle ph);
int delete_grid(
	void* geom,				//the geometry on which the grid is
	int type				//indicate the type of the geom: curve,face,domain
							//2: curve; 3: face; 4: domain
	);
int delete_grid_point(GPntHandle p_hd);
int delete_grid_elem(ElemHandle ele_hd);
int update_grid_data();		//update the grid data according the refcount

//extern SurfBndBKGMesh * g_bkgmesh;
//extern std::map<GBLoop *, SurfBndBKGMesh *> g_surfBndBKGMeshs;

//extern int* ibnd;
//extern int ibnd_size;
/* -------------------------------------------------------------------------
 * 离散数据类声明 declaration of class DiscretizationData
 * 设计意图: 
 *	   1. 它保存了几何体离散后所获得的离散信息 
 *     2. 它只对GeometryDiscretization类及其集成类开放其写接口；
          对其它类而言，只能读取数据
       3. 缺省的数据包括一个单元数组和一个节点数组
       4. 针对GeometryDiscretization的派生类，用户可以根据需求
	      设计DiscretizationData的派生类
 * Design Considerations:
 *     1. It saves the discretization data for a geometry
 *     2. It is readable only for GeometryDiscretization object and its 
 *        children objects.
 *     3. Default memebers include an element array and a grid point array
 *     4. Derived classes could be designed for child classes of 
          GeometryDiscretization
 * -------------------------------------------------------------------------*/
class DiscretizationData
{
public:
	DiscretizationData()
		:elem_handle_iterator(&elem_handle_array), 
		gpnt_handle_iterator(&gpnt_handle_array)
	{
		meshPolicyType = -1;
		elem_handle_array.clear();
		gpnt_handle_array.clear();
	}
	~DiscretizationData() {}
	
	/* 遍历所有单元 travel all elements */
	ElemArrayIterator& get_elem_iterator() {return elem_handle_iterator;}
	
	/* 遍历所有网格节点 travel all nodes */
	GPntArrayIterator& get_gpnt_iterator() {return gpnt_handle_iterator;}
	
	/* 得到网格点索引数组 */
	GPntHandleArray* get_gpnt_array(){return &gpnt_handle_array;}

	GPntHandleArray& get_gpnt_array_ref(){return gpnt_handle_array;}
	
	/* 得到网格单元索引数组 */
	ElemHandleArray* get_elem_array(){return &elem_handle_array;}

	ElemHandleArray& get_elem_array_ref(){return elem_handle_array;}
	
	
	int number_of_elems() { return elem_handle_array.size(); }
	int number_of_gpnts() { return gpnt_handle_array.size(); }

	
	virtual void set_mesh_policy(int mtype){meshPolicyType = mtype;}

	virtual int get_mesh_policy(){return meshPolicyType;}

	void set_size(double sz) {size = sz;} // added [1/28/2009 leon]
	double get_size() {return size;} // added [1/28/2009 leon]

protected:
	
	/* 网格化使用的算法 */
	int meshPolicyType;

	double size; 
	// size of the discretized object. e.g. the spacing value for meshing surface. [1/28/2009 leon]

private:
/* 几何体离散化后得到的单元数组 
	* array of generated elements */
	ElemHandleArray elem_handle_array;
	/* 几何体离散化后得到的节点数组 */
	GPntHandleArray gpnt_handle_array;
	
	ElemArrayIterator elem_handle_iterator;
	GPntArrayIterator gpnt_handle_iterator;
	
	friend class GeometryDiscretization;
};

class GeometryDiscretization
{
public:
	GeometryDiscretization(GlobalSpacing *glob_sp = NULL)
	{
		geom_dim = DIM_ZERO;
		geom_handle = GEOMETRY_NULL;
		global_spacing_enabled = false;
		per_curve_enabled = false;

		global_space = glob_sp;
		d_data = NULL;
	}
	~GeometryDiscretization() 
	{
	}
	
	int get_dimensions() {return geom_dim;}
	/* 
	   设置几何体的维度
	   set the dimensions
	   注意：一个子类往往其维度是确定的，或
	   仅仅可以进行部分设定，因此，将函数设定
	   为虚函数，当返回false时表示不接受用户设定。
	   Notice: Not all discretization child objects allow the users to set
	   the dimensions. This function could be overloaded by child objects 
	   to return false when an unaccepted calling is performed by the users.
	 */
	virtual bool set_dimensions(int dim) {geom_dim = dim; return true;}

	GeometryHandle get_geometry() {return geom_handle;}

	virtual bool set_geometry(GeometryHandle gh) = 0;
        //annotation by zhvliu
    //virtual int discretization() = 0;
	
	/* 检查密度控制开关
	   check the spacing control switches */
	bool is_global_spacing_enabled() {return global_spacing_enabled;}
	bool is_per_curve_enabled() {return per_curve_enabled;}

	/* 
	   设置密度控制开关
	   set the spacing control switches 
	   注意：并不是没一类离散化方法都允许用户设定部分或全部密度控制方式，
	   因此，将其设定为虚函数，当返回false时表示不接受用户设定。
	   Notice: Not all discretization child objects allow the users to set
	   some parts or all of spacing control switches. These functions could 
	   be overloaded by child objects to return false when an unaccepted calling
	   is performed by the users.
	 */
	virtual bool enable_global_spacing(bool is_on) 
	{
		global_spacing_enabled = is_on; 
		return true;
	}

	virtual bool enable_per_curve(bool is_on)
	{
		per_curve_enabled = is_on;
		return true;
	}

	DiscretizationData *get_discretization_data() {return d_data;}
	
protected:
	ElemHandleArray *get_elem_array() 
	{
		if (d_data)
		{
			return &d_data->elem_handle_array;
		}
		return NULL;
	}
	GPntHandleArray *get_gpnt_array()
	{
		if (d_data)
		{
			return &d_data->gpnt_handle_array;
		}
		return NULL;
	}

protected:
	/* 问题维度 problem dimensions */
	int geom_dim;
	
	/* 几何体句柄 geometry handle*/
	GeometryHandle geom_handle; 
		
	/* 是否启动全局密度控制尺寸，缺省值为true 
	 * toggle switch for global spacing, default value is true */
	bool global_spacing_enabled;

	/* 逐个曲线指定密度控制策略的方式是否启动，缺省值为false 
	 * toggle switch for per curve segmentation way, default value is false */
	bool per_curve_enabled;

	/* 全局密度控制信息 */
	GlobalSpacing *global_space;
	
	/* 结果数据　result data */
	DiscretizationData *d_data;
	

};

class PointDiscretization: public GeometryDiscretization
{
public:
	PointDiscretization(GlobalSpacing *glob_sp = NULL);
	~PointDiscretization()
	{
		if (geom_handle != POINT_NULL)
		{
			::del_point_handle(geom_handle);
		}
	}

	virtual bool set_geometry(GeometryHandle gh);
    //annotation by zhvliu
//	virtual int discretization();
protected:
};

class CurveDiscretization: public GeometryDiscretization
{
public:
	CurveDiscretization(GlobalSpacing *glob_sp = NULL);
	~CurveDiscretization()
	{
		if (geom_handle != CURVE_NULL)
		{
			::del_curve_handle(geom_handle);
		}
	}

	virtual bool set_geometry(GeometryHandle gh);
	virtual int discretization();

	//added 2012/6/20
	int discretization_sample();
    int discretization_adaptive();
	// added [12/26/2008 ly]
	void setAdaptive(bool adaptive) {
		this->adaptive = adaptive;
	}
	// added [1/9/2009 ly]
	void setSampling(bool sampling) {
		this->sampling = sampling;
	}
	bool smoothSamples(); // added [2/8/2009 leon]
protected:
private:
	// added [12/26/2008 ly]
	bool adaptive;
	// added [1/9/2009 ly]
	bool sampling;
};

class FaceDiscretization: public GeometryDiscretization
{
public:
	FaceDiscretization(GlobalSpacing *glob_sp = NULL);
	~FaceDiscretization() 		
	{
		if (geom_handle != FACE_NULL)
		{
			::del_face_handle(geom_handle);
		}
	}

	virtual bool set_geometry(GeometryHandle gh);
	virtual int discretization();
	//zlj 2012/6/19
	// 采样离散
	int discretization_sample_2d();
	//重新离散曲线
	int discretization_adaptive_2d();
	//zlj 2012/9/2
	int discretization_sample_3d();
	int discretization_adaptive_3d();
     
	int get_boundary_info(int& nbp, int &nbe);
	
    //virtual int construct_boundary_info(GBLoop* face, int& nbp,double*& conp, int*& lpnp);//annotation by zhvliu
	// added [12/26/2008 ly]
	void setAdaptive(bool adaptive) {
		this->adaptive = adaptive;
	}

	// added [1/9/2009 ly]
	void setSampling(bool sampling) {
		this->sampling = sampling;
	}

	bool smoothSamples(); // added [2/8/2009 leon]
protected:
/* -------------------------------------------------------------------
* 在边界网格生成之后，将其节点信息加入利索数据对象的节点数组中
* 切记：
*     不能假定边界网格的节点在数组中是连续排列的，更不能假定
*     它就在刚生成的网格节点段中
* 返回值：
*     错误代码
* -----------------------------------------------------------------*/
	virtual int attach_boundary_nodes();
//private:
public:
	// Extracted [1/9/2009 ly]
	void prepareLocalBound(GBLoop * loop, 
		int mapsize, int * g_to_l, int * l_to_g1, 
		int prob_dim, int *ibnd, double *bcoord, int *nbp, int *nbe);
	
	// Extracted [1/8/2009 ly]
	int discretizeCurvesAroundFace(GBFace * face, CurveDiscretization &cv_disc, int *max_nbe, int *max_nbp, int sample);

private:
	// added [12/26/2008 ly]
	int adaptive;
	// added [1/9/2009 ly]
	int sampling;
};

class DomainDiscretization: public GeometryDiscretization
{
public:
	DomainDiscretization(GlobalSpacing *glob_sp = NULL);
	~DomainDiscretization()
	{
		if (geom_handle != DOMAIN_NULL)
		{
			::del_domain_handle(geom_handle);
		}
	}

	virtual bool set_geometry(GeometryHandle gh);
    //annotation by zhvliu
//	virtual int discretization();
	//确定某个区域中的面的引用的正负值
	int bndFaceOrientation(int*& surf_ort);
	
protected:
	virtual int attach_boundary_nodes();
};

/* -----------------------------------------------------------------------------
 * 一些辅助函数 some help functions
 * ----------------------------------------------------------------------------*/
/* 创建GBPoint对象(如果没有则它)
 * create a GBPoint object, if no such a object */
GBPoint *create_point(int iv);

int create_grid_element(int conn[], int nconn,void* pgeom);
int create_grid_element(int gpntA, int gpntB, int gpntC, int gpntD,void* pface);

/* 创建网格点，将其加入网格点数组，返回数组索引
 * create a grid point, insert it into the grid point arrya,
   and return its index
 */
int create_grid_point(double x, double y, double z,int tp);

/* ---------------------------------------------------------------------------------
* 将系统中定义的2D密度控制信息以基本类型的数组形式获取，供网格化函数调用
* 参数：
*		nbmn				背景网格顶点数目
*	    bmn_info			背景网格顶点信息，纬度为3：[x y spacing]
*		nbme				背景网格单元数目
*		bm_topu				背景网格单元拓扑信息，前3维为顶点编号，后3维为相邻关系
*		nsrc[2]				点源数目, 线源数目
*		src_info			点源、线源信息，每个点源包含5维信息：[x y spacing r1 r2]; 
*		                    每个线源包含2个点源
* 返回值: 
*      错误代码
*
* 描述：
*
* -------------------------------------------------------------------------------*/
int construct_primal_spacing_2D(
			double *gvalue,		/* 全局密度尺寸 */
			int *nbmn,			/* 背景网格顶点数目 */
			double **bmn_info,	/* 背景网格顶点信息，纬度为3：[x y spacing] */
			int *nbme,			/* 背景网格单元数目 */
			int **bm_topu,		/* 背景网格单元拓扑信息，前3维为顶点编号，后3维为相邻关系 */
			int nsrc[2],		/* 点源数目, 线源数目 */
			double **src_info	/* 点源、线源信息，每个点源包含5维信息：
									[x y spacing r1 r2]; 每个线源包含2个点源 */
			);

/* ---------------------------------------------------------------------------------
* 将系统中定义的3D密度控制信息以基本类型的数组形式获取，供网格化函数调用
* 参数：
*		nbmn				背景网格顶点数目
*	    bmn_info			背景网格顶点信息，纬度为4：[x y z spacing]
*		nbme				背景网格单元数目
*		bm_topu				背景网格单元拓扑信息，前4维为顶点编号，后4维为相邻关系
*		nsrc[3]				点源数目, 线源数目, 面源数目
*		src_info			点源、线源、面源信息，每个点源包含6维信息：
*							[x y z spacing r1 r2]; 每个线源包含2个点源;
*							每个面源包含3个点源;
* 返回值: 
*      错误代码
*
* 描述：
*
* -------------------------------------------------------------------------------*/
int construct_primal_spacing_3D(
			double *gvalue,		/* 全局密度尺寸 */
			int *nbmn,			/* 背景网格顶点数目 */
			double **bmn_info,	/* 背景网格顶点信息，纬度为4：[x y z spacing] */
			int *nbme,			/* 背景网格单元数目 */
			int **bm_topu,		/* 背景网格单元拓扑信息，前4维为顶点编号，后4维为相邻关系 */
			int nsrc[3],		/* 点源数目, 线源数目,  面源数目*/
			double **src_info	/* 点源、线源、面源信息，每个点源包含6维信息：
									[x y z spacing r1 r2]; 每个线源包含2个点源;
									每个面源包含3个点源; */
			);

//typedef std::vector<CurveDiscretization>		CurveDiscretizationArray;
//typedef std::vector<CurveDiscretization>::iterator		CurveDiscretizationArrayIterator;
//extern CurveDiscretizationArray					curveDiscretizationArray;
//extern CurveDiscretizationArrayIterator			curve_discretization_iterator;
//bool readFr2(const char* fname, int *nbp, double **lcoord, int **ibnd);
//bool readBa2(const char* fname, int *nbmn, double **nd_info, int *nbme, int **bm_topu, int nsrc[2], double **src_info);
//bool writePl2(const char* fname, int nbp, int ngp, double *g_coord, int nel, int *g_topu, int *ibnd, int *prt_arr);

/* ---------------------------------------------------------------------------------
* 计算出当前各类网格单元的数量，包括曲线网格单元、面网格单元和体网格单元
*
* 返回值: 
*      错误代码
*
* 描述：
*
* -------------------------------------------------------------------------------*/
int calculate_grid_number();


//for debug
void write_out_info(int& nbp, double*& lcoord, int& ibnd_size, int*& ibnd);

#endif /* __eemas_geometrydiscretization_h__ */
