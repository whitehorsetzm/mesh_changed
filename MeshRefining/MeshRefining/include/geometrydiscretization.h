/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2007��10��11��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 11, 2007
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953165
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��geometrydiscretization.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ�����0ά~3ά���������ɢ��������ɢ���ݽṹ
 *     Define a set of data structures & discretization 
 *     algorithms for zero ~ three dimensions geometries.
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
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
 * ��ɢ���������� declaration of class DiscretizationData
 * �����ͼ: 
 *	   1. �������˼�������ɢ������õ���ɢ��Ϣ 
 *     2. ��ֻ��GeometryDiscretization�༰�伯���࿪����д�ӿڣ�
          ����������ԣ�ֻ�ܶ�ȡ����
       3. ȱʡ�����ݰ���һ����Ԫ�����һ���ڵ�����
       4. ���GeometryDiscretization�������࣬�û����Ը�������
	      ���DiscretizationData��������
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
	
	/* �������е�Ԫ travel all elements */
	ElemArrayIterator& get_elem_iterator() {return elem_handle_iterator;}
	
	/* ������������ڵ� travel all nodes */
	GPntArrayIterator& get_gpnt_iterator() {return gpnt_handle_iterator;}
	
	/* �õ�������������� */
	GPntHandleArray* get_gpnt_array(){return &gpnt_handle_array;}

	GPntHandleArray& get_gpnt_array_ref(){return gpnt_handle_array;}
	
	/* �õ�����Ԫ�������� */
	ElemHandleArray* get_elem_array(){return &elem_handle_array;}

	ElemHandleArray& get_elem_array_ref(){return elem_handle_array;}
	
	
	int number_of_elems() { return elem_handle_array.size(); }
	int number_of_gpnts() { return gpnt_handle_array.size(); }

	
	virtual void set_mesh_policy(int mtype){meshPolicyType = mtype;}

	virtual int get_mesh_policy(){return meshPolicyType;}

	void set_size(double sz) {size = sz;} // added [1/28/2009 leon]
	double get_size() {return size;} // added [1/28/2009 leon]

protected:
	
	/* ����ʹ�õ��㷨 */
	int meshPolicyType;

	double size; 
	// size of the discretized object. e.g. the spacing value for meshing surface. [1/28/2009 leon]

private:
/* ��������ɢ����õ��ĵ�Ԫ���� 
	* array of generated elements */
	ElemHandleArray elem_handle_array;
	/* ��������ɢ����õ��Ľڵ����� */
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
	   ���ü������ά��
	   set the dimensions
	   ע�⣺һ������������ά����ȷ���ģ���
	   �������Խ��в����趨����ˣ��������趨
	   Ϊ�麯����������falseʱ��ʾ�������û��趨��
	   Notice: Not all discretization child objects allow the users to set
	   the dimensions. This function could be overloaded by child objects 
	   to return false when an unaccepted calling is performed by the users.
	 */
	virtual bool set_dimensions(int dim) {geom_dim = dim; return true;}

	GeometryHandle get_geometry() {return geom_handle;}

	virtual bool set_geometry(GeometryHandle gh) = 0;
        //annotation by zhvliu
    //virtual int discretization() = 0;
	
	/* ����ܶȿ��ƿ���
	   check the spacing control switches */
	bool is_global_spacing_enabled() {return global_spacing_enabled;}
	bool is_per_curve_enabled() {return per_curve_enabled;}

	/* 
	   �����ܶȿ��ƿ���
	   set the spacing control switches 
	   ע�⣺������ûһ����ɢ�������������û��趨���ֻ�ȫ���ܶȿ��Ʒ�ʽ��
	   ��ˣ������趨Ϊ�麯����������falseʱ��ʾ�������û��趨��
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
	/* ����ά�� problem dimensions */
	int geom_dim;
	
	/* �������� geometry handle*/
	GeometryHandle geom_handle; 
		
	/* �Ƿ�����ȫ���ܶȿ��Ƴߴ磬ȱʡֵΪtrue 
	 * toggle switch for global spacing, default value is true */
	bool global_spacing_enabled;

	/* �������ָ���ܶȿ��Ʋ��Եķ�ʽ�Ƿ�������ȱʡֵΪfalse 
	 * toggle switch for per curve segmentation way, default value is false */
	bool per_curve_enabled;

	/* ȫ���ܶȿ�����Ϣ */
	GlobalSpacing *global_space;
	
	/* ������ݡ�result data */
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
	// ������ɢ
	int discretization_sample_2d();
	//������ɢ����
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
* �ڱ߽���������֮�󣬽���ڵ���Ϣ�����������ݶ���Ľڵ�������
* �мǣ�
*     ���ܼٶ��߽�����Ľڵ������������������еģ������ܼٶ�
*     �����ڸ����ɵ�����ڵ����
* ����ֵ��
*     �������
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
	//ȷ��ĳ�������е�������õ�����ֵ
	int bndFaceOrientation(int*& surf_ort);
	
protected:
	virtual int attach_boundary_nodes();
};

/* -----------------------------------------------------------------------------
 * һЩ�������� some help functions
 * ----------------------------------------------------------------------------*/
/* ����GBPoint����(���û������)
 * create a GBPoint object, if no such a object */
GBPoint *create_point(int iv);

int create_grid_element(int conn[], int nconn,void* pgeom);
int create_grid_element(int gpntA, int gpntB, int gpntC, int gpntD,void* pface);

/* ��������㣬���������������飬������������
 * create a grid point, insert it into the grid point arrya,
   and return its index
 */
int create_grid_point(double x, double y, double z,int tp);

/* ---------------------------------------------------------------------------------
* ��ϵͳ�ж����2D�ܶȿ�����Ϣ�Ի������͵�������ʽ��ȡ�������񻯺�������
* ������
*		nbmn				�������񶥵���Ŀ
*	    bmn_info			�������񶥵���Ϣ��γ��Ϊ3��[x y spacing]
*		nbme				��������Ԫ��Ŀ
*		bm_topu				��������Ԫ������Ϣ��ǰ3άΪ�����ţ���3άΪ���ڹ�ϵ
*		nsrc[2]				��Դ��Ŀ, ��Դ��Ŀ
*		src_info			��Դ����Դ��Ϣ��ÿ����Դ����5ά��Ϣ��[x y spacing r1 r2]; 
*		                    ÿ����Դ����2����Դ
* ����ֵ: 
*      �������
*
* ������
*
* -------------------------------------------------------------------------------*/
int construct_primal_spacing_2D(
			double *gvalue,		/* ȫ���ܶȳߴ� */
			int *nbmn,			/* �������񶥵���Ŀ */
			double **bmn_info,	/* �������񶥵���Ϣ��γ��Ϊ3��[x y spacing] */
			int *nbme,			/* ��������Ԫ��Ŀ */
			int **bm_topu,		/* ��������Ԫ������Ϣ��ǰ3άΪ�����ţ���3άΪ���ڹ�ϵ */
			int nsrc[2],		/* ��Դ��Ŀ, ��Դ��Ŀ */
			double **src_info	/* ��Դ����Դ��Ϣ��ÿ����Դ����5ά��Ϣ��
									[x y spacing r1 r2]; ÿ����Դ����2����Դ */
			);

/* ---------------------------------------------------------------------------------
* ��ϵͳ�ж����3D�ܶȿ�����Ϣ�Ի������͵�������ʽ��ȡ�������񻯺�������
* ������
*		nbmn				�������񶥵���Ŀ
*	    bmn_info			�������񶥵���Ϣ��γ��Ϊ4��[x y z spacing]
*		nbme				��������Ԫ��Ŀ
*		bm_topu				��������Ԫ������Ϣ��ǰ4άΪ�����ţ���4άΪ���ڹ�ϵ
*		nsrc[3]				��Դ��Ŀ, ��Դ��Ŀ, ��Դ��Ŀ
*		src_info			��Դ����Դ����Դ��Ϣ��ÿ����Դ����6ά��Ϣ��
*							[x y z spacing r1 r2]; ÿ����Դ����2����Դ;
*							ÿ����Դ����3����Դ;
* ����ֵ: 
*      �������
*
* ������
*
* -------------------------------------------------------------------------------*/
int construct_primal_spacing_3D(
			double *gvalue,		/* ȫ���ܶȳߴ� */
			int *nbmn,			/* �������񶥵���Ŀ */
			double **bmn_info,	/* �������񶥵���Ϣ��γ��Ϊ4��[x y z spacing] */
			int *nbme,			/* ��������Ԫ��Ŀ */
			int **bm_topu,		/* ��������Ԫ������Ϣ��ǰ4άΪ�����ţ���4άΪ���ڹ�ϵ */
			int nsrc[3],		/* ��Դ��Ŀ, ��Դ��Ŀ,  ��Դ��Ŀ*/
			double **src_info	/* ��Դ����Դ����Դ��Ϣ��ÿ����Դ����6ά��Ϣ��
									[x y z spacing r1 r2]; ÿ����Դ����2����Դ;
									ÿ����Դ����3����Դ; */
			);

//typedef std::vector<CurveDiscretization>		CurveDiscretizationArray;
//typedef std::vector<CurveDiscretization>::iterator		CurveDiscretizationArrayIterator;
//extern CurveDiscretizationArray					curveDiscretizationArray;
//extern CurveDiscretizationArrayIterator			curve_discretization_iterator;
//bool readFr2(const char* fname, int *nbp, double **lcoord, int **ibnd);
//bool readBa2(const char* fname, int *nbmn, double **nd_info, int *nbme, int **bm_topu, int nsrc[2], double **src_info);
//bool writePl2(const char* fname, int nbp, int ngp, double *g_coord, int nel, int *g_topu, int *ibnd, int *prt_arr);

/* ---------------------------------------------------------------------------------
* �������ǰ��������Ԫ��������������������Ԫ��������Ԫ��������Ԫ
*
* ����ֵ: 
*      �������
*
* ������
*
* -------------------------------------------------------------------------------*/
int calculate_grid_number();


//for debug
void write_out_info(int& nbp, double*& lcoord, int& ibnd_size, int*& ibnd);

#endif /* __eemas_geometrydiscretization_h__ */
