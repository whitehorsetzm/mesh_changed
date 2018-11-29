/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * ��ѧ��Ӧ��ģ��ĸ��ܻ���
 * Enabling Environment for Multi-displinary Application Simulations
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2007��10��15��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, Oct. 15, 2007
 * 
 * ��ϵ��ʽ (For further information, please conctact)
 *   �绰 (Tel)��+86-571-87953165
 *   ���� (Fax)��+86-571-87953167
 *   ���� (Mail)��chenjj@zju.edu.cn
 *
 * �ļ����� (File Name)��spacing.h
 * ��ʼ�汾 (Initial Version): V1.0
 * ���ܽ��� (Function Introduction��
 *     ������һ���ܶȿ��ƻ���
 *     Define a set of element spacing controlling scheme.
 * 
 * 
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):			����ʤ
 * �޸����� (Revision Date):	2007��12��22��
 * ��ǰ�汾 (Current Version):	
 * �޸Ľ��� (Revision Introduction):
 *		��GlobalSpacing���������˻�ñ�����������Դȫ����С�ߴ�Ĺ��ܣ������ʵ��
 *		��Ҫ��ͨ���ڸ��������������get_min_spacing()����
 *
 * ------------------------------------------------------------------------------
 * 
 *
 * -----------------------------�޸ļ�¼ (Revision Record)------------------------
 * �޸��� (Revisor):
 * �޸����� (Revision Date):
 * ��ǰ�汾 (Current Version):
 * �޸Ľ��� (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

# ifndef __eemas_spacing_h__
# define __eemas_spacing_h__

#include "entityiterator.h"
#include   "ferguson_curve.h"

/* ���Ƕ���ȫ���ܶȿ�����Ϣ����ߴ���Ŀ���ֲ��;
   ǰ12λ���������ڴ洢��𣬺�24λ������洢����
   We redefine the global spacing data structures to
   improve the resuing performance of these codes 
 */
#define GLO_VAL 0x1
#define BKG_MSH 0x2
#define GRID_SRC 0x4

/* Category GLO_VAL */
#define GLO_VAL_NORMAL ((1<<12) & GLO_VAL)

/* Category BKG_MESH */
#define BKG_MSH_NORMAL ((1<<12) & BKG_MSH)

/* Category GRID_SRC */
#define PNT_SRC ((1<<3) & GRID_SRC)
#define LIN_SRC ((1<<3) & GRID_SRC)
#define TRI_SRC ((1<<3) & GRID_SRC)

#define MAX_SPACING_VALUE 1e12

#define MAX_GLOBAL_SPACING 1.e12

#define minmin(a,b) (((a) < (b)) ? (a) : (b))

enum SpatialType
{
	TWODIMENSION = 2,
		THREEDIMENSION
};


class GlobalSpacingObject
{
public:
	GlobalSpacingObject() { obj_enabled = true; }
	~GlobalSpacingObject() {}

	virtual int get_type() = 0;
	virtual double get_pnt_spacing(double x, double y, double z = 0.0) = 0;
 	virtual double get_pnt_BCS( int &elem, double x, double y, double z = 0.0, int hints = 0) = 0;
	virtual bool enable(bool is_on) {obj_enabled = is_on; return true;}
	bool is_enabled() {return obj_enabled;}

	virtual double get_min_spacing() = 0;//�����õ����������ļ��е�һ����С�ߴ���Ϣ����psue�з�����ɢ����ʱҪ��
protected:
	bool obj_enabled;
};

typedef EntityHandle GSOHandle;
typedef ArrayIterator GSOHandleIterator;
typedef std::vector<GSOHandle> GSOHandleArray;

GlobalSpacingObject *get_spacing_object(GSOHandle gh);

class GlobalSpacing
{
public:
	GlobalSpacing():spa_obj_iterator(&spacing_objects) {}
	~GlobalSpacing() 
	{
		GSOHandleIterator& iterator = get_iterator();
		iterator.restart();
		while (!iterator.is_end())
		{
			GlobalSpacingObject *obj = 
				get_spacing_object(*iterator);
			if (obj)
			{
				delete obj;
			}
			++iterator;
		}
	}
	
	double minimal_spacing(int &elem, double x, double y, double z = 0.0,int hints = 0);

	GSOHandleIterator& get_iterator();
	unsigned int get_object_count();

	int add_GSO_object(GlobalSpacingObject *obj);
	int del_GSO_object(GlobalSpacingObject *obj);
	
	//added by ����ʤ -- 2007_11_30
	double get_min_spacing()//�����õ����������ļ��е�һ����С�ߴ���Ϣ��ʹ��psue��������ɢ����ʱҪ��
	{
		double min_s = MAX_GLOBAL_SPACING, s;

		GSOHandleIterator& iterator = get_iterator();
		iterator.restart();
		while (!iterator.is_end())
		{
			GlobalSpacingObject *obj = 
				get_spacing_object(*iterator);
			if (obj)
			{
				s = obj->get_min_spacing();
				if (s < min_s)
					min_s = s;
			}
			++iterator;
		}
		return min_s;
	}

protected:
	GSOHandleArray spacing_objects;
	GSOHandleIterator spa_obj_iterator;

};//class GlobalSpacing

typedef struct DirSpacing
{
	double axis_x, axis_y, axis_z;
	double spacing;
} DirSpacing;

class BKGNode
{
public:
	BKGNode() { pos_x = pos_y = pos_z = 0.0; }
	BKGNode(double cx, double cy, double cz) {pos_x = cx; pos_y = cy; pos_z = cz;}
	~BKGNode() {}

	double get_x() { return pos_x; }
	double get_y() { return pos_y; }
	double get_z() { return pos_z; }
	virtual double get_spacing(double x = 1.0, double y = 0.0, double z = 0.0) = 0;
	virtual void set_spacing(double sp) = 0;
	virtual double get_min_spacing() = 0;
	virtual BKGNode *copy() = 0;

protected:
	double pos_x, pos_y, pos_z;
};

class IsoBKGNode: public BKGNode
{
public:
	IsoBKGNode() { spacing = 0.0; }
	IsoBKGNode(double cx, double cy, double cz, double sp)
		:BKGNode(cx, cy, cz)
	{
		spacing = sp;
	}
	~IsoBKGNode() {}

	virtual double get_spacing(double x = 1.0, double y = 0.0, double z = 0.0)
	{
		return spacing;
	}
	virtual void set_spacing(double sp)
	{
		spacing = sp;
	}

	virtual double get_min_spacing()
	{
		return spacing;
	}
	virtual BKGNode *copy()
	{
		return new IsoBKGNode(pos_x, pos_y, pos_z, spacing);
	}

protected:
	double spacing;
};

#if 0
class AniBKGNode: public BKGNode
{
public:
	AniBKGNode() {dir_spacing_array = NULL; array_size = 0;}
	AniBKGNode(double cx, double cy, double cz,
		DirSpacing *sp_array, unsigned int size)
		:BKGNode(cx, cy, cz)
	{
		dir_spacing_array = sp_array;
		array_size = size;
	}
	~AniBKGNode() { if (dir_spacing_array) delete[] dir_spacing_array; }
	
	virtual double get_spacing(double x = 1.0, double y = 0.0, double z = 0.0);
	virtual void set_spacing(double sp);
	virtual double get_min_spacing()
	{
		double min_s = MAX_GLOBAL_SPACING,s;

		for(int i = 0; i<array_size; i++)
		{
			s = dir_spacing_array[i].spacing;
			if(s < min_s)
				min_s = s;
		}
		return min_s;
	}

	virtual BKGNode *copy()
	{
		unsigned int size = array_size;
		DirSpacing *sp_array = NULL;
		int i;

		if (size > 0)
		{
			sp_array = new DirSpacing[size];
			for (i = 0; i < size; i++)
			{
				sp_array[i].axis_x = dir_spacing_array[i].axis_x;
				sp_array[i].axis_y = dir_spacing_array[i].axis_y;
				sp_array[i].axis_z = dir_spacing_array[i].axis_z;

				sp_array[i].spacing = dir_spacing_array[i].spacing;
			}
		}

		return new AniBKGNode(pos_x, pos_y, pos_z, sp_array, size);
	}
protected:
	DirSpacing *dir_spacing_array;
	unsigned int array_size;
};
#endif

enum {LIE_IN, LIE_ON, LIE_OUT};

class BKGMesh: public GlobalSpacingObject
{
public:
	BKGMesh() { nodes_per_elem = 0; }
	~BKGMesh() 
	{
		int i;
		for (i = 0; i < node_array.size(); i++)
		{
			if (node_array[i])
				delete node_array[i];
		}
		if (visited)
		{
			delete []visited;
		}
		if (activeDirect)
		{
			delete[]activeDirect;
		}
	}
	
	bool set_elem_indices(int *indices, int size, int nconn)
	{
		int i;
		if (indices && size > 0 && nconn > 0 && size%nconn == 0)
		{
			elem_indices.resize(size);
			//add at 2012/12/7  visited and activeDirect
			visited = new char[size/nconn];
			activeDirect = new char[size/nconn][3];
			
			for (i = 0; i < size; i++)
			{
				elem_indices[i] = indices[i];
			}
			
			nodes_per_elem = nconn;

			return true;
		}

		return false;
	}

	bool set_elem_neighs(int *neighs, int size)
	{
		int i;
		if (neighs && size > 0)
		{
			elem_neighs.resize(size);
			for (i = 0; i < size; i++)
			{
				elem_neighs[i] = neighs[i];
			}
			
			return true;
		}
		
		return false;
	}
	bool set_elem_valid(int *valid, int size)
	{
		int i;
		if (valid && size > 0)
		{
			elem_valid.resize(size);
			for (i = 0; i < size; i++)
			{
				elem_valid[i] = valid[i];
			}
			
			return true;
		}
		return false;

	}

	bool set_node_array(BKGNode **array, int size)
	{
		int i;
		if (array && size > 0)
		{
			for (i = 0; i < node_array.size(); i++)
			{
				if (node_array[i])
					delete node_array[i];
			}
			node_array.resize(size);
			for (i = 0; i < size; i++)
			{
				node_array[i] = array[i];
			}

			return true;
		}

		return false;
	}
	bool cal_min_spacing()
	{
	   mimi_spacing = MAX_GLOBAL_SPACING;
	   double s;
	 	for(int i=0; i < node_array.size(); i++)
		{
			s = node_array[i]->get_min_spacing();
			if(s < mimi_spacing)
	 			mimi_spacing = s;	
		}
		return true;

	}

	virtual int get_type() {return BKG_MSH_NORMAL;}
	virtual double get_pnt_spacing(double x, double y, double z = 0.0);
	virtual double get_pnt_BCS( int &elem, double x, double y, double z = 0.0, int hints= 0);
	virtual double search(int hints,int* elem,double x, double y, double z = 0.0 );
	virtual double get_min_spacing()
	{
		double min_s = MAX_GLOBAL_SPACING,s;
		
// 		for(int i=0; i < node_array.size(); i++)
// 		{
// 			s = node_array[i]->get_min_spacing();
// 			if(s < min_s)
// 				min_s = s;
// 		}
		if (min_s < mimi_spacing)
			 min_s = mimi_spacing;
		return min_s;
	}

public:
	std::vector<int> elem_indices;
	std::vector<int> elem_neighs;
	std::vector<int> elem_valid;
	std::vector<BKGNode*> node_array;
	int nodes_per_elem;
	char *visited;
	char (*activeDirect)[3];
	double mimi_spacing;

};
class SurfBndBKGMesh: public BKGMesh
{
public:
	virtual ~SurfBndBKGMesh();
    //annotation by zhvliu
//	virtual double get_pnt_spacing(double x, double y, double z = 0.0);
	void setSurface(EEMAS::FergusonSurface *fs) {
		this->fs = fs;
	}
	void set_cen_rad_arrays(int nel, double *cen, double *rad);
	
protected:
	EEMAS::FergusonSurface *fs;
	std::vector<BKGNode *> cen_array;
	std::vector<double> rad_array;
};


class GlobalValue: public GlobalSpacingObject
{
public:
	GlobalValue(double sp = MAX_SPACING_VALUE) {spacing = sp;}
	~GlobalValue() {} 

	virtual int get_type() {return GLO_VAL_NORMAL;} 
	virtual double get_pnt_spacing(double x, double y, double z = 0.0) 
	{ return spacing;}
	//get_pnt_BCS need to be changed
	virtual double get_pnt_BCS( int &elem, double x, double y, double z = 0.0, int hints= 0)
	{
		return spacing;
	}

	virtual double get_min_spacing()
	{
		return spacing;
	}

protected:
	double spacing;
};

class PointSource: public GlobalSpacingObject
{
public:
	PointSource() { bkg_node = NULL; inner_radius = outer_radius = 0.0; }
	PointSource(const PointSource& src)
	{
		inner_radius = src.inner_radius;
		outer_radius = src.outer_radius;
		if (src.bkg_node)
		{
			bkg_node = src.bkg_node->copy();
		}
		else
		{
			bkg_node = NULL;
		}
	}
	~PointSource() { if (bkg_node) delete bkg_node; }

	BKGNode *get_BKG_node() {return bkg_node;}
	double get_inner_radius() {return inner_radius;}
	double get_outer_radius() {return outer_radius;}

	bool set_radii(double in, double ou) 
	{
		if (ou > in && in >= 0.0 && ou > 0.0)
		{
			inner_radius = in;
			outer_radius = ou;
			return true;
		}
		return false;
	}
	bool attach_BKG_node(BKGNode *node)
	{
		if (node)
		{
			if (bkg_node)
				delete bkg_node;
			bkg_node = node;

			return true;
		}
		return false;
	}

	virtual int get_type() { return PNT_SRC; }
	virtual double get_pnt_spacing(double x, double y, double z = 0.0);
	//need to be changed
    virtual double get_pnt_BCS( int &elem, double x, double y, double z = 0.0, int hints= 0)
	{
		return 10000000;
	}
	virtual double get_min_spacing()
	{
		double min_s = MAX_GLOBAL_SPACING;
		if(bkg_node)
			min_s = bkg_node->get_min_spacing();
		return min_s;
	}

protected:
	BKGNode *bkg_node;
	double inner_radius, outer_radius;
};

class LineSource: public GlobalSpacingObject
{
public:
	LineSource() 
	{
		pnt_src1 = pnt_src2 = NULL;
	}
	~LineSource() 
	{
		if (pnt_src1) delete pnt_src1;
		if (pnt_src2) delete pnt_src2;
	}

	PointSource *get_pnt_src1() {return pnt_src1;}
	PointSource *get_pnt_src2() {return pnt_src2;}
	
	bool attach_pnt_src1(PointSource *src)
	{
		if (src)
		{
			if (pnt_src1) delete pnt_src1;
			pnt_src1 = src;

			return true;
		}
		return false;
	}
	bool attach_pnt_src2(PointSource *src)
	{
		if (src)
		{
			if (pnt_src2) delete pnt_src2;
			pnt_src2 = src;

			return true;
		}
		return false;
	}

	virtual int get_type() { return LIN_SRC; }
	virtual double get_pnt_spacing(double x, double y, double z = 0.0);
   	virtual double get_pnt_BCS( int &elem, double x, double y, double z = 0.0, int hints= 0)
	{
		return 10000000000000000;
	}
	virtual double get_min_spacing()
	{
		double s1 = MAX_GLOBAL_SPACING,s2 = MAX_GLOBAL_SPACING;
		if(pnt_src1)
			s1 = pnt_src1->get_min_spacing();
		if(pnt_src2)
			s2 = pnt_src2->get_min_spacing();
		return s1 < s2? s1:s2;
	}


protected:
	PointSource *pnt_src1, *pnt_src2;
};

class TriangleSource: public GlobalSpacingObject 
{
public:
	TriangleSource() { pnt_src1 = pnt_src2 = pnt_src3 = NULL; }
	~TriangleSource() 
	{
		if (pnt_src1) delete pnt_src1;
		if (pnt_src2) delete pnt_src2;
		if (pnt_src3) delete pnt_src3;
	}

	PointSource *get_pnt_src1() {return pnt_src1;}
	PointSource *get_pnt_src2() {return pnt_src2;}
	PointSource *get_pnt_src3() {return pnt_src3;}

	bool attach_pnt_src1(PointSource *src)
	{
		if (src)
		{
			if (pnt_src1) delete pnt_src1;
			pnt_src1 = src;

			return true;
		}
		return false;
	}
	bool attach_pnt_src2(PointSource *src)
	{
		if (src)
		{
			if (pnt_src2) delete pnt_src2;
			pnt_src2 = src;

			return true;
		}
		return false;
	}
	bool attach_pnt_src3(PointSource *src)
	{
		if (src)
		{
			if (pnt_src3) delete pnt_src3;
			pnt_src3 = src;

			return true;
		}
		return false;
	}

	virtual int get_type() { return TRI_SRC; }
	virtual double get_pnt_spacing(double x, double y, double z = 0.0);
	//need to be changed
	virtual double get_pnt_BCS( int &elem, double x, double y, double z = 0.0, int hints= 0)
	{
		return 10000000000000000;
	}
	virtual double get_min_spacing()
	{
		double s1 = MAX_GLOBAL_SPACING,s2 = MAX_GLOBAL_SPACING,s3 = MAX_GLOBAL_SPACING;
		if(pnt_src1)
			s1 = pnt_src1->get_min_spacing();
		if(pnt_src2)
			s2 = pnt_src2->get_min_spacing();
		if(pnt_src3)
			s3 = pnt_src3->get_min_spacing();
		return minmin(minmin(s1, s2), s3);
	}

protected:
	PointSource *pnt_src1, *pnt_src2, *pnt_src3;
};

bool calc_point_coord_2d(double x,  double y,  double z,
						 double x1, double y1, double z1,
						 double *w1, double epsilon = 1.0e-18);

bool calc_point_coord_3d(double x,  double y,  double z,
						 double x1, double y1, double z1,
						 double *w1, double epsilon = 1.0e-18);
/* 
 * �����(x, y)�ĳ�������
 * calculate length coord of a point, (x, y)
 */
bool calc_length_coord_2d(double x,  double y, 
					    double x1, double y1,
					    double x2, double y2, 
					    double *w1, double *w2,  
						double epsilon = 1.0e-18);

/* 
 * p(x, y, z); p1(x1, y1, z1); p2(x2, y2, z2)
 * �����p���߶�p1p2��ͶӰ��ĳ�������
 * calculate length coord of a point, the projection of p, in the plane, p1, p2, p3
 */
bool calc_length_coord_3d(double x,  double y,  double z, 
					    double x1, double y1, double z1,
					    double x2, double y2, double z2,
					    double *w1, double *w2, 
						double epsilon = 1.0e-18);

/* 
 * �����(x, y)���������
 * calculate area coord of a point, (x, y)
 */
bool calc_area_coord_2d(double x,  double y, 
					    double x1, double y1,
					    double x2, double y2, 
					    double x3, double y3,
					    double *w1, double *w2, double *w3, 
						double epsilon = 1.0e-18);

/* 
 * p(x, y, z); p1(x1, y1, z1); p2(x2, y2, z2); p3(x3, y3, z3)
 * �����p��ƽ��p1p2p3��ͶӰ����������
 * calculate area coord of a point, the projection of p, in the plane, p1, p2, p3
 */
bool calc_area_coord_3d(double x,  double y,  double z, 
					    double x1, double y1, double z1,
					    double x2, double y2, double z2,
					    double x3, double y3, double z3,
					    double *w1, double *w2, double *w3, 
						double epsilon = 1.0e-18);

bool calc_volume_coord(double x,  double y,  double z,
				   	   double x1, double y1, double z1,
					   double x2, double y2, double z2,
					   double x3, double y3, double z3,
					   double x4, double y4, double z4,
					   double *w1, double *w2, double *w3, double *w4,
					   double epsilon = 1.0e-18);

double calc_volume(double x1, double y1, double z1,
				   double x2, double y2, double z2,
				   double x3, double y3, double z3,
				   double x4, double y4, double z4);

# endif /* __eemas_spacing_h__ */
