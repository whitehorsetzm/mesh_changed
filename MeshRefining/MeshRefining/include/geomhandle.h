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
 * 文件名称 (File Name)：geomhandle.h
 * 初始版本 (Initial Version): V1.0
 * 功能介绍 (Function Introduction：
 *     定义了一套针对数据引用机制, 使得通过GeometryHandle来引用相应的
 *     几何信息变得更健壮和更轻松
 *     Define a set of data reference structures, to make it robust and easy to
 *     get geometry information with GeometryHandle objects
 * 
 *
 * -----------------------------修改记录 (Revision Record)------------------------
 * 修改者 (Revisor):
 * 修改日期 (Revision Date):
 * 当前版本 (Current Version):
 * 修改介绍 (Revision Introduction):
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __eemas_geomhandle_h__
#define __eemas_geomhandle_h__

#ifdef EEMAS_USE_OLD_STD_HEADERS
# include <vector.h>
#else
# include <vector>
#endif
# include "entityiterator.h"
#include "global.h"


class DiscretizationData;

typedef struct GeometryRef
{
	void* pointer;
	unsigned int index;
	unsigned int ref_count;
	DiscretizationData *d_data;
} GeometryRef, PointRef, CurveRef, LoopRef, FaceRef, DomainRef;

void free_geom_ref(GeometryRef* ref);

typedef EntityHandle GeometryHandle;
typedef GeometryHandle PointHandle;
typedef GeometryHandle CurveHandle;
typedef GeometryHandle LoopHandle;
typedef GeometryHandle FaceHandle;
typedef GeometryHandle DomainHandle;

#define GEOMETRY_NULL NULL
#define POINT_NULL GEOMETRY_NULL
#define CURVE_NULL GEOMETRY_NULL
#define LOOP_NULL GEOMETRY_NULL
#define FACE_NULL GEOMETRY_NULL
#define DOMAIN_NULL GEOMETRY_NULL

typedef std::list<PointHandle> PointRefContainer;
typedef std::list<CurveHandle> CurveRefContainer;
typedef std::list<LoopHandle> LoopRefContainer;
typedef std::list<FaceHandle> FaceRefContainer;
typedef std::list<DomainHandle> DomainRefContainer;

typedef ListIterator PointRefIterator;
typedef ListIterator CurveRefIterator;
typedef ListIterator LoopRefIterator;
typedef ListIterator FaceRefIterator;
typedef ListIterator DomainRefIterator;

extern PointRefContainer point_ref_container;
extern CurveRefContainer curve_ref_container;
extern LoopRefContainer loop_ref_container;
extern FaceRefContainer face_ref_container;
extern DomainRefContainer domain_ref_container;

extern PointRefIterator point_ref_iterator;
extern CurveRefIterator curve_ref_iterator;
extern LoopRefIterator loop_ref_iterator;
extern FaceRefIterator face_ref_iterator;
extern DomainRefIterator domain_ref_iterator;


bool set_point_ref_null(GBPoint *p);
bool set_curve_ref_null(GBCurve *c);
bool set_loop_ref_null(GBLoop *l);
bool set_face_ref_null(GBFace *f);


bool add_point_ref(GBPoint *p, unsigned int idx);
bool add_curve_ref(GBCurve *c, unsigned int idx);
bool add_loop_ref(GBLoop *l, unsigned int idx);
bool add_face_ref(GBFace *f, unsigned int idx);


PointHandle map_point_handle(GBPoint *p);
PointHandle map_point_handle(unsigned int idx);
CurveHandle map_curve_handle(GBCurve *c);
CurveHandle map_curve_handle(unsigned int idx);
LoopHandle map_loop_handle(GBLoop *l);
LoopHandle map_loop_handle(unsigned int idx);
FaceHandle map_face_handle(GBFace *f);
FaceHandle map_face_handle(unsigned int idx);

DomainHandle map_domain_handle(unsigned int idx);

void del_point_handle(PointHandle ph);
void del_curve_handle(CurveHandle ch);
void del_loop_handle(LoopHandle lh);
void del_face_handle(FaceHandle fh);
void del_domain_handle(DomainHandle dh);

void assign_point_handle(PointHandle ph1, PointHandle &ph2);
void assign_curve_handle(CurveHandle ch1, CurveHandle &ch2);
void assign_loop_handle(LoopHandle lh1, LoopHandle &lh2);
void assign_face_handle(FaceHandle fh1, FaceHandle &fh2);
void assign_domain_handle(DomainHandle dh1, DomainHandle &dh2);

GBPoint *get_point(PointHandle ph);
GBCurve *get_curve(CurveHandle ch);
GBLoop *get_loop(LoopHandle lh);
GBFace *get_face(FaceHandle fh);


unsigned int get_point_index(PointHandle ph);
unsigned int get_curve_index(CurveHandle ch);
unsigned int get_loop_index(LoopHandle lh);
unsigned int get_face_index(FaceHandle fh);
unsigned int get_domain_index(DomainHandle dh);

DiscretizationData *get_point_ddata(PointHandle ph);
DiscretizationData *get_curve_ddata(CurveHandle ch);
DiscretizationData *get_loop_ddata(LoopHandle lh);
DiscretizationData *get_face_ddata(FaceHandle fh);
DiscretizationData *get_domain_ddata(DomainHandle dh);

void attach_point_ddata(PointHandle ph, DiscretizationData *data);
void attach_curve_ddata(PointHandle ch, DiscretizationData *data);
void attach_loop_ddata(PointHandle lh, DiscretizationData *data);
void attach_face_ddata(PointHandle fh, DiscretizationData *data);
void attach_domain_ddata(PointHandle dh, DiscretizationData *data);

void add_ref_count(GeometryHandle gh);
void subtract_ref_count(GeometryHandle gh);
bool is_valid_handle(GeometryHandle gh);
#endif /* __eemas_geomhandle_h__ */