#ifndef GLOBAL_H
#define GLOBAL_H

#include<vector>
#include<unordered_map>
#include"vector.h"
#include<set>
#include "config.h"
#include "string.h"
//#include<iostream>
using namespace std;
//using namespace stdext;
#define ADD_BOX 0

//#define test;

extern double eps_PP;//eps of point and point

extern double eps_dotProduct;

extern double eps_global;//global epsilon

extern double eps_PE;

extern double eps_PT;

extern double eps_EE;

//added by sunlisten
/* 几何体的维度 geometric dimensions*/
enum 
{
	DIM_ZERO = 0,	/* 二维点 2D point*/
	DIM_HALF,		/* 三维点 3D point*/
	DIM_ONE,		/* 二维曲线 2D curve*/
	DIM_ONE_HALF,	/* 三维曲线 3D curve*/
	DIM_TWO,		/* 二维曲面 2D surface (plane) */
	DIM_TWO_HALF,	/* 三维曲面 3D surfaces*/
	DIM_THREE		/* 三维点 3D volume*/
};

//the four intersection situation of curves
#define SITUATION0 0 //meaning intersection(intersection exists)
#define SITUATION1 1 //meaning intersection case 1
#define SITUATION2 2 //meaning intersection case 2
#define SITUATION3 3 //meaning intersection case 3
#define SITUATION4 4 //meaning intersection case 4


//mesh policy for curve
#define CURVE_POLICY_PRO_END  1		//mesh the curve proportionally, indicate the ratio between two end segments
#define CURVE_POLICY_PRO_NER  2		//mesh the curve proportionally, indicate the ratio between neighboring segments
#define CURVE_POLICY_UNI  3			//mesh the curve uniformly
#define CURVE_POLICY_BKG  4			//mesh the curve according to background mesh
// added [1/9/2009 ly]
#define CURVE_POLICY_SAMPLE 5		//sample the curve


//mesh policy for plane
#define PLANE_POLICY_STRTRI_MAP  50		//generate structured triangles using transfinite-mapping
#define PLANE_POLICY_STRQUA_MAP  51		//generate structured quadrilaterals using transfinite-mapping
#define PLANE_POLICY_UNSTRTRI_DEL  52		//generate unstructured triangles using 2d-delaunay method
#define PLANE_POLICY_UNSTRQUA_AUTOSUB  53	//generate unstructured quadrilaterals using auto-multi-subdomain method
#define PLANE_POLICY_UNSTRQUA_MANUSUB  54	//generate unstructured quadrilaterals using manual-multi-subdomain method

//mesh policy for surface
#define SURF_POLICY_STRTRI_STR  100		//generate structured triangles using structured method
#define SURF_POLICY_STRQUA_STR	101		//generate structured quadrilaterals using structured method
#define SURF_POLICY_UNSTRTRI_MAP  102		//generate unstructured triangles using mapping based method
#define SURF_POLICY_SAMPLE 103	//sample the surface // added [1/28/2009 leon]

//mesh policy for volume
#define VOL_POLICY_UNSTRTETR_DEL  200 		//generate unstructured tetrahedras using 3d-delaunay method
#define VOL_POLICY_STRHEX_SWP  201 		//generate structured hexahedrals using sweep method

//grid element/point type
#define GRIDNONE	0	//none
#define GRIDPNT		1	//point element
#define GRIDLIN		2	//line element
#define GRIDTRI		3	//triangle element
#define GRIDTET		4	//tetrahedral element
#define GRIDQUA		5	//quadrilateral element
#define GRIDHEX		6	//hexahedral element
#define CONSTRAINLINE 10

#define BOXPNT      0
#define NV			0
#define BOXSIZE     500

#ifndef maxmax
#define maxmax(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#define CHECK_FREE(a) if(a){ delete []a; a = NULL;}

#ifndef minmin
#define minmin(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#if 0
typedef double coordT;
typedef double vectT;
#else
typedef double coordT;
typedef double vectT;
#endif
typedef int ptf;

typedef struct
{
  coordT x, y, z;
} GBCoords;

typedef struct
{
	vectT x, y, z;
} GBVect;

//input data structure











/* Basic Geometry */

class GBVertex
{
public:
    struct GBPoint		*pp_t;
	GBCoords		coord;
    int index;

    int referenceID;


};


typedef struct GBPoint
{
	int 			index;
	int			v_i;
	struct GBPoint		*next;
} GBPoint;
class Bedge// boundary edges
{
public:
	int id;// id of Bedge
	int p1;//one point of the edge
	int p2;//another point of the edge 
	int triID;//the triangle id which the edge belongs to
	int faceID;//the face which the edge belongs to
	int curveID;//the curve which the edge belongs to 
public:
	Bedge()
	{
		id=-1;
		p1=-1;
		p2=-1;
		triID=-1;
		faceID=-1;
		curveID=-1;//represent the consistent curve is null
    }
	~Bedge(){}
	Bedge(const Bedge &b)
	{
		this->id=b.id;
		this->p1=b.p1;
		this->p2=b.p2;
		this->triID=b.triID;
		this->curveID=b.curveID;
		this->faceID=b.faceID;
	
	}
	int anotherpnt(int pnt)
	{
		if(pnt==p1)
			return p2;
		if(pnt==p2)
			return p1;

		return -1;
	
	}

};
class GBCurve
{
public:
    int         preIndex;
	int 		index;
	int			type;
	int			ncp;
	int			*v_i;
    Vector      *d1_vector;		/* 型值点切向量数组 */
	double      *arc_length;		/* 弧长数组, 前data_point_num-1个成员为各个分段曲线的弧长， 最后一个成员为总弧长 */

	int ngpt;			//曲线端点对应的网格点的个数
	int endgpt[2];		//曲线端点对应的网格点编号
    set<int> discretePointsID;
    set<int> discreteEdgeIDs;
    int situation;
    int NumVertex;//Number of vertex(discrete mesh vertex) the curve contains

    double u[2];//parameter length of the curve  initial is 0-1

    bool enable;
    int  referenceID;//if enable==false which means the curves is unvalid, this integer denotes the new merged Curve ID;
//note: enable =false and referenceID=-1 means the curve has splited several new curves
 //enable =false and referenceID!=-1 means the curve has been merged and the new ID is referenceID
    GBCurve()
    {
        discretePointsID.clear();
        discreteEdgeIDs.clear();
        situation = -1;
        NumVertex=2;
        ncp=-1;
        u[0]=-1;
        u[1]=ncp;
        preIndex=-1;
        enable=true;
        referenceID=-1;

    }
    ~GBCurve()
    {

        if(v_i!=nullptr)
        {
            delete []v_i;
        }
        if(d1_vector!=nullptr)
        {
            delete []d1_vector;
        }
        if(arc_length!=nullptr)
        {
            delete []arc_length;
        }

    }

    GBCurve(const GBCurve &c)
    {
        index=c.index;
        type=c.type;
        ncp=c.ncp;
        v_i=new int[ncp];
        for(int i=0;i<ncp;i++)
        {
            v_i[i]=c.v_i[i];
        }
        d1_vector=new Vector[ncp];
        for(int i=0;i<ncp;i++)
        {
            d1_vector[i]=c.d1_vector[i];
        }
        arc_length=new double[ncp];
        for(int i=0;i<ncp;i++)
        {
            arc_length[i]=c.arc_length[i];
        }
        ngpt=c.ngpt;
        endgpt[0]=c.endgpt[0];
        endgpt[1]=c.endgpt[1];
        discreteEdgeIDs=c.discreteEdgeIDs;
        discretePointsID=c.discretePointsID;
        situation=c.situation;
        NumVertex=c.NumVertex;
        u[0]=c.u[0];
        u[1]=c.u[1];
        preIndex=c.preIndex;
        enable=c.enable;
        referenceID=c.referenceID;

    }
    bool operator()(const GBCurve&c)
    {
        index=c.index;
        type=c.type;
        ncp=c.ncp;
        v_i=new int[ncp];
        for(int i=0;i<ncp;i++)
        {
            v_i[i]=c.v_i[i];
        }
        for(int i=0;i<ncp;i++)
        {
            d1_vector[i]=c.d1_vector[i];
        }
        for(int i=0;i<ncp;i++)
        {
            arc_length[i]=c.arc_length[i];
        }
        ngpt=c.ngpt;
        endgpt[0]=c.endgpt[0];
        endgpt[1]=c.endgpt[1];
        discreteEdgeIDs=c.discreteEdgeIDs;
        discretePointsID=c.discretePointsID;
        situation=c.situation;
        NumVertex=c.NumVertex;
        u[0]=c.u[0];
        u[1]=c.u[1];
        preIndex=c.preIndex;
        enable=c.enable;
        referenceID=c.referenceID;

        return true;

    }
    bool operator <(const GBCurve&curveA)const
    {
        if(discretePointsID<curveA.discretePointsID)
        {
            return true;
        }

        else
            return false;


    }
    bool operator ==(const GBCurve&curveA)const
    {

        if(discretePointsID==curveA.discretePointsID)
        {
            return true;
        }
        else
            return false;
    }




};
class GBFace
{
public:
    int 			index;
    int				type;
    int             picked;
    int				nsp1;			// number of points in a curve       direction u
    int				nsp2;			// number of curves in the face      direction v
    //int				*nsp1[3];        //three points in the face
    int				**v_i;
    GBVect			**normal;//not used
    int				nsl;
    struct GBLoop	**lp_t;

    int              ancientFace;  //donating the old face

    int domainID; //donating the domain ID
    //vector<Bedge*> boundaryedges;
    //vector<Bedge*> *boundaryedges;

    //enable == false and referenceID ==-1 means the face has been splited ; enable == true and referenceID!=-1, means the face has been merged,
    //the referecnceID denotes the merged face ID

    bool enable;//denoting whether valid

    int referenceID;//if the value is not -1, it means the face has been merged, and its value denote the new faceID;

    set<int>facets;

    set<  set<int> > classifyFacets;


    GBFace()
    {
        index=-1;
        type=0;
        picked=-1;
        nsp1=0;
        nsp2=0;
        v_i=nullptr;
        normal=nullptr;
        nsl=0;
        lp_t=nullptr;
        domainID=-1;
        enable=true;
        referenceID=-1;
        facets.clear();
        ancientFace=-1;
        classifyFacets.clear();

    }
    GBFace(const GBFace &face)
:classifyFacets(face.classifyFacets),
    facets(face.facets)
    {

        index=face.index;
        type=face.type;
        picked=face.picked;
        nsp1=face.nsp1;
        nsp2=face.nsp2;
        if(nsp2>0)
        {
            v_i=new int*[nsp2];
        }
        for(int i=0;i<nsp2;i++)
        {
            v_i[i]=new int[nsp1];
            for(int j=0;j<nsp1;j++)
            {
                v_i[i][j]=face.v_i[i][j];


            }


        }
        nsl=face.nsl;
        lp_t=new GBLoop*[nsl];
        for(int i=0;i<nsl;i++)
        {

           lp_t[i]=face.lp_t[i];
        }
        domainID=face.domainID;
        ancientFace=face.ancientFace;
        enable=face.enable;
        referenceID=face.referenceID;

        normal=nullptr;



    }
    bool operator<(const GBFace&face)const
    {
        if(facets<face.facets)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    bool operator==(const GBFace&face)
    {
        if(facets==face.facets)
        {
            return true;
        }
        else
        {
            return false;
        }

    }
    bool operator()(const GBFace&face)
    {
        index=face.index;
        type=face.type;
        picked=face.picked;
        nsp1=face.nsp1;
        nsp2=face.nsp2;
        if(nsp2>0)
        {
            v_i=new int*[nsp2];
        }
        for(int i=0;i<nsp2;i++)
        {
            v_i[i]=new int[nsp1];
            for(int j=0;j<nsp1;j++)
            {
                v_i[i][j]=face.v_i[i][j];


            }


        }
        nsl=face.nsl;
        lp_t=new GBLoop*[nsl];
        for(int i=0;i<nsl;i++)
        {

           lp_t[i]=face.lp_t[i];
        }
        enable=face.enable;
        referenceID=face.referenceID;
        facets=face.facets;
        classifyFacets=face.classifyFacets;
        domainID=face.domainID;
        ancientFace=face.ancientFace;

        normal=nullptr;




    }
    ~GBFace()
    {
        if(v_i!=nullptr)
        {
            for(int i=0;i<nsp2;i++)
            {

                delete []v_i[i];
            }
            delete []v_i;

            v_i=nullptr;

        }
        if(normal!=nullptr)
        {
            for(int i=0;i<nsp2;i++)
            {

                delete []normal[i];
            }
            delete []normal;
            normal=nullptr;

        }

        if(lp_t!=nullptr)
        {
            delete []lp_t;

        }



    }

};

class GBLoop
{
public:
	int 		index;
	int			type;
	int			nlc;
    GBCurve		**cp_t;    // curve pointer
	int         nhl; //hardLine 的数量  //curve 可以截断，hardLine是不变的
	struct GBHardLine **hl_t;
    int			*refer;// refer to curve number  start from 0 For example: number 0 curve  is  GBcurve cpt[0]
    GBFace		*fp_t;//face pointer
	int 			f_i;
    GBLoop *next;
    bool enable;
    GBLoop()
    {
        index=-1;
        type=-1;
        nlc=0;
        cp_t=nullptr;
        nhl=0;
        hl_t=nullptr;
        refer=nullptr;
        fp_t=nullptr;
        f_i=-1;
        next=nullptr;
        enable=true;
    }
    GBLoop(const GBLoop &loop)
    {
        index=loop.index;
        type=loop.type;
        nlc=loop.nlc;
        cp_t=new GBCurve*[nlc];
        for(int i=0;i<nlc;i++)
        {
            cp_t[i]=loop.cp_t[i];
        }
        nhl=loop.nhl;
        hl_t=new GBHardLine*[nhl];
        for(int i=0;i<nhl;i++)
        {
            hl_t[i]=loop.hl_t[i];
        }
        refer=new int[nlc];
        for(int i=0;i<loop.nlc;i++)
        {
            refer[i]=loop.refer[i];
        }
        fp_t=loop.fp_t;
        f_i=loop.f_i;
        next=loop.next;
        enable=loop.enable;

    }
    ~GBLoop()
    {

        if(refer!=nullptr)
        {
            delete []refer;
            refer=nullptr;
        }
        if(cp_t!=nullptr)
        {
            delete []cp_t;
            cp_t=nullptr;
        }
        if(hl_t!=nullptr)
        {
            delete []hl_t;
            hl_t=nullptr;
        }
    }

};
class GBDomain
{
public:
    int index;
    vector<int>surfaceID;


    GBDomain()
    {
        index=-1;
        surfaceID.clear();
    }
};


//constraint line->hard line
typedef struct GBHardLine
{
	int                  nc; //number of curves in a hardline
	struct GBCurve		**cp_t;
	struct GBLoop       *lp_t;
}GBHardLine;


//store the parameter coord in the curve
class CoordU
{
public:
	int pntidex;//grid point
	double u;//the parameter coord in specific curve


};


typedef struct	
{
 	double			rad_1;
 	double			rad_2;
	double			den; 
} GBSVal;


typedef struct
{
	double			x;
	double			y;
	double			z;
} GBBkgVal;


typedef struct GBBkgPoint
{
	int			in_flag;
	int 			index;
	int			marked;
	GBBkgVal		val;
	GBCoords		coord;
	struct GBBkgPoint	*next;
} GBBkgPoint;


typedef struct GBBkgElement
{
	int index;
	struct GBBkgPoint	*connp_t[3];
	int neig[3];
	GBVect			normal;
	int flag;
	struct GBBkgElement	*next;
} GBBkgElement;


typedef struct GBPSor
{
	int			in_flag;
	GBSVal			sor;
	GBCoords		coord;
	struct GBPSor		*next;
} GBPSor;


typedef struct GBLSor
{
	int			p_picked[2];
	int			in_flag[2];
	GBSVal			sor[2];
	GBCoords		coord[2];
	struct GBLSor		*next;
} GBLSor;


typedef struct GBTSor
{
	int			p_picked[3];
	int			in_flag[3];
	GBSVal			sor[3];
	GBCoords		coord[3];
	struct GBTSor		*next;
} GBTSor;


/* Boundary Grid / Result Grid */
class DiscretPoint
{
public:
    int index;
    double x;
    double y;
    double z;
    set<int> linkedPoints;
    set<int> linkedFacets;
    bool vertex;
    bool curvePoint;  //if the curve Point
    set<int> geometryCurveID;

    int NumGeometry;
    set<int> geometryID;

    DiscretPoint()
    {
        index=-1;
        x=0;
        y=0;
        z=0;
        linkedPoints.clear();
        linkedFacets.clear();
        geometryCurveID.clear();
        curvePoint=false;
        vertex=false;
        NumGeometry=0;
        geometryID.clear();
    }
    bool operator=(const DiscretPoint&copy)
    {
        index=copy.index;
        x=copy.x;
        y=copy.y;
        z=copy.z;
        linkedPoints=copy.linkedPoints;
        linkedFacets=copy.linkedFacets;
        vertex=copy.vertex;
        curvePoint=copy.curvePoint;
        geometryCurveID=copy.geometryCurveID;
        NumGeometry=copy.NumGeometry;
        geometryID=copy.geometryID;


        return true;
    }
};
class DiscretVertex
{
public:
    int index;
    int pointID;
    set<int> geometryID;
    int uniformedGeomtryID;
    int NumGeometry;

    DiscretVertex()
    {

        index=-1;
        pointID=-1;
        NumGeometry=-1;
        uniformedGeomtryID=-1;
        geometryID.clear();

    }


};

class DiscretEdge
{
public:
    int index;
    int startPoint;
    int endPoint;
    int facetHead;//indicate the first facetID that contains the edge
    bool line;//whether it is a line which means a boundary edge
    set<int>geometryID;
    int NumGeometry;
    bool valid;
    set<int>interSectionPoints;
    set<int>newIDs;//if valid is false , this variable is enabled


    DiscretEdge()
    {
        index=-1;
        startPoint=-1;
        endPoint=-1;
        facetHead=-1;
        line=false;
        NumGeometry=0;
        geometryID.clear();
        valid=true;
        interSectionPoints.clear();
        newIDs.clear();

    }
    DiscretEdge(const DiscretEdge&edge)
    {

        index=edge.index;
        startPoint=edge.startPoint;
        endPoint=edge.endPoint;
        facetHead=edge.facetHead;
        line=edge.line;
        geometryID=edge.geometryID;
        NumGeometry=edge.NumGeometry;
        valid=edge.valid;
        interSectionPoints=edge.interSectionPoints;
        newIDs=edge.newIDs;
    }
    int anotherPointID(int pointID)
    {
        if(pointID==startPoint)
        {
            return endPoint;
        }
        else if(pointID=endPoint)
        {
            return startPoint;
        }
        else
        {
            return -1;
        }
    }
    bool operator=(const DiscretEdge&edge)
    {
        index=edge.index;
        startPoint=edge.startPoint;
        endPoint=edge.endPoint;
        facetHead=edge.facetHead;
        line=edge.line;
        geometryID=edge.geometryID;
        NumGeometry=edge.NumGeometry;
        valid=edge.valid;
        interSectionPoints=edge.interSectionPoints;
        newIDs=edge.newIDs;

        return true;
    }

    bool operator< (const DiscretEdge &ad) const
    {

        if(startPoint!=ad.startPoint||startPoint!=startPoint)
        {
            if(startPoint<ad.startPoint)
                return true;
            else
                return false;

        }


        else
            return false;
    }



};
class DiscretFacet
{
public:
    int index;
    int points[3];
    int edges[3];
    int neighborFacets[3];// the index corresponded the index of edges ; the value indicate the next neighborhood facet that shared the specific edge
    set<int>geometryID;
    int NumGeometry;
    bool valid;
    set<int>intersectPoints;
    vector<DiscretEdge>intersectEdges;
    DiscretFacet()
    {
        index=-1;
        points[0]=-1;
        points[1]=-1;
        points[2]=-1;
        edges[0]=-1;
        edges[1]=-1;
        edges[2]=-1;
        neighborFacets[0]=-1;
        neighborFacets[1]=-1;
        neighborFacets[2]=-1;
        NumGeometry=-1;
        geometryID.clear();
        valid=true;
        intersectPoints.clear();
        intersectEdges.clear();

    }
    bool operator=(const DiscretFacet&other)
    {
        index=other.index;
        points[0]=other.points[0];
        points[1]=other.points[1];
        points[2]=other.points[2];

        edges[0]=other.edges[0];
        edges[1]=other.edges[1];
        edges[2]=other.edges[2];

        neighborFacets[0]=other.neighborFacets[0];
        neighborFacets[1]=other.neighborFacets[1];
        neighborFacets[2]=other.neighborFacets[2];

        geometryID=other.geometryID;

        NumGeometry=other.NumGeometry;

        valid=other.valid;

        intersectPoints=other.intersectPoints;
        intersectEdges=other.intersectEdges;
        return true;
    }
    int thirdPointID(int a,int b)
    {
        bool flag[3]={false,false,false};
        for(int i=0;i<3;i++)
        {
            if(points[i]==a||points[i]==b)
            {
                flag[i]=true;
            }
        }

        for(int i=0;i<3;i++)
        {
            if(flag[i]==false)
                return points[i];
        }
       // return c;
    }

};

class DiscretBox
{
public:
    DiscretBox()
    {
        xmin=xmax=ymin=ymax=zmin=zmax=0;
    }
    bool operator=(const DiscretBox&box)
    {
        xmin=box.xmin;
        xmax=box.xmax;
        ymin=box.ymin;
        ymax=box.ymax;
        zmin=box.zmin;
        zmax=box.zmax;

        return true;

    }
public:
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
};
class DiscretSolid
{
public:
    int index;
    int NumPoints;
    int NumVertices;
    int NumLines;
    int NumFacets;
    int NumAllEdegs;
    char solidName[256];
    DiscretPoint*discretPoints;
    DiscretVertex*discretVertices;
    DiscretEdge*discreteLines;
    DiscretFacet*discreteFacets;
    DiscretBox discretBox;
    vector<DiscretEdge>allEdges;

    map<int,int>faceID_to_bodyID;

    unordered_map<string,int> *hashEdges;

    vector<DiscretFacet>newFacets;
    vector<DiscretEdge>newLines;
    vector<DiscretPoint>newPoints;


    DiscretSolid(DiscretSolid& other)
        : allEdges(other.allEdges)
    {
        discretPoints = other.discretPoints;
        other.discretPoints = nullptr;
        discretVertices = other.discretVertices;
        other.discretVertices = nullptr;
        discreteLines = other.discreteLines;
        other.discreteLines = nullptr;
        discreteFacets = other.discreteFacets;
        other.discreteFacets = nullptr;
        discretBox=other.discretBox;
        hashEdges = other.hashEdges;
        other.hashEdges = nullptr;
        index = other.index;
        NumPoints = other.NumPoints;
        NumVertices = other.NumVertices;
        NumLines = other.NumLines;
        NumFacets = other.NumFacets;
        NumAllEdegs = other.NumAllEdegs;
        strcpy(solidName, other.solidName);
        faceID_to_bodyID=other.faceID_to_bodyID;
        newFacets=other.newFacets;
        newLines=other.newLines;
        newPoints=other.newPoints;

    }

    DiscretSolid()
    {
        index=-1;
        NumPoints=-1;
        NumLines=-1;
        NumAllEdegs=-1;
        NumFacets=-1;
        NumVertices=-1;
        strcpy(solidName, "default");
        discretPoints=nullptr;
        discreteLines=nullptr;
        discreteFacets=nullptr;
        discretVertices=nullptr;
        hashEdges=nullptr;
        allEdges.clear();
        faceID_to_bodyID.clear();
        newFacets.clear();
        newLines.clear();
        newPoints.clear();

    }
    ~DiscretSolid()
    {

        if(discretPoints!=nullptr)
        {
            delete []discretPoints;
            discretPoints=nullptr;
        }
        if(discreteLines!=nullptr)
        {
            delete []discreteLines;
            discreteLines=nullptr;
        }
        if(discretVertices!=nullptr)
        {
            delete []discretVertices;
            discretVertices=nullptr;
        }
        if(discreteFacets!=nullptr)
        {
            delete []discreteFacets;
            discreteFacets=nullptr;

        }
        if(hashEdges!=nullptr)
        {
            hashEdges->clear();
            delete hashEdges;
            hashEdges=nullptr;
        }
    }
};

class GBSolid
{
public:
    int nvertices_G; //点de数量
    int ncurves_G; //曲线的数量
    int nloops_G;//环的数量
    int nfaces_G;//面数量
    int ndomain_G;//number of  bodies

    char solidName[256];
    int index; //ID

    GBVertex *vertex_G;
    GBCurve *curve_G;
    GBFace *face_G;
    GBLoop *loop_G;
    GBDomain*domian_G;


    vector<GBCurve>newCurves;
    vector<GBFace>newFaces;

    vector<GBVertex>newVertex;

    vector<GBLoop>newLoops;// the old loop is not used after the new loop

    GBSolid()
    {

        strcpy(solidName, "default");
        index=-1;
        nvertices_G=-1;
        ncurves_G=-1;
        nloops_G=-1;
        nfaces_G=-1;
        newCurves.clear();
        newFaces.clear();
        newVertex.clear();
        newLoops.clear();

        vertex_G=nullptr;
        curve_G=nullptr;
        face_G=nullptr;
        loop_G=nullptr;

    }
    int clear()
    {
        strcpy(solidName, "default");
        index=-1;
        nvertices_G=-1;
        ncurves_G=-1;
        nloops_G=-1;
        nfaces_G=-1;
        newCurves.clear();
        newFaces.clear();
        newVertex.clear();
        newLoops.clear();


        if(vertex_G!=nullptr)
        {
            delete []vertex_G;
            vertex_G=nullptr;
        }
        if(curve_G!=nullptr)
        {
            delete []curve_G;
            curve_G=nullptr;
        }
        if(face_G!=nullptr)
        {
            delete[]face_G;
            face_G=nullptr;
        }
        if(loop_G!=nullptr)
        {
            delete []loop_G;
            loop_G=nullptr;
        }
        return 1;
    }
    ~GBSolid()
    {
      clear();
    }

};



//extern vector<Bedge>  boundaryedges;//boundary edges

extern GBPoint *point_hd;
extern GBCurve *curve_G;
extern GBVertex *vertex_G;
extern GBFace *face_G;
extern GBHardLine *hardLine_G;
extern GBLoop *loop_G;
extern int npoints_G;
extern int nvertices_G; //点得数量
extern int ncurves_G; //曲线的数量
extern int nloops_G;//环的数量
extern int nfaces_G;//面数量
extern int nhardline_G;
extern int save_hardLines;
extern int nboxPoint; //外包围盒的点数


extern GBPSor *ptsor_G;
extern GBLSor *lnsor_G;
extern int nlnsor_G; //线源
extern int nptsor_G; //点源

extern GBBkgPoint *bkgpt_G;
extern GBBkgElement *bkgel_G;
extern double *coord;
extern int nbkgpt_G;//网格点的数量
extern int nbkgel_G;//网格单元的数量

extern int minimum_spacing;


extern double		tolg;	/* tolerance for geometry,decided by max size 
								of the geometry points--default:1.0e-6 */
extern double		GBCOINTOL; /* tolg*tolg  Tolerance for coincident vertices 
								--dafault:1.0e-4 */
extern int  nboud_G;


#endif
