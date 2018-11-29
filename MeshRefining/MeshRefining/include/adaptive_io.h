#ifndef ADA_IO_H
#define ADA_IO_H
#include <stdio.h>
#include "math.h"
#include "global.h"
#include "ferguson_curve.h"
#include "ferguson_surface.h"
#include<iostream>
#include "data_io.h"
//#include "cork.h"


//int vtk_to_corkmesh(VTKInput&vtk,CorkTriMesh&mesh1,CorkTriMesh&mesh2,int body1,int body2);
//int writeVTK_corkmesh(char*filename,CorkTriMesh&meshout);


//3d
int read_gm3(const char *filename);
//bulid the relationship based on points
int buildRelationshipByPoint(DiscretSolid *discreteSolid);

int buildFacetRelationshipByEdge(DiscretSolid *discreteSolid);
//attaching the Lines into Curves
int buildBoundaryInformation(DiscretSolid*discreteSolid, GBSolid *gbsolid);
//attaching the vertex and the line information

int buildBoundaryInformationOnDiscreteSolid(DiscretSolid &discreteSolid);

int rebuildTopologySolid(DiscretSolid*discreteSolid, GBSolid *gbsolid);

int rebulidTopologyVertex(DiscretSolid*discreteSolid, GBSolid *gbsolid);
int rebuildTopologyCurve(DiscretSolid*discreteSolid, GBSolid *gbsolid);
int rebuildTopologyFace(DiscretSolid*discreteSolid, GBSolid *gbsolid);

int mergeCurve(DiscretSolid&discreteSolid,GBSolid&gbsolid);

int mergeFace(DiscretSolid&discretSolid, GBSolid&gbsolid);

int anotherPointInLines(DiscretSolid&discreteSolid,set<int>&lines,int pointID,int &linesID);

int splitCurve(DiscretSolid&discreteSolid, GBSolid&gbsolid, int NumNewCurve, multimap<double, int> &parameterss, int CurveID);

int constructNewFace(DiscretSolid&discreteSolid,GBSolid&gbsolid,int faceID);

int updateFaceLoop(DiscretSolid&discreteSolid,GBSolid&gbsolid);

int MergeGBVertex(GBSolid&gbsolid,int precision=6);

int AvoidDuplicateVertexInCurve(GBSolid&gbsolid);

int prepareWriteGM3(DiscretSolid &discretSolid, GBSolid&gbsolid);


int gbsolidToGM3DATA(GBSolid&gbsolid,GM3Data &geomData);

string IntToString(int m);

int StringToInt(string m);

int read_vtk(const char *filename, DiscretSolid *discretsolid);

int VTKInputToDiscretSolid(VTKInput&vtkinput,DiscretSolid&discretsolid);

int mergeMap_VTKpoints(VTKInput&vtkinput, map<int,int>&ID_referenceID, int precision=6);

int merge_vtk_points(VTKInput&vtkinput,int precision=6);

int read_vtk(const char*filename, VTKInput&vtkinput, int precision=6);

int read_vtk_surface_imprint(const char*filename, VTKInput&vtkinput);//read vtkfile for surface imprint

int read_gm3(const char*filename,GBSolid*gbsolid);

int checkConnectivity(GBSolid&gbsolid,DiscretSolid&discreteSolid);
int faceConnectivityProcess(GBSolid&gbsolid,DiscretSolid &discreteSolid,int faceID);

void  createFergusonCurve(GBSolid *gbsolid, GBCurve *c, EEMAS::FergusonCurve *fc);

void  createFergusonCurve_General(Vector pointIDs[],int NumPoints,EEMAS::FergusonCurve *fc);

void  createFergusonFace(GBFace c,EEMAS::FergusonSurface* fer_face);

//void save_geom_3d_gm3(char* filename);
//void find_basicgeom_index(void);

int AttributesModification(VTKInput &vtk,const map<int,int>&oldID_to_newID_curve, const map<int,int>&oldID_to_newID_face);

int writeVTK(char*filename , VTKInput &vtk);
int WriteGM3(char *filename, GBSolid&gbsolid);

int WritePLS(char*filename,VTKInput&vtk);

int WritePLS_valid(char*filename,VTKInput&vtk,const map<int,int>&oldID_to_newID_face);

//assign map at this function;
int WriteFliVlm_valid(char*filename, GM3Data &gm3data, map<int,int>&oldID_to_newID_curve, map<int,int>&oldID_to_newID_face);

int WriteGM3_valid(char*filename, GM3Data &gm3data, const map<int,int>&oldID_to_newID_curve, const map<int,int>&oldID_to_newID_face);

int WriteFliVlm(char*filename,GM3Data &gm3data);

//int WriteVLM(char*filename,GM3Data &gm3data);

int DiscretSolidToVTKInput(VTKInput&vtk,DiscretSolid&discretsolid);

typedef struct
{
	double x;
	double y;
	double z;
}Vec3f;
class Comp3d
{
public:
    Comp3d(Vector start)
    {
        origin.x=start.x;
        origin.y=start.y;
        origin.z=start.z;
    }
    bool operator() ( const Vector &a,  const Vector &b) const
    {
        //double d1= (a.x-origin.x)*(a.x-origin.x)+
        double d1=(a.x-origin.x)*(a.x-origin.x)+(a.y-origin.y)*(a.y-origin.y)+(a.z-origin.z)*(a.z-origin.z);
        double d2=(b.x-origin.x)*(b.x-origin.x)+(b.y-origin.y)*(b.y-origin.y)+(b.z-origin.z)*(b.z-origin.z);
        if(d1<d2)
            return true;
        else
            return false;
    }
public:
    Vector origin;
};

class CmpVec
{
public:

    CmpVec(double m=6)
    {
        precision=m;

    }

    bool operator()( const Vec3f& _v0, const Vec3f& _v1 ) const
    {
        /*
        if(0)//old version
        {

            if (fabs(_v0.x - _v1.x) <= eps_)
            {
                if (fabs(_v0.y - _v1.y) <= eps_)
                {
                    return (_v0.z < _v1.z - eps_);//old version

                }
                else return (_v0.y < _v1.y - eps_);

            }
            else return (_v0.x < _v1.x - eps_);
        }
        */


        long long disx=(long long)(_v0.x*pow(10.0,precision)) - (long long)(_v1.x*pow(10.0,precision));
        long long disy=(long long)(_v0.y*pow(10.0,precision)) - (long long)(_v1.y*pow(10.0,precision));
        long long disz=(long long)(_v0.z*pow(10.0,precision)) - (long long)(_v1.z*pow(10.0,precision));
        if(disx==0)
        {
            if(disy==0)
            {
                return disz<0;

            }
            else
                return disy<0;
        }
        else
            return disx<0;

    }

private:
    int precision;
};

#endif
