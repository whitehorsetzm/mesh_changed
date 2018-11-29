#ifndef TABLE_H
#define TABLE_H
#include<map>
#include<vector>
#include<set>
#include<queue>
#include"ferguson_curve.h"
#include"ferguson_surface.h"
#include"adaptive_io.h"
#include"dataclass.h"
class HYBRID_MESH;
class HybridMesh;
using namespace std;
using namespace EEMAS;
typedef int64_t idx_t;
extern HybridMesh mesh;
class FACET_temp
{
public:
    FACET_temp();
    idx_t iSurface;
    idx_t index;
    int vertices[3];
};




class Node_temp
{
public:
    Node_temp();
    idx_t index;//global index
    Vector coord;
    int flag;
};

class HybridMesh
{
public:
    HybridMesh();
    ~HybridMesh();
public:
    int NumNodes;
    int NumTris;
    int flag;

    Node_temp *nodes;

    FACET_temp *pTris; // triangular boundary facets
};

class temp{
public:
    int patch_id;
    int curve_id;
    int node_id;
    temp():patch_id(-1),curve_id(-1),node_id(-1){/*cout<<"sucess"<<endl;*/}
};


class Refletion
{
public:
    Refletion():has_initial(0){}
    Refletion(Refletion *a){
        this->subject_table=a->subject_table;
        this->curves=a->curves;
        this->sufaces=a->sufaces;
    }
    void read_gm(const char *filename);

    void initial(int NumNodes,double *nodecoord,int Numtris,int *vertices);

    int GetPatch_ID(int face_ID);

    void attachFace(int face_ID,int patch_ID);

    int  detachFace(int face_ID);

    void reflection(int patch_ID_1,int patch_ID_2,double &x,double &y,double &z);


//    void edgeToFace(HYBRID_MESH &mesh,int &face_id_1, int &face_id_2,string edge_name);

//    void initialEdgeTable(HYBRID_MESH &mesh);

//    Vector subject_test(int patch_ID_1,int patch_ID_2,Vector coord);

//    int GetCurve_ID(int patch_id_1,int patch_id_2,double x,double y,double z);

//    void project_to_surface(int patch_id,double x,double y,double z);

//    void project_to_curve(int patch_id_1,double x,double y,double z);

//    Vector subject_face_id(int face_ID_1,int face_ID_2,Vector coord);


    map<int,int> subject_table;
    bool has_initial;
private:
    void datainitail(HybridMesh& mesh);

    Vector subjectPatchId(int patch_ID_1,int patch_ID_2,Vector coord);

    FergusonCurve* findcurve(int patch_id_1,int patch_id_2,Vector coord);

    int findsurface(int curve_id_1,int curve_id_2,int curve_id_3,Vector vector);

    temp *ref_table;
    GBSolid gbsolid;
    HybridMesh mesh;

    FergusonSurface *sufaces;
    FergusonCurve   *curves;
//    vector<string> lines;

};

class vertex{
public:
    vertex(){
        neighbor.clear();
        flag=-1;
        surface=nullptr;
        curve=nullptr;
        index=-1;
        face.clear();
        on_curve=false;
    }
  int index;
  set<int> neighbor;
  set<int> face;  //face contain this pointer
  FergusonSurface *surface;
  FergusonCurve   *curve;
  int flag;
  bool on_curve;
};


#endif // TABLE_H
