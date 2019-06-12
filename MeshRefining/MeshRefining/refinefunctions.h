#ifndef REFINEFUNCTIONS_H
#define REFINEFUNCTIONS_H
#include "dataclass.h"
#include "globaldefine.h"
#include <map>
#include <string>
#include "mpi.h"
#include <vector>
#include "mshgen3d_def.h"
#include "reflection.h"
using namespace std;
class Refletion;
class newNode
{
public:
    newNode& operator =(const newNode&node)
    {
        coord=node.coord;
        localID=node.localID;

        return *this;
    }
public:
    Vector coord;
    int localID;

};



string IntToString(int m);

int StringToInt(string m);



int partition(HYBRID_MESH&tetrasfile, HYBRID_MESH *tetrasPart, int nparts,int stride=1,int offset=0,int type=0);

int partition_test(HYBRID_MESH&tetrasfile, HYBRID_MESH *tetrasPart, int nparts,int stride=1,int offset=0,int type=0);


int sortPointsID(HYBRID_MESH&mesh);


int meshRefining(HYBRID_MESH &tetrasfile,HYBRID_MESH &newTetrasfile,int partMarker,Refletion &_table);

int meshRefining(HYBRID_MESH &tetrasfile,HYBRID_MESH &newTetrasfile,int partMarker);

int sixNodesPattern(HYBRID_MESH&oldmesh, HYBRID_MESH &newmesh, map<string, newNode> &edgeHash,Refletion &ref);

int sixNodesPattern(HYBRID_MESH&oldmesh, HYBRID_MESH &newmesh, map<string, newNode> &edgeHash);

int splitFacet(HYBRID_MESH&oldmesh, HYBRID_MESH &newmesh);

int constructFacets(HYBRID_MESH&mesh, HYBRID_MESH &globalMesh , map<string, int64_t> &tri_globalID);//construct inter facets

int constructFacets_test(HYBRID_MESH&mesh, HYBRID_MESH& globalMesh,map<string,int64_t>&tri_globalID);

int constructOneTriangle(int *vertices,TRI& triangle);//construct one triangle

int updateTriIndex(HYBRID_MESH &mesh, map<string, int64_t> &tri_globalID);

int managerProc(int rank, int nparts, HYBRID_MESH *tetrasPart, HYBRID_MESH& construcTetras,int &refineSize,vector<string>&bcstring,MPI_Comm comm=MPI_COMM_WORLD);

int managerProc(int rank, int nparts, HYBRID_MESH *tetrasPart, HYBRID_MESH& construcTetras,Refletion ref ,Refletion construcref, int &refineSize,vector<string>&bcstring,MPI_Comm comm=MPI_COMM_WORLD);

int workerProc(int rank, HYBRID_MESH&construcTetras,Refletion construcref, int &refineSize,vector<string>&bcstring,MPI_Comm comm=MPI_COMM_WORLD);

int workerProc(int rank, HYBRID_MESH&construcTetras, int &refineSize,vector<string>&bcstring,MPI_Comm comm=MPI_COMM_WORLD);

int meshImproved(HYBRID_MESH&mesh);

#endif // REFINEFUNCTIONS_H
