#ifndef DATACLASS_H
#define DATACLASS_H

#include <cstdlib>
#include <set>
#include "vector.h"
#include <set>

using std::set;

#if __cplusplus < 201103L
#define nullptr NULL
#endif

typedef int64_t idx_t;

/* if define macro _INTERNAL_FACE_FIRST_ and define macro _SURFACE_FACE_FIRST_
 * mark global index of faces in order:
 *      internal face,
 *      physical boundayr face,
 *      interface face.
 * if define macro _INTERNAL_FACE_FIRST_ but not define macro _SURFACE_FACE_FIRST_
 * mark global index of faces in order:
 *      internal face,
 *      interface face,
 *      physical boundayr face.
 * if not define the two macroes
 * mark global index of faces in order:
 *      boundary face (physical boundayr face and interface face)
 *      internal face
 */
#define _INTERNAL_FACE_FIRST_

using std::set;

template <int NoVertices=4,int NoFacets=4>
class ELEMENT
{
public:
    ELEMENT();
    ELEMENT &operator=(const ELEMENT &element);


    ~ELEMENT();
public:
    idx_t index;
    int localID;
    int partMarker;
    int vertices[NoVertices];
    int neighbors[NoFacets];
    double neighborsmark[NoFacets];
    static const int NumVertices = NoVertices;
    static const int NumNeighbors = NoFacets;
};

template <int NoVertices,int NoEdges>
class FACET
{
public:
    FACET()
    {
        index=-1;
        iCell=-1;
        iSurf=-1;
        iOppoProc=-1;
        partMarker=-1;
        iSurface=-1;
        for(int i=0;i<NoEdges;i++)
        {
            addedNodes[i]=-1;
        }
        for(int i=0;i<NoVertices;i++)
        {
            addedNodes[i]=-1;
        }
    }
    ~FACET() { }
    FACET& operator =(const FACET &facet);

    idx_t index;
    idx_t iSurface;
    // int iCell; // the id of the element which contains this facet
    int iCell; // the local id of the element which contains this facet
    // bool bInterface; // true for interface, false for boundary facet.      //not used
    int iSurf; // the geometric surface index
    int iOppoProc; // the processor index of its neighbour
    // int iOppoIndex; // the local index of the facet on its neighboring processor  //not used
    int addedNodes[NoEdges];
    int vertices[NoVertices];
    int partMarker; //the part marker of  the triangle
};


template <int NoVertices,int NoEdges>
FACET<NoVertices,NoEdges>& FACET<NoVertices,NoEdges>::operator =(const FACET<NoVertices,NoEdges>&facet)
{
    index=facet.index;
    iCell=facet.iCell;
    iSurf=facet.iSurf;
    iSurface=facet.iSurface;
    iOppoProc=facet.iOppoProc;
    partMarker=facet.partMarker;
    for(int i=0;i<NoVertices;i++)
    {
        vertices[i]=facet.vertices[i];
    }
    for(int i=0;i<NoEdges;i++)
    {
        addedNodes[i]=facet.addedNodes[i];
    }
    return *this;
}

template <int NoVertices,int NoFacets>
ELEMENT<NoVertices,NoFacets>::ELEMENT()
{
    //if()
    index=-1;
    localID=-1;
    partMarker=-1;
    // NumVertices=NoVertices;
    // NumNeighbors=NoFacets;
    // vertices=new int[NumVertices];
    // neighbors=new int[NumNeighbors];

    for(int i=0;i<NumNeighbors;i++)
    {
        neighbors[i]=-1;
    }
    for(int i=0;i<NumNeighbors;i++)
    {
        neighborsmark[i]=-1;
    }
}

template <int NoVertices,int NoFacets>
ELEMENT<NoVertices,NoFacets> &ELEMENT<NoVertices,NoFacets>::operator =(const ELEMENT<NoVertices,NoFacets>&element)
{
    //if()
    index=element.index;
    localID=element.localID;
    partMarker=element.partMarker;
    // NumVertices=element.NumVertices;
    // NumNeighbors=element.NumNeighbors;
    // vertices=new int[NumVertices];
    for(int i=0;i<NumVertices;i++)
    {
        vertices[i]=element.vertices[i];
    }
    // neighbors=new int[NumNeighbors];
    for(int i=0;i<NumVertices;i++)
    {
        neighbors[i]=element.neighbors[i];
    }
    return *this;
}

template <int NoVertices,int NoFacets>
ELEMENT<NoVertices,NoFacets>::~ELEMENT()
{
    /*
    if(vertices!=nullptr)
    {
        delete []vertices;

        vertices=nullptr;
    }
    if(neighbors!=nullptr)
    {
        delete []neighbors;

        neighbors=nullptr;
    }
    */
}


typedef ELEMENT<4,4> TETRAS;
typedef ELEMENT<8,6> HEX;
typedef ELEMENT<6,5> PRISM;
typedef ELEMENT<5,5> PYRAMID;
typedef FACET<3, 3> TRI;
typedef FACET<4, 4> QUAD;
//typedef FACET<3,3> TRI;

class Node
{
public:
    Node& operator=(const Node &node)
    {
        index=node.index;
        partMarker=node.partMarker;
        coord=node.coord;
        localID=node.localID;
        procs=node.procs;

        NumProcs=node.NumProcs;

     //   procs = node.procs;
        // iOppoProc = node.iOppoProc;
        // iOppoIndex = node.iOppoIndex;
        return *this;
    }
    Node();
public:
    idx_t index;//global index
    int partMarker;
    Vector coord;
    int localID;
    int flag;

    set<int> procs;
    int NumProcs;

    // int iOppoProc; // the processor index of its neighbour
    // int iOppoIndex; // the local index of the node on its neighboring processor
};

//template <int NoVertices,int NoFacets>
class HYBRID_MESH
{
public:
    HYBRID_MESH();
    ~HYBRID_MESH();
    int clear();
    int numOfCells() const { return NumTetras + NumHexes + NumPrsm + NumPyra; }
    idx_t numOfFacets() const { return NumTris + NumQuads; }
public:
    int NumNodes;
    int NumTetras;
    int NumHexes;
    int NumPyra;
    int NumPrsm;

    idx_t NumUniqueSurfFacets;
    idx_t NumUniqueInterfFacets;

    int NumTris;
    int NumQuads;

    Node *nodes;

    TETRAS *pTetras;
    PYRAMID *pPyras;
    PRISM *pPrisms;
    HEX *pHexes;

    TRI *pTris; // triangular boundary facets
    QUAD *pQuads; // quadangular boundary facets
};


#endif // DATACLASS_H
