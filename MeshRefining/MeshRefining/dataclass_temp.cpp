#include "globaldefine.h"
#include "dataclass.h"


Node::Node()
{
    index=-1;
    localID=-1;
    coord.x=0;
    coord.y=0;
    coord.z=0;
    flag=0;
    procs.clear();
    NumProcs=procs.size();
    partMarker=-1;
}

HYBRID_MESH::HYBRID_MESH()
{
    NumNodes=0;
    NumTetras=0;
    NumPrsm=0;
    NumPyra=0;
    NumHexes=0;
	NumUniqueInterfFacets=0;
	NumUniqueSurfFacets=0;
	NumTris=0;
    NumQuads=0;
    nodes=nullptr;
    pTetras=nullptr;
    pHexes=nullptr;
    pTris=nullptr;
    pPyras = nullptr;
    pPrisms = nullptr;
    pQuads = nullptr;
}
int HYBRID_MESH::clear()
{
    if(nodes!=nullptr)
    {
        delete []nodes;
        nodes=nullptr;
    }
    if(pTetras!=nullptr)
    {
        delete []pTetras;
        pTetras=nullptr;
    }
    if(pHexes!=nullptr)
    {
        delete []pHexes;
        pHexes=nullptr;
    }
    if (pPyras != nullptr)
    {
        delete [] pPyras;
        pPyras = nullptr;
    }
    if (pPrisms != nullptr)
    {
        delete [] pPrisms;
        pPrisms = nullptr;
    }
    if(pTris!=nullptr)
    {
        delete []pTris;
        pTris=nullptr;
    }
    if (pQuads != nullptr)
    {
        delete [] pQuads;
        pQuads = nullptr;
    }
    NumNodes=0;
    NumTetras=0;
    NumHexes=0;
    NumPrsm = 0;
    NumPyra = 0;
    NumTris = 0;
    NumQuads = 0;
	NumUniqueInterfFacets=0;
	NumUniqueSurfFacets=0;
	return 1;
}
HYBRID_MESH::~HYBRID_MESH()
{
    clear();
}
