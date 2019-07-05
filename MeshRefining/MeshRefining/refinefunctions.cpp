#include "refinefunctions.h"
#include "metis.h"
#include <map>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <set>
#include "dataio.h"
#include <map>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "algorithm"
#include <unordered_map>
#include <fstream>
using namespace std;


//string IntToString(int m)
//{
//    stringstream stream;
//    stream<<m;
//    string temp;
//    stream>>temp;
//    return temp;

//}

//int StringToInt(string m)
//{
//    stringstream stream;
//    stream<<m;
//    int temp;
//    stream>>temp;
//    return temp;

//}
/*
 *METIS_API(int) METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
                  idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
                  idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options,
                  idx_t *edgecut, idx_t *part);
 */
int sortPointsID(HYBRID_MESH&mesh)
{
    int nPoint=mesh.NumNodes;
    map<idx_t,int>global_local;
    map<int,int>local_newLocal;
    idx_t gloablID=-1;
    int localID=-1;
    int newLocalID=-1;
    map<idx_t,int>::iterator it;
    for(int i=0;i<nPoint;i++)
    {
        gloablID=mesh.nodes[i].index;
        localID=mesh.nodes[i].localID;
        global_local[gloablID]=localID;
    }
    int count=0;
    for(it=global_local.begin();it!=global_local.end();it++)
    {
        localID=it->second;
        local_newLocal[localID]=count;
        count++;
    }
    assert(nPoint==count);
    //update the points
    Node*nodes_new=new Node[nPoint];
    for(int i=0;i<nPoint;i++)
    {
        int newID=local_newLocal[i];
        nodes_new[newID]=mesh.nodes[i];
        nodes_new[newID].localID=newID;

    }
    if(mesh.nodes!=nullptr)
        delete[]mesh.nodes;

    mesh.nodes=nodes_new;


    //update the newLocalID
    int temp=-1;
    for(int i=0;i<mesh.NumTris;i++)
    {
        for(int j=0;j<3;j++)
        {
            temp=mesh.pTris[i].vertices[j];
            mesh.pTris[i].vertices[j]=local_newLocal[temp];
        }
    }
    for(int i=0;i<mesh.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            temp=mesh.pTetras[i].vertices[j];
            mesh.pTetras[i].vertices[j]=local_newLocal[temp];
        }

    }

    //
    return 1;
}
int partition_test(HYBRID_MESH&tetrasfile, HYBRID_MESH *tetrasPart, int nparts, int stride, int offset, int type)
{


    TETRAS *ptetras=tetrasfile.pTetras;

    PRISM *pprisms=tetrasfile.pPrisms;

    idx_t NumVertices=tetrasfile.NumNodes;

    idx_t NumElement=tetrasfile.numOfCells();

    idx_t NumPrsm=tetrasfile.NumPrsm;

    idx_t NumHexes=tetrasfile.NumHexes;

    idx_t NumTetras=tetrasfile.NumTetras;

    idx_t ncommon=3;

    idx_t NumParts=nparts;

    idx_t *eptr=new idx_t[NumElement+1];    //存的eind下标,so  need to index to get range

    idx_t *eind=new idx_t[NumPrsm*6+NumTetras*4+NumHexes*8];

    idx_t *vwgt=new idx_t[NumElement];
    for(int i=0;i<NumPrsm+1;i++)
    {
        eptr[i]=i*6;
    }
    for(int i=NumPrsm+1;i<NumPrsm+NumTetras+1;i++)
    {
        eptr[i]=NumPrsm*6+(i-NumPrsm)*4;
    }        //important!!!!!!!!

    for(int i=0;i<NumPrsm;i++)
    {
        for(int j=0;j<6;j++)
        {
            eind[i*6+j]=pprisms[i].vertices[j];
        }
    }
    for(int i=0;i<NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            eind[NumPrsm*6+i*4+j]=ptetras[i].vertices[j];
        }
    }
    for(int i=0;i<NumPrsm;i++)
    {
        vwgt[i]=6;
    }

    for(int i=NumPrsm;i<NumPrsm+NumTetras;i++)
    {
        vwgt[i]=4;
    }

    idx_t *epart=new idx_t[NumElement];

    idx_t *npart=new idx_t[NumVertices];

    idx_t objvalue=0;

    idx_t options[40];

    METIS_SetDefaultOptions(options);

    options[METIS_OPTION_NUMBERING]=0;

    if(nparts>1)
    {

        int results=METIS_PartMeshDual(&NumElement,&NumVertices,eptr,eind,vwgt,nullptr,&ncommon,&NumParts,nullptr,
                                       nullptr,&objvalue ,epart,npart);

    }
    else
    {
        for(int i=0;i<NumElement;i++)
        {
            epart[i]=0;
        }
        for(int i=0;i<NumVertices;i++)
        {
            npart[i]=0;
        }
    }


    for(int i=0;i<NumTetras;i++)
    {
      //  tetrasfile.pTetras[i].partMarker=epart[i]*stride+offset;//shift partMarker
        tetrasfile.pTetras[i].partMarker=epart[NumPrsm+i];
        tetrasfile.pTetras[i].index=NumPrsm+i;//reset the index
    }


    for(int i=0;i<NumPrsm;i++)
    {
      //  tetrasfile.pTetras[i].partMarker=epart[i]*stride+offset;//shift partMarker
        tetrasfile.pPrisms[i].partMarker=epart[i];
        tetrasfile.pPrisms[i].index=i;//reset the index
    }

    for(int i=0;i<NumVertices;i++)
    {
        tetrasfile.nodes[i].partMarker=npart[i];
    }

    if(eptr!=nullptr)
        delete []eptr;
    eptr=nullptr;
    if(eind!=nullptr)
        delete []eind;
    eind=nullptr;
    int tetrasID=-1;


    for(int i=0;i<tetrasfile.NumTris;i++)
    {
        tetrasID=tetrasfile.pTris[i].iCell;
        if(tetrasID<NumPrsm)
        tetrasfile.pTris[i].partMarker=tetrasfile.pPrisms[tetrasID].partMarker;
        else if(tetrasID<NumPrsm+NumTetras)
        tetrasfile.pTris[i].partMarker=tetrasfile.pTetras[tetrasID-NumPrsm].partMarker;
    }



//set the origin partMarker

    //insert part information to node
    for(int i=0;i<tetrasfile.NumPrsm;i++)
    {
        int pointID=-1;


        int procID=tetrasfile.pPrisms[i].partMarker;

        for(int j=0;j<6;j++)
        {
            pointID=tetrasfile.pPrisms[i].vertices[j];
            tetrasfile.nodes[pointID].procs.insert(procID);
        }
     //   pointID=tetrasfile.pTetras[i].vertices[0];
    }

    //insert part information to node
    for(int i=0;i<tetrasfile.NumTetras;i++)
    {
        int pointID=-1;


        int procID=tetrasfile.pTetras[i].partMarker;

        for(int j=0;j<4;j++)
        {
            pointID=tetrasfile.pTetras[i].vertices[j];
            tetrasfile.nodes[pointID].procs.insert(procID);
        }
     //   pointID=tetrasfile.pTetras[i].vertices[0];
    }


    //prepare the data to deliver



    //set<>
    set<int> *elementPart=new set<int>[nparts];

    set<int> *tempTerasPart=new set<int>[nparts];

    set<int> *tempPrismPart=new set<int>[nparts];

    set<int> *facetPart=new set<int>[nparts];

    int partID=-1;
    for(int i=0;i<NumElement;i++)
    {
        //   tetrasfile.elements[i].partMarker=epart[i];
        partID=epart[i];
        elementPart[partID].insert(i);
        if(vwgt[i]==6){
            tempPrismPart[partID].insert(i);
        }
        else if(vwgt[i]==4){

            tempTerasPart[partID].insert(i);
        }

    }

    for(int i=0;i<tetrasfile.NumTris;i++)
    {
        partID=tetrasfile.pTris[i].partMarker;
        facetPart[partID].insert(i);

    }



    int tempSize=-1;
    int tempID=-1;
    int tempMID=-1;

    //set<int>pointID;
    map<int,int>global_to_local;

    map<int,int>local_to_global;

    set<int>::iterator setIter;
    map<int,int>::iterator mapIter;
    int count=0;
    set<int> trin;
    set<int> tetn;
    for(int i=0;i<nparts;i++)
    {
        tempSize=tempTerasPart[i].size();
        tetrasPart[i].NumTetras=tempSize;
        tetrasPart[i].pTetras=new TETRAS[tempSize];
        tempSize=tempPrismPart[i].size();
        tetrasPart[i].NumPrsm=tempSize;
        tetrasPart[i].pPrisms=new PRISM[tempSize];

        count=0;
        //insert pointID
        global_to_local.clear();
        local_to_global.clear();
        for(setIter=tempPrismPart[i].begin();setIter!=tempPrismPart[i].end();setIter++)
        {
            tempID=*setIter;

            for(int j=0;j<6;j++)
            {
                tempMID=tetrasfile.pPrisms[tempID].vertices[j];
                tetn.insert(tempMID);
                mapIter=global_to_local.find(tempMID);
                if(mapIter!=global_to_local.end())
                {
                    //find it, do nothing
                }
                else
                {
                    //not found
                    global_to_local[tempMID]=count;
                    local_to_global[count]=tempMID;
                    count++;
                }
            }
        }
//        ofstream a;
//        a.open("fffff.txt");
//        for(int i=0;i<tetrasfile.NumTetras;++i){
//            for(int j=0;j<4;++j)
//            a<<tetrasfile.pTetras[i].vertices[j]<<" ";
//            a<<endl;
//        }
//        a.close();
        for(setIter=tempTerasPart[i].begin();setIter!=tempTerasPart[i].end();setIter++)
        {
            tempID=*setIter-NumPrsm;
            for(int j=0;j<4;j++)
            {

                tempMID=tetrasfile.pTetras[tempID].vertices[j];
                tetn.insert(tempMID);
                mapIter=global_to_local.find(tempMID);
                if(mapIter!=global_to_local.end())
                {
                    //find it, do nothing
                }
                else
                {
                    //not found
                    global_to_local[tempMID]=count;
                    local_to_global[count]=tempMID;
                //   cout<<tempMID<<endl;
                    count++;
                }
            }
        }
//        for(int i=0;i<global_to_local.size();++i){
//            cout<<global_to_local[i]<<endl;
//        }
        tetrasPart[i].NumNodes=count;
        tetrasPart[i].nodes=new Node[tetrasPart[i].NumNodes];
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            mapIter=local_to_global.find(j);
            tempMID=mapIter->second;
            tetrasPart[i].nodes[j]=tetrasfile.nodes[tempMID];
            tetrasPart[i].nodes[j].localID=j;
        }
        count=0;
        //insert elementID
        for(setIter=tempPrismPart[i].begin();setIter!=tempPrismPart[i].end();setIter++)
        {
            tempID=*setIter;
            tetrasPart[i].pPrisms[count]=tetrasfile.pPrisms[tempID];

            tetrasPart[i].pPrisms[count].localID=count;

          //  tetrasPart[i].pPrisms[count].index=

            int tempPartMarker=tetrasPart[i].pPrisms[count].partMarker;

            tetrasPart[i].pPrisms[count].partMarker=tempPartMarker*stride+offset;//////////////shift the partMarker


            for(int j=0;j<6;j++)
            {
                tempMID=tetrasPart[i].pPrisms[count].vertices[j];
                tetrasPart[i].pPrisms[count].vertices[j]=global_to_local[tempMID];
            }
            count++;

        }

        for(setIter=tempTerasPart[i].begin();setIter!=tempTerasPart[i].end();setIter++)
        {
            tempID=*setIter;
          tetrasPart[i].pTetras[count-tetrasPart[i].NumPrsm]=tetrasfile.pTetras[tempID-tetrasfile.NumPrsm];
          tetrasPart[i].pTetras[count-tetrasPart[i].NumPrsm].localID=count;

            //tetrasPart[i].pTetras[count].index=

            int tempPartMarker=tetrasPart[i].pTetras[count-tetrasPart[i].NumPrsm].partMarker;

            tetrasPart[i].pTetras[count-tetrasPart[i].NumPrsm].partMarker=tempPartMarker*stride+offset;//////////////shift the partMarker

            for(int j=0;j<4;j++)
            {

                tempMID=tetrasPart[i].pTetras[count-tetrasPart[i].NumPrsm].vertices[j];

                tetrasPart[i].pTetras[count-tetrasPart[i].NumPrsm].vertices[j]=global_to_local[tempMID];
            }
            count++;

        }
//        cout<<"tempSize: "<<tempSize<<" count:"<<count<<endl;
        assert(count==tempSize);
        //insert facets
        tempSize=facetPart[i].size();
        tetrasPart[i].NumTris=tempSize;
        tetrasPart[i].pTris=new TRI[tetrasPart[i].NumTris];
        count=0;
        for(setIter=facetPart[i].begin();setIter!=facetPart[i].end();setIter++)
        {
            tempID=*setIter;

            tetrasPart[i].pTris[count]=tetrasfile.pTris[tempID];
      //      tetrasPart[i].pTris[count].index=count;///////////////////////////////annotated by zhvliu

            int tempPartMarker=tetrasPart[i].pTris[count].partMarker;

            tetrasPart[i].pTris[count].partMarker=tempPartMarker*stride+offset;//////////////shift the partMarker

            for(int j=0;j<3;j++)
            {

                tempMID=tetrasPart[i].pTris[count].vertices[j];
//                                if(tetn.find(tempMID)==tetn.end())
//                                {
//                                    cout<<"errorororo"<<endl;
//                                    cout<<global_to_local[tempMID]<<endl;
//                                }
//                cout<<tempMID<<"  "<<global_to_local[tempMID]<<endl;
                tetrasPart[i].pTris[count].vertices[j]=global_to_local[tempMID];//here have problem

            }
            count++;
        }
    }

    if(elementPart!=nullptr)
    {
        delete []elementPart;
    }
    if(facetPart!=nullptr)
    {
        delete []facetPart;
    }


    for(int i=0;i<NumElement;i++)
    {
        if(i<NumPrsm)
        tetrasfile.pPrisms[i].partMarker=epart[i]*stride+offset;//shift partMarker
        else
        tetrasfile.pTetras[i-NumPrsm].partMarker=epart[i]*stride+offset;//shift partMarker

      //  tetrasfile.pTetras[i].localID=i;
      //  tetrasfile.pTetras[i].partMarker=epart[i];
    }

/*
    char testName[256];
    sprintf(testName,"%s%d%s","second",offset,".vtk");
    writeVTKFile(testName,tetrasfile);


    sprintf(testName,"%s%d%s","secondTriangle",offset,".vtk");
    writeTriangleVTKFile(testName,tetrasfile);

*/
    if(epart!=nullptr)
        delete []epart;
    epart=nullptr;
    if(npart!=nullptr)
        delete []npart;
    npart=nullptr;

    map<string,int64_t> tri_globalID;
    for(int i=0;i<nparts;i++)
    {
        constructFacets_test(tetrasPart[i],tetrasfile,tri_globalID);
        char name[256];
        sprintf(name,"%s%d%s","cubeInter",i,".vtk");
       // writeTriangleVTKFile(name,tetrasPart[i]);
    }


    if(type==1)//which means the second partition
    {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
       // map<string,int>tri_globalID;

        cout<<"Entering type==1"<<endl;

        int originFacets=tetrasfile.NumTris;
        int sumFacets=0;
        for(int i=0;i<nparts;i++)
        {
            sumFacets+=tetrasPart[i].NumTris;
        }
        cout<<"originFacets: "<<originFacets<<" sumFacets: "<<sumFacets<<endl;
        int nNewLocalFacets=(sumFacets-originFacets)/2;
        cout<<"nNewLocalFacets: "<<nNewLocalFacets<<" tri_gloablID.size= "<<tri_globalID.size()<<endl;
        assert(nNewLocalFacets==tri_globalID.size());



        int64_t offset = nNewLocalFacets;
        int64_t nNewGlobalFacets = nNewLocalFacets;
        if (comm_size > 1)
        {
            MPI_Scan(&(static_cast<const int64_t&>(nNewLocalFacets)), &offset, 1,
                     MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&nNewLocalFacets, &nNewGlobalFacets, 1,
                          MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        }
        offset -= nNewLocalFacets;
        offset += tetrasfile.NumUniqueSurfFacets + tetrasfile.NumUniqueInterfFacets;
        for (int i = 0; i < nparts; ++i)
        {
            tetrasPart[i].NumUniqueSurfFacets = tetrasfile.NumUniqueSurfFacets;
            tetrasPart[i].NumUniqueInterfFacets = tetrasfile.NumUniqueInterfFacets + nNewGlobalFacets;
        }

        map<string,int64_t>::iterator mapIter;
        count=0;
        for(mapIter=tri_globalID.begin();mapIter!=tri_globalID.end();mapIter++)
        {
            mapIter->second=count+offset;
            count++;
        }


        for(int i=0;i<nparts;i++)
        {
            updateTriIndex(tetrasPart[i],tri_globalID);
     //       findiCellFast(tetrasPart[i]);///////find icell
        }



    }

    /*
    const int nUniqueFacets = refinedMesh.numOfFacets() - nDuplicatedFacets;
    int offset = nUniqueFacets;
    MPI_Scan(&nUniqueFacets, &offset, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(const_cast<int*>(&nUniqueFacets), &refinedMesh.NumUniqueFacets,
                  1, MPI_INT, MPI_SUM, comm);

    if (rank == 0)
        printf("##: nTotalUniqueFacets = %d\n", refinedMesh.NumUniqueFacets);

    offset -= nUniqueFacets;

    */
    return 1;
}

int partition(HYBRID_MESH&tetrasfile, HYBRID_MESH *tetrasPart, int nparts, int stride, int offset, int type)
{

    TETRAS *ptetras=tetrasfile.pTetras;

    idx_t NumVertices=tetrasfile.NumNodes;

    idx_t NumElement=tetrasfile.NumTetras;

    idx_t ncommon=3;

    idx_t NumParts=nparts;

    idx_t *eptr=new idx_t[NumElement+1];

    idx_t *eind=new idx_t[NumElement*4];

    for(int i=0;i<NumElement+1;i++)
    {
        eptr[i]=i*4;
    }

    for(int i=0;i<NumElement;i++)
    {
        for(int j=0;j<4;j++)
        {
            eind[i*4+j]=ptetras[i].vertices[j];
        }
    }

    idx_t *epart=new idx_t[NumElement];

    idx_t *npart=new idx_t[NumVertices];

    idx_t objvalue=0;

    idx_t options[40];

    METIS_SetDefaultOptions(options);

    options[METIS_OPTION_NUMBERING]=0;


    if(nparts>1)
    {
        int results=METIS_PartMeshDual(&NumElement,&NumVertices,eptr,eind,nullptr,nullptr,&ncommon,&NumParts,nullptr,
                                       nullptr,&objvalue ,epart,npart);
    }
    else
    {
        for(int i=0;i<NumElement;i++)
        {
            epart[i]=0;
        }
        for(int i=0;i<NumVertices;i++)
        {
            npart[i]=0;
        }
    }

    for(int i=0;i<NumElement;i++)
    {
      //  tetrasfile.pTetras[i].partMarker=epart[i]*stride+offset;//shift partMarker
        tetrasfile.pTetras[i].partMarker=epart[i];
        tetrasfile.pTetras[i].index=i;//reset the index
    }
    for(int i=0;i<NumVertices;i++)
    {
        tetrasfile.nodes[i].partMarker=npart[i];
    }

    if(eptr!=nullptr)
        delete []eptr;
    eptr=nullptr;
    if(eind!=nullptr)
        delete []eind;
    eind=nullptr;


    int tetrasID=-1;
    for(int i=0;i<tetrasfile.NumTris;i++)
    {
        tetrasID=tetrasfile.pTris[i].iCell;
        tetrasfile.pTris[i].partMarker=tetrasfile.pTetras[tetrasID].partMarker;
    }

//set the origin partMarker

    //insert part information to node
    for(int i=0;i<tetrasfile.NumTetras;i++)
    {
        int pointID=-1;


        int procID=tetrasfile.pTetras[i].partMarker;

        for(int j=0;j<4;j++)
        {
            pointID=tetrasfile.pTetras[i].vertices[j];
            tetrasfile.nodes[pointID].procs.insert(procID);
        }
     //   pointID=tetrasfile.pTetras[i].vertices[0];
    }

    //prepare the data to deliver



    //set<>

    set<int> *elementPart=new set<int>[nparts];

    set<int> *facetPart=new set<int>[nparts];

    int partID=-1;
    for(int i=0;i<NumElement;i++)
    {
        //   tetrasfile.elements[i].partMarker=epart[i];
        partID=epart[i];
        elementPart[partID].insert(i);

    }
    for(int i=0;i<tetrasfile.NumTris;i++)
    {
        partID=tetrasfile.pTris[i].partMarker;
        facetPart[partID].insert(i);
    }



    int tempSize=-1;
    int tempID=-1;
    int tempMID=-1;

    //set<int>pointID;
    map<int,int>global_to_local;

    map<int,int>local_to_global;

    set<int>::iterator setIter;
    map<int,int>::iterator mapIter;
    int count=0;
    for(int i=0;i<nparts;i++)
    {
        tempSize=elementPart[i].size();
        tetrasPart[i].NumTetras=tempSize;
        tetrasPart[i].pTetras=new TETRAS[tempSize];

        count=0;
        //insert pointID
        global_to_local.clear();
        local_to_global.clear();
        for(setIter=elementPart[i].begin();setIter!=elementPart[i].end();setIter++)
        {
            tempID=*setIter;

            for(int j=0;j<4;j++)
            {
                tempMID=tetrasfile.pTetras[tempID].vertices[j];
                mapIter=global_to_local.find(tempMID);
                if(mapIter!=global_to_local.end())
                {
                    //find it, do nothing
                }
                else
                {
                    //not found
                    global_to_local[tempMID]=count;
                    local_to_global[count]=tempMID;
                    count++;
                }
            }
        }
        tetrasPart[i].NumNodes=count;
        tetrasPart[i].nodes=new Node[tetrasPart[i].NumNodes];

        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            mapIter=local_to_global.find(j);
            tempMID=mapIter->second;
            tetrasPart[i].nodes[j]=tetrasfile.nodes[tempMID];
            tetrasPart[i].nodes[j].localID=j;
        }

        count=0;
        //insert elementID
        for(setIter=elementPart[i].begin();setIter!=elementPart[i].end();setIter++)
        {
            tempID=*setIter;
            tetrasPart[i].pTetras[count]=tetrasfile.pTetras[tempID];

            tetrasPart[i].pTetras[count].localID=count;

            //tetrasPart[i].pTetras[count].index=

            int tempPartMarker=tetrasPart[i].pTetras[count].partMarker;

            tetrasPart[i].pTetras[count].partMarker=tempPartMarker*stride+offset;//////////////shift the partMarker


            for(int j=0;j<4;j++)
            {
                tempMID=tetrasPart[i].pTetras[count].vertices[j];
                tetrasPart[i].pTetras[count].vertices[j]=global_to_local[tempMID];
            }
            count++;

        }
   //     cout<<"tempSize: "<<tempSize<<" count:"<<count<<endl;
        assert(count==tempSize);
        //insert facets
        tempSize=facetPart[i].size();
        tetrasPart[i].NumTris=tempSize;
        tetrasPart[i].pTris=new TRI[tetrasPart[i].NumTris];
        count=0;
        for(setIter=facetPart[i].begin();setIter!=facetPart[i].end();setIter++)
        {
            tempID=*setIter;
            tetrasPart[i].pTris[count]=tetrasfile.pTris[tempID];
      //      tetrasPart[i].pTris[count].index=count;///////////////////////////////annotated by zhvliu

            int tempPartMarker=tetrasPart[i].pTris[count].partMarker;

            tetrasPart[i].pTris[count].partMarker=tempPartMarker*stride+offset;//////////////shift the partMarker

            for(int j=0;j<3;j++)
            {
                tempMID=tetrasPart[i].pTris[count].vertices[j];
                tetrasPart[i].pTris[count].vertices[j]=global_to_local[tempMID];
            }
            count++;
        }
    }

    if(elementPart!=nullptr)
    {
        delete []elementPart;
    }
    if(facetPart!=nullptr)
    {
        delete []facetPart;
    }


    for(int i=0;i<NumElement;i++)
    {
        tetrasfile.pTetras[i].partMarker=epart[i]*stride+offset;//shift partMarker

      //  tetrasfile.pTetras[i].localID=i;
      //  tetrasfile.pTetras[i].partMarker=epart[i];
    }

/*
    char testName[256];
    sprintf(testName,"%s%d%s","second",offset,".vtk");
    writeVTKFile(testName,tetrasfile);


    sprintf(testName,"%s%d%s","secondTriangle",offset,".vtk");
    writeTriangleVTKFile(testName,tetrasfile);

*/
    if(epart!=nullptr)
        delete []epart;
    epart=nullptr;
    if(npart!=nullptr)
        delete []npart;
    npart=nullptr;

    map<string,int64_t> tri_globalID;
    for(int i=0;i<nparts;i++)
    {
        constructFacets(tetrasPart[i],tetrasfile,tri_globalID);
        char name[256];
        sprintf(name,"%s%d%s","cubeInter",i,".vtk");
       // writeTriangleVTKFile(name,tetrasPart[i]);
    }


    if(type==1)//which means the second partition
    {
		int comm_size;
		MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
       // map<string,int>tri_globalID;

        cout<<"Entering type==1"<<endl;

        int originFacets=tetrasfile.NumTris;
        int sumFacets=0;
        for(int i=0;i<nparts;i++)
        {
            sumFacets+=tetrasPart[i].NumTris;
        }
        cout<<"originFacets: "<<originFacets<<" sumFacets: "<<sumFacets<<endl;
        int nNewLocalFacets=(sumFacets-originFacets)/2;
        cout<<"nNewLocalFacets: "<<nNewLocalFacets<<" tri_gloablID.size= "<<tri_globalID.size()<<endl;
        assert(nNewLocalFacets==tri_globalID.size());



        int64_t offset = nNewLocalFacets;
        int64_t nNewGlobalFacets = nNewLocalFacets;
		if (comm_size > 1)
		{
			MPI_Scan(&(static_cast<const int64_t&>(nNewLocalFacets)), &offset, 1,
					 MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&nNewLocalFacets, &nNewGlobalFacets, 1,
						  MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
		}
        offset -= nNewLocalFacets;
		offset += tetrasfile.NumUniqueSurfFacets + tetrasfile.NumUniqueInterfFacets;
        for (int i = 0; i < nparts; ++i)
		{
			tetrasPart[i].NumUniqueSurfFacets = tetrasfile.NumUniqueSurfFacets;
			tetrasPart[i].NumUniqueInterfFacets = tetrasfile.NumUniqueInterfFacets + nNewGlobalFacets;
		}

        map<string,int64_t>::iterator mapIter;
        count=0;
        for(mapIter=tri_globalID.begin();mapIter!=tri_globalID.end();mapIter++)
        {
            mapIter->second=count+offset;
            count++;
        }


        for(int i=0;i<nparts;i++)
        {
            updateTriIndex(tetrasPart[i],tri_globalID);
     //       findiCellFast(tetrasPart[i]);///////find icell
        }



    }

    /*
    const int nUniqueFacets = refinedMesh.numOfFacets() - nDuplicatedFacets;
    int offset = nUniqueFacets;
    MPI_Scan(&nUniqueFacets, &offset, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(const_cast<int*>(&nUniqueFacets), &refinedMesh.NumUniqueFacets,
                  1, MPI_INT, MPI_SUM, comm);

    if (rank == 0)
        printf("##: nTotalUniqueFacets = %d\n", refinedMesh.NumUniqueFacets);

    offset -= nUniqueFacets;

    */
    return 1;
}

int updateTriIndex(HYBRID_MESH &mesh,map<string,int64_t>&tri_globalID)
{
    int vertices[3];
    map<string,int64_t>::iterator iter;
    int count=0;
    for(int i=0;i<mesh.NumTris;i++)
    {
        for(int j=0;j<3;j++)
        {
            vertices[j]=mesh.pTris[i].vertices[j];
            vertices[j]=mesh.nodes[vertices[j]].index;
        }
        sort(vertices,vertices+3);
        string temp=IntToString(vertices[0])+"_"+IntToString(vertices[1])+"_"+IntToString(vertices[2]);
        iter=tri_globalID.find(temp);
        if(iter==tri_globalID.end())
        {
            //not found , which means it is not the new interface
        }
        else
        {
            //found it
            mesh.pTris[i].index=iter->second;
            count++;
        }

    }
    return 1;
}

int meshRefining_test(HYBRID_MESH &tetrasfile,HYBRID_MESH &newTetrasfile,int partMarker)
{
#ifdef DEBUG
    set<int> sum;
//    cout<<tetrasfile.NumTetras<<endl;
    for(int i=0;i<tetrasfile.NumPrsm;++i){
        for(int j=0;j<6;++j)
        sum.insert(tetrasfile.pPrisms[i].vertices[j]);
    }
    for(int i=0;i<tetrasfile.NumTetras;++i){
        for(int j=0;j<4;++j)
        sum.insert(tetrasfile.pTetras[i].vertices[j]);

    }
    cout<<"sum  cell===="<<sum.size()<<endl;
    for(int i=0;i<tetrasfile.NumTris;++i){
        for(int j=0;j<3;++j)
        sum.insert(tetrasfile.pTris[i].vertices[j]);
    }
    cout<<"sum  all====="<<sum.size()<<endl;
#endif
    //HYBRID_MESH newTetrasfile;
    //delete the self Procs
    set<int>::iterator setIter;
    for(int i=0;i<tetrasfile.NumNodes;i++)
    {
     //   set<int> &procs = tetrasfile.nodes[i].procs;
       // for(setIter=tetrasfile.nodes[i].procs.begin();setIter!=tetrasfile.nodes[i].procs.end();)
        tetrasfile.nodes[i].procs.erase(partMarker);
        tetrasfile.nodes[i].NumProcs=tetrasfile.nodes[i].procs.size();
        // assert(tetrasfile.nodes[i].procs.size() < 2);
    }


    map<string,newNode>edgeHash;
    map<string,newNode>faceHash;

    edgeHash.clear();
    int NumNodes=tetrasfile.NumNodes;
    int NumNewNodes=NumNodes;
    int NumElements=tetrasfile.numOfCells();

    int edgeID[6][2]=
    {
        0,1,
        0,2,
        0,3,
        1,2,
        1,3,
        2,3
    };

    int prism_edgeID[9][2]=
    {
        0,1,//0
        0,2,//1
        0,3,//2
        1,2,//3
        1,4,//4
        2,5,//5
        3,4,//6
        3,5,//7
        4,5//8
    };
    int prism_faceID[3][4]=
    {
      0,3,5,2,//0
      0,1,4,3,//1
      1,2,5,4 //2
    };

    int Estart=-1;
    int Eend=-1;
    int vec[4];
    int startID=-1;
    int endID=-1;
    int mid=-1;
    string temp;
    pair<string,newNode> edgePair;
    pair<string,newNode> facePair;
    for(int i=0;i<NumElements;i++)
    {
        if(i<tetrasfile.NumPrsm){
            for(int j=0;j<3;j++){
                for(int k=0;k<4;++k){
                    vec[k]=tetrasfile.pPrisms[i].vertices[prism_faceID[j][k]];
                }
                sort(vec,vec+4);
                temp=IntToString(vec[0])+"_"+IntToString(vec[1])+"_"+IntToString(vec[2])+"_"+IntToString(vec[3]);
                if(faceHash.find(temp)!=faceHash.end())
                {
                    //find it, do nothing
                }
                else
                {
                    //not found
                    facePair.first=temp;
                    newNode nodetemp;
                    nodetemp.localID=NumNewNodes;
                    NumNewNodes++;
                    nodetemp.coord=tetrasfile.nodes[vec[0]].coord+tetrasfile.nodes[vec[1]].coord+tetrasfile.nodes[vec[2]].coord+tetrasfile.nodes[vec[3]].coord;
                    nodetemp.coord=nodetemp.coord/4;
                    facePair.second=nodetemp;

                    faceHash.insert(facePair);
                }
            }
            for(int j=0;j<9;j++)
            {
                Estart=prism_edgeID[j][0];
                Eend=prism_edgeID[j][1];
                startID=tetrasfile.pPrisms[i].vertices[Estart];
                endID=tetrasfile.pPrisms[i].vertices[Eend];
                if(startID>endID)
                {
                    mid=startID;
                    startID=endID;
                    endID=mid;
                }
                temp=IntToString(startID)+"_"+IntToString(endID);
                if(edgeHash.find(temp)!=edgeHash.end())
                {
                    //find it, do nothing
                }
                else
                {
                    //not found
                    edgePair.first=temp;
                    newNode nodetemp;
                    nodetemp.localID=NumNewNodes;
                    NumNewNodes++;
                    nodetemp.coord=tetrasfile.nodes[startID].coord+tetrasfile.nodes[endID].coord;
                    nodetemp.coord=nodetemp.coord/2;
                    edgePair.second=nodetemp;
                    edgeHash.insert(edgePair);
                }
            }
        }
       else{
        for(int j=0;j<6;j++)
        {
            Estart=edgeID[j][0];
            Eend=edgeID[j][1];
            startID=tetrasfile.pTetras[i-tetrasfile.NumPrsm].vertices[Estart];
            endID=tetrasfile.pTetras[i-tetrasfile.NumPrsm].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            if(edgeHash.find(temp)!=edgeHash.end())
            {
                //find it, do nothing
            }
            else
            {
                //not found
                edgePair.first=temp;
                newNode nodetemp;
                nodetemp.localID=NumNewNodes;
                NumNewNodes++;
                nodetemp.coord=tetrasfile.nodes[startID].coord+tetrasfile.nodes[endID].coord;
                nodetemp.coord=nodetemp.coord/2;
                edgePair.second=nodetemp;

                edgeHash.insert(edgePair);
            }
        }
}
    }

    newTetrasfile.NumNodes=NumNewNodes;
    newTetrasfile.nodes=new Node[NumNewNodes];

    //---------------//锟斤拷锟斤拷锟铰憋拷锟斤拷锟斤拷锟斤拷锟侥斤拷锟斤拷锟斤拷息 //---------------///
    for(int i=0;i<NumNodes;i++)
    {
        newTetrasfile.nodes[i]=tetrasfile.nodes[i];//local ID
    }
    map<string,newNode>::iterator iter;

    newNode midtemp;
    int indextemp=-1;
    for(iter=edgeHash.begin();iter!=edgeHash.end();iter++)
    {
        midtemp=iter->second;
        indextemp=midtemp.localID;
        newTetrasfile.nodes[indextemp].index=-1;
        newTetrasfile.nodes[indextemp].localID=indextemp;
        newTetrasfile.nodes[indextemp].coord=midtemp.coord;
    }

    for(iter=faceHash.begin();iter!=faceHash.end();iter++)
    {
        midtemp=iter->second;
        indextemp=midtemp.localID;
        newTetrasfile.nodes[indextemp].index=-1;
        newTetrasfile.nodes[indextemp].localID=indextemp;
        newTetrasfile.nodes[indextemp].coord=midtemp.coord;
    }



    //---------------//锟斤拷锟斤拷锟铰憋拷锟斤拷锟斤拷锟斤拷锟侥碉拷元锟斤拷息 //---------------///

    int NumNewElements=NumElements*8;//one element are divide into 8 parts
    newTetrasfile.NumTetras=tetrasfile.NumTetras*8;
    newTetrasfile.NumPrsm=tetrasfile.NumPrsm*8;
    newTetrasfile.pTetras=new TETRAS[tetrasfile.NumTetras*8];
    newTetrasfile.pPrisms=new PRISM[tetrasfile.NumPrsm*8];




//    sixNodesPattern(tetrasfile,newTetrasfile,edgeHash);
    sixNodesPattern_test(tetrasfile,newTetrasfile,edgeHash,faceHash);
    printf("oldelems=%d;newelems=%d \n",NumElements,NumNewElements);

    //find the iCell of the triangles
    setupCellNeig_test(newTetrasfile.NumNodes,newTetrasfile.NumTetras,&newTetrasfile);//local

//    findiCellFast(newTetrasfile);





    return 1;
}


int meshRefining(HYBRID_MESH &tetrasfile,HYBRID_MESH &newTetrasfile,int partMarker)
{
    //HYBRID_MESH newTetrasfile;
    //delete the self Procs
    set<int>::iterator setIter;
    for(int i=0;i<tetrasfile.NumNodes;i++)
    {
     //   set<int> &procs = tetrasfile.nodes[i].procs;
       // for(setIter=tetrasfile.nodes[i].procs.begin();setIter!=tetrasfile.nodes[i].procs.end();)
        tetrasfile.nodes[i].procs.erase(partMarker);
        tetrasfile.nodes[i].NumProcs=tetrasfile.nodes[i].procs.size();
        // assert(tetrasfile.nodes[i].procs.size() < 2);
    }


    map<string,newNode>edgeHash;

    edgeHash.clear();
    int NumNodes=tetrasfile.NumNodes;
    int NumNewNodes=NumNodes;
    int NumElements=tetrasfile.NumTetras;

    int edgeID[6][2]=
    {
        0,1,
        0,2,
        0,3,
        1,2,
        1,3,
        2,3
    };

    int Estart=-1;
    int Eend=-1;
    int startID=-1;
    int endID=-1;
    int mid=-1;
    string temp;
    pair<string,newNode> edgePair;
    for(int i=0;i<NumElements;i++)
    {
        for(int j=0;j<6;j++)
        {
            Estart=edgeID[j][0];
            Eend=edgeID[j][1];
            startID=tetrasfile.pTetras[i].vertices[Estart];
            endID=tetrasfile.pTetras[i].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            if(edgeHash.find(temp)!=edgeHash.end())
            {
                //find it, do nothing
            }
            else
            {
                //not found
                edgePair.first=temp;
                newNode nodetemp;
                nodetemp.localID=NumNewNodes;
                NumNewNodes++;
                nodetemp.coord=tetrasfile.nodes[startID].coord+tetrasfile.nodes[endID].coord;
                nodetemp.coord=nodetemp.coord/2;
                edgePair.second=nodetemp;

                edgeHash.insert(edgePair);
            }
        }

    }

    newTetrasfile.NumNodes=NumNewNodes;
    newTetrasfile.nodes=new Node[NumNewNodes];

    //---------------//锟斤拷锟斤拷锟铰憋拷锟斤拷锟斤拷锟斤拷锟侥斤拷锟斤拷锟斤拷息 //---------------///
    for(int i=0;i<NumNodes;i++)
    {
        newTetrasfile.nodes[i]=tetrasfile.nodes[i];//local ID
    }
    map<string,newNode>::iterator iter;

    newNode midtemp;
    int indextemp=-1;
    for(iter=edgeHash.begin();iter!=edgeHash.end();iter++)
    {
        midtemp=iter->second;
        indextemp=midtemp.localID;
        newTetrasfile.nodes[indextemp].index=-1;
        newTetrasfile.nodes[indextemp].localID=indextemp;
        newTetrasfile.nodes[indextemp].coord=midtemp.coord;
    }



    //---------------//锟斤拷锟斤拷锟铰憋拷锟斤拷锟斤拷锟斤拷锟侥碉拷元锟斤拷息 //---------------///

    int NumNewElements=NumElements*8;//one element are divide into 8 parts
    newTetrasfile.NumTetras=NumNewElements;
    newTetrasfile.pTetras=new TETRAS[NumNewElements];





    sixNodesPattern(tetrasfile,newTetrasfile,edgeHash);


    printf("oldelems=%d;newelems=%d \n",NumElements,NumNewElements);


    //find the iCell of the triangles
    setupCellNeig(newTetrasfile.NumNodes,newTetrasfile.NumTetras,newTetrasfile.pTetras);//local

//    findiCellFast(newTetrasfile);





    return 1;
}
int meshRefining(HYBRID_MESH &tetrasfile,HYBRID_MESH &newTetrasfile,int partMarker,Refletion &ref)
{
//    for(int i=0;i<tetrasfile.NumTris;++i)
//    cout<<tetrasfile.pTris[i].iSurface<<endl;
    //delete the self Procs
//    for(int i=0;i<tetrasfile.NumTris;++i){
//       cout<< tetrasfile.pTris[i].iSurface<<endl;
//    }
    cout<<"meshrefining start!!"<<endl;
    map<int,int> reset;
    int k=0;
    for(auto i=ref.subject_table.begin();i!=ref.subject_table.end();++i){
        reset[k++]=i->second;
   }
    ref.subject_table=reset;


    for(int i=0;i<tetrasfile.NumNodes;i++)
    {
     //   set<int> &procs = tetrasfile.nodes[i].procs;
       // for(setIter=tetrasfile.nodes[i].procs.begin();setIter!=tetrasfile.nodes[i].procs.end();)
        tetrasfile.nodes[i].procs.erase(partMarker);
        tetrasfile.nodes[i].NumProcs=tetrasfile.nodes[i].procs.size();
        // assert(tetrasfile.nodes[i].procs.size() < 2);
    }
    map<string,newNode>edgeHash;

    edgeHash.clear();
    int NumNodes=tetrasfile.NumNodes;
    int NumNewNodes=NumNodes;
    int NumElements=tetrasfile.NumTetras;

    int edgeID[6][2]=
    {
        0,1,
        0,2,
        0,3,
        1,2,
        1,3,
        2,3
    };
    vector<string> lines;
    unordered_map<string,vector<int>> lines_map;
    string line_1;
    string line_2;
    string line_3;
   int v[3];

   lines.clear();
   lines_map.clear();

   for(int i=0;i<tetrasfile.NumTris;i++)
   {
       v[0]=tetrasfile.pTris[i].vertices[0];
       v[1]=tetrasfile.pTris[i].vertices[1];
       v[2]=tetrasfile.pTris[i].vertices[2];
       sort(v,v+3);
       line_1=IntToString(v[0])+"_"+IntToString(v[1]);
       line_2=IntToString(v[1])+"_"+IntToString(v[2]);
       line_3=IntToString(v[0])+"_"+IntToString(v[2]);
       lines.push_back(line_1);
       lines.push_back(line_2);
       lines.push_back(line_3);
       lines_map[line_1].push_back(i);
       lines_map[line_2].push_back(i);
       lines_map[line_3].push_back(i);
   }
 cout<<"map established!!!!!!!"<<endl;
    int Estart=-1;
    int Eend=-1;
    int startID=-1;
    int endID=-1;
    int mid=-1;
    string temp;
    pair<string,newNode> edgePair;
    for(int i=0;i<NumElements;i++)
    {
        for(int j=0;j<6;j++)
        {
            Estart=edgeID[j][0];
            Eend=edgeID[j][1];
            startID=tetrasfile.pTetras[i].vertices[Estart];
            endID=tetrasfile.pTetras[i].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            if(edgeHash.find(temp)!=edgeHash.end())
            {
                //find it, do nothing
            }
            else
            {
                //not found

                edgePair.first=temp;
                newNode nodetemp;
                nodetemp.localID=NumNewNodes;
                NumNewNodes++;
                nodetemp.coord=tetrasfile.nodes[startID].coord+tetrasfile.nodes[endID].coord;
                nodetemp.coord=nodetemp.coord/2;
                int face_id_1=-1;
                int face_id_2=-1;
                int patch_id_1=-1;
                int patch_id_2=-1;
                double x,y,z;


                vector<string>::iterator strIter;
                vector<int> face_id;
//                strIter=lines.begin();
//                   while(strIter!=lines.end()){
//                       strIter=find(strIter,lines.end(),temp);
//                       if(strIter!=lines.end()){
//                           face_id.push_back((strIter-lines.begin())/3);//problem?strIter
//                           strIter++;
//                       }
//                   }//#add    效率极低
                   face_id=lines_map[temp];


            if(face_id.size()>1){
            for(int k=0;k<face_id.size()-1;++k){
                if(ref.subject_table[face_id[k]]!=ref.subject_table[face_id[k+1]]){
                    face_id_1=face_id[k];
                    face_id_2=face_id[k+1];
                    break;
                }
                face_id_1=face_id[k];
                face_id_2=face_id[k];
            }
            }    //simple judge
         //#add
                if(face_id_1!=-1&&face_id_2!=-1){
                    patch_id_1=ref.GetPatch_ID(face_id_1);
                    patch_id_2=ref.GetPatch_ID(face_id_2);

                }
                x=nodetemp.coord.x;
                y=nodetemp.coord.y;
                z=nodetemp.coord.z;
                if(face_id_1!=-1&&face_id_2!=-1){
                patch_id_1=tetrasfile.pTris[face_id_1].iSurface;
                patch_id_2=tetrasfile.pTris[face_id_2].iSurface;
                }
                if(patch_id_2!=-1&&patch_id_1!=-1){
//                nodetemp.coord=ref.subject_patch_id(patch_id_2,patch_id_1,nodetemp.coord);//#add
                  ref.reflection(patch_id_1,patch_id_2,x,y,z);
                    nodetemp.coord.x=x;
                    nodetemp.coord.y=y;
                    nodetemp.coord.z=z;
                }
                face_id_1=-1;
                face_id_2=-1;
                patch_id_1=-1;
                patch_id_2=-1;
                edgePair.second=nodetemp;
                edgeHash.insert(edgePair);
            }
        }
   //     cout<<i<<endl;
     }

     newTetrasfile.NumNodes=NumNewNodes;
    newTetrasfile.nodes=new Node[NumNewNodes];
   // cout<<"test here"<<endl;
    //---------------//锟斤拷锟斤拷锟铰憋拷锟斤拷锟斤拷锟斤拷锟侥斤拷锟斤拷锟斤拷息 //---------------///
    for(int i=0;i<NumNodes;i++)
    {
        newTetrasfile.nodes[i]=tetrasfile.nodes[i];//local ID
    }
    map<string,newNode>::iterator _iter;

    newNode midtemp;
    int indextemp=-1;
    for(_iter=edgeHash.begin();_iter!=edgeHash.end();_iter++)
    {
        midtemp=_iter->second;
        indextemp=midtemp.localID;
        newTetrasfile.nodes[indextemp].index=-1;
        newTetrasfile.nodes[indextemp].localID=indextemp;
        newTetrasfile.nodes[indextemp].coord=midtemp.coord;
    }



    //---------------//锟斤拷锟斤拷锟铰憋拷锟斤拷锟斤拷锟斤拷锟侥碉拷元锟斤拷息 //---------------///

    int NumNewElements=NumElements*8;//one element are divide into 8 parts
    newTetrasfile.NumTetras=NumNewElements;
    newTetrasfile.pTetras=new TETRAS[NumNewElements];





    sixNodesPattern(tetrasfile,newTetrasfile,edgeHash, ref);

    printf("oldelems=%d;newelems=%d \n",NumElements,NumNewElements);


    //find the iCell of the triangles
    setupCellNeig(newTetrasfile.NumNodes,newTetrasfile.NumTetras,newTetrasfile.pTetras);//local

 //  findiCellFast(newTetrasfile);
    writeVTKFile("TEST.vtk", newTetrasfile);

    return 1;
}

int splitFacet(HYBRID_MESH&oldmesh, HYBRID_MESH &newmesh)
{
    return 1;
}

int sixNodesPattern(HYBRID_MESH &oldmesh, HYBRID_MESH &newmesh, map<string, newNode> &edgeHash,Refletion &ref)
{
    int edgeID[6][2]=
    {
        0,1,//0
        0,2,//1
        0,3,//2
        1,2,//3
        1,3,//4
        2,3//5
    };
    int newPointID[6];//corresponded with edgeID

    int Estart=-1;
    int Eend=-1;
    int startID=-1;
    int endID=-1;

    int count=0;

    TETRAS tetrasTemp;

    map<string,newNode>::iterator iter;
    newNode Nodetemp;
    int mid=-1;
    string temp;
    for(int i=0;i<oldmesh.NumTetras;i++)
    {
        //get the six new points
        for(int j=0;j<6;j++)
        {
            Estart=edgeID[j][0];
            Eend=edgeID[j][1];
            startID=oldmesh.pTetras[i].vertices[Estart];
            endID=oldmesh.pTetras[i].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
                cout<<"Error in finding the new point!"<<endl;
                exit(1);
            }
            else
            {
                Nodetemp=iter->second;
                newPointID[j]=Nodetemp.localID;
            }
        }
        //construct new 8 tetras


        //No1
        newmesh.pTetras[count].vertices[0]=oldmesh.pTetras[i].vertices[0];
        newmesh.pTetras[count].vertices[1]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[2]=newPointID[1];//0-2
        newmesh.pTetras[count].vertices[3]=newPointID[2];//0-3
        count++;

        //No2
        newmesh.pTetras[count].vertices[0]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[1]=oldmesh.pTetras[i].vertices[1];
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;


        //No3
        newmesh.pTetras[count].vertices[0]=newPointID[1];//0-2
        newmesh.pTetras[count].vertices[1]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[2]=oldmesh.pTetras[i].vertices[2];
        newmesh.pTetras[count].vertices[3]=newPointID[5];//2-3
        count++;


        //No4
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[4];//1-3
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=oldmesh.pTetras[i].vertices[3];
        count++;


        //No5
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;


        //No6
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;

        //No7
        newmesh.pTetras[count].vertices[0]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[1]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[1];//0-2
        count++;


        //No8
        newmesh.pTetras[count].vertices[0]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[1]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=newPointID[1];//0-2
        count++;

    }
    assert(count==newmesh.NumTetras);
    //split the facets
    newNode midtemp;
    for(int i=0;i<oldmesh.NumTris;i++)
    {
        for(int j=0;j<3;j++)
        {
            startID=oldmesh.pTris[i].vertices[j%3];
            endID=oldmesh.pTris[i].vertices[(j+1)%3];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            string temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
                cout<<"Error in finding addedNodes!"<<endl;
            }
            else
            {
                midtemp=iter->second;
                oldmesh.pTris[i].addedNodes[(j+2)%3]=midtemp.localID;
            }
        }
    }
    newmesh.NumTris=oldmesh.NumTris*4;
    newmesh.pTris=new TRI[newmesh.NumTris];
    count=0;
    for(int i=0;i<oldmesh.NumTris;i++)
    {
        //NO1
        int patch_id;
//        if(oldmesh.pTris[i].iOppoProc!=-1)
        patch_id=ref.detachFace(i);
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[0];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[2];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[1];
//        if(oldmesh.pTris[i].iOppoProc!=-1)
        ref.attachFace(count+oldmesh.NumTris,patch_id);
        count++;
        //NO2
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[1];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[0];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[2];
//        if(oldmesh.pTris[i].iOppoProc!=-1)
        ref.attachFace(count+oldmesh.NumTris,patch_id);
        count++;
        //NO3
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[2];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[1];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[0];
//        if(oldmesh.pTris[i].iOppoProc!=-1)
        ref.attachFace(count+oldmesh.NumTris,patch_id);
        count++;
        //NO4
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].addedNodes[0];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[1];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[2];
//        if(oldmesh.pTris[i].iOppoProc!=-1)
        ref.attachFace(count+oldmesh.NumTris,patch_id);
        count++;
    }

    assert(count==newmesh.NumTris);

    return 1;

}
int sixNodesPattern_test(HYBRID_MESH &oldmesh, HYBRID_MESH &newmesh, map<string, newNode> &edgeHash,map<string, newNode> &faceHash)
{
    int edgeID_prism[9][2]=
    {
        0,1,//0
        0,2,//1
        0,3,//2
        1,2,//3
        1,4,//4
        2,5,//5
        3,4,//6
        3,5,//7
        4,5,//8
    };
    int prism_faceID[3][4]=
    {
      0,3,5,2,//0
      0,1,4,3,//1
      1,2,5,4 //2
    };

    int newPointID_prism[12];//corresponded with edgeID
    int Estart=-1;
    int Eend=-1;
    int startID=-1;
    int endID=-1;
    int count=0;

    PRISM prismsTemp;

    map<string,newNode>::iterator iter;
    newNode Nodetemp;
    int mid=-1;
    string temp;
    int vec[4];
    for(int i=0;i<oldmesh.NumPrsm;i++)
    {
        //get the eleven new points
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<4;++k){
                vec[k]=oldmesh.pPrisms[i].vertices[prism_faceID[j][k]];
            }
            sort(vec,vec+4);
            temp=IntToString(vec[0])+"_"+IntToString(vec[1])+"_"+IntToString(vec[2])+"_"+IntToString(vec[3]);
            iter=faceHash.find(temp);
            if(iter==faceHash.end())
            {
                cout<<"Error in finding the new point!"<<endl;
                exit(1);
            }
            else
            {
                Nodetemp=iter->second;
                newPointID_prism[j+9]=Nodetemp.localID;
            }

        }

        for(int j=0;j<9;j++)
        {
            Estart=edgeID_prism[j][0];
            Eend=edgeID_prism[j][1];
            startID=oldmesh.pPrisms[i].vertices[Estart];
            endID=oldmesh.pPrisms[i].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
                cout<<"Error in finding the new point!"<<endl;
                exit(1);
            }
            else
            {
                Nodetemp=iter->second;
                newPointID_prism[j]=Nodetemp.localID;
            }
        }
        for(int k=0;k<12;++k){
            if(newPointID_prism[k]<0||newPointID_prism[k]>9000000)
                cout<<newPointID_prism[k]<<endl;
        }

        //construct new 8 tetras
        //No1
        newmesh.pPrisms[count].vertices[0]=oldmesh.pPrisms[i].vertices[0];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[0];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[1];
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[2];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[10];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[9];
        count++;

        //No2
        newmesh.pPrisms[count].vertices[0]=oldmesh.pPrisms[i].vertices[1];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[3];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[0];
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[4];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[11];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[10];
        count++;

        //No3
        newmesh.pPrisms[count].vertices[0]=oldmesh.pPrisms[i].vertices[2];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[1];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[3];
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[5];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[9];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[11];
        count++;


        //No4
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[2];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[10];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[9];
        newmesh.pPrisms[count].vertices[0]=oldmesh.pPrisms[i].vertices[3];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[6];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[7];
        count++;


        //No5
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[4];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[11];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[10];
        newmesh.pPrisms[count].vertices[0]=oldmesh.pPrisms[i].vertices[4];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[8];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[6];
        count++;


        //No6
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[5];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[9];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[11];
        newmesh.pPrisms[count].vertices[0]=oldmesh.pPrisms[i].vertices[5];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[7];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[8];
        count++;

        //No7
        newmesh.pPrisms[count].vertices[0]=newPointID_prism[0];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[3];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[1];
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[10];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[11];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[9];
        count++;


        //No8
        newmesh.pPrisms[count].vertices[0]=newPointID_prism[6];
        newmesh.pPrisms[count].vertices[1]=newPointID_prism[8];
        newmesh.pPrisms[count].vertices[2]=newPointID_prism[7];
        newmesh.pPrisms[count].vertices[3]=newPointID_prism[10];
        newmesh.pPrisms[count].vertices[4]=newPointID_prism[11];
        newmesh.pPrisms[count].vertices[5]=newPointID_prism[9];
        count++;
    }
    assert(count==newmesh.NumPrsm);

    int edgeID[6][2]=
    {
        0,1,//0
        0,2,//1
        0,3,//2
        1,2,//3
        1,3,//4
        2,3//5
    };
    int newPointID[6];//corresponded with edgeID

     Estart=-1;
     Eend=-1;
     startID=-1;
     endID=-1;

     count=0;

    TETRAS tetrasTemp;

    mid=-1;
    for(int i=0;i<oldmesh.NumTetras;i++)
    {
        //get the six new points
        for(int j=0;j<6;j++)
        {
            Estart=edgeID[j][0];
            Eend=edgeID[j][1];
            startID=oldmesh.pTetras[i].vertices[Estart];
            endID=oldmesh.pTetras[i].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
                cout<<"Error in finding the new point!"<<endl;
                exit(1);
            }
            else
            {
                Nodetemp=iter->second;
                newPointID[j]=Nodetemp.localID;
            }
        }
        //construct new 8 tetras


        //No1
        newmesh.pTetras[count].vertices[0]=oldmesh.pTetras[i].vertices[0];
        newmesh.pTetras[count].vertices[1]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[2]=newPointID[1];//0-2
        newmesh.pTetras[count].vertices[3]=newPointID[2];//0-3
        count++;

        //No2
        newmesh.pTetras[count].vertices[0]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[1]=oldmesh.pTetras[i].vertices[1];
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;


        //No3
        newmesh.pTetras[count].vertices[0]=newPointID[1];//0-2
        newmesh.pTetras[count].vertices[1]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[2]=oldmesh.pTetras[i].vertices[2];
        newmesh.pTetras[count].vertices[3]=newPointID[5];//2-3
        count++;


        //No4
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[4];//1-3
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=oldmesh.pTetras[i].vertices[3];
        count++;


        //No5
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;


        //No6
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;

        //No7
        newmesh.pTetras[count].vertices[0]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[1]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[1];//0-2
        count++;


        //No8
        newmesh.pTetras[count].vertices[0]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[1]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=newPointID[1];//0-2
        count++;

    }
    assert(count==newmesh.NumTetras);
    //split the facets
    newNode midtemp;
    for(int i=0;i<oldmesh.NumTris;i++)
    {
        for(int j=0;j<3;j++)
        {
            startID=oldmesh.pTris[i].vertices[j%3];
            endID=oldmesh.pTris[i].vertices[(j+1)%3];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            string temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
         //       cout<<"Error in finding addedNodes!"<<endl;
            }
            else
            {
                midtemp=iter->second;
                oldmesh.pTris[i].addedNodes[(j+2)%3]=midtemp.localID;
            }
        }
    }
    newmesh.NumTris=oldmesh.NumTris*4;
    newmesh.pTris=new TRI[newmesh.NumTris];
    count=0;
    for(int i=0;i<oldmesh.NumTris;i++)
    {
        //NO1
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[0];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[2];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[1];
        count++;
        //NO2
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[1];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[0];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[2];

        count++;
        //NO3
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[2];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[1];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[0];
        count++;
        //NO4
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].addedNodes[0];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[1];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[2];
        count++;
    }

    assert(count==newmesh.NumTris);

    return 1;
}

int sixNodesPattern(HYBRID_MESH &oldmesh, HYBRID_MESH &newmesh, map<string, newNode> &edgeHash)
{
    int edgeID[6][2]=
    {
        0,1,//0
        0,2,//1
        0,3,//2
        1,2,//3
        1,3,//4
        2,3//5
    };
    int newPointID[6];//corresponded with edgeID

    int Estart=-1;
    int Eend=-1;
    int startID=-1;
    int endID=-1;

    int count=0;

    TETRAS tetrasTemp;

    map<string,newNode>::iterator iter;
    newNode Nodetemp;
    int mid=-1;
    string temp;
    for(int i=0;i<oldmesh.NumTetras;i++)
    {
        //get the six new points
        for(int j=0;j<6;j++)
        {
            Estart=edgeID[j][0];
            Eend=edgeID[j][1];
            startID=oldmesh.pTetras[i].vertices[Estart];
            endID=oldmesh.pTetras[i].vertices[Eend];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
                cout<<"Error in finding the new point!"<<endl;
                exit(1);
            }
            else
            {
                Nodetemp=iter->second;
                newPointID[j]=Nodetemp.localID;
            }
        }
        //construct new 8 tetras


        //No1
        newmesh.pTetras[count].vertices[0]=oldmesh.pTetras[i].vertices[0];
        newmesh.pTetras[count].vertices[1]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[2]=newPointID[1];//0-2
        newmesh.pTetras[count].vertices[3]=newPointID[2];//0-3
        count++;

        //No2
        newmesh.pTetras[count].vertices[0]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[1]=oldmesh.pTetras[i].vertices[1];
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;


        //No3
        newmesh.pTetras[count].vertices[0]=newPointID[1];//0-2
        newmesh.pTetras[count].vertices[1]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[2]=oldmesh.pTetras[i].vertices[2];
        newmesh.pTetras[count].vertices[3]=newPointID[5];//2-3
        count++;


        //No4
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[4];//1-3
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=oldmesh.pTetras[i].vertices[3];
        count++;


        //No5
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;


        //No6
        newmesh.pTetras[count].vertices[0]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[1]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=newPointID[4];//1-3
        count++;

        //No7
        newmesh.pTetras[count].vertices[0]=newPointID[0];//0-1
        newmesh.pTetras[count].vertices[1]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[2]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[3]=newPointID[1];//0-2
        count++;


        //No8
        newmesh.pTetras[count].vertices[0]=newPointID[3];//1-2
        newmesh.pTetras[count].vertices[1]=newPointID[2];//0-3
        newmesh.pTetras[count].vertices[2]=newPointID[5];//2-3
        newmesh.pTetras[count].vertices[3]=newPointID[1];//0-2
        count++;

    }
    assert(count==newmesh.NumTetras);
    //split the facets
    newNode midtemp;
    for(int i=0;i<oldmesh.NumTris;i++)
    {
        for(int j=0;j<3;j++)
        {
            startID=oldmesh.pTris[i].vertices[j%3];
            endID=oldmesh.pTris[i].vertices[(j+1)%3];
            if(startID>endID)
            {
                mid=startID;
                startID=endID;
                endID=mid;
            }
            string temp=IntToString(startID)+"_"+IntToString(endID);
            iter=edgeHash.find(temp);
            if(iter==edgeHash.end())
            {
         //       cout<<"Error in finding addedNodes!"<<endl;
            }
            else
            {
                midtemp=iter->second;
                oldmesh.pTris[i].addedNodes[(j+2)%3]=midtemp.localID;
            }
        }
    }
    newmesh.NumTris=oldmesh.NumTris*4;
    newmesh.pTris=new TRI[newmesh.NumTris];
    count=0;
    for(int i=0;i<oldmesh.NumTris;i++)
    {
        //NO1
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[0];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[2];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[1];
        count++;
        //NO2
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[1];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[0];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[2];

        count++;
        //NO3
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].vertices[2];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[1];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[0];
        count++;
        //NO4
        newmesh.pTris[count]=oldmesh.pTris[i];
        newmesh.pTris[count].vertices[0]=oldmesh.pTris[i].addedNodes[0];
        newmesh.pTris[count].vertices[1]=oldmesh.pTris[i].addedNodes[1];
        newmesh.pTris[count].vertices[2]=oldmesh.pTris[i].addedNodes[2];
        count++;
    }

    assert(count==newmesh.NumTris);

    return 1;

}

int constructFacets_test(HYBRID_MESH&mesh, HYBRID_MESH& globalMesh,map<string,int64_t>&tri_globalID)
{

    int originNum=mesh.NumTris;
    for(int i=0;i<mesh.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            mesh.pTetras[i].neighbors[j]=-1;
        }
    }
    for(int i=0;i<mesh.NumPrsm;i++)
    {
        for(int j=0;j<5;j++)
        {
            mesh.pPrisms[i].neighbors[j]=-1;
        }
    }
    cout<<"Tetras Number is: "<<mesh.NumTetras<<endl;
    cout<<"Prism Number is: "<<mesh.NumPrsm<<endl;
    setupCellNeig_test(mesh.NumNodes,mesh.NumTetras,&mesh);
    map<string,bool>triMap;

    int vertices[3];
    for(int i=0;i<mesh.NumTris;i++)
    {
        vertices[0]=mesh.pTris[i].vertices[0];
        vertices[1]=mesh.pTris[i].vertices[1];
        vertices[2]=mesh.pTris[i].vertices[2];
        sort(vertices,vertices+3);
        string temp=to_string(vertices[0])+"_";
        temp+=to_string(vertices[1]);
        temp+="_";
        temp+=to_string(vertices[2]);
//        for(int j=0;j<3;++j)
//          cout<<mesh.pTris[i].vertices[j]<<"   ";
//        cout<<endl;
        triMap[temp]=true;
    }



    int triVID[4][3]={
        1,3,2 ,
        0,2,3,
        0,3,1,
        0,1,2
    };


    int* pris[5];
    int pris1[3] = {0, 1, 2};
    int pris2[3] = {3, 4, 5};
    int pris3[4] = {0, 1, 4, 3};
    int pris4[4] = {1, 2, 5, 4};
    int pris5[4] = {0, 2, 5, 3};
    pris[0] = pris1;
    pris[1] = pris2;
    pris[2] = pris3;
    pris[3] = pris4;
    pris[4] = pris5;

//    ofstream a;
//    a.open("triange.txt");
//    for(int i=0;i<mesh.NumTetras;++i){
//        for(int j=0;j<4;++j)
//        a<<mesh.pTetras[i].neighbors[j]<<"   ";
//        a<<endl;
//    }
//    a.close();

             ofstream a;
             a.open("vtk.vtk");
             a<<"# vtk DataFile Version 2.0"<<endl;
             a<<"boundary layer mesh"<<endl;
             a<<"ASCII"<<endl;
             a<<"DATASET UNSTRUCTURED_GRID"<<endl;
             a<<"POINTS "<<mesh.NumNodes<<" float"<<endl;
             for(int i=0;i<mesh.NumNodes;++i){
                 a<<mesh.nodes[i].coord.x<<" "<<mesh.nodes[i].coord.y<<" "<<mesh.nodes[i].coord.z<<endl;
             }
             a<<"CELLS"<<" "<<mesh.numOfCells()<<" "<<mesh.NumPrsm*7+mesh.NumTetras*5<<endl;

             for(int i=0;i<mesh.NumPrsm;++i){
                 a<<6;
                 for(int j=0;j<6;++j){
                     a<<" "<<mesh.pPrisms[i].vertices[j];
                 }
                 a<<endl;
             }
             for(int i=0;i<mesh.NumTetras;++i){
                 a<<4;
                 for(int j=0;j<4;++j){
                     a<<" "<<mesh.pTetras[i].vertices[j];
                 }
                 a<<endl;
             }
             a<<"CELL_TYPES "<<mesh.numOfCells()<<endl;
             for(int i=0;i<mesh.NumPrsm;++i){
                 a<<13<<endl;
             }
             for(int i=0;i<mesh.NumTetras;++i){
                 a<<10<<endl;
             }
             a.close();

             a.open("table.txt");
             for(int i=0;i<mesh.NumPrsm;++i){
                 a<<mesh.pPrisms[i].index<<" ";
                 for(int j=0;j<5;++j)
                 a<<mesh.pPrisms[i].neighbors[j]<<"  ";
                 a<<"vertices  ";
                 for(int j=0;j<6;++j)
                 a<<mesh.pPrisms[i].vertices[j]<<"  ";
                 a<<endl;
             }
             for(int i=0;i<mesh.NumTetras;++i){
                 a<<mesh.pTetras[i].index<<" ";
                 for(int j=0;j<4;++j)
                 a<<mesh.pTetras[i].neighbors[j]<<"  ";
                 a<<"vertices  ";
                 for(int j=0;j<4;++j)
                 a<<mesh.pTetras[i].vertices[j]<<"  ";
                 a<<endl;
             }
             a.close();

    vector<TRI> interFacets;
    int count=0;
    int Numtri=0;
    for(int i=0;i<mesh.NumTetras;i++)
    {
        int forth=-1;
        for(int j=0;j<4;j++)
        {
            if(mesh.pTetras[i].neighbors[j]==-1)
            {

                Numtri++;
                forth=j;
                count=0;
                for(int k=0;k<3;k++)
                {
                    vertices[k]=mesh.pTetras[i].vertices[triVID[j][k]];
                }

                /*
                for(int k=0;k<4;k++)
                {
                    if(k!=forth)
                    {
                        vertices[count]=mesh.pTetras[i].vertices[k];
                        count++;

                    }
                }*/
                TRI triangle;
                constructOneTriangle(vertices,triangle);
                triangle.iCell=i+mesh.NumPrsm;
                int globalID=mesh.pTetras[i].index;


                if(globalID<globalMesh.NumPrsm||globalID>=globalMesh.NumPrsm+globalMesh.NumTetras)
                {
//                    cout<<i<<endl;
//                      cout<<globalID<<endl;
                    cout<<"Error in finding the gloabl index!"<<endl;
                }
                globalID=globalMesh.pTetras[globalID-globalMesh.NumPrsm].neighbors[forth];
       //        cout<<mesh.pTetras[i].neighbors[forth]<<endl;
                if(globalID==-1)
                {
                    triangle.iOppoProc=-1;
                //      cout<<triangle.iOppoProc<<endl;

                }
                else if(globalID<globalMesh.NumPrsm)
                {
                    triangle.iOppoProc=globalMesh.pPrisms[globalID].partMarker;
                    //cout<<"forth:  "<<temp<<endl;
                   // cout<<triangle.iOppoProc<<endl;

                }
                else
                {
                    triangle.iOppoProc=globalMesh.pTetras[globalID-globalMesh.NumPrsm].partMarker;
                    //cout<<"forth:  "<<temp<<endl;
                   // cout<<triangle.iOppoProc<<endl;

                }
                interFacets.push_back(triangle);
            }
            else
            {
                //nothing
            }
        }
    }


    for(int i=0;i<mesh.NumPrsm;i++)
    {
        int forth=-1;
        for(int j=0;j<2;j++)
        {
            if(mesh.pPrisms[i].neighbors[j]==-1)
            {

                Numtri++;
                forth=j;
                count=0;
                for(int k=0;k<3;k++)
                {
                    vertices[k]=mesh.pPrisms[i].vertices[pris[j][k]];
                }

                /*
                for(int k=0;k<4;k++)
                {
                    if(k!=forth)
                    {
                        vertices[count]=mesh.pTetras[i].vertices[k];
                        count++;

                    }
                }*/
                TRI triangle;
                constructOneTriangle(vertices,triangle);
                triangle.iCell=i;
                int globalID=mesh.pPrisms[i].index;


                if(globalID>globalMesh.NumPrsm)
                {
//                    cout<<i<<endl;
//                      cout<<globalID<<endl;
                    cout<<"Error in finding the gloabl index!"<<endl;
                }
                globalID=globalMesh.pPrisms[globalID].neighbors[forth];
       //        cout<<mesh.pTetras[i].neighbors[forth]<<endl;
                if(globalID==-1)
                {
                    triangle.iOppoProc=-1;
                //      cout<<triangle.iOppoProc<<endl;

                }
                else if(globalID<globalMesh.NumPrsm)
                {
                    triangle.iOppoProc=globalMesh.pPrisms[globalID].partMarker;
                    //cout<<"forth:  "<<temp<<endl;
                   // cout<<triangle.iOppoProc<<endl;

                }
                else
                {
                    triangle.iOppoProc=globalMesh.pTetras[globalID-globalMesh.NumPrsm].partMarker;
                    //cout<<"forth:  "<<temp<<endl;
                   // cout<<triangle.iOppoProc<<endl;

                }
                interFacets.push_back(triangle);
            }
            else
            {
                //nothing
            }
        }
    }





    TRI * facets= new TRI[mesh.NumTris];

    for(int i=0;i<mesh.NumTris;i++)
    {
        facets[i]=mesh.pTris[i];
    }
    count=mesh.NumTris;
  //  cout<<interFacets.size()<<"!!!!!!!!!!!!!!!!!!!!"<<count<<endl;    //butong  xuyao xiugai

    mesh.NumTris=interFacets.size();

    if(mesh.pTris!=nullptr)
    {
        delete []mesh.pTris;
        mesh.pTris=nullptr;
    }
    mesh.pTris=new TRI[mesh.NumTris];
    for(int i=0;i<count;i++)
    {
        mesh.pTris[i]=facets[i];
    }
    if(facets!=nullptr)
    {
        delete []facets;
        facets=nullptr;
    }

    //   mesh.pTris=facets;
    //facets=nullptr;
 //   cout<<"No error ..........................."<<endl;

    //  count=;
    int  sum1=0;
    map<string,bool>triMap2;
    for(int i=0;i<interFacets.size();i++){
            vertices[0]=interFacets[i].vertices[0];
            vertices[1]=interFacets[i].vertices[1];
            vertices[2]=interFacets[i].vertices[2];
            sort(vertices,vertices+3);
            string temp=to_string(vertices[0])+"_";
            temp+=to_string(vertices[1]);
            temp+="_";
            temp+=to_string(vertices[2]);
            triMap2[temp]=true;

    }
    for(auto i: triMap){
        string temp=i.first;
        if(triMap2.find(temp)!=triMap2.end()){
            ;
        }
        else
            ;
     //       cout<<temp<<endl;
    }
    for(int i=0;i<interFacets.size();i++)
    {

        vertices[0]=interFacets[i].vertices[0];
        vertices[1]=interFacets[i].vertices[1];
        vertices[2]=interFacets[i].vertices[2];
        sort(vertices,vertices+3);
        string temp=to_string(vertices[0])+"_";
        temp+=to_string(vertices[1]);
        temp+="_";
        temp+=to_string(vertices[2]);

        if(triMap.find(temp)!=triMap.end())
        {
            sum1++;
            //case : surface facet or interface created in the last partition

        }
        else
        {
            //case : new interface
                            //cout<<"test"<<endl;
            int globalVertices[3];
            globalVertices[0]=mesh.nodes[vertices[0]].index;
            globalVertices[1]=mesh.nodes[vertices[1]].index;
            globalVertices[2]=mesh.nodes[vertices[2]].index;
            sort(globalVertices,globalVertices+3);
            string gTemp=to_string(globalVertices[0])+"_"+to_string(globalVertices[1])+"_"+to_string(globalVertices[2]);

            tri_globalID[gTemp]=-1;//temp value == -1////need gloabl index;

            mesh.pTris[count]=interFacets[i];
            count++;
        }
    }
#ifdef DEBUG
    cout<<interFacets.size()<<"!!!!!!!!!!!!!!!!!!!!"<<count<<endl;    //butong  xuyao xiugai
    cout<<"originNum =="<<originNum<<endl;
    cout<<"sum1=="<<sum1<<endl;
    cout<<"Counter of the facets is: "<<Numtri<<endl;
    cout<<"Surface triangle Number is "<<originNum<<endl;
    cout<<"All triangles Number is:"<<interFacets.size()<<endl;
#endif
//    assert(count==mesh.NumTris);

    cout<<"Counter of the facets is: "<<Numtri<<endl;
    cout<<"Surface triangle Number is "<<originNum<<endl;
    cout<<"All triangles Number is:"<<interFacets.size()<<endl;



    return 1;
}


int constructFacets(HYBRID_MESH&mesh, HYBRID_MESH& globalMesh,map<string,int64_t>&tri_globalID)
{

    int originNum=mesh.NumTris;
    for(int i=0;i<mesh.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            mesh.pTetras[i].neighbors[j]=-1;
        }
    }
    cout<<"Tetras Number is: "<<mesh.NumTetras<<endl;

    setupCellNeig(mesh.NumNodes,mesh.NumTetras,mesh.pTetras);

    map<string,bool>triMap;

    int vertices[3];
    for(int i=0;i<mesh.NumTris;i++)
    {
        vertices[0]=mesh.pTris[i].vertices[0];
        vertices[1]=mesh.pTris[i].vertices[1];
        vertices[2]=mesh.pTris[i].vertices[2];
        sort(vertices,vertices+3);
        string temp=IntToString(vertices[0])+"_";
        temp+=IntToString(vertices[1]);
        temp+="_";
        temp+=IntToString(vertices[2]);
        triMap[temp]=true;
    }



    int triVID[4][3]={
        1,3,2 ,
        0,2,3,
        0,3,1,
        0,1,2
    };


    vector<TRI> interFacets;
    int count=0;
    int Numtri=0;
    for(int i=0;i<mesh.NumTetras;i++)
    {
        int forth=-1;
        for(int j=0;j<4;j++)
        {
            if(mesh.pTetras[i].neighbors[j]==-1)
            {

                Numtri++;
                forth=j;
                count=0;
                for(int k=0;k<3;k++)
                {
                    vertices[k]=mesh.pTetras[i].vertices[triVID[j][k]];
                }

                /*
                for(int k=0;k<4;k++)
                {
                    if(k!=forth)
                    {
                        vertices[count]=mesh.pTetras[i].vertices[k];
                        count++;

                    }
                }*/
                TRI triangle;
                constructOneTriangle(vertices,triangle);
                triangle.iCell=i;
                int globalID=mesh.pTetras[i].index;///////////////////
               // cout<<globalID<<endl;
                if(globalID<0||globalID>=globalMesh.NumTetras)
                {
                    cout<<"Error in finding the gloabl index!"<<endl;
                }
                globalID=globalMesh.pTetras[globalID].neighbors[forth];
       //        cout<<mesh.pTetras[i].neighbors[forth]<<endl;
                if(globalID==-1)
                {
                    triangle.iOppoProc=-1;
                //      cout<<triangle.iOppoProc<<endl;

                }
                else
                {
                    triangle.iOppoProc=globalMesh.pTetras[globalID].partMarker;
                    //cout<<"forth:  "<<temp<<endl;
                   // cout<<triangle.iOppoProc<<endl;

                }
                interFacets.push_back(triangle);
            }
            else
            {
                //nothing
            }
        }
    }

   // cout<<"No error ..........................."<<endl;

    TRI * facets= new TRI[mesh.NumTris];

    for(int i=0;i<mesh.NumTris;i++)
    {
        facets[i]=mesh.pTris[i];
    }
    count=mesh.NumTris;

    mesh.NumTris=interFacets.size();

    if(mesh.pTris!=nullptr)
    {
        delete []mesh.pTris;
        mesh.pTris=nullptr;
    }
    mesh.pTris=new TRI[mesh.NumTris];
    for(int i=0;i<count;i++)
    {
        mesh.pTris[i]=facets[i];
    }
    if(facets!=nullptr)
    {
        delete []facets;
        facets=nullptr;
    }

    //   mesh.pTris=facets;
    //facets=nullptr;

    //  count=;
    for(int i=0;i<interFacets.size();i++)
    {
        vertices[0]=interFacets[i].vertices[0];
        vertices[1]=interFacets[i].vertices[1];
        vertices[2]=interFacets[i].vertices[2];
        sort(vertices,vertices+3);
        string temp=IntToString(vertices[0])+"_";
        temp+=IntToString(vertices[1]);
        temp+="_";
        temp+=IntToString(vertices[2]);
        if(triMap.find(temp)!=triMap.end())
        {
            //case : surface facet or interface created in the last partition


        }
        else
        {
            //case : new interface
                            //cout<<"test"<<endl;
            int globalVertices[3];
            globalVertices[0]=mesh.nodes[vertices[0]].index;
            globalVertices[1]=mesh.nodes[vertices[1]].index;
            globalVertices[2]=mesh.nodes[vertices[2]].index;
            sort(globalVertices,globalVertices+3);
            string gTemp=IntToString(globalVertices[0])+"_"+IntToString(globalVertices[1])+"_"+IntToString(globalVertices[2]);

            tri_globalID[gTemp]=-1;//temp value == -1////need gloabl index;

            mesh.pTris[count]=interFacets[i];
            count++;
        }
    }

    assert(count==mesh.NumTris);

    cout<<"Counter of the facets is: "<<Numtri<<endl;
    cout<<"Surface triangle Number is "<<originNum<<endl;
    cout<<"All triangles Number is:"<<interFacets.size()<<endl;



    return 1;
}

int constructOneTriangle(int *vertices,TRI& triangle)
{
    //TRIANGLE triangle;

    triangle.vertices[0]=vertices[0];
    triangle.vertices[1]=vertices[1];
    triangle.vertices[2]=vertices[2];

    // triangle.

    return 1;
}




int managerProc(int rank, int nparts, HYBRID_MESH*tetrasPart, HYBRID_MESH& construcTetras,Refletion ref ,Refletion construcref, int &refineSize,vector<string>&bcstring, MPI_Comm comm)
{
    //prepare nodes data to send

    int count=0;
    int *sendcounts=new int[nparts];

    int *displs=new int[nparts];
    int offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumNodes*3;
        displs[i]=offset;
        offset+=tetrasPart[i].NumNodes*3;//need to check
        count=count+tetrasPart[i].NumNodes;
    }


    int tablesize = ref.subject_table.size();
    int *sendtable=new int[122];
    for(int i=0;i<122;++i){
        sendtable[i]=ref.subject_table[i];
    }
    MPI_Bcast(&tablesize,1,MPI_INT,0,comm);
    MPI_Bcast(sendtable,tablesize,MPI_INT,0,comm);

    double *coordata=new double[count*3];

    cout<<"count size: "<<count*3<<endl;

    int countwo=0;

    // int *nodeGlobabIndex=nullptr;

    int *nodeGlobabIndex=new int[count];

    int *procSize=new int[count];

    int count3=0;

    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            coordata[countwo]=tetrasPart[i].nodes[j].coord.x;
            countwo++;
            coordata[countwo]=tetrasPart[i].nodes[j].coord.y;
            countwo++;
            coordata[countwo]=tetrasPart[i].nodes[j].coord.z;
            countwo++;

            nodeGlobabIndex[count3]=tetrasPart[i].nodes[j].index;

            procSize[count3]=tetrasPart[i].nodes[j].procs.size();

            count3++;

        }
    }
    assert(count*3==countwo);

    assert(count3==count);


    int root=0;

    int *countbuf=new int[nparts];

    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumNodes*3;
    }

    int recvNum=0;

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the node number



    /*
    int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
            */

    cout<<"RecvNum: "<<recvNum<<endl;

    int recvcount=recvNum;

    double *recvbuf=new double[recvcount];

    MPI_Scatterv(coordata,sendcounts,displs,MPI_DOUBLE,recvbuf,recvcount,MPI_DOUBLE,rank,comm);//send the coord ;

    if(coordata!=nullptr)
        delete []coordata;

    coordata=nullptr;

    construcTetras.NumNodes=recvNum/3;
    construcTetras.nodes=new Node[construcTetras.NumNodes];

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].coord.x=recvbuf[i*3+0];
        construcTetras.nodes[i].coord.y=recvbuf[i*3+1];
        construcTetras.nodes[i].coord.z=recvbuf[i*3+2];
    }
    if(recvbuf!=nullptr)
        delete []recvbuf;
    recvbuf=nullptr;

    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumNodes;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumNodes;
        // count=count+tetrasPart[i].NumNodes;
    }

    int NumNodes=recvcount/3;

    int *nodeGlobalID=new int[NumNodes];
    int *recvProcSize=new int[NumNodes];


    //for(int i=0;i<)

    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,nodeGlobalID,NumNodes,MPI_INT,rank,comm);//send the global index of the node ;

    MPI_Scatterv(procSize,sendcounts,displs,MPI_INT,recvProcSize,NumNodes,MPI_INT,rank,comm);//send the procs size of every node ;


    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].NumProcs=recvProcSize[i];
    }

    if(procSize!=nullptr)
    {
        delete []procSize;
        procSize=nullptr;
    }

    int allSetSize=0;
    int *sumProcs=new int[nparts];
    for(int i=0;i<nparts;i++)
    {
        sumProcs[i]=0;
    }

    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            sumProcs[i]+=tetrasPart[i].nodes[j].procs.size();
           // allSetSize+=tetrasPart[i].nodes[j].procs.size();
        }
        allSetSize+=sumProcs[i];
    }
    int setSizePart=-1;
    MPI_Scatter(sumProcs,1,MPI_INT,&setSizePart,1,MPI_INT,rank,comm);//send the size of procs of every part

    int *allSetValue=new int[allSetSize];
    set<int>::iterator setIter;
    count=0;
    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            for(setIter=tetrasPart[i].nodes[j].procs.begin();setIter!=tetrasPart[i].nodes[j].procs.end();setIter++)
            {
                int setValue=*setIter;
                allSetValue[count]=setValue;
                count++;

            }
        }
    }

    offset=0;
    for(int i=0;i<nparts;i++)
    {
        displs[i]=offset;//need to check
        offset+=sumProcs[i];
    }
    int *recvSetValue=new int[setSizePart];
    MPI_Scatterv(allSetValue,sumProcs,displs,MPI_INT,recvSetValue,setSizePart,MPI_INT,rank,comm);//send the value of the procs



    int allCount=0;
    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        int NumProcs=construcTetras.nodes[i].NumProcs;
       // count=0;
        for(int j=0;j<NumProcs;j++)
        {
            construcTetras.nodes[i].procs.insert(recvSetValue[allCount]);
            allCount++;
        }
    }

    assert(allCount==setSizePart);

    if(recvSetValue!=nullptr)
        delete []recvSetValue;
    recvSetValue=nullptr;


    if(nodeGlobabIndex!=nullptr)
        delete []nodeGlobabIndex;
    nodeGlobabIndex=nullptr;

    if(procSize!=nullptr)
        delete []procSize;
    procSize=nullptr;

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].index=nodeGlobalID[i];
       // construcTetras.nodes[i].procs.
    }

    if(nodeGlobalID!=nullptr)
        delete []nodeGlobalID;

    nodeGlobalID=nullptr;

    int sumTetras=0;
    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumTetras*4;
        sumTetras+=tetrasPart[i].NumTetras*4;
    }

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the tetras number

    construcTetras.NumTetras=recvNum/4;

    cout<<"TetrasNum: "<<recvNum/4<<endl;

    construcTetras.pTetras=new TETRAS[construcTetras.NumTetras];

    int *tetrasVertices=new int[sumTetras];

    count=0;

    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumTetras;j++)
        {
            tetrasVertices[count*4+0]=tetrasPart[i].pTetras[j].vertices[0];
            tetrasVertices[count*4+1]=tetrasPart[i].pTetras[j].vertices[1];
            tetrasVertices[count*4+2]=tetrasPart[i].pTetras[j].vertices[2];
            tetrasVertices[count*4+3]=tetrasPart[i].pTetras[j].vertices[3];
            count++;

        }
    }


    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumTetras*4;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumTetras*4;
        // count=count+tetrasPart[i].NumNodes;
    }

    int *vTetras=new int[recvNum];

    MPI_Scatterv(tetrasVertices,sendcounts,displs,MPI_INT,vTetras,recvNum,MPI_INT,rank,comm);//send the tetras  ;


    if(tetrasVertices!=nullptr)
        delete []tetrasVertices;
    tetrasVertices=nullptr;

    for(int i=0;i<construcTetras.NumTetras;i++)
    {
        construcTetras.pTetras[i].partMarker=rank;
        construcTetras.pTetras[i].vertices[0]=vTetras[i*4+0];
        construcTetras.pTetras[i].vertices[1]=vTetras[i*4+1];
        construcTetras.pTetras[i].vertices[2]=vTetras[i*4+2];
        construcTetras.pTetras[i].vertices[3]=vTetras[i*4+3];

    }
    if(vTetras!=nullptr)
        delete []vTetras;
    vTetras=nullptr;

    int sumTris=0;
    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumTris*3;
        sumTris+=tetrasPart[i].NumTris*3;
    }

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the tris number

    construcTetras.NumTris=recvNum/3;

    cout<<"TrisNum: "<<recvNum/3<<endl;

    construcTetras.pTris=new TRI[construcTetras.NumTris];

    int *triVertices=new int[sumTris];

    count=0;
    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumTris;j++)
        {
            triVertices[count*3+0]=tetrasPart[i].pTris[j].vertices[0];
            triVertices[count*3+1]=tetrasPart[i].pTris[j].vertices[1];
            triVertices[count*3+2]=tetrasPart[i].pTris[j].vertices[2];
            count++;
        }
    }
    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumTris*3;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumTris*3;
        // count=count+tetrasPart[i].NumNodes;
    }

    int *vTris=new int[recvNum];

    MPI_Scatterv(triVertices,sendcounts,displs,MPI_INT,vTris,recvNum,MPI_INT,rank,comm);//send the tris  ;

    if(triVertices!=nullptr)
        delete []triVertices;
    triVertices=nullptr;

    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].partMarker=rank;
        construcTetras.pTris[i].vertices[0]=vTris[i*3+0];
        construcTetras.pTris[i].vertices[1]=vTris[i*3+1];
        construcTetras.pTris[i].vertices[2]=vTris[i*3+2];

    }

    if(vTris!=nullptr)
        delete []vTris;
    vTris=nullptr;
    count=0;

    sumTris=0;
    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumTris*3;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumTris*3;
        sumTris+=tetrasPart[i].NumTris*3;
        // count=count+tetrasPart[i].NumNodes;
    }


    int *attributs=new int[sumTris];// the first value is iSurf , the second value is iOppoProc
    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumTris;j++)
        {
            attributs[count*3+0]=tetrasPart[i].pTris[j].iSurf;
            attributs[count*3+1]=tetrasPart[i].pTris[j].iOppoProc;
            attributs[count*3+2]=tetrasPart[i].pTris[j].iSurface;
            count++;
        }
    }
    assert(count==sumTris/3);

    cout<<"Triangle Number is: "<<sumTris/3<<endl;

    int *recvAttr=new int[3*construcTetras.NumTris];

    recvNum=construcTetras.NumTris*3;

    MPI_Scatterv(attributs,sendcounts,displs,MPI_INT,recvAttr,recvNum,MPI_INT,rank,comm);//send the tris  ;//send the tris attributes


    if(attributs!=nullptr)
        delete []attributs;
    attributs=nullptr;
    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].iSurf=recvAttr[i*3+0];
        construcTetras.pTris[i].iOppoProc=recvAttr[i*3+1];
        construcTetras.pTris[i].iSurface=recvAttr[i*3+2];

    }

    if(recvAttr!=nullptr)
        delete []recvAttr;
    recvAttr=nullptr;


    int *refines=new int[nparts];

    for(int i=0;i<nparts;i++)
    {
        refines[i]=refineSize;
    }
    int recvRefine=-1;
    MPI_Scatter(refines,1,MPI_INT,&recvRefine,1,MPI_INT,rank,comm);

    if(refines!=nullptr)
        delete []refines;
    refines=nullptr;

    stringstream ss;
    for (int i = 0; i < bcstring.size(); ++i)
        ss << bcstring[i] << "\n";
    const string& bufStr = ss.str();
    cout << bufStr << endl;
    int stringSize=bufStr.size();

    //int *strSize=new int[];

   // MPI_Scatter()
    MPI_Bcast(&stringSize,1,MPI_INT,rank,comm);
    MPI_Bcast(const_cast<char*>(bufStr.c_str()), stringSize, MPI_CHAR,rank,comm);

    cout<<"MPI starting..."<<endl;
    return 1;

}

int workerProc(int rank, HYBRID_MESH&construcTetras,Refletion construcref, int &refineSize, vector<string> &bcstring, MPI_Comm comm)
{

    int tablesize=0;
   MPI_Bcast(&tablesize,1,MPI_INT,0,comm);
       cout<<tablesize<<endl;
    int *sendtable=new int[tablesize];
     MPI_Bcast(sendtable,tablesize,MPI_INT,0,comm);
     for(int i=0;i<tablesize;++i){
         cout<<sendtable[i]<<endl;
     }
    int countbuf[2];
    int recvNum=0;
    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);

    cout<<"RecvNum: "<<recvNum<<endl;

    double *recvbuf=new double[recvNum];

    double coordata[2];

    int sendcounts[2];

    int displs[2];

    MPI_Scatterv(coordata,sendcounts,displs,MPI_DOUBLE,recvbuf,recvNum,MPI_DOUBLE,0,comm);//receive the coord ;

    construcTetras.NumNodes=recvNum/3;
    construcTetras.nodes=new Node[construcTetras.NumNodes];

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].localID=i;
        construcTetras.nodes[i].coord.x=recvbuf[i*3+0];
        construcTetras.nodes[i].coord.y=recvbuf[i*3+1];
        construcTetras.nodes[i].coord.z=recvbuf[i*3+2];
    }

    if(recvbuf!=nullptr)
    {
        delete []recvbuf;
        recvbuf=nullptr;
    }

    int NumNodes=recvNum/3;



    int *nodeGlobalID=new int[NumNodes];

    cout<<"Numnodes: "<<NumNodes<<endl;

    int *nodeGlobabIndex=new int[2];

    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,nodeGlobalID,NumNodes,MPI_INT,0,comm);//receive the global index of the node ;

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].index=nodeGlobalID[i];
    }

    if(nodeGlobalID!=nullptr)
    {
        delete []nodeGlobalID;
        nodeGlobalID=nullptr;
    }

    int *recvProcSize=new int[NumNodes];
    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,recvProcSize,NumNodes,MPI_INT,0,comm);//receive the procs size of every node ;


    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].NumProcs=recvProcSize[i];
    }

    if(recvProcSize!=nullptr)
        delete []recvProcSize;
    recvProcSize=nullptr;

    int setSizePart=-1;
    MPI_Scatter(countbuf,1,MPI_INT,&setSizePart,1,MPI_INT,0,comm);//receive the size of procs of every part

    int *recvSetValue=new int[setSizePart];
    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,recvSetValue,setSizePart,MPI_INT,0,comm);//receive the value of the procs



    int allCount=0;
    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        int NumProcs=construcTetras.nodes[i].NumProcs;
 //       count=0;
        for(int j=0;j<NumProcs;j++)
        {
            construcTetras.nodes[i].procs.insert(recvSetValue[allCount]);
            allCount++;
        }
    }

    assert(allCount==setSizePart);

    if(recvSetValue!=nullptr)
        delete []recvSetValue;
    recvSetValue=nullptr;


    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);//receive the tetras number

    construcTetras.NumTetras=recvNum/4;

    cout<<"TetrasNum: "<<recvNum/4<<endl;

    construcTetras.pTetras=new TETRAS[construcTetras.NumTetras];

    int *vTetras=new int[recvNum];

    int *tetrasVertices=new int[2];

    MPI_Scatterv(tetrasVertices,sendcounts,displs,MPI_INT,vTetras,recvNum,MPI_INT,0,comm);//receive the tetras  ;


    for(int i=0;i<construcTetras.NumTetras;i++)
    {
        construcTetras.pTetras[i].partMarker=rank;
        construcTetras.pTetras[i].localID=i;
        construcTetras.pTetras[i].vertices[0]=vTetras[i*4+0];
        construcTetras.pTetras[i].vertices[1]=vTetras[i*4+1];
        construcTetras.pTetras[i].vertices[2]=vTetras[i*4+2];
        construcTetras.pTetras[i].vertices[3]=vTetras[i*4+3];

    }

    if(vTetras!=nullptr)
    {
        delete []vTetras;
        vTetras=nullptr;
    }



    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);//recv the tris number

    construcTetras.NumTris=recvNum/3;

    cout<<"TrisNum: "<<recvNum/3<<endl;

    construcTetras.pTris=new TRI[construcTetras.NumTris];

    int triVertices[2];

    int *vTris=new int[construcTetras.NumTris*3];

  //  cout<<"worker zhvliu..............................."<<endl;

    MPI_Scatterv(triVertices,sendcounts,displs,MPI_INT,vTris,recvNum,MPI_INT,0,comm);//recv the tris  ;

    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].partMarker=rank;
        construcTetras.pTris[i].vertices[0]=vTris[i*3+0];
        construcTetras.pTris[i].vertices[1]=vTris[i*3+1];
        construcTetras.pTris[i].vertices[2]=vTris[i*3+2];

    }

    if(vTris!=nullptr)
    {
        delete []vTris;
        vTris=nullptr;
    }

    int attributs[3];


    recvNum=construcTetras.NumTris*3;

    int *recvAttr=new int[recvNum];



    MPI_Scatterv(attributs,sendcounts,displs,MPI_INT,recvAttr,recvNum,MPI_INT,0,comm);//recv send the tris attributes


    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].iSurf=recvAttr[i*3+0];
        construcTetras.pTris[i].iOppoProc=recvAttr[i*3+1];
        construcTetras.pTris[i].iSurface=recvAttr[i*3+2];

    }


    MPI_Scatter(attributs,1,MPI_INT,&refineSize,1,MPI_INT,0,comm);

    if(recvAttr!=nullptr)
    {
        delete []recvAttr;
        recvAttr=nullptr;
    }

    // cout << bufStr << endl;
    int stringSize;
    MPI_Bcast(&stringSize,1,MPI_INT,0,comm);

    cout<<"StringSize == "<<stringSize<<endl;

    char *buf = new char[stringSize + 1];
    MPI_Bcast(buf, stringSize, MPI_CHAR,0,comm);
    buf[stringSize] = '\0';

    stringstream ss;
    ss << string(buf);
    string bufStr;
    while (getline(ss, bufStr))
        bcstring.push_back(bufStr);

    // cout << string(buf) << "Received!\n";



    delete [] buf;

    return 1;
}
int managerProc(int rank, int nparts, HYBRID_MESH*tetrasPart, HYBRID_MESH& construcTetras, int &refineSize,vector<string>&bcstring, MPI_Comm comm)
{
    //prepare nodes data to send

    int count=0;
    int *sendcounts=new int[nparts];

    int *displs=new int[nparts];

    int offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumNodes*3;
        displs[i]=offset;
        offset+=tetrasPart[i].NumNodes*3;//need to check
        count=count+tetrasPart[i].NumNodes;
    }

    double *coordata=new double[count*3];

    cout<<"count size: "<<count*3<<endl;

    int countwo=0;

    // int *nodeGlobabIndex=nullptr;

    int *nodeGlobabIndex=new int[count];

    int *procSize=new int[count];

    int count3=0;

    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            coordata[countwo]=tetrasPart[i].nodes[j].coord.x;
            countwo++;
            coordata[countwo]=tetrasPart[i].nodes[j].coord.y;
            countwo++;
            coordata[countwo]=tetrasPart[i].nodes[j].coord.z;
            countwo++;

            nodeGlobabIndex[count3]=tetrasPart[i].nodes[j].index;

            procSize[count3]=tetrasPart[i].nodes[j].procs.size();

            count3++;

        }
    }
    assert(count*3==countwo);

    assert(count3==count);


    int root=0;

    int *countbuf=new int[nparts];

    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumNodes*3;
    }

    int recvNum=0;

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the node number



    /*
    int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
            */

    cout<<"RecvNum: "<<recvNum<<endl;

    int recvcount=recvNum;

    double *recvbuf=new double[recvcount];

    MPI_Scatterv(coordata,sendcounts,displs,MPI_DOUBLE,recvbuf,recvcount,MPI_DOUBLE,rank,comm);//send the coord ;

    if(coordata!=nullptr)
        delete []coordata;

    coordata=nullptr;

    construcTetras.NumNodes=recvNum/3;
    construcTetras.nodes=new Node[construcTetras.NumNodes];

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].coord.x=recvbuf[i*3+0];
        construcTetras.nodes[i].coord.y=recvbuf[i*3+1];
        construcTetras.nodes[i].coord.z=recvbuf[i*3+2];
    }
    if(recvbuf!=nullptr)
        delete []recvbuf;
    recvbuf=nullptr;

    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumNodes;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumNodes;
        // count=count+tetrasPart[i].NumNodes;
    }

    int NumNodes=recvcount/3;

    int *nodeGlobalID=new int[NumNodes];
    int *recvProcSize=new int[NumNodes];


    //for(int i=0;i<)

    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,nodeGlobalID,NumNodes,MPI_INT,rank,comm);//send the global index of the node ;

    MPI_Scatterv(procSize,sendcounts,displs,MPI_INT,recvProcSize,NumNodes,MPI_INT,rank,comm);//send the procs size of every node ;


    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].NumProcs=recvProcSize[i];
    }

    if(procSize!=nullptr)
    {
        delete []procSize;
        procSize=nullptr;
    }

    int allSetSize=0;
    int *sumProcs=new int[nparts];
    for(int i=0;i<nparts;i++)
    {
        sumProcs[i]=0;
    }

    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            sumProcs[i]+=tetrasPart[i].nodes[j].procs.size();
           // allSetSize+=tetrasPart[i].nodes[j].procs.size();
        }
        allSetSize+=sumProcs[i];
    }
    int setSizePart=-1;
    MPI_Scatter(sumProcs,1,MPI_INT,&setSizePart,1,MPI_INT,rank,comm);//send the size of procs of every part

    int *allSetValue=new int[allSetSize];
    set<int>::iterator setIter;
    count=0;
    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumNodes;j++)
        {
            for(setIter=tetrasPart[i].nodes[j].procs.begin();setIter!=tetrasPart[i].nodes[j].procs.end();setIter++)
            {
                int setValue=*setIter;
                allSetValue[count]=setValue;
                count++;

            }
        }
    }

    offset=0;
    for(int i=0;i<nparts;i++)
    {
        displs[i]=offset;//need to check
        offset+=sumProcs[i];
    }
    int *recvSetValue=new int[setSizePart];
    MPI_Scatterv(allSetValue,sumProcs,displs,MPI_INT,recvSetValue,setSizePart,MPI_INT,rank,comm);//send the value of the procs



    int allCount=0;
    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        int NumProcs=construcTetras.nodes[i].NumProcs;
       // count=0;
        for(int j=0;j<NumProcs;j++)
        {
            construcTetras.nodes[i].procs.insert(recvSetValue[allCount]);
            allCount++;
        }
    }

    assert(allCount==setSizePart);

    if(recvSetValue!=nullptr)
        delete []recvSetValue;
    recvSetValue=nullptr;


    if(nodeGlobabIndex!=nullptr)
        delete []nodeGlobabIndex;
    nodeGlobabIndex=nullptr;

    if(procSize!=nullptr)
        delete []procSize;
    procSize=nullptr;

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].index=nodeGlobalID[i];
       // construcTetras.nodes[i].procs.
    }

    if(nodeGlobalID!=nullptr)
        delete []nodeGlobalID;

    nodeGlobalID=nullptr;


    int sumPrism=0;
    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumPrsm*6;
        sumPrism+=tetrasPart[i].NumPrsm*6;
    }

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the prism number
    //add
        construcTetras.NumPrsm=recvNum/6;
 cout<<"PrismNum: "<<recvNum/6<<endl;
     construcTetras.pPrisms=new PRISM[construcTetras.NumPrsm];
   int *prismsVertices=new int[sumPrism];

   count=0;

   for(int i=0;i<nparts;i++)
   {
       for(int j=0;j<tetrasPart[i].NumPrsm;j++)
       {
//           for(int k=0;k<6;++k){
//               cout<<tetrasPart[i].pPrisms[j].vertices[k]<<"  ";

//           }
 //          cout<<endl;
           prismsVertices[count*6+0]=tetrasPart[i].pPrisms[j].vertices[0];
           prismsVertices[count*6+1]=tetrasPart[i].pPrisms[j].vertices[1];
           prismsVertices[count*6+2]=tetrasPart[i].pPrisms[j].vertices[2];
           prismsVertices[count*6+3]=tetrasPart[i].pPrisms[j].vertices[3];
           prismsVertices[count*6+4]=tetrasPart[i].pPrisms[j].vertices[4];
           prismsVertices[count*6+5]=tetrasPart[i].pPrisms[j].vertices[5];
           count++;

       }
   }

   offset=0;
   for(int i=0;i<nparts;i++)
   {
       sendcounts[i]=tetrasPart[i].NumPrsm*6;
       displs[i]=offset;//need to check
       offset+=tetrasPart[i].NumPrsm*6;
       // count=count+tetrasPart[i].NumNodes;
   }

   int *vPrisms=new int[recvNum];

   MPI_Scatterv(prismsVertices,sendcounts,displs,MPI_INT,vPrisms,recvNum,MPI_INT,rank,comm);//send the prisms  ;

   if(prismsVertices!=nullptr)
       delete []prismsVertices;
   prismsVertices=nullptr;

   for(int i=0;i<construcTetras.NumPrsm;i++)
   {
       construcTetras.pPrisms[i].partMarker=rank;
       construcTetras.pPrisms[i].vertices[0]=vPrisms[i*6+0];
       construcTetras.pPrisms[i].vertices[1]=vPrisms[i*6+1];
       construcTetras.pPrisms[i].vertices[2]=vPrisms[i*6+2];
       construcTetras.pPrisms[i].vertices[3]=vPrisms[i*6+3];
       construcTetras.pPrisms[i].vertices[4]=vPrisms[i*6+4];
       construcTetras.pPrisms[i].vertices[5]=vPrisms[i*6+5];
}
   if(vPrisms!=nullptr)
       delete []vPrisms;
   vPrisms=nullptr;



    int sumTetras=0;
    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumTetras*4;
        sumTetras+=tetrasPart[i].NumTetras*4;
    }

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the tetras number

    construcTetras.NumTetras=recvNum/4;

    cout<<"TetrasNum: "<<recvNum/4<<endl;

    construcTetras.pTetras=new TETRAS[construcTetras.NumTetras];

    int *tetrasVertices=new int[sumTetras];

    count=0;

    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumTetras;j++)
        {
            tetrasVertices[count*4+0]=tetrasPart[i].pTetras[j].vertices[0];
            tetrasVertices[count*4+1]=tetrasPart[i].pTetras[j].vertices[1];
            tetrasVertices[count*4+2]=tetrasPart[i].pTetras[j].vertices[2];
            tetrasVertices[count*4+3]=tetrasPart[i].pTetras[j].vertices[3];
            count++;

        }
    }


    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumTetras*4;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumTetras*4;
        // count=count+tetrasPart[i].NumNodes;
    }

    int *vTetras=new int[recvNum];

    MPI_Scatterv(tetrasVertices,sendcounts,displs,MPI_INT,vTetras,recvNum,MPI_INT,rank,comm);//send the tetras  ;


    if(tetrasVertices!=nullptr)
        delete []tetrasVertices;
    tetrasVertices=nullptr;

    for(int i=0;i<construcTetras.NumTetras;i++)
    {
        construcTetras.pTetras[i].partMarker=rank;
        construcTetras.pTetras[i].vertices[0]=vTetras[i*4+0];
        construcTetras.pTetras[i].vertices[1]=vTetras[i*4+1];
        construcTetras.pTetras[i].vertices[2]=vTetras[i*4+2];
        construcTetras.pTetras[i].vertices[3]=vTetras[i*4+3];

    }
    if(vTetras!=nullptr)
        delete []vTetras;
    vTetras=nullptr;

    int sumTris=0;
    for(int i=0;i<nparts;i++)
    {
        countbuf[i]=tetrasPart[i].NumTris*3;
        sumTris+=tetrasPart[i].NumTris*3;
    }

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,rank,comm);//send the tris number

    construcTetras.NumTris=recvNum/3;

    cout<<"TrisNum: "<<recvNum/3<<endl;

    construcTetras.pTris=new TRI[construcTetras.NumTris];

    int *triVertices=new int[sumTris];

    count=0;
    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumTris;j++)
        {
            triVertices[count*3+0]=tetrasPart[i].pTris[j].vertices[0];
            triVertices[count*3+1]=tetrasPart[i].pTris[j].vertices[1];
            triVertices[count*3+2]=tetrasPart[i].pTris[j].vertices[2];
            count++;
        }
    }
    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumTris*3;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumTris*3;
        // count=count+tetrasPart[i].NumNodes;
    }

    int *vTris=new int[recvNum];

    MPI_Scatterv(triVertices,sendcounts,displs,MPI_INT,vTris,recvNum,MPI_INT,rank,comm);//send the tris  ;

    if(triVertices!=nullptr)
        delete []triVertices;
    triVertices=nullptr;

    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].partMarker=rank;
        construcTetras.pTris[i].vertices[0]=vTris[i*3+0];
        construcTetras.pTris[i].vertices[1]=vTris[i*3+1];
        construcTetras.pTris[i].vertices[2]=vTris[i*3+2];

    }

    if(vTris!=nullptr)
        delete []vTris;
    vTris=nullptr;
    count=0;

    sumTris=0;
    offset=0;
    for(int i=0;i<nparts;i++)
    {
        sendcounts[i]=tetrasPart[i].NumTris*3;
        displs[i]=offset;//need to check
        offset+=tetrasPart[i].NumTris*3;
        sumTris+=tetrasPart[i].NumTris*3;
        // count=count+tetrasPart[i].NumNodes;
    }


    int *attributs=new int[sumTris];// the first value is iSurf , the second value is iOppoProc
    for(int i=0;i<nparts;i++)
    {
        for(int j=0;j<tetrasPart[i].NumTris;j++)
        {
            attributs[count*3+0]=tetrasPart[i].pTris[j].iSurf;
            attributs[count*3+1]=tetrasPart[i].pTris[j].iOppoProc;
            attributs[count*3+2]=tetrasPart[i].pTris[j].iSurface;
            count++;
        }
    }
    assert(count==sumTris/3);

    cout<<"Triangle Number is: "<<sumTris/3<<endl;

    int *recvAttr=new int[3*construcTetras.NumTris];

    recvNum=construcTetras.NumTris*3;

    MPI_Scatterv(attributs,sendcounts,displs,MPI_INT,recvAttr,recvNum,MPI_INT,rank,comm);//send the tris  ;//send the tris attributes


    if(attributs!=nullptr)
        delete []attributs;
    attributs=nullptr;
    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].iSurf=recvAttr[i*3+0];
        construcTetras.pTris[i].iOppoProc=recvAttr[i*3+1];
        construcTetras.pTris[i].iSurface=recvAttr[i*3+2];

    }

    if(recvAttr!=nullptr)
        delete []recvAttr;
    recvAttr=nullptr;


    int *refines=new int[nparts];

    for(int i=0;i<nparts;i++)
    {
        refines[i]=refineSize;
    }
    int recvRefine=-1;
    MPI_Scatter(refines,1,MPI_INT,&recvRefine,1,MPI_INT,rank,comm);

    if(refines!=nullptr)
        delete []refines;
    refines=nullptr;

    stringstream ss;
    for (int i = 0; i < bcstring.size(); ++i)
        ss << bcstring[i] << "\n";
    const string& bufStr = ss.str();
    cout << bufStr << endl;
    int stringSize=bufStr.size();

    //int *strSize=new int[];

   // MPI_Scatter()
    MPI_Bcast(&stringSize,1,MPI_INT,rank,comm);
    MPI_Bcast(const_cast<char*>(bufStr.c_str()), stringSize, MPI_CHAR,rank,comm);

    cout<<"MPI starting..."<<endl;
    return 1;

}

int workerProc(int rank, HYBRID_MESH&construcTetras, int &refineSize, vector<string> &bcstring, MPI_Comm comm)
{

    int countbuf[2];

    int recvNum=0;

    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);

    cout<<"RecvNum: "<<recvNum<<endl;

    double *recvbuf=new double[recvNum];

    double coordata[2];

    int sendcounts[2];

    int displs[2];

    MPI_Scatterv(coordata,sendcounts,displs,MPI_DOUBLE,recvbuf,recvNum,MPI_DOUBLE,0,comm);//receive the coord ;

    construcTetras.NumNodes=recvNum/3;
    construcTetras.nodes=new Node[construcTetras.NumNodes];

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].localID=i;
        construcTetras.nodes[i].coord.x=recvbuf[i*3+0];
        construcTetras.nodes[i].coord.y=recvbuf[i*3+1];
        construcTetras.nodes[i].coord.z=recvbuf[i*3+2];
    }

    if(recvbuf!=nullptr)
    {
        delete []recvbuf;
        recvbuf=nullptr;
    }

    int NumNodes=recvNum/3;



    int *nodeGlobalID=new int[NumNodes];

    cout<<"Numnodes: "<<NumNodes<<endl;

    int *nodeGlobabIndex=new int[2];

    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,nodeGlobalID,NumNodes,MPI_INT,0,comm);//receive the global index of the node ;

    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].index=nodeGlobalID[i];
    }

    if(nodeGlobalID!=nullptr)
    {
        delete []nodeGlobalID;
        nodeGlobalID=nullptr;
    }

    int *recvProcSize=new int[NumNodes];
    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,recvProcSize,NumNodes,MPI_INT,0,comm);//receive the procs size of every node ;


    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        construcTetras.nodes[i].NumProcs=recvProcSize[i];
    }

    if(recvProcSize!=nullptr)
        delete []recvProcSize;
    recvProcSize=nullptr;

    int setSizePart=-1;
    MPI_Scatter(countbuf,1,MPI_INT,&setSizePart,1,MPI_INT,0,comm);//receive the size of procs of every part

    int *recvSetValue=new int[setSizePart];
    MPI_Scatterv(nodeGlobabIndex,sendcounts,displs,MPI_INT,recvSetValue,setSizePart,MPI_INT,0,comm);//receive the value of the procs



    int allCount=0;
    for(int i=0;i<construcTetras.NumNodes;i++)
    {
        int NumProcs=construcTetras.nodes[i].NumProcs;
 //       count=0;
        for(int j=0;j<NumProcs;j++)
        {
            construcTetras.nodes[i].procs.insert(recvSetValue[allCount]);
            allCount++;
        }
    }

    assert(allCount==setSizePart);

    if(recvSetValue!=nullptr)
        delete []recvSetValue;
    recvSetValue=nullptr;


    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);//receive the prisms number

    construcTetras.NumPrsm=recvNum/6;

    cout<<"PrismNum: "<<recvNum/6<<endl;

    construcTetras.pPrisms=new PRISM[construcTetras.NumPrsm];

    int *vPrism=new int[recvNum];

    int *prismsVertices=new int[2];

    MPI_Scatterv(prismsVertices,sendcounts,displs,MPI_INT,vPrism,recvNum,MPI_INT,0,comm);//receive the prisms  ;


    for(int i=0;i<construcTetras.NumPrsm;i++)
    {
        construcTetras.pPrisms[i].partMarker=rank;
        construcTetras.pPrisms[i].localID=i;
        construcTetras.pPrisms[i].vertices[0]=vPrism[i*6+0];
        construcTetras.pPrisms[i].vertices[1]=vPrism[i*6+1];
        construcTetras.pPrisms[i].vertices[2]=vPrism[i*6+2];
        construcTetras.pPrisms[i].vertices[3]=vPrism[i*6+3];
        construcTetras.pPrisms[i].vertices[4]=vPrism[i*6+4];
        construcTetras.pPrisms[i].vertices[5]=vPrism[i*6+5];

    }

    if(vPrism!=nullptr)
    {
        delete []vPrism;
        vPrism=nullptr;
    }


    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);//receive the tetras number

    construcTetras.NumTetras=recvNum/4;

    cout<<"TetrasNum: "<<recvNum/4<<endl;

    construcTetras.pTetras=new TETRAS[construcTetras.NumTetras];

    int *vTetras=new int[recvNum];

    int *tetrasVertices=new int[2];

    MPI_Scatterv(tetrasVertices,sendcounts,displs,MPI_INT,vTetras,recvNum,MPI_INT,0,comm);//receive the tetras  ;


    for(int i=0;i<construcTetras.NumTetras;i++)
    {
        construcTetras.pTetras[i].partMarker=rank;
        construcTetras.pTetras[i].localID=i;
        construcTetras.pTetras[i].vertices[0]=vTetras[i*4+0];
        construcTetras.pTetras[i].vertices[1]=vTetras[i*4+1];
        construcTetras.pTetras[i].vertices[2]=vTetras[i*4+2];
        construcTetras.pTetras[i].vertices[3]=vTetras[i*4+3];

    }

    if(vTetras!=nullptr)
    {
        delete []vTetras;
        vTetras=nullptr;
    }



    MPI_Scatter(countbuf,1,MPI_INT,&recvNum,1,MPI_INT,0,comm);//recv the tris number

    construcTetras.NumTris=recvNum/3;

    cout<<"TrisNum: "<<recvNum/3<<endl;

    construcTetras.pTris=new TRI[construcTetras.NumTris];

    int triVertices[2];

    int *vTris=new int[construcTetras.NumTris*3];

  //  cout<<"worker zhvliu..............................."<<endl;

    MPI_Scatterv(triVertices,sendcounts,displs,MPI_INT,vTris,recvNum,MPI_INT,0,comm);//recv the tris  ;

    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].partMarker=rank;
        construcTetras.pTris[i].vertices[0]=vTris[i*3+0];
        construcTetras.pTris[i].vertices[1]=vTris[i*3+1];
        construcTetras.pTris[i].vertices[2]=vTris[i*3+2];

    }

    if(vTris!=nullptr)
    {
        delete []vTris;
        vTris=nullptr;
    }

    int attributs[3];


    recvNum=construcTetras.NumTris*3;

    int *recvAttr=new int[recvNum];



    MPI_Scatterv(attributs,sendcounts,displs,MPI_INT,recvAttr,recvNum,MPI_INT,0,comm);//recv send the tris attributes


    for(int i=0;i<construcTetras.NumTris;i++)
    {
        construcTetras.pTris[i].iSurf=recvAttr[i*3+0];
        construcTetras.pTris[i].iOppoProc=recvAttr[i*3+1];
        construcTetras.pTris[i].iSurface=recvAttr[i*3+2];

    }


    MPI_Scatter(attributs,1,MPI_INT,&refineSize,1,MPI_INT,0,comm);

    if(recvAttr!=nullptr)
    {
        delete []recvAttr;
        recvAttr=nullptr;
    }

    // cout << bufStr << endl;
    int stringSize;
    MPI_Bcast(&stringSize,1,MPI_INT,0,comm);

    cout<<"StringSize == "<<stringSize<<endl;

    char *buf = new char[stringSize + 1];
    MPI_Bcast(buf, stringSize, MPI_CHAR,0,comm);
    buf[stringSize] = '\0';

    stringstream ss;
    ss << string(buf);
    string bufStr;
    while (getline(ss, bufStr))
        bcstring.push_back(bufStr);

    // cout << string(buf) << "Received!\n";



    delete [] buf;

    return 1;
}
int meshImproved(HYBRID_MESH&mesh)
{

    /*
    //meshGen3D_memo();
    //meshGen3D_memo();//need to modified
    //meshGen3D_memo();
    int nNodes=mesh.NumNodes;
    int nTetras=mesh.NumTetras;
    int nFacets=mesh.NumTris;
    double *x=new double[nNodes];
    double *y=new double[nNodes];
    double *z=new double[nNodes];
    int *tetras=new int[nTetras*4];
    int *facets=new int[nFacets*3];

    double *outx=nullptr;
    double *outy=nullptr;
    double *outz=nullptr;
    int outNodes=0;
    int outNTetras=0;
    int *ouTetras=nullptr;
    meshImproved(nNodes,x,y,z,nTetras,tetras,nFacets,facets,outNodes,outx,outy,outz,outNTetras,ouTetras);
    */
    return 1;
}










