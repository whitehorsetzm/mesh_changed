#include "dataio.h"
#include "stdio.h"
#include "ctype.h"
#include <string.h>
#include "stdlib.h"
#include "assert.h"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include "cgnslib.h"

using namespace std;
bool isNullOrComment(char* chLine)
{
    char *str = chLine;
    int numOfValidC = 0;
    while (*str != '\0' && numOfValidC == 0)
    {
        if (!isspace(*str))
        {
            if (numOfValidC == 0 && *str == '#')
                return true;
            numOfValidC++;
        }
        str++;
    }

    return numOfValidC == 0;
}


int readValidLine(FILE *fp, char *chLine, int len)
{
    if (fp && !feof(fp))
    {
        do
        {
            fgets(chLine, len, fp);
            if (!isNullOrComment(chLine))
                return 0; // succeed
        }
        while (!feof(fp));
    }

    return 1; // fail
}

int eraseWhiteSpace(char *chLine)
{
    char *strHd = chLine;
    char *strTl = chLine + strlen(chLine) - 1;
    int i = 0;

    while (strHd < strTl)
    {
        /* 锟接猴拷锟斤拷锟斤拷 */
        if (!isspace(*strHd))
            break;
        strHd++;
    }

    while (strHd < strTl)
    {
        /* 锟斤拷前锟斤拷锟斤拷 */
        if (!isspace(*strTl))
            break;
        strTl--;
    }

    if (strTl > strHd)
    {
        if (strHd > chLine)
        {
            i = 0;
            do
            {
                chLine[i++] = *(strHd++);
            }
            while (strHd <= strTl);
            chLine[i] = '\0';
        }
        else
            *(strTl+1) = '\0';
    }

    return 1; // fail
}

int setupCellNeig(int nNodes, int nElems, TETRAS *pBKGElems)
{

    int nAllocFaceSize = nNodes * 10, nAllocNodeFaceSize = 20384;
    int nFaceSize = 0;
    int *vecRefIntFHash = NULL;
    InterFace * vecInterFaces = NULL, *vecInterFaces_Temp = NULL;
    int *ndIFaces = NULL, *ndIFaces_Temp = NULL;
    int errCode = 0;
    int cellIdx, faceAdIdx, lftCell, rgtCell;
    int facNdIdx1, facNdIdx2, facNdIdx3, minFacNdIdx, faceIt;
    int ndIFaceSize = 0;
    int i, j, k, nCommon;
    InterFace faceAd;
    static int cf[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    int ndSize, clSize;
    //BKGElem *pElem = NULL, *pLft = NULL, *pRgt = NULL;

    TETRAS *pElem = nullptr, *pLft=nullptr, *pRgt=nullptr;

    vecRefIntFHash = (int *) malloc(sizeof(int)*nNodes);
    vecInterFaces = (InterFace *) malloc(sizeof(InterFace)*nAllocFaceSize);
    ndIFaces = (int *) malloc(sizeof(int)*nAllocNodeFaceSize);
    if (!vecRefIntFHash || !vecInterFaces || !ndIFaces)
    {
        errCode = -1;
        goto FAIL;
    }

    ndSize = nNodes;
    clSize = nElems;

    for (i = 0; i < ndSize; i++)
        vecRefIntFHash[i] = -1;
    for (cellIdx = 0; cellIdx < clSize; cellIdx++)
    {
        pElem = &pBKGElems[cellIdx];
        for (i = 0; i <= BKG_MESH_DIM; i++)
            pElem->neighbors[i] = -1;
    }

    for (cellIdx = 0; cellIdx < clSize; cellIdx++)
    {
        if (cellIdx % 1000000 == 0 || cellIdx == clSize - 1)
            printf("%%%.2f.\n", 100.0*(cellIdx)/clSize);

        pElem = &pBKGElems[cellIdx];

        for (i = 0; i <= BKG_MESH_DIM; i++) {
            facNdIdx1 = pElem->vertices[cf[i][0]];
            facNdIdx2 = pElem->vertices[cf[i][1]];
            facNdIdx3 = pElem->vertices[cf[i][2]];

            minFacNdIdx = facNdIdx1 < facNdIdx2 ? facNdIdx1 : facNdIdx2;
            minFacNdIdx = minFacNdIdx < facNdIdx3 ? minFacNdIdx : facNdIdx3;
            //		nodeInterFace(minFacNdIdx, ndIFaces, &ndIFaceSize, vecInterFaces, vecRefIntFHash);

            j = 0;
            faceIt = vecRefIntFHash[minFacNdIdx];
            while (faceIt >= 0)
            {
                if (j >= nAllocNodeFaceSize)
                {
                    nAllocNodeFaceSize += MAX_VALUE(100, (int)(nAllocNodeFaceSize * 0.1));
                    ndIFaces_Temp = (int*)realloc(ndIFaces, sizeof(int)*nAllocNodeFaceSize);
                    if (!ndIFaces_Temp)
                    {
                        errCode = -1;
                        goto FAIL;
                    }
                    ndIFaces = ndIFaces_Temp;
                }
                ndIFaces[j++] = faceIt;
                faceIt= vecInterFaces[faceIt].hashNxt;
            }
            ndIFaceSize = j;


            nCommon = 0;
            for (j = 0;	j < ndIFaceSize; j++) {
                nCommon = 0;
                for (k = 0; k < 3; k++)
                    if (vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx1 ||
                            vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx2 ||
                            vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx3)
                        nCommon++;
                if (nCommon >= 3)
                    break;
            }

            if (nCommon < 3) {
                /* 没锟斤拷锟揭碉拷 */
                if (nFaceSize >= nAllocFaceSize)
                {
                    nAllocFaceSize += MAX_VALUE(100, (int)(nAllocFaceSize * 0.1));
                    vecInterFaces_Temp = (InterFace*)realloc(vecInterFaces, sizeof(InterFace)*nAllocFaceSize);
                    if (!vecInterFaces_Temp)
                    {
                        errCode = -1;
                        goto FAIL;
                    }
                    vecInterFaces = vecInterFaces_Temp;
                }

                faceAd.conn[0] = facNdIdx1;
                faceAd.conn[1] = facNdIdx2;
                faceAd.conn[2] = facNdIdx3;
                faceAd.lftCell = cellIdx;
                faceAd.rgtCell = -1;

                faceAdIdx = nFaceSize;
                faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];
                vecRefIntFHash[minFacNdIdx] = faceAdIdx;
                vecInterFaces[nFaceSize++] = faceAd;
            }
            else {

                if(!(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                     vecInterFaces[ndIFaces[j]].rgtCell < 0))
                    cout<<cellIdx<<endl;

                assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                        vecInterFaces[ndIFaces[j]].rgtCell < 0);


                vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;

                lftCell = vecInterFaces[ndIFaces[j]].lftCell;
                rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;
                pLft = &pBKGElems[lftCell];
                pRgt = &pBKGElems[rgtCell];
                for (k = 0; k < 4; k++)
                    if (pLft->vertices[k] != facNdIdx1 &&
                            pLft->vertices[k] != facNdIdx2 &&
                            pLft->vertices[k] != facNdIdx3)
                        break;
                assert(k < 4);
                pRgt->neighbors[i] = lftCell;
                pLft->neighbors[k] = rgtCell;
            }
        }
    }

    goto END;
FAIL:
END:
    if (vecRefIntFHash)
        free(vecRefIntFHash);
    if (vecInterFaces)
        free(vecInterFaces);
    if (ndIFaces)
        free(ndIFaces);

    return errCode; /* S_OK */
}


int readVTKPLSFile(char *fname,HYBRID_MESH& file)
{

    FILE *fp = NULL;
    Node *pBKGNodes = nullptr;
    TETRAS *pBKGElems = nullptr;
    Node *pNode=nullptr;
    TETRAS *pElem=nullptr;

    // ELEMENT<4,4> ele;

    //ELEMENT<4,4>*p = new ELEMENT<4,4>();

    int nNodeNum = 0, nElemNum = 0;
    char chLine[MAX_FILE_LINE];
    int i, j, conn, nValueNum;
    int errCode = 0;


    char vtkfile[256];
    sprintf(vtkfile, "%s%s",fname,".vtk");

    if (!vtkfile)
    {
        printf("No VTK file name is given.\n");
        return 2;
    }
    fp = fopen(vtkfile, "r");
    if (!(fp = fopen(vtkfile, "r")))
    {
        printf("Cannot open file %s.\n", vtkfile);
        return 2;
    }

    /* points */
    do
    {
        readValidLine(fp, chLine, MAX_FILE_LINE);
        eraseWhiteSpace(chLine);
        if (strstr(chLine, "POINTS") != NULL)
            break;
    }
    while (!feof(fp));

    if (feof(fp))
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    sscanf(chLine, "POINTS %d", &nNodeNum);

    file.NumNodes=nNodeNum;

    if (nNodeNum <= 0)
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    if (!(file.nodes = new Node[nNodeNum]))
    {
        errCode = -1; /* not enough memory */
        goto FAIL;
    }

    pBKGNodes=file.nodes;

    for (i = 0; i < nNodeNum; i++)
    {
        pNode = &pBKGNodes[i];
        pNode->index=i;
        fscanf(fp, "%lf", &pNode->coord.x);
        fscanf(fp, "%lf", &pNode->coord.y);
        fscanf(fp, "%lf", &pNode->coord.z);
    }

    /* elements */
    do
    {
        readValidLine(fp, chLine, MAX_FILE_LINE);
        eraseWhiteSpace(chLine);
        if (strstr(chLine, "CELLS") != NULL)
            break;
    }
    while (!feof(fp));

    if (feof(fp))
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    sscanf(chLine, "CELLS %d", &nElemNum);

    file.NumTetras=nElemNum;

    if (nNodeNum <= 0)
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    if (!(file.pTetras = new TETRAS[nElemNum]))
    {
        errCode = -1; /* not enough memory */
        goto FAIL;
    }

    pBKGElems=file.pTetras;

    for (i = 0; i < nElemNum; i++)
    {
        fscanf(fp, "%d", &nValueNum);
        if (nValueNum != BKG_MESH_DIM + 1)
        {
            errCode = 3; /* Invalid file */
            goto FAIL;
        }

        pElem = &pBKGElems[i];
        pElem->index=i;
        for (j = 0; j <= BKG_MESH_DIM; j++)
        {
            fscanf(fp, "%d", &conn);
            //if (conn <= 0 || conn > nNodeNum)
            if (conn < 0 || conn >= nNodeNum)
            {
                errCode = 3; /* Invalid file */
                goto FAIL;
            }
            pElem->vertices[j] = conn;/*conn - 1;*/
            pElem->neighbors[j] = -1;
        }
    }

    /* reconstruct neighbouring indices */
    errCode = setupCellNeig(nNodeNum, nElemNum, pBKGElems);
    if (errCode != 0 && errCode != 1)
        goto FAIL;


    goto END;
FAIL:

END:
    fclose(fp);


//read pls file
    char plsfile[256];
    sprintf(plsfile, "%s%s",fname,".pls");
    fp = fopen(plsfile, "r");
    if (!plsfile)
    {
        printf("Error: cann't open file %s!\n", plsfile);
        return 0;
    }

    int NumTris=0; int NumPoints=0;  int tmp=-1;

    fgets(chLine, MAX_FILE_LINE, fp);
    sscanf(chLine, "%d %d", &NumTris, &NumPoints);
    //fscanf(fp, "%d %d %d %d %d %d\n", &NumTris, &NumPoints, &tmp, &tmp, &tmp, &tmp);

   // *nbele = nel;
   // (*bele) = new int[nel * 4];
	file.NumTris=NumTris;
	file.pTris=new TRI[NumTris];

    float ftmp=0;
    for (int i=0; i<NumPoints; i++)
    {
        fscanf(fp, "%d %f %f %f\n", &tmp, &ftmp, &ftmp, &ftmp);
    }

    for (int i=0; i<NumTris; i++)
    {
        int id1, id2, id3, fidx;
        fscanf(fp, "%d %d %d %d %d\n", &tmp, &id1, &id2, &id3, &fidx);
        file.pTris[i].vertices[0]=id1-1;
        file.pTris[i].vertices[1]=id2-1;
        file.pTris[i].vertices[2]=id3-1;
        file.pTris[i].iSurf=fidx-1;
        file.pTris[i].index=i;
    }

    fclose(fp);
    fp = nullptr;

    findiCellFast(file);


    return errCode;
}


int findiCellFast(HYBRID_MESH &file)
{
    int neigbor=-2;
    int count=0;


    set<int>cellIDs;

    map<string,int> trimap;

    int v[3];
    for(int i=0;i<file.NumTris;i++)
    {
        v[0]=file.pTris[i].vertices[0];
        v[1]=file.pTris[i].vertices[1];
        v[2]=file.pTris[i].vertices[2];
        sort(v,v+3);
        string temp=IntToString(v[0])+"_"+IntToString(v[1])+"_"+IntToString(v[2]);
        trimap[temp]=i;
    }


    map<string,int>::iterator mapIter;
    for(int i=0;i<file.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            neigbor=file.pTetras[i].neighbors[j];

            if(neigbor==-1)
            {
               // count++;

                cellIDs.insert(i);
                count=0;
                for(int k=0;k<4;k++)
                {
                    if(k!=j)
                    {
                        v[count]=file.pTetras[i].vertices[k];
                        count++;

                    }
                }
                sort(v,v+3);
                string temp=IntToString(v[0])+"_"+IntToString(v[1])+"_"+IntToString(v[2]);
                mapIter=trimap.find(temp);
                if(mapIter==trimap.end())
                {
                    cout<<"Error in finding iCell!"<<endl;
                }
                else
                {
                    int triID=mapIter->second;
                    file.pTris[triID].iCell=i;
                }
                assert(count==3);
            }
        }
    }

    cout<<"Finding iCell finished!"<<endl;
    return 1;
}

int findiCell(HYBRID_MESH&file)
{
    int neigbor=-2;
    int count=0;


    set<int>cellIDs;


    for(int i=0;i<file.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            neigbor=file.pTetras[i].neighbors[j];
            if(neigbor==-1)
            {
                count++;
                cellIDs.insert(i);
            }
        }
    }



    set<int>::iterator iter;

    //if(iter==cellIDs.end())
    //{
    //    cout<<"Error "<<endl;
    //}

    cout<<"NumTri :" <<file.NumTris<<endl;
   // assert(count==file.NumTri);


    count=0;
    int tempID=-1;
    int triVer[3];
    int tetVer[4];
    int triCompare[3];
    for(int i=0;i<file.NumTris;i++)
    {
        triVer[0]=file.pTris[i].vertices[0];
        triVer[1]=file.pTris[i].vertices[1];
        triVer[2]=file.pTris[i].vertices[2];



        bool same=false;
        for(iter=cellIDs.begin();iter!=cellIDs.end();iter++)
        {
            tempID=*iter;

           triCompare[0]=file.pTetras[tempID].vertices[0];
           triCompare[1]=file.pTetras[tempID].vertices[1];
           triCompare[2]=file.pTetras[tempID].vertices[2];
           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }

           triCompare[0]=file.pTetras[tempID].vertices[0];
           triCompare[1]=file.pTetras[tempID].vertices[1];
           triCompare[2]=file.pTetras[tempID].vertices[3];

           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }
           triCompare[0]=file.pTetras[tempID].vertices[0];
           triCompare[1]=file.pTetras[tempID].vertices[2];
           triCompare[2]=file.pTetras[tempID].vertices[3];

           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }
           triCompare[0]=file.pTetras[tempID].vertices[1];
           triCompare[1]=file.pTetras[tempID].vertices[2];
           triCompare[2]=file.pTetras[tempID].vertices[3];

           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }

        }
        if(same==false)
        {
            cout<<"Error in finding the owner tetras! Tri v is "<<triVer[0]<<" "
               <<triVer[1]<<" "<<triVer[2]<<endl;
            exit(0);
        }
    }
    cout<<"Finding iCell finished!"<<endl;

    return 1;
}

bool tetrasContainTriangle(int *tri, int *tetras)
{
    /*
    sort(tri,tri+3);
    sort(tetras,tetras+4);

    if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
    else if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
    else if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
    else if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
	*/
    return true;
    
}
bool isSameTriangle(int *tri,int *triCom)
{
    sort(tri,tri+3);
    sort(triCom,triCom+3);
    if(tri[0]==triCom[0]&&tri[1]==triCom[1]&&tri[2]==triCom[2])
    {
        return true;
    }
    else
        return false;
}
int writeTriangleVTKFile(char *filename,HYBRID_MESH&file)
{
    int i, j = 0, k = 0, iFac = 0, ii1 = 0, ii2 = 0, ii3 = 0;
    int faceNumber = 0;
   // MYPOINT pnt;
  //  DesFacet *pDF = NULL;
   // bool *visited = NULL;

	int elemSize = file.NumTris, nodeSize = file.NumNodes;

    FILE *fout = fopen(filename,"w");
    fprintf(fout,"# vtk DataFile Version 2.0\n");
    fprintf(fout,"background mesh\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout,"POINTS %d double\n", nodeSize);


    for (i = 0; i < nodeSize; i++)
    {
        int index = i + 1;
        double x = file.nodes[i].coord.x;
        double y = file.nodes[i].coord.y;
        double z = file.nodes[i].coord.z;
        fprintf(fout,"%lf %lf %lf\n",x, y, z);

    }
    fprintf(fout,"CELLS %d %d \n", elemSize, (elemSize) * 4);

    for (i = 0; i < elemSize; i++)
    {

        int index = i + 1;
        int form0 = file.pTris[i].vertices[0];///////////+1
        int form1 = file.pTris[i].vertices[1];
        int form2 = file.pTris[i].vertices[2];
        //int form3 = file.pTetras[i].vertices[3];

        fprintf(fout,"3 %d %d %d\n",form0,form1,form2);

    }
    fprintf(fout,"CELL_TYPES %d \n", elemSize);
    for (i = 0; i < elemSize; i++)
    {
        fprintf(fout,"%d\n", 5);
    }

	fprintf(fout,"\nPOINT_DATA %d \n", nodeSize);

#if 0
    fprintf(fout,"SCALARS procs_begin int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < nodeSize; i++)
    {
        //double size = m_pNodes[i].space/m_scale;
       // fprintf(fout,"%d\n",file.nodes[i].index);
        fprintf(fout,"%d\n",*file.nodes[i].procs.begin());
    }
#endif


    fprintf(fout,"SCALARS nodeID int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < nodeSize; i++)
    {
        //double size = m_pNodes[i].space/m_scale;
        fprintf(fout,"%d\n",file.nodes[i].index);
        // fprintf(fout,"%d\n",file.nodes[i].procs.size());
    }

	fprintf(fout,"\nCELL_DATA %d \n", elemSize);
    fprintf(fout,"SCALARS iOppoProc int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < elemSize; i++)
    {

       // double quality =  m_pElems[i].minAngle;
        fprintf(fout,"%d\n",file.pTris[i].iOppoProc);
      //  fprintf(fout,"%d\n",file.pTris[i].iCell);
    }


    fprintf(fout,"SCALARS iCell int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < elemSize; i++)
    {
        fprintf(fout,"%d\n",file.pTris[i].iCell);
    }

	fprintf(fout,"\nSCALARS facetID int\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i = 0; i < elemSize; i++)
	{
	   // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%lld\n",file.pTris[i].index);
	}

	fprintf(fout,"\nSCALARS iSurf int\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i = 0; i < elemSize; i++)
	{
	   // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%d\n",file.pTris[i].iSurf);
	}

	fclose(fout);
    return 1;
}

int writeVTKFile(char *filename,HYBRID_MESH&file)
{

    int i, j = 0, k = 0, iFac = 0, ii1 = 0, ii2 = 0, ii3 = 0;
    int faceNumber = 0;
   // MYPOINT pnt;
  //  DesFacet *pDF = NULL;
   // bool *visited = NULL;

    int elemSize = file.NumTetras, nodeSize = file.NumNodes;

    FILE *fout = fopen(filename,"w");
    fprintf(fout,"# vtk DataFile Version 2.0\n");
    fprintf(fout,"background mesh\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout,"POINTS %d double\n", nodeSize);


    for (i = 0; i < nodeSize; i++)
    {
        int index = i + 1;
        double x = file.nodes[i].coord.x;
        double y = file.nodes[i].coord.y;
        double z = file.nodes[i].coord.z;
        fprintf(fout,"%lf %lf %lf\n",x, y, z);

    }
    fprintf(fout,"CELLS %d %d \n", elemSize, (elemSize) * 5);

    for (i = 0; i < elemSize; i++)
    {

        int index = i + 1;
        int form0 = file.pTetras[i].vertices[0];///////////+1
        int form1 = file.pTetras[i].vertices[1];
        int form2 = file.pTetras[i].vertices[2];
        int form3 = file.pTetras[i].vertices[3];

        fprintf(fout,"4 %d %d %d %d\n",form0,form1,form2,form3);

    }
    fprintf(fout,"CELL_TYPES %d \n", elemSize);
    for (i = 0; i < elemSize; i++)
    {
        fprintf(fout,"%d\n", 10);
    }

	fprintf(fout,"\nPOINT_DATA %d \n", nodeSize);
	fprintf(fout,"SCALARS nodeID int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < nodeSize; i++)
    {
        //double size = m_pNodes[i].space/m_scale;
        fprintf(fout,"%d\n",file.nodes[i].index);
    }

	fprintf(fout,"\nCELL_DATA %d \n", elemSize);
    fprintf(fout,"SCALARS part int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < elemSize; i++)
    {
       // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%d\n",file.pTetras[i].partMarker);
	}


	fprintf(fout,"\nSCALARS cellID int\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i = 0; i < elemSize; i++)
	{
	   // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%d\n",file.pTetras[i].index);
	}


    fclose(fout);
    return 1;


    //return 1;
}
#if 1
int readCGNS(char*filename,HYBRID_MESH&mesh,vector<string>&bcstring)
{
  //  float x[21*17*9],y[21*17*9],z[21*17*9];
    cgsize_t isize[3][1];//ielem[20*16*8][8];
    int index_file,index_base,index_zone;
    cgsize_t irmin,irmax,istart,iend;
    int nsections,index_sect,nbndry,iparent_flag;
    cgsize_t iparentdata;
    char zonename[33],sectionname[33];
    CGNS_ENUMT(ElementType_t) itype;

/* READ X, Y, Z GRID POINTS FROM CGNS FILE */
/* open CGNS file for read-only */
    int filetype=-1;
    if (!cg_is_cgns(filename,&filetype))
        cout<<"CGNS FILE: type="<<filetype<<endl;
    else
        cout<<"NOT CGNS FILE."<<endl;

    cg_set_file_type(filetype);
    if (cg_open(filename,CG_MODE_READ,&index_file))
        cg_error_exit();
    cout<<"Open CGNS file"<<endl;
/* we know there is only one base (real working code would check!) */
    int nbases=-1;
    cg_nbases(index_file,&nbases);
    if(nbases!=1)
    {
        cout<<"The bases is more than 1!---"<<nbases<<endl;
    }
    cout<<"Get base number"<<endl;
    index_base=1;
/* we know there is only one zone (real working code would check!) */
    index_zone=1;
    int nzones=-1;
    cg_nzones(index_file,index_base,&nzones);
    if(nzones!=1)
    {
        cout<<"The zones is more than 1!---"<<nzones<<endl;
    }
/* get zone size (and name - although not needed here) */
    cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);
    cout<<"read Zone"<<endl;
/* lower range index */
    irmin=1;
/* upper range index of vertices */
    irmax=isize[0][0];
    cout<<"irmax = "<< irmax <<endl;

    int nNodes=irmax;

    cgsize_t nTetras=isize[1][0];
    cgsize_t *tetras=new cgsize_t[nTetras*4];

/* read grid coordinates */
    cout<<"Starting to read grid coordinates"<<endl;
    int ncoord=-1;
    cg_ncoords(index_file,index_base,index_zone,&ncoord);

    float *x=new float[irmax];
    float *y=new float[irmax];
    float *z=new float[irmax];

    cout<<"Read grid coordinates"<<endl;
    cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,x);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,y);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,z);

    cout<<"The number of vertices is: "<<irmax<<" "<<ncoord<<endl;
/* find out how many sections */
    cg_nsections(index_file,index_base,index_zone,&nsections);
    printf("\nnumber of sections=%i\n",nsections);
/* read element connectivity */

    int nTriSum=0;

    //cgsize_t *triSection[3];
    cgsize_t **triSection=nullptr;
    triSection=new cgsize_t*[nsections-1];
    int *nTriSection=new int[nsections-1];
    //vector<vector<cgsize_t>>triSection;
    int count=0;

    for (index_sect=1; index_sect <= nsections; index_sect++)
    {
      cg_section_read(index_file,index_base,index_zone,index_sect,sectionname,
                      &itype,&istart,&iend,&nbndry,&iparent_flag);


      printf("\nReading section data...\n");
      printf("   section name=%s\n",sectionname);
      printf("   section type=%s\n",ElementTypeName[itype]);
      printf("   istart,iend=%i, %i\n",(int)istart,(int)iend);
      if (itype == CGNS_ENUMV(TETRA_4))
      {
        printf("   reading element data for this Tetras\n");
        cg_elements_read(index_file,index_base,index_zone,index_sect,tetras, \
                         &iparentdata);
      }

      if (itype == CGNS_ENUMV(TRI_3))
      {
          int nTri=iend-istart+1;
          nTriSection[count]=nTri;
          nTriSum+=nTri;
          triSection[count]=new cgsize_t[nTri*3];
          string strName=sectionname;
          bcstring.push_back(strName);
          printf("   reading element data for this element\n");
          cg_elements_read(index_file,index_base,index_zone,index_sect,triSection[count], \
                  &iparentdata);
          count++;
      }

    }
/* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully read unstructured grid from file grid_c.cgns\n");

    mesh.NumNodes=nNodes;
    mesh.nodes=new Node[nNodes];
    mesh.NumTetras=nTetras;
    mesh.pTetras=new TETRAS[nTetras];
    mesh.NumTris=nTriSum;
    mesh.pTris=new TRI[nTriSum];

    for(int i=0;i<nNodes;i++)
    {
        mesh.nodes[i].index=i;
        mesh.nodes[i].coord.x=x[i];
        mesh.nodes[i].coord.y=y[i];
        mesh.nodes[i].coord.z=z[i];
    }
    for(int i=0;i<nTetras;i++)
    {
        mesh.pTetras[i].index=i;
        mesh.pTetras[i].vertices[0]=tetras[i*4+0]-1;
        mesh.pTetras[i].vertices[1]=tetras[i*4+1]-1;
        mesh.pTetras[i].vertices[2]=tetras[i*4+2]-1;
        mesh.pTetras[i].vertices[3]=tetras[i*4+3]-1;
    }

    count=0;
   // mesh.pTris[i].index=i;
    for(int j=0;j<nsections-1;j++)
    {
        for(int k=0;k<nTriSection[j];k++)
        {
            mesh.pTris[count].index=count;
            mesh.pTris[count].vertices[0]=triSection[j][k*3+0]-1;
            mesh.pTris[count].vertices[1]=triSection[j][k*3+1]-1;
            mesh.pTris[count].vertices[2]=triSection[j][k*3+2]-1;
            mesh.pTris[count].iSurf=j;
            count++;
        }
    }

    assert(count==nTriSum);


    return 1;
}
#endif
