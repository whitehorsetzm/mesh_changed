#include "openfoamfile.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cassert>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>

#include <mpi.h>

using namespace std;


OpenFoamFile::OpenFoamFile(int iRank, const HYBRID_MESH *pMesh, idx_t nFace, idx_t nCell, const BoundaryCondition *bc, string path, string name, int nProc, idx_t nAllInFace)
    :pHybridMesh(pMesh), nFaceBeg(nFace), nCellBeg(nCell), pBC(bc), nAllInternalFace(nAllInFace)
{
    filename = name;
    homepath = path + "/" + name;
    nProcessor = nProc;
    iProc = iRank;

    vecHead.push_back("/*--------------------------------*- C++ -*----------------------------------*\\");
    vecHead.push_back("| =========                 |                                                 |");
    vecHead.push_back("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |");
    vecHead.push_back("|  \\\\    /   O peration     | Version:  xxxxxxx                               |");
    vecHead.push_back("|   \\\\  /    A nd           | Revision: exported                              |");
    vecHead.push_back("|    \\\\/     M anipulation  | Web:      http://www.OpenFOAM.org               |");
    vecHead.push_back("\\*---------------------------------------------------------------------------*/");
    vecHead.push_back("FoamFile");
    vecHead.push_back("{");
    vecHead.push_back("    version     2.0;");
    vecHead.push_back("    format      ascii;");

    lablist      = "    class       labelList;";
    vecfield     = "    class       vectorField;";
    bdymesh      = "    class       polyBoundaryMesh;";
    facelist     = "    class       faceList;";
    location     = "    location    \"constant/polyMesh\";";
    objectPoint  = "    object      points;";
    objectPtPAd  = "    object      pointProcAddressing;";
    objectBdy    = "    object      boundary;";
    objectBdyAd  = "    object      boundaryProcAddressing;";
    objectCellAd = "    object      cellProcAddressing;";
    objectFaceAd = "    object      faceProcAddressing;";
    objectFace   = "    object      faces;";
    objectNeib   = "    object      neighbour;";
    objectOwner  = "    object      owner;";

}

#if 1
/************************************************************************/
/* 创建OpenFoam文件格式输出的文件夹树，并初始化各文件名                       */
/************************************************************************/
int OpenFoamFile::creatFile()
{
    int i = 0, i_err = 0;
    string path;
    string pathProc;

    /* 创建根目录 */
    path = homepath + "/";
    i_err = createDictionary( path );

    if(i_err == 0)
	{
	 cout<< "Create home directionay in Proc " << iProc << "successfuly!"<<endl;
	}
#ifdef DEBUG
    cout<<path<<":"<<i_err<<endl;
#endif DEBUG

    /* 创建总体网格0目录 */
//    path = homepath + "\\0";
//    i_err = createDictionary( path );
//    cout<<path<<":"<<i_err<<endl;
//    field = path;

    /* 创建总体网格system目录 */
//    path = homepath + "\\system";
//    i_err = createDictionary( path );
//    cout<<path<<":"<<i_err<<endl;
//    system = path;

    /* 创建总体网格constant\polyMesh目录
       并初始化总体网格文件名
    */
	if (nProcessor == 1)
	{
		path = homepath + "/constant/polyMesh";
		i_err = createDictionary( path );
		if(i_err != 0)
			cout<<path<<":"<<i_err<<endl;
#ifdef DEBUG
		cout<<path<<":"<<i_err<<endl;
#endif DEBUG
		points = path + "/points";
		faces = path + "/faces";
		neighbour = path + "/neighbour";
		owner = path + "/owner";
		boundary = path + "/boundary";
	}
	else
	{
		/* 创建分块网格目录，并初始化分块网格文件名 */
		pathProc = homepath + "/processor" + intToString(iProc) + "/";
		i_err = createDictionary( pathProc );
		if(i_err != 0)
			cout<<pathProc<<":"<<i_err<<endl;
		
#ifdef DEBUG
		cout<<pathProc<<":"<<i_err<<endl;
#endif DEBUG

		path = pathProc + "/0/";
		i_err = createDictionary( path );
#ifdef DEBUG
		cout<<path<<":"<<i_err<<endl;
#endif DEBUG
		fieldProc = path;

		path = pathProc + "/constant//polyMesh/";
		i_err = createDictionary( path );
#ifdef DEBUG
		cout<<path<<":"<<i_err<<endl;
#endif DEBUG

		pointsProc     = path + "/points";
		pointAddreProc = path + "/pointProcAddressing";
		facesProc      = path + "/faces";
		faceAddreProc  = path + "/faceProcAddressing";
		neighbourProc  = path + "/neighbour";
		ownerProc      = path + "/owner";
		cellAddreProc  = path + "/cellProcAddressing";
		boundaryProc   = path + "/boundary";
		bdyAddreProc   = path + "/boundaryProcAddressing";
	}

    return 0;
}
#else
int OpenFoamFile::creatFile()
{
    int i = 0, i_err = 0;
    string path;
    string pathProc;

    /* 创建根目录 */
    path = homepath;
    i_err = createDictionary( path );
    cout<<path<<":"<<i_err<<endl;

    /* 创建总体网格0目录 */
    path = homepath + "\\0";
    i_err = createDictionary( path );
    cout<<path<<":"<<i_err<<endl;
    field = path;

    /* 创建总体网格system目录 */
    path = homepath + "\\system";
    i_err = createDictionary( path );
    cout<<path<<":"<<i_err<<endl;
    system = path;

    /* 创建总体网格constant\polyMesh目录
       并初始化总体网格文件名
    */
    path = homepath + "\\constant\\polyMesh";
    i_err = createDictionary( path );
    cout<<path<<":"<<i_err<<endl;

    points = path + "\\points";
    faces = path + "\\faces";
    neighbour = path + "\\neighbour";
    owner = path + "\\owner";
    boundary = path + "\\boundary";

    /* 创建分块网格目录，并初始化分块网格文件名 */
    pointsProc     = new string[nProcessor];
    pointAddreProc = new string[nProcessor];
    facesProc      = new string[nProcessor];
    faceAddreProc  = new string[nProcessor];
    neighbourProc  = new string[nProcessor];
    ownerProc      = new string[nProcessor];
    cellAddreProc  = new string[nProcessor];
    boundaryProc   = new string[nProcessor];
    bdyAddreProc   = new string[nProcessor];
    fieldProc      = new string[nProcessor];

    for (i = 0; i < nProcessor; i++)
    {
        string str_num;
        stringstream stream;
        stream << i;
        stream >> str_num;
        pathProc = homepath + "\\processor" + str_num;
        i_err = createDictionary( pathProc );
        cout<<pathProc<<":"<<i_err<<endl;

        path = pathProc + "\\0";
        i_err = createDictionary( path );
        cout<<path<<":"<<i_err<<endl;
        fieldProc[i] = path;

        path = pathProc + "\\constant\\polyMesh";
        i_err = createDictionary( path );
        cout<<path<<":"<<i_err<<endl;

        pointsProc[i]     = path + "\\points";
        pointAddreProc[i] = path + "\\pointProcAddressing";
        facesProc[i]      = path + "\\faces";
        faceAddreProc[i]  = path + "\\faceProcAddressing";
        neighbourProc[i]  = path + "\\neighbour";
        ownerProc[i]      = path + "\\owner";
        cellAddreProc[i]  = path + "\\cellProcAddressing";
        boundaryProc[i]   = path + "\\boundary";
        bdyAddreProc[i]   = path + "\\boundaryProcAddressing";
    }
    return 0;
}
#endif

int OpenFoamFile::writeMesh_OpenFoam_Points()
{
    int i;
    int nPoints = pHybridMesh->NumNodes;
	string fname;
	if (nProcessor == 1)
		fname = points;
	else if (nProcessor > 1)
		fname = pointsProc;

    FILE *fp = fopen(fname.c_str(), "w");
    if (!fp)
    {
        printf("Cannot open file %s.\n", fname.c_str());
        return 1;	/* E_INVALID_ARG */
    }

    writePointsHead(fp);

    // 输出分区网格节点数目
    fprintf(fp,  "%d\n(\n", nPoints);

    // 网格节点坐标
    for(i = 0; i < nPoints; ++i)
    {
        fprintf(fp, "(%f %f %f)\n", pHybridMesh->nodes[i].coord.x,
                pHybridMesh->nodes[i].coord.y, pHybridMesh->nodes[i].coord.z);
    }

    fprintf(fp, ")");

    writeTail(fp);

    if (fp)
    {
        fclose(fp);
    }

    return 0;

}

/* 输出文件points头信息 */
int OpenFoamFile::writePointsHead(FILE *fp)
{
    int i;
    if (!fp) return 1;

    for (i = 0; i < vecHead.size() ; i++)
    {
        fprintf(fp,"%s\n",vecHead[i].c_str());
    }
    fprintf(fp,"%s\n",vecfield.c_str());
    fprintf(fp,"%s\n",location.c_str());
    fprintf(fp,"%s\n}\n",objectPoint.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");

    return 0;
}

/* 输出分区网格pointProcAddressing文件 */
int OpenFoamFile::writeProcMesh_OpenFoam_pointProcAddressing()
{
    int i;
    int nPoints = pHybridMesh->NumNodes;
    FILE *fp = fopen(pointAddreProc.c_str(), "w");

    if (!fp)
    {
        printf("Cannot open file %s.\n", pointAddreProc.c_str());
        return 1;	/* E_INVALID_ARG */
    }

    writePointPAdHead(fp);

    fprintf(fp, "%d\n(\n", nPoints);

    /* 输出point全局编号 */
    for(i = 0; i < nPoints; ++i)
    {
        fprintf(fp,"%lld\n", pHybridMesh->nodes[i].index);
    }

    fprintf(fp,")\n");

    writeTail(fp);

    if (fp)
    {
        fclose(fp);
    }
    return 0;
}

/* 输出文件pointProcAddressing头信息 */
int OpenFoamFile::writePointPAdHead(FILE *fp)
{
    int i;
    if (!fp)	return 1;
    for (i = 0; i < vecHead.size(); ++i)
    {
        fprintf(fp,"%s\n",vecHead[i].c_str());
    }
    fprintf(fp,"%s\n",lablist.c_str());
    fprintf(fp,"%s\n",location.c_str());
    fprintf(fp,"%s\n}\n",objectPtPAd.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
    return 0;
}

/* 输出分区网格面片相关文件,包括faces, faceProcAddressing, owner, neighbor */
int OpenFoamFile::writeMesh_OpenFoam_AllFaceFiles()
{
	int rtn = 0;
    int i, j, iCell, iOwner, iFace, iNeig, iBdry, iNode;
	int nPoints = 0;
	int nTetr, nPrsm, nPyra, nHexa, nElem;
	int nTria, nQuad;
	int nBdyFaces, nInternalFaces, nFaces; // 边界面片数, 内部面片数, 总面片数 

	TETRAS *pTetras;
	HEX *pHexes;

	TRI *pTri;
	QUAD *pQuad;

	Node *pNodes = pHybridMesh->nodes;

    const int faceTriaNodes[][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}}; // 前三个点右手法则背离第四个点
    //const int faceTriaNodes[][3] = {{1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};  // 前三个点右手法则指向第四个点

	FILE *fileFaces = nullptr, *fileOwner = nullptr, *fileNeigh = nullptr, *fileFaceAd =nullptr;	

    double pntCoord[4][3];

	nPoints = pHybridMesh->NumNodes;

	nTetr = pHybridMesh->NumTetras;
	nPrsm = pHybridMesh->NumPrsm;
	nPyra = pHybridMesh->NumPyra;
	nHexa = pHybridMesh->NumHexes;
	nElem = nTetr + nPrsm + nPyra + nHexa;

	nTria = pHybridMesh->NumTris;
	nQuad = pHybridMesh->NumQuads;

	nBdyFaces = nTria + nQuad;
	nInternalFaces = (4*nTetr + 5*nPrsm + 5*nPyra + 6*nHexa - nBdyFaces)/2;
	nFaces = nBdyFaces + nInternalFaces;

	pTetras = pHybridMesh->pTetras;
	pHexes = pHybridMesh->pHexes;

	pTri = pHybridMesh->pTris;
	pQuad = pHybridMesh->pQuads;

	string fnameFaces, fnameOwner, fnameNeigh;
	if (nProcessor == 1)
	{
		fnameFaces = faces;
		fnameOwner = owner;
		fnameNeigh = neighbour;
	}
	else if (nProcessor > 1)
	{
		fnameFaces = facesProc;
		fnameOwner = ownerProc;
		fnameNeigh = neighbourProc;
	}
	fileFaces = fopen(fnameFaces.c_str(), "w");
	fileOwner = fopen(fnameOwner.c_str(), "w");
	fileNeigh = fopen(fnameNeigh.c_str(), "w");
	if (nProcessor > 1)
	{
		fileFaceAd = fopen(faceAddreProc.c_str(), "w");
	}
	if (!fileFaces)
	{
		printf("Cannot open file %s.\n", fnameFaces.c_str());
		rtn = 1;/* E_INVALID_ARG */
	}
	if (!fileOwner)
	{
		printf("Cannot open file %s.\n", fnameOwner.c_str());
		rtn = 1;/* E_INVALID_ARG */
	}
	if (!fileNeigh)
	{
		printf("Cannot open file %s.\n", fnameNeigh.c_str());
		rtn = 1;/* E_INVALID_ARG */
	}
	if (nProcessor > 1)
	{
		if (!fileFaceAd)
		{
			printf("Cannot open file %s.\n", faceAddreProc.c_str());
			rtn = 1;/* E_INVALID_ARG */
		}
	}
	if (rtn == 1)
	{
		if (fileFaces)
			fclose(fileFaces);
		if (fileOwner)
			fclose(fileOwner);
		if (fileNeigh)
			fclose(fileNeigh);
		if (nProcessor > 1)
		{		
			if (fileFaceAd) fclose(fileFaceAd);
		}
		return rtn;
	}

	writeFacesHead(fileFaces);
	writeOwnerHead(fileOwner, nPoints, nElem, nFaces, nInternalFaces);
	writeNeighbHead(fileNeigh, nPoints, nElem, nFaces, nInternalFaces);
	fprintf(fileFaces,"%d\n(\n", nFaces);
	fprintf(fileOwner,"%d\n(\n",nFaces);
	fprintf(fileNeigh,"%d\n(\n",nInternalFaces);
	if (nProcessor > 1)
	{
		writeFacePAdHead(fileFaceAd);
		fprintf(fileFaceAd,"%d\n(\n",nFaces);
	}

	/************************************************************************/
	/* 先统计本分区起始全局面片编号                                         */
	/************************************************************************/
//	int nFaceTotal = pHybridMesh->NumUniqueFacets;
//	int nFaceBeg;
//	int *nInternalFacesProc;
//	nInternalFacesProc = new int[nProcessor](); // 记录每个分区内部面片数目,初始化为0

//	MPI_Allgather(&nInternalFaces, 1, MPI_INT,
//		nInternalFacesProc, 1, MPI_INT, MPI_COMM_WORLD);

//	nFaceBeg = nFaceTotal;
//	for (i = 0; i < iProc; ++i)
//	{
//		nFaceBeg += nInternalFacesProc[i];
//	}

    /*
	内部面：这些面连接着两个单元（不会多于两个单元）。对于每一个内部面，
		面的法线方向朝向单元号较大的单元内部，比如，对于 2 单元和 5 单元来说，法
		线方向指向 5 单元；
	边界面：自从它们被分配在边界域上时， 这些面就属于一个单元。因此一个
		边界面只属于一个单元和边界类别。这些面的法向指向计算域外面。
	*/

	/* 先输出内部面 */
#ifdef DEBUG
	cout << "Proc " << iProc << " Write Internal Faces Begin" << endl;

    for(iCell = 0; iCell < nTetr; ++iCell)
    {
        for (j = 0; j < 4; ++j)
        {
            pntCoord[j][0] = pNodes[pTetras[iCell].vertices[j]].coord.x;
            pntCoord[j][1] = pNodes[pTetras[iCell].vertices[j]].coord.y;
            pntCoord[j][2] = pNodes[pTetras[iCell].vertices[j]].coord.z;
        }
        if (volumeTetrSign(pntCoord[0], pntCoord[1], pntCoord[2], pntCoord[3]) == false)
        {
            cout << "iCell = " << iCell << " is wrong direction!" << endl;
        }
        else
        {
            //cout << "iCell = " << iCell << " is right direction!" << endl;
        }
    }
#endif // DEBUG

	int idx = 0;
    map<int, int> mapNeibLocN;
    map<int, int>::iterator itMapIntInt;
	/* 四面体 */    
	/* pHybridMesh中四面体节点的顺序：前三个点右手法则指向第四个点 */
	for (iCell = 0; iCell < nTetr; ++iCell)
	{
        mapNeibLocN.clear();
		for (i = 0; i < 4; ++i)
		{
			iNeig = pTetras[iCell].neighbors[i]; 
			if (iNeig >= 0) // 邻居单元不为负的即为内部面
			{
				if (iCell < iNeig) // 面片右手法则指向编号大的单元
				{
                    mapNeibLocN.insert(make_pair(iNeig,i));
				}
			}
		}

        for(itMapIntInt = mapNeibLocN.begin(); itMapIntInt != mapNeibLocN.end(); ++itMapIntInt)
        {
            i = itMapIntInt->second;
            iNeig = pTetras[iCell].neighbors[i];
            // 输出内部三角形面片节点的局部编号
            fprintf(fileFaces,"3(%d %d %d)\n", pTetras[iCell].vertices[faceTriaNodes[i][0]],
                pTetras[iCell].vertices[faceTriaNodes[i][1]], pTetras[iCell].vertices[faceTriaNodes[i][2]]);
            // 输出内部三角面片的全局编号
            if (nProcessor > 1)
            {
                fprintf(fileFaceAd, "%lld\n", nFaceBeg+idx+1); // * 注意：从1开始
                idx++;
            }

            // 输出内部三角面片的所属单元的局部编号
            iOwner = iCell;
            fprintf(fileOwner,"%d\n",iOwner);
            // 输出内部三角面片的邻居单元的局部编号
            fprintf(fileNeigh,"%d\n",iNeig);
        }
	}
	/* 三棱柱 */
	for (i = 0; i < nPrsm; ++i)
	{
		// 
	}
	/* 金字塔 */
	for (i = 0; i < nPyra; ++i)
	{
		// 
	}
	/* 六面体 */
	for (i = 0; i < nHexa; ++i)
	{
		// 
	}
#ifdef DEBUG
	cout << "Proc " << iProc << " Write Internal Faces Finished" << endl;
#endif // DEBUG
	/************************************************************************/
	/*  输出边界面                                                          */
	/************************************************************************/
#ifdef DEBUG
	cout << "Proc " << iProc << " Write Boundary Faces Begin" << endl;
#endif // DEBUG
	vector<int>::iterator result;
	vector<int> vecSurfId;
	vector<int> vecFaceId;
	vector<vector<int> > vecBCFaceIds;  // 用于记录每个边界条件中包含的边界面片编号
	multimap<string,int> mapBCFaceIds;
	string sProcBCName;
	int nBcTria, nBcQuad; 
	int nBcFace;
	TRI tri;
	QUAD quad;
	TETRAS tetr;

	int iTriPnt[3];

	/************************************************************************/
	/* 遍历所有边界面片
	   按照边界条件顺序以及处理器交界面边界条件顺序
	   将面片编号放入指定容器 */
	/************************************************************************/
	vecBCFaceIds.resize(pBC->bcInfo.vecBC.size());
	/* 三角形 */
	for (iFace = 0; iFace < nTria; ++iFace)
	{	
		i = 0; 
		vecSurfId.clear();
		result = vecSurfId.end();
		while (result == vecSurfId.end() && i < pBC->bcInfo.vecBC.size())
		{
			vecSurfId = pBC->bcInfo.vecBC[i].faceIDs;
			result = find(vecSurfId.begin(), vecSurfId.end(), pTri[iFace].iSurf+1); //边界条件中几何面片编号从1开始， pTri中的几何面片编号从0开始
			if (result == vecSurfId.end())
				++i;
		}

		if (result != vecSurfId.end()) // 找到对应的物理边界条件
		{
			vecBCFaceIds[i].push_back(iFace);	
		}
		else
		{
			assert(pTri[iFace].iOppoProc >= 0 && pTri[iFace].iOppoProc != iProc);
			sProcBCName = "procBoundary" + intToString(iProc) + "to" + intToString(pTri[iFace].iOppoProc);
			mapBCFaceIds.insert(map<string,int>::value_type(sProcBCName,iFace));
		}
	}

	/* 四边形 */
	for (iFace = 0; iFace < nQuad; ++iFace)
	{	
		i = 0; 
		vecSurfId.clear();
		while (result == vecSurfId.end() && i < pBC->bcInfo.vecBC.size())
		{
			vecSurfId = pBC->bcInfo.vecBC[i++].faceIDs;
			result = find(vecSurfId.begin(), vecSurfId.end(), pQuad[iFace].iSurf); 
		}

		if (result != vecSurfId.end()) // 找到对应的物理边界条件
		{
			vecBCFaceIds[i].push_back(iFace+nTria); // 注意,这里为了区分四边形的编号加了三角形边界面片数目		
		}
		else
		{
			assert(pQuad[iFace].iOppoProc >= 0 && pQuad[iFace].iOppoProc != iProc);
			sProcBCName = "procBoundary" + intToString(iProc) + "to" + intToString(pQuad[iFace].iOppoProc);
			mapBCFaceIds.insert(map<string,int>::value_type(sProcBCName,iFace+nTria));// 注意,这里为了区分四边形的编号加了三角形边界面片数目
		}
	}

	/* 输出物理边界面片 */
	nBcFace = nInternalFaces;
	int nCount; 
	BoundaryOF bcOF;
	for (i = 0; i < vecBCFaceIds.size(); ++i)
	{
		/* 记录此边界条件信息,用于输出boundary文件 */		
		bcOF.myProcNo = iProc;
		bcOF.name = pBC->bcInfo.vecBC[i].name;
		bcOF.type = "patch";
		bcOF.startFace = nBcFace;
		
		nCount = 0;
		for (iFace = 0; iFace < vecBCFaceIds[i].size(); ++iFace)
		{
			if (iFace < nTria)
			{
				tri = pTri[vecBCFaceIds[i][iFace]];
                assert(tri.iCell < nTetr);
#ifdef DEBUG
                //cout << "iFace = " << iFace << endl;
				iTriPnt[0] = tri.vertices[0];
				iTriPnt[1] = tri.vertices[1];
				iTriPnt[2] = tri.vertices[2];
				// 先检查面片方向（要求右手法则法向指向区域外部）

				tetr = pTetras[tri.iCell];
				for (j = 0; j < 4;)
				{
					if (findIntNum(tetr.vertices[j], iTriPnt, 3) == false)
						break;
					++j;
				}
				assert(j < 4);
				pntCoord[3][0]  = pNodes[tetr.vertices[j]].coord.x;
				pntCoord[3][1]  = pNodes[tetr.vertices[j]].coord.y;
				pntCoord[3][2]  = pNodes[tetr.vertices[j]].coord.z;
				for (j = 0; j < 3; ++j)
				{
					pntCoord[j][0] = pNodes[tri.vertices[0]].coord.x;
					pntCoord[j][1] = pNodes[tri.vertices[0]].coord.y;
					pntCoord[j][2] = pNodes[tri.vertices[0]].coord.z;
				}	
				if (volumeTetrSign(pntCoord[0], pntCoord[1], pntCoord[2], pntCoord[3]))
				{
					fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[0],
						tri.vertices[2], tri.vertices[1]);
				}
				else
				{
                    cout << "boundary face" << endl;
					fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[0],
						tri.vertices[1], tri.vertices[2]);
				}
#else
                fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[0],
                    tri.vertices[2], tri.vertices[1]);
#endif //DEBUG

				if (nProcessor > 1)
				{
#ifdef _INTERNAL_FACE_FIRST_
#ifdef _SURFACE_FACE_FIRST_
                    fprintf(fileFaceAd,"%lld\n", nAllInternalFace + tri.index+1); // phycail boundary face gIdx < interface face gIdx
#else
                    fprintf(fileFaceAd,"%lld\n", nAllInternalFace + tri.index + pHybridMesh->NumUniqueInterfFacets +1); // phycail boundary face gIdx > interface face gIdx
#endif //_SURFACE_FACE_FIRST_
#else
                    fprintf(fileFaceAd,"%lld\n", tri.index+1); // * 注意：从1开始
#endif //_INTERNAL_FACE_FIRST_
				}				

				fprintf(fileOwner, "%d\n", tri.iCell);
				nCount++;
			}
			else
			{
				quad = pQuad[vecBCFaceIds[i][iFace-nTria]];

				// 先检查面片方向（要求右手法则法向指向区域外部）
				printf("Check face direction! (Need to be done)\n");
				fprintf(fileFaces, "4(%d %d %d %d)\n", quad.vertices[0],
					quad.vertices[1], quad.vertices[2], quad.vertices[3]);

				if (nProcessor > 1)
				{
#ifdef _INTERNAL_FACE_FIRST_
#ifdef _SURFACE_FACE_FIRST_
                    fprintf(fileFaceAd,"%lld\n", nAllInternalFace + quad.index+1); // phycail boundary face gIdx < interface face gIdx
#else
                    fprintf(fileFaceAd,"%lld\n", nAllInternalFace + quad.index + pHybridMesh->NumUniqueInterfFacets +1); // phycail boundary face gIdx > interface face gIdx
#endif //_SURFACE_FACE_FIRST_
#else
                    fprintf(fileFaceAd, "%lld\n", quad.index+1); // * 注意：从1开始
#endif //_INTERNAL_FACE_FIRST_
				}				

				fprintf(fileOwner, "%d\n", quad.iCell);
				nCount++;
			}
		}
		bcOF.nFace = nCount;		
		vecBdyOF.push_back(bcOF);
		nBcFace += nCount;
	}
#ifdef DEBUG
	cout << "Proc " << iProc << " Write Boundary Faces Finished" << endl;
	if (fileFaces)
		fflush(fileFaces);
	if (fileOwner)
		fflush(fileOwner);
	if (fileNeigh)
		fflush(fileNeigh);
	if (nProcessor > 1)
	{		
		if (fileFaceAd) fflush(fileFaceAd);
	}
#endif // DEBUG
	/* 输出交界面面片 */
#ifdef DEBUG
	cout << "Proc " << iProc << " Write Interface Faces Begin" << endl;
#endif // DEBUG
	multimap<string, int>::size_type sizeCount = 0;
	multimap<string, int>::iterator itMap;
	itMap = mapBCFaceIds.begin();
	while (itMap != mapBCFaceIds.end())
	{
		sProcBCName = itMap->first;
		sizeCount = mapBCFaceIds.count(sProcBCName);

		bcOF.myProcNo = iProc;
		bcOF.name = sProcBCName;
		bcOF.type = "processor"; // 交界面边界条件类型
		bcOF.startFace = nBcFace;
        bcOF.matchTolerance =  0.0001;
		
		nCount = 0;

        /* sort write order by global face index */

        /* the global index of the first node of the face must match the coupled equvalent in neibghbor part mesh
            if not, there will be error like:
            **Error in coupled point location: 14 faces have their 0th or consecutive vertex not opposite their coupled equivalent. Average mismatch 0.*/
        map<int, int> mapWriteTriaOrder;
        map<int, int> mapWriteQuadOrder;
        int gIdx[4];
        for (multimap<string, int>::size_type cnt = 0; cnt != sizeCount; ++cnt,++itMap)
        {
            iFace = itMap->second;
            if(iFace < nTria)
            {
                tri = pTri[iFace];
                mapWriteTriaOrder.insert(make_pair(tri.index,iFace));
            }
            else
            {
                // quad
            }
        }
        /* write triangle */
        for (itMapIntInt = mapWriteTriaOrder.begin(); itMapIntInt != mapWriteTriaOrder.end(); ++itMapIntInt)
        {
            iFace = itMapIntInt->second;
            tri = pTri[iFace];
            iNode = 0;
            for(i = 0; i < 3; ++i)
            {
                gIdx[i] = pNodes[tri.vertices[i]].index;
            }

            iNode = findIdxMinInt(gIdx, 3);
            assert(iNode >= 0 && iNode <= 2);
            fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[iNode],
                tri.vertices[(iNode+2)%3], tri.vertices[(iNode+1)%3]);

            if (nProcessor > 1)
            {
                if (tri.iOppoProc > iProc) // 表示面方向指向邻居分区中的单元
                {
#ifdef _INTERNAL_FACE_FIRST_
#ifdef _SURFACE_FACE_FIRST_
                    fprintf(fileFaceAd,"%lld\n", nAllInternalFace + tri.index+1); // phycail boundary face gIdx < interface face gIdx
#else
                    fprintf(fileFaceAd,"%lld\n", nAllInternalFace + tri.index - pHybridMesh->NumUniqueSurfFacets +1); // phycail boundary face gIdx > interface face gIdx
#endif //_SURFACE_FACE_FIRST_
#else
                    fprintf(fileFaceAd,"%lld\n", tri.index+1); // * 注意：从1开始
#endif //_INTERNAL_FACE_FIRST_
                }
                else
                {
#ifdef _INTERNAL_FACE_FIRST_
#ifdef _SURFACE_FACE_FIRST_
                    fprintf(fileFaceAd,"%lld\n", -1*(nAllInternalFace + tri.index+1)); // phycail boundary face gIdx < interface face gIdx
#else
                    fprintf(fileFaceAd,"%lld\n", -1*(nAllInternalFace + tri.index - pHybridMesh->NumUniqueSurfFacets +1)); // phycail boundary face gIdx > interface face gIdx
#endif //_SURFACE_FACE_FIRST_
#else
                    fprintf(fileFaceAd,"%lld\n", -1*(tri.index+1)); // * 注意：从1开始
#endif //_INTERNAL_FACE_FIRST_
                }
            }

            fprintf(fileOwner, "%d\n", tri.iCell);
            nCount++;
            iNeig = tri.iOppoProc;
        }
        /* write quad */
        for (itMapIntInt = mapWriteQuadOrder.begin(); itMapIntInt != mapWriteQuadOrder.end(); ++itMapIntInt)
        {
            //				fprintf(fileOwner, "%d\n", quad.iCell);
            //				nCount++;
            //				iNeig = quad.iOppoProc;
        }


//		for (multimap<string, int>::size_type cnt = 0; cnt != sizeCount; ++cnt,++itMap)
//		{
//			iFace = itMap->second;
//			if (iFace < nTria)
//			{
//				tri = pTri[iFace];

//#ifdef CHECKMESH
//				iTriPnt[0] = tri.vertices[0];
//				iTriPnt[1] = tri.vertices[1];
//				iTriPnt[2] = tri.vertices[2];
//				// 先检查面片方向（要求右手法则法向指向区域外部）
//				tetr = pTetras[tri.iCell];
//				for (j = 0; j < 4;)
//				{
//					if (findIntNum(tetr.vertices[j], iTriPnt, 3) == false)
//						break;
//					++j;
//				}
//				assert(j < 4);
//				pntCoord[3][0]  = pNodes[tetr.vertices[j]].coord.x;
//				pntCoord[3][1]  = pNodes[tetr.vertices[j]].coord.y;
//				pntCoord[3][2]  = pNodes[tetr.vertices[j]].coord.z;
//				for (j = 0; j < 3; ++j)
//				{
//					pntCoord[j][0] = pNodes[tri.vertices[0]].coord.x;
//					pntCoord[j][1] = pNodes[tri.vertices[0]].coord.y;
//					pntCoord[j][2] = pNodes[tri.vertices[0]].coord.z;
//				}
//				if (volumeTetrSign(pntCoord[0], pntCoord[1], pntCoord[2], pntCoord[3]))
//				{
//					fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[0],
//						tri.vertices[2], tri.vertices[1]);
//				}
//				else
//				{
//					fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[0],
//						tri.vertices[1], tri.vertices[2]);
//				}
//#else
//                fprintf(fileFaces, "3(%d %d %d)\n", tri.vertices[0],
//                    tri.vertices[2], tri.vertices[1]);
//#endif

//				if (nProcessor > 1)
//				{
//					if (tri.iOppoProc > iProc) // 表示面方向指向邻居分区中的单元
//					{
//						fprintf(fileFaceAd,"%d\n", tri.index+1); // * 注意：从1开始
//					}
//					else
//					{
//						fprintf(fileFaceAd,"%d\n", -1*(tri.index+1)); // * 注意：从1开始
//					}
//				}

//				fprintf(fileOwner, "%d\n", tri.iCell);
//				nCount++;
//				iNeig = tri.iOppoProc;
//			}
//			else
//			{
//				quad = pQuad[iFace-nTria];
//				// 先检查面片方向（要求右手法则法向指向区域外部）
//				printf("Check face direction! (Need to be done)\n");
//				fprintf(fileFaces, "4(%d %d %d %d)\n", quad.vertices[0],
//					quad.vertices[1], quad.vertices[2], quad.vertices[3]);

//				if (nProcessor > 1)
//				{
//					if (tri.iOppoProc > iProc) // 表示面方向指向邻居分区中的单元
//					{
//						fprintf(fileFaceAd,"%d\n", quad.index+1); // * 注意：从1开始
//					}
//					else
//					{
//						fprintf(fileFaceAd,"%d\n", -1*(quad.index+1)); // * 注意：从1开始
//					}
//				}

//				fprintf(fileOwner, "%d\n", quad.iCell);
//				nCount++;
//				iNeig = quad.iOppoProc;
//			}
//		}
		nBcFace += nCount;
		bcOF.nFace = nCount;
		bcOF.myProcNo = iProc;
		bcOF.neighbProcNo = iNeig;
		vecBdyOF.push_back(bcOF);
    }


#ifdef DEBUG
	cout << "Proc " << iProc << " Write Interface Faces Finished" << endl;
	if (fileFaces)
		fflush(fileFaces);
	if (fileOwner)
		fflush(fileOwner);
	if (fileNeigh)
		fflush(fileNeigh);
	if (nProcessor > 1)
	{		
		if (fileFaceAd) fflush(fileFaceAd);
	}
#endif // DEBUG


	//nBcFace = nInternalFaces;
	//for (i = 0; i < pBC->bcInfo.vecBC.size(); ++i)
	//{
	//	vecSurfId = pBC->bcInfo.vecBC[i].faceIDs;
	//	nBcTria = 0;
	//	for (iFace = 0; iFace < nTria; ++iFace)
	//	{
	//		result = find(vecSurfId.begin(), vecSurfId.end(), pTri[iFace].iSurf); 
	//		if (result != vecSurfId.end()) // 找到相应的几何面片
	//		{
	//			vecFaceId.push_back(iFace);
	//			nBcTria++;
	//		}
	//	}
	//	nBcQuad = 0;
	//	for (iFace = 0; iFace < nQuad; ++iFace)
	//	{
	//		result = find(vecSurfId.begin(), vecSurfId.end(), pQuad[iFace].iSurf);
	//		if (result != vecSurfId.end())
	//		{
	//			vecFaceId.push_back(iFace);
	//			nBcQuad++;
	//		}
	//	}
	//	/* 记录此边界条件信息,用于输出boundary文件 */
	//	BoundaryOF bcOF;
	//	bcOF.myProcNo = iProc;
	//	bcOF.name = pBC->bcInfo.vecBC[i].name;
	//	bcOF.type = "patch";
	//	bcOF.startFace = nBcFace;
	//	nBcFace += nBcTria + nBcQuad;
	//	bcOF.nFace = nBcTria + nBcQuad;	
	//	vecBdyOF.push_back(bcOF);


	//	for (iFace = 0; iFace < nBcTria; ++iFace)
	//	{
	//		tri = pTri[vecFaceId[iFace]];
	//		fprintf(fileFaces, "3（%d %d %d)\n", tri.vertices[0],
	//			tri.vertices[1], tri.vertices[2]);

	//		fprintf(fileFaceAd,"%d\n", tri.index+1); // * 注意：从1开始

	//		fprintf(fileOwner, "%d\n", tri.iCell[0]);			
	//	}
	//	for (iFace = nBcTria; iFace < nBcTria + nBcQuad; ++iFace)
	//	{
	//		quad = pQuad[vecFaceId[iFace]];
	//		fprintf(fileFaces, "4(%d %d %d %d)\n", quad.vertices[0],
	//			quad.vertices[1], quad.vertices[2], quad.vertices[3]);

	//		fprintf(fileFaceAd, "%d\n", quad.index+1); // * 注意：从1开始

	//		fprintf(fileOwner, "%d\n", quad.iCell[0]);
	//	}		
	//}


	fprintf(fileFaces,")");	
	fprintf(fileOwner,")");
	fprintf(fileNeigh,")");
	writeTail(fileFaces);
	writeTail(fileOwner);
	writeTail(fileNeigh);
	if (nProcessor > 1)
	{
		fprintf(fileFaceAd,")");
		writeTail(fileFaceAd);
	}

	if (fileFaces)
		fclose(fileFaces);
	if (fileOwner)
		fclose(fileOwner);
	if (fileNeigh)
		fclose(fileNeigh);
	if (nProcessor > 1)
	{		
		if (fileFaceAd) fclose(fileFaceAd);
	}

	return rtn;
}


/* 输出faces文件头信息 */
int OpenFoamFile::writeFacesHead(FILE *fp)
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); ++i)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",facelist.c_str());
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectFace.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/************************************************************************/
/* 输出faceProcAddressing文件头信息                                     */
/************************************************************************/
int OpenFoamFile::writeFacePAdHead(FILE *fp)
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); ++i)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",lablist.c_str());
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectFaceAd.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/************************************************************************/
/* 输出owner文件头信息                                                  */
/************************************************************************/
int OpenFoamFile::writeOwnerHead(FILE *fp, int nPoints, 
								 int nCells, int nFaces, int nInternalFaces)
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); ++i)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",lablist.c_str());
	/*    note        "nPoints: 620 nCells: 2397 nFaces: 5223 nInternalFaces: 4365"*/
	fprintf(fp,"    note        \"nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d\"\n",
		nPoints,nCells,nFaces,nInternalFaces);
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectOwner.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/************************************************************************/
/* 输出neighbour文件头信息                                              */
/************************************************************************/
int OpenFoamFile::writeNeighbHead(FILE *fp, int nPoints,
					int nCells, int nFaces, int nInternalFaces )
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); i ++)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",lablist.c_str());
	/*    note        "nPoints: 620 nCells: 2397 nFaces: 5223 nInternalFaces: 4365"*/
	fprintf(fp,"    note        \"nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d\"\n",
		nPoints,nCells,nFaces,nInternalFaces);
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectNeib.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/* 输出分区网格cellProcAddressing文件 */
int OpenFoamFile::writeProcMesh_OpenFoam_cellProcAddressing()
{
    idx_t i;
	TETRAS *pTet = pHybridMesh->pTetras;
	PYRAMID *pPyra = pHybridMesh->pPyras;
	PRISM *pPrsm = pHybridMesh->pPrisms;
	HEX *pHexa = pHybridMesh->pHexes;
	
    idx_t nTetr = pHybridMesh->NumTetras;
    idx_t nPrsm = pHybridMesh->NumPrsm;
    idx_t nPyra = pHybridMesh->NumPyra;
    idx_t nHexa = pHybridMesh->NumHexes;
    int nCell = nTetr + nPrsm + nPyra + nHexa;


//	/************************************************************************/
//    /* 先统计本分区起始全局单元编号                                         */
//	/************************************************************************/
//	int nCellBeg;
//	int *pCellProc;
//	pCellProc = new int[nProcessor](); // 记录每个分区单元数,初始化为0

//	MPI_Allgather(&nCell, 1, MPI_INT,
//		pCellProc, 1, MPI_INT, MPI_COMM_WORLD);

//	nCellBeg = 0;
//	for (i = 0; i < iProc; ++i)
//	{
//		nCellBeg += pCellProc[i];
//	}
	
	FILE *fp = fopen(cellAddreProc.c_str(), "w");
	if (!fp)
	{
		printf("Cannot open file %s.\n", cellAddreProc.c_str());
		return 1;;/* E_INVALID_ARG */
	}

	writeCellPAdHead(fp);

	fprintf(fp, "%d\n(\n", nCell);
	for (i = 0; i < nTetr; ++i)
	{
        fprintf(fp,"%lld\n", i+nCellBeg);
	}
	for (i = 0; i < nPyra; ++i)
	{
        fprintf(fp,"%lld\n", i+nTetr+nCellBeg);
	}
	for (i = 0; i < nPrsm; ++i)
	{
        fprintf(fp, "%lld\n", i+nTetr+nPyra+nCellBeg);
	}
	for (i = 0; i < nHexa; ++i)
	{
        fprintf(fp, "%lld\n", i+nTetr+nPyra+nPrsm+nCellBeg);
	}

	fprintf(fp,")");
	writeTail(fp);
	if (fp)
	{
		fclose(fp);
	}
	return 0;
}

int OpenFoamFile::writeCellPAdHead(FILE *fp)
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); i ++)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",lablist.c_str());
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectCellAd.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/* 输出分区网格boundary文件 */
int OpenFoamFile::writeMesh_OpenFoam_boundary()
{
	int i;
	string fname;
	if (nProcessor == 1)
		fname = boundary;
	else if(nProcessor > 1)
		fname = boundaryProc;

	FILE *fp = fopen(fname.c_str(), "w");
	if (!fp)
	{
		printf("Cannot open file %s.\n", fname.c_str());
		return 2;	/* E_INVALID_ARG */
	}

	writeBoundaryHead(fp);

	fprintf(fp,"%d\n(\n",vecBdyOF.size());
	for (i = 0; i < vecBdyOF.size(); i++)
	{
		fprintf(fp,"    %s\n",vecBdyOF[i].name.c_str());
		fprintf(fp,"    {\n");
		fprintf(fp,"        type            %s;\n",vecBdyOF[i].type.c_str());
		fprintf(fp,"        nFaces          %d;\n",vecBdyOF[i].nFace);
		fprintf(fp,"        startFace       %d;\n",vecBdyOF[i].startFace);

		if ( i >= pBC->bcInfo.vecBC.size() )
		{/* 交界面边界条件 */
            fprintf(fp,"        matchTolerance  %lf;\n",vecBdyOF[i].matchTolerance);
			fprintf(fp,"        myProcNo        %d;\n",vecBdyOF[i].myProcNo);
			fprintf(fp,"        neighbProcNo    %d;\n",vecBdyOF[i].neighbProcNo);
		}
		fprintf(fp,"    }\n");
	}
	fprintf(fp,")");
	writeTail(fp);
	if (fp)
	{
		fclose(fp);
	}
	return 0;

}

/************************************************************************/
/* 输出文件boundary头信息                                               */
/************************************************************************/
int OpenFoamFile::writeBoundaryHead(FILE *fp)
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); i ++)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",bdymesh.c_str());
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectBdy.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/* 输出分区网格boundaryProcAddressing文件 */
int OpenFoamFile::writeProcMesh_OpenFoam_boundaryProcAddressing()
{
	int i;
	FILE *fp = fopen(bdyAddreProc.c_str(), "w");

	if (!fp)
	{
		printf("Cannot open file %s.\n", bdyAddreProc.c_str());
		return 2;	/* E_INVALID_ARG */
	}

	writeBoundaryPAdHead(fp);
	//3(0 1 -1)
	fprintf(fp,"%d(",vecBdyOF.size());
	for (i = 0; i < vecBdyOF.size(); i++)
	{
		if ( i < pBC->bcInfo.vecBC.size() )
		{/* 非交界面边界条件 */
			fprintf(fp,"%d ",i);
		}
		else
		{/* 交界面边界条件 */
			if ( i == vecBdyOF.size()-1 )
				fprintf(fp,"-1)");
			else
				fprintf(fp,"-1 ");
		}
	}
	writeTail(fp);
	if (fp)
	{
		fclose(fp);
	}

	return 0;
}

/************************************************************************/
/* 输出boundaryProcAddressing文件头信息                                 */
/************************************************************************/
int OpenFoamFile::writeBoundaryPAdHead(FILE *fp)
{
	int i;
	if (!fp)	return 1;	
	for (i = 0; i < vecHead.size(); i ++)
	{
		fprintf(fp,"%s\n",vecHead[i].c_str());
	}
	fprintf(fp,"%s\n",lablist.c_str());
	fprintf(fp,"%s\n",location.c_str());
	fprintf(fp,"%s\n}\n",objectBdyAd.c_str());
    fprintf(fp,"%s\n\n\n","// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");
	return 0;
}

/* 输出文件尾 */
int OpenFoamFile::writeTail(FILE *fp)
{
    if (!fp)	return 1;
    fprintf(fp,"\n\n%s","// ************************************************************************* //");
    return 0;
}

/* 输出并行分块网格 */
int OpenFoamFile::write_OpenFoam_Mesh()
{
	if (nProcessor == 1)
	{
		cout << "******* Begin to output OpenFOAM Mesh in serial model, please wait! *******"<< endl;      
	}
	else
	{
		assert(nProcessor > 1);
#if 0
		if (iProc == 0)
		{
            cout << "******* Proc 0: Begin to output OpenFOAM Mesh in parallel model, please wait! *******"<< endl;
		}
#else
            cout << "******* Proc "<< iProc <<": Begin to output OpenFOAM Mesh in parallel model, please wait! *******"<< endl;
#endif
	}

	creatFile();

	writeMesh_OpenFoam_Points();	
	writeMesh_OpenFoam_AllFaceFiles();
	writeMesh_OpenFoam_boundary(); // 必须在writeMesh_OpenFoam_AllFaceFiles函数之后调用
	if (nProcessor > 1)
	{
		writeProcMesh_OpenFoam_pointProcAddressing();
		writeProcMesh_OpenFoam_boundaryProcAddressing();
		writeProcMesh_OpenFoam_cellProcAddressing();
	}

	if (nProcessor == 1)
	{
        cout << "******* Output OpenFOAM Mesh in serial model succeed!               *******"<< endl;
	}
	else
	{
		assert(nProcessor > 1);
#if 0
		if (iProc == 0)
		{
            cout << "******* Proc 0: Output OpenFOAM Mesh in parallel model succeed!               *******"<< endl;
		}
#else
            cout << "******* Proc "<< iProc << ": Output OpenFOAM Mesh in parallel model succeed!               *******"<< endl;
#endif
	}
	
	return 0;
}
