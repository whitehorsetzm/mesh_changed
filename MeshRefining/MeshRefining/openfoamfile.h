#ifndef OPENFOAMFILE_H
#define OPENFOAMFILE_H

#include "dataclass.h"
#include "boundarycondition.h"

#include <string>
#include <vector>
using namespace std;

/* 子区域所包含的边界条件,在生成区域网格的时候确定 */
typedef struct BoundaryConditionOpenFoam
{
	string	   name;
	string     type;
	int	       nFace;
	int        startFace;
    double     matchTolerance;
	int 	   myProcNo;
	int        neighbProcNo;
}BoundaryOF;

class OpenFoamFile
{
public:
    OpenFoamFile(int iRank, const HYBRID_MESH *mesh, idx_t nFace, idx_t nCell, const BoundaryCondition *bc, string path = "./", string name = "OpenFoamFile", int nProc = 0, idx_t nAllInFace = 0);

    int creatFile();

    /* 路径、分块数、文件夹名称 */
    string homepath;
    int nProcessor;
    int iProc;
    string filename;
    const idx_t nFaceBeg;
    const idx_t nCellBeg;
    const idx_t nAllInternalFace; // useless when define macro _INTERNAL_FACE_FIRST_

    /* 总体网格文件名 */
    string points;
    string faces;
    string neighbour;
    string owner;
    string boundary;
    string field;
    string system;

    /* 分块网格文件名 */
//    string *pointsProc;
//    string *pointAddreProc;
//    string *facesProc;
//    string *faceAddreProc;
//    string *neighbourProc;
//    string *ownerProc;
//    string *cellAddreProc;
//    string *boundaryProc;
//    string *bdyAddreProc;
//    string *fieldProc;
    string pointsProc;
    string pointAddreProc;
    string facesProc;
    string faceAddreProc;
    string neighbourProc;
    string ownerProc;
    string cellAddreProc;
    string boundaryProc;
    string bdyAddreProc;
    string fieldProc;

private:
    /* 输出分区网格points文件 */
    int writeMesh_OpenFoam_Points();

    /* 输出文件points头信息 */
    int writePointsHead(FILE *fp);

    /* 输出分区网格pointProcAddressing文件 */
    int writeProcMesh_OpenFoam_pointProcAddressing();

    /* 输出文件pointProcAddressing头信息 */
    int writePointPAdHead(FILE *fp);

	/* 输出分区网格面片相关文件,包括faces, owner, neighbor, faceProcAddressing */
	int writeMesh_OpenFoam_AllFaceFiles();

	/* 输出faces文件头信息 */
	int writeFacesHead(FILE *fp);

	/* 输出faceProcAddressing文件头信息 */
	int writeFacePAdHead(FILE *fp);

	/* 输出owner文件头信息 */
	int writeOwnerHead(FILE *fp, int nPoints, 
		int nCells, int nFaces, int nInternalFaces);

	/* 输出neighbour文件头信息 */
	int writeNeighbHead(FILE *fp, int nPoints, 
		int nCells, int nFaces, int nInternalFaces);

    ///* 输出分区网格faces文件 */
    //int writeProcMesh_OpenFoam_Faces();

    ///* 输出分区网格faceProcAddressing文件 */
    //int writeProcMesh_OpenFoam_faceProcAddressing();

    ///* 输出分区网格owner文件 */
    //int writeProcMesh_OpenFoam_owner();

    ///* 输出分区网格neighbor文件 */
    //int writeProcMesh_OpenFoam_neighbor();

    /* 输出分区网格cellProcAddressing文件 */
    int writeProcMesh_OpenFoam_cellProcAddressing();

	/* 输出cellProcAddressing文件头信息 */
	int writeCellPAdHead(FILE *fp);

    /* 输出分区网格boundary文件 */
    int writeMesh_OpenFoam_boundary();

	/* 输出文件boundary头信息 */
	int writeBoundaryHead(FILE *fp);

	/* 输出分区网格boundaryProcAddressing文件 */
	int writeProcMesh_OpenFoam_boundaryProcAddressing();

	/* 输出boundaryProcAddressing文件头信息 */
	int writeBoundaryPAdHead(FILE *fp);

    /* 输出文件尾 */
    int writeTail(FILE *fp);

public:
	/* 输出网格 */
	int write_OpenFoam_Mesh();

	/* 输出整体网格 */



private:
    const HYBRID_MESH* pHybridMesh;


    vector<string> vecHead;
    string lablist;
    string vecfield;
    string bdymesh;
    string facelist;
    string location;
    string objectPoint;
    string objectPtPAd;
    string objectBdy;
    string objectBdyAd;
    string objectCellAd;
    string objectFaceAd;
    string objectFace;
    string objectNeib;
    string objectOwner;

	const BoundaryCondition *pBC;
	vector<BoundaryOF> vecBdyOF;

};

#endif // OPENFOAMFILE_H
