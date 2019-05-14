#ifndef __tetra_mesh_h__
#define __tetra_mesh_h__

#include "iso3d_define.h"
#include <vector>
// 
// typedef struct TetraCell {
// 	int form[4];
// 	int neig[4];
//  } TetraCell;
#define MAX_NODE_INTERIOR_FACE_SIZE 20384//5096
#define MAX_NODE_CELL_SIZE 20384//5096

class TetraMesh {
	enum {REAL_MESH, TOPU_MESH};
	typedef struct InterFace {
		int form[3];
		int lftCell, rgtCell;
		int hashNxt;
	} InterFace;

public:
	TetraMesh();
	~TetraMesh();
	
	virtual int initialize(REAL nodes[], int ndSize, int clNodes[], int clSize, bool autoNeig = true);
	virtual int initialize(REAL nodes[], int ndSize, int clNodes[], int clNeigs[], int clSize);
	
	int nodeSize();
	int nodeCoords(int iNode, REAL coords[]);
//	int addNode(double coords[]);

	int cellSize();
	int cellNodes(int iCell, int clNodes[]);

	int cellNeigs(int iCell, int clNeigs[]);
	int cellNeig(int iCell, int iCode);

	int nodeCells(int node, int cells[], int *size);
	int edgeCells(int node1, int node2, int cells[], int *size);
	int faceCells(int node1, int node2, int node3, int *lftCell, int *rgtCell);
	static bool isLeftCell(int code1, int code2, int code3);

protected:
	int initialNodeCellHash(int ndSize);
	int setupCellNeig(int ndSize, int clSize);
	int nodeInterFace(int minFacNdIdx, int ndIFaces[], int *ndIFaceSize, 
		std::vector<InterFace> &vecInterFaces, std::vector<int>& vecRefIntFHash);
	int isCellNode(int iCell, int iNode);

protected:
	std::vector<REAL>		m_vecNodeCords;
	std::vector<INTEGER>	m_vecCellNodes;
	std::vector<INTEGER>	m_vecCellNeigs;
	std::vector<INTEGER>    m_vecNodClHash;
	int m_nMaxNodeIdx;
//	int m_meshType;
};

#endif /* __tetra_mesh_h__ */