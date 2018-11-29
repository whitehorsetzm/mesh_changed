#ifndef DATAIO_H
#define DATAIO_H

#include "dataclass.h"
#include "globaldefine.h"
#include "stdio.h"
#include <vector>
#include "refinefunctions.h"

#define MAX_FILE_LINE   1024

#define BKG_MESH_DIM    3

#define MAX_VALUE(x, y) (x) > (y) ? (x) : (y)
using namespace std;

class InterFace
{
public:
    int conn[BKG_MESH_DIM];
    int lftCell, rgtCell;
    int hashNxt;
};

bool isNullOrComment(char* chLine);

int eraseWhiteSpace(char *chLine);

int readValidLine(FILE *fp, char *chLine, int len);

int setupCellNeig(int nNodes, int nElems, TETRAS *pBKGElems);

int readVTKPLSFile(char *fname, HYBRID_MESH&file);

int writeVTKFile(char *filename,HYBRID_MESH&file);

int writeTriangleVTKFile(char *filename,HYBRID_MESH&file);

bool tetrasContainTriangle(int *tri, int *tetras);

bool isSameTriangle(int *tri1,int *tri2);

int findiCell(HYBRID_MESH &file);

int findiCellFast(HYBRID_MESH &file);

//int readVTKPLSFile(char*fname,TETRASFILE&file);

int readCGNS(char*filename, HYBRID_MESH&mesh, vector<string> &bcstring);

#endif // DATAIO_H
