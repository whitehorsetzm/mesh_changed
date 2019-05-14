#ifndef GR_Classes
#define GR_Classes 1

/* Alternating digital tree for efficient searching, either for verts or
   bounding boxes of regions. */
class ADT;

/* Bdry face classes */
class BFace;
#define pBFInvalidBFace (static_cast<BFace*>(NULL))

/* 2D Bdry faces */
class BdryEdgeBase;
class BdryEdge;
class IntBdryEdge;

/* 3D Bdry faces */
class TriBFaceBase;
class TriBFace;
class IntTriBFace;

class QuadBFaceBase;
class QuadBFace;
class IntQuadBFace;

/* Bdry faces with cell-vert connectivity */
class BFaceCV;
class EdgeBFaceCV;
class TriBFaceCV;
class QuadBFaceCV;

/* Ultimate base classes for boundary reps, in both 2D and 3D. */
class BdryPatch;
class BdryRep;

/* Classes for 2D boundary rep */
class Bdry2D;
class BdryPatch2D;
class BdryArc;
class BdryCubicParam2D;
class BdrySeg;

/* Classes for 3D boundary / surface rep */
class Bdry3D;
class BdryPatch3D;
class BdryPolygon;

/* Cell classes */
class CellSkel;
class Cell;
#define pCInvalidCell (static_cast<Cell*>(NULL))

class TriCell;
class QuadCell;

class TetCell;
class PyrCell;
class PrismCell;
class HexCell;

class SimplexCell;

/* Cell classes with cell-vertex connectivity */
class CellCV;
class TriCellCV;
class QuadCellCV;
class TetCellCV;
class PyrCellCV;
class PrismCellCV;
class HexCellCV;

class BadEdge;
class ConstrainEdge;
class ConstrainFace;

/* Face classes */
class Face;
#define pFInvalidFace (static_cast<Face*>(NULL))

class EdgeFace;
class TriFace;
class QuadFace;

class PreserveEdge;

/* Classes for priority queue-based point insertion */
class InsertionQueue;
class InsertionQueueEntry;
class WatsonInfo;

/* List template */
template <class T> class List;

/* Mesh classes */
class Mesh;
class Mesh2D;
class SurfMesh;
class VolMesh;

/* Mesh quality assessment */
class Quality;

/* Vertices */
class Vert;
#define pVInvalidVert (static_cast<Vert*>(NULL))
class VertConnect;

/* Front entry for anisotropic meshing */
class AnisoFront;
class AnisoFrontEntry;
#define pAFEInvalid (static_cast<AnisoFrontEntry*>(NULL))

/* Stuff for global topology optimization. */
class MeshFrag;
class MeshFragContainer;
class FrontFace;
class FeasData;

#endif
