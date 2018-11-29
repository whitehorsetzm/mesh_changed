#ifndef SURFACEIMPRINT_H
#define SURFACEIMPRINT_H


#include "global.h"
#include "math.h"
#define REAL double
typedef REAL *point;

// Labels that signify the result of triangle-triangle intersection test.
enum interresult {DISJOINT, INTERSECT, SHAREVERT, SHAREEDGE, SHAREFACE,
                  TOUCHEDGE, TOUCHFACE, ACROSSVERT, ACROSSEDGE, ACROSSFACE,
                  COLLISIONFACE, ACROSSSEG, ACROSSSUB};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Robust Geometric predicates                                               //
//                                                                           //
// Geometric predicates are simple tests of spatial relations of a set of d- //
// dimensional points, such as the orientation test and the point-in-sphere  //
// test. Each of these tests is performed by evaluating the sign of a deter- //
// minant of a matrix whose entries are the coordinates of these points.  If //
// the computation is performed by using the floating-point numbers, e.g.,   //
// the single or double precision numbers in C/C++, roundoff error may cause //
// an incorrect result. This may either lead to a wrong result or eventually //
// lead to a failure of the program.  Computing the predicates exactly will  //
// avoid the error and make the program robust.                              //
//                                                                           //
// The following routines are the robust geometric predicates for 3D orient- //
// ation test and point-in-sphere test.  They were implemented by Shewchuk.  //
// The source code are generously provided by him in the public domain,      //
// http://www.cs.cmu.edu/~quake/robust.html. predicates.cxx is a C++ version //
// of the original C code.                                                   //
//                                                                           //
// The original predicates of Shewchuk only use "dynamic filters", i.e., it  //
// computes the error at run time step by step. TetGen first adds a "static  //
// filter" in each predicate. It estimates the maximal possible error in all //
// cases.  So it can safely and quickly answer many easy cases.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
extern REAL orient3dfast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
//REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
//

extern REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_test()    Triangle-edge intersection test.                       //
//                                                                           //
// This routine takes a triangle T (with vertices A, B, C) and an edge E (P, //
// Q) in 3D, and tests if they intersect each other.                         //
//                                                                           //
// If the point 'R' is not NULL, it lies strictly above the plane defined by //
// A, B, C. It is used in test when T and E are coplanar.                    //
//                                                                           //
// If T and E intersect each other, they may intersect in different ways. If //
// 'level' > 0, their intersection type will be reported 'types' and 'pos'.  //
//                                                                           //
// The return value indicates one of the following cases:                    //
//   - 0, T and E are disjoint.                                              //
//   - 1, T and E intersect each other.                                      //
//   - 2, T and E are not coplanar. They intersect at a single point.        //
//   - 4, T and E are coplanar. They intersect at a single point or a line   //
//        segment (if types[1] != DISJOINT).                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#define SETVECTOR3(V, a0, a1, a2) (V)[0] = (a0); (V)[1] = (a1); (V)[2] = (a2)

#define SWAP2(a0, a1, tmp) (tmp) = (a0); (a0) = (a1); (a1) = (tmp)

REAL dot(REAL* v1, REAL* v2);

REAL distance(REAL* p1, REAL* p2);


void cross(REAL* v1, REAL* v2, REAL* n);

void facenormal(double * pa, double * pb, double * pc, REAL *n, int pivot,
                            REAL* lav);

int tri_edge_tail(double* A,double * B,double * C,double* P,double * Q,double* R,
                              REAL sP,REAL sQ,int level,int *types,int *pos);

int tri_edge_inter_tail(double*, double*,double*, double*, double*, REAL, REAL);
int tri_tri_inter(double *p1, double *B, double*C, double*O, double *P, double*Q);

int interecursive(DiscretSolid&solid, DiscretPoint *pointsArray, DiscretFacet*facetsArray, int arraySize, int axis, double bxmin, double bxmax, double bymin, double bymax,
              double bzmin, double bzmax, int &internum, multimap<int,int>&entiy_inter_entity, unordered_map<string, bool> &visit, int type);

//type = 1 points intersect lines;type =2 points intersection triangle ; type =3 lines intersectioin lines
int detectFacetIntersection(DiscretSolid&discreteSolid,multimap<int,int>&entity_inter_entity,int type);

bool intersectionProcess(int facet1,int facet2,DiscretSolid&solid);


int process_point_inter_edge_remesh(DiscretSolid&discreteSolid,multimap<int,int>&entity_inter_entity);

int process_point_inter_edge(DiscretSolid&discreteSolid,multimap<int,int>&entity_inter_entity);

int process_edge_inter_edge(DiscretSolid&discreteSolid,multimap<int,int>&entity_inter_entity);

int process_point_inter_triangle(DiscretSolid&discreteSolid,multimap<int,int>&entity_inter_entity);

int process_point_inter_triangle_remsh(DiscretSolid&discreteSolid,multimap<int,int>&entity_inter_entity);

bool checkInterscetion(Vector c, Vector a,Vector b);//point intersection edge

double getDistance(Vector c, Vector a, Vector b);//a,b is the endpoints on line, c is the point being tested

Vector getProjectPointOnEdge(Vector c, Vector a, Vector b);

Vector getProjectPointOnTrianglePlane(Vector p, Vector A, Vector B, Vector C);

bool checkIntersection_P_T(Vector P, Vector A, Vector B, Vector C);//point intersection triangle

bool determinePointInTriangle(Vector P, Vector A, Vector B, Vector C);

bool checkIntersection_E_E(Vector edgeP1,Vector edgeP2,Vector edgeP3,Vector edgeP4);//edge intersection edge

int resolve_point_inter_edges(DiscretSolid&solid);

int resolve_point_inter_triangles(DiscretSolid&solid);

int resolve_edge_inter_edge(DiscretSolid&solid);

int discreteSolid_to_validDiscreteSolid(DiscretSolid&solid);

int clearEdgeInformation(DiscretSolid&solid);

bool IsEqual(double d1, double d2);

double DistanceSegmentToSegment(double x1, double y1, double z1,
                                            double x2, double y2, double z2,
                                            double x3, double y3, double z3,
                                            double x4, double y4, double z4);

Vector getIntersectPointOnEdge_Edge(double x1, double y1, double z1,
                                  double x2, double y2, double z2,
                                  double x3, double y3, double z3,
                                  double x4, double y4, double z4);

int splitEdgeAndTriangle(DiscretSolid&discreteSolid,int pointID,int EdgeID);//the point is on the edge , the function of the function is to split the edge
                                                                            //edge and divide the triangle
int mergePoints(DiscretSolid&solid);//merge points by the distance

int clearTriangle(DiscretSolid&solid);

int formConnecteDomain(const DiscretSolid &solid, int faceID, const set<int> &facets, vector<set<int> > &connecteDomain);

#endif // SURFACEIMPRINT_H
