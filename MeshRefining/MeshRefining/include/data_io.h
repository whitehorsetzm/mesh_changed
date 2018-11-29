#ifndef DATA_IO_H
#define DATA_IO_H
#include<vector>
#include<set>
#include<unordered_map>
using namespace std;
class PointInput
{
public:
    int index;
    double coord[3];
public:
    PointInput();

};
class VertexInput
{
public:
    int PointID;
    bool valid;
    set<int>attributes;
public:
    VertexInput();



};
class LineInput
{
public:
    int pointID1;
    int pointID2;
    bool valid;
    set<int>attributes;
public:
    LineInput();

};
class TriangleInput
{
public:
    int points[3];
    bool valid;
    set<int>attributes;
public:
    TriangleInput();

};

class VTKInput
{
public:
    int NumPoints;
    int NumLines;
    int NumVertex;
    int NumValidVertex;//used to avoid overlap
    int NumValidLines;//used to avoid overlap
    int NumTraingles;
    int NumValidTriangles;//used to avoid overlap
    PointInput*points;
    VertexInput*vertices;
    LineInput*lines;
    TriangleInput*triangles;

    unordered_map<int,bool>validVertex;//used to avoid overlap
    unordered_map<string,bool>validLines;//used to avoid overlap
    unordered_map<string,bool>validTraingles;//used to avoid overlap
public:
    VTKInput();

    int clear();


    ~VTKInput();

};

class GM3Point
{
public:
    GM3Point();
public:
    double coord[3];
    int index;
    bool vertex;
    set<int>parents;//curve
};
class GM3Curve
{
public:
    GM3Curve();
    ~GM3Curve();
public:
    int index;
    int NumPoints;
    int *pointIDs;
    bool valid;
    set<int>parents;//face
    set<int>mesh;//unvalid, because the vtkFile does not save all the edges

};
class GM3Face
{
public:
    GM3Face();
    ~GM3Face();
public:
    int index;
    int NumU;
    int NumV;
    int **pointIDs;
    int LoopID;
    bool valid;
    set<int>parents;//body
    set<int>mesh;
};
class GM3Loop
{
public:
    GM3Loop();
    ~GM3Loop();
public:
    int index;
    vector<int>curveIDs;
    bool valid;
};
class GM3Body
{
public:
    GM3Body();
    ~GM3Body();
public:
    int index;
    vector<int>faceIDs;
};
class GM3Data
{
public:
    GM3Data();
    ~GM3Data();
    int clear();
public:
    int NumPoints;
    int NumCurves;
    int NumFaces;
    int NumLoops;
    int NumBodies;
    GM3Point *points;
    GM3Curve *curves;
    GM3Face *faces;
    GM3Loop *loops;
    GM3Body *bodies;

};


#endif // DATA_IO_H
