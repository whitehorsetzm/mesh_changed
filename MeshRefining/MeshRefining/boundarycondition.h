#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <cstdio>
#include <string>
#include <vector>

using namespace std;


typedef struct BC
{
    string name;
    string type;
    vector<int> faceIDs;
}BC;

typedef struct BCInfo
{
    int numOfFaces;
    vector<string> types;
    vector<BC> vecBC;
}BCInfo;

class BoundaryCondition
{
public:
    BoundaryCondition(const char *filename);
    BoundaryCondition(const vector<string> vecBcName);

    int readFileBc(const char *filename);

public:

    BCInfo bcInfo;

};

#endif // BOUNDARYCONDITION_H
