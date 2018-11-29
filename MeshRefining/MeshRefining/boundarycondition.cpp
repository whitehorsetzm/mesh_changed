#include "boundarycondition.h"


BoundaryCondition::BoundaryCondition(const char *filename)
{
    readFileBc(filename);
}

BoundaryCondition::BoundaryCondition(const vector<string> vecBcName)
{
    int numBc = vecBcName.size();
    bcInfo.numOfFaces = numBc; // here, one Boundary condition represent one face;
	bcInfo.vecBC.resize(numBc);
    for( int i = 0; i < numBc; ++i )
    {
        bcInfo.vecBC[i].name = vecBcName[i];
        bcInfo.vecBC[i].type = "";
        bcInfo.vecBC[i].faceIDs.push_back(i+1); // start form 1, as the same as read from a bc file
    }
}



int BoundaryCondition::readFileBc(const char *filename)
{
    FILE *fin = NULL;
    int faceNum, faceID, numOfBCs;

    if ((fin = fopen(filename, "r")) == 0)
    {
        fprintf(stderr, "Can't open file BC (%s).\n", filename);
        fclose(fin);
        return 1;
    }

    fscanf(fin,"%d %d\n", &bcInfo.numOfFaces, &numOfBCs);

    if (bcInfo.numOfFaces <= 0 || numOfBCs <= 0)
    {
        fprintf(stderr, "Unexpected numOfFaces (%d) or # of BCs (%d) defined in layer data file (%s).\n",
            bcInfo.numOfFaces, numOfBCs, filename);
        return 2;   /* E_WRONG_FILE_FORMAT */
    }

    bcInfo.vecBC.resize(numOfBCs);

    /* bc names */
    string strFullName;
    char cFullName[1024];
    string::size_type pos;
    for(int i = 0; i<numOfBCs;i++)
    {
        fscanf(fin,"%s\n",&cFullName);
        strFullName = string(cFullName);
        pos = strFullName.find_first_of("_",0);
        if(pos != string::npos)
        {
            string strType(strFullName,0,pos);
            string strName(strFullName,pos+1,strFullName.length()-pos);
            bcInfo.vecBC[i].type = strType;
            bcInfo.vecBC[i].name = strName;
        }
        else
        {
            bcInfo.vecBC[i].type = "";
            bcInfo.vecBC[i].name = strFullName;
        }

    }

    for(int i = 0;i<numOfBCs; i++)
    {
        fscanf(fin, "%d\n",&faceNum);
        if(faceNum > 0 && faceNum <= bcInfo.numOfFaces)
        {
            for(int j = 0;j<faceNum;j++)
            {
                fscanf(fin,"%d\n",&faceID);
                if(faceID > 0 && faceID <= bcInfo.numOfFaces)
                {
                   // sprintf(key,"Face%d",faceID);
                    bcInfo.vecBC[i].faceIDs.push_back(faceID);
                }
            }
        }
    }
    if(fin)
        fclose(fin);

    return 0;
}
