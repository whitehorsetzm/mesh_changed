/*
********struct of the config to define the property**********
**
**
*/
#ifndef CONFIG_H
#define CONFIG_H
#pragma once
#include <map>
#include<vector>


#define FILENAMESIZE 120
#define DEFAULTCFNAME "TRAN_CONFIG"

//char config[FILENAMESIZE];

class COORDS
{
public:
	double cood[3];

};

typedef struct {
	char filenam[FILENAMESIZE];		// the project name
    int NumFile;                    //the number of files
	char command[FILENAMESIZE];
	double sTolerance;				//set sTolerance, default value 1.0e-8
	double gTolerance;				//set gTolerance, default value 1.0e-8
	bool   del_small_edges;			//set del_small_edges, default value false
	double spacing;
	int geotype;
	int unum;
	int vnum;
	int outflow;                   //option for generating a outflow field
	double L;
	double xmin,ymin,zmin;
	double xmax,ymax,zmax;
	double center_coord[3];
	double a;                        //half long axis of elipsoid
	double b;                        //half short axis of elipoid
	char rotate;                     // the rotation axis;

	int step;
	char method[256];

	double global_spacing;
    double min_spacing;
	double inner_spacing;
	double conV;
	double conR;

	int hull;
	//double hullcoord[3];

	std::vector<COORDS> hullcoord;

    int  proximity_num;
    double curvature_angle;
    double expand_ratio;
	int outputM;

	std::map<int, double> mapfcbku;
	std::map<int, double> mapfcbkv;

	std::map<int,double> mapFaceSize;
}ConfigArgc;



extern ConfigArgc cf;

#endif
