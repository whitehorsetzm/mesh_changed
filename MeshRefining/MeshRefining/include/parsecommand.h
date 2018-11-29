#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "config.h"
#include "global.h"

using namespace std;

#define ERRORCONFIG -1
#define DEFAULTCONFIG 1
#define INPUTCONFIG 0
#define NODEFAULTCONFIG -2

char config[FILENAMESIZE];

int GetConfigName(int argc,char ** argv)
{
	int i=argc;
	if(i>2)
	{
		printf("Error : Please type in the config filename.\n");
		exit(ERRORCONFIG);
	}
	else if(i==1)
	{
		printf("Note : Use the default config.\n");
		strcpy(config,DEFAULTCFNAME);
		return DEFAULTCONFIG;
	}
	else
	{
		strcpy(config,argv[1]);
		return INPUTCONFIG;
	}
}

void SetDefaultConfig(ConfigArgc &cf)
{
	strcpy(cf.filenam, "automesh");
    cf.NumFile=-1;
	cf.gTolerance = 1.0e-6;
	cf.sTolerance = 1.0e-6;
	cf.del_small_edges = false;
	cf.geotype = 1;
	cf.unum = 4;
	cf.vnum = 4;
	cf.outflow = 0;
	cf.L=24;
	cf.xmin=0;
	cf.xmax=0;
	cf.ymax=0;
	cf.ymin=0;
	cf.zmax=0;
	cf.zmin=0;
	cf.center_coord[0]=0;
	cf.center_coord[1]=0;
	cf.center_coord[2]=0;

	cf.proximity_num=1;
    cf.curvature_angle=10;
    cf.expand_ratio=1.2;
	cf.step=500;
	strcpy(cf.method,"ma86");
	cf.hull=0;
	cf.outputM=0;
	cf.inner_spacing=BOXSIZE;
	cf.conR=0;
	cf.conV=0;
	cf.hullcoord.clear();
	//cf.method="ma57";
}

int parsecommand(int argc, char ** argv, ConfigArgc &cf)
{
	char word[100];      /* variable name from config file */
    int value=0;         /* variable value from config file */
    double fvalue=0.0; /* variable for floats in config file */
    char str[100];
    char fstr[100];
	char temp[100];
    
    /* temp stuff for line reading */
    int nbytes = 1000;
    char line[1000];
    int numassigned;
	double coord[3];
	char rotationaxis;
    int NumFile=-1;

	GetConfigName(argc,argv);
	FILE *fp=fopen(config,"a");
	if(!fp)
	{
		printf("Config File does not exist.\n");
		exit(NODEFAULTCONFIG);
	}
	
	SetDefaultConfig(cf);

	fclose(fp);
	fp=fopen(config,"r");
	if(!fp)
	{
		printf("Config File does not exist.\n");
		exit(NODEFAULTCONFIG);
	}

	
	while (fgets(line, nbytes, fp) != NULL)
	{
		/* attempt to fetch a variable name and value from the config file */
		numassigned = sscanf(line, "%s %s", word, temp);
		fvalue = atof(fstr);
		/* check if this is a comment */
		if (word[0] == '#' || word[0] == '\n' || numassigned < 2) continue;

		if (!strcmp(word, "filename")) {
            sscanf(line, "%s %s %d", word, str,&NumFile);
			strcpy(cf.filenam, str);
            cf.NumFile=NumFile;
		}
		if (!strcmp(word, "command")) {
			sscanf(line, "%s %s", word, str);
			strcpy(cf.command, str);
		}
		if (!strcmp(word, "sTolerance")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.sTolerance = fvalue;
		}
		if (!strcmp(word, "gTolerance")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.gTolerance = fvalue;
		}
		if (!strcmp(word, "del_small_edges")) {
			sscanf(line, "%s %d", word, &value);
			cf.del_small_edges = value;
		}
		if (!strcmp(word, "hull")) {
			sscanf(line, "%s %d", word, &value);
			cf.hull = value;
		}

		if (!strcmp(word, "outputM")) {
			sscanf(line, "%s %d", word, &value);
			//cf.hull = value;
			cf.outputM=value;
		}


		if (!strcmp(word, "hull_coord")) {
			sscanf(line, "%s %lf %lf %lf", word, &coord[0],&coord[1],&coord[2]);
			COORDS temp;
			temp.cood[0]=coord[0];
			temp.cood[1]=coord[1];
			temp.cood[2]=coord[2];
			/*
			cf.hullcoord[0]=coord[0];
			cf.hullcoord[1]=coord[1];
			cf.hullcoord[2]=coord[2];
			*/
			cf.hullcoord.push_back(temp);
		}
		
		if (!strcmp(word, "global_spacing")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.spacing = fvalue;
			cf.global_spacing=fvalue;
		}

		if (!strcmp(word, "inner_spacing")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.inner_spacing = fvalue;
			//cf.global_spacing=fvalue;
		}

		if (!strcmp(word, "conR")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.conR = fvalue;
			//cf.global_spacing=fvalue;
		}
		if (!strcmp(word, "conV")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.conV = fvalue;
			//cf.global_spacing=fvalue;
		}


		if (!strcmp(word, "min_spacing")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.min_spacing=fvalue;
		}
		if (!strcmp(word, "curvature_angle")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.curvature_angle=fvalue;
		}
		if (!strcmp(word, "expand_ratio")) {
			sscanf(line, "%s %lf", word, &fvalue);
			cf.expand_ratio=fvalue;
		}
		if (!strcmp(word, "geometry_type")) {
			sscanf(line, "%s %d", word, &value);
			cf.geotype = value;
		}
		if (!strcmp(word, "unum")) {
			sscanf(line, "%s %d", word, &value);
			cf.unum = value;
		}
		if (!strcmp(word, "step")) {
			sscanf(line, "%s %d", word, &value);
			cf.step = value;
		}
		if (!strcmp(word, "proximity_num")) {
			sscanf(line, "%s %d", word, &value);
			cf.proximity_num = value;
		}
		if (!strcmp(word, "vnum")) {
			sscanf(line, "%s %d", word, &value);
			cf.vnum = value;
		}
		if (!strcmp(word, "bkfaceu")) {
			int iface;
			sscanf(line, "%s %d %lf", word, &iface, &fvalue);
			cf.mapfcbku.insert(std::make_pair(iface-1, fvalue));
		}
		if (!strcmp(word, "bkfacev")) {
			int iface;
			sscanf(line, "%s %d %lf", word, &iface, &fvalue);
			cf.mapfcbkv.insert(std::make_pair(iface-1, fvalue));
		}
		if (!strcmp(word, "face")) {
			int iface;
			sscanf(line, "%s %d %lf", word, &iface, &fvalue);
			cf.mapFaceSize.insert(std::make_pair(iface-1, fvalue));
		}
	}
	cf.gTolerance = cf.sTolerance;
	fclose(fp);

	return 0;
}




