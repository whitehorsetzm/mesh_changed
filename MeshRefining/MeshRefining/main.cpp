#define main_1
#ifdef main_
#include "smooth.h"
#endif
#ifdef main_1
#include <iostream>
#include "dataio.h"
#include "string.h"
#include "refinefunctions.h"
#include "mpi.h"
#include "assert.h"
#include "openfoamfile.h"
#include "boundary.h"
#include "dataclass.h"
#include "fstream"
//#include "smooth.h"
#define VTK_OUTPUT;
#define DEBUG
using namespace std;
int writepl3(char *filename, HYBRID_MESH &file)
{
    int npt,nface,ncell,i;
    FILE *fout = fopen(filename,"w");
    npt=file.NumNodes;
    nface=file.NumTris;
    ncell=file.NumTetras;
     fprintf(fout,"%d %d %d\n",ncell,npt,nface);
    for (i = 0; i < npt; i++)
    {
        int index = i + 1;
        double x = file.nodes[i].coord.x;
        double y = file.nodes[i].coord.y;
        double z = file.nodes[i].coord.z;
        fprintf(fout,"%d %lf %lf %lf\n",index,x, y, z);

    }

    for (i = 0; i < ncell; i++)
    {

        int index = i + 1;
        int form0 = file.pTetras[i].vertices[0]+1;///////////+1
        int form1 = file.pTetras[i].vertices[1]+1;
        int form2 = file.pTetras[i].vertices[2]+1;
        int form3 = file.pTetras[i].vertices[3]+1;

        fprintf(fout,"%d %d %d %d %d\n",index,form0,form1,form2,form3);

    }

    for (i = 0; i < nface; i++)
    {

        int index = i + 1;
        int face0 = file.pTris[i].vertices[0]+1;
        int face1 = file.pTris[i].vertices[1]+1;
        int face2 = file.pTris[i].vertices[2]+1;
        int parent= file.pTris[i].iCell+1;
   //     cout<<parent<<endl;
        int iPatch= file.pTris[i].iSurf;
//        cout<<iPatch+1<<" "<<iPatch<<endl;
        fprintf(fout,"%d %d %d %d %d %d\n",index,face0,face1,face2,parent,iPatch);

    }
   fclose(fout);


}
struct Arguments
{
    Arguments(int nOutputParts)
        : caseName(nullptr), expectedVol(0), refineTimes(1),
          // nWorkingProcs(num_working_procs)
          nOutputParts(nOutputParts)
    { }

    const char *caseName;
    const char *gmName;
    int expectedVol;
    int refineTimes;
    int nOutputParts;
    // int nWorkingProcs;
};

static void printHint()
{
    printf("MeshRefining - A automatic mesh refining tool.\n");
    printf("Usage: mpirun|mpiexec -np <nprocs> MeshRefining <case_name>\n"
           "\t\t[-e|--expected-volume <expected_volume>]\n"
           "\t\t[-o|--output-parts <output_parts> | \n"
           // "\t\t[-p|--num-working-procs <num_working_procs> | \n"
           "\t\t -r|--refine-times <refine_times>]\n");
    printf("\t<case_name>: file name of the case (no suffix included).\n");
    printf("\t<expected_volume>: expected number of cells\n"
           "\t\t(default: refining for only once).\n");
    printf("\t<output_parts>: number of output parts\n"
           "\t\t(default: equals to the number of MPI processes).\n");
    // printf("\t<num_working_procs>: number of working processes\n"
    // 	   "\t\t(default: equals to the number of MPI processes).\n");
}

static bool checkArgs(const int argc, const char *argv[], Arguments& arguments)
{
    if (argc < 2)
    {
        printf("Error: too few arguments!\n");
        return false;
    }
    arguments.caseName = argv[1];
    for (int i = 2; i < argc; ++i)
    {
        bool error = false;
        if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--expected-volume"))
        {
            int count = sscanf(argv[++i], "%d", &arguments.expectedVol);
            if (count != 1)
                error = true;
        }/* else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--num-working-procs"))
        {
            int count = sscanf(argv[++i], "%d", &arguments.nWorkingProcs);
            if (count != 1)
                error = true;
        }*/ else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--refine-times"))
        {
            int count = sscanf(argv[++i], "%d", &arguments.refineTimes);
            if (count != 1)
                error = true;
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output-parts"))
        {
            int count = sscanf(argv[++i], "%d", &arguments.nOutputParts);
            if (count != 1)
                error = true;
        } else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--read-geometry"))
        {
            arguments.gmName=argv[++i];
        }

        if (error)
        {
            printf("Error: invalid argument \"%s\"!\n", argv[i]);
            return false;
        }
    }
    return true;
}

/*
static MPI_Comm buildWorkComm(MPI_Comm commWorld, int nWorkingProcs)
{
    int rank, commWorldSize;
    MPI_Comm_rank(commWorld, &rank);
    MPI_Comm_size(commWorld, &commWorldSize);
    int color = rank < nWorkingProcs ? 0 : MPI_UNDEFINED;
    int key = rank;
    MPI_Comm workComm;
    MPI_Comm_split(commWorld, color, key, &workComm);
    return workComm;
}

static MPI_Comm buildRescatterComm(MPI_Comm commWorld, int nWorkingProcs)
{
    int rank, commWorldSize;
    MPI_Comm_rank(commWorld, &rank);
    MPI_Comm_size(commWorld, &commWorldSize);
    int color = rank % nWorkingProcs;
    int key = rank;
    MPI_Comm rescatterComm;
    MPI_Comm_split(commWorld, color, key, &rescatterComm);
    return rescatterComm;
}
*/

int main1(int argc, char *argv[])
{
    int rank=-1, size=-1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPITypeAssistant *mpiTypeAssistant = new MPITypeAssistant(); // Used for MPI type constructing and destructing.
    HYBRID_MESH newTetrasfile;
    HYBRID_MESH construcTetras;
    Refletion ref;
 //   Refletion construcref;
    Arguments arguments(size);
    if (!checkArgs(argc, const_cast<const char**>(argv), arguments))
    {
        printHint();
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    char filename[256];
    char filenameBC[256];

    sprintf(filename, "%s", arguments.caseName);



    int refineTimes = arguments.refineTimes;
    int expectVolume = arguments.expectedVol;
    int nOutputParts = arguments.nOutputParts;


    // cout<<expectVolume<<endl;
    // cout<<filename<<endl;
    strcpy(filenameBC, filename);
    strcat(filenameBC, ".bc");
    vector<string>bcstring;
    HYBRID_MESH *meshParts = nullptr;
    HYBRID_MESH tetrasfile;
    /*
    MPI_Comm workComm = buildWorkComm(MPI_COMM_WORLD, arguments.nWorkingProcs);
    MPI_Comm rescatterComm = buildRescatterComm(MPI_Comm, arguments.nWorkingProcs);
    int rescatterCommSize = -1;
    int rescatterCommRank = -1;
    MPI_Comm_size(workComm, &rescatterCommSize);
    MPI_Comm_rank(workComm, &rescatterCommRank);
    if (workComm != MPI_COMM_NULL)
    {
        int workCommSize = -1;
        int workCommRank = -1;
        MPI_Comm_size(workComm, &workCommSize);
        MPI_Comm_rank(workComm, &workCommRank);
    */

//  size=2;
//    char CGNSfile[256];
//   char GM3file[256];
//    strcpy(CGNSfile, filename);
//    strcat(CGNSfile, ".cgns");
//    cout<<CGNSfile<<endl;
//    readCGNS(CGNSfile,tetrasfile,bcstring);
//    if(arguments.gmName!=nullptr){
//        char gm3name[256];
//    sprintf(gm3name, "%s", arguments.gmName);
//        strcpy(GM3file, gm3name);
//        strcat(GM3file, ".gm3");
//        cout<<GM3file<<endl;
//        double coord[tetrasfile.NumNodes*3];
//        int    vertices[tetrasfile.NumTris*3];
//        for(int i=0;i<tetrasfile.NumNodes;++i){
//            coord[i*3+0]=tetrasfile.nodes[i].coord.x;
//            coord[i*3+1]=tetrasfile.nodes[i].coord.y;
//            coord[i*3+2]=tetrasfile.nodes[i].coord.z;
//        }
//        for(int i=0;i<tetrasfile.NumTris;++i){
//            vertices[i*3+0]=tetrasfile.pTris[i].vertices[0];
//            vertices[i*3+1]=tetrasfile.pTris[i].vertices[1];
//            vertices[i*3+2]=tetrasfile.pTris[i].vertices[2];
//        }
//        cout<<"test1"<<endl;
//        ref.initial(GM3file,tetrasfile.NumNodes,coord,tetrasfile.NumTris,vertices);
//        for(int i=0;i<tetrasfile.NumTris;++i){
//            tetrasfile.pTris[i].iSurface=ref.subject_table[i];
//        }
    if(arguments.gmName!=nullptr){
        cout<<"read gm3!!!!"<<endl;
        char GM3file[256];
        char gm3name[256];
    sprintf(gm3name, "%s", arguments.gmName);
        strcpy(GM3file, gm3name);
        strcat(GM3file, ".gm3");
     ref.read_gm(GM3file);
    }

    if (rank == 0)
    {
        char CGNSfile[256];
        char GM3file[256];
        strcpy(CGNSfile, filename);
        strcat(CGNSfile, ".cgns");
         cout<<CGNSfile<<endl;
           cout<<"rank =========="<<rank<<endl;
         readCGNS_temp(CGNSfile,tetrasfile,bcstring);
         if(arguments.gmName!=nullptr){                                 //??????????prolem
             char gm3name[256];
         sprintf(gm3name, "%s", arguments.gmName);
             strcpy(GM3file, gm3name);
             strcat(GM3file, ".gm3");
             cout<<GM3file<<endl;
             double coord[tetrasfile.NumNodes*3];
             int    vertices[tetrasfile.NumTris*3];
             for(int i=0;i<tetrasfile.NumNodes;++i){
                 coord[i*3+0]=tetrasfile.nodes[i].coord.x;
                 coord[i*3+1]=tetrasfile.nodes[i].coord.y;
                 coord[i*3+2]=tetrasfile.nodes[i].coord.z;
             }
             for(int i=0;i<tetrasfile.NumTris;++i){
                 vertices[i*3+0]=tetrasfile.pTris[i].vertices[0];
                 vertices[i*3+1]=tetrasfile.pTris[i].vertices[1];
                 vertices[i*3+2]=tetrasfile.pTris[i].vertices[2];
             }
//             GBSolid gbsolid;
//             read_gm3(GM3file,&gbsolid);
             ref.initial(tetrasfile.NumNodes,coord,tetrasfile.NumTris,vertices);
             for(int i=0;i<tetrasfile.NumTris;++i){
                 tetrasfile.pTris[i].iSurface=ref.subject_table[i];
             }
         }
         setupCellNeig(tetrasfile.NumNodes,tetrasfile.NumTetras,tetrasfile.pTetras);
//         ofstream neighbor;
//         neighbor.open("neighbor.txt");
//          for(int i=0;i<tetrasfile.NumTetras;++i){
//              neighbor<<"neighbor size : "<<i<<endl;
//               for(auto a:tetrasfile.pTetras[i].neighbors)      //运行有不确定结果   ok
//                neighbor<<a<<" ";
//               neighbor<<endl;
//          }
         findiCellFast(tetrasfile);
        writepl3("test.pl3", tetrasfile);
       // exit(1);
    //    readVTKPLSFile(filename,tetrasfile);
        tetrasfile.NumUniqueSurfFacets=tetrasfile.NumTris;
        int nparts=size;
        if (expectVolume != 0)
            refineTimes = expectVolume/tetrasfile.NumTetras/8;

        meshParts = new HYBRID_MESH [nparts];        
        partition(tetrasfile,meshParts,nparts);


#ifdef VTK_OUTPUT
        writeVTKFile("cube_partition.vtk", tetrasfile);
        writeTriangleVTKFile("cube_surface.vtk", tetrasfile);
#endif

   //     exit(1);
        managerProc(rank,nparts,meshParts,construcTetras,refineTimes,bcstring,MPI_COMM_WORLD);
        delete []meshParts;
        meshParts=nullptr;

    }
    else
        workerProc(rank,construcTetras,refineTimes,bcstring,MPI_COMM_WORLD);

    cout<<"refine times "<<refineTimes<<endl;


    for (int i = 0; i < bcstring.size(); ++i)
        cout << bcstring[i];
    cout <<rank<<"  Received!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";


    HYBRID_MESH *originalMesh=nullptr;
    HYBRID_MESH *refinedMesh=nullptr;
    if(size==1)
    {
        originalMesh = &tetrasfile;
        refinedMesh = &newTetrasfile;
    }
    else
    {
        originalMesh = &construcTetras;
        refinedMesh = &newTetrasfile;
    }

    for (int i = 0; i < refineTimes; i++)
    {

        if(refinedMesh!=nullptr)
            refinedMesh->clear();
        cout<<"originalMesh->NumNodes"<<originalMesh->NumNodes<<endl;
        if(arguments.gmName!=nullptr)
        meshRefining(*originalMesh, *refinedMesh, rank,ref);
        else
        meshRefining(*originalMesh, *refinedMesh, rank);


        unifyBoundaries(*originalMesh, *refinedMesh, MPI_COMM_WORLD);
        if (rank == 0)
            printf("The %d-th refining finished.\n", i + 1);

        swap(originalMesh, refinedMesh);
    }

    swap(originalMesh, refinedMesh);

    //optimization procedure

    //
    if(refineTimes==0)
    {

    }

    if(originalMesh!=nullptr)
    {
        originalMesh->clear();
    }
    //tetra  optimalization


#ifdef VTK_OUTPUT
    char outputName[256];
    // /*

    sprintf(outputName,"%s%d%s",filename,rank,".vtk");
//	writeVTKFile(outputName, *originalMesh);

    sprintf(outputName,"%s%d%s",filename,rank,"new.vtk");
    writeVTKFile(outputName, *refinedMesh);

    sprintf(outputName,"%s%d%s",filename,rank,"surface_origin.vtk");
//	writeTriangleVTKFile(outputName, *originalMesh);

    sprintf(outputName,"%s%d%s",filename,rank,"surface_new.vtk");
    writeTriangleVTKFile(outputName, *refinedMesh);

#endif


    if (meshParts != nullptr)
        delete [] meshParts;
    //mesh quality improving


    //end improving


    int nParts = nOutputParts / size;
    if (rank < nOutputParts % size)
        ++nParts;



    cout<<"New nParts: "<<nParts<<endl;



    if(size==nOutputParts)
    {
        meshParts=refinedMesh;
        refinedMesh=nullptr;
        //refinedMesh
        cout<<"Skip the second partition!"<<endl;//nothing
    }

    else
    {
        meshParts = new HYBRID_MESH [nParts];
        cout<<"Second Partition!"<<endl;

        partition(*refinedMesh, meshParts, nParts, size, rank,1);

        cout<<"Second Partition finished!"<<endl;
        refinedMesh->clear();
        unifyRepartitionedFacets(meshParts, nParts, MPI_COMM_WORLD);


    }
    //sort the index of points
    for(int i=0;i<nParts;i++)
    {
        sortPointsID(meshParts[i]);

    }

    //
#ifdef VTK_OUTPUT
    for (int i = 0; i < nParts; ++i)
    {

        sprintf(outputName,"%s%d%s",filename,i * size + rank,"final.vtk");
        cout<<filename<<endl;
        writeVTKFile(outputName, meshParts[i]);
        sprintf(outputName,"%s%d%s",filename,i * size + rank, "surface_final.vtk");
        writeTriangleVTKFile(outputName, meshParts[i]);
        cout<<"output vtk..."<<endl;
    }
#endif

   // exit(0);

    /************************************************************************/
    /* 输出OpenFoam格式网格                                                   */
    /************************************************************************/
#ifdef DEBUG
    cout << "rank = " << rank << ": Begin OpenFoam Module" <<endl;

#endif //DEBUG
    /* calc the number of internal face & calc the number of cell in one part mesh */

    int *pNumInternalFacesSubMesh = new int[nParts]();  // the total internal face number in one part mesh.
    int *pNumInternalFacesProcMesh = new int[size]();  // the total internal face number in One Processor.
    int nInternalFacesProcMesh = 0;
    idx_t nFace = 0, nInternalFaceBegIdx = 0, nCountFace = 0;

    int *pNumCellsSubMesh = new int[nParts]();  // the total cell number in one part mesh.
    int *pNumCellsProcMesh = new int[size]();  // the total cell number in One Processor.
    int nTetr, nPrsm, nPyra, nHexa, nTria, nQuad;
    int nCellProcMesh = 0;
    idx_t nCell = 0, nCellBegIdx = 0, nCountCell = 0;

    for(int i = 0; i < nParts; ++i)
    {
        nTetr = meshParts[i].NumTetras;
        nPrsm = meshParts[i].NumPrsm;
        nPyra = meshParts[i].NumPyra;
        nHexa = meshParts[i].NumHexes;

        nTria = meshParts[i].NumTris;
        nQuad = meshParts[i].NumQuads;

        pNumInternalFacesSubMesh[i] = (4*nTetr + 5*nPrsm + 5*nPyra + 6*nHexa - nTria - nQuad)/2;
        nInternalFacesProcMesh += pNumInternalFacesSubMesh[i];

        pNumCellsSubMesh[i] = nTetr + nPrsm + nPyra + nHexa;
        nCellProcMesh += pNumCellsSubMesh[i];
    }
    //cout << "Proc: " << rank << " nCellProcMesh = " << nCellProcMesh << endl;

    MPI_Allgather(&nInternalFacesProcMesh, 1, MPI_INT,
        pNumInternalFacesProcMesh, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&nCellProcMesh, 1, MPI_INT,
        pNumCellsProcMesh, 1, MPI_INT, MPI_COMM_WORLD);

    nFace = 0;
    nCell = 0;
    for (int i = 0; i < rank; ++i)
    {
        nFace += pNumInternalFacesProcMesh[i];
        nCell += pNumCellsProcMesh[i];
    }

    int nAllInternalFace = nFace;
    for (int i = rank; i < size; ++i)
    {
        nAllInternalFace += pNumInternalFacesProcMesh[i];
    }
    //cout << "Proc: " << rank << " nCell = " << nCell << endl;

    /* set boundary condition */
#ifndef _BC_FILE_
    BoundaryCondition bc(bcstring);
#else
    /*  read boundary condition file */
    BoundaryCondition bc(filenameBC);
#endif //_CGNS_

    nCountFace = 0;
    nCountCell = 0;
    for (int i = 0; i < nParts; ++i)
    {
#ifdef DEBUG
        cout << "rank = " << rank << ": Begin to write " << i*size+rank << "subMesh" <<endl;
#endif //DEBUG
#ifdef _INTERNAL_FACE_FIRST_
        nInternalFaceBegIdx = nFace + nCountFace;
#else
        //nInternalFaceBegIdx = meshParts[i].NumUniqueFacets + nFace + nCountFace;
        nInternalFaceBegIdx = meshParts[i].NumUniqueSurfFacets + meshParts[i].NumUniqueInterfFacets + nFace + nCountFace;
#endif
        nCountFace += pNumInternalFacesSubMesh[i];

        nCellBegIdx = nCell + nCountCell;
//        if(rank == 1)
//            cout << "Proc: " << rank << " , part " << i << " ncell = " << nCell << ", nCountCell = " << nCountCell << endl;
        nCountCell += pNumCellsSubMesh[i];

//        cout << "Proc: " << rank << " , part " << i << " nCellBegIdx = " << nCellBegIdx << endl;

#if 0
        if(rank == 1 && i == 1)
        {
            TRI *pTri_tmp = meshParts[i].pTris;
            for(int j = 0; j < meshParts[i].NumTris; ++j)
            {
                cout << "rank = " << rank << " Tri[" << j <<"].iOppoProc = " <<
                        pTri_tmp[j].iOppoProc << endl;
            }
        }
#endif
	
        OpenFoamFile outputOF(i*size+rank, &meshParts[i], nInternalFaceBegIdx, nCellBegIdx, &bc, "./", filename, nOutputParts, nAllInternalFace);
        outputOF.write_OpenFoam_Mesh();
    }

    delete [] pNumCellsProcMesh;
    delete [] pNumCellsSubMesh;
    delete [] pNumInternalFacesProcMesh;
    delete [] pNumInternalFacesSubMesh;
    if(size==nOutputParts)
    {
        meshParts=nullptr;
    }
    {
        delete [] meshParts;
    }
    delete mpiTypeAssistant;


    MPI_Finalize();

    if (rank == 0)
        cout<<"Finished!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    return 0;
}
#endif
int main(int argc, char *argv[]){
#ifdef main_1
    main1(argc, argv);
#endif

#ifdef main_2
    main2(argc, argv);
#endif
}
