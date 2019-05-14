///* ----------------------------------------------------------------------------
// * ----------------------------------------------------------------------------
// *
// * ��ά����ͬ��Delaunay���������� (�汾�ţ�1.0Belta)
// * 3D Isotropic Delaunay Mesh Generation (Version 1.0Belta)
// *
// * �½��� �й� �㽭��ѧ��������ѧ�����о�����
// * ��Ȩ����	  2012��12��30��
// * Chen Jianjun  Center for Engineering & Scientific Computation,
// * Zhejiang University, P. R. China
// * Copyright reserved, 2012, 12, 30
// *
// * ��ϵ��ʽ
// *   �绰��+86-571-87951883
// *   ���棺+86-571-87953167
// *   ���䣺chenjj@zju.edu.cn
// * For further information, please conctact
// *  Tel: +86-571-87951883
// *  Fax: +86-571-87953167
// * Mail: chenjj@zju.edu.cn
// *
// * ------------------------------------------------------------------------------
// * ------------------------------------------------------------------------------*/
//#include <stdio.h>
//#include <assert.h>
//#include "iso3d_error.h"
//#include "iso3d_utility.h"
//#include "iso3d.h"
//#ifdef _TIMING_PERFORMANCE
//#include "spr.h"
//#include "time_perform.h"
//#endif
//#include "locsmoothing.h"
//#include<iostream>
//#include"smooth.h"
//void printHead();
//void printTimingProfile();
//void printPrompt();
//int dtiso3d_init(int argc, char **argv);
//int dtiso3d_readVlmMesh(char chExampleName[], DTIso3D& iso3d);
//int dtiso3d_improving(char chExampleName[], DTIso3D& iso3d);
//int dtiso3d_meshing(char chExampleName[], DTIso3D& iso3d);

//static bool g_bMeshingEnabled = false;  /* �������ɹ����Ƿ������� -M */
//static bool g_bQuaImprEnabled = false;	/* �����Ż������Ƿ������� -Q */

//#define MAX_FILE_LEN 1024


//int dtiso3d_meshing(char chExampleName[], DTIso3D& iso3d)
//{
//#ifdef _ERROR_CHK
//    int nErrCd;
//#endif /* _ERROR_CHK */
//    char strPLS[MAX_FILE_LEN], strBA3[MAX_FILE_LEN], strPL3[MAX_FILE_LEN];
//#ifdef _NGB
//    char strNGB[MAX_FILE_LEN];
//#endif

//    strcpy(strPLS, chExampleName);
//    strcat(strPLS, ".pls");
//    strcpy(strBA3, chExampleName);
//    strcat(strBA3, ".ba3");
//    strcpy(strPL3, chExampleName);
//    strcat(strPL3, ".pl3");
//#ifdef _NGB
//    strcpy(strNGB, chExampleName);
//    strcat(strNGB, ".init.ngb");
//#endif /* NGB */

//    printf("Example Name: %s\n", chExampleName);

//#ifdef _TIMING_PERFORMANCE
//    time_start[READ_AND_INIT_TIME] = SPRlogTime();
//#endif
//    /*
//     * initialization
//     */
//    if (iso3d.readBA3(strBA3) <= 0)
//    {
//        printf("fail to read BA3 file %s\n", strBA3);
//    }

//    if (iso3d.readPLS(strPLS) <= 0)
//    {
//        printf("fail to read PLS file %s\n", strPLS);
//        exit(1);
//    }

////	if (iso3d.isSelfIntersect())
////	{
////		printf("Self-intersection found in the surface. \n");
////		return 1;
////	}

//#ifdef _TIMING_PERFORMANCE
//    time_final[READ_AND_INIT_TIME] = time_start[BOUND_POINT_INST_TIME] = SPRlogTime();
//#endif
//    /*
//     * insert boundary nodes
//     */
//    iso3d.bndPntInst();

////	iso3d.readTetgenELE("mohne-a_cf-intmediate.ele", true);
////	iso3d.readTetgenELE("mohne-half-part-intmediate.ele", true);

////	if (iso3d.isSelfIntersect())
////	{
////		printf("Self-intersection found in the surface. \n");
////		return 1;
////	}

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData();
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK*/


//#ifdef _TIMING_PERFORMANCE
//    time_final[BOUND_POINT_INST_TIME] = SPRlogTime();
//#ifdef _RECOVERY_FINAL
//    time_start[INNER_POINT_INST_TIME] = time_final[BOUND_POINT_INST_TIME];
//#endif
//#endif

//#ifndef _NO_FIELD_POINTS
//#ifdef _RECOVERY_FINAL
//    /*
//     * insert inner points
//     */
//    iso3d.innerPntInst();
//#endif /* _RECOVERY_FINAL */

//#ifdef _TIMING_PERFORMANCE
//#ifdef _RECOVERY_FINAL
//    time_final[INNER_POINT_INST_TIME] = SPRlogTime();
//#endif /* _RECOVERY_FINAL */
//#endif /* _TIMING_PERFORMANCE */
//#endif /* #ifndef _NO_FIELD_POINTS */

//#ifdef _RECOVER_WITHOUT_SPS

//#ifdef _TIMING_PERFORMANCE
//    time_start[LSSP_RECV_EDGE_TIME] = SPRlogTime();
//#endif /* _TIMING_PERFORMANCE */

//    for (int i = 0; i < 0; i++)
//    {
//        iso3d.checkEdgeSP();
//        iso3d.recvEdges_LessStnPnt();
//    }

////	iso3d.writePL3("initmesh.pl3");
////	iso3d.writeTetgenELE("mohne-half-part-intmediate.1.ele");
////	exit(1);

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData();
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK*/

//#ifdef _TIMING_PERFORMANCE
//    time_final[LSSP_RECV_EDGE_TIME] = time_start[NOSP_RECV_EDGE_TIME] = SPRlogTime();
//#endif /* _TIMING_PERFORMANCE */

////	iso3d.checkEdgeSP();
////	exit(1);

//#if 1

//    iso3d.setBndProtectionFlag(true);
//    iso3d.abstNodFirEle();	/* ��ȡ�׵�Ԫ��������ϣ�� */
//    for (int i = 0; i < 1; i++)
//    {
//        iso3d.recvBndEdge_NoStnPnt_V2();
//        iso3d.recvBndFace_NoStnPnt_V2();
//    }

//    for (int i = 0; i < 0; i++)
//    {
//        iso3d.classifyBndEdges(NULL);
//        iso3d.recvBndEdge_NoStnPnt();
//        iso3d.checkEdgeSP();
//    }
////	iso3d.writeTetgenELE("cami1a-intmediate.1.ele");
////	exit(1);

//#if 0
//    iso3d.setBndProtectionFlag(true);
//    iso3d.classifyBndEdges(NULL);
//    for (int i = 0; i < 0; i++)
//    {
//        iso3d.recvBndEdge_NoStnPnt();
//        iso3d.checkEdgeSP();
//    }
//#endif

////	iso3d.writeTetgenELE("mohne_a-intmediate.1.ele");
//#else
//    iso3d.recvEdges_NoStnPnt();
//#endif
////	exit(1);

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData();
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK*/

//#ifdef _TIMING_PERFORMANCE
//    time_final[NOSP_RECV_EDGE_TIME] = time_start[NOSP_RECV_FACE_TIME] = SPRlogTime();
//#endif /* _TIMING_PERFORMANCE */

//#if 0
//    iso3d.checkEdgeSP();

//    iso3d.recvFaces_NoStnPnt();


//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData();
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK*/
//#endif

//#ifdef _TIMING_PERFORMANCE
//    time_final[NOSP_RECV_FACE_TIME] = SPRlogTime();
//#endif /* _TIMING_PERFORMANCE */
//#endif /* _RECOVER_WITHOUT_SPS */


//#ifdef _TIMING_PERFORMANCE
//    time_start[CONF_RECV_EDGE_TIME] = SPRlogTime();
//#endif /* _TIMING_PERFORMANCE */

////	iso3d.checkEdgeSP();
////	iso3d.writeTetgenELE("mohne-half-part-intmediate.ele");
//    //exit(1);

//    /*
//     * recover boundary edges
//     */
//    iso3d.recvEdges();

////	if (iso3d.checkConfEdges() != 0)
////		exit(1);

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData();
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK*/

//#ifdef _TIMING_PERFORMANCE
//    time_final[CONF_RECV_EDGE_TIME] = time_start[CONF_RECV_FACE_TIME] = SPRlogTime();
//#endif

//    /*
//     * recovery boudary faces
//     */
//    iso3d.recvFaces();

// // if (iso3d.checkOuIn(true) != 0)
// // {
// //		exit(1);
// //	}

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData();
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK*/

//#ifdef _TIMING_PERFORMANCE
//    time_final[CONF_RECV_FACE_TIME] = time_start[TYPE_DEL_ELEM_TIME] = SPRlogTime();
//#endif

//    /*
//     * remove outer elements
//     */
//    iso3d.typeEles();

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//#ifdef _TIMING_PERFORMANCE
//    time_final[TYPE_DEL_ELEM_TIME] = SPRlogTime();
//#endif

////	iso3d.printConformalRecvInfo();

//#ifdef _CONSTRAINED_RECOVERY

//#ifdef _TIMING_PERFORMANCE
//    time_start[CONS_RECV_EDGE_TIME] = SPRlogTime();
//#endif
//    iso3d.prepareRecvDesBnds();
//    iso3d.recvDesEdges();

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//#ifdef _TIMING_PERFORMANCE
//    time_final[CONS_RECV_EDGE_TIME] = time_start[CONS_RECV_FACE_TIME] = SPRlogTime();
//#endif

//    iso3d.recvDesFacets();

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//#ifdef _TIMING_PERFORMANCE
//    time_final[CONS_RECV_FACE_TIME] = time_start[CONS_RECV_PNT_RMV_TIME] = SPRlogTime();
//#endif

//    for (int i = 0; i < 1; i++)
//    {
//        iso3d.removeInnerSteinerPoint_SPR();	/* ���ÿ�ǻ�����㷨����Steiner�� */
//        iso3d.removeInnerSteinerPoint_Flip();
//    }

// //	iso3d.removeEdgeSteinerPoint_SPR();
// //	iso3d.removeFacSteinerPoint_SPR();
// //	iso3d.removeAllLeftNodes();
//    iso3d.printSteinerPntInfo();

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//#ifdef _TIMING_PERFORMANCE
//    time_final[CONS_RECV_PNT_RMV_TIME] = SPRlogTime();
//#endif

//#endif /* _CONSTRAINED_RECOVERY */

//    iso3d.typeDelEles();

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//#ifndef _NO_FIELD_POINTS
//#ifndef _RECOVERY_FINAL
//#ifdef _TIMING_PERFORMANCE
//    time_start[INNER_POINT_INST_TIME] = SPRlogTime();
//#endif

//    iso3d.innerPntInst();

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//#ifdef _CONSTRAINED_RECOVERY
//#ifdef _TIMING_PERFORMANCE
//    time_final[INNER_POINT_INST_TIME] = SPRlogTime();
//#endif

//    /* �ڲ����ڲ����󣬿ɳ�����ɾ������Steiner�� */
//    iso3d.prepareRmvSteinerPntAftInstInnNode();
//    iso3d.removeInnerSteinerPoint_SPR();	/* ���ÿ�ǻ�����㷨����Steiner�� */
//    iso3d.removeEdgeSteinerPoint_SPR(false);
//    iso3d.removeFacSteinerPoint_SPR(false);
// //	iso3d.removeAllLeftNodes();
//    iso3d.printSteinerPntInfo();

//#ifdef _TIMING_PERFORMANCE
//    time_final[CONS_RECV_PNT_RMV_TIME] += SPRlogTime() - time_final[INNER_POINT_INST_TIME];
//#endif
//#endif /* _CONSTRAINED_RECOVERY */

//#endif /* #ifndef _RECOVERY_FINAL */
//#endif/* _NO_FIELD_POINTS */


//#ifdef _TIMING_PERFORMANCE
//    time_start[WRIT_AND_FINA_TIME] = SPRlogTime();
//#endif

//    iso3d.rmvNodsAndEles();
//#ifdef _CONSTRAINED_RECOVERY
//    iso3d.fillParents(false);
//#else
//    iso3d.fillParents(true);
//#endif

//    strcpy(strPL3, chExampleName);
//    strcat(strPL3, ".init.pl3");
//    iso3d.writePL3(strPL3);
////	iso3d.writeNGB(strNGB);
////	exit(1);

//#ifdef _NGB
//    printf("\nPrepare write NGB %s, please wait ...\n", strNGB);
//    iso3d.writeNGB(strNGB);
//#endif /* _NGB */
//#ifdef _TIMING_PERFORMANCE
//    time_final[WRIT_AND_FINA_TIME] = SPRlogTime();
//#endif

//    return 1;
//}

//int dtiso3d_improving(char chExampleName[], DTIso3D& iso3d)
//{
//#ifdef _ERROR_CHK
//    int nErrCd;
//#endif /* _ERROR_CHK */
//    char strPL3[MAX_FILE_LEN], strQual[MAX_FILE_LEN];
//#ifdef _NGB
//    char strNGB[MAX_FILE_LEN];
//#endif

//#ifdef _TIMING_PERFORMANCE
//    time_start[QUAL_IMPRV_TIME] = SPRlogTime();
//#endif
//    /*added*/
//    //iso3d.printResult();
//    //printTimingProfile();
//#if 0
//    //�Ż�
//    iso3d.updatebdnpnts();
//    iso3d.initoptimize();
////	iso3d.setsmoothpar(argv[1]);
//#endif

//    strcpy(strQual, chExampleName);
//    strcat(strQual, ".init.allangles.qual");
//    iso3d.printElemQuality(strQual, QUALALLANGLES, 0.0, false);

//    iso3d.meshopt();

//#ifdef _TIMING_PERFORMANCE
//    time_final[QUAL_IMPRV_TIME] = SPRlogTime();
//#endif

//    strcpy(strQual, chExampleName);
//    strcat(strQual, ".allangles.qual");
//    iso3d.printElemQuality(strQual, QUALALLANGLES, 0.0, false);

//    strcpy(strQual, chExampleName);
//    strcat(strQual, ".sineangle.qual");
//    iso3d.printElemQuality(strQual, QUALSINEANGLE, 30, false);

//    strcpy(strQual, chExampleName);
//    strcat(strQual, ".voledge.qual");
//    iso3d.printElemQuality(strQual, QUALVLRMS3RATIO, 0.3, false);
////
//// 	iso3d.shelltransform(sin(ANGLE(30.0)));
//// 	//iso3d.tranverseST(1);
//// 	iso3d.smooth(8);
////
//// 	iso3d.checkbadcells();
//// 	iso3d.smooth(8);
////
//// 	iso3d.shelltransform(sin(ANGLE(30.0)));
//// 	//iso3d.tranverseST(1);
//// 	iso3d.smooth(8);

//// 	iso3d.checkbadcells();
//// 	iso3d.smooth(8);

//    strcpy(strQual, chExampleName);
//    strcat(strQual, ".bad_voledge.pl3");
//    iso3d.printBadCells(strQual, QUALVLRMS3RATIO, 0.1);

//    iso3d.rmvNodsAndEles();
//    iso3d.fillParents();

////	iso3d.evaluateQuality();
////	iso3d.printbadcells();

//    strcpy(strPL3, chExampleName);
////	strcat(strPL3, ".pl3");
//    strcat(strPL3, ".improved.pl3");
//    iso3d.writePL3(strPL3);

//    return 1;
//}

//int dtiso3d_readVlmMesh(char chExampleName[], DTIso3D& iso3d)
//{
//    char strPL3[MAX_FILE_LEN];

//    strcpy(strPL3, chExampleName);
//    strcat(strPL3, ".pl3");

//    return iso3d.readPL3(strPL3);
//}

////#ifndef _ERROR_CHK
////#define _ERROR_CHK
////#endif
//int main2(int argc, char **argv)
//{
//    DTIso3D iso3d;
//#ifdef _ERROR_CHK
//    int nErrCd;
//#endif /* _ERROR_CHK */
//    char strPLS[MAX_FILE_LEN], strBA3[MAX_FILE_LEN], strPL3[MAX_FILE_LEN], strQual[MAX_FILE_LEN];
//#ifdef _NGB
//    char strNGB[MAX_FILE_LEN];
//#endif

//    /*
//     * ��ȡ��������
//     */
//    GEOM_FUNC::exactinit();
//    /*
//     * ����SPR�����Ĳ���
//     */
//    spr_time_limit = 1.0; /* ʱ�����ƣ�=0.0��ʾ������ */

//#if 0
//    int iNod, iElem;
//    int elems[1][4] = {{0, 1, 2, 3}};
//    REAL qual, quals[1], qualgrads[6][3], volume, volumegrad[3], qualmeasure = 5;

//    num_vertices = 4;
//    vertices[0][0] = 0.021635;
//    vertices[0][1] = -0.036032;
//    vertices[0][2] = 0.003098;
//    vertices[1][0] = 0.021751;
//    vertices[1][1] = -0.036169;
//    vertices[1][2] = 0.003115;
//    vertices[2][0] = 0.021656;
//    vertices[2][1] = -0.035976;
//    vertices[2][2] = 0.003141;
//    vertices[3][0] = 0.021518;
//    vertices[3][1] = -0.035895;
//    vertices[3][2] = 0.003080;
//    quals[0] = tetquality(vertices[0], vertices[1], vertices[3], vertices[2], qualmeasure);

//    ::setMesh_LocSmoothing(elems, 1, quals);
//    getoptinfo(0, 0, &qual, qualgrads, &volume, volumegrad, qualmeasure);
//#endif

//    printHead();

//    if (dtiso3d_init(argc, argv) == 0)
//    {
//        printf("Cannot initialize. Terminate\n");
//        exit(1);
//    }

//    if (strlen(argv[1]) + 5 > MAX_FILE_LEN)
//    {
//        printf("too long file name %s", argv[1]);
//        exit(1);
//    }

// //   SPRlogTime_Init();
//#ifdef _TIMING_PERFORMANCE
//    time_start[TOTAL_EXEC_TIME] =  SPRlogTime();
//#endif

//    if (g_bMeshingEnabled)
//        dtiso3d_meshing(argv[1], iso3d);

//    if (g_bQuaImprEnabled)
//    {
//        if (!g_bMeshingEnabled)
//        {
//#ifdef _TIMING_PERFORMANCE
//            time_start[READ_AND_INIT_TIME] = SPRlogTime();
//#endif
//            if (dtiso3d_readVlmMesh(argv[1], iso3d) != 1)
//            {
//                printf("Cannot read the volume mesh.\n");
//                exit(1);
//            }
//#ifdef _TIMING_PERFORMANCE
//            time_final[READ_AND_INIT_TIME] = SPRlogTime();
//#endif
//            /* ������ϣ����ע�⣺Ŀǰֻ���Ƕ����������⣻�Ƕ���������������Ҫ��չ */
//            iso3d.abstNodFirEle();
//        }
//        dtiso3d_improving(argv[1], iso3d);
//    }


//#ifdef _TIMING_PERFORMANCE
//    time_final[TOTAL_EXEC_TIME] = SPRlogTime();
//#endif

//#ifdef _ERROR_CHK
//    nErrCd = iso3d.checkGlobalData(false);
//    if (nErrCd != 0)
//    {
//        printError(nErrCd);
//        exit(1);
//    }
//#endif /* _ERROR_CHK */

//    iso3d.printResult();
//    printTimingProfile();
//    return 0;
//}

///*
// * ����������ͷ��Ϣ
// * print the head info. of the program (PDMG)
// */
//void printHead()
//{
//    printf(" ********************************************************************\n");
//    printf("               3D Isotropic Delaunay Mesh Generation \n");
//    printf("                           Version 1.0Belta\n\n");
//    printf("         Center for Engineering & Scientific Computation\n");
//    printf("                Zhejiang University, P. R. China\n");
//    printf("	      CHEN Jianjun   Copyright reserved, 27/11/2014\n\n");
//    printf("               For further information, please contact Dr. CHEN JJ\n");
//    printf("                         Tel: +86-571-87951883\n");
//    printf("                         Fax: +86-571-87953167\n");
//    printf("                        Mail: chenjj@zju.edu.cn\n");
//    printf(" ********************************************************************\n\n\n");
//}

///*
// * ������Ҫ�Ľ�����Ϣ
// * print the resulting info. of the program (PDMG)
// */
//void printTimingProfile()
//{
//#ifdef _TIMING_PERFORMANCE
//    double totalTime = 0.0;
//    double percentFactor = 0.0, computingRatio;
//    int totalECCallings, succECCallings = 0, totalECCallingsSub = 0, succECCallingsSub = 0, i;

//    for (i = TIME_START_IDX; i <= TYPE_DEL_ELEM_TIME; i++)
//        time_elaps[i] = time_final[i] - time_start[i];
//    totalTime = time_elaps[TOTAL_EXEC_TIME];

////	if (totalTime <= 0.0)
////		return;

//    percentFactor = 100.0/totalTime;

//    printf(" ********************************************************************\n");
//    printf("                            Timing Profile\n");
//    printf(" ********************************************************************\n");
//    printf("%-24s%-12s%-12s\n", "Step Name", "Time(s)", "Percentage(%)");
//    printf("%-24s%-12.2f%-12.2f\n", "Read & Init Input", time_elaps[READ_AND_INIT_TIME], percentFactor*time_elaps[READ_AND_INIT_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Insert Boun. Pnts", time_elaps[BOUND_POINT_INST_TIME], percentFactor*time_elaps[BOUND_POINT_INST_TIME]);
//#ifndef _NO_FIELD_POINTS
//#ifdef _RECOVERY_FINAL
//    printf("%-24s%-12.2f%-12.2f\n", "Insert Innr. Pnts", time_elaps[INNER_POINT_INST_TIME], percentFactor*time_elaps[INNER_POINT_INST_TIME]);
//#endif
//#endif /* _NO_FIELD_POINTS */
//#ifdef _RECOVER_WITHOUT_SPS
//    printf("%-24s%-12.2f%-12.2f\n", "LsSP. Recv. Edges", time_elaps[LSSP_RECV_EDGE_TIME], percentFactor*time_elaps[LSSP_RECV_EDGE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "NoSP. Recv. Edges", time_elaps[NOSP_RECV_EDGE_TIME], percentFactor*time_elaps[NOSP_RECV_EDGE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "NoSP. Recv. Faces", time_elaps[NOSP_RECV_FACE_TIME], percentFactor*time_elaps[NOSP_RECV_FACE_TIME]);
//#endif /* __RECOVER_WITHOUT_SPS */
//    printf("%-24s%-12.2f%-12.2f\n", "Conf. Recv. Edges", time_elaps[CONF_RECV_EDGE_TIME], percentFactor*time_elaps[CONF_RECV_EDGE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Conf. Recv. Faces", time_elaps[CONF_RECV_FACE_TIME], percentFactor*time_elaps[CONF_RECV_FACE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Dele. Outer Elems.", time_elaps[TYPE_DEL_ELEM_TIME], percentFactor*time_elaps[TYPE_DEL_ELEM_TIME]);
//#ifdef _CONSTRAINED_RECOVERY
//    printf("%-24s%-12.2f%-12.2f\n", "Cons. Recv. Edges", time_elaps[CONS_RECV_EDGE_TIME], percentFactor*time_elaps[CONS_RECV_EDGE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Cons. Recv. Faces", time_elaps[CONS_RECV_FACE_TIME], percentFactor*time_elaps[CONS_RECV_FACE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Rmv. Steiner Pnts", time_elaps[CONS_RECV_PNT_RMV_TIME], percentFactor*time_elaps[CONS_RECV_PNT_RMV_TIME]);
//#endif
//#ifndef _NO_FIELD_POINTS
//#ifndef _RECOVERY_FINAL
//    printf("%-24s%-12.2f%-12.2f\n", "Insert Innr. Pnts", time_elaps[INNER_POINT_INST_TIME], percentFactor*time_elaps[INNER_POINT_INST_TIME]);
//#endif
//#endif /* _NO_FIELD_POINTS */
//    printf("%-24s%-12.2f%-12.2f\n", "Finish and Output", time_elaps[WRIT_AND_FINA_TIME], percentFactor*time_elaps[WRIT_AND_FINA_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Recover Edge Loop", time_elaps[RECV_EDG_LOOP_TIME], percentFactor*time_elaps[RECV_EDG_LOOP_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Recover Face Loop", time_elaps[RECV_FAC_LOOP_TIME], percentFactor*time_elaps[RECV_FAC_LOOP_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Remesh EDG. Shell", time_elaps[REMESH_ASHELL_TIME], percentFactor*time_elaps[REMESH_ASHELL_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Remove Edge Recv.", time_elaps[REMOV_AN_EDGE_TIME], percentFactor*time_elaps[REMOV_AN_EDGE_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Reflex Edge Check", time_elaps[REFLEX_EDG_CH_TIME], percentFactor*time_elaps[REFLEX_EDG_CH_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Search ALL Shells", time_elaps[FIND_ALLSHELL_TIME], percentFactor*time_elaps[FIND_ALLSHELL_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Dyn. Rm. Sh. Main", time_elaps[DYN_RSHL_MAIN_TIME], percentFactor*time_elaps[DYN_RSHL_MAIN_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Chk. ED. Recovery", time_elaps[CHK_EDGE_RECV_TIME], percentFactor*time_elaps[CHK_EDGE_RECV_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Chk. FC. Recovery", time_elaps[CHK_FACE_RECV_TIME], percentFactor*time_elaps[CHK_FACE_RECV_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Qual. Improvement", time_elaps[QUAL_IMPRV_TIME], percentFactor*time_elaps[QUAL_IMPRV_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Qual. Recurs. ST.", time_elaps[QUAL_RECV_SH_TRANS_TIME], percentFactor*time_elaps[QUAL_RECV_SH_TRANS_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Opt. Mesh Smthing", time_elaps[MESH_SMOOTHING_TIME], percentFactor*time_elaps[MESH_SMOOTHING_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Edge Contraction ", time_elaps[EDGE_CONTRACT_TIME], percentFactor*time_elaps[EDGE_CONTRACT_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Edge Splitting   ", time_elaps[EDGE_SPLITTING_TIME], percentFactor*time_elaps[EDGE_SPLITTING_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Qual. Evaluation ", time_elaps[QUALITY_EVAL_TIME], percentFactor*time_elaps[QUALITY_EVAL_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Total Execu. Time", time_elaps[TOTAL_EXEC_TIME], 100.0);
//#endif /* _TIMING_PERFORMANCE */
//#ifdef _VERBOSE
//#ifdef _ROBUST_BW_KERNEL
//    printf("%-24s%-12.2f%-12.2f\n", "Form Init. Cavity", time_elaps[FORM_INIT_CAVITY_TIME], percentFactor*time_elaps[FORM_INIT_CAVITY_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Modify the Cavity", time_elaps[MODIFY_CAVITY_TIME], percentFactor*time_elaps[MODIFY_CAVITY_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Update the Trian.", time_elaps[UPDATE_TRIANGULATION_TIME], percentFactor*time_elaps[UPDATE_TRIANGULATION_TIME]);
//    printf("%-24s%-12.2f%-12.2f\n", "Create Inner Pnt.", time_elaps[CREATE_INNER_POINT_TIME], percentFactor*time_elaps[CREATE_INNER_POINT_TIME]);
//#ifdef _INSIDE_TIMING_CALC
//    printf("%-24s%-12.2f%-12.2f\n", "Calc. Elem. Para.", time_elaps[CREATE_ELEM_PARA_TIME], percentFactor*time_elaps[CREATE_ELEM_PARA_TIME]);
//#endif /* _INSIDE_TIMING_CALC */
//#endif /* _ROBUST_BW_KERNEL */
//#endif /* _VERBOSE */
//    printf("%-24s%-12.2f%-12.2f\n", "Call Poly. Recon.", time_elaps[SPR_CALLING_TIME], percentFactor*time_elaps[SPR_CALLING_TIME]);
//#ifdef _VERBOSE
//    printf("%-24s%-12.2f%-12.2f\t%fsec. per node\n", "Cln. af. Suc. Ad.", time_elaps[ADD_SUCC_CLEAN_TIME],
//        percentFactor*time_elaps[ADD_SUCC_CLEAN_TIME],time_elaps[ADD_SUCC_CLEAN_TIME]/(double)g_nInnInstPoint);
//    printf("%-24s%-12.2f%-12.2f\t%fsec. per node\n", "Cln. af. FAI. Ad.", time_elaps[ADD_FAIL_CLEAN_TIME],
//        percentFactor*time_elaps[ADD_FAIL_CLEAN_TIME],time_elaps[ADD_FAIL_CLEAN_TIME]/(double)(g_nInnTryPoints - g_nInnInstPoint));
//#endif /* _VERBOSE */

//    printf("\nNo. of interior points that are attempted to be inserted and of which really inserted:: (%d,%d).\n", g_nInnTryPoints, g_nInnInstPoint);
//    printf("Average cavity size for the really inserted points: %f.\n", g_nInnInstPoint > 0 ? g_nInnTotalCaviElem /(double)g_nInnInstPoint : 0.0);
//    printf("#Callings of the SPR routine and of which the successful ones: (%d,%d).\n", g_nSPRCallings, g_nSPRCallSucc);
//    printf("Average Size of the reconnected polyhedra  & Average time cost per SPR calling. #face: %f, Time cost(s):%f.\n",
//        (double)g_nSPRPolySize/g_nSPRCallings, time_elaps[SPR_CALLING_TIME]/g_nSPRCallings);
//    if (g_bQuaImprEnabled)
//    {
//        printf("#Callings of the smoothing routine in the smoothing passes: %d\n", g_nGlobalSmoothingCallings);
//        printf("#Successful Callings: %d\n", g_nGlobalSmoothingCallSucc);
//        printf("#Callings of the smoothing routine in the edge contraction passes: %d\n", g_nEdContSmoothingCallings);
//        printf("#Successful Callings: %d\n", g_nEdContSmoothingCallSucc);
//        printf("#Callings of the smoothing routine in the edge splitting passes: %d\n", g_nEdSpltSmoothingCallings);
//        printf("#Successful Callings: %d\n", g_nEdSpltSmoothingCallSucc);

//        totalECCallings = succECCallings = 0;
//        for (i = 0; i < MAX_EC_QUALITY_RATIO_DIVIDE; i++)
//        {
//            succECCallings += g_nECQualityRatio_SuccStatistics[i];
//            totalECCallings += g_nECQualityRatio_TotalStatistics[i];
//        }
//        printf("Distribution of quality ratios for successful edge contraction callings.\n");
//        printf("#Succ. Callings/#Total Callings:%d/%d = %f%%.\n", succECCallings, totalECCallings,
//            totalECCallings > 0 ? 100.0*succECCallings/totalECCallings : 0.0);
//        printf("Min = %f; Max = %f.\n", g_dECQualityRatio_MinRatio, g_dECQualityRatio_MaxRatio);
//        printf("Distribution of the ratios.\n");
//        totalECCallingsSub = 0;
//        succECCallingsSub = 0;
//        for (i = 0; i < MAX_EC_QUALITY_RATIO_DIVIDE; i++)
//        {
//            totalECCallingsSub += g_nECQualityRatio_TotalStatistics[i];
//            succECCallingsSub += g_nECQualityRatio_SuccStatistics[i];
//            printf(" [%f--%f]:%d/%d\t= %4.1f%%; Score VS Load: %4.1f%% : %4.1f%%.\n",
//                g_dECQualityRatio_Divide[i], g_dECQualityRatio_Divide[i+1],
//                g_nECQualityRatio_SuccStatistics[i], g_nECQualityRatio_TotalStatistics[i],
//                g_nECQualityRatio_TotalStatistics[i] > 0 ?
//                100.0*g_nECQualityRatio_SuccStatistics[i]/g_nECQualityRatio_TotalStatistics[i] :0.0,
//                succECCallings > 0 ? 100.0*succECCallingsSub/succECCallings : 0.0,
//                totalECCallings > 0 ? 100.0*totalECCallingsSub/totalECCallings : 0.0);
//        }
//    }
//}

///*
// * ��������������������ʽҪ��
// * print the requirment for the command line
// */
//void printPrompt()
//{
//    printf("Please check input arguments!\n");
//    printf("Four options:\n");
//    printf("\n\tOption 1: DTIso3D example_name");
//    printf("\n\tOption 2: DTIso3D example_name -MQ\n");
//    printf("\tDTIso3D reads a surface mesh (example_name.pls) &\n\ta mesh size map (optional, example_name.ba3)\n");
//    printf("\tand then generates a volume mesh, and also improves \n\tthe volume mesh.\n");
//    printf("\n\tOption 3: DTIso3D example_name -M\n");
//    printf("\tDTIso3D reads a surface mesh (example_name.pls) &\n\ta mesh size map (optional, example_name.ba3)\n");
//    printf("\tand then generates a volume mesh, but does not \n\timprove the volume mesh.\n");
//    printf("\n\tOption 4: DTIso3D example_name -Q\n");
//    printf("\tDTIso3D reads example_name.pl3 (volume mesh) and \n\tthen improves the volume mesh.\n\n");
//}

//int dtiso3d_init(int argc, char **argv)
//{
//    int lenOfString = 0, iChar;
//    char* pChar = NULL, ch;

//    if (argc != 2 && argc != 3)
//    {
//        printPrompt();
//        return 0;
//    }

//    if (argc == 3)
//        {
//        pChar = argv[2];
//        lenOfString = strlen(argv[2]);
//        if (lenOfString <= 1 && *pChar != '-')
//        {
//            printPrompt();
//            return 0;
//        }
//        pChar++;

//        while (*pChar != '\0')
//        {
//            switch (*pChar)
//            {
//            case 'M':
//                g_bMeshingEnabled = true;
//                break;
//            case 'Q':
//                g_bQuaImprEnabled = true;
//                break;
//            default:
//                printPrompt();
//                return 0;
//            }
//            pChar++;
//        }
//    }
//    else
//    {
//        g_bMeshingEnabled = g_bQuaImprEnabled = true;
//    }

//    return 1;
//}
