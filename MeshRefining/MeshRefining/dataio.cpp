#include "dataio.h"
#include "stdio.h"
#include "ctype.h"
#include <string.h>
#include "stdlib.h"
#include "assert.h"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include "cgnslib.h"
#include <fstream>

using namespace std;
bool isNullOrComment(char* chLine)
{
    char *str = chLine;
    int numOfValidC = 0;
    while (*str != '\0' && numOfValidC == 0)
    {
        if (!isspace(*str))
        {
            if (numOfValidC == 0 && *str == '#')
                return true;
            numOfValidC++;
        }
        str++;
    }

    return numOfValidC == 0;
}


int readValidLine(FILE *fp, char *chLine, int len)
{
    if (fp && !feof(fp))
    {
        do
        {
            fgets(chLine, len, fp);
            if (!isNullOrComment(chLine))
                return 0; // succeed
        }
        while (!feof(fp));
    }

    return 1; // fail
}

int eraseWhiteSpace(char *chLine)
{
    char *strHd = chLine;
    char *strTl = chLine + strlen(chLine) - 1;
    int i = 0;

    while (strHd < strTl)
    {
        /* 锟接猴拷锟斤拷锟斤拷 */
        if (!isspace(*strHd))
            break;
        strHd++;
    }

    while (strHd < strTl)
    {
        /* 锟斤拷前锟斤拷锟斤拷 */
        if (!isspace(*strTl))
            break;
        strTl--;
    }

    if (strTl > strHd)
    {
        if (strHd > chLine)
        {
            i = 0;
            do
            {
                chLine[i++] = *(strHd++);
            }
            while (strHd <= strTl);
            chLine[i] = '\0';
        }
        else
            *(strTl+1) = '\0';
    }

    return 1; // fail
}

int setupCellNeig_test(int nNodes, int nElems, HYBRID_MESH *mesh)
{
    HEX *pHexes=mesh->pHexes;
    PRISM *pPrism=mesh->pPrisms;
    TETRAS *pTetras=mesh->pTetras;

    int   nHexes=mesh->NumHexes;
    int   nPrism=mesh->NumPrsm;
    int   nTetras=mesh->NumTetras;

    int nAllocFaceSize = nNodes * 10, nAllocNodeFaceSize = 20384;
    int nFaceSize = 0;
    int *vecRefIntFHash = NULL;
    InterFace * vecInterFaces = NULL, *vecInterFaces_Temp = NULL;
    int *ndIFaces = NULL, *ndIFaces_Temp = NULL;
    int errCode = 0;
    int cellIdx, faceAdIdx, lftCell, rgtCell;
    int facNdIdx1, facNdIdx2, facNdIdx3,facNdIdx4, minFacNdIdx, faceIt;
    int ndIFaceSize = 0;
    int i, j, k, nCommon;
    InterFace faceAd;

    int* pris[5];
    int pris1[3] = {0, 1, 2};
    int pris2[3] = {3, 4, 5};
    int pris3[4] = {0, 1, 4, 3};
    int pris4[4] = {1, 2, 5, 4};
    int pris5[4] = {0, 2, 5, 3};
    pris[0] = pris1;
    pris[1] = pris2;
    pris[2] = pris3;
    pris[3] = pris4;
    pris[4] = pris5;
    static int tet[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    static int hex[6][4] = {{0, 1, 5, 4}, {3, 2, 6 ,7}, {5, 6, 7, 4}, {1, 2, 3, 0}, {1, 2, 5, 6}, {0, 3, 7, 4}};


    int ndSize, clSize;
    //BKGElem *pElem = NULL, *pLft = NULL, *pRgt = NULL;

    HEX *pHexes_t=nullptr;
    PRISM *pPrism_t=nullptr;
    TETRAS *pTetras_t=nullptr;

    vecRefIntFHash = (int *) malloc(sizeof(int)*nNodes);                     //以最小点的值为下标，用面ID为值不断更新，每次更新都记录到face的next，故可以形成链状，最后记录的是链表的起点。
    vecInterFaces = (InterFace *) malloc(sizeof(InterFace)*nAllocFaceSize);  //储存着每个面的信息
    ndIFaces = (int *) malloc(sizeof(int)*nAllocNodeFaceSize);               //储存着面ID

    if (!vecRefIntFHash || !vecInterFaces || !ndIFaces)
    {
        errCode = -1;                                                //申请内存空间并验证是否成功
        goto FAIL;
    }

    ndSize = nNodes;
    clSize = nHexes+nPrism+nTetras;
    for (i = 0; i < ndSize; i++)
        vecRefIntFHash[i] = -1;
    for (cellIdx = 0; cellIdx < nPrism; cellIdx++)
    {

        pPrism_t = &pPrism[cellIdx];
        for (i = 0; i <= 4; i++)                 //初始化neighbors值为-1  BKG_MESH_DIM == 3
            pPrism_t->neighbors[i] = -1;
    }
        cout<<"first"<<endl;

    for (cellIdx = 0; cellIdx < nHexes; cellIdx++)
    {

        pHexes_t = &pHexes[cellIdx];
        for (i = 0; i <= 4; i++)                 //初始化neighbors值为-1  BKG_MESH_DIM == 3
            pHexes_t->neighbors[i] = -1;
    }
        cout<<"first"<<endl;
    for (cellIdx = 0; cellIdx <nTetras; cellIdx++)
    {

        pTetras_t = &pTetras[cellIdx];
        for (i = 0; i <= 4; i++)                 //初始化neighbors值为-1  BKG_MESH_DIM == 3
            pTetras_t->neighbors[i] = -1;
    }
        cout<<"first"<<endl;

    for (cellIdx = 0; cellIdx < clSize; cellIdx++)
    {
        if (cellIdx % 1000000 == 0 || cellIdx == clSize - 1)
            printf("%%%.2f.\n", 100.0*(cellIdx)/clSize);             //输出工作百分比（仅显示0%和99.99%）这个for内为主要工作循环
    if(cellIdx<nPrism){
        PRISM  *pLft=nullptr, *pRgt=nullptr;
        pPrism_t = &pPrism[cellIdx];
        for (i = 0; i <= 4; i++) {
            if(sizeof(pris[i])/sizeof(int)==4){
            facNdIdx1 = pPrism_t->vertices[pris[i][0]];
            facNdIdx2 = pPrism_t->vertices[pris[i][1]];                //BKG_MESH_DIM == 3  new BKG_MESH_DIM 应该为5
            facNdIdx3 = pPrism_t->vertices[pris[i][2]];
            facNdIdx4 = pPrism_t->vertices[pris[i][3]];

       int facevec[4]={facNdIdx1,facNdIdx2,facNdIdx3,facNdIdx4};
       sort(facevec,facevec+4);

            int min=nNodes+10;
            for(int index=0;index<4;++index) {
                 if(min>pPrism_t->vertices[pris[i][index]])
                 min=pPrism_t->vertices[pris[i][index]];
            }
            minFacNdIdx = min;


            j = 0;
            faceIt = vecRefIntFHash[minFacNdIdx];      //该点对应的neighbors为faceIt
            while (faceIt >= 0)                        //若该点(最小点)周围已经有面被录入容器
            {
                if (j >= nAllocNodeFaceSize)           //j是什么？     ndIFaces空间不足 重新分配
                {
                    nAllocNodeFaceSize += MAX_VALUE(100, (int)(nAllocNodeFaceSize * 0.1));
                    ndIFaces_Temp = (int*)realloc(ndIFaces, sizeof(int)*nAllocNodeFaceSize);
                    if (!ndIFaces_Temp)
                    {
                        errCode = -1;
                        goto FAIL;
                    }
                    ndIFaces = ndIFaces_Temp;
                }
                ndIFaces[j++] = faceIt;                //ndIFaces空间足够，则按顺序将已经录入的面id加入数组ndIFaces
                                                                                          //class InterFace
                                                                                          //{
                                                                                          //public:
                                                                                          //    int conn[BKG_MESH_DIM];  //面的点ID
                                                                                          //    int lftCell, rgtCell;    //左右邻接体
                                                                                          //    int hashNxt;             //
                                                                                          //};

                faceIt= vecInterFaces[faceIt].hashNxt;
            }
            ndIFaceSize = j;                      //j为已录入的数量
                                                           //j以及vecRefIntFHash的作用应该是为了加速查找重合面，否则遍历vecInterFaces即可
                                                           //vecRefIntFHash中存的面必然和当前面有联系
                                                           //

            nCommon = 0;
            for (j = 0;	j < ndIFaceSize; j++) {                 // j为暂时最小点周围已经加入容器的面的数量（非最终总数）对faceIt 也就是一个点周围一圈的面进行遍历找出重合面
                nCommon = 0;
                for (k = 0; k < 4; k++)
                    if (vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx1 ||
                            vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx2 ||
                            vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx3||
                        vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx4)
                        nCommon++;

                if (nCommon >= 4)
                {
            //        cout<<"find it !!!!!!!!!"<<endl;
                    break;
                }
            }

            if (nCommon < 4) {//未找到重合面，循环结束退出，说明该面尚未放入容器，该面尚未登记过lftCell。则进行登记操作
                /* 没锟斤拷锟揭碉拷 */
                if (nFaceSize >= nAllocFaceSize)                      //若空间不足则扩大空间
                {
                    nAllocFaceSize += MAX_VALUE(100, (int)(nAllocFaceSize * 0.1));
                    vecInterFaces_Temp = (InterFace*)realloc(vecInterFaces, sizeof(InterFace)*nAllocFaceSize);
                    if (!vecInterFaces_Temp)
                    {
                        errCode = -1;
                        goto FAIL;
                    }
                    vecInterFaces = vecInterFaces_Temp;
                }

                faceAd.conn[0] = facNdIdx1;          //面的conn点
                faceAd.conn[1] = facNdIdx2;
                faceAd.conn[2] = facNdIdx3;
                faceAd.conn[3] = facNdIdx4;
                faceAd.lftCell = cellIdx;           //当前体中的面毫无疑问有一个邻接体是当前体
                faceAd.rgtCell = -1;

                faceAdIdx = nFaceSize;
                faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];     //上一个faceAdIdx被放入此处，第一个faceAD的hasNxt为-1
                vecRefIntFHash[minFacNdIdx] = faceAdIdx;          //faceAdIdx被放入此处，给下一个faceAd作为next。
                vecInterFaces[nFaceSize++] = faceAd;              //faceAd顺序放入容器
            }



            else {

                if(!(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                     vecInterFaces[ndIFaces[j]].rgtCell < 0))
                    cout<<cellIdx<<endl;

                assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                        vecInterFaces[ndIFaces[j]].rgtCell < 0);


                vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;              //既然有重合面，说明该面已经登记过，说明该面已有lftCell
                                                                           //那么当前体就是该面的rgtCell。

                lftCell = vecInterFaces[ndIFaces[j]].lftCell;
                rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;
                pLft = &pPrism[lftCell];
                pRgt = &pPrism[rgtCell];                             //

             for(int k=0;k<6;++k){
                   if(pLft->neighbors[k]==-1)
                {
                    pLft->neighbors[k] = rgtCell;
                    pLft->neighborsmark[k]=double((facevec[0]+facevec[1]*2+facevec[2]*3+facevec[3]*4))/10;
                    break;
                }
             }

                pRgt->neighbors[i] = lftCell;
                pRgt->neighborsmark[i]=double((facevec[0]+facevec[1]*2+facevec[2]*3+facevec[3]*4))/10;

                                                            //在两个相邻的体的neighbors中互相储存对方。（下标为什么为i和k？）
                                              //i在此循环中代表该体的当前面，k则代表当前面的对点编号
            }                                                         //由cf[4][3]可以看出，此时k恒等于i
                                                                    //此时对于rgtCell来讲，自然是能保证neighbors[0]不被覆盖
                                                                    //lftCell如何保证？
        }

            else
            {


                facNdIdx1 = pPrism_t->vertices[pris[i][0]];
                facNdIdx2 = pPrism_t->vertices[pris[i][1]];                //BKG_MESH_DIM == 3  new BKG_MESH_DIM 应该为5
                facNdIdx3 = pPrism_t->vertices[pris[i][2]];

                int facevec[3]={facNdIdx1,facNdIdx2,facNdIdx3};
                sort(facevec,facevec+3);

                     int min=nNodes+10;
                     for(int index=0;index<3;++index) {
                          if(min>pPrism_t->vertices[pris[i][index]])
                          min=pPrism_t->vertices[pris[i][index]];
                     }
                     minFacNdIdx = min;

                j = 0;
                faceIt = vecRefIntFHash[minFacNdIdx];
                while (faceIt >= 0)                        //若该点(最小点)已经有对应的Neighbors（面的neighbors储存在对角点）
                {
                    if (j >= nAllocNodeFaceSize)           //j是什么？     ndIFaces空间不足 重新分配
                    {
                        nAllocNodeFaceSize += MAX_VALUE(100, (int)(nAllocNodeFaceSize * 0.1));
                        ndIFaces_Temp = (int*)realloc(ndIFaces, sizeof(int)*nAllocNodeFaceSize);
                        if (!ndIFaces_Temp)
                        {
                            errCode = -1;
                            goto FAIL;
                        }
                        ndIFaces = ndIFaces_Temp;
                    }
                    ndIFaces[j++] = faceIt;                //ndIFaces空间足够，则按顺序将
                                                                                              //class InterFace
                                                                                              //{
                                                                                              //public:
                                                                                              //    int conn[BKG_MESH_DIM];  //面的点ID
                                                                                              //    int lftCell, rgtCell;    //左右邻接体
                                                                                              //    int hashNxt;             //
                                                                                              //};

                    faceIt= vecInterFaces[faceIt].hashNxt;  //neighbors的hasnext的意义？
                }
                ndIFaceSize = j;
                                                               //j以及vecRefIntFHash的作用应该是为了加速查找重合面，否则遍历vecInterFaces即可
                                                               //vecRefIntFHash中存的面必然和当前面有联系
                                                               //

                nCommon = 0;
                for (j = 0;	j < ndIFaceSize; j++) {                 //对faceIt 也就是一个点周围一圈的面进行遍历找出重合面
                    nCommon = 0;
                    for (k = 0; k < 3; k++)
                        if (vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx1 ||
                                vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx2 ||
                                vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx3)
                            nCommon++;
                    if (nCommon >= 3)
                    {
          //              cout<<"find it !!!!!!!!!"<<endl;
                        break;
                    }
                }

                if (nCommon < 3) {//未找到重合面，循环结束退出，说明该面尚未放入容器，该面尚未登记过lftCell。则进行登记操作
                    if (nFaceSize >= nAllocFaceSize)                      //若空间不足则扩大空间
                    {
                        nAllocFaceSize += MAX_VALUE(100, (int)(nAllocFaceSize * 0.1));
                        vecInterFaces_Temp = (InterFace*)realloc(vecInterFaces, sizeof(InterFace)*nAllocFaceSize);
                        if (!vecInterFaces_Temp)
                        {
                            errCode = -1;
                            goto FAIL;
                        }
                        vecInterFaces = vecInterFaces_Temp;
                    }

                    faceAd.conn[0] = facNdIdx1;          //面的conn点
                    faceAd.conn[1] = facNdIdx2;
                    faceAd.conn[2] = facNdIdx3;
                    faceAd.lftCell = cellIdx;           //当前体中的面毫无疑问有一个邻接体是当前体
                    faceAd.rgtCell = -1;

                    faceAdIdx = nFaceSize;
                    faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];     //上一个faceAdIdx被放入此处，第一个faceAD的hasNxt为-1
                    vecRefIntFHash[minFacNdIdx] = faceAdIdx;          //faceAdIdx被放入此处，给下一个faceAd作为next。
                    vecInterFaces[nFaceSize++] = faceAd;              //faceAd顺序放入容器
                }


                                                                               //找到重合面由break退出
                else {

                    if(!(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                         vecInterFaces[ndIFaces[j]].rgtCell < 0))
                    {
                        cout<<ndIFaces[j]<<endl;
                    }

                    assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                            vecInterFaces[ndIFaces[j]].rgtCell < 0);


                    vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;              //既然有重合面，说明该面已经登记过，说明该面已有lftCell
                                                                               //那么当前体就是该面的rgtCell。

                    lftCell = vecInterFaces[ndIFaces[j]].lftCell;
                    rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;
                    pLft = &pPrism[lftCell];
                    pRgt = &pPrism[rgtCell];                             //

                 for(int k=0;k<6;++k){

                       if(pLft->neighbors[k]==-1)
                    {
                        pLft->neighbors[k] = rgtCell;
                        pLft->neighborsmark[k]=double((facevec[0]+facevec[1]*2+facevec[2]*3))/10;
                        break;
                    }
                 }

                    pRgt->neighbors[i] = lftCell;
                    pRgt->neighborsmark[i]=double((facevec[0]+facevec[1]*2+facevec[2]*3))/10;

                                                                //在两个相邻的体的neighbors中互相储存对方。（下标为什么为i和k？）
                                                  //i在此循环中代表该体的当前面，k则代表当前面的对点编号
                }                                                         //由cf[4][3]可以看出，此时k恒等于i
                                                                        //此时对于rgtCell来讲，自然是能保证neighbors[0]不被覆盖
                                                                        //lftCell如何保证？

            }



    }

    }
    else if(cellIdx<nPrism+nHexes){
         pHexes_t = &pHexes[cellIdx];
         HEX  *pLft=nullptr, *pRgt=nullptr;
        for (i = 0; i <= 5; i++) {
                    facNdIdx1 = pHexes_t->vertices[hex[i][0]];
                    facNdIdx2 = pHexes_t->vertices[hex[i][1]];                //BKG_MESH_DIM == 3  new BKG_MESH_DIM 应该为5
                    facNdIdx3 = pHexes_t->vertices[hex[i][2]];
                    facNdIdx4 = pHexes_t->vertices[hex[i][3]];

        int facevec[4]={facNdIdx1,facNdIdx2,facNdIdx3,facNdIdx4};
               sort(facevec,facevec+4);

                    int min=nNodes+10;
                    for(int index=0;index<4;++index) {
                         if(min>pPrism_t->vertices[hex[i][index]])
                         min=pPrism_t->vertices[hex[i][index]];
                    }
                    minFacNdIdx = min;
         //      cout<<facNdIdx1<<" "<<facNdIdx2<<" "<<facNdIdx3<<" "<<facNdIdx4<<" "<<min<<endl;
              //      cout<<min<<endl;
                           //找出最小序号点
                    //		nodeInterFace(minFacNdIdx, ndIFaces, &ndIFaceSize, vecInterFaces, vecRefIntFHash);


                    j = 0;
                    faceIt = vecRefIntFHash[minFacNdIdx];      //该点对应的neighbors为faceIt
                    while (faceIt >= 0)                        //若该点(最小点)已经有对应的Neighbors（面的neighbors储存在对角点）
                    {
                        if (j >= nAllocNodeFaceSize)           //j是什么？     ndIFaces空间不足 重新分配
                        {
                            nAllocNodeFaceSize += MAX_VALUE(100, (int)(nAllocNodeFaceSize * 0.1));
                            ndIFaces_Temp = (int*)realloc(ndIFaces, sizeof(int)*nAllocNodeFaceSize);
                            if (!ndIFaces_Temp)
                            {
                                errCode = -1;
                                goto FAIL;
                            }
                            ndIFaces = ndIFaces_Temp;
                        }
                        ndIFaces[j++] = faceIt;                //ndIFaces空间足够，则按顺序将
                                                                                                  //class InterFace
                                                                                                  //{
                                                                                                  //public:
                                                                                                  //    int conn[BKG_MESH_DIM];  //面的点ID
                                                                                                  //    int lftCell, rgtCell;    //左右邻接体
                                                                                                  //    int hashNxt;             //
                                                                                                  //};

                        faceIt= vecInterFaces[faceIt].hashNxt;  //neighbors的hasnext的意义？
                    }
                    ndIFaceSize = j;
              //      cout<<"j   "<<j<<endl;
                                                                   //j以及vecRefIntFHash的作用应该是为了加速查找重合面，否则遍历vecInterFaces即可
                                                                   //vecRefIntFHash中存的面必然和当前面有联系
                                                                   //

                    nCommon = 0;
                    for (j = 0;	j < ndIFaceSize; j++) {                 //对faceIt 也就是一个点周围一圈的面进行遍历找出重合面
                        nCommon = 0;
                        for (k = 0; k < 4; k++)
                            if (vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx1 ||
                                    vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx2 ||
                                    vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx3||
                                vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx4)
                                nCommon++;
                        if (nCommon >= 4)
                        {
          //                  cout<<"find it !!!!!!!!!"<<endl;
                            break;
                        }
                    }
               //     cout<<"nCommon"<<nCommon<<endl;

                    if (nCommon < 4) {//未找到重合面，循环结束退出，说明该面尚未放入容器，该面尚未登记过lftCell。则进行登记操作
                        /* 没锟斤拷锟揭碉拷 */
                        if (nFaceSize >= nAllocFaceSize)                      //若空间不足则扩大空间
                        {
                            nAllocFaceSize += MAX_VALUE(100, (int)(nAllocFaceSize * 0.1));
                            vecInterFaces_Temp = (InterFace*)realloc(vecInterFaces, sizeof(InterFace)*nAllocFaceSize);
                            if (!vecInterFaces_Temp)
                            {
                                errCode = -1;
                                goto FAIL;
                            }
                            vecInterFaces = vecInterFaces_Temp;
                        }

                        faceAd.conn[0] = facNdIdx1;          //面的conn点
                        faceAd.conn[1] = facNdIdx2;
                        faceAd.conn[2] = facNdIdx3;
                        faceAd.conn[3] = facNdIdx4;
                        faceAd.lftCell = cellIdx;           //当前体中的面毫无疑问有一个邻接体是当前体
                        faceAd.rgtCell = -1;

                        faceAdIdx = nFaceSize;
                        faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];     //上一个faceAdIdx被放入此处，第一个faceAD的hasNxt为-1
                        vecRefIntFHash[minFacNdIdx] = faceAdIdx;          //faceAdIdx被放入此处，给下一个faceAd作为next。
                        vecInterFaces[nFaceSize++] = faceAd;              //faceAd顺序放入容器
                    //    cout<<nFaceSize<<endl;
                    }


                                                                                   //找到重合面由break退出
                    else {

                        if(!(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                             vecInterFaces[ndIFaces[j]].rgtCell < 0))
                            cout<<cellIdx<<endl;

                        assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                                vecInterFaces[ndIFaces[j]].rgtCell < 0);


                        vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;              //既然有重合面，说明该面已经登记过，说明该面已有lftCell
                                                                                   //那么当前体就是该面的rgtCell。

                        lftCell = vecInterFaces[ndIFaces[j]].lftCell;
                        rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;
                        pLft = &pHexes[lftCell];
                        pRgt = &pHexes[rgtCell];                             //

                     for(int k=0;k<6;++k){
          //               cout<<"start"<<endl;
                           if(pLft->neighbors[k]==-1)
                        {
                            pLft->neighbors[k] = rgtCell;
                            pLft->neighborsmark[k]=double((facevec[0]+facevec[1]*2+facevec[2]*3+facevec[3]*4))/10;
                         //   cout<<pLft->neighborsmark[k]<<endl;
                            break;
                        }
                     }

                        pRgt->neighbors[i] = lftCell;
                        pRgt->neighborsmark[i]=double((facevec[0]+facevec[1]*2+facevec[2]*3+facevec[3]*4))/10;

                                                                    //在两个相邻的体的neighbors中互相储存对方。（下标为什么为i和k？）
                                                      //i在此循环中代表该体的当前面，k则代表当前面的对点编号
                    }                                                         //由cf[4][3]可以看出，此时k恒等于i
                                                                            //此时对于rgtCell来讲，自然是能保证neighbors[0]不被覆盖
                                                                            //lftCell如何保证？
                }
            }

    else if(cellIdx<nPrism+nHexes+nTetras){
       pTetras_t = &pTetras[cellIdx];
       TETRAS  *pLft=nullptr, *pRgt=nullptr;
       for (i = 0; i <= 3; i++) {
                   facNdIdx1 = pTetras_t->vertices[tet[i][0]];
                   facNdIdx2 = pTetras_t->vertices[tet[i][1]];                //BKG_MESH_DIM == 3  new BKG_MESH_DIM 应该为5
                   facNdIdx3 = pTetras_t->vertices[tet[i][2]];

       int facevec[3]={facNdIdx1,facNdIdx2,facNdIdx3};
              sort(facevec,facevec+3);

                   int min=nNodes+10;
                   for(int index=0;index<3;++index) {
                        if(min>pTetras_t->vertices[tet[i][index]])
                        min=pTetras_t->vertices[tet[i][index]];
                   }
                   minFacNdIdx = min;
        //      cout<<facNdIdx1<<" "<<facNdIdx2<<" "<<facNdIdx3<<" "<<facNdIdx4<<" "<<min<<endl;
             //      cout<<min<<endl;
                          //找出最小序号点
                   //		nodeInterFace(minFacNdIdx, ndIFaces, &ndIFaceSize, vecInterFaces, vecRefIntFHash);


                   j = 0;
                   faceIt = vecRefIntFHash[minFacNdIdx];      //该点对应的neighbors为faceIt
                   while (faceIt >= 0)                        //若该点(最小点)已经有对应的Neighbors（面的neighbors储存在对角点）
                   {
                       if (j >= nAllocNodeFaceSize)           //j是什么？     ndIFaces空间不足 重新分配
                       {
                           nAllocNodeFaceSize += MAX_VALUE(100, (int)(nAllocNodeFaceSize * 0.1));
                           ndIFaces_Temp = (int*)realloc(ndIFaces, sizeof(int)*nAllocNodeFaceSize);
                           if (!ndIFaces_Temp)
                           {
                               errCode = -1;
                               goto FAIL;
                           }
                           ndIFaces = ndIFaces_Temp;
                       }
                       ndIFaces[j++] = faceIt;                //ndIFaces空间足够，则按顺序将
                                                                                                 //class InterFace
                                                                                                 //{
                                                                                                 //public:
                                                                                                 //    int conn[BKG_MESH_DIM];  //面的点ID
                                                                                                 //    int lftCell, rgtCell;    //左右邻接体
                                                                                                 //    int hashNxt;             //
                                                                                                 //};

                       faceIt= vecInterFaces[faceIt].hashNxt;  //neighbors的hasnext的意义？
                   }
                   ndIFaceSize = j;
             //      cout<<"j   "<<j<<endl;
                                                                  //j以及vecRefIntFHash的作用应该是为了加速查找重合面，否则遍历vecInterFaces即可
                                                                  //vecRefIntFHash中存的面必然和当前面有联系
                                                                  //

                   nCommon = 0;
                   for (j = 0;	j < ndIFaceSize; j++) {                 //对faceIt 也就是一个点周围一圈的面进行遍历找出重合面
                       nCommon = 0;
                       for (k = 0; k < 3; k++)
                           if (vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx1 ||
                                   vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx2 ||
                                   vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx3)
                               nCommon++;
                       if (nCommon >= 3)
                       {
         //                  cout<<"find it !!!!!!!!!"<<endl;
                           break;
                       }
                   }
              //     cout<<"nCommon"<<nCommon<<endl;

                   if (nCommon < 3) {//未找到重合面，循环结束退出，说明该面尚未放入容器，该面尚未登记过lftCell。则进行登记操作
                       /* 没锟斤拷锟揭碉拷 */
                       if (nFaceSize >= nAllocFaceSize)                      //若空间不足则扩大空间
                       {
                           nAllocFaceSize += MAX_VALUE(100, (int)(nAllocFaceSize * 0.1));
                           vecInterFaces_Temp = (InterFace*)realloc(vecInterFaces, sizeof(InterFace)*nAllocFaceSize);
                           if (!vecInterFaces_Temp)
                           {
                               errCode = -1;
                               goto FAIL;
                           }
                           vecInterFaces = vecInterFaces_Temp;
                       }

                       faceAd.conn[0] = facNdIdx1;          //面的conn点
                       faceAd.conn[1] = facNdIdx2;
                       faceAd.conn[2] = facNdIdx3;
                       faceAd.lftCell = cellIdx;           //当前体中的面毫无疑问有一个邻接体是当前体
                       faceAd.rgtCell = -1;

                       faceAdIdx = nFaceSize;
                       faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];     //上一个faceAdIdx被放入此处，第一个faceAD的hasNxt为-1
                       vecRefIntFHash[minFacNdIdx] = faceAdIdx;          //faceAdIdx被放入此处，给下一个faceAd作为next。
                       vecInterFaces[nFaceSize++] = faceAd;              //faceAd顺序放入容器
                   //    cout<<nFaceSize<<endl;
                   }


                                                                                  //找到重合面由break退出
                   else {

                       if(!(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                            vecInterFaces[ndIFaces[j]].rgtCell < 0))
                           cout<<cellIdx<<endl;

                       assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                               vecInterFaces[ndIFaces[j]].rgtCell < 0);


                       vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;              //既然有重合面，说明该面已经登记过，说明该面已有lftCell
                                                                                  //那么当前体就是该面的rgtCell。

                       lftCell = vecInterFaces[ndIFaces[j]].lftCell;
                       rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;
                       pLft = &pTetras[lftCell];
                       pRgt = &pTetras[rgtCell];                             //

                    for(int k=0;k<4;++k){
         //               cout<<"start"<<endl;
                          if(pLft->neighbors[k]==-1)
                       {
                           pLft->neighbors[k] = rgtCell;
                           pLft->neighborsmark[k]=double((facevec[0]+facevec[1]*2+facevec[2]*3))/10;
                        //   cout<<pLft->neighborsmark[k]<<endl;
                           break;
                       }
                    }

                       pRgt->neighbors[i] = lftCell;
                       pRgt->neighborsmark[i]=double((facevec[0]+facevec[1]*2+facevec[2]*3))/10;

                                                                   //在两个相邻的体的neighbors中互相储存对方。（下标为什么为i和k？）
                                                     //i在此循环中代表该体的当前面，k则代表当前面的对点编号
                   }                                                         //由cf[4][3]可以看出，此时k恒等于i
                                                                           //此时对于rgtCell来讲，自然是能保证neighbors[0]不被覆盖
                                                                           //lftCell如何保证？
               }
    }
    }

    goto END;
FAIL:
END:
    if (vecRefIntFHash)
        free(vecRefIntFHash);
    if (vecInterFaces)
        free(vecInterFaces);
    if (ndIFaces)
        free(ndIFaces);

    return errCode; /* S_OK */
}


int setupCellNeig(int nNodes, int nElems, TETRAS *pBKGElems)
{

    int nAllocFaceSize = nNodes * 10, nAllocNodeFaceSize = 20384;
    int nFaceSize = 0;
    int *vecRefIntFHash = NULL;
    InterFace * vecInterFaces = NULL, *vecInterFaces_Temp = NULL;
    int *ndIFaces = NULL, *ndIFaces_Temp = NULL;
    int errCode = 0;
    int cellIdx, faceAdIdx, lftCell, rgtCell;
    int facNdIdx1, facNdIdx2, facNdIdx3, minFacNdIdx, faceIt;
    int ndIFaceSize = 0;
    int i, j, k, nCommon;
    InterFace faceAd;
    static int cf[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    int ndSize, clSize;
    //BKGElem *pElem = NULL, *pLft = NULL, *pRgt = NULL;

    TETRAS *pElem = nullptr, *pLft=nullptr, *pRgt=nullptr;

    vecRefIntFHash = (int *) malloc(sizeof(int)*nNodes);
    vecInterFaces = (InterFace *) malloc(sizeof(InterFace)*nAllocFaceSize);
    ndIFaces = (int *) malloc(sizeof(int)*nAllocNodeFaceSize);
    if (!vecRefIntFHash || !vecInterFaces || !ndIFaces)
    {
        errCode = -1;
        goto FAIL;
    }

    ndSize = nNodes;
    clSize = nElems;

    for (i = 0; i < ndSize; i++)
        vecRefIntFHash[i] = -1;
    for (cellIdx = 0; cellIdx < clSize; cellIdx++)
    {
        pElem = &pBKGElems[cellIdx];
        for (i = 0; i <= BKG_MESH_DIM; i++)
            pElem->neighbors[i] = -1;
    }

    for (cellIdx = 0; cellIdx < clSize; cellIdx++)
    {
        if (cellIdx % 1000000 == 0 || cellIdx == clSize - 1)
            printf("%%%.2f.\n", 100.0*(cellIdx)/clSize);

        pElem = &pBKGElems[cellIdx];

        for (i = 0; i <= BKG_MESH_DIM; i++) {
            facNdIdx1 = pElem->vertices[cf[i][0]];
            facNdIdx2 = pElem->vertices[cf[i][1]];
            facNdIdx3 = pElem->vertices[cf[i][2]];

            minFacNdIdx = facNdIdx1 < facNdIdx2 ? facNdIdx1 : facNdIdx2;
            minFacNdIdx = minFacNdIdx < facNdIdx3 ? minFacNdIdx : facNdIdx3;
            //		nodeInterFace(minFacNdIdx, ndIFaces, &ndIFaceSize, vecInterFaces, vecRefIntFHash);

            j = 0;
            faceIt = vecRefIntFHash[minFacNdIdx];
            while (faceIt >= 0)
            {
                if (j >= nAllocNodeFaceSize)
                {
                    nAllocNodeFaceSize += MAX_VALUE(100, (int)(nAllocNodeFaceSize * 0.1));
                    ndIFaces_Temp = (int*)realloc(ndIFaces, sizeof(int)*nAllocNodeFaceSize);
                    if (!ndIFaces_Temp)
                    {
                        errCode = -1;
                        goto FAIL;
                    }
                    ndIFaces = ndIFaces_Temp;
                }
                ndIFaces[j++] = faceIt;
                faceIt= vecInterFaces[faceIt].hashNxt;
            }
            ndIFaceSize = j;


            nCommon = 0;
            for (j = 0;	j < ndIFaceSize; j++) {
                nCommon = 0;
                for (k = 0; k < 3; k++)
                    if (vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx1 ||
                            vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx2 ||
                            vecInterFaces[ndIFaces[j]].conn[k] == facNdIdx3)
                        nCommon++;
                if (nCommon >= 3)
                    break;
            }

            if (nCommon < 3) {
                /* 没锟斤拷锟揭碉拷 */
                if (nFaceSize >= nAllocFaceSize)
                {
                    nAllocFaceSize += MAX_VALUE(100, (int)(nAllocFaceSize * 0.1));
                    vecInterFaces_Temp = (InterFace*)realloc(vecInterFaces, sizeof(InterFace)*nAllocFaceSize);
                    if (!vecInterFaces_Temp)
                    {
                        errCode = -1;
                        goto FAIL;
                    }
                    vecInterFaces = vecInterFaces_Temp;
                }

                faceAd.conn[0] = facNdIdx1;
                faceAd.conn[1] = facNdIdx2;
                faceAd.conn[2] = facNdIdx3;
                faceAd.lftCell = cellIdx;
                faceAd.rgtCell = -1;

                faceAdIdx = nFaceSize;
                faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];
                vecRefIntFHash[minFacNdIdx] = faceAdIdx;
                vecInterFaces[nFaceSize++] = faceAd;
            }
            else {

                if(!(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                     vecInterFaces[ndIFaces[j]].rgtCell < 0))
                    cout<<cellIdx<<endl;

                assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
                        vecInterFaces[ndIFaces[j]].rgtCell < 0);


                vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;

                lftCell = vecInterFaces[ndIFaces[j]].lftCell;
                rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;
                pLft = &pBKGElems[lftCell];
                pRgt = &pBKGElems[rgtCell];
                for (k = 0; k < 4; k++)
                    if (pLft->vertices[k] != facNdIdx1 &&
                            pLft->vertices[k] != facNdIdx2 &&
                            pLft->vertices[k] != facNdIdx3)
                        break;
                assert(k < 4);
                pRgt->neighbors[i] = lftCell;
      //          cout<<i<<"  "<<k <<endl;
                pLft->neighbors[k] = rgtCell;
            }
        }
    }

    goto END;
FAIL:
END:
    if (vecRefIntFHash)
        free(vecRefIntFHash);
    if (vecInterFaces)
        free(vecInterFaces);
    if (ndIFaces)
        free(ndIFaces);

    return errCode; /* S_OK */
}


int readVTKPLSFile(char *fname,HYBRID_MESH& file)
{

    FILE *fp = NULL;
    Node *pBKGNodes = nullptr;
    TETRAS *pBKGElems = nullptr;
    Node *pNode=nullptr;
    TETRAS *pElem=nullptr;

    // ELEMENT<4,4> ele;

    //ELEMENT<4,4>*p = new ELEMENT<4,4>();

    int nNodeNum = 0, nElemNum = 0;
    char chLine[MAX_FILE_LINE];
    int i, j, conn, nValueNum;
    int errCode = 0;


    char vtkfile[256];
    sprintf(vtkfile, "%s%s",fname,".vtk");

    if (!vtkfile)
    {
        printf("No VTK file name is given.\n");
        return 2;
    }
    fp = fopen(vtkfile, "r");
    if (!(fp = fopen(vtkfile, "r")))
    {
        printf("Cannot open file %s.\n", vtkfile);
        return 2;
    }

    /* points */
    do
    {
        readValidLine(fp, chLine, MAX_FILE_LINE);
        eraseWhiteSpace(chLine);
        if (strstr(chLine, "POINTS") != NULL)
            break;
    }
    while (!feof(fp));

    if (feof(fp))
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    sscanf(chLine, "POINTS %d", &nNodeNum);

    file.NumNodes=nNodeNum;

    if (nNodeNum <= 0)
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    if (!(file.nodes = new Node[nNodeNum]))
    {
        errCode = -1; /* not enough memory */
        goto FAIL;
    }

    pBKGNodes=file.nodes;

    for (i = 0; i < nNodeNum; i++)
    {
        pNode = &pBKGNodes[i];
        pNode->index=i;
        fscanf(fp, "%lf", &pNode->coord.x);
        fscanf(fp, "%lf", &pNode->coord.y);
        fscanf(fp, "%lf", &pNode->coord.z);
    }

    /* elements */
    do
    {
        readValidLine(fp, chLine, MAX_FILE_LINE);
        eraseWhiteSpace(chLine);
        if (strstr(chLine, "CELLS") != NULL)
            break;
    }
    while (!feof(fp));

    if (feof(fp))
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    sscanf(chLine, "CELLS %d", &nElemNum);

    file.NumTetras=nElemNum;

    if (nNodeNum <= 0)
    {
        errCode = 3; /* Invalid file */
        goto FAIL;
    }

    if (!(file.pTetras = new TETRAS[nElemNum]))
    {
        errCode = -1; /* not enough memory */
        goto FAIL;
    }

    pBKGElems=file.pTetras;

    for (i = 0; i < nElemNum; i++)
    {
        fscanf(fp, "%d", &nValueNum);
        if (nValueNum != BKG_MESH_DIM + 1)
        {
            errCode = 3; /* Invalid file */
            goto FAIL;
        }

        pElem = &pBKGElems[i];
        pElem->index=i;
        for (j = 0; j <= BKG_MESH_DIM; j++)
        {
            fscanf(fp, "%d", &conn);
            //if (conn <= 0 || conn > nNodeNum)
            if (conn < 0 || conn >= nNodeNum)
            {
                errCode = 3; /* Invalid file */
                goto FAIL;
            }
            pElem->vertices[j] = conn;/*conn - 1;*/
            pElem->neighbors[j] = -1;
        }
    }

    /* reconstruct neighbouring indices */
    errCode = setupCellNeig(nNodeNum, nElemNum, pBKGElems);
    if (errCode != 0 && errCode != 1)
        goto FAIL;


    goto END;
FAIL:

END:
    fclose(fp);


//read pls file
    char plsfile[256];
    sprintf(plsfile, "%s%s",fname,".pls");
    fp = fopen(plsfile, "r");
    if (!plsfile)
    {
        printf("Error: cann't open file %s!\n", plsfile);
        return 0;
    }

    int NumTris=0; int NumPoints=0;  int tmp=-1;

    fgets(chLine, MAX_FILE_LINE, fp);
    sscanf(chLine, "%d %d", &NumTris, &NumPoints);
    //fscanf(fp, "%d %d %d %d %d %d\n", &NumTris, &NumPoints, &tmp, &tmp, &tmp, &tmp);

   // *nbele = nel;
   // (*bele) = new int[nel * 4];
	file.NumTris=NumTris;
	file.pTris=new TRI[NumTris];

    float ftmp=0;
    for (int i=0; i<NumPoints; i++)
    {
        fscanf(fp, "%d %f %f %f\n", &tmp, &ftmp, &ftmp, &ftmp);
    }

    for (int i=0; i<NumTris; i++)
    {
        int id1, id2, id3, fidx;
        fscanf(fp, "%d %d %d %d %d\n", &tmp, &id1, &id2, &id3, &fidx);
        file.pTris[i].vertices[0]=id1-1;
        file.pTris[i].vertices[1]=id2-1;
        file.pTris[i].vertices[2]=id3-1;
        file.pTris[i].iSurf=fidx-1;
        file.pTris[i].index=i;
    }

    fclose(fp);
    fp = nullptr;

    //findiCellFast(file);


    return errCode;
}


int findiCellFast(HYBRID_MESH &file)
{
    int neigbor=-2;
    int count=0;


    set<int>cellIDs;

    map<string,int> trimap;

    int v[3];
    for(int i=0;i<file.NumTris;i++)
    {
        v[0]=file.pTris[i].vertices[0];
        v[1]=file.pTris[i].vertices[1];
        v[2]=file.pTris[i].vertices[2];
        sort(v,v+3);
        string temp=IntToString(v[0])+"_"+IntToString(v[1])+"_"+IntToString(v[2]);
        trimap[temp]=i;
    }


    map<string,int>::iterator mapIter;
    for(int i=0;i<file.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            neigbor=file.pTetras[i].neighbors[j];

            if(neigbor==-1)
            {
               // count++;

                cellIDs.insert(i);
                count=0;
                for(int k=0;k<4;k++)
                {
                    if(k!=j)
                    {
                        v[count]=file.pTetras[i].vertices[k];
                        count++;

                    }
                }
                sort(v,v+3);
                string temp=IntToString(v[0])+"_"+IntToString(v[1])+"_"+IntToString(v[2]);
                mapIter=trimap.find(temp);
                if(mapIter==trimap.end())
                {
                    cout<<"Error in finding iCell!"<<endl;
                }
                else
                {
                    int triID=mapIter->second;
                    file.pTris[triID].iCell=i;
                 //   cout<<i<<endl;
                }
                assert(count==3);
            }
        }
    }

    cout<<"Finding iCell finished!"<<endl;
    return 1;
}

int findiCell(HYBRID_MESH&file)
{
    int neigbor=-2;
    int count=0;


    set<int>cellIDs;


    for(int i=0;i<file.NumTetras;i++)
    {
        for(int j=0;j<4;j++)
        {
            neigbor=file.pTetras[i].neighbors[j];
            if(neigbor==-1)
            {
                count++;
                cellIDs.insert(i);
            }
        }
    }



    set<int>::iterator iter;

    //if(iter==cellIDs.end())
    //{
    //    cout<<"Error "<<endl;
    //}

    cout<<"NumTri :" <<file.NumTris<<endl;
   // assert(count==file.NumTri);


    count=0;
    int tempID=-1;
    int triVer[3];
    int tetVer[4];
    int triCompare[3];
    for(int i=0;i<file.NumTris;i++)
    {
        triVer[0]=file.pTris[i].vertices[0];
        triVer[1]=file.pTris[i].vertices[1];
        triVer[2]=file.pTris[i].vertices[2];



        bool same=false;
        for(iter=cellIDs.begin();iter!=cellIDs.end();iter++)
        {
            tempID=*iter;

           triCompare[0]=file.pTetras[tempID].vertices[0];
           triCompare[1]=file.pTetras[tempID].vertices[1];
           triCompare[2]=file.pTetras[tempID].vertices[2];
           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }

           triCompare[0]=file.pTetras[tempID].vertices[0];
           triCompare[1]=file.pTetras[tempID].vertices[1];
           triCompare[2]=file.pTetras[tempID].vertices[3];

           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }
           triCompare[0]=file.pTetras[tempID].vertices[0];
           triCompare[1]=file.pTetras[tempID].vertices[2];
           triCompare[2]=file.pTetras[tempID].vertices[3];

           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }
           triCompare[0]=file.pTetras[tempID].vertices[1];
           triCompare[1]=file.pTetras[tempID].vertices[2];
           triCompare[2]=file.pTetras[tempID].vertices[3];

           same=isSameTriangle(triVer,triCompare);
           if(same==true)
           {
               file.pTris[i].iCell=tempID;
               break;
           }

        }
        if(same==false)
        {
            cout<<"Error in finding the owner tetras! Tri v is "<<triVer[0]<<" "
               <<triVer[1]<<" "<<triVer[2]<<endl;
            exit(0);
        }
    }
    cout<<"Finding iCell finished!"<<endl;

    return 1;
}

bool tetrasContainTriangle(int *tri, int *tetras)
{
    /*
    sort(tri,tri+3);
    sort(tetras,tetras+4);

    if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
    else if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
    else if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
    else if(tri[0]==tetras[0]&&tri[0]==tetras[0]&&tri[0]==tetras[0])
        return true;
	*/
    return true;
    
}
bool isSameTriangle(int *tri,int *triCom)
{
    sort(tri,tri+3);
    sort(triCom,triCom+3);
    if(tri[0]==triCom[0]&&tri[1]==triCom[1]&&tri[2]==triCom[2])
    {
        return true;
    }
    else
        return false;
}
int writeTriangleVTKFile(char *filename,HYBRID_MESH&file)
{
    int i, j = 0, k = 0, iFac = 0, ii1 = 0, ii2 = 0, ii3 = 0;
    int faceNumber = 0;
   // MYPOINT pnt;
  //  DesFacet *pDF = NULL;
   // bool *visited = NULL;

	int elemSize = file.NumTris, nodeSize = file.NumNodes;

    FILE *fout = fopen(filename,"w");
    fprintf(fout,"# vtk DataFile Version 2.0\n");
    fprintf(fout,"background mesh\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout,"POINTS %d double\n", nodeSize);


    for (i = 0; i < nodeSize; i++)
    {
        int index = i + 1;
        double x = file.nodes[i].coord.x;
        double y = file.nodes[i].coord.y;
        double z = file.nodes[i].coord.z;
        fprintf(fout,"%lf %lf %lf\n",x, y, z);

    }
    fprintf(fout,"CELLS %d %d \n", elemSize, (elemSize) * 4);

    for (i = 0; i < elemSize; i++)
    {

        int index = i + 1;
        int form0 = file.pTris[i].vertices[0];///////////+1
        int form1 = file.pTris[i].vertices[1];
        int form2 = file.pTris[i].vertices[2];
        //int form3 = file.pTetras[i].vertices[3];

        fprintf(fout,"3 %d %d %d\n",form0,form1,form2);

    }
    fprintf(fout,"CELL_TYPES %d \n", elemSize);
    for (i = 0; i < elemSize; i++)
    {
        fprintf(fout,"%d\n", 5);
    }

	fprintf(fout,"\nPOINT_DATA %d \n", nodeSize);

#if 0
    fprintf(fout,"SCALARS procs_begin int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < nodeSize; i++)
    {
        //double size = m_pNodes[i].space/m_scale;
       // fprintf(fout,"%d\n",file.nodes[i].index);
        fprintf(fout,"%d\n",*file.nodes[i].procs.begin());
    }
#endif


    fprintf(fout,"SCALARS nodeID int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < nodeSize; i++)
    {
        //double size = m_pNodes[i].space/m_scale;
        fprintf(fout,"%d\n",file.nodes[i].index);
        // fprintf(fout,"%d\n",file.nodes[i].procs.size());
    }

	fprintf(fout,"\nCELL_DATA %d \n", elemSize);
    fprintf(fout,"SCALARS iOppoProc int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < elemSize; i++)
    {

       // double quality =  m_pElems[i].minAngle;
        fprintf(fout,"%d\n",file.pTris[i].iOppoProc);
      //  fprintf(fout,"%d\n",file.pTris[i].iCell);
    }


    fprintf(fout,"SCALARS iCell int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < elemSize; i++)
    {
        fprintf(fout,"%d\n",file.pTris[i].iCell);
    }

	fprintf(fout,"\nSCALARS facetID int\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i = 0; i < elemSize; i++)
	{
	   // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%lld\n",file.pTris[i].index);
	}

	fprintf(fout,"\nSCALARS iSurf int\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i = 0; i < elemSize; i++)
	{
	   // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%d\n",file.pTris[i].iSurf);
	}

	fclose(fout);
    return 1;
}

int writeVTKFile(char *filename,HYBRID_MESH&file)
{

    int i, j = 0, k = 0, iFac = 0, ii1 = 0, ii2 = 0, ii3 = 0;
    int faceNumber = 0;
   // MYPOINT pnt;
  //  DesFacet *pDF = NULL;
   // bool *visited = NULL;

    int elemSize = file.NumTetras, nodeSize = file.NumNodes;

    FILE *fout = fopen(filename,"w");
    fprintf(fout,"# vtk DataFile Version 2.0\n");
    fprintf(fout,"background mesh\n");
    fprintf(fout,"ASCII\n");
    fprintf(fout,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout,"POINTS %d double\n", nodeSize);


    for (i = 0; i < nodeSize; i++)
    {
        int index = i + 1;
        double x = file.nodes[i].coord.x;
        double y = file.nodes[i].coord.y;
        double z = file.nodes[i].coord.z;
        fprintf(fout,"%lf %lf %lf\n",x, y, z);

    }
    fprintf(fout,"CELLS %d %d \n", elemSize, (elemSize) * 5);

    for (i = 0; i < elemSize; i++)
    {

        int index = i + 1;
        int form0 = file.pTetras[i].vertices[0];///////////+1
        int form1 = file.pTetras[i].vertices[1];
        int form2 = file.pTetras[i].vertices[2];
        int form3 = file.pTetras[i].vertices[3];

        fprintf(fout,"4 %d %d %d %d\n",form0,form1,form2,form3);

    }
    fprintf(fout,"CELL_TYPES %d \n", elemSize);
    for (i = 0; i < elemSize; i++)
    {
        fprintf(fout,"%d\n", 10);
    }

	fprintf(fout,"\nPOINT_DATA %d \n", nodeSize);
	fprintf(fout,"SCALARS nodeID int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < nodeSize; i++)
    {
        //double size = m_pNodes[i].space/m_scale;
        fprintf(fout,"%d\n",file.nodes[i].index);
    }

	fprintf(fout,"\nCELL_DATA %d \n", elemSize);
    fprintf(fout,"SCALARS part int\n");
    fprintf(fout,"LOOKUP_TABLE default\n");
    for (i = 0; i < elemSize; i++)
    {
       // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%d\n",file.pTetras[i].partMarker);
	}


	fprintf(fout,"\nSCALARS cellID int\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i = 0; i < elemSize; i++)
	{
	   // double quality =  m_pElems[i].minAngle;
		fprintf(fout,"%d\n",file.pTetras[i].index);
	}


    fclose(fout);
    return 1;


    //return 1;
}
#if 1
int readCGNS_temp(char*filename,HYBRID_MESH&mesh,vector<string>&bcstring)
{
  //  float x[21*17*9],y[21*17*9],z[21*17*9];
    cgsize_t isize[3][1];//ielem[20*16*8][8];
    int index_file,index_base,index_zone;
    cgsize_t irmin,irmax,istart,iend;
    int nsections,index_sect,nbndry,iparent_flag;
    cgsize_t iparentdata;
    char zonename[33],sectionname[33];
    CGNS_ENUMT(ElementType_t) itype;

/* READ X, Y, Z GRID POINTS FROM CGNS FILE */
/* open CGNS file for read-only */
    int filetype=-1;
    if (!cg_is_cgns(filename,&filetype))
        cout<<"CGNS FILE: type="<<filetype<<endl;
    else
        cout<<"NOT CGNS FILE."<<endl;

    cg_set_file_type(filetype);
    if (cg_open(filename,CG_MODE_READ,&index_file))
        cg_error_exit();
    cout<<"Open CGNS file"<<endl;
/* we know there is only one base (real working code would check!) */
    int nbases=-1;
    cg_nbases(index_file,&nbases);
    if(nbases!=1)
    {
        cout<<"The bases is more than 1!---"<<nbases<<endl;
    }
    cout<<"Get base number"<<endl;
    index_base=1;
/* we know there is only one zone (real working code would check!) */
    index_zone=1;
    int nzones=-1;
    cg_nzones(index_file,index_base,&nzones);
    if(nzones!=1)
    {
        cout<<"The zones is more than 1!---"<<nzones<<endl;
    }
/* get zone size (and name - although not needed here) */
    cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);
    cout<<"read Zone"<<endl;
/* lower range index */
    irmin=1;
/* upper range index of vertices */
    irmax=isize[0][0];
    cout<<"irmax = "<< irmax <<endl;

    int nNodes=irmax;

    cgsize_t nprism;
    cgsize_t ntetras;
    cgsize_t *prism;
    cgsize_t *tetras;


/* read grid coordinates */
    cout<<"Starting to read grid coordinates"<<endl;
    int ncoord=-1;
    cg_ncoords(index_file,index_base,index_zone,&ncoord);

    float *x=new float[irmax];
    float *y=new float[irmax];
    float *z=new float[irmax];

    cout<<"Read grid coordinates"<<endl;
    cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,x);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,y);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,z);

    cout<<"The number of vertices is: "<<irmax<<" "<<ncoord<<endl;
/* find out how many sections */
    cg_nsections(index_file,index_base,index_zone,&nsections);
    printf("\nnumber of sections=%i\n",nsections);
/* read element connectivity */

    int nTriSum=0;

    //cgsize_t *triSection[3];
    cgsize_t **triSection=nullptr;
    cgsize_t **quadSection=nullptr;
    triSection=new cgsize_t*[nsections-2];
    quadSection=new cgsize_t*[nsections-1];
    int *nTriSection=new int[nsections-2];
    int *nQuadSection=new int[nsections-1];
    //vector<vector<cgsize_t>>triSection;
    int tri_count=0;
    int quad_count=0;

    for (index_sect=1; index_sect <= nsections; index_sect++)
    {
      cg_section_read(index_file,index_base,index_zone,index_sect,sectionname,
                      &itype,&istart,&iend,&nbndry,&iparent_flag);


      printf("\nReading section data...\n");
      printf("   section name=%s\n",sectionname);
      printf("   section type=%s\n",ElementTypeName[itype]);
      printf("   istart,iend=%i, %i\n",(int)istart,(int)iend);
      if (itype == CGNS_ENUMV(TETRA_4))
      {
        printf("   reading element data for this Tetras\n");
        ntetras=(iend-istart+1);
       tetras=new cgsize_t[ntetras*4];
        cg_elements_read(index_file,index_base,index_zone,index_sect,tetras, \
                         &iparentdata);
      }
      if (itype == CGNS_ENUMV(PENTA_6))
      {
         nprism=(iend-istart+1);
         prism=new cgsize_t[nprism*6];
        printf("   reading element data for this Prism\n");
        cg_elements_read(index_file,index_base,index_zone,index_sect,prism, \
                         &iparentdata);
      }

      if (itype == CGNS_ENUMV(TRI_3))
      {
          int nTri=iend-istart+1;
          nTriSection[tri_count]=nTri;
          nTriSum+=nTri;
          triSection[tri_count]=new cgsize_t[nTri*3];
          string strName=sectionname;
          bcstring.push_back(strName);
          printf("   reading element data for this TRIS_3\n");
          cg_elements_read(index_file,index_base,index_zone,index_sect,triSection[tri_count], \
                  &iparentdata);
          tri_count++;
      }
      if (itype == CGNS_ENUMV(QUAD_4))
      {
          int nQuad=iend-istart+1;
          nQuadSection[quad_count]=nQuad;
          nTriSum+=nQuad;
          quadSection[quad_count]=new cgsize_t[nQuad*4];
          string strName=sectionname;
          bcstring.push_back(strName);
          printf("   reading element data for this QUAD_4\n");
          cg_elements_read(index_file,index_base,index_zone,index_sect,triSection[quad_count], \
                  &iparentdata);
          quad_count++;
      }
    }
/* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully read unstructured grid from file grid_c.cgns\n");

    mesh.NumNodes=nNodes;
    mesh.nodes=new Node[nNodes];
    mesh.NumPrsm=nprism;
    mesh.pPrisms=new PRISM[nprism];
    mesh.NumTetras=ntetras;
    mesh.pTetras=new TETRAS[ntetras];
    mesh.NumTris=nTriSum;
    mesh.pTris=new TRI[nTriSum];

    for(int i=0;i<nNodes;i++)
    {
        mesh.nodes[i].index=i;
        mesh.nodes[i].coord.x=x[i];
        mesh.nodes[i].coord.y=y[i];
        mesh.nodes[i].coord.z=z[i];
    }

    for(int i=0;i<nprism;i++)
    {

        mesh.pTetras[i].index=i;
        mesh.pTetras[i].vertices[0]=prism[i*6+0]-1;
        mesh.pTetras[i].vertices[1]=prism[i*6+1]-1;
        mesh.pTetras[i].vertices[2]=prism[i*6+2]-1;
        mesh.pTetras[i].vertices[3]=prism[i*6+3]-1;
        mesh.pTetras[i].vertices[4]=prism[i*6+4]-1;
        mesh.pTetras[i].vertices[5]=prism[i*6+5]-1;
    }

    for(int i=0;i<ntetras;i++)
    {
        mesh.pTetras[i].index=i;
        mesh.pTetras[i].vertices[0]=tetras[i*4+0]-1;
        mesh.pTetras[i].vertices[1]=tetras[i*4+1]-1;
        mesh.pTetras[i].vertices[2]=tetras[i*4+2]-1;
        mesh.pTetras[i].vertices[3]=tetras[i*4+3]-1;
    }


    tri_count=0;
   // mesh.pTris[i].index=i;
    for(int j=0;j<nsections-2;j++)
    {
        for(int k=0;k<nTriSection[j];k++)
        {
            mesh.pTris[tri_count].index=tri_count;
            mesh.pTris[tri_count].vertices[0]=triSection[j][k*3+0]-1;
            mesh.pTris[tri_count].vertices[1]=triSection[j][k*3+1]-1;
            mesh.pTris[tri_count].vertices[2]=triSection[j][k*3+2]-1;
            mesh.pTris[tri_count].iSurf=j;
            tri_count++;
        }
    }
 cout<<"succesfull"<<endl;
    assert(count==nTriSum);

 cout<<"succesfull"<<endl;
    return 1;
}

int readCGNS(char*filename,HYBRID_MESH&mesh,vector<string>&bcstring)
{
  //  float x[21*17*9],y[21*17*9],z[21*17*9];
    cgsize_t isize[3][1];//ielem[20*16*8][8];
    int index_file,index_base,index_zone;
    cgsize_t irmin,irmax,istart,iend;
    int nsections,index_sect,nbndry,iparent_flag;
    cgsize_t iparentdata;
    char zonename[33],sectionname[33];
    CGNS_ENUMT(ElementType_t) itype;

/* READ X, Y, Z GRID POINTS FROM CGNS FILE */
/* open CGNS file for read-only */
    int filetype=-1;
    if (!cg_is_cgns(filename,&filetype))
        cout<<"CGNS FILE: type="<<filetype<<endl;
    else
        cout<<"NOT CGNS FILE."<<endl;

    cg_set_file_type(filetype);
    if (cg_open(filename,CG_MODE_READ,&index_file))
        cg_error_exit();
    cout<<"Open CGNS file"<<endl;
/* we know there is only one base (real working code would check!) */
    int nbases=-1;
    cg_nbases(index_file,&nbases);
    if(nbases!=1)
    {
        cout<<"The bases is more than 1!---"<<nbases<<endl;
    }
    cout<<"Get base number"<<endl;
    index_base=1;
/* we know there is only one zone (real working code would check!) */
    index_zone=1;
    int nzones=-1;
    cg_nzones(index_file,index_base,&nzones);
    if(nzones!=1)
    {
        cout<<"The zones is more than 1!---"<<nzones<<endl;
    }
/* get zone size (and name - although not needed here) */
    cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);
    cout<<"read Zone"<<endl;
/* lower range index */
    irmin=1;
/* upper range index of vertices */
    irmax=isize[0][0];
    cout<<"irmax = "<< irmax <<endl;

    int nNodes=irmax;

    cgsize_t nTetras=isize[1][0];
    cgsize_t *tetras=new cgsize_t[nTetras*4];

/* read grid coordinates */
    cout<<"Starting to read grid coordinates"<<endl;
    int ncoord=-1;
    cg_ncoords(index_file,index_base,index_zone,&ncoord);

    float *x=new float[irmax];
    float *y=new float[irmax];
    float *z=new float[irmax];

    cout<<"Read grid coordinates"<<endl;
    cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,x);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,y);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
                  CGNS_ENUMV(RealSingle),&irmin,&irmax,z);

    cout<<"The number of vertices is: "<<irmax<<" "<<ncoord<<endl;
/* find out how many sections */
    cg_nsections(index_file,index_base,index_zone,&nsections);
    printf("\nnumber of sections=%i\n",nsections);
/* read element connectivity */

    int nTriSum=0;

    //cgsize_t *triSection[3];
    cgsize_t **triSection=nullptr;
    triSection=new cgsize_t*[nsections-1];
    int *nTriSection=new int[nsections-1];
    //vector<vector<cgsize_t>>triSection;
    int count=0;

    for (index_sect=1; index_sect <= nsections; index_sect++)
    {
      cg_section_read(index_file,index_base,index_zone,index_sect,sectionname,
                      &itype,&istart,&iend,&nbndry,&iparent_flag);


      printf("\nReading section data...\n");
      printf("   section name=%s\n",sectionname);
      printf("   section type=%s\n",ElementTypeName[itype]);
      printf("   istart,iend=%i, %i\n",(int)istart,(int)iend);
      if (itype == CGNS_ENUMV(TETRA_4))
      {
        printf("   reading element data for this Tetras\n");
        cg_elements_read(index_file,index_base,index_zone,index_sect,tetras, \
                         &iparentdata);
      }

      if (itype == CGNS_ENUMV(TRI_3))
      {
          int nTri=iend-istart+1;
          nTriSection[count]=nTri;
          nTriSum+=nTri;
          triSection[count]=new cgsize_t[nTri*3];
          string strName=sectionname;
          bcstring.push_back(strName);
          printf("   reading element data for this element\n");
          cg_elements_read(index_file,index_base,index_zone,index_sect,triSection[count], \
                  &iparentdata);
          count++;
      }

    }
/* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully read unstructured grid from file grid_c.cgns\n");

    mesh.NumNodes=nNodes;
    mesh.nodes=new Node[nNodes];
    mesh.NumTetras=nTetras;
    mesh.pTetras=new TETRAS[nTetras];
    mesh.NumTris=nTriSum;
    mesh.pTris=new TRI[nTriSum];

    for(int i=0;i<nNodes;i++)
    {
        mesh.nodes[i].index=i;
        mesh.nodes[i].coord.x=x[i];
        mesh.nodes[i].coord.y=y[i];
        mesh.nodes[i].coord.z=z[i];
    }
    for(int i=0;i<nTetras;i++)
    {
        mesh.pTetras[i].index=i;
        mesh.pTetras[i].vertices[0]=tetras[i*4+0]-1;
        mesh.pTetras[i].vertices[1]=tetras[i*4+1]-1;
        mesh.pTetras[i].vertices[2]=tetras[i*4+2]-1;
        mesh.pTetras[i].vertices[3]=tetras[i*4+3]-1;
    }

    count=0;
   // mesh.pTris[i].index=i;
    for(int j=0;j<nsections-1;j++)
    {
        for(int k=0;k<nTriSection[j];k++)
        {
            mesh.pTris[count].index=count;
            mesh.pTris[count].vertices[0]=triSection[j][k*3+0]-1;
            mesh.pTris[count].vertices[1]=triSection[j][k*3+1]-1;
            mesh.pTris[count].vertices[2]=triSection[j][k*3+2]-1;
            mesh.pTris[count].iSurf=j;
            count++;
        }
    }

    assert(count==nTriSum);


    return 1;
}
#endif
