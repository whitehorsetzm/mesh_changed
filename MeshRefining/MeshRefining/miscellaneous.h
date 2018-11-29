#ifndef MISCELLANEOUS_H
#define MISCELLANEOUS_H

#include <string>
#include <vector>

/************************************************************************/
/*  文件及文件夹相关操作                                                   */
/************************************************************************/

/* 递归创建文件夹 */
// 返回值： 0 成功，1 失败
int createDictionary(std::string strDic);

/* 判断文件夹是否存在 */
bool isFileExist( std::string &path );

/* 创建文件夹 */
int makeDir( std::string strDic );



/************************************************************************/
/*  字符串相关操作                                                        */
/************************************************************************/

/**/
std::vector<std::string> splitStr( std::string &str, char div );

void trimLeftRight( std::string &str );

void insert_symbol( std::string &str, std::string letter );

std::string intToString(int value);


/************************************************************************/
/* 数字相关                                                             */
/************************************************************************/

/* iRegion[n]中是否存在iSearch. 存在返回true, 否则false */
bool findIntNum(int iSearch, const int* iRegion, int n);

/* find the index of minimize integer */
int findIdxMinInt(const int *iArray, int n);

/************************************************************************/
/* 几何相关                                                             */
/************************************************************************/

/* 根据4个点坐标判断体积正负 */
bool volumeTetrSign(const double *pnt1, const double *pnt2, const double *pnt3, const double *pnt4);

#endif // MISCELLANEOUS_H
