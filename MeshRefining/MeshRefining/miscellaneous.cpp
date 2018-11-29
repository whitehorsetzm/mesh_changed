#include "miscellaneous.h"

#include <vector>
#include <sstream>
#include <cstdlib>
#include <iostream>

#include <string>
#include <stdint.h>
#ifdef WIN32
//#include <windows.h>
#include <io.h>

//
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#endif

using namespace std;
#define MAX_PATH_LEN 2048

/************************************************************************/
/*  文件及文件夹相关操作                                                   */
/************************************************************************/
#if 1
/* 递归创建文件夹 */
// 返回值： 0 成功，1 失败
#ifdef WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

// 从左到右依次判断文件夹是否存在,不存在就创建
// example: /home/root/mkdir/1/2/3/4/
// 注意:最后一个如果是文件夹的话,需要加上 '\' 或者 '/'
int createDictionary(std::string directoryPath)
{
    uint32_t dirPathLen = directoryPath.length();
    if (dirPathLen > MAX_PATH_LEN)
    {
        return -1;
    }
    char tmpDirPath[MAX_PATH_LEN] = { 0 };
    for (uint32_t i = 0; i < dirPathLen; ++i)
    {
        tmpDirPath[i] = directoryPath[i];
        if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
        {
            if (ACCESS(tmpDirPath, 0) != 0)
            {
                int32_t ret = MKDIR(tmpDirPath);
                if (ret != 0)
                {
                    return ret;
                }
            }
        }
    }
    return 0;
}
#else
int createDictionary(std::string strDic)
{
	int rtn = 0;
    if( isFileExist( strDic ) )
        return 0;
#if 1
    std::string dir = "\"" + strDic + "\"";
    rtn = makeDir( dir );
#else
    std::vector<std::string> str_vec = splitStr( strDic, '\\' );
    int i = 0;
    std::string str_tmp = str_vec[0];
    for(i=1;i<str_vec.size();i++)
    {
        str_tmp += std::string("\\") + str_vec[i];
        /* 非零表示成功，零表示失败 */
        if(CreateDirectory(str_tmp.c_str(),NULL) == 0)
            return 1; // fail
    }
#endif
    return rtn;
}
#endif

/* 判断文件夹是否存在 */
bool isFileExist( std::string &path )
{
    int error = 0;
    //error = _access( path.c_str(), 0 );
    error = access( path.c_str(), 0 );
    if( error == 0 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

/* 创建文件夹 */
int makeDir( std::string strDic )
{
    int rtn;
    std::string mkdirCmd;
    std::string cmd;
#ifdef WIN32
    mkdirCmd = "mkdir";
#else
    mkdirCmd = "mkdir -p";
#endif

    cmd = mkdirCmd + " " + strDic;

    rtn = system( cmd.c_str() );

    return rtn;
}

/************************************************************************/
/*  字符串相关操作                                                        */
/************************************************************************/


std::vector<std::string> splitStr( std::string &str, char div )
{
    trimLeftRight( str );
    char tmp_1[2];
    tmp_1[0] = div;
    tmp_1[1] = '\0';
    insert_symbol( str, std::string( tmp_1 ) );
    int i = 0;
    int length = str.length();
    std::vector<std::string> result;
    result.clear();
    int posL = 0;
    for(i=0;i<length;i++)
    {
        if( str[i] == div )
        {
            if( i == 0 ) continue;
            result.push_back( str.substr( posL + 1, i - posL - 1 ) );
            posL = i;
        }
    }
    return result;
}

void trimLeftRight( std::string &str )
{
    int pos = 0;
    int lengths = str.length();
    if( lengths <= 0 ) return;
    int countSpace = 0;
    char spacedw = ' ';
    int i = -1;
    while( str[++i] == spacedw )
    {
        countSpace ++;
    }
    str.erase( pos, countSpace );

    lengths = str.length();
    if( lengths <= 0 ) return;
    i = lengths;
    countSpace = 0;
    while( str[--i] == spacedw )
    {
        countSpace ++;
    }
    str.erase( i+1, countSpace );
}

void insert_symbol( std::string &str, std::string letter )
{
    if( str == "" ) return;
    str.insert( 0, letter );
    std::string tmp = letter;
    str += tmp;
}


std::string intToString(int value)
{
	string str;
	stringstream stream;
	stream << value;
	stream >> str;

	return str;
}


bool findIntNum(int iSearch, const int* iRegion, int n)
{
	for (int i = 0; i < n; ++i)
	{
		if (iSearch == iRegion[i])
		{
			return true;
		}
	}
	return false;
}

/* find the index of minimize integer */
int findIdxMinInt(const int *iArray, int n)
{
    if(n <= 0)
        return -1;

    int iMin = 0;
    for(int i = 1; i < n; ++i)
    {
        if(iArray[iMin]>iArray[i])
            iMin = i;
    }

    return iMin;

}

/* 根据4个点坐标求得四面体体积 */
bool volumeTetrSign(const double *pnt1, const double *pnt2, const double *pnt3, const double *pnt4)
{
	double x21,y21,z21,x31,y31,z31,x41,y41,z41,vol;
	x21 = pnt2[0] - pnt1[0];
	y21 = pnt2[1] - pnt1[1];
	z21 = pnt2[2] - pnt1[2];
	x31 = pnt3[0] - pnt1[0];
	y31 = pnt3[1] - pnt1[1];
	z31 = pnt3[2] - pnt1[2];
	x41 = pnt4[0] - pnt1[0];
	y41 = pnt4[1] - pnt1[1];
	z41 = pnt4[2] - pnt1[2];
	vol = x21*y31*z41 + y21*z31*x41 + z21*x31*y41 - z21*y31*x41 - y21*x31*z41 - x21*z31*y41;

	if(vol<0.000)
		return false;
	return true;
}
