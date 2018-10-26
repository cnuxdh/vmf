/***************************************************************************
* Copyright(c) 2003 , TianYuan 3D Technology
* All rights reversed 
*
* 文件名称 : Global.h
* 文件标识 :
* 摘    要 : 常量、数据结构的定义
*
* 当前版本 : 1.0
* 作    者 : 李仁举
* 完成日期 : 
*
* 取代版本 : 
* 原 作 者 : 
* 完成日期 : 2003.12.10
/***************************************************************************/

#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <windows.h>
// 常量定义
#define PI 3.1415926

#define MAX_CIRCLE_NUM  200

/*
// 标定块
#define MAX_CIRCLE_NUM			200    // 标志点的最大个数
#define ROW						9      // 定标块上点的行数
#define COLUMN					11     // 定标块上点的列数
#define SPACING					30     // 定标块上标志点之间的X和Y方向间距
#define R_MIN                   6      // 标定块上小圆半径
#define R_MAX                   10     // 标定块上大圆半径

// CCD 
#define WIDTH					780    // 图像宽度
#define HEIGHT					582    // 图像高度
#define SIZE_X                  6.474  // CCD大小(x)
#define SIZE_Y                  4.8306  // CCD大小(y)
#define FOCUS                   12 

#define MAXARRAYBOUND			780*582+1  // 每图像象素个数+1
*/
// 测量点距
#define POINT_SPACING           0.5 
// 数据结构
// 2004.11.26 系统参数
struct SystemInfo
{
	unsigned int nWidth  ;  // 图像宽度
	unsigned int nHeight ;  // 图像高度

	double   dSizeX      ;  // CCD大小(x)
	double   dSizeY      ;  // CCD大小(x)
	double   dFocus      ;  // 镜头焦距大小

	double   dPointSpacing ;  // 测量点距大小

	unsigned int  nCalibType ; // 标定块类型，1：9行11列；2：7行9列
	unsigned int  nRow   ;  // 标定块上点的行数
	unsigned int  nColumn ; // 标定块上点的列数
	
	double   dCalibPointSpacing ; // 标定块上标志点的点距

	double   dRMin       ;  // 标定块上小圆半径
	double   dRMax       ;  // 标定块上大圆半径

	double   dRRefPoint  ;  // 标志点拼接中用的标志点大小
};


// 每个点的结构，在ANN, PostProc, Reconstruct 中用到
struct TVertex
{
	float x ;
	float y ;
	float z ;
	bool  bDelete ;    // 此点是否被删除
	bool  bBoundary ;  // 此点是否在两个点云的相交区域内，如果在，在哪个区域，不在－1

	// 法向量(nx, ny, nz ) , 法向量是否存在
	float nx ;
	float ny ;
	float nz ; 
	bool  bNormal ; 

	// 颜色(r,b,b), 颜色是否存在
	unsigned char r ;
	unsigned char g ;
	unsigned char b ;
	bool  bColor ;

	// 曲率，曲率是否存在
	float  fCurvature ;
	bool   bCurvature;
	// 对应于二维图象上的坐标
	// 右图象，一般只用右图象
	unsigned short u ;
	unsigned short v ;

    bool  bShow; // 该点是否显示
    unsigned short nRowNum ; // 该点的行号，用于在按比例均匀显示时行的控制
   
    unsigned char nGray ;  // 该点的灰度值
    bool          bGray ;  // 该点是否有灰度值
} ;        

// 三维点的结构
struct Point3D {
	double x ;
	double y ;
	double z ;
};

// 椭圆信息
struct CEllipseInfo{
	double cx ; // 椭圆中心x坐标
	double cy  ; // 椭圆中心y坐标
	double a  ; // 椭圆半长轴
	double b  ; // 椭圆半短轴
	double alpha ;      // 长轴的倾角
};

// 椭圆方程系数: x*x + b*x*y +c*y*y + d*x + e*y + f = 0 
struct CEllipseCoeff{
	double b ; // 椭圆中心x坐标
	double c  ; // 椭圆中心y坐标
	double d  ; // 椭圆半长轴
	double e  ; // 椭圆半短轴
	double f ;      // 长轴的倾角
};

// 两组图象上匹配点的结构
struct CMatchPoint{
	double lx ;  // 左图象， x坐标
	double ly ;  // 左图象， y坐标
	double rx ;  // 
	double ry ; 
	double lPhase; //  左图象上此点的相位
	double rPhase ; // 右图象上此点的相位
	double X  ;
	double Y  ; 
	double Z  ; 
};

// 网格点检测结构
struct CPoints
{
	float x;
	float y;
	int   col;
	int   row;
};
#endif // _GLOBAL_H