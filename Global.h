/***************************************************************************
* Copyright(c) 2003 , TianYuan 3D Technology
* All rights reversed 
*
* �ļ����� : Global.h
* �ļ���ʶ :
* ժ    Ҫ : ���������ݽṹ�Ķ���
*
* ��ǰ�汾 : 1.0
* ��    �� : ���ʾ�
* ������� : 
*
* ȡ���汾 : 
* ԭ �� �� : 
* ������� : 2003.12.10
/***************************************************************************/

#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <windows.h>
// ��������
#define PI 3.1415926

#define MAX_CIRCLE_NUM  200

/*
// �궨��
#define MAX_CIRCLE_NUM			200    // ��־���������
#define ROW						9      // ������ϵ������
#define COLUMN					11     // ������ϵ������
#define SPACING					30     // ������ϱ�־��֮���X��Y������
#define R_MIN                   6      // �궨����СԲ�뾶
#define R_MAX                   10     // �궨���ϴ�Բ�뾶

// CCD 
#define WIDTH					780    // ͼ����
#define HEIGHT					582    // ͼ��߶�
#define SIZE_X                  6.474  // CCD��С(x)
#define SIZE_Y                  4.8306  // CCD��С(y)
#define FOCUS                   12 

#define MAXARRAYBOUND			780*582+1  // ÿͼ�����ظ���+1
*/
// �������
#define POINT_SPACING           0.5 
// ���ݽṹ
// 2004.11.26 ϵͳ����
struct SystemInfo
{
	unsigned int nWidth  ;  // ͼ����
	unsigned int nHeight ;  // ͼ��߶�

	double   dSizeX      ;  // CCD��С(x)
	double   dSizeY      ;  // CCD��С(x)
	double   dFocus      ;  // ��ͷ�����С

	double   dPointSpacing ;  // ��������С

	unsigned int  nCalibType ; // �궨�����ͣ�1��9��11�У�2��7��9��
	unsigned int  nRow   ;  // �궨���ϵ������
	unsigned int  nColumn ; // �궨���ϵ������
	
	double   dCalibPointSpacing ; // �궨���ϱ�־��ĵ��

	double   dRMin       ;  // �궨����СԲ�뾶
	double   dRMax       ;  // �궨���ϴ�Բ�뾶

	double   dRRefPoint  ;  // ��־��ƴ�����õı�־���С
};


// ÿ����Ľṹ����ANN, PostProc, Reconstruct ���õ�
struct TVertex
{
	float x ;
	float y ;
	float z ;
	bool  bDelete ;    // �˵��Ƿ�ɾ��
	bool  bBoundary ;  // �˵��Ƿ����������Ƶ��ཻ�����ڣ�����ڣ����ĸ����򣬲��ڣ�1

	// ������(nx, ny, nz ) , �������Ƿ����
	float nx ;
	float ny ;
	float nz ; 
	bool  bNormal ; 

	// ��ɫ(r,b,b), ��ɫ�Ƿ����
	unsigned char r ;
	unsigned char g ;
	unsigned char b ;
	bool  bColor ;

	// ���ʣ������Ƿ����
	float  fCurvature ;
	bool   bCurvature;
	// ��Ӧ�ڶ�άͼ���ϵ�����
	// ��ͼ��һ��ֻ����ͼ��
	unsigned short u ;
	unsigned short v ;

    bool  bShow; // �õ��Ƿ���ʾ
    unsigned short nRowNum ; // �õ���кţ������ڰ�����������ʾʱ�еĿ���
   
    unsigned char nGray ;  // �õ�ĻҶ�ֵ
    bool          bGray ;  // �õ��Ƿ��лҶ�ֵ
} ;        

// ��ά��Ľṹ
struct Point3D {
	double x ;
	double y ;
	double z ;
};

// ��Բ��Ϣ
struct CEllipseInfo{
	double cx ; // ��Բ����x����
	double cy  ; // ��Բ����y����
	double a  ; // ��Բ�볤��
	double b  ; // ��Բ�����
	double alpha ;      // ��������
};

// ��Բ����ϵ��: x*x + b*x*y +c*y*y + d*x + e*y + f = 0 
struct CEllipseCoeff{
	double b ; // ��Բ����x����
	double c  ; // ��Բ����y����
	double d  ; // ��Բ�볤��
	double e  ; // ��Բ�����
	double f ;      // ��������
};

// ����ͼ����ƥ���Ľṹ
struct CMatchPoint{
	double lx ;  // ��ͼ�� x����
	double ly ;  // ��ͼ�� y����
	double rx ;  // 
	double ry ; 
	double lPhase; //  ��ͼ���ϴ˵����λ
	double rPhase ; // ��ͼ���ϴ˵����λ
	double X  ;
	double Y  ; 
	double Z  ; 
};

// �������ṹ
struct CPoints
{
	float x;
	float y;
	int   col;
	int   row;
};
#endif // _GLOBAL_H