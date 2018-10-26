#ifndef ABS_ORI_HPP
#define ABS_ORI_HPP

#include<vector>
using namespace std;


//affine similarity macro
#define SIMILARITY_X(x,y,p) p[0]*x-p[1]*y+p[2]
#define SIMILARITY_Y(x,y,p) p[1]*x+p[0]*y+p[3]

typedef struct struct_GeoInfo
{
    double left, top;        //ground coordinate of top-left corner 
    double minlon, maxlat;  //(lon,lat) of top-left corner
    double dx, dy;           //resolution;
    int    wd, ht;           //dimension of image
    int    nband;
    int    zoneNumber;      //zone number 
    int    type;            //data type
    const char* projectRef; //
}stGeoInfo;

typedef struct stLon2Pix
{
	double lon, lat, height;
	int x, y;
}Lon2Pix;

typedef struct stPOINT3
{
	double x, y, z;
	bool   bIsOutSide;
	double f;
}POINT3, POINT3D;

//parameters for absolute pose estimation
typedef struct ABS_POS_PARAMS
{
	double scale;
	double R[9];
	double T[3];
}stAbsPOS;

int RansanAbsOri(vector<POINT3D> srcPts, vector<POINT3D> dstPts,
	stAbsPOS& absPosParams);
int AbsOriP3(vector<POINT3D> freePts, vector<POINT3D> grdPts,
	stAbsPOS& absPosParams);

int  SimRansac(vector<POINT3D>& srcPt, vector<POINT3D>& desPt, double *ap);

void SimilarityTransform(vector<POINT3D>& oriPt, vector<POINT3D>& desPt, double* p);


int  GdalWriteImageColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b,
    int ht, int wd, stGeoInfo geoInfo);


#endif // !ABS_ORI_HPP



