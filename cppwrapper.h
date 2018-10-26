
#ifndef WRAPPER_H  
#define WRAPPER_H  

#include <opencv\cv.h>
#include <opencv\highgui.h>

#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/imgproc/imgproc.hpp> 
#include <opencv2/features2d/features2d.hpp>



#include<vector>
using namespace std;

using namespace cv;

#ifdef _WIN32
# define DLL_EXPORT __declspec( dllexport )
#else
# define DLL_EXPORT
#endif


#define HIST_LEN 128

//typedef struct structHist
//{
//	double h[HIST_LEN];
//	int    sum;
//}stHIST;

class stHIST
{
public:
	stHIST(){ memset(h, 0, sizeof(double)*HIST_LEN); }
	~stHIST(){}
	
	double h[HIST_LEN];
	int    sum;

private:
};


//
//#ifdef __cplusplus  
//extern "C" {
//#endif  



DLL_EXPORT void    PanoVideoMatchAndWrap(IplImage* pPanoPatch, IplImage* pFrame, double* hp);
DLL_EXPORT int     CalculateGradientHist(IplImage* pSrc, double* pHist, int nHistLen);
DLL_EXPORT int     HistMatch(stHIST* pHists, int nHistNum, stHIST* pFrameHist, double* cov);
DLL_EXPORT int     CalcPanoHoG(IplImage* pPano, stHIST* pHists, int nHistNum, IplImage** pImagePatchs);
DLL_EXPORT double  Cov(double* v1, double *v2, int len);

void GenerateBlockHist1(Mat src, stHIST* pHistSrc, int bht, int bwd, int blocksize);
void GenerateBlockHist(Mat src, vector<stHIST>& histSrc, int bht, int bwd, int blocksize);

int     HistMatch(vector<stHIST>& pHists, int nHistNum, stHIST& pFrameHist, double* cov);



void    BenchToPano(IplImage* pPano, IplImage* pBench, double* sp);

//#ifdef __cplusplus  
//};
//#endif   
#endif