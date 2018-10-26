
#include"changedetect.h"
#include"cppwrapper.h"
#include"mymatrix.h"



#include <opencv\cv.h>
#include <opencv\highgui.h>

using namespace std;
using namespace cv;


typedef struct stMyPointF
{
	double x, y;
}MyPointF;

//affine transform macro
#define AFFINE_X(x,y,p) p[0]*x+p[1]*y+p[2]
#define AFFINE_Y(x,y,p) p[3]*x+p[4]*y+p[5]

//similarity transform
#define SIMILARITY_X(x,y,p) p[0]*x-p[1]*y+p[2]
#define SIMILARITY_Y(x,y,p) p[1]*x+p[0]*y+p[3]

//homography
#define HOMOGRAPHY_X(x,y,p) (p[0]*x+p[1]*y+p[2]) / (p[6]*x+p[7]*y+p[8]) 
#define HOMOGRAPHY_Y(x,y,p) (p[3]*x+p[4]*y+p[5]) / (p[6]*x+p[7]*y+p[8])


/*function from 5point lib to select index randomly
input:
n: the length of arr
k: number of random
arr: array saving the random values
*/
void ChooseRandIndex(int n, int k, int *arr)
{
	int i;

	if (k > n)
	{
		fprintf(stderr, "[choose] Error: k > n\n");
		return;
	}

	for (i = 0; i < k; i++)
	{
		while (1)
		{
			int idx = rand() % n;
			int j, redo = 0;
			for (j = 0; j < i; j++)
			{
				if (idx == arr[j])
				{
					redo = 1;
					break;
				}
			}
			if (!redo)
			{
				arr[i] = idx;
				break;
			}
		}
	}
}


/* 
x' = ax - by + c
y' = bx + ay + d;
*/
void SimilarityTransform(MyPointF* oriPt, MyPointF* desPt, int num, double* p)
{
	double* A = NULL;
	double* At = NULL;
	double* AtA = NULL;
	double* L = NULL;
	double* AtL = NULL;
	double* res = NULL;
	int  paraNum = 4;
	int  i, j;

	A = (double*)malloc(num * 2 * paraNum*sizeof(double));
	At = (double*)malloc(num * 2 * paraNum*sizeof(double));
	L = (double*)malloc(num * 2 * sizeof(double));
	AtA = (double*)malloc(paraNum*paraNum*sizeof(double));
	AtL = (double*)malloc(paraNum*sizeof(double));
	res = (double*)malloc(paraNum*sizeof(double));

	for (i = 0; i<num; i++)
	{
		A[2 * i*paraNum] = oriPt[i].x;
		A[2 * i*paraNum + 1] = -oriPt[i].y;
		A[2 * i*paraNum + 2] = 1;
		A[2 * i*paraNum + 3] = 0;

		A[(2 * i + 1)*paraNum] = oriPt[i].y;
		A[(2 * i + 1)*paraNum + 1] = oriPt[i].x;
		A[(2 * i + 1)*paraNum + 2] = 0;
		A[(2 * i + 1)*paraNum + 3] = 1;

		L[2 * i] = desPt[i].x;
		L[2 * i + 1] = desPt[i].y;
	}

	mytranspose(A, At, num * 2, paraNum);
	mymult(At, A, AtA, paraNum, num * 2, paraNum);
	mymult(At, L, AtL, paraNum, num * 2, 1);
	myinvers_matrix(AtA, paraNum);
	mymult(AtA, AtL, res, paraNum, paraNum, 1);

	for (i = 0; i<paraNum; i++)
	{
		p[i] = res[i];
	}

	free(A);
	free(At);
	free(L);
	free(AtA);
	free(AtL);
	free(res);
}

/* 仿射变换
x' = ax+by+c;
y' = dx+ey+f;
输入参数：
op：变换前的坐标
dp：变换后的坐标
num：点的数量
输出：
p：  仿射系数
*/
void AffineTransformFloat(MyPointF* op, MyPointF* dp, int num, float* p)
{
	int i, j, k;
	double res[6];
	double A[2][6];
	double L[2];
	double At[6][2];
	double AtA[6][6];
	double AtL[6];
	double sAtA[6][6];
	double sAtL[6];

	for (j = 0; j<6; j++)
	{
		for (i = 0; i<6; i++)
		{
			sAtA[j][i] = 0;
		}
		sAtL[j] = 0;
	}

	for (k = 0; k<num; k++)
	{
		A[0][0] = op[k].x;
		A[0][1] = op[k].y;
		A[0][2] = 1;
		A[0][3] = 0;
		A[0][4] = 0;
		A[0][5] = 0;

		A[1][0] = 0;
		A[1][1] = 0;
		A[1][2] = 0;
		A[1][3] = op[k].x;
		A[1][4] = op[k].y;
		A[1][5] = 1;

		L[0] = dp[k].x;
		L[1] = dp[k].y;

		mytranspose(*A, *At, 2, 6);
		mymult(*At, *A, *AtA, 6, 2, 6);
		mymult(*At, L, AtL, 6, 2, 1);

		for (j = 0; j<6; j++)
		{
			for (i = 0; i<6; i++)
			{
				sAtA[j][i] += AtA[j][i];
			}
			sAtL[j] += AtL[j];
		}
	}

	myinvers_matrix(*sAtA, 6);
	mymult(*sAtA, sAtL, res, 6, 6, 1);

	for (i = 0; i<6; i++)
		p[i] = res[i];
}

int  AffineRansac(MyPointF* desPt, MyPointF* srcPt, int npt, double *ap)
{
	int ransac_rounds = 500;
	int maxInliers = 0;
	int DistanceThrshold = 16;

#define RANSAC_PT_NUM 4

	for (int round = 0; round < ransac_rounds; round++)
	{
		int indices[RANSAC_PT_NUM];
		ChooseRandIndex(npt, RANSAC_PT_NUM, indices);

		MyPointF dpts[RANSAC_PT_NUM];
		MyPointF spts[RANSAC_PT_NUM];
		for (int k = 0; k < RANSAC_PT_NUM; k++)
		{
			int randi = indices[k];
			dpts[k].x = desPt[randi].x;
			dpts[k].y = desPt[randi].y;
			spts[k].x = srcPt[randi].x;
			spts[k].y = srcPt[randi].y;
		}

		double p[6];
		SimilarityTransform(spts, dpts, RANSAC_PT_NUM, p);
		//float p[6];
		//AffineTransformFloat(spts, dpts, 3, p);

		//error
		double error = 0;
		int inliers = 0;
		//vInliers.clear();
		for (int i = 0; i < npt; i++)
		{
			double tx, ty;
			double dx = desPt[i].x;
			double dy = desPt[i].y;
			double ox = srcPt[i].x;
			double oy = srcPt[i].y;

			//tx = AFFINE_X(ox, oy, p);
			//ty = AFFINE_Y(ox, oy, p);

			tx = SIMILARITY_X(ox, oy, p);
			ty = SIMILARITY_Y(ox, oy, p);

			double len = sqrt((tx - dx)*(tx - dx) + (ty - dy)*(ty - dy));
			if (len < DistanceThrshold)
			{
				//vInliers.push_back(i);
				inliers++;
			}
			error += len;
		}
		error /= double(npt);

		if (maxInliers < inliers)
		{
			maxInliers = inliers;

			for (int k = 0; k < 6; k++) {
				ap[k] = p[k];
			}
		}
	}

	return maxInliers;

}

//find homography using opencv
Mat  HomographyRansac(MyPointF* desPt, MyPointF* srcPt, int npt, double *hp){

	//-- Localize the object
	std::vector<Point2f> despts;
	std::vector<Point2f> srcpts;

	for (int i = 0; i < npt; i++)
	{
		Point2f dp, sp;

		dp.x = desPt[i].x;
		dp.y = desPt[i].y;
		sp.x = srcPt[i].x;
		sp.y = srcPt[i].y;
		
		despts.push_back(dp);
		srcpts.push_back(sp);
	}

	Mat H = findHomography(srcpts, despts, CV_RANSAC, 8);
	int index = 0;
	for (int j = 0; j < 3; j++){
		for (int i = 0; i < 3; i++){
			hp[index] = H.at<double>(j, i);
			index++;
		}
	}		

	return H;
}


int CalcPanoHists(IplImage* pPano, stHIST* pHists, int nHistNum){

	IplImage* pGray = pPano;
	if (pPano->nChannels == 3){
		pGray = cvCreateImage(cvGetSize(pPano), IPL_DEPTH_8U, 1);
		cvCvtColor(pPano, pGray, CV_BGR2GRAY);//cvCvtColor(src,des,CV_BGR2GRAY
	}

	int ht = pGray->height;
	int wd = pGray->width;
	int swd = pGray->widthStep;

	if (ht > wd)
		return -1;

	int step = (wd - ht) / (nHistNum - 1);

	int subSize = ht*0.5;

	for (int i = 0; i < nHistNum; i++){

		int cx = subSize + i*step;
		int cy = subSize;
		int l = max(0, cx - subSize);
		int r = min(wd - 1, cx + subSize);
		int t = 0;
		int b = ht - 1;

		double h[HIST_LEN];
		memset(h, 0, sizeof(double)*HIST_LEN);
		int interval = 256 / HIST_LEN;
		for (int m = t; m < b; m++)
			for (int n = l; n < r; n++){
				unsigned char pv = (unsigned char)(pGray->imageData[m*swd + n]);
				int index = (double)(pv) / (double)(interval)+0.5;
				index = min(HIST_LEN - 1, index);
				h[index]++;
			}
		//normalize the hist
		double sum = 0;
		for (int k = 0; k < HIST_LEN; k++){
			sum += h[k];
		}
		for (int k = 0; k < HIST_LEN; k++){
			h[k] = h[k] / sum;
		}

		memcpy(pHists[i].h, h, sizeof(double)*HIST_LEN);
	}

	cvReleaseImage(&pGray);

	return 0;
}

void CalcFrameHist(IplImage* pFrame, stHIST* pHist){

	IplImage* pGray = pFrame;
	if (pFrame->nChannels == 3){
		pGray = cvCreateImage(cvGetSize(pFrame), IPL_DEPTH_8U, 1);
		cvCvtColor(pFrame, pGray, CV_BGR2GRAY);//cvCvtColor(src,des,CV_BGR2GRAY
	}

	int ht = pGray->height;
	int wd = pGray->width;
	int swd = pGray->widthStep;

	double h[HIST_LEN];
	memset(h, 0, sizeof(double)*HIST_LEN);
	int interval = 256 / HIST_LEN;
	for (int j = 0; j < ht; j++){
		for (int i = 0; i < wd; i++){
			unsigned char pv = (unsigned char)(pGray->imageData[j*swd + i]);
			int index = (double)(pv) / (double)(interval)+0.5;
			index = min(HIST_LEN - 1, index);
			h[index]++;
		}
	}

	//normalize the hist
	double sum = 0;
	for (int k = 0; k < HIST_LEN; k++){
		sum += h[k];
	}
	for (int k = 0; k < HIST_LEN; k++){
		h[k] = h[k] / sum;
	}
	memcpy(pHist->h, h, sizeof(double)*HIST_LEN);

	cvReleaseImage(&pGray);
} 

double  Cov(double* v1, double *v2, int len)
{
	int i, j;
	float mx;
	float my;
	float dx2, dy2, dxy;
	float fcov;

	mx = 0;
	my = 0;
	for (i = 0; i<len; i++)
	{
		mx += v1[i];
		my += v2[i];
	}
	mx = mx / len;
	my = my / len;

	dx2 = 0;
	dy2 = 0;
	dxy = 0;
	for (i = 0; i<len; i++)
	{
		dx2 += (v1[i] - mx)*(v1[i] - mx);
		dy2 += (v2[i] - my)*(v2[i] - my);
		dxy += (v1[i] - mx)*(v2[i] - my);
	}
	fcov = dxy / (sqrt(dx2)*sqrt(dy2));
	return fcov;
}

int HistMatch(vector<stHIST>& pHists, int nHistNum, stHIST& pFrameHist, double* cov)
{
	int nIndex = 0;

	//printf("hist match.... \n");
	double maxcov = -1;
	for (int i = 0; i < nHistNum; i++){

		double cor = Cov(pHists[i].h, pFrameHist.h, HIST_LEN);
		//printf("%lf ", cor);
		if (maxcov < cor){
			maxcov = cor;
			nIndex = i;
		}
	}
	//printf(" \n");

	*cov = maxcov;

	return nIndex;
}

int HistMatch(stHIST* pHists, int nHistNum, stHIST* pFrameHist, double* cov){

	int nIndex = 0;

	//printf("hist match.... \n");
	double maxcov = -1;
	for (int i = 0; i < nHistNum; i++){

		double cor = Cov(pHists[i].h, pFrameHist->h, HIST_LEN);
		//printf("%lf ", cor);
		if (maxcov < cor){
			maxcov = cor;
			nIndex = i;
		}
	}
	printf(" \n");

	*cov = maxcov;

	return nIndex;
}

#define PI 3.1415926

//!computes the angle in degrees (0..360) of the vector (x,y), extracted from OpenCV
float cvFastAtan2(double y, double x)
{
	double a, x2 = (double)x*x, y2 = (double)y*y;
	if (y2 <= x2)
	{
		a = (180. / PI)*x*y*(x2 + 0.43157974*y2) / (x2*x2 + y2*(0.76443945*x2 + 0.05831938*y2) + DBL_EPSILON);
		return (float)(x < 0 ? a + 180 : y >= 0 ? a : 360 + a);
	}
	a = (180. / PI)*x*y*(y2 + 0.43157974*x2) / (y2*y2 + x2*(0.76443945*y2 + 0.05831938*x2) + DBL_EPSILON);

	return (float)(y >= 0 ? 90 - a : 270 - a);
}

int CalculateGradientHist(IplImage* pSrc, double* pHist, int nHistLen)
{
	IplImage* pImage = cvCloneImage(pSrc);
	if (pSrc->nChannels == 3){
		cvReleaseImage(&pImage);
		pImage = cvCreateImage(cvGetSize(pSrc), IPL_DEPTH_8U, 1);
		cvCvtColor(pSrc, pImage, CV_BGR2GRAY);//cvCvtColor(src,des,CV_BGR2GRAY
	}

	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;

	//2. sobel gradient
	IplImage* pDx = cvCreateImage(cvGetSize(pImage), IPL_DEPTH_16S, 1);
	cvSobel(pImage, pDx, 0, 1, 3);

	IplImage* pDy = cvCreateImage(cvGetSize(pImage), IPL_DEPTH_16S, 1);
	cvSobel(pImage, pDy, 1, 0, 3);

	//cvSaveImage("c:\\temp\\dx.jpg", pDx);
	//cvSaveImage("c:\\temp\\dy.jpg", pDy);

	//3. gradient histogram
	double sum = 0;
	for (int j = 0; j < ht; j++)
		for (int i = 0; i < wd; i++)
		{
			double dx = pDx->imageData[j*scanwd + i];
			double dy = pDy->imageData[j*scanwd + i];

			if (sqrt(dx*dx + dy*dy) < 24)
				continue;

			double angle = cvFastAtan2(dy, dx);
			//if (angle > 180) angle = angle - 180;

			double fIndex = angle / 3.0;
			int    index = (int)(fIndex);
			double w = fIndex - index;

			pHist[index] += 1 - w;
			pHist[index + 1] += w;

			sum++;
		}

	for (int i = 0; i < nHistLen; i++)
	{
		pHist[i] /= sum;
	}

	

	cvReleaseImage(&pImage);
	cvReleaseImage(&pDx);
	cvReleaseImage(&pDy);

	return int(sum);
}


int CalcPanoHoG(IplImage* pPano, stHIST* pHists, int nHistNum, IplImage** pImagePatchs){

	IplImage* pGray = cvCloneImage(pPano);
	if (pPano->nChannels == 3){
		cvReleaseImage(&pGray);
		pGray = cvCreateImage(cvGetSize(pPano), IPL_DEPTH_8U, 1);
		cvCvtColor(pPano, pGray, CV_BGR2GRAY);//cvCvtColor(src,des,CV_BGR2GRAY
	}

	int ht = pGray->height;
	int wd = pGray->width;
	int swd = pGray->widthStep;

	if (ht > wd)
		return -1;

	int step = (wd - ht) / (nHistNum - 1);

	int subSize = ht*0.5*1.7;

	for (int i = 0; i < nHistNum; i++){

		int cx = subSize + i*step;
		int cy = subSize;
		int l = max(0, cx - subSize);
		int r = min(wd - 1, cx + subSize);
		int t = 0;
		int b = ht - 1;

		//crop image
		CvSize size = cvSize(r - l, b - t);
		cvSetImageROI(pGray, cvRect(l, t, size.width, size.height));
		IplImage* pDest = cvCreateImage(size, pGray->depth, pGray->nChannels);
		cvCopy(pGray, pDest, NULL);
		cvResetImageROI(pGray);

		//save the image patchs
		cvSetImageROI(pPano, cvRect(l, t, size.width, size.height));
		IplImage* pDestColor = cvCreateImage(size, pPano->depth, pPano->nChannels);
		cvCopy(pPano, pDestColor, NULL);
		cvResetImageROI(pPano);
		pImagePatchs[i] = cvCloneImage(pDestColor);
		cvReleaseImage(&pDestColor);


		//
		CalculateGradientHist(pDest, pHists[i].h, 62);

		//char file[256];
		//sprintf(file, "c:\\temp\\crop_%d.jpg", i);
		//cvSaveImage(file, pImagePatchs[i], NULL);

		//memcpy(pHists[i].h, h, sizeof(double)*HIST_LEN);

		cvReleaseImage(&pDest);
	}

	cvReleaseImage(&pGray);
}


void CalculateGradientAngle(Mat image, double* pGradientAngle)
{
	IplImage src = image;   //cvCloneImage(pSrc);

	IplImage* pImage = cvCloneImage(&src);
	if (src.nChannels == 3)
	{	
		cvReleaseImage(&pImage);
		pImage = cvCreateImage(cvGetSize(&src), IPL_DEPTH_8U, 1);
		cvCvtColor(&src, pImage, CV_BGR2GRAY);//cvCvtColor(src,des,CV_BGR2GRAY
	}

	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;

	//2. sobel gradient
	IplImage* pDx = cvCreateImage(cvGetSize(pImage), IPL_DEPTH_16S, 1);
	cvSobel(pImage, pDx, 0, 1, 3);

	IplImage* pDy = cvCreateImage(cvGetSize(pImage), IPL_DEPTH_16S, 1);
	cvSobel(pImage, pDy, 1, 0, 3);

	//3. gradient histogram
	//double sum = 0;
	for (int j = 0; j < ht; j++)
		for (int i = 0; i < wd; i++)
		{
			double dx = pDx->imageData[j*scanwd + i];
			double dy = pDy->imageData[j*scanwd + i];

			if (sqrt(dx*dx + dy*dy) < 16)
				continue;

			double angle = cvFastAtan2(dy, dx);
			//if (angle > 180) angle = angle - 180;

			//double fIndex = angle / 3.0;
			//int    index = (int)(fIndex);
			//double w = fIndex - index;
			//pHist[index] += 1 - w;
			//pHist[index + 1] += w;
			//sum++;
			
			pGradientAngle[j*wd + i] = angle;
		}

	cvReleaseImage(&pImage);
}
 
//firstly calculate the angle, then hist
void GenerateBlockHist(Mat src, stHIST* pHistSrc, int bht, int bwd, int blocksize)
{
	int ht = src.rows;
	int wd = src.cols;

	double* pAnglesSrc = new double[ht*wd];
	for (int i = 0; i < ht*wd; i++)
	{
		pAnglesSrc[i] = -1;
	}
	CalculateGradientAngle(src, pAnglesSrc);
	
	for (int j = 0; j < bht; j++)
	{
		for (int i = 0; i < bwd; i++)
		{
			/*Rect r;
			r.x = i*blocksize;
			r.y = j*blocksize;
			r.width  = blocksize;
			r.height = blocksize;*/

			int l = i*blocksize;
			int r = min(wd, l + blocksize);
			int t = j*blocksize;
			int b = min(ht, t + blocksize);

			double sum = 0;
			for (int jj = t; jj < b; jj++)
				for (int ii = l; ii < r; ii++)
				{
					double angle = pAnglesSrc[jj*wd + ii];
					if (angle < 0)
						continue;

					double fIndex = angle / 3.0;
					int    index = (int)(fIndex);
					double w = fIndex - index;

					pHistSrc[j*bwd + i].h[index] += 1 - w;
					pHistSrc[j*bwd + i].h[index + 1] += w;
					sum++;
				}

			for (int k = 0; k < HIST_LEN; k++)
			{
				pHistSrc[j*bwd + i].h[k] /= sum;
			}
			pHistSrc[j*bwd + i].sum = sum;

			//Mat srcBlock = src(r);
			//Mat dstBlock = dst(r);
			//imwrite("c:\\temp\\block-src.jpg", srcBlock);
			//imwrite("c:\\temp\\block-dst.jpg", dstBlock);
		}
	}

	delete[] pAnglesSrc;
}

//firstly seperate the image into blocks, and then calculate the histogram of gradient
void GenerateBlockHist(Mat src, vector<stHIST>& histSrc, int bht, int bwd, int blocksize)
{
	int ht = src.rows;
	int wd = src.cols;

	for (int j = 0; j < bht; j++)
	{
		for (int i = 0; i < bwd; i++)
		{
			Rect r;
			r.x = i*blocksize;
			r.y = j*blocksize;
			r.width = blocksize;
			r.height = blocksize;

			Mat srcBlock = src(r);
			IplImage psrcblock = srcBlock;
			int texPixelSum = CalculateGradientHist(&psrcblock, histSrc[j*bwd + i].h, HIST_LEN);
			histSrc[j*bwd + i].sum = texPixelSum;
		}
	}
}


void GenerateBlockHist1(Mat src, stHIST* pHistSrc, int bht, int bwd, int blocksize)
{
	int ht = src.rows;
	int wd = src.cols;

	for (int j = 0; j < bht; j++)
	{
		for (int i = 0; i < bwd; i++)
		{
			Rect r;
			r.x = i*blocksize;
			r.y = j*blocksize;
			r.width  = blocksize;
			r.height = blocksize;

			Mat srcBlock = src(r);
			IplImage psrcblock = srcBlock;
			int texPixelSum = CalculateGradientHist(&psrcblock, pHistSrc[j*bwd+i].h, HIST_LEN);
			pHistSrc[j*bwd + i].sum = texPixelSum;
		}
	}
}

void BenchToPano(IplImage* pPano, IplImage* pBench, double* sp)
{
	ORB orb;
	vector<KeyPoint> kpPanoPatch, kpFrame; //keypoints
	Mat dpPanoPatch, dpFrame; //descrptors
	Mat mPanoPatch(pPano);
	Mat mFrame(pBench);
	
	orb(mPanoPatch, Mat(), kpPanoPatch, dpPanoPatch);
	orb(mFrame, Mat(), kpFrame, dpFrame);

	BFMatcher matcher(NORM_HAMMING);
	vector<DMatch> matches;
	matcher.match(dpPanoPatch, dpFrame, matches);

	double max_dist = 0; double min_dist = 100;
	//-- Quick calculation of max and min distances between keypoints
	for (int i = 0; i < dpPanoPatch.rows; i++)
	{
		double dist = matches[i].distance;
		if (dist < min_dist) min_dist = dist;
		if (dist > max_dist) max_dist = dist;
	}
	//printf("-- Max dist : %f \n", max_dist);
	//printf("-- Min dist : %f \n", min_dist);
	//-- Draw only "good" matches (i.e. whose distance is less than 0.6*max_dist )
	//-- PS.- radiusMatch can also be used here.
	std::vector< DMatch > good_matches;
	for (int i = 0; i < dpPanoPatch.rows; i++)
	{
		if (matches[i].distance < 0.6*max_dist)
		{
			good_matches.push_back(matches[i]);
		}
	}
		
	//calculate the affine parameters based on RANSAC
	int nmatch = good_matches.size();
	MyPointF* pSrcPts = new MyPointF[nmatch];
	MyPointF* pDstPts = new MyPointF[nmatch];
	for (int i = 0; i < nmatch; i++){
		int si = good_matches[i].queryIdx;
		int di = good_matches[i].trainIdx;
		pSrcPts[i].x = kpPanoPatch[si].pt.x;
		pSrcPts[i].y = kpPanoPatch[si].pt.y;
		pDstPts[i].x = kpFrame[di].pt.x;
		pDstPts[i].y = kpFrame[di].pt.y;
	}

	AffineRansac(pSrcPts, pDstPts, nmatch, sp);
}


void PanoVideoMatchAndWrap(IplImage* pPanoPatch, IplImage* pFrame, double* hp)
{
	//int   nFeatures = 5000;   // fSettings["ORBextractor.nFeatures"];
	//float fScaleFactor = 1.2; // fSettings["ORBextractor.scaleFactor"];
	//int   nLevels = 3;        // fSettings["ORBextractor.nLevels"];
	//int   fIniThFAST = 20;    // fSettings["ORBextractor.iniThFAST"];
	//int   fMinThFAST = 7;     // fSettings["ORBextractor.minThFAST"];
	//ORBextractor* pIniORBextractor = new ORBextractor(nFeatures,
	//	fScaleFactor, nLevels, fIniThFAST, fMinThFAST);


	ORB orb;
	
	vector<KeyPoint> kpPanoPatch, kpFrame; //keypoints
	Mat dpPanoPatch, dpFrame; //descrptors
	
	Mat mPanoPatch(pPanoPatch);
	Mat mFrame(pFrame);
		
	//resize the frame
	//double scale = (double)(mPanoPatch.rows) / (double)( mFrame.rows);
	//Size dsize = Size(mFrame.cols*scale, mFrame.rows*scale);
	//Mat resizeFrame = Mat(dsize, CV_8U);
	//resize(mFrame, resizeFrame, dsize);

	//imwrite("c:\\temp\\src.jpg", resizeFrame);
	//imwrite("c:\\temp\\dst.jpg", mFrame);
	
	//orb.detect(panoPatch, kpPanoPatch );
	//orb.detect(frame, kpFrame);
	orb(mPanoPatch, Mat(), kpPanoPatch, dpPanoPatch);
	orb(mFrame, Mat(), kpFrame, dpFrame);

	//BruteForceMatcher<HammingLUT> matcher;
	BFMatcher matcher(NORM_HAMMING);
	vector<DMatch> matches;
	matcher.match(dpPanoPatch, dpFrame, matches);

	double max_dist = 0; double min_dist = 100;
	//-- Quick calculation of max and min distances between keypoints
	for (int i = 0; i < dpPanoPatch.rows; i++)
	{
		double dist = matches[i].distance;
		if (dist < min_dist) min_dist = dist;
		if (dist > max_dist) max_dist = dist;
	}
	//printf("-- Max dist : %f \n", max_dist);
	//printf("-- Min dist : %f \n", min_dist);
	//-- Draw only "good" matches (i.e. whose distance is less than 0.6*max_dist )
	//-- PS.- radiusMatch can also be used here.
	std::vector< DMatch > good_matches;
	for (int i = 0; i < dpPanoPatch.rows; i++)
	{
		if (matches[i].distance < 0.6*max_dist)
		{
			good_matches.push_back(matches[i]);
		}
	}

	/*
	Mat img_matches;
	drawMatches(mPanoPatch, kpPanoPatch, mFrame, kpFrame,
		good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
		vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	*/
	//imwrite("c:\\temp\\orb-match.jpg", img_matches);
	
	//calculate the affine parameters based on RANSAC
	int nmatch = good_matches.size();
	MyPointF* pSrcPts = new MyPointF[nmatch];
	MyPointF* pDstPts = new MyPointF[nmatch];
	for (int i = 0; i < nmatch; i++){
		int si = good_matches[i].queryIdx;
		int di = good_matches[i].trainIdx;
		pSrcPts[i].x = kpPanoPatch[si].pt.x;
		pSrcPts[i].y = kpPanoPatch[si].pt.y;
		pDstPts[i].x = kpFrame[di].pt.x;
		pDstPts[i].y = kpFrame[di].pt.y;
	}

	//double ap[9];
	//AffineRansac(pDstPts, pSrcPts, nmatch, ap);
	Mat hm = HomographyRansac(pSrcPts, pDstPts, nmatch, hp);

	//warp image
	//Mat warp = mPanoPatch.clone();
	//warp = create(Size(mPanoPatch.cols, mPanoPatch.rows), mPanoPatch.type(), CV_RGB(0, 0, 0));
	

	/*
	for (int j = 0; j < ht; j++)
	{
		uchar* sp = warp.ptr<uchar>(j);
		for (int i = 0; i < wd; i++)
		{
			//int dx = AFFINE_X(i, j, ap);
			//int dy = AFFINE_Y(i, j, ap);

			//int dx = SIMILARITY_X(i, j, ap);
			//int dy = SIMILARITY_Y(i, j, ap);

			double dx = HOMOGRAPHY_X(i, j, ap);
			double dy = HOMOGRAPHY_Y(i, j, ap);

			int ix = int(dx);
			int iy = int(dy);

			double deltax = dx - ix;
			double deltay = dy - iy;
			
			if ( ix < 0 || dx >= (resizeFrame.cols-1) || dy < 0 || dy >= (resizeFrame.rows-1) )
				continue;

			//warpPerspective

			uchar* dp = resizeFrame.ptr<uchar>(iy);

			sp[i * 3]   = dp[ix * 3];
			sp[i * 3+1] = dp[ix * 3+1];
			sp[i * 3+2] = dp[ix * 3+2];
		}
	}
	*/
	
	/*Mat warpimage(mPanoPatch.rows, mPanoPatch.cols, mPanoPatch.type(), CV_RGB(0, 0, 0));
	int ht = warpimage.rows;
	int wd = warpimage.cols;
	warpPerspective(resizeFrame, warpimage, hm, warpimage.size(), CV_INTER_NN);
    */

	//warp.create(warpimage.size(), warpimage.type());
	//warpimage.copyTo(warp);
	//warp = warpimage.clone();

	//imwrite("c:\\temp\\warp-src.jpg", mPanoPatch);
	//imwrite("c:\\temp\\warp-dst.jpg", warp);
	//namedWindow("warp-src");
	//namedWindow("warp-dst");
	//imshow("warp-src", mPanoPatch);
	//imshow("warp-dst", warp);	
	//change detection
	//ChangeDetection(mPanoPatch, warp);
	//return warp;

	//return hm;
}

