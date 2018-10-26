#include"absori.h"
#include"Matrix.h"

#include<vector>
using namespace std;


#include<opencv2/opencv.hpp>

using namespace cv;

#include"gdal_priv.h"


//similarity transform
#define SIMILARITY_X(x,y,p) p[0]*x-p[1]*y+p[2]
#define SIMILARITY_Y(x,y,p) p[1]*x+p[0]*y+p[3]



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


//
int RansanAbsOri(vector<POINT3D> srcPts, vector<POINT3D> dstPts,
	stAbsPOS& absPosParams)
{

	//RANSAC absolute orientation
	int ransac_rounds = 1000;
	int minWrongNumber = 10000000;
	int DistanceThrshold = 16;
	for (int round = 0; round < ransac_rounds; round++)
	{
		int indices[3];
		ChooseRandIndex(srcPts.size(), 3, indices);

		vector<POINT3D> freePts;
		vector<POINT3D> grdPts;
		freePts.resize(3);
		grdPts.resize(3);
		for (int i = 0; i<3; i++)
		{
			int ri = indices[i];

			freePts[i].x = srcPts[ri].x;
			freePts[i].y = srcPts[ri].y;
			freePts[i].z = srcPts[ri].z;

			grdPts[i].x = dstPts[ri].x;
			grdPts[i].y = dstPts[ri].y;
			grdPts[i].z = dstPts[ri].z;
		}

		stAbsPOS absPara;
		AbsOriP3(freePts, grdPts, absPara);

		//calculate the errors
		int nWrongNumber = 0;
		for (int i = 0; i<srcPts.size(); i++)
		{
			double fP[3];
			fP[0] = srcPts[i].x;
			fP[1] = srcPts[i].y;
			fP[2] = srcPts[i].z;
            //printf("free camera center: %lf %lf %lf \n", fP[0], fP[1], fP[2]);

			double tp[3];
			matrix_product(3, 3, 3, 1, absPara.R, fP, tp);
			for (int k = 0; k<3; k++)
				fP[k] = absPara.scale*tp[k] + absPara.T[k];

            //printf("transformed camera center: %lf %lf %lf \n", fP[0], fP[1], fP[2]);


			double gP[3];
			gP[0] = dstPts[i].x;
			gP[1] = dstPts[i].y;
			gP[2] = dstPts[i].z;

            //printf("gps: %lf %lf %lf \n", gP[0], gP[1], gP[2]);


			double len = 0;
			for (int k = 0; k<3; k++)
				len += (fP[k] - gP[k])*(fP[k] - gP[k]);
			len = sqrt(len);
			if (len>DistanceThrshold)
				nWrongNumber++;
		}

		if (minWrongNumber>nWrongNumber)
		{
			minWrongNumber = nWrongNumber;
			absPosParams = absPara;
		}
	}
	printf("minimux wrong number: %d \n", minWrongNumber);

    printf("abs scale: %lf \n", absPosParams.scale);
    printf("abs translation ... ");
    for (int k = 0; k < 3; k++)
    {
        printf("%lf ", absPosParams.T[k]);
    }
    printf("\n");

	return 0;
}



/*
using only three points to calculate the absolute orientation,
written by Donghai,Xie, 2015.10.10
the function can be invoked in RANSAC
inputs:
freePts,grdPts: free and ground control points
*/
int AbsOriP3(vector<POINT3D> freePts, vector<POINT3D> grdPts,
	stAbsPOS& absPosParams)
{
	if (freePts.size() != grdPts.size())
	{
		return -1;
	}
	if (freePts.size() != 3)
	{
		printf("size must be 3 ! \n");
		return -1;
	}

	//calculate the centroid
	double fcx = 0, fcy = 0, fcz = 0;
	double gcx = 0, gcy = 0, gcz = 0;
	for (int j = 0; j<freePts.size(); j++)
	{
		double xj = freePts[j].x;
		double yj = freePts[j].y;
		double zj = freePts[j].z;
		double gxj = grdPts[j].x;
		double gyj = grdPts[j].y;
		double gzj = grdPts[j].z;

		fcx += xj;		fcy += yj;		fcz += zj;
		gcx += gxj;		gcy += gyj;		gcz += gzj;
	}
	POINT3D freeCenterPt, grdCenterPt;
	freeCenterPt.x = fcx / double(freePts.size());
	freeCenterPt.y = fcy / double(freePts.size());
	freeCenterPt.z = fcz / double(freePts.size());
	grdCenterPt.x = gcx / double(freePts.size());
	grdCenterPt.y = gcy / double(freePts.size());
	grdCenterPt.z = gcz / double(freePts.size());

	//
	for (int j = 0; j<freePts.size(); j++)
	{
		freePts[j].x = freePts[j].x - freeCenterPt.x;
		freePts[j].y = freePts[j].y - freeCenterPt.y;
		freePts[j].z = freePts[j].z - freeCenterPt.z;
		grdPts[j].x = grdPts[j].x - grdCenterPt.x;
		grdPts[j].y = grdPts[j].y - grdCenterPt.y;
		grdPts[j].z = grdPts[j].z - grdCenterPt.z;
	}

	//calculate the scale
	double sumScale = 0;
	int    ns = 0;
	for (int j = 0; j<freePts.size(); j++)
	{
		double fx = freePts[j].x;
		double fy = freePts[j].y;
		double fz = freePts[j].z;
		double gx = grdPts[j].x;
		double gy = grdPts[j].y;
		double gz = grdPts[j].z;
		double gLen = sqrt(gx*gx + gy * gy + gz * gz);
		double fLen = sqrt(fx*fx + fy * fy + fz * fz);

		if (fLen != 0)
		{
			sumScale += gLen / fLen;
			ns++;
		}
	}
	absPosParams.scale = sumScale / (double)(ns);

	//calculate the R
	double sumC[9];
	memset(sumC, 0, sizeof(double) * 9);
	for (int j = 0; j<freePts.size(); j++)
	{
		double gp[3]; //ground point
		double fp[3]; //free point

		fp[0] = freePts[j].x;
		fp[1] = freePts[j].y;
		fp[2] = freePts[j].z;
		gp[0] = grdPts[j].x;
		gp[1] = grdPts[j].y;
		gp[2] = grdPts[j].z;

		double mc[9];
		matrix_product(3, 1, 1, 3, gp, fp, mc);
		for (int k = 0; k<9; k++)
			sumC[k] += mc[k];
	}
	//svd
	Mat mc, mw, mu, mvt;
	mc.create(3, 3, CV_32FC1);
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			mc.at<float>(j, i) = sumC[j * 3 + i];
		}
	}

	SVD::compute(mc, mw, mu, mvt);
	
	double U[9];
	double VT[9];
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			U[j * 3 + i] = mu.at<float>(j, i);
			VT[j * 3 + i] = mvt.at<float>(j, i);
		}
	}

	matrix_product(3, 3, 3, 3, U, VT, absPosParams.R);

	//calculate the T
	double gpc[3]; //center of the ground points
	double fpc[3]; //center of the free points
	gpc[0] = grdCenterPt.x;	    gpc[1] = grdCenterPt.y;	    gpc[2] = grdCenterPt.z;
	fpc[0] = freeCenterPt.x;	fpc[1] = freeCenterPt.y;	fpc[2] = freeCenterPt.z;
	double tc[3];
	matrix_product(3, 3, 3, 1, absPosParams.R, fpc, tc);
	for (int i = 0; i<3; i++)
		absPosParams.T[i] = gpc[i] - absPosParams.scale*tc[i];

	return 0;
}





/*
x' = ax - by + c
y' = bx + ay + d;
*/
void SimilarityTransform(vector<POINT3D>& oriPt, vector<POINT3D>& desPt, double* p)
{
    int num = oriPt.size();

    double* A = NULL;
    double* At = NULL;
    double* AtA = NULL;
    double* L = NULL;
    double* AtL = NULL;
    double* res = NULL;
    int  paraNum = 4;
    int  i, j;

    A = (double*)malloc(num * 2 * paraNum * sizeof(double));
    At = (double*)malloc(num * 2 * paraNum * sizeof(double));
    L = (double*)malloc(num * 2 * sizeof(double));
    AtA = (double*)malloc(paraNum*paraNum * sizeof(double));
    AtL = (double*)malloc(paraNum * sizeof(double));
    res = (double*)malloc(paraNum * sizeof(double));

    for (i = 0; i<num; i++)
    {
        A[2 * i*paraNum]     = oriPt[i].x;
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

    transpose(A, At, num * 2, paraNum);
    mult(At, A, AtA, paraNum, num * 2, paraNum);
    mult(At, L, AtL, paraNum, num * 2, 1);
    invers_matrix(AtA, paraNum);
    mult(AtA, AtL, res, paraNum, paraNum, 1);

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


int  SimRansac(vector<POINT3D>& srcPt, vector<POINT3D>& desPt, double *ap)
{
    int npt = srcPt.size();
    int ransac_rounds = 100;
    int maxInliers = 0;
    int DistanceThrshold = 16;

#define RANSAC_PT_NUM 4 

    for (int round = 0; round < ransac_rounds; round++)
    {
        int indices[RANSAC_PT_NUM];
        ChooseRandIndex(npt, RANSAC_PT_NUM, indices);

        vector<POINT3D> dpts(RANSAC_PT_NUM);
        vector<POINT3D> spts(RANSAC_PT_NUM);
        for (int k = 0; k < RANSAC_PT_NUM; k++)
        {
            int randi = indices[k];
            dpts[k].x = desPt[randi].x;
            dpts[k].y = desPt[randi].y;
            spts[k].x = srcPt[randi].x;
            spts[k].y = srcPt[randi].y;
        }

        double p[6];
        SimilarityTransform(spts, dpts, p);
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


int  GdalWriteImageColor(char* savePath, unsigned char* r, unsigned char* g, unsigned char* b,
    int ht, int wd, stGeoInfo geoInfo)
{
    GDALAllRegister();
    GDALDriver* poDriver = NULL;
    GDALDataset *poDataset = NULL;   //GDALÊý¾Ý¼¯
    GDALRasterBand *poBand = NULL;
    char **papszOptions = NULL;
    //GDALRasterBand *poBand;
    poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
    if (poDriver == NULL)
    {
        printf("Gdal Open file failed ! \n");
        return 0;
    }

    poDataset = poDriver->Create(savePath, wd, ht, 3, GDT_Byte, papszOptions);

    double  geoTransform[6];
    memset(geoTransform, 0, sizeof(double) * 6);
    geoTransform[0] = geoInfo.left;
    geoTransform[3] = geoInfo.top;
    geoTransform[1] = geoInfo.dx;
    geoTransform[5] = geoInfo.dy;
    poDataset->SetGeoTransform(geoTransform);
    poDataset->SetProjection(geoInfo.projectRef);

    poBand = poDataset->GetRasterBand(1);
    poBand->RasterIO(GF_Write, 0, 0, wd, ht, r, wd, ht, GDT_Byte, 0, 0);
    poBand = poDataset->GetRasterBand(2);
    poBand->RasterIO(GF_Write, 0, 0, wd, ht, g, wd, ht, GDT_Byte, 0, 0);
    poBand = poDataset->GetRasterBand(3);
    poBand->RasterIO(GF_Write, 0, 0, wd, ht, b, wd, ht, GDT_Byte, 0, 0);

	GDALSetRasterNoDataValue(poDataset->GetRasterBand(1), 0);
	GDALSetRasterNoDataValue(poDataset->GetRasterBand(2), 0);
	GDALSetRasterNoDataValue(poDataset->GetRasterBand(3), 0);
    //close
    GDALClose((GDALDatasetH)poDataset);
    return 1;
}