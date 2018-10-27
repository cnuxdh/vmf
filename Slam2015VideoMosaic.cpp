// Slam2015VideoMosaic.cpp : 定义控制台应用程序的入口点。
//

//#include "stdafx.h"
#ifdef WIN32
#include"windows.h"
#include "atlimage.h"
#endif
#include<iostream>
#include<algorithm>
#include<fstream>
#include<chrono>

#include<direct.h>

#include<opencv2/core/core.hpp>

#include<System.h>
#include"Frame.h"

#include"UavInfo.h"
using namespace std;
using namespace cv;
using namespace ORB_SLAM2;

void readsubtxt(char* file, std::vector<GpsInfo> & m_gpsList) {
	m_gpsList.clear();
	FILE* fp = fopen(file, "r");
	char str[100];
	if (fp)
	{
		char charlist[50][50] = { "" };/*指定分隔后子字符串存储的位置，这里定义二维字符串数组*/
		while (!feof(fp))
		{
			memset(str, 0, 100);
			if (fgets(str, sizeof(str), fp))
			{
				GpsInfo m_gps;
				sscanf(str, "%lf %lf %f %f", &m_gps.lat, &m_gps.lon, &m_gps.height, &m_gps.heading);
				m_gpsList.push_back(m_gps);

			}
		}
	}
}
inline int ShowMat(cv::Mat img, HWND hWndDisplay)
{
	if (img.channels()<3)
	{
		return -1;
	}
	try
	{
		RECT rect;
		GetClientRect(hWndDisplay, &rect);
		cv::Mat imgShow(abs(rect.top - rect.bottom), abs(rect.right - rect.left), CV_8UC3);
		resize(img, imgShow, imgShow.size());

		ATL::CImage CI;
		int w = imgShow.cols;//宽
		int h = imgShow.rows;//高
		int channels = imgShow.channels();//通道数

		CI.Create(w, h, 8 * channels);
		uchar *pS;
		uchar *pImg = (uchar *)CI.GetBits();//得到CImage数据区地址
		int step = CI.GetPitch();
		for (int i = 0; i<h; i++)
		{
			pS = imgShow.ptr<uchar>(i);
			for (int j = 0; j<w; j++)
			{
				for (int k = 0; k<3; k++)
					*(pImg + i*step + j * 3 + k) = pS[j * 3 + k];
				//注意到这里的step不用乘以3
			}
		}

		HDC dc;
		dc = GetDC(hWndDisplay);
		CI.Draw(dc, 0, 0);
		ReleaseDC(hWndDisplay, dc);
		CI.Destroy();
	}
	catch (Exception* e)
	{
		int a = 0;
	}
	return 0;
}
#define MAX_PATH 256

//loading video file, added by xdh, 2018.7.28
//C:\Work\Data\turkin-video\rgb-video\sub.mp4  ORBvoc.bin  cfg.yaml 0 0
//C:\Work\Data\turkin-video\00060\test1.mp4 ORBvoc.bin  cfg.yaml 0 0
int  main(int argc, char **argv)
{
	if (argc != 6)
	{
		cerr << endl << "Usage: ./mono path_to_video binfile yamlfile realhandle mosaichandle" << endl;
		return -1;
	}
	std::vector<GpsInfo>  m_gpsList;
	readsubtxt("sub.txt", m_gpsList);
	//get the path containing the exe 
	char path[MAX_PATH];
	_getcwd(path, MAX_PATH);

	//valcabulary library file
	char vocfile[MAX_PATH];
	strcpy(vocfile, argv[2]);
	/*strcat(vocfile, "\\ORBvoc.bin");*/
	printf("%s \n", vocfile);

	//setting file
	char settingfile[MAX_PATH];
	strcpy(settingfile, argv[3]);
	//strcat(settingfile, "\\cfg.yaml");
	printf("%s \n", settingfile);
	//input handler
	int64 wndhandle = atoi(argv[4]);
	printf("%d \n", wndhandle);
	ORB_SLAM2::System SLAM(vocfile, settingfile, wndhandle, ORB_SLAM2::System::MONOCULAR, false);

	//input video file
	printf("video file: %s \n", argv[1]);
	int64 capturehandle = atoi(argv[5]);
	printf("capture hwnd:%d \n", capturehandle);
	HWND hwnddisplay = NULL;
	if (capturehandle > 0) {
		hwnddisplay = (HWND)capturehandle;
	}
	VideoCapture cap;
	cap.open(argv[1]);

	int nFrameRatio = 25;

	if (!cap.isOpened()) {
		printf("failed to open video file ! \n");
		return -1;
	}
	int frameCount = 0;
	frameCount = cap.get(CV_CAP_PROP_FRAME_COUNT);

	//namedWindow("video");

	//seconds for neigboring frames
	double timefortrack = 1.0 / double(nFrameRatio);
	double tframe = 0;

	int nKeyFrame = 0;

	Mat frame;
	GpsInfo info;
    int frameIndex = 0;
	while (1) {

        frameIndex++;

        if (frameIndex % 2 == 0)
            continue;

		tframe += timefortrack;
		cap >> frame;
		if (frame.empty())
			break;

		//resize the image
		Mat simg;
		//double ratio = 640.0 / (double)(frame.cols);
        int sht = SLAM.GetFixedHt();// frame.rows*ratio;
        int swd = SLAM.GetFixedWd();// frame.cols*ratio;
		resize(frame, simg, Size(swd, sht), 0, 0, INTER_LINEAR);

		//printf("%d %d \n", sht, swd);
		if (hwnddisplay) {
			ShowMat(simg, hwnddisplay);
			cv::waitKey(30);
		}
		else
		{
			cv::imshow("video", simg);
			cv::waitKey(30);
		}

		// Pass the image to the SLAM system
		//printf("%lf \n", tframe);
		if (!m_gpsList.empty()) {
			
            int frame_num = cap.get(CV_CAP_PROP_POS_FRAMES);
            //printf("frame num: %d \n", frame_num);

			int count = m_gpsList.size();
			info = m_gpsList[int(1.0*frame_num / frameCount*count)];

            //printf("%lf %lf \n", info.lat, info.lon);

		}

		SLAM.TrackMonocular(simg, tframe, info);

	}

    SLAM.FinalMosaic();
    //SLAM.SaveMapPoints("c:\\temp\\map.txt");

	cap.release();


	return 0;
}


