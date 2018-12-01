/******************************************************************************

  This file is part of Map2DFusion.

  Copyright 2016 (c)  Yong Zhao <zd5945@126.com> http://www.zhaoyong.adv-ci.com

  ----------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/
#ifdef WIN32
#include"windows.h"
#include "atlimage.h"
#endif

#include "MultiBandMap2D.h"
#include "types.h"
#include "SPtr.h"
#include "svar.h"

#include"System.h"
#include"Converter.h"
#include"FitObject.h"
#include"Matrix.h"
#include"LatLong-UTMconversion.h"

#include"gdal_priv.h"
#include"ogr_spatialref.h"


//#include <gui/gl/glHelper.h>
//#include <GL/gl.h>
//#include <gui/gl/SignalHandle.h>
//#include <base/Svar/Svar.h>
//#include <base/time/Global_Timer.h>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/stitching/stitcher.hpp>

#include"smooth.h"


//#define HAS_GOOGLEMAP
#ifdef HAS_GOOGLEMAP
#include <hardware/Gps/utils_GPS.h>
#include <base/Svar/Scommand.h>
#endif

using namespace std;

/**

  __________max
  |    |    |
  |____|____|
  |    |    |
  |____|____|
 min
 */



bool Map2DPrepare::prepare(const pi::SE3d& plane, const PinHoleParameters& camera,
	const std::deque<std::pair<cv::Mat, pi::SE3d> >& frames)
{
	if (frames.size() == 0 || camera.w <= 0 || camera.h <= 0 || camera.fx == 0 || camera.fy == 0)
	{
		cerr << "Map2D::prepare:Not valid prepare!\n";
		return false;
	}
	_camera = camera; _fxinv = 1. / camera.fx; _fyinv = 1. / camera.fy;
	_plane  = plane;
	_frames = frames;
	
    //deleted by xdh, 2018.9.13
    //   for (std::deque<std::pair<cv::Mat, pi::SE3d> >::iterator it = _frames.begin(); it != _frames.end(); it++)
	//{
	//	pi::SE3d& pose = it->second;
	//	pose = plane.inverse()*pose;//plane coordinate
	//}
	
    return true;
}

MultiBandMap2DCPU::MultiBandMap2DCPUEle::~MultiBandMap2DCPUEle()
{

    //if(texName) pi::gl::Signal_Handle::instance().delete_texture(texName);
}

bool MultiBandMap2DCPU::MultiBandMap2DCPUEle::normalizeUsingWeightMap(const cv::Mat& weight, cv::Mat& src)
{
    if(!(src.type()==CV_32FC3&&weight.type()==CV_32FC1)) return false;
    pi::Point3f* srcP=(pi::Point3f*)src.data;
    float*    weightP=(float*)weight.data;
    for(float* Pend=weightP+weight.cols*weight.rows;weightP!=Pend;weightP++,srcP++)
        *srcP=(*srcP)/(float)(*weightP+1e-5);
    return true;
}

bool MultiBandMap2DCPU::MultiBandMap2DCPUEle::mulWeightMap(const cv::Mat& weight, cv::Mat& src)
{
    if(!(src.type()==CV_32FC3&&weight.type()==CV_32FC1)) return false;
    pi::Point3f* srcP=(pi::Point3f*)src.data;
    float*    weightP=(float*)weight.data;
    for(float* Pend=weightP+weight.cols*weight.rows;weightP!=Pend;weightP++,srcP++)
        *srcP=(*srcP)*(*weightP);
    return true;
}

cv::Mat MultiBandMap2DCPU::MultiBandMap2DCPUEle::blend(const std::vector<SPtr<MultiBandMap2DCPUEle> >& neighbors)
{
    if(!pyr_laplace.size()) return cv::Mat();
    if(neighbors.size()==9)
    {
        //blend with neighbors, this obtains better visualization
        int flag=0;
        for(int i=0;i<neighbors.size();i++)
        {
            flag<<=1;
            if(neighbors[i].get()&&neighbors[i]->pyr_laplace.size())
                flag|=1;
        }
        switch (flag) {
        case 0X01FF:
        {
            vector<cv::Mat> pyr_laplaceClone(pyr_laplace.size());
            for(int i=0;i<pyr_laplace.size();i++)
            {
                int borderSize=1<<(pyr_laplace.size()-i-1);
                int srcrows=pyr_laplace[i].rows;
                int dstrows=srcrows+(borderSize<<1);
                pyr_laplaceClone[i]=cv::Mat(dstrows,dstrows,pyr_laplace[i].type());

                for(int y=0;y<3;y++)
                    for(int x=0;x<3;x++)
                {
                    const SPtr<MultiBandMap2DCPUEle>& ele=neighbors[3*y+x];
                    //pi::ReadMutex lock(ele->mutexData);
                    if(ele->pyr_laplace[i].empty())
                        continue;
                    cv::Rect      src,dst;
                    src.width =dst.width =(x==1)?srcrows:borderSize;
                    src.height=dst.height=(y==1)?srcrows:borderSize;
                    src.x=(x==0)?(srcrows-borderSize):0;
                    src.y=(y==0)?(srcrows-borderSize):0;
                    dst.x=(x==0)?0:((x==1)?borderSize:(dstrows-borderSize));
                    dst.y=(y==0)?0:((y==1)?borderSize:(dstrows-borderSize));
                    ele->pyr_laplace[i](src).copyTo(pyr_laplaceClone[i](dst));
                }
            }

            cv::detail::restoreImageFromLaplacePyr(pyr_laplaceClone);

            {
                cv::Mat result;
                int borderSize=1<<(pyr_laplace.size()-1);
                pyr_laplaceClone[0](cv::Rect(borderSize,borderSize,ELE_PIXELS,ELE_PIXELS)).copyTo(result);
                return  result.setTo(cv::Scalar::all(0),weights[0]==0);
            }
        }
            break;
        default:
            break;
        }
    }

    {
        //blend by self
        vector<cv::Mat> pyr_laplaceClone(pyr_laplace.size());
        for(int i=0;i<pyr_laplace.size();i++)
        {
            pyr_laplaceClone[i]=pyr_laplace[i].clone();
        }

        cv::detail::restoreImageFromLaplacePyr(pyr_laplaceClone);

        return  pyr_laplaceClone[0].setTo(cv::Scalar::all(0),weights[0]==0);
    }
}

// this is a bad idea, just for test
bool MultiBandMap2DCPU::MultiBandMap2DCPUEle::updateTexture(const std::vector<SPtr<MultiBandMap2DCPUEle> >& neighbors)
{
    cv::Mat tmp=blend(neighbors);
    unsigned int type=0;
    
	/*
	if(tmp.empty()) return false;
	else if(tmp.type()==CV_16SC3)
    {
        tmp.convertTo(tmp,CV_8UC3);
        type=GL_UNSIGNED_BYTE;
    }
    else if(tmp.type()==CV_32FC3)
        type=GL_FLOAT;
    if(!type) return false;

    {
        if(texName==0)// create texture
        {
            glGenTextures(1, &texName);
            glBindTexture(GL_TEXTURE_2D,texName);
            glTexImage2D(GL_TEXTURE_2D, 0,
                         GL_RGB,tmp.cols,tmp.rows, 0,
						 GL_RGB, type, tmp.data);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
        }
        else
        {
            glBindTexture(GL_TEXTURE_2D,texName);
            glTexImage2D(GL_TEXTURE_2D, 0,
                         GL_RGB,tmp.cols,tmp.rows, 0,
						 GL_RGB, type, tmp.data);
        }
    }

    SvarWithType<cv::Mat>::instance()["LastTexMat"]=tmp;
    SvarWithType<cv::Mat>::instance()["LastTexMatWeight"]=weights[0].clone();

    Ischanged=false;
	*/

    return true;
}

MultiBandMap2DCPU::MultiBandMap2DCPUData::MultiBandMap2DCPUData(double eleSize_,double lengthPixel_,pi::Point3d max_,pi::Point3d min_,
             int w_,int h_,const std::vector<SPtr<MultiBandMap2DCPUEle> >& d_)
    :_eleSize(eleSize_),_eleSizeInv(1./eleSize_),
      _lengthPixel(lengthPixel_),_lengthPixelInv(1./lengthPixel_),
      _min(min_),_max(max_),_w(w_),_h(h_),_data(d_)
{
    //_gpsOrigin=svar.get_var("GPS.Origin",_gpsOrigin);
}

bool MultiBandMap2DCPU::MultiBandMap2DCPUData::prepare(SPtr<MultiBandMap2DCPUPrepare> prepared)
{
    if(_w||_h) return false;//already prepared
    {
        _max=pi::Point3d(-1e10,-1e10,-1e10);
        _min=-_max;
        for(std::deque<std::pair<cv::Mat,pi::SE3d> >::iterator it=prepared->_frames.begin();
            it!=prepared->_frames.end();it++)
        {
            pi::SE3d& pose=it->second;
            pi::Point3d& t=pose.get_translation();
            _max.x=t.x>_max.x?t.x:_max.x;
            _max.y=t.y>_max.y?t.y:_max.y;
            _max.z=t.z>_max.z?t.z:_max.z;
            _min.x=t.x<_min.x?t.x:_min.x;
            _min.y=t.y<_min.y?t.y:_min.y;
            _min.z=t.z<_min.z?t.z:_min.z;
        }
        if(_min.z*_max.z<=0) return false;
        cout<<"Box:Min:"<<_min<<",Max:"<<_max<<endl;
    }
    //estimate w,h and bonding box
    {
        double maxh;
        if(_max.z>0) maxh=_max.z;
        else maxh=-_min.z;
        pi::Point3d line=prepared->UnProject(pi::Point2d(prepared->_camera.w,prepared->_camera.h))
                -prepared->UnProject(pi::Point2d(0,0));
        double radius=0.5*maxh*sqrt((line.x*line.x+line.y*line.y));
        
		_lengthPixel=svar.GetDouble("Map2D.Resolution",0);
        if(!_lengthPixel)
        {
            cout<<"Auto resolution from max height "<<maxh<<"m.\n";
            _lengthPixel=2*radius/sqrt(prepared->_camera.w*prepared->_camera.w
                                       +prepared->_camera.h*prepared->_camera.h);

            _lengthPixel/=svar.GetDouble("Map2D.Scale",1);
        }
        cout<<"Map2D.Resolution="<<_lengthPixel<<endl;
        _lengthPixelInv=1./_lengthPixel;
        _min=_min-pi::Point3d(radius,radius,0);
        _max=_max+pi::Point3d(radius,radius,0);
        pi::Point3d center=0.5*(_min+_max);
        _min=2*_min-center;_max=2*_max-center;
        _eleSize=ELE_PIXELS*_lengthPixel;
        _eleSizeInv=1./_eleSize;
        {
            _w=ceil((_max.x-_min.x)/_eleSize);
            _h=ceil((_max.y-_min.y)/_eleSize);
            _max.x=_min.x+_eleSize*_w;
            _max.y=_min.y+_eleSize*_h;
            _data.resize(_w*_h);

            printf("data size: %d %d \n", _w, _h);
        }
    }
    _gpsOrigin=svar.get_var("GPS.Origin",_gpsOrigin);
    return true;
}

MultiBandMap2DCPU::MultiBandMap2DCPU(Map* pMap, const string &strSettingPath, CUDPClient* messageSendClient, int64 wndDisplay)
	:mbFinishRequested(false), mbFinished(true), mbStopped(true), mbStopRequested(false), hWndDisplay(NULL), messageSendClient(NULL)
{
	_thread = false;
	alpha = 0;
	_highQualityShow = 1;
	_bandNum = 5;
	_bandNum = min(_bandNum, static_cast<int>(ceil(log(ELE_PIXELS) / log(2.0))));

	//reading settings
	cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
	mFx = fSettings["Camera.fx"];
	mFy = fSettings["Camera.fy"];
	mCx = fSettings["Camera.cx"];
	mCy = fSettings["Camera.cy"];
	mHt = fSettings["Camera.ht"];
	mWd = fSettings["Camera.wd"];
	printf("camra paras: %lf %lf %lf %lf %lf %lf \n", mFx, mFy, mCx, mCy, mHt, mWd);
	if (wndDisplay != 0){
		hWndDisplay = HWND(wndDisplay);
	}
	mpMap = pMap;

	mnCurrentKeyNumber = 0;
	this->messageSendClient = messageSendClient;
	//namedWindow("mosaic");

    //for plane fitting
    mPlanePos = pi::SE3d(0, 0, 0, 0, 0, 0, 1);

}

MultiBandMap2DCPU::MultiBandMap2DCPU(int64 wndDisplay,bool thread)
    :alpha(svar.GetInt("Map2D.Alpha",0)),
     _valid(false),_thread(thread),
     _bandNum(svar.GetInt("MultiBandMap2DCPU.BandNumber",5)),
	 _highQualityShow(svar.GetInt("MultiBandMap2DCPU.HighQualityShow", 1)), hWndDisplay(NULL)
{
    _bandNum=min(_bandNum, static_cast<int>(ceil(log(ELE_PIXELS) / log(2.0))));
	if (wndDisplay != 0){
		hWndDisplay = HWND(wndDisplay);
	}
}

bool MultiBandMap2DCPU::prepare(const pi::SE3d& plane,const PinHoleParameters& camera,
                const std::deque<std::pair<cv::Mat,pi::SE3d> >& frames)
{
	printf("plane: %lf %lf %lf \n", plane.get_translation().x,
		plane.get_translation().y, plane.get_translation().z);

    //insert frames
    SPtr<MultiBandMap2DCPUPrepare> p(new MultiBandMap2DCPUPrepare);
    SPtr<MultiBandMap2DCPUData>    d(new MultiBandMap2DCPUData);

    if(p->prepare(plane,camera,frames))
        if(d->prepare(p))
        {
            //pi::WriteMutex lock(mutex);
            prepared=p;
            data=d;
            weightImage.release();

            //if(_thread&&!isRunning())
            //    start();
            
			_valid=true;
            return true;
        }
    return false;
}

bool MultiBandMap2DCPU::feed(cv::Mat img,const pi::SE3d& pose)
{
    if(!_valid) return false;
    SPtr<MultiBandMap2DCPUPrepare> p;
    SPtr<MultiBandMap2DCPUData>    d;
    {
        //pi::ReadMutex lock(mutex);
        p=prepared;d=data;
    }
    std::pair<cv::Mat,pi::SE3d> frame(img,p->_plane.inverse()*pose);
    if(_thread)
    {
        //pi::WriteMutex lock(p->mutexFrames);
        p->_frames.push_back(frame);
        if(p->_frames.size()>20) p->_frames.pop_front();
        return true;
    }
    else
    {
        return renderFrame(frame);
    }
}

bool MultiBandMap2DCPU::renderFrame(const std::pair<cv::Mat,pi::SE3d>& frame)
{
    SPtr<MultiBandMap2DCPUPrepare> p;
    SPtr<MultiBandMap2DCPUData>    d;
    {
        //pi::ReadMutex lock(mutex);
        p=prepared;d=data;
    }

    if(frame.first.cols!=p->_camera.w||frame.first.rows!=p->_camera.h||frame.first.type()!=CV_8UC3)
    {
		printf("%d %d   %lf %lf \n ", frame.first.cols, frame.first.rows, 
			p->_camera.w, p->_camera.h);

        cerr<<"MultiBandMap2DCPU::renderFrame: frame.first.cols!=p->_camera.w||frame.first.rows!=p->_camera.h||frame.first.type()!=CV_8UC3\n";
        return false;
    }

    // 1. pose->pts
    std::vector<pi::Point2d>          imgPts;
    {
        imgPts.reserve(4);
        imgPts.push_back(pi::Point2d(0,0));
        imgPts.push_back(pi::Point2d(p->_camera.w,0));
        imgPts.push_back(pi::Point2d(0,p->_camera.h));
        imgPts.push_back(pi::Point2d(p->_camera.w,p->_camera.h));
    }
    vector<pi::Point2d> pts;
    pts.reserve(imgPts.size());
    pi::Point3d downLook(0,0,-1);
    if(frame.second.get_translation().z<0) downLook=pi::Point3d(0,0,1);
    for(int i=0;i<imgPts.size();i++)
    {
        pi::Point3d axis=frame.second.get_rotation()*p->UnProject(imgPts[i]);
        if(axis.dot(downLook)<0.4)
        {
            return false;
        }
        axis=frame.second.get_translation()
                -axis*(frame.second.get_translation().z/axis.z);

        pts.push_back(pi::Point2d(axis.x,axis.y));
        //printf("%lf %lf \n", axis.x, axis.y);
    }
    // 2. dest location?
    double xmin=pts[0].x;
    double xmax=xmin;
    double ymin=pts[0].y;
    double ymax=ymin;
    for(int i=1;i<pts.size();i++)
    {
        if(pts[i].x<xmin) xmin=pts[i].x;
        if(pts[i].y<ymin) ymin=pts[i].y;
        if(pts[i].x>xmax) xmax=pts[i].x;
        if(pts[i].y>ymax) ymax=pts[i].y;
    }
    if(xmin<d->min().x||xmax>d->max().x||ymin<d->min().y||ymax>d->max().y)
    {
        if(p!=prepared)//what if prepare called?
        {
            return false;
        }
        if(!spreadMap(xmin,ymin,xmax,ymax))
        {
            return false;
        }
        else
        {
            //pi::ReadMutex lock(mutex);
            if(p!=prepared)//what if prepare called?
            {
                return false;
            }
            d=data;//new data
        }
    }
    int xminInt=floor((xmin-d->min().x)*d->eleSizeInv());
    int yminInt=floor((ymin-d->min().y)*d->eleSizeInv());
    int xmaxInt= ceil((xmax-d->min().x)*d->eleSizeInv());
    int ymaxInt= ceil((ymax-d->min().y)*d->eleSizeInv());
    if(xminInt<0||yminInt<0||xmaxInt>d->w()||ymaxInt>d->h()||xminInt>=xmaxInt||yminInt>=ymaxInt)
    {
        cerr<<"MultiBandMap2DCPU::renderFrame:should never happen!\n";
        return false;
    }
    {
        xmin=d->min().x+d->eleSize()*xminInt;
        ymin=d->min().y+d->eleSize()*yminInt;
        xmax=d->min().x+d->eleSize()*xmaxInt;
        ymax=d->min().y+d->eleSize()*ymaxInt;
    }

    // 3.prepare weight and warp images
	//printf("prepare weight and warp images \n");
    cv::Mat weight_src;
    if(weightImage.empty()||weightImage.cols!=frame.first.cols||weightImage.rows!=frame.first.rows)
    {
        //pi::WriteMutex lock(mutex);
        int w=frame.first.cols;
        int h=frame.first.rows;
        weightImage.create(h,w,CV_32FC1);
        float *p=(float*)weightImage.data;
        float x_center=w/2;
        float y_center=h/2;
        float dis_max=sqrt(x_center*x_center+y_center*y_center);
        int weightType=svar.GetInt("Map2D.WeightType",0);
        for(int i=0;i<h;i++)
            for(int j=0;j<w;j++)
            {
                float dis=(i-y_center)*(i-y_center)+(j-x_center)*(j-x_center);
                dis=1-sqrt(dis)/dis_max;
                if(0==weightType)
                    *p=dis;
                else *p=dis*dis;
                if(*p<=1e-5) *p=1e-5;
                p++;
            }
        weight_src=weightImage.clone();
    }
    else
    {
        //pi::ReadMutex lock(mutex);
        weight_src=weightImage.clone();
    }

    std::vector<cv::Point2f>          imgPtsCV;
    {
        imgPtsCV.reserve(imgPts.size());
        for(int i=0;i<imgPts.size();i++)
            imgPtsCV.push_back(cv::Point2f(imgPts[i].x,imgPts[i].y));
    }
    std::vector<cv::Point2f> destPoints;
    destPoints.reserve(imgPtsCV.size());
    for(int i=0;i<imgPtsCV.size();i++)
    {
        destPoints.push_back(cv::Point2f((pts[i].x-xmin)*d->lengthPixelInv(),
                             (pts[i].y-ymin)*d->lengthPixelInv()));
        
        //destPoints.push_back(cv::Point2f((pts[i].x-xmin)*d->lengthPixelInv(),
        //                     (ymax-pts[i].y)*d->lengthPixelInv()));
         
    }

    cv::Mat transmtx = cv::getPerspectiveTransform(imgPtsCV, destPoints);

    cv::Mat img_src;
    if(svar.GetInt("MultiBandMap2DCPU.ForceFloat",0))
        frame.first.convertTo(img_src,CV_32FC3,1./255.);
    else
        frame.first.convertTo(img_src,CV_16SC3);

    cv::Mat weight_warped((ymaxInt-yminInt)*ELE_PIXELS,(xmaxInt-xminInt)*ELE_PIXELS,CV_32FC1);
    cv::Mat image_warped((ymaxInt-yminInt)*ELE_PIXELS,(xmaxInt-xminInt)*ELE_PIXELS,img_src.type());
    cv::warpPerspective(img_src, image_warped, transmtx, image_warped.size(),cv::INTER_LINEAR,cv::BORDER_REFLECT);
    cv::warpPerspective(weight_src, weight_warped, transmtx, weight_warped.size(),cv::INTER_NEAREST);

    /*if(svar.GetInt("ShowWarped",0))
    {
        cv::imshow("image_warped",image_warped);
        cv::imshow("weight_warped",weight_warped);
        if(svar.GetInt("SaveImageWarped"))
        {
            cout<<"Saving warped image.\n";
            cv::imwrite("image_warped.png",image_warped);
            cv::imwrite("weight_warped.png",weight_warped);
        }
        cv::waitKey(0);
    }*/

    // 4. blender dst to eles
	//printf("blending...\n");
	std::vector<cv::Mat> pyr_laplace;
    cv::detail::createLaplacePyr(image_warped, _bandNum, pyr_laplace);

    std::vector<cv::Mat> pyr_weights(_bandNum+1);
    pyr_weights[0]=weight_warped;
    for (int i = 0; i < _bandNum; ++i)
        cv::pyrDown(pyr_weights[i], pyr_weights[i + 1]);

    //pi::timer.enter("MultiBandMap2DCPU::Apply");
	//printf("Applying ... \n");
    std::vector<SPtr<MultiBandMap2DCPUEle> > dataCopy=d->data();
    for(int x=xminInt;x<xmaxInt;x++)
        for(int y=yminInt;y<ymaxInt;y++)
        {
            SPtr<MultiBandMap2DCPUEle> ele=dataCopy[y*d->w()+x];
            if(!ele.get())
            {
                ele=d->ele(y*d->w()+x);
            }
            {
                //pi::WriteMutex lock(ele->mutexData);
                if(!ele->pyr_laplace.size())
                {
                    ele->pyr_laplace.resize(_bandNum+1);
                    ele->weights.resize(_bandNum+1);
                }

                int width=ELE_PIXELS,height=ELE_PIXELS;

                for (int i = 0; i <= _bandNum; ++i)
                {
                    if(ele->pyr_laplace[i].empty())
                    {
                        //fresh
                        cv::Rect rect(width*(x-xminInt),height*(y-yminInt),width,height);
                        pyr_laplace[i](rect).copyTo(ele->pyr_laplace[i]);
                        pyr_weights[i](rect).copyTo(ele->weights[i]);
                    }
                    else
                    {
                        if(pyr_laplace[i].type()==CV_32FC3)
                        {
                            int org =(x-xminInt)*width+(y-yminInt)*height*pyr_laplace[i].cols;
                            int skip=pyr_laplace[i].cols-ele->pyr_laplace[i].cols;

                            pi::Point3f *srcL=((pi::Point3f*)pyr_laplace[i].data)+org;
                            float       *srcW=((float*)pyr_weights[i].data)+org;

                            pi::Point3f *dstL=(pi::Point3f*)ele->pyr_laplace[i].data;
                            float       *dstW=(float*)ele->weights[i].data;

                            for(int eleY=0;eleY<height;eleY++,srcL+=skip,srcW+=skip)
                                for(int eleX=0;eleX<width;eleX++,srcL++,dstL++,srcW++,dstW++)
                                {
                                    if((*srcW)>=(*dstW))
                                    {
                                        *dstL=(*srcL);
                                        *dstW=*srcW;
                                    }
                                }
                        }
                        else if(pyr_laplace[i].type()==CV_16SC3)
                        {
                            int org =(x-xminInt)*width+(y-yminInt)*height*pyr_laplace[i].cols;
                            int skip=pyr_laplace[i].cols-ele->pyr_laplace[i].cols;

                            pi::Point3_<short> *srcL=((pi::Point3_<short>*)pyr_laplace[i].data)+org;
                            float       *srcW=((float*)pyr_weights[i].data)+org;

                            pi::Point3_<short> *dstL=(pi::Point3_<short>*)ele->pyr_laplace[i].data;
                            float       *dstW=(float*)ele->weights[i].data;

                            for(int eleY=0;eleY<height;eleY++,srcL+=skip,srcW+=skip)
                                for(int eleX=0;eleX<width;eleX++,srcL++,dstL++,srcW++,dstW++)
                                {
                                    if((*srcW)>=(*dstW))
                                    {
                                        *dstL=(*srcL);
                                        *dstW=*srcW;
                                    }
                                }
                        }
                    }
                    width/=2;height/=2;
                }
                ele->Ischanged=true;
            }
        }
    
	//pi::timer.leave("MultiBandMap2DCPU::Apply");

    return true;
}


bool MultiBandMap2DCPU::spreadMap(double xmin,double ymin,double xmax,double ymax)
{
    //pi::timer.enter("MultiBandMap2DCPU::spreadMap");
	printf("spread map...\n");
    SPtr<MultiBandMap2DCPUData> d;
    {
        //pi::ReadMutex lock(mutex);
        d=data;
    }
    int xminInt=floor((xmin-d->min().x)*d->eleSizeInv());
    int yminInt=floor((ymin-d->min().y)*d->eleSizeInv());
    int xmaxInt= ceil((xmax-d->min().x)*d->eleSizeInv());
    int ymaxInt= ceil((ymax-d->min().y)*d->eleSizeInv());
    xminInt=min(xminInt,0); yminInt=min(yminInt,0);
    xmaxInt=max(xmaxInt,d->w()); ymaxInt=max(ymaxInt,d->h());
    int w=xmaxInt-xminInt;
    int h=ymaxInt-yminInt;
    pi::Point2d min,max;
    {
        min.x=d->min().x+d->eleSize()*xminInt;
        min.y=d->min().y+d->eleSize()*yminInt;
        max.x=min.x+w*d->eleSize();
        max.y=min.y+h*d->eleSize();
    }
    std::vector<SPtr<MultiBandMap2DCPUEle> > dataOld=d->data();
    std::vector<SPtr<MultiBandMap2DCPUEle> > dataCopy;
    dataCopy.resize(w*h);
    {
        for(int x=0,xend=d->w();x<xend;x++)
            for(int y=0,yend=d->h();y<yend;y++)
            {
                dataCopy[x-xminInt+(y-yminInt)*w]=dataOld[y*d->w()+x];
            }
    }
    //apply
    {
        //pi::WriteMutex lock(mutex);
        data=SPtr<MultiBandMap2DCPUData>(new MultiBandMap2DCPUData(d->eleSize(),d->lengthPixel(),
                                                 pi::Point3d(max.x,max.y,d->max().z),
                                                 pi::Point3d(min.x,min.y,d->min().z),
                                                 w,h,dataCopy));
    }
    //pi::timer.leave("MultiBandMap2DCPU::spreadMap");
    return true;
}

bool MultiBandMap2DCPU::getFrame(std::pair<cv::Mat,pi::SE3d>& frame)
{
    //pi::ReadMutex lock(mutex);
    //pi::ReadMutex lock1(prepared->mutexFrames);
    if(prepared->_frames.size())
    {
        frame=prepared->_frames.front();
        prepared->_frames.pop_front();
        return true;
    }
    else 
		return false;
}


// 向量归一化
int Normalize(double *invect, int num)
{
	double t = 0.0;
	int i = 0;
	for (i = 0; i<num; i++)
	{
		t += (invect[i] * invect[i]);
	}
	t = sqrt(t);
	for (i = 0; i<num; i++)
	{
		invect[i] = invect[i] / t;
	}

	return 1;
}


void FitPlane1(double* px, double* py, double* pz, int np, double* R)
{
	int i, j, k;

	// 拟合平面，旋转到与XZ平面平行
	double *errorinfo = new double[np];
	double *errorpoints = new double[np];
	double *tmpotx = new double[np];
	double *tmpoty = new double[np];
	double *tmpotz = new double[np];
	double para[100];
	int paranum = 100;
	int infonum = 0;
	FILE* hp = NULL;
	//char  tmppath[256] = "d:\\R.txt";

	FitObject(2, px, py, pz, np,
		para, paranum,
		errorinfo, errorpoints,
		tmpotx, tmpoty, tmpotz, infonum);

	/*
	//写文件
	hp = fopen(tmppath,"wt");
	for (i=0;i<4;i++)
	{
	fprintf(hp, "%d %lf %lf %lf\n", i, tmpotx[i],tmpoty[i],tmpotz[i]);
	}
	fprintf(hp, "%lf %lf %lf %lf\n",para[0],para[1],para[2],para[3]);
	*/

	int pairarray[100];
	pairarray[0] = 31;
	int objectinfo[3] = { 1, 0, 0 };
	Normalize(para, 3);
	double tx[1] = { para[0] };
	double ty[1] = { para[1] };
	double tz[1] = { para[2] };
	ComputeParameter(tx, ty, tz, pairarray, objectinfo, 1, para, paranum);
	//fprintf(hp, "%lf %lf %lf\n  %lf %lf %lf\n %lf %lf %lf\n", 
	//	para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
	//fclose(hp);

	//rotation matrix from original data to horizontal plane
	for (i = 0; i<9; i++)
		R[i] = para[i];

	delete[]errorinfo;
	delete[]errorpoints;
	delete[]tmpotx;
	delete[]tmpoty;
	delete[]tmpotz;
}

void PlaneFit(double* px, double* py, double* pz, int np, 
	double& mx, double& my, double& mz, 
	double& wx, double& wy, double& wz, double& w)
{
	printf("plane fitting... \n");
	printf("point number: %d \n", np);

	double R[9];
	FitPlane1(px, py, pz, np, R);

	for (int i = 0; i < np; i++)
	{
		double fpt[3];
		fpt[0] = px[i];
		fpt[1] = py[i];
		fpt[2] = pz[i];
		double res[3];
		mult(R, fpt, res, 3, 3, 1);
		px[i] = res[0];
		py[i] = res[1];
		pz[i] = res[2];
	}
	//WritePly("c:\\temp\\rotatemap.ply", npt - 100, px, py, pz, pr, pg, pb);
	mx = 0;
	my = 0;
	mz = 0;
	for (int i = 0; i < np; i++)
	{
		mx += px[i];
		my += py[i];
		mz += pz[i];
	}
	mx /= np;
	my /= np;
	mz /= np;

	printf("translation: %lf %lf %lf \n", mx, my, mz);

	invers_matrix(R, 3);

	//from rotation to quartion
	Eigen::Matrix<double, 3, 3> M;
	M << R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8];
	Eigen::Quaterniond q(M);
	wx = q.x();
	wy = q.y();
	wz = q.z();
	w = q.w();

	printf("quartion: %lf %lf %lf %lf \n", wx, wy, wz, w);
}

inline int MultiBandMap2DCPU::ShowMat(cv::Mat img, HWND hWndDisplay)
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

void FrameAbsOri(double* absR, double* absT, double absScale, 
	Mat R, Mat T, double &x, double &y, double &z,
	double &wx, double &wy, double &wz, double &w)
{
	double fR[9];
	double fT[3];
	int index = 0;
	for (int jj = 0; jj < 3; jj++) {
		fT[jj] = T.at<float>(jj);
		for (int ii = 0; ii < 3; ii++) {
			fR[index] = R.at<float>(jj, ii);
		}
	}

	double invAbsR[9];
	memcpy(invAbsR, absR, sizeof(double) * 9);
	invers_matrix(invAbsR, 3);

	//new R
	double newR[9];
	matrix_product(3, 3, 3, 3, fR, invAbsR, newR);
	
	//new T
	double newT[3];
	matrix_product(3, 3, 3, 1, newR, absT, newT);
	for (int k = 0; k<3; k++)
		newT[k] = absScale *fT[k] - newT[k];

	x = newT[0];
	y = newT[1];
	z = newT[2];
	for (int jj = 0; jj < 3; jj++) {
		for (int ii = 0; ii < 3; ii++) {
			R.at<float>(jj, ii) = newR[index];
		}
	}

	vector<float> q = Converter::toQuaternion(R);
	wx = q[0];
	wy = q[1];
	wz = q[2];
	w = q[3];

}

void MultiBandMap2DCPU::planeFit(bool bIsFit)
{
    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    int nKeyFrames = vpKFs.size();

    const vector<MapPoint*> &vpMPs = mpMap->GetAllMapPoints();
    int nmappt = vpMPs.size();    
    
    //remove the wrong map points 
    double* px = new double[nmappt];
    double* py = new double[nmappt];
    double* pz = new double[nmappt];
    int nvalidpt = 0;
    for (size_t i = 0, iend = nmappt; i < iend; i++)
    {
        if (vpMPs[i]->isBad())
            continue;
        cv::Mat pos = vpMPs[i]->GetWorldPos();
        px[nvalidpt] = pos.at<float>(0);
        py[nvalidpt] = pos.at<float>(1);
        pz[nvalidpt] = pos.at<float>(2);
        nvalidpt++;
    }
    
    int nSmoothPt = PointCloudFilter(px, py, pz, nvalidpt, 64, 0.6);
        
    double mx = 0;
    double my = 0;
    double mz = 0;
    for (size_t i = 0, iend = nSmoothPt; i < iend; i++)
    {
        mx += px[i];
        my += py[i];
        mz += pz[i];
    }
    mx /= double(nSmoothPt);
    my /= double(nSmoothPt);
    mz /= double(nSmoothPt);

    if (bIsFit)
    {
        double wx, wy, wz, w;
        PlaneFit(px, py, pz, nSmoothPt, mx, my, mz, wx, wy, wz, w);
        mPlanePos = pi::SE3d(mx, my, mz, wx, wy, wz, w);
    }
    else
    {
        mPlanePos = pi::SE3d(mx, my, mz, 0, 0, 0, 1);
    }

    /*
    FILE* fp = fopen("c:\\temp\\slampt.xyz", "w");
    for (int i = 0; i < nKeyFrames; i++)
    {
        KeyFrame* pKF = vpKFs[i];
        cv::Mat pos = pKF->GetCameraCenter();
        double xs = pos.at<float>(0);
        double ys = pos.at<float>(1);
        double zs = pos.at<float>(2);       
        fprintf(fp, "%lf %lf %lf %d %d %d \n", xs, ys, zs, 0, 255, 0);     
    }
    for (size_t i = 0, iend = nSmoothPt; i < iend; i++)
    {        
        fprintf(fp, "%lf %lf %lf %d %d %d \n", px[i], py[i], pz[i], 255, 0, 0);
    }
    fclose(fp);
    */

    delete[] px;
    delete[] py;
    delete[] pz;

    /*
    //////////////////  plane fitting ///////////////////////////////
    const vector<MapPoint*> &vpMPs = mpMap->GetAllMapPoints();
    int nmappt = vpMPs.size();
    double* px = new double[nmappt];
    double* py = new double[nmappt];
    double* pz = new double[nmappt];
    int npt = 0;

    for (size_t i = 0, iend = vpMPs.size(); i < iend; i++)
    {
        if (vpMPs[i]->isBad())
            continue;
        cv::Mat pos = vpMPs[i]->GetWorldPos();

        double fp[3];
        fp[0] = pos.at<float>(0);
        fp[1] = pos.at<float>(1);
        fp[2] = pos.at<float>(2);
        //printf("free map point: %lf %lf %lf \n", fp[0], fp[1], fp[2]);

        //fprintf(freeout, "%lf %lf %lf %d %d %d \n", fp[0], fp[1], fp[2], 1,0,0);

        //if (bIsAbsOri)
        //{
        //    //absolute orientation
        //    double tp[3];
        //    matrix_product(3, 3, 3, 1, mAbsParas.R, fp, tp);
        //    for (int k = 0; k < 3; k++)
        //        fp[k] = mAbsParas.scale*tp[k] + mAbsParas.T[k];
        //    
        //    printf("%lf %lf %lf -> %lf %lf %lf \n", tp[0], tp[1], tp[2],
        //        fp[0], fp[1], fp[2]);

        //    //output the points after abs orientation
        //    //FILE* fout = fopen("c:\\temp\\init.xyz", "a+");
        //    //fprintf(fout, "%lf %lf %lf %d %d %d \n", tp[0], tp[1], tp[2], 1, 1, 1);
        //    //fclose(fout);
        //}

        px[npt] = fp[0];
        py[npt] = fp[1];
        pz[npt] = fp[2];

        npt++;
    }
    //fclose(freeout);

    double mx, my, mz, wx, wy, wz, w;
    PlaneFit(px, py, pz, npt, mx, my, mz, wx, wy, wz, w);
    mPlanePos = pi::SE3d(mx, my, mz, wx, wy, wz, w);
    delete[] px;
    delete[] py;
    delete[] pz;

    //reorientation using the results of plane fitting
    FILE* initfp = fopen("c:\\temp\\init.xyz", "w");
    FILE* fitfp = fopen("c:\\temp\\fit.xyz", "w");
    for (int i = 0; i < nKeyFrames; i++)
    {
        KeyFrame* pKF = vpKFs[i];
        cv::Mat pos = pKF->GetCameraCenter();
        //cv::Mat pos = pKF->GetTranslation();

        double xs = pos.at<float>(0);
        double ys = pos.at<float>(1);
        double zs = pos.at<float>(2);
        pi::SE3d nse = mPlanePos.inverse()*pi::SE3d(xs, ys, zs, 0, 0, 0, 1);

        fprintf(initfp, "%lf %lf %lf %d %d %d \n", xs, ys, zs, 0, 255, 0);
        fprintf(fitfp, "%lf %lf %lf %d %d %d \n",
            nse.get_translation().x,
            nse.get_translation().y,
            nse.get_translation().z, 0, 255, 0);
    }

    for (size_t i = 0, iend = vpMPs.size(); i < iend; i++)
    {
        if (vpMPs[i]->isBad())
            continue;
        cv::Mat pos = vpMPs[i]->GetWorldPos();

        double x, y, z;
        x = pos.at<float>(0);
        y = pos.at<float>(1);
        z = pos.at<float>(2);
        fprintf(initfp, "%lf %lf %lf %d %d %d \n", x, y, z, 255, 0, 0);

        pi::SE3d nse = mPlanePos.inverse()*pi::SE3d(x, y, z, 0, 0, 0, 1);
        double nx = nse.get_translation().x;
        double ny = nse.get_translation().y;
        double nz = nse.get_translation().z;
        fprintf(fitfp, "%lf %lf %lf %d %d %d \n", nx, ny, nz, 255, 0, 0);
    }
    fclose(initfp);
    fclose(fitfp);
    */
}



void MultiBandMap2DCPU::run()
{
	mbFinished = false;
	mbStopped = false;

    bool bIsAbsOri = false;
	
	//cv::resizeWindow("mosaic", 1024, 640);
	int prepareMosaicCount = PREPARE_MOSAICFRAMES;
	while (1)
	{
		if (Stop())
		{
			while (isStopped())
			{
				usleep(3000);
			}
		}

		if (CheckFinish())
			break;

		//usleep(3000);
		//initialize the data 
		vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
		int nKeyFrames = vpKFs.size();

		if (nKeyFrames > mnCurrentKeyNumber){
			mnCurrentKeyNumber = nKeyFrames;
			printf("current key number: %d \n", mnCurrentKeyNumber);
		}
		else{
			//printf("key frames: %d \n", nKeyFrames);
			continue;
		}
		
        //************************ I: plane fitting *******************************
		if (nKeyFrames == prepareMosaicCount || (nKeyFrames>0 && nKeyFrames % 20 == 0))
        //if (nKeyFrames == prepareMosaicCount )
		{
            sort(vpKFs.begin(), vpKFs.end(), KeyFrame::lId);

            //camera parameters
            PinHoleParameters camParas;
            camParas.fx = mFx;
            camParas.fy = mFy;
            camParas.cx = mCx;
            camParas.cy = mCy;
            camParas.w = mWd;
            camParas.h = mHt;

            //1. plane fitting
            if (nKeyFrames == prepareMosaicCount)
            {
                planeFit(false);
            }
            else
            {
                planeFit(true);
            }

			//2. transform based on plane fitting
			printf("generating frames and rendering... \n");
			deque<std::pair<cv::Mat, pi::SE3d> > frames;	            
			for (int i = 0; i < nKeyFrames; i++)
			{
                KeyFrame* pKF = vpKFs[i];
				//printf("%d ", i);
				std::pair<cv::Mat, pi::SE3d> frame;
				frame.first = pKF->mFrameImage;

				double x, y, z, wx, wy, wz, w;
				cv::Mat R = pKF->GetRotation().t();
				cv::Mat t = pKF->GetCameraCenter();
                //cv::Mat t = pKF->GetTranslation();				
                {				
                    x = t.at<float>(0);
				    y = t.at<float>(1);
				    z = t.at<float>(2);
				    vector<float> q = Converter::toQuaternion(R);
				    wx = q[0];
				    wy = q[1];
				    wz = q[2];
				    w  = q[3];
                }
				printf("pos: %lf %lf %lf %lf %lf %lf %lf \n", x, y, z, wx, wy, wz, w);
				frame.second = mPlanePos.inverse()*pi::SE3d(x,y,z,wx,wy,wz,w);
                //frame.second = pi::SE3d(x, y, z, wx, wy, wz, w);
                frames.push_back(frame);           
			}
                       
            //prepare
            printf("\n prepare  ... \n");
            prepare(mPlanePos, camParas, frames);
            for (int i = 0; i < nKeyFrames; i++)
            {
                renderFrame(frames[i]);
            }            			
		}

        //****************** II: mosaic and geotransformation ******************
		if (nKeyFrames > prepareMosaicCount)
		{
			sort(vpKFs.begin(), vpKFs.end(), KeyFrame::lId);

            ///////////////// 3. rendering the last keyframe /////////////
			KeyFrame* pKF = vpKFs[nKeyFrames-1];
			Mat image = pKF->mFrameImage;
			double x, y, z, wx, wy, wz, w;
			cv::Mat R = pKF->GetRotation().t();
			vector<float> q = Converter::toQuaternion(R);
			cv::Mat t = pKF->GetCameraCenter();
            //cv::Mat t = pKF->GetTranslation();
            {
			    x = t.at<float>(0);
			    y = t.at<float>(1);
			    z = t.at<float>(2);
			    wx = q[0];
			    wy = q[1];
			    wz = q[2];
			    w = q[3];
            }            			
			printf("pos: %lf %lf %lf %lf %lf %lf %lf \n", x, y, z, wx, wy, wz, w);
			pi::SE3d pos = pi::SE3d(x, y, z, wx, wy, wz, w);
			printf("feed image ...\n");
			//feed(image, pos);			
            std::pair<cv::Mat, pi::SE3d> frame(image, prepared->_plane.inverse()*pos);
            //std::pair<cv::Mat, pi::SE3d> frame(image, pos);
            renderFrame(frame);
			
					
            ////////////////// 4. similarity transform /////////////////
            Mat result;
            int xoff, yoff;
            bool bIsOk = mosaic(result, xoff, yoff);
            if (!bIsOk)
                continue;

            //draw camera locations
            vector<POINT3D> srcPts;
            vector<POINT3D> dstPts;
            int ht = result.rows;
            int wd = result.cols;
            double minx = data->min().x + xoff * data->eleSize();
            double miny = data->min().y + yoff * data->eleSize();
            double slon = vpKFs[0]->mGpsInfo.lon;
            mZoneNumber = int((slon + 180) / 6) + 1;
            for (int k = 0; k < nKeyFrames; k++)
            {
                KeyFrame* pKF = vpKFs[k];

                double x, y, z, wx, wy, wz, w;
                cv::Mat R = pKF->GetRotation().t();
                vector<float> q = Converter::toQuaternion(R);
                cv::Mat t = pKF->GetCameraCenter();
                //cv::Mat t = pKF->GetTranslation();
                {
                    x = t.at<float>(0);
                    y = t.at<float>(1);
                    z = t.at<float>(2);
                    wx = q[0];
                    wy = q[1];
                    wz = q[2];
                    w = q[3];
                }
                pi::SE3d pos = pi::SE3d(x, y, z, wx, wy, wz, w);
                std::pair<cv::Mat, pi::SE3d> frame(image, prepared->_plane.inverse()*pos);

                double cx = frame.second.get_translation().x;
                double cy = frame.second.get_translation().y;
                double ix = (cx - minx)*data->lengthPixelInv();
                double iy = (cy - miny)*data->lengthPixelInv();
                printf("camera location: %lf %lf \n", ix, iy);
                //cv::circle(result, Point(ix, iy), 3, CV_RGB(255, 0, 0));
                POINT3D sp;
                sp.x = ix - wd * 0.5;
                sp.y = ht * 0.5 - iy;
                srcPts.push_back(sp);

                //convert from lon,lat to utm
                POINT3D dstpt;
                double lat = pKF->mGpsInfo.lat;
                double lon = pKF->mGpsInfo.lon;
                printf("lat,lon: %lf %lf \n", lat, lon);
                double gz = pKF->mGpsInfo.height;
                double gx, gy;
                LLtoUTM(23, lat, lon, gy, gx, mZoneNumber);
                POINT3D dp;
                dp.x = gx;
                dp.y = gy;
                dp.z = gz;
                double la, lo;
                UTMtoLL(23, gy, gx, mZoneNumber, la, lo);
                printf("%lf %lf  %lf %lf \n", gy, gx, la, lo);
                dstPts.push_back(dp);
            }
            double sp[4];
            //SimilarityTransform(srcPts, dstPts, sp);
            SimRansac(srcPts, dstPts, sp);
            
            for (int k = 0; k < srcPts.size(); k++)
            {
                double dx = srcPts[k].x*sp[0] - srcPts[k].y*sp[1] + sp[2];
                double dy = srcPts[k].x*sp[1] + srcPts[k].y*sp[0] + sp[3];
                printf("%lf %lf %lf %lf \n", dx, dy, dstPts[k].x, dstPts[k].y);
            }
            for (int k = 0; k < 4; k++)
            {
                mSP[k] = sp[k];
            }
                        
            /////// 5. determine the scope of transformed image //////
            std::vector<cv::Point2f> imgPtsCV;
            imgPtsCV.resize(4);
            imgPtsCV[0].x = -wd * 0.5;   imgPtsCV[0].y = ht * 0.5;
            imgPtsCV[1].x = wd * 0.5;    imgPtsCV[1].y = ht * 0.5;
            imgPtsCV[2].x = -wd * 0.5;   imgPtsCV[2].y = -ht * 0.5;
            imgPtsCV[3].x = wd * 0.5;    imgPtsCV[3].y = -ht * 0.5;
            double nsp[4];
            double dnorm = sqrt(sp[0] * sp[0] + sp[1] * sp[1]);
            nsp[0] = sp[0] / dnorm;
            nsp[1] = sp[1] / dnorm;
            nsp[2] = 0;
            nsp[3] = 0;
            int minx1 = 100000;
            int maxx1 = -100000;
            int miny1 = 100000;
            int maxy1 = -100000;
            double gx_min = 100000000;
            double gx_max = -100000000;
            double gy_min = 100000000;
            double gy_max = -100000000;
            for (int k = 0; k < 4; k++)
            {
                int dx = SIMILARITY_X(imgPtsCV[k].x, imgPtsCV[k].y, nsp);
                int dy = SIMILARITY_Y(imgPtsCV[k].x, imgPtsCV[k].y, nsp);
                if (minx1 > dx) minx1 = dx;
                if (maxx1 < dx) maxx1 = dx;
                if (miny1 > dy) miny1 = dy;
                if (maxy1 < dy) maxy1 = dy;

                double tx = SIMILARITY_X(imgPtsCV[k].x, imgPtsCV[k].y, mSP);
                double ty = SIMILARITY_Y(imgPtsCV[k].x, imgPtsCV[k].y, mSP);
                if (gx_min > tx) gx_min = tx;
                if (gx_max < tx) gx_max = tx;
                if (gy_min > ty) gy_min = ty;
                if (gy_max < ty) gy_max = ty;
            } 
            printf("%d %d %d %d \n", minx1, maxx1, miny1, maxy1);
            printf("%lf %lf %lf %lf \n", gx_min, gx_max, gy_min, gy_max);
            
            //////////// 6. resample the image based on transform //////////
            cv::Mat warp_res( maxy1 - miny1 + 1, maxx1 - minx1 + 1, CV_8UC3);
            int warpht = warp_res.rows;
            int warpwd = warp_res.cols;
            int row = 0;
            int col = 0;
            for (int j = miny1; j < maxy1; j++)
            {
                col = 0;
                uchar* pdst = warp_res.ptr<uchar>(warpht-row-1);
                for (int i = minx1; i < maxx1; i++)
                {
                    int sx = nsp[0]*i  + nsp[1]*j; // SIMILARITY_X(i, j, nisp);
                    int sy = -nsp[1]*i + nsp[0]*j; // SIMILARITY_Y(i, j, nisp);
                    sx = sx + wd*0.5;
                    sy = ht*0.5 - sy;
                    if (sx < 0)   sx = 0;
                    if (sx >= wd) sx = wd - 1;
                    if (sy < 0)   sy = 0;
                    if (sy >= ht) sy = ht - 1;
                    uchar* psrc = result.ptr<uchar>(sy); 
                    pdst[col * 3]     = psrc[sx * 3];
                    pdst[col * 3 + 1] = psrc[sx * 3 + 1];
                    pdst[col * 3 + 2] = psrc[sx * 3 + 2];
                    col++;
                }
                row++;
            }  

            /////////////// 7. add geoinfomation /////////////////
            double resolution = dnorm;
            int rht = warp_res.rows;
            int rwd = warp_res.cols;
            unsigned char* pr = new unsigned char[rht*rwd];
            unsigned char* pg = new unsigned char[rht*rwd];
            unsigned char* pb = new unsigned char[rht*rwd];
            int index = 0;
            for (int m = 0; m < rht; m++)
            {
                uchar *p = warp_res.ptr<uchar>(m);
                for (int n = 0; n < rwd; n++)
                {
                    pr[index] = p[n * 3];
                    pg[index] = p[n * 3 + 1];
                    pb[index] = p[n * 3 + 2];
                    index++;
                }
            }
            OGRSpatialReference oSRS;
            oSRS.SetUTM(mZoneNumber);
            oSRS.SetWellKnownGeogCS("WGS84");
            char    *pszWKT = NULL;
            oSRS.exportToWkt(&pszWKT);
            stGeoInfo geoinfo;
            geoinfo.left = gx_min;
            geoinfo.top = gy_max;
            geoinfo.dx = resolution;
            geoinfo.dy = -resolution;
            geoinfo.projectRef = pszWKT;
            GdalWriteImageColor("c:\\temp\\mosaic.tif", pb, pg, pr,
                rht, rwd, geoinfo);
            delete[] pr;
            delete[] pg;
            delete[] pb;

            //////////////////////////////////////////////////////////////////
			Mat imagecontours=imread("c:\\temp\\mosaic.tif");
			//result.copyTo(imagecontours);
			cv::cvtColor(imagecontours, imagecontours, CV_BGR2GRAY);
			cv::threshold(imagecontours, imagecontours, 5, 255, CV_THRESH_BINARY);
			cv::Mat element7(7, 7, CV_8U, cv::Scalar(1));
			cv::morphologyEx(imagecontours, imagecontours, cv::MORPH_OPEN, element7);
			cv::morphologyEx(imagecontours, imagecontours, cv::MORPH_CLOSE, element7);
			vector<Mat> contours;
			cv::findContours(imagecontours, contours, CV_RETR_EXTERNAL,//检索外部轮廓
				CV_CHAIN_APPROX_SIMPLE);//每个轮廓的拐点像素
			std::vector<cv::Point> vPoints;
			for (size_t i = 0; i < contours.size(); i++)
			{
				int nbrOfPoints = contours[i].rows;
				int cols = contours[i].cols;
				for (int j = 0; j < nbrOfPoints; ++j)
				{ 
					cv::Point temp;
					temp.x = contours[i].at<int>(j, 0);
					temp.y = contours[i].at<int>(j, 1);
					vPoints.push_back(temp);
				}
			}
			//计算轮廓像素对应经纬度坐标
			double  geoTransform[6];
			memset(geoTransform, 0, sizeof(double) * 6);
			geoTransform[0] = geoinfo.left;
			geoTransform[3] = geoinfo.top;
			geoTransform[1] = geoinfo.dx;
			geoTransform[5] = geoinfo.dy;
			std::vector<Lon2Pix> vLon2Pixs;
			int nXSize = rwd;// #列数
			int nYSize = rht; //#行数
			for (cv::Point i : vPoints) {
				cv::Point d = i;
				Lon2Pix lon2pix;
				lon2pix.x = d.x;
				lon2pix.y = d.y;
				lon2pix.height = 0;
				double UTMNorthing = geoTransform[0] + d.x*geoTransform[1] + d.y*geoTransform[2];
				double	UTMEasting = geoTransform[3] + d.x*geoTransform[4] + d.y*geoTransform[5];
				/*OGRSpatialReference *polatlon=oSRS.CloneGeogCS();
				OGRCoordinateTransformation *transformation=OGRCreateCoordinateTransformation(&oSRS,polatlon);
				if (transformation) {
				transformation->Transform(1, &UTMNorthing, &UTMEasting);
				}*/
				UTMtoLL(23, UTMEasting, UTMNorthing, mZoneNumber, lon2pix.lat, lon2pix.lon);
				vLon2Pixs.push_back(lon2pix);
			}

			if (hWndDisplay) {
				/*printf_s("save mosaic... \n");
				imwrite("c:\\temp\\mosaic.jpg", result);*/
				//通信无效 
				if (messageSendClient) {
					char vocfile[MAX_PATH];
					_getcwd(vocfile, MAX_PATH);
					strcat(vocfile, "\\send.xml");
					FileStorage file(vocfile, FileStorage::WRITE | FileStorage::FORMAT_XML);// | FileStorage::MEMORY);
					file << "MatchPoints" << "[";
					for (Lon2Pix i : vLon2Pixs) {
						Lon2Pix d = i;
						file << "{";
						file << "lon" << d.lon;
						file << "lat" << d.lat;
						file << "height" << d.height;
						file << "x" << d.x;
						file << "y" << d.y;
						file << "}";
					}
					file << "]";
					file << "MosaicFilePath" << "c:\\temp\\mosaic.tif";
					file << "ImageWidth" << warp_res.cols;
					file << "ImageHeight" << warp_res.rows;
					file.release();
					/*string str = file.releaseAndGetString();*/
					messageSendClient->SendMessageToServer((char*)vocfile, strlen(vocfile));
					//cv::waitKey(2);

				}
				ShowMat(warp_res, hWndDisplay);
			}
			else
			{
				/*printf_s("save mosaic... \n");
				imwrite("c:\\temp\\mosaic.jpg", result);*/
				if (messageSendClient) {
					char vocfile[MAX_PATH];
					_getcwd(vocfile, MAX_PATH);
					strcat(vocfile, "\\send.xml");
					FileStorage file(vocfile, FileStorage::WRITE | FileStorage::FORMAT_XML);// | FileStorage::MEMORY);
					file << "MatchPoints" << "[";
					for (Lon2Pix i : vLon2Pixs) {
						Lon2Pix d = i;
						file << "{";
						file << "lon" << d.lon;
						file << "lat" << d.lat;
						file << "height" << d.height;
						file << "x" << d.x;
						file << "y" << d.y;
						file << "}";
					}
					file << "]";
					file << "MosaicFilePath" << "c:\\temp\\mosaic.tif";
					file << "ImageWidth" << warp_res.cols;
					file << "ImageHeight" << warp_res.rows;
					file.release();
					/*string str = file.releaseAndGetString();*/
					messageSendClient->SendMessageToServer((char*)vocfile, strlen(vocfile));
				}
				int ht = warp_res.rows;
				int wd = warp_res.cols;
				int sht = 640;
				double scale = (double)(sht) / (double)(ht);
				int swd = wd * scale;
				cv::namedWindow("mosaic", WINDOW_NORMAL);
				cv::resizeWindow("mosaic", swd, sht);
				cv::imshow("mosaic", warp_res);
				cv::waitKey(100);
				/*printf("save mosaic... \n");
				imwrite("c:\\temp\\mosaic.jpg", result);*/
			}          			
		}		
	}	
	SetFinish();
}
int  MultiBandMap2DCPU::finalMosaic()
{
    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    int nKeyFrames = vpKFs.size();

    //camera parameters
    PinHoleParameters camParas;
    camParas.fx = mFx;
    camParas.fy = mFy;
    camParas.cx = mCx;
    camParas.cy = mCy;
    camParas.w = mWd;
    camParas.h = mHt;

    //1. plane fitting
    planeFit(true);

    //2. transform based on plane fitting
    printf("generating frames and rendering... \n");
    deque<std::pair<cv::Mat, pi::SE3d> > frames;
    for (int i = 0; i < nKeyFrames; i++)
    {
        KeyFrame* pKF = vpKFs[i];
        //printf("%d ", i);
        std::pair<cv::Mat, pi::SE3d> frame;
        frame.first = pKF->mFrameImage;

        double x, y, z, wx, wy, wz, w;
        cv::Mat R = pKF->GetRotation().t();
        cv::Mat t = pKF->GetCameraCenter();
        //cv::Mat t = pKF->GetTranslation();				
        {
            x = t.at<float>(0);
            y = t.at<float>(1);
            z = t.at<float>(2);
            vector<float> q = Converter::toQuaternion(R);
            wx = q[0];
            wy = q[1];
            wz = q[2];
            w = q[3];
        }
        printf("pos: %lf %lf %lf %lf %lf %lf %lf \n", x, y, z, wx, wy, wz, w);
        frame.second = mPlanePos.inverse()*pi::SE3d(x, y, z, wx, wy, wz, w);
        //frame.second = pi::SE3d(x, y, z, wx, wy, wz, w);
        frames.push_back(frame);
    }

    ////////////////// 3. prepare and rendering all frames //////
    printf("\n prepare  ... \n");
    prepare(mPlanePos, camParas, frames);
    for (int i = 0; i < nKeyFrames; i++)
    {
        renderFrame(frames[i]);
    }
    Mat result;
    int xoff, yoff;
    mosaic(result, xoff, yoff);

    ////////////////// 4. similarity transform /////////////////
    //camera locations
    vector<POINT3D> srcPts;
    vector<POINT3D> dstPts;
    int ht = result.rows;
    int wd = result.cols;
    double minx = data->min().x + xoff * data->eleSize();
    double miny = data->min().y + yoff * data->eleSize();
    double slon = vpKFs[0]->mGpsInfo.lon;
    mZoneNumber = int((slon + 180) / 6) + 1;
    for (int k = 0; k < nKeyFrames; k++)
    {
        KeyFrame* pKF = vpKFs[k];

        double x, y, z, wx, wy, wz, w;
        cv::Mat R = pKF->GetRotation().t();
        vector<float> q = Converter::toQuaternion(R);
        cv::Mat t = pKF->GetCameraCenter();
        //cv::Mat t = pKF->GetTranslation();
        {
            x = t.at<float>(0);
            y = t.at<float>(1);
            z = t.at<float>(2);
            wx = q[0];
            wy = q[1];
            wz = q[2];
            w = q[3];
        }
        pi::SE3d pos = pi::SE3d(x, y, z, wx, wy, wz, w);
        //std::pair<cv::Mat, pi::SE3d> frame(image, prepared->_plane.inverse()*pos);
        pi::SE3d tp = prepared->_plane.inverse()*pos;

        double cx = tp.get_translation().x;
        double cy = tp.get_translation().y;
        double ix = (cx - minx)*data->lengthPixelInv();
        double iy = (cy - miny)*data->lengthPixelInv();
        printf("camera location: %lf %lf \n", ix, iy);
        //cv::circle(result, Point(ix, iy), 3, CV_RGB(255, 0, 0));
        POINT3D sp;
        sp.x = ix - wd * 0.5;
        sp.y = ht * 0.5 - iy;
        srcPts.push_back(sp);

        //convert from lon,lat to utm
        POINT3D dstpt;
        double lat = pKF->mGpsInfo.lat;
        double lon = pKF->mGpsInfo.lon;
        printf("lat,lon: %lf %lf \n", lat, lon);
        double gz = pKF->mGpsInfo.height;
        double gx, gy;
        LLtoUTM(23, lat, lon, gy, gx, mZoneNumber);
        POINT3D dp;
        dp.x = gx;
        dp.y = gy;
        dp.z = gz;
        double la, lo;
        UTMtoLL(23, gy, gx, mZoneNumber, la, lo);
        printf("%lf %lf  %lf %lf \n", gy, gx, la, lo);
        dstPts.push_back(dp);
    }
    double sp[4];
    //SimilarityTransform(srcPts, dstPts, sp);
    SimRansac(srcPts, dstPts, sp);
    for (int k = 0; k < 4; k++)
    {
        mSP[k] = sp[k];
    }

    /////// 5. determine the scope of transformed image //////
    std::vector<cv::Point2f> imgPtsCV;
    imgPtsCV.resize(4);
    imgPtsCV[0].x = -wd * 0.5;   imgPtsCV[0].y = ht * 0.5;
    imgPtsCV[1].x = wd * 0.5;    imgPtsCV[1].y = ht * 0.5;
    imgPtsCV[2].x = -wd * 0.5;   imgPtsCV[2].y = -ht * 0.5;
    imgPtsCV[3].x = wd * 0.5;    imgPtsCV[3].y = -ht * 0.5;
    double nsp[4];
    double dnorm = sqrt(sp[0] * sp[0] + sp[1] * sp[1]);
    nsp[0] = sp[0] / dnorm;
    nsp[1] = sp[1] / dnorm;
    nsp[2] = 0;
    nsp[3] = 0;
    int minx1 = 100000;
    int maxx1 = -100000;
    int miny1 = 100000;
    int maxy1 = -100000;
    double gx_min = 100000000;
    double gx_max = -100000000;
    double gy_min = 100000000;
    double gy_max = -100000000;
    for (int k = 0; k < 4; k++)
    {
        int dx = SIMILARITY_X(imgPtsCV[k].x, imgPtsCV[k].y, nsp);
        int dy = SIMILARITY_Y(imgPtsCV[k].x, imgPtsCV[k].y, nsp);
        if (minx1 > dx) minx1 = dx;
        if (maxx1 < dx) maxx1 = dx;
        if (miny1 > dy) miny1 = dy;
        if (maxy1 < dy) maxy1 = dy;

        double tx = SIMILARITY_X(imgPtsCV[k].x, imgPtsCV[k].y, mSP);
        double ty = SIMILARITY_Y(imgPtsCV[k].x, imgPtsCV[k].y, mSP);
        if (gx_min > tx) gx_min = tx;
        if (gx_max < tx) gx_max = tx;
        if (gy_min > ty) gy_min = ty;
        if (gy_max < ty) gy_max = ty;
    }
    printf("%d %d %d %d \n", minx1, maxx1, miny1, maxy1);
    printf("%lf %lf %lf %lf \n", gx_min, gx_max, gy_min, gy_max);

    //////////// 6. resample the image based on transform //////////
    cv::Mat warp_res(maxy1 - miny1 + 1, maxx1 - minx1 + 1, CV_8UC3);
    int warpht = warp_res.rows;
    int warpwd = warp_res.cols;
    int row = 0;
    int col = 0;
    for (int j = miny1; j < maxy1; j++)
    {
        col = 0;
        uchar* pdst = warp_res.ptr<uchar>(warpht - row - 1);
        for (int i = minx1; i < maxx1; i++)
        {
            int sx = nsp[0] * i + nsp[1] * j; // SIMILARITY_X(i, j, nisp);
            int sy = -nsp[1] * i + nsp[0] * j; // SIMILARITY_Y(i, j, nisp);
            sx = sx + wd * 0.5;
            sy = ht * 0.5 - sy;
            if (sx < 0)   sx = 0;
            if (sx >= wd) sx = wd - 1;
            if (sy < 0)   sy = 0;
            if (sy >= ht) sy = ht - 1;
            uchar* psrc = result.ptr<uchar>(sy);
            pdst[col * 3] = psrc[sx * 3];
            pdst[col * 3 + 1] = psrc[sx * 3 + 1];
            pdst[col * 3 + 2] = psrc[sx * 3 + 2];
            col++;
        }
        row++;
    }

    /////////////// 7. add geoinfomation /////////////////
    double resolution = dnorm;
    int rht = warp_res.rows;
    int rwd = warp_res.cols;
    unsigned char* pr = new unsigned char[rht*rwd];
    unsigned char* pg = new unsigned char[rht*rwd];
    unsigned char* pb = new unsigned char[rht*rwd];
    int index = 0;
    for (int m = 0; m < rht; m++)
    {
        uchar *p = warp_res.ptr<uchar>(m);
        for (int n = 0; n < rwd; n++)
        {
            pr[index] = p[n * 3];
            pg[index] = p[n * 3 + 1];
            pb[index] = p[n * 3 + 2];
            index++;
        }
    }
    OGRSpatialReference oSRS;
    oSRS.SetUTM(mZoneNumber);
    oSRS.SetWellKnownGeogCS("WGS84");
    char    *pszWKT = NULL;
    oSRS.exportToWkt(&pszWKT);
    stGeoInfo geoinfo;
    geoinfo.left = gx_min;
    geoinfo.top = gy_max;
    geoinfo.dx = resolution;
    geoinfo.dy = -resolution;
    geoinfo.projectRef = pszWKT;
    GdalWriteImageColor("c:\\temp\\mosaic.tif", pb, pg, pr,
        rht, rwd, geoinfo);
    delete[] pr;
    delete[] pg;
    delete[] pb;

    printf("Final mosaic is finished ! \n");

    return 0;
}

void MultiBandMap2DCPU::draw()
{
	/*
    if(!_valid) return;

    SPtr<MultiBandMap2DCPUPrepare> p;
    SPtr<MultiBandMap2DCPUData>    d;
    {
        pi::ReadMutex lock(mutex);
        p=prepared;d=data;
    }
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrix(p->_plane);
    //draw deque frames
    pi::TicTac ticTac;
    ticTac.Tic();
    {
        std::deque<std::pair<cv::Mat,pi::SE3d> > frames=p->getFrames();
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        for(std::deque<std::pair<cv::Mat,pi::SE3d> >::iterator it=frames.begin();it!=frames.end();it++)
        {
            pi::SE3d& pose=it->second;
            glColor3ub(255,0,0);
            glVertex(pose.get_translation());
            glVertex(pose*pi::Point3d(1,0,0));
            glColor3ub(0,255,0);
            glVertex(pose.get_translation());
            glVertex(pose*pi::Point3d(0,1,0));
            glColor3ub(0,0,255);
            glVertex(pose.get_translation());
            glVertex(pose*pi::Point3d(0,0,1));
        }
        glEnd();
    }
    //draw global area
    if(svar.GetInt("Map2D.DrawArea"))
    {
        pi::Point3d _min=d->min();
        pi::Point3d _max=d->max();
        glColor3ub(255,0,0);
        glBegin(GL_LINES);
        glVertex3d(_min.x,_min.y,0);
        glVertex3d(_min.x,_max.y,0);
        glVertex3d(_min.x,_min.y,0);
        glVertex3d(_max.x,_min.y,0);
        glVertex3d(_max.x,_min.y,0);
        glVertex3d(_max.x,_max.y,0);
        glVertex3d(_min.x,_max.y,0);
        glVertex3d(_max.x,_max.y,0);
        glEnd();
    }

    //draw textures
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
//    glEnable(GL_LIGHTING);
    if(alpha)
    {
        glEnable(GL_ALPHA_TEST);
        glAlphaFunc(GL_GREATER, 0.1f);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    }
    GLint last_texture_ID;
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture_ID);
    std::vector<SPtr<MultiBandMap2DCPUEle> > dataCopy=d->data();
    int wCopy=d->w(),hCopy=d->h();
    glColor3ub(255,255,255);
    for(int x=0;x<wCopy;x++)
        for(int y=0;y<hCopy;y++)
        {
            int idxData=y*wCopy+x;
            float x0=d->min().x+x*d->eleSize();
            float y0=d->min().y+y*d->eleSize();
            float x1=x0+d->eleSize();
            float y1=y0+d->eleSize();
            SPtr<MultiBandMap2DCPUEle> ele=dataCopy[idxData];
            if(!ele.get())  continue;
            {
                {
                    pi::ReadMutex lock(ele->mutexData);
                    if(!(ele->pyr_laplace.size()&&ele->weights.size()
                         &&ele->pyr_laplace.size()==ele->weights.size())) continue;
                    if(ele->Ischanged)
                    {
                        pi::timer.enter("MultiBandMap2DCPU::updateTexture");
                        bool updated=false,inborder=false;
                        if(_highQualityShow)
                        {
                            vector<SPtr<MultiBandMap2DCPUEle> > neighbors;
                            neighbors.reserve(9);
                            for(int yi=y-1;yi<=y+1;yi++)
                                for(int xi=x-1;xi<=x+1;xi++)
                                {
                                    if(yi<0||yi>=hCopy||xi<0||xi>=wCopy)
                                    {
                                        neighbors.push_back(SPtr<MultiBandMap2DCPUEle>());
                                        inborder=true;
                                    }
                                    else neighbors.push_back(dataCopy[yi*wCopy+xi]);
                                }
                            updated=ele->updateTexture(neighbors);
                        }
                        else
                            updated=ele->updateTexture();
                        pi::timer.leave("MultiBandMap2DCPU::updateTexture");

                        if(updated&&!inborder&&svar.GetInt("Fuse2Google"))
                        {
                            pi::timer.enter("MultiBandMap2DCPU::fuseGoogle");
                            stringstream cmd;
                            pi::Point3d  worldTl=p->_plane*pi::Point3d(x0,y0,0);
                            pi::Point3d  worldBr=p->_plane*pi::Point3d(x1,y1,0);
                            pi::Point3d  gpsTl,gpsBr;
                            pi::calcLngLatFromDistance(d->gpsOrigin().x,d->gpsOrigin().y,worldTl.x,worldTl.y,gpsTl.x,gpsTl.y);
                            pi::calcLngLatFromDistance(d->gpsOrigin().x,d->gpsOrigin().y,worldBr.x,worldBr.y,gpsBr.x,gpsBr.y);
//                            cout<<"world:"<<worldBr<<"origin:"<<d->gpsOrigin()<<endl;
                            cmd<<"Map2DUpdate LastTexMat "<< setiosflags(ios::fixed)
                              << setprecision(9)<<gpsTl<<" "<<gpsBr;
//                            cout<<cmd.str()<<endl;
                            scommand.Call("MapWidget",cmd.str());
                            pi::timer.leave("MultiBandMap2DCPU::fuseGoogle");

                        }
                    }
                }
                if(ele->texName)
                {
                    glBindTexture(GL_TEXTURE_2D,ele->texName);
                    glBegin(GL_QUADS);
                    glTexCoord2f(0.0f, 0.0f); glVertex3f(x0,y0,0);
                    glTexCoord2f(0.0f, 1.0f); glVertex3f(x0,y1,0);
                    glTexCoord2f(1.0f, 1.0f); glVertex3f(x1,y1,0);
                    glTexCoord2f(1.0f, 0.0f); glVertex3f(x1,y0,0);
                    glEnd();
                }
            }
        }
    glBindTexture(GL_TEXTURE_2D, last_texture_ID);
    glPopMatrix();
	*/
}

bool MultiBandMap2DCPU::mosaic(Mat& result, int& xoff, int& yoff){

	// determin minmax
	SPtr<MultiBandMap2DCPUPrepare> p;
	SPtr<MultiBandMap2DCPUData>    d;
	{
		//pi::ReadMutex lock(mutex);
		p = prepared;
		d = data;
	}
	if (d->w() == 0 || d->h() == 0) 
		return false;

	pi::Point2i minInt(1e6, 1e6), maxInt(-1e6, -1e6);
	int contentCount = 0;
	for (int x = 0; x<d->w(); x++)
		for (int y = 0; y<d->h(); y++)
		{
			SPtr<MultiBandMap2DCPUEle> ele = d->data()[x + y*d->w()];
			if (!ele.get()) continue;
			{
				//pi::ReadMutex lock(ele->mutexData);
				if (!ele->pyr_laplace.size()) continue;
			}
			contentCount++;
			minInt.x = min(minInt.x, x); minInt.y = min(minInt.y, y);
			maxInt.x = max(maxInt.x, x); maxInt.y = max(maxInt.y, y);
		}

	printf("%d \n", _bandNum);
	maxInt = maxInt + pi::Point2i(1, 1);

    //added by xiedonghai, 2018.9.14
    xoff = minInt.x;
    yoff = minInt.y;

    
	pi::Point2i wh = maxInt - minInt;

    if (wh.x < 1 || wh.y < 1)
        return false;

	printf("debugtestwhx:%d why:%d\n", wh.x,wh.y);
	vector<cv::Mat> pyr_laplace(_bandNum + 1);
	vector<cv::Mat> pyr_weights(_bandNum + 1);
	for (int i = 0; i <= 0; i++)
		pyr_weights[i] = cv::Mat::zeros(wh.y*ELE_PIXELS, wh.x*ELE_PIXELS, CV_32FC1);
	for (int x = minInt.x; x<maxInt.x; x++)
		for (int y = minInt.y; y<maxInt.y; y++)
		{
			//printf("%d %d \n", x, y);
			SPtr<MultiBandMap2DCPUEle> ele = d->data()[x + y*d->w()];
			if (!ele.get()) continue;
			{
				//pi::ReadMutex lock(ele->mutexData);
				if (!ele->pyr_laplace.size()) continue;
				int width = ELE_PIXELS, height = ELE_PIXELS;

				for (int i = 0; i <= _bandNum; ++i)
				{
					cv::Rect rect(width*(x - minInt.x), height*(y - minInt.y), width, height);
					if (pyr_laplace[i].empty())
						pyr_laplace[i] = cv::Mat::zeros(wh.y*height, wh.x*width, ele->pyr_laplace[i].type());
					ele->pyr_laplace[i].copyTo(pyr_laplace[i](rect));
					if (i == 0)
						ele->weights[i].copyTo(pyr_weights[i](rect));
					height >>= 1; width >>= 1;
				}
			}
		}
	
	//printf("laplace..\n");
	//printf("%d %d \n", pyr_laplace.)
	cv::detail::restoreImageFromLaplacePyr(pyr_laplace);

	result = pyr_laplace[0].clone();
	if (result.type() == CV_16SC3) 
		result.convertTo(result, CV_8UC3);
	result.setTo(cv::Scalar::all(0), pyr_weights[0] == 0);
	
	printf("mosaic finished! \n");
	return true;
}

int MultiBandMap2DCPU::saveAsTiff(string filename)
{



    return 0;
}

bool MultiBandMap2DCPU::save(const std::string& filename)
{
    // determin minmax
    SPtr<MultiBandMap2DCPUPrepare> p;
    SPtr<MultiBandMap2DCPUData>    d;
    {
        //pi::ReadMutex lock(mutex);
        p=prepared;
		d=data;
    }
    if(d->w()==0||d->h()==0) return false;

    pi::Point2i minInt(1e6,1e6),maxInt(-1e6,-1e6);
    int contentCount=0;
    for(int x=0;x<d->w();x++)
        for(int y=0;y<d->h();y++)
        {
            SPtr<MultiBandMap2DCPUEle> ele=d->data()[x+y*d->w()];
            if(!ele.get()) continue;
            {
                //pi::ReadMutex lock(ele->mutexData);
                if(!ele->pyr_laplace.size()) continue;
            }
            contentCount++;
            minInt.x=min(minInt.x,x); minInt.y=min(minInt.y,y);
            maxInt.x=max(maxInt.x,x); maxInt.y=max(maxInt.y,y);
        }

    maxInt=maxInt+pi::Point2i(1,1);
    pi::Point2i wh=maxInt-minInt;
    vector<cv::Mat> pyr_laplace(_bandNum+1);
    vector<cv::Mat> pyr_weights(_bandNum+1);
    for(int i=0;i<=0;i++)
        pyr_weights[i]=cv::Mat::zeros(wh.y*ELE_PIXELS,wh.x*ELE_PIXELS,CV_32FC1);

    for(int x=minInt.x;x<maxInt.x;x++)
        for(int y=minInt.y;y<maxInt.y;y++)
        {
            SPtr<MultiBandMap2DCPUEle> ele=d->data()[x+y*d->w()];
            if(!ele.get()) continue;
            {
                //pi::ReadMutex lock(ele->mutexData);
                if(!ele->pyr_laplace.size()) continue;
                int width=ELE_PIXELS,height=ELE_PIXELS;

                for (int i = 0; i <= _bandNum; ++i)
                {
                    cv::Rect rect(width*(x-minInt.x),height*(y-minInt.y),width,height);
                    if(pyr_laplace[i].empty())
                        pyr_laplace[i]=cv::Mat::zeros(wh.y*height,wh.x*width,ele->pyr_laplace[i].type());
                    ele->pyr_laplace[i].copyTo(pyr_laplace[i](rect));
                    if(i==0)
                        ele->weights[i].copyTo(pyr_weights[i](rect));
                    height>>=1;width>>=1;
                }
            }
        }

    cv::detail::restoreImageFromLaplacePyr(pyr_laplace);

    cv::Mat result=pyr_laplace[0];
    if(result.type()==CV_16SC3) result.convertTo(result,CV_8UC3);
    result.setTo(cv::Scalar::all(svar.GetInt("Result.BackGroundColor")),pyr_weights[0]==0);
    cv::imwrite(filename,result);
    cout<<"Resolution:["<<result.cols<<" "<<result.rows<<"]";
    if(svar.exist("GPS.Origin"))
          cout<<",_lengthPixel:"<<d->lengthPixel()
       <<",Area:"<<contentCount*d->eleSize()*d->eleSize()<<endl;
    return true;
}


void MultiBandMap2DCPU::RequestFinish()
{
	unique_lock<mutex> lock(mMutexFinish);
	mbFinishRequested = true;
}

bool MultiBandMap2DCPU::CheckFinish()
{
	unique_lock<mutex> lock(mMutexFinish);
	return mbFinishRequested;
}

void MultiBandMap2DCPU::SetFinish()
{
	unique_lock<mutex> lock(mMutexFinish);
	mbFinished = true;
}

bool MultiBandMap2DCPU::isFinished()
{
	unique_lock<mutex> lock(mMutexFinish);
	return mbFinished;
}

void MultiBandMap2DCPU::RequestStop()
{
	unique_lock<mutex> lock(mMutexStop);
	if (!mbStopped)
		mbStopRequested = true;
}

bool MultiBandMap2DCPU::isStopped()
{
	unique_lock<mutex> lock(mMutexStop);
	return mbStopped;
}

bool MultiBandMap2DCPU::Stop()
{
	unique_lock<mutex> lock(mMutexStop);
	unique_lock<mutex> lock2(mMutexFinish);

	if (mbFinishRequested)
		return false;
	else if (mbStopRequested)
	{
		mbStopped = true;
		mbStopRequested = false;
		return true;
	}

	return false;

}

void MultiBandMap2DCPU::Release()
{
	unique_lock<mutex> lock(mMutexStop);
	mbStopped = false;
}


//****** absolute orientation, added by xiedonghai, 2018.8.31 **********
//if (bIsAbsOri)
//{
//    double slon = vpKFs[0]->mGpsInfo.lon;
//    int zoneNumber = int((slon + 180) / 6) + 1;
//    vector<POINT3D> srcpts;
//    vector<POINT3D> dstpts;
//    for (int i = 0; i < nKeyFrames; i++)
//    {
//        KeyFrame* pKF = vpKFs[i];
//        double lat = pKF->mGpsInfo.lat;
//        double lon = pKF->mGpsInfo.lon;
//        printf("lat,lon: %lf %lf \n", lat, lon);

//        //convert from lon,lat to utm
//        double gz = pKF->mGpsInfo.height;
//        double gx, gy;
//        LLtoUTM(23, lat, lon, gy, gx, zoneNumber);
//        POINT3D dp;
//        dp.x = gx;
//        dp.y = gy;
//        dp.z = gz;
//        dstpts.push_back(dp);

//        POINT3D sp;
//        cv::Mat t = pKF->GetCameraCenter();
//        sp.x = t.at<float>(0);
//        sp.y = t.at<float>(1);
//        sp.z = t.at<float>(2);
//        srcpts.push_back(sp);

//        printf("grd pt:  %lf %lf %lf \n", dp.x, dp.y, dp.z);
//        printf("free pt: %lf %lf %lf \n", sp.x, sp.y, sp.z);
//    }

//    //stAbsPOS absparas;
//    RansanAbsOri(srcpts, dstpts, mAbsParas);

//    FILE* fp = fopen("c:\\temp\\traj.txt", "w");
//    for (int k = 0; k < srcpts.size(); k++)
//    {
//        double fP[3];
//        fP[0] = srcpts[k].x;
//        fP[1] = srcpts[k].y;
//        fP[2] = srcpts[k].z;
//        double tp[3];
//        matrix_product(3, 3, 3, 1, mAbsParas.R, fP, tp);
//        for (int k = 0; k<3; k++)
//            fP[k] = mAbsParas.scale*tp[k] + mAbsParas.T[k];
//        
//        fprintf(fp, "%lf %lf %lf   %lf %lf %lf  %lf %lf %lf  \n",
//            srcpts[k].x, srcpts[k].y, srcpts[k].z,
//            dstpts[k].x, dstpts[k].y, dstpts[k].z,
//            fp[0], fp[1], fp[2]);
//    }
//    fclose(fp);
//}


