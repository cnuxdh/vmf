
#include "smooth.h"


//pcl library
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/statistical_outlier_removal.h>

int PointCloudFilter(double* px, double* py, double* pz, int nPt,
	int numNeibor, double stddev_mult)
{
	//prepare the 3D points
	//pcl::PointCloud<pcl::PointXYZ> ptCloud;
	pcl::PointCloud<pcl::PointXYZ>::Ptr ptCloud(new pcl::PointCloud<pcl::PointXYZ>);

	//pcl::PointXYZ pt;
	ptCloud->height = 1;
	ptCloud->width = nPt;
	ptCloud->points.resize(nPt);
	for (int i = 0; i<nPt; i++)
	{
		ptCloud->points[i].x = px[i];
		ptCloud->points[i].y = py[i];
		ptCloud->points[i].z = pz[i];
	}

	//reading the file into the cloud structure
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);

	std::cerr << "Cloud before filtering: " << std::endl;
	std::cerr << *ptCloud << std::endl;

	// Create the filtering object
	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
	sor.setInputCloud(ptCloud);
	sor.setMeanK(numNeibor);
	sor.setStddevMulThresh(stddev_mult);
	sor.filter(*cloud_filtered);

	std::cerr << "Cloud after filtering: " << std::endl;
	std::cerr << *cloud_filtered << std::endl;

	//output the smoothed points
	int nSmoothPt = 0;
	nSmoothPt = cloud_filtered->size();
	for (int i = 0; i<nSmoothPt; i++)
	{
		px[i] = cloud_filtered->points[i].x;
		py[i] = cloud_filtered->points[i].y;
		pz[i] = cloud_filtered->points[i].z;
	}

	return nSmoothPt;
}


//
//int PointCloudFilter(double* px, double* py, double* pz, int nPt,
//	int numNeibor, double stddev_mult)
//{
//	//prepare the 3D points
//	//pcl::PointCloud<pcl::PointXYZ> ptCloud;
//	pcl::PointCloud<pcl::PointXYZ>::Ptr ptCloud(new pcl::PointCloud<pcl::PointXYZ>);
//
//	//pcl::PointXYZ pt;
//	ptCloud->height = 1;
//	ptCloud->width = nPt;
//	ptCloud->points.resize(nPt);
//	for(int i=0; i<nPt; i++)
//	{
//		ptCloud->points[i].x = px[i];
//		ptCloud->points[i].y = py[i];
//		ptCloud->points[i].z = pz[i];
//	}
//
//	//reading the file into the cloud structure
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
//
//	std::cerr << "Cloud before filtering: " << std::endl;
//	std::cerr << *ptCloud << std::endl;
//
//	// Create the filtering object
//	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
//	sor.setInputCloud(ptCloud);
//	sor.setMeanK(numNeibor);
//	sor.setStddevMulThresh(stddev_mult);
//	sor.filter(*cloud_filtered);
//
//	std::cerr << "Cloud after filtering: " << std::endl;
//	std::cerr << *cloud_filtered << std::endl;
//
//
//	//pcl::PCDWriter writer;
//	//writer.write<pcl::PointXYZ> ("d:\\table_scene_lms400_inliers.pcd", *cloud_filtered, false);
//
//	//sor.setNegative (true);
//	//sor.filter (*cloud_filtered);
//	//writer.write<pcl::PointXYZ> ("d:\\table_scene_lms400_outliers.pcd", *cloud_filtered, false);
//
//	//output the smoothed points
//	int nSmoothPt = 0;
//	nSmoothPt = cloud_filtered->size();
//	for(int i=0; i<nSmoothPt; i++)
//	{
//		px[i] = cloud_filtered->points[i].x;
//		py[i] = cloud_filtered->points[i].y;
//		pz[i] = cloud_filtered->points[i].z;
//	}
//
//	return nSmoothPt;
//}
