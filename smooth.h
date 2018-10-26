#ifndef  PC_SMOOTH_H
#define  PC_SMOOTH_H

//  numNeibor=60,  stddev_mult=0.6
int PointCloudFilter(double* px, double* py, double* pz, int nPt,
	int numNeibor, double stddev_mult);


#endif

