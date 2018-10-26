#ifndef FIT_OBJECT_H
#define FIT_OBJECT_H



double FitPlane(double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror);
double FitObject(int mode, double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror, double *outx, double *outy, double *outz, int &infopnum);
bool   ComputeParameter(double *inx, double *iny, double *inz, int *pairarray, int *objectinfo, int innum, double *para, int &paranum);



#endif