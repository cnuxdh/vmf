

//#include "API_FitObject.h"
//#include "API_MathFunc.h"
//#include "API_Fitting.h"

#include "Global.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "FitObject.h"



double FitLine(double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror)
{
	// Z!=0 ( X - para[1] ) / para[0] = ( Y - para[3] ) / para[2] = Z / 1
	// Z==0  para[0]*X + Y + para[1] = 0
	int i=0;
	double A[4]={0.0}, B[4]={0.0}, C[4]={0.0};
	double zsum=0.0;
	double err = 0.0;
	double maxerr=0.0;
	// 判断 Z 是否全部为0
	for (i=0;i<innum;i++)
	{
		zsum += fabs(inz[i]);
	}
	
	if (zsum == 0.0)
	{
		for (i=0;i<innum;i++)
		{
			A[0] += inx[i]*inx[i];
			A[1] += inx[i];
			A[2] += inx[i];
			A[3] += 1;
			
			B[0] -= (inx[i]*iny[i]);
			B[1] -= iny[i];
		}
		//brinv(A,2);
		//brmul(A,B,2,2,1,C);
		invers_matrix(A,2);
		mult(A,B,C,2,2,1);

		para[0] = C[0];
		para[1] = C[1];
		
		for (i=0;i<innum;i++)
		{
			// para[2]*X + para[0]*Y - 2*para[0]*para[2]*Z - (para[2]*para[1] + para[0]*para[3])
			double terr = fabs(para[0]*inx[i]+iny[i]+para[1])
				/sqrt(para[0]*para[0]+1);
			(maxerr<terr)?maxerr = terr:maxerr;
			err += terr;
			pointerror[i] = terr;
		}
		paranum = 2;
		err = err/innum;
	}
	else
	{
		for (i=0;i<innum;i++)
		{
			A[0] += inx[i]*inz[i];
			A[1] += inx[i];
			A[2] += iny[i]*inz[i];
			A[3] += iny[i];
			
			B[0] += inz[i]*inz[i];
			B[1] += inz[i];
			B[2] += inz[i];
			B[3]  = innum;
		}
		//brinv(B,2);
		//brmul(A,B,2,2,2,C);
		
		invers_matrix(B,2);
		mult(A,B,C,2,2,2);

		para[0] = C[0];
		para[1] = C[1];
		para[2] = C[2];
		para[3] = C[3];
		
		for (i=0;i<innum;i++)
		{
			// para[2]*X + para[0]*Y - 2*para[0]*para[2]*Z - (para[2]*para[1] + para[0]*para[3])
			double terr = fabs(para[2]*inx[i]+para[0]*iny[i]-2.0*para[0]*para[2]*inz[i]-(para[2]*para[1]+para[0]*para[3]))
				/sqrt(para[2]*para[2]+para[0]*para[0]+4.0*para[0]*para[0]*para[2]*para[2]);
			(maxerr<terr)?maxerr = terr:maxerr;
			err += terr;
			pointerror[i] = terr;
		}
		err = err/innum;
		paranum = 4;
	}
	outerrinfo[0] = maxerr;
	outerrinfo[1] = err;		
	return err;
}


double FitPlane(double *inx, double *iny, double *inz, 
				int innum, double *para, int &paranum, 
				double *outerrinfo, double *pointerror)
{
	int i=0;
	double A[9]={0.0}, B[3]={0.0}, C[3]={0.0};
	
	double xsum=0.0, ysum=0.0, zsum=0.0;
	for (i=0;i<innum;i++)
	{
		xsum += fabs(inx[i]);
		ysum += fabs(iny[i]);
		zsum += fabs(inz[i]);
	}
	if (xsum == 0.0 && ysum != 0.0 && zsum != 0.0)
	{
		para[0] = 1.0;
		para[1] = 0.0;
		para[2] = 0.0;
		para[3] = 0.0;
	}
	else if (xsum != 0.0 && ysum == 0.0 && zsum != 0.0)
	{
		para[0] = 0.0;
		para[1] = 1.0;
		para[2] = 0.0;
		para[3] = 0.0;
	}
	else if (xsum != 0.0 && ysum != 0.0 && zsum == 0.0)
	{
		para[0] = 0.0;
		para[1] = 0.0;
		para[2] = 1.0;
		para[3] = 0.0;
	}
	else if(xsum != 0.0 && ysum != 0.0 && zsum != 0.0)
	{
		for (i=0;i<innum;i++)
		{
			A[0] += inx[i]*inx[i];
			A[1] += inx[i]*iny[i];
			A[2] += inx[i];
			A[3] += inx[i]*iny[i];
			A[4] += iny[i]*iny[i];
			A[5] += iny[i];
			A[6] += inx[i];
			A[7] += iny[i];
			A[8] += 1.0;
			
			
			B[0] -= (inx[i]*inz[i]);
			B[1] -= (iny[i]*inz[i]);
			B[2] -= inz[i];
		}
		
		//brinv(A,3);
		//brmul(A,B,3,3,1,C);
		invers_matrix(A,3);
		mult(A,B,C,3,3,1);
		
		double totcom = sqrt(C[0]*C[0]+C[1]*C[1]+1.0);
		
		para[0] = C[0]/totcom;
		para[1] = C[1]/totcom;
		para[2] = 1.0/totcom;
		para[3] = C[2]/totcom;
		
		
	}
	else
	{
		para[0] = 0.0;
		para[1] = 0.0;
		para[2] = 0.0;
		para[3] = 0.0;
	}
	
	
	double err=0.0;
	double maxerr=0.0;
	for (i=0;i<innum;i++)
	{
		double terr = fabs(para[0]*inx[i]+para[1]*iny[i]+para[2]*inz[i]+para[3])
			/sqrt(para[0]*para[0]+para[1]*para[1]+para[2]*para[2]);
		(maxerr<terr)?maxerr=terr:maxerr;
		err += terr;
		pointerror[i] = terr;
	}
	err = err/innum;
	outerrinfo[0] = maxerr;
	outerrinfo[1] = err;
	paranum = 4;
		
	return err;
}


double FitCircle(double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror)
{
	double *intpx = new double[innum];
	double *intpy = new double[innum];
	double *intpz = new double[innum];
	Point3D ctr,vect;
	bool flag[1] = {true};
	double rad=0.0;
	double *errarray = new double[innum];
	int i=0;
	

	
	for (i=0;i<innum;i++)
	{
		intpx[i] = inx[i];
		intpy[i] = iny[i];
		intpz[i] = inz[i];
	}


	// 需要重写
//	Circle3DFitting(intpx,intpy,intpz,flag,innum,ctr,rad,vect,errarray);
//////////////////////////////////////////////////////////////////////////

	double *planepara = new double[20];
	double *tmpouterr = new double[innum];
	double *tmpdataerr = new double[innum];
	FitPlane(intpx,intpy,intpz,innum,planepara,paranum,tmpouterr,tmpdataerr);
	double ppx,ppy,ppz;
	if (planepara[2]<0.0)
	{
		ppx = -1.0*planepara[0];
		ppy = -1.0*planepara[1];
		ppz = -1.0*planepara[2];
	}
	else
	{
		ppx = planepara[0];
		ppy = planepara[1];
		ppz = planepara[2];
	}
	vect.x = ppx;
	vect.y = ppy;
	vect.z = ppz;


	// 旋转数据使其与XY面平行
	double R[9],RT[9];
	R[0] = 1 - ppx*ppx/(1+ppz) ;
	R[1] = -ppx*ppy/(1+ppz) ;
	R[2] = -ppx ;
	R[3] = -ppx*ppy/(1+ppz) ;
	R[4] = 1 - ppy*ppy/(1+ppz) ;
	R[5] = -ppy ;
	R[6] = ppx ;
	R[7] = ppy ;
	R[8] = ppz ;

	//Trans( R, 3, 3, RT ) ;
	transpose(R,RT,3,3);

	// 旋转后的数据
	double *pX = new double [innum] ;
	double *pY = new double [innum] ;
	double *pZ = new double [innum] ;
	
	for( i=0; i<innum ; i++ )
	{
		double temp1[3] ; // 原始数据
		double temp2[3] ; // 旋转后数据
		temp1[0] = intpx[i];
		temp1[1] = intpy[i] ;
		temp1[2] = intpz[i];
		
		//brmul( R, temp1, 3, 3 ,1 , temp2 ) ;
		mult(R,temp1,temp2,3,3,1);
		
		pX[i] = temp2[0];
		pY[i] = temp2[1] ;
		pZ[i] = temp2[2] ;
	}

	// 在二维空间内拟合椭圆
	double cvx, cvy  ; 
	
	double cvz = 0.0 ;
	for( i=0 ; i<innum ; i++ )
	{
		cvz+= pZ[i]/innum ;
	}
	double *pDis = new double [innum];
	//Circle2DFitting( pX , pY , innum , cvx, cvy , rad ,  pDis ) ;
	double center[3] ;
	center[0] = cvx ;
	center[1] = cvy ;
	center[2] = cvz ;
	
	double center_org[3];
	//brmul( RT, center, 3 , 3 , 1 , center_org ) ;
	mult(RT,center,center_org,3,3,1);
	
	// 输出结果
	ctr.x = center_org[0];
	ctr.y = center_org[1] ;
	ctr.z = center_org[2] ;


	delete planepara;
	delete tmpouterr;
	delete tmpdataerr;
	delete pX;
	delete pY;
	delete pZ;
	delete pDis;

//////////////////////////////////////////////////////////////////////////




	double err = 0.0, maxerr = 0.0;
	double thisplane[4];
	// 由圆心点, 先求出 m*X+n*Y+p*Z+q = 0 的平面
	thisplane[0] = vect.x;
	thisplane[1] = vect.y;
	thisplane[2] = vect.z;
	thisplane[3] = -1.0*(vect.x*ctr.x+vect.y*ctr.y+vect.z*ctr.z);



	for (i=0;i<innum;i++)
	{
		double toplanedist = fabs(thisplane[0]*intpx[i]+thisplane[1]*intpy[i]+thisplane[2]*intpz[i]+thisplane[3])
			/sqrt(thisplane[0]*thisplane[0]+thisplane[1]*thisplane[1]+thisplane[2]*thisplane[2]);
		double onpx,onpy,onpz;
		onpz = thisplane[0]*vect.x*intpz[i]-thisplane[0]*vect.z*intpx[i]+thisplane[1]*vect.y*intpz[i]-thisplane[1]*vect.z*intpy[i]-thisplane[3]*vect.z;
		onpz = onpz/(thisplane[0]*vect.x+thisplane[1]*vect.y+thisplane[2]*vect.z);

		onpx = vect.x*(onpz-intpz[i])/vect.z+intpx[i];
		onpy = vect.y*(onpz-intpz[i])/vect.z+intpy[i];

		double onplanedist = sqrt((onpx-ctr.x)*(onpx-ctr.x)+(onpy-ctr.y)*(onpy-ctr.y)+(onpz-ctr.z)*(onpz-ctr.z));
		onplanedist = fabs(onplanedist-rad);

		double terr = sqrt(onplanedist*onplanedist+toplanedist*toplanedist);

		(maxerr<terr)?maxerr=terr:maxerr;
		err += terr;
		pointerror[i] = terr;
	}
	err = err/innum;
	outerrinfo[0] = maxerr;
	outerrinfo[1] = err;



	para[0] = ctr.x;
	para[1] = ctr.y;
	para[2] = ctr.z;
	para[3] = vect.x;
	para[4] = vect.y;
	para[5] = vect.z;
	para[6] = rad;
	paranum = 7;


	delete errarray;
	delete intpx;
	delete intpy;
	delete intpz;
	errarray = NULL;
	intpx = NULL;
	intpy = NULL;
	intpz = NULL;

	return err;
}



BOOL Ellipse2DFitting(double *EdgeX, double *EdgeY, int nPoint,  CEllipseCoeff *pCoeff)
{
// 椭圆方程系数: x*x + b*x*y +c*y*y + d*x + e*y + f = 0 
	int i,j;
	double px,py;
	//用椭圆拟合来找到中心
	static	double A[5][5],b[5];

	if(nPoint<5)
		return FALSE;
	int n=3;
	for(i=0;i<5;i++)
		{
			b[i]=0;
			for(j=0;j<5;j++)
				A[i][j]=0;
		}

	for(i=0;i<nPoint;i++)
	{
		px=EdgeX[i];
		py=EdgeY[i];
		A[0][0]+=px*px*py*py;
		A[0][1]+=px*py*py*py;
		A[0][2]+=px*px*py;
		A[0][3]+=px*py*py;
		A[0][4]+=px*py;

		A[1][0]+=px*py*py*py;
		A[1][1]+=py*py*py*py;
		A[1][2]+=px*py*py;
		A[1][3]+=py*py*py;
		A[1][4]+=py*py;

		A[2][0]+=px*px*py;
		A[2][1]+=px*py*py;
		A[2][2]+=px*px;
		A[2][3]+=px*py;
		A[2][4]+=px;

		A[3][0]+=px*py*py;
		A[3][1]+=py*py*py;
		A[3][2]+=px*py;
		A[3][3]+=py*py;
		A[3][4]+=py;

		A[4][0]+=px*py;
		A[4][1]+=py*py;
		A[4][2]+=px;
		A[4][3]+=py;
		A[4][4]+=1;

		b[0]-=px*px*px*py;
		b[1]-=px*px*py*py;
		b[2]-=px*px*px;
		b[3]-=px*px*py;
		b[4]-=px*px;
	}

	//computing the five coefficients  of ellipse 
	double CenterX,CenterY;
	double h, c, d , e , f ;
	//if(agaus((double *)A,b,5)!=0)
	{
		h=b[0];
		c=b[1];
		d=b[2];
		e=b[3];
		f=b[4];
		pCoeff->b = h ;
		pCoeff->c = c ;
		pCoeff->d = d ;
		pCoeff->e = e ;
		pCoeff->f = f ;
		CenterX = (h*e-2*c*d)/(4*c-h*h);
	    CenterY = (h*d-2*e)/(4*c-h*h);
	}
	
	//计算椭圆的半长轴和半短轴
	double a1,c1,f1;
	if(h >= 0)
		a1 = (1+c+sqrt((1-c)*(1-c)+h*h))/2;
	else
		a1 = (1+c-sqrt((1-c)*(1-c)+h*h))/2;
	
	c1 = (1+c)-a1;
	f1 = f+(CenterX*d+CenterY*e)/2;

	// 再加上边界点到拟合距离判断

	if((-f1/a1 >0) && (-f1/c1>0))
	{
		return TRUE;
	}

	return FALSE;
}


// 从椭圆方程系数得到椭圆参数
BOOL GetEllipseInfoFromEllipseCoefficient( CEllipseCoeff *pCoeff , CEllipseInfo *pEllipseInfo )
{
	double h, c , d , e , f ;
	double CenterX , CenterY ,HalfLong, HalfShort , Angle;

	h = pCoeff->b ;
	c = pCoeff->c ;
	d = pCoeff->d ;
	e = pCoeff->e ;
	f = pCoeff->f ;

	CenterX = (h*e-2*c*d)/(4*c-h*h);
	CenterY = (h*d-2*e)/(4*c-h*h);

	// 倾角计算
	// tg(2*Alpha) = b/(1-c);
	Angle = atan( h/(1-c) )/2;
	
	double a2, c2 ;
	a2 = cos(Angle)*cos(Angle) + h*sin(Angle)*cos(Angle) + c*sin(Angle)*sin(Angle) ;
	c2 = sin(Angle)*sin(Angle) - h*sin(Angle)*cos(Angle) + c*cos(Angle)*cos(Angle) ;

	//计算椭圆的半长轴和半短轴
	double a1,c1,f1;
	a1 = cos(Angle)*cos(Angle) + h*sin(Angle)*cos(Angle) + c*sin(Angle)*sin(Angle) ;
	c1 = sin(Angle)*sin(Angle) - h*sin(Angle)*cos(Angle) + c*cos(Angle)*cos(Angle) ;
	Angle = Angle*180/PI;

/*
	if(h >= 0)
		a1 = ( 1+c+sqrt((1-c)*(1-c)+h*h) )/2;
	else
		a1 = ( 1+c-sqrt((1-c)*(1-c)+h*h) )/2;
*/	
//	c1 = (1+c)-a1;
	double f2 ;
	f2 = f+(CenterX*d+CenterY*e)/2;
	f1 = f+(CenterX*CenterX + h*CenterX*CenterY + c*CenterY*CenterY + d*CenterX + e*CenterY ) ;
	// 再加上边界点到拟合距离判断

	if((-f1/a1 >0) && (-f1/c1>0))
	{
		HalfLong = sqrt(-f1/a1);
		HalfShort = sqrt(-f1/c1);

		pEllipseInfo->cx = CenterX ;
		pEllipseInfo->cy = CenterY ;
		pEllipseInfo->a  = HalfLong ;
		pEllipseInfo->b  = HalfShort ;
		pEllipseInfo->alpha = Angle ;
		return TRUE;
	}

	return FALSE;
}



double FitElipse(double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror)
{
	double *intpx = new double[innum];
	double *intpy = new double[innum];
	double *intpz = new double[innum];
	Point3D ctr,vect;
	double axis_a, axis_b;
	double *errarray = new double[innum];
	int i=0;
	
	for (i=0;i<innum;i++)
	{
		intpx[i] = inx[i];
		intpy[i] = iny[i];
		intpz[i] = inz[i];
	}

//////////////////////////////////////////////////////////////////////////
	double *planepara = new double[20];
	double *tmpouterr = new double[innum];
	double *tmpdataerr = new double[innum];
	FitPlane(intpx,intpy,intpz,innum,planepara,paranum,tmpouterr,tmpdataerr);
	double ppx,ppy,ppz;
	if (planepara[2]<0.0)
	{
		ppx = -1.0*planepara[0];
		ppy = -1.0*planepara[1];
		ppz = -1.0*planepara[2];
	}
	else
	{
		ppx = planepara[0];
		ppy = planepara[1];
		ppz = planepara[2];
	}
	vect.x = ppx;
	vect.y = ppy;
	vect.z = ppz;


	// 旋转数据使其与XY面平行
	double R[9],RT[9];
	R[0] = 1 - ppx*ppx/(1+ppz) ;
	R[1] = -ppx*ppy/(1+ppz) ;
	R[2] = -ppx ;
	R[3] = -ppx*ppy/(1+ppz) ;
	R[4] = 1 - ppy*ppy/(1+ppz) ;
	R[5] = -ppy ;
	R[6] = ppx ;
	R[7] = ppy ;
	R[8] = ppz ;
	//Trans( R, 3, 3, RT ) ;
	transpose(R,RT,3,3);

	// 旋转后的数据
	double *pX = new double [innum] ;
	double *pY = new double [innum] ;
	double *pZ = new double [innum] ;
	
	for( i=0; i<innum ; i++ )
	{
		double temp1[3] ; // 原始数据
		double temp2[3] ; // 旋转后数据
		temp1[0] = intpx[i];
		temp1[1] = intpy[i] ;
		temp1[2] = intpz[i];
		
		//brmul( R, temp1, 3, 3 ,1 , temp2 ) ;
		mult(R,temp1,temp2,3,3,1);
		
		pX[i] = temp2[0];
		pY[i] = temp2[1] ;
		pZ[i] = temp2[2] ;
	}


	CEllipseCoeff Coeff;
	BOOL bEllipse = Ellipse2DFitting(pX, pY, innum,  &Coeff);
	
	CEllipseInfo EllipseInfo;
	BOOL bEllipseInfo = GetEllipseInfoFromEllipseCoefficient( &Coeff , &EllipseInfo );




	// 在二维空间内拟合椭圆
	double cvx=0.0, cvy=0.0  ; 
	
	double cvz = 0.0 ;
	for( i=0 ; i<innum ; i++ )
	{
		cvz+= pZ[i]/innum ;
	}

	double center[3] ;
	center[0] = EllipseInfo.cx ;
	center[1] = EllipseInfo.cy ;
	center[2] = cvz ;
	
	double center_org[3];
	//brmul( RT, center, 3 , 3 , 1 , center_org ) ;
	mult(RT,center, center_org,3,3,1);
	
	// 输出结果
	ctr.x = center_org[0];
	ctr.y = center_org[1];
	ctr.z = center_org[2];
	
	axis_a = EllipseInfo.a;
	axis_b = EllipseInfo.b;




	delete planepara;
	delete tmpouterr;
	delete tmpdataerr;
	delete pX;
	delete pY;
	delete pZ;


//////////////////////////////////////////////////////////////////////////




	double err = 0.0, maxerr = 0.0;
	double thisplane[4];
	// 由圆心点, 先求出 m*X+n*Y+p*Z+q = 0 的平面
	thisplane[0] = vect.x;
	thisplane[1] = vect.y;
	thisplane[2] = vect.z;
	thisplane[3] = -1.0*(vect.x*ctr.x+vect.y*ctr.y+vect.z*ctr.z);


	for (i=0;i<innum;i++)
	{
		double toplanedist = fabs(thisplane[0]*intpx[i]+thisplane[1]*intpy[i]+thisplane[2]*intpz[i]+thisplane[3])
			/sqrt(thisplane[0]*thisplane[0]+thisplane[1]*thisplane[1]+thisplane[2]*thisplane[2]);
		double onpx,onpy,onpz;
		onpz = thisplane[0]*vect.x*intpz[i]-thisplane[0]*vect.z*intpx[i]+thisplane[1]*vect.y*intpz[i]-thisplane[1]*vect.z*intpy[i]-thisplane[3]*vect.z;
		onpz = onpz/(thisplane[0]*vect.x+thisplane[1]*vect.y+thisplane[2]*vect.z);

		onpx = vect.x*(onpz-intpz[i])/vect.z+intpx[i];
		onpy = vect.y*(onpz-intpz[i])/vect.z+intpy[i];


		double onplanedist = axis_a*100.0;
		double tmpxyzin[3], tmpxyzout[3], tmpxyzcenter[3];
		tmpxyzin[0] = onpx;
		tmpxyzin[1] = onpy;
		tmpxyzin[2] = onpz;
		//brmul(R,tmpxyzin,3,3,1,tmpxyzout);
		tmpxyzin[0] = ctr.x;
		tmpxyzin[1] = ctr.y;
		tmpxyzin[2] = ctr.z;
		//brmul(R,tmpxyzin,3,3,1,tmpxyzcenter);
		for (double j=0.0;j<PI*2.0;j+=0.01)
		{
			double ptx = axis_a*cos(j)+tmpxyzcenter[0];
			double pty = axis_b*sin(j)+tmpxyzcenter[1];
			double ptz = tmpxyzcenter[2];
			double tmpplanedist = sqrt((tmpxyzout[0]-ptx)*(tmpxyzout[0]-ptx)+(tmpxyzout[1]-pty)*(tmpxyzout[1]-pty)+(tmpxyzout[2]-ptz)*(tmpxyzout[2]-ptz));
			if (tmpplanedist<onplanedist)
			{
				onplanedist = tmpplanedist;
			}
		}



		double terr = sqrt(onplanedist*onplanedist+toplanedist*toplanedist);

		(maxerr<terr)?maxerr=terr:maxerr;
		err += terr;
		pointerror[i] = terr;
	}
	err = err/innum;
	outerrinfo[0] = maxerr;
	outerrinfo[1] = err;


	para[0] = ctr.x;
	para[1] = ctr.y;
	para[2] = ctr.z;
	para[3] = vect.x;
	para[4] = vect.y;
	para[5] = vect.z;
	para[6] = axis_a;
	para[7] = axis_b;
	paranum = 8;


	delete errarray;
	delete intpx;
	delete intpy;
	delete intpz;
	errarray = NULL;
	intpx = NULL;
	intpy = NULL;
	intpz = NULL;


	return err;
}

double FitSphere(double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror)
{
	// x*x + y*y + z*z + para[0]*x + para[1]*y + para[2]*z + para[3] = 0
	int i=0;
	double A[16]={0.0}, B[4]={0.0}, C[4]={0.0};

	for (i=0;i<innum;i++)
	{
		A[0] += inx[i]*inx[i];
		A[1] += inx[i]*iny[i];
		A[2] += inx[i]*inz[i];
		A[3] += inx[i];

		A[4] += inx[i]*iny[i];
		A[5] += iny[i]*iny[i];
		A[6] += iny[i]*inz[i];
		A[7] += iny[i];
		
		A[8] += inx[i]*inz[i];
		A[9] += iny[i]*inz[i];
		A[10] += inz[i]*inz[i];
		A[11] += inz[i];

		A[12] += inx[i];
		A[13] += iny[i];
		A[14] += inz[i];
		A[15] += 1.0;
		
		B[0] -= (inx[i]*(inx[i]*inx[i] + iny[i]*iny[i] + inz[i]*inz[i]));
		B[1] -= (iny[i]*(inx[i]*inx[i] + iny[i]*iny[i] + inz[i]*inz[i]));
		B[2] -= (inz[i]*(inx[i]*inx[i] + iny[i]*iny[i] + inz[i]*inz[i]));
		B[3] -= (inx[i]*inx[i] + iny[i]*iny[i] + inz[i]*inz[i]);
	}
	
	//brinv(A,4);
	//brmul(A,B,4,4,1,C);
	invers_matrix(A,4);
	mult(A,B,C,4,4,1);

	para[0] = C[0];
	para[1] = C[1];
	para[2] = C[2];
	para[3] = C[3];
	paranum = 4;


	double err=0.0;
	double maxerr=0.0;
	double crx,cry,crz,crr;
	crx = -para[0]/2.0;
	cry = -para[1]/2.0;
	crz = -para[2]/2.0;
	crr = sqrt(crx*crx+cry*cry+crz*crz-para[3]);
	
	for (i=0;i<innum;i++)
	{
		double terr = fabs( sqrt( (crx-inx[i])*(crx-inx[i]) + (cry-iny[i])*(cry-iny[i]) + (crz-inz[i])*(crz-inz[i]) ) - crr );
		(maxerr<terr)?maxerr=terr:maxerr;
		err += terr;
		pointerror[i] = terr;
	}
	err = err/innum;
	outerrinfo[0] = maxerr;
	outerrinfo[1] = err;


	return err;
}

double FitCylinder(double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror)
{
	int i=0;
	double Rad[1],errt[1];

	Point3D *tp3d = new Point3D[innum];
	Point3D ctr[1], vect[1];
	for (i=0;i<innum;i++)
	{
		tp3d[i].x = inx[i];
		tp3d[i].y = iny[i];
		tp3d[i].z = inz[i];
	}



	double *dis = new double[innum];
	//Cylinder3DFitting(tp3d,innum,ctr,vect,Rad,errt,dis);
	delete dis;
	dis = NULL;
	// 误差分析
	double maxerr = 0.0, err = 0.0, cyheight=0.0;


	// (X - Cx)/Vx = (Y - Cy)/Vy = (Z - Cz)/Vz，点到该直线距离与半径的差为误差
	for (i=0;i<innum;i++)
	{


// 		para[0] = vect[0].x/vect[0].z;
// 		para[1] = -1.0*para[0]*ctr[0].z+ctr[0].x;
// 		para[2] = vect[0].y/vect[0].z;
// 		para[3] = -1.0*para[2]*ctr[0].z+ctr[0].y;
// 
// 		// para[2]*X + para[0]*Y - 2*para[0]*para[2]*Z - (para[2]*para[1] + para[0]*para[3])
// 		double terr = fabs(para[2]*inx[i]+para[0]*iny[i]-2.0*para[0]*para[2]*inz[i]-(para[2]*para[1]+para[0]*para[3]))
// 			/sqrt(para[2]*para[2]+para[0]*para[0]+4.0*para[0]*para[0]*para[2]*para[2]);
// 		terr = fabs(terr-Rad[0]);
// 		(maxerr<terr)?maxerr=terr:maxerr;
// 		pointerror[i] = terr;
// 		err += terr;



		// 新的2010.2.24, SBlW，直接求点到中间截面的距离
		
		double vx = vect[0].x, vy = vect[0].y, vz = vect[0].z;
		double cx = ctr[0].x, cy = ctr[0].y, cz = ctr[0].z;
		double vd = -1.0*(vx*cx+vy*cy+vz*cz);
// 		double chuz = (para[0]*inx[i]+para[2]*iny[i]+inz[i]-para[0]*para[1]-para[2]*para[3])
// 					/(para[0]*para[0]+para[1]*para[1]+1);
// 		double chux = para[0]*chuz+para[1];
// 		double chuy = para[2]*chuz+para[3];

		double chuz = vz*vz*inz[i]+vy*vy*cy+vx*vx*cx-vy*vz*(cy-iny[i])-vx*vz*(cx-inx[i]);
		double chux = vx/vz*(chuz-cx)+cx;
		double chuy = vy/vz*(chuz-cy)+cy;
		double distxyz = sqrt((chux-ctr[0].x)*(chux-ctr[0].x)+(chuy-ctr[0].y)*(chuy-ctr[0].y)+(chuz-ctr[0].z)*(chuz-ctr[0].z));
	
		distxyz = fabs(vx*inx[i]+vy*iny[i]+vz*inz[i]+vd)/sqrt(vx*vx+vy*vy+vz*vz);
		
		cyheight<distxyz?cyheight = distxyz:cyheight;


		// 误差的求取，已知和中心点的距离，点到中间截面的距离，直角三角形

		double tocenterdist = sqrt((inx[i]-cx)*(inx[i]-cx)+(iny[i]-cy)*(iny[i]-cy)+(inz[i]-cz)*(inz[i]-cz));
		double terr = sqrt(tocenterdist*tocenterdist-distxyz*distxyz);
 		terr = fabs(terr-Rad[0]);
		(maxerr<terr)?maxerr=terr:maxerr;
		pointerror[i] = terr;
 		err += terr;
	}

	err = err/innum;
	outerrinfo[0] = maxerr;
	outerrinfo[1] = err;


	para[0] = ctr[0].x;
	para[1] = ctr[0].y;
	para[2] = ctr[0].z;
	para[3] = vect[0].x;
	para[4] = vect[0].y;
	para[5] = vect[0].z;
	para[6] = Rad[0];
	para[7] = cyheight;
	paranum = 8;

	delete tp3d;
	tp3d = NULL;

	return errt[0];
}



bool GetInfoPointsfromFitLine(double *inx, double *iny, double *inz, int innum, double *para, int mode, double *outx, double *outy, double *outz)
{
	int i=0;
	double *otx = (double *)calloc(sizeof(double),innum);
	double *oty = (double *)calloc(sizeof(double),innum);
	double *otz = (double *)calloc(sizeof(double),innum);

	if (mode==1)
	{
		// 平面上, para[0]*X + Y + para[1] = 0
		for (i=0;i<innum;i++)
		{
			double p2 = iny[i] - inx[i]/para[0];
			otx[i] = (-p2-para[1])/(para[0]+1.0/para[0]);
			oty[i] = -para[0]*otx[i]-para[1];
			otz[i] = 0.0;
		}
	}
	else if (mode==2)
	{
		// Z!=0
		for (i=0;i<innum;i++)
		{
			// X = para[0]*Z + para[1]
			// Y = para[2]*Z + para[3]
			// para[0]*a + para[2]*b + 1 = 0;
			double p11 = para[0]*inx[i] + para[2]*iny[i] + inz[i] - para[0]*para[1] - para[2]*para[3];
			double p22 = para[0]*para[0] + para[2]*para[2] + 1.0;

			otz[i] = p11/p22;
			otx[i] = para[0]*otz[i] + para[1];
			oty[i] = para[2]*otz[i] + para[3];

		}
	}

	double minx=otx[0], miny=oty[0], minz=otz[0], maxx=otx[0], maxy=oty[0], maxz=otz[0];
	for (i=1;i<innum;i++)
	{
		if (minx>otx[i])
		{
			minx = otx[i];
			miny = oty[i];
			minz = otz[i];
		}
		else if (minx==otx[i])
		{
			if(miny>oty[i])
			{
				minx = otx[i];
				miny = oty[i];
				minz = otz[i];
			}
			else if (miny==oty[i])
			{
				if(minz>otz[i])
				{
					minx = otx[i];
					miny = oty[i];
					minz = otz[i];
				}
			}
		}
	}
	for (i=0;i<innum;i++)
	{
		if (maxx<otx[i])
		{
			maxx = otx[i];
			maxy = oty[i];
			maxz = otz[i];
		}
		else if (maxx==otx[i])
		{
			if(maxy>oty[i])
			{
				maxx = otx[i];
				maxy = oty[i];
				maxz = otz[i];
			}
			else if (maxy==oty[i])
			{
				if(maxz>otz[i])
				{
					maxx = otx[i];
					maxy = oty[i];
					maxz = otz[i];
				}
			}
		}
	}

	outx[0] = minx;
	outy[0] = miny;
	outz[0] = minz;
	outx[1] = maxx;
	outy[1] = maxy;
	outz[1] = maxz;
	free(otx);
	free(oty);
	free(otz);


	return 1;
}


bool GetInfoPointsfromFitPlane(double *inx, double *iny, double *inz, int innum, double *para, int mode, double *outx, double *outy, double *outz)
{
	// 垂足为：z=( (p0*p0+p1*p1)*inz - p0*p2*inx - p1*p2*iny - p2*p3 )/(p0*p0 + p1*p1 + p2*p2)
	int i;
	double *otx = (double *)calloc(innum,sizeof(double));
	double *oty = (double *)calloc(innum,sizeof(double));
	double *otz = (double *)calloc(innum,sizeof(double));
	for (i=0;i<innum;i++)
	{
		otz[i] = ((para[0]*para[0] + para[1]*para[1])*inz[i]
			- para[0]*para[2]*inx[i] - para[1]*para[2]*iny[i] - para[2]*para[3])
			/(para[0]*para[0] + para[1]*para[1] + para[2]*para[2]);
		otx[i] = para[0]/para[2]*(otz[i] - inz[i]) + inx[i];
		oty[i] = para[1]/para[2]*(otz[i] - inz[i]) + iny[i];
	}
	
	double minx=otx[0], miny=oty[0], minz=otz[0], maxx=otx[0], maxy=oty[0], maxz=otz[0];
	for (i=1;i<innum;i++)
	{
		if (minx>otx[i])
		{
			minx = otx[i];
		}
		if (maxx<otx[i])
		{
			maxx = otx[i];
		}
	}
	for (i=1;i<innum;i++)
	{
		if(miny>oty[i])
		{
			miny = oty[i];
		}
		if(maxy<oty[i])
		{
			maxy = oty[i];
		}
	}
	outx[0] = minx;
	outy[0] = miny;
	outz[0] = ( para[0]*outx[0] + para[1]*outy[0] + para[3] )/(0.0-para[2]);
	
	outx[1] = minx;
	outy[1] = maxy;
	outz[1] = ( para[0]*outx[1] + para[1]*outy[1] + para[3] )/(0.0-para[2]);
	
	outx[2] = maxx;
	outy[2] = maxy;
	outz[2] = ( para[0]*outx[2] + para[1]*outy[2] + para[3] )/(0.0-para[2]);
	
	outx[3] = maxx;
	outy[3] = miny;
	outz[3] = ( para[0]*outx[3] + para[1]*outy[3] + para[3] )/(0.0-para[2]);
	
	return 1;
}



bool GetInfoPointsfromFitCircle(double *inx, double *iny, double *inz, int innum, double *para, int mode, double *outx, double *outy, double *outz )
{
	outx[0] = para[0];
	outy[0] = para[1];
	outz[0] = para[2];

	return 1;
}


bool GetInfoPointsfromFitElipse(double *inx, double *iny, double *inz, int innum, double *para, int mode, double *outx, double *outy, double *outz )
{
	outx[0] = para[0];
	outy[0] = para[1];
	outz[0] = para[2];
	
	return 1;
}


bool GetInfoPointsfromFitSphere(double *inx, double *iny, double *inz, int innum, double *para, int mode, double *outx, double *outy, double *outz )
{
	double crx,cry,crz,crr;
	crx = -para[0]/2.0;
	cry = -para[1]/2.0;
	crz = -para[2]/2.0;
	crr = sqrt(crx*crx+cry*cry+crz*crz-para[3]);

	outx[0] = crx;
	outy[0] = cry;
	outz[0] = crz;
	outx[1] = crr;
	outy[1] = crr;
	outz[1] = crr;

	return 1;
}

bool GetInfoPointsfromFitCylinder(double *inx, double *iny, double *inz, int innum, double *para, int mode, double *outx, double *outy, double *outz )
{
	// 圆柱的高度如何获得
	// para,012 CxCyCz, 345 Vectxyz, 6 Radius 7 Height
	outx[0] = para[0];
	outy[0] = para[1];
	outz[0] = para[2];			// 圆柱中心
	// 圆柱底面 Vx*Cx + Vy*Cy + Vz*Cz + D = 0
	double D1 = -1.0*(para[0]*para[3] + para[1]*para[4] + para[2]*para[5]);
	double tx1,tx2,ty1,ty2,tz,ta,tb,tc,tt;
	tz = para[2]+1.0;								// attention，这个选择要注意
	ta = para[3]*para[3]+para[4]*para[4];
	tb = 2.0*((para[0]+para[5]*para[2]+D1)*para[4]-para[1]*para[3]*para[3]);
	tc = (para[0]+para[2]*para[5]+D1)*(para[0]+para[2]*para[5]+D1)
		+para[1]*para[1]*para[3]*para[3]-para[3]*para[3]*para[6]*para[6];
	tt = tb*tb-4.0*ta*tc;
	if (tt>0)
	{
		ty1 = (-tb+sqrt(tt))/2.0*ta;
		ty2 = (-tb-sqrt(tt))/2.0*ta;
	}
	tx1 = -1.0*(tz*para[5]+D1+para[4]*ty1)/para[3];
	tx2 = -1.0*(tz*para[5]+D1+para[4]*ty2)/para[3];

	double Cx,Cy,Cz,VX,VY,VZ,CR,thx1,thx2,thx3,thx4,thy1,thy2,thy3,thy4,thz1,thz2,thz3,thz4;
	// 以(tx1, ty1, tz) 为起始点，设为0度，选取底面圆上4点
	// 第二个点是180度处点
	//(X-Cx)/(thx1-Cx) = (Y-Cy)/(thy1-Cy) = (Z-Cz)/(thz1-Cz)
	Cx = para[0];
	Cy = para[1];
	Cz = para[2];
	VX = para[3];
	VY = para[4];
	VZ = para[5];
	CR = para[6];
	thx1 = tx1;
	thy1 = ty1;
	thz1 = tz;

	thz2 = (thz1-Cz)*(thz1-Cz)*CR*CR
		/((thx1-Cx)*(thx1-Cx)+(thy1-Cy)*(thy1-Cy)+(thz1-Cz)*(thz1-Cz));
	thz2 = sqrt(thz2)+Cz;
	thx2 = (thz2-Cz)/(thz1-Cz)*(thx1-Cx)+Cx;
	thy2 = (thz2-Cz)/(thz1-Cz)*(thy1-Cy)+Cy;

	// 以下正确
	double vx0,vy0=1.0,vz0=2.0;
	vx0 = -1.0*(VY*vy0+VZ*vz0)/VX;

	thz1 = sqrt(CR*CR*vz0*vz0/(vx0*vx0+vy0*vy0+vz0*vz0))+Cz;
	thx1 = vx0/vz0*(thz1-Cz)+Cx;
	thy1 = vy0/vz0*(thz1-Cz)+Cy;
		
	thz2 = -1.0*sqrt(CR*CR*vz0*vz0/(vx0*vx0+vy0*vy0+vz0*vz0))+Cz;
	thx2 = vx0/vz0*(thz2-Cz)+Cx;
	thy2 = vy0/vz0*(thz2-Cz)+Cy;

	// 两点
	double vx1,vy1,vz1=1.0;
	vy1 = (VX*(thz1-thz2) + VZ*(thx2-thx1))
		/(VX*(thy2-thy1) + VY*(thx1-thx2));
	vx1 = -1.0*(VZ+VY*vy1)/VX;
	// (X-Cx)/vx1 = (Y-Cy)/vy1 = (Z-Cz)/vz1
	thz3 = sqrt(CR*CR*vz1*vz1/(vx1*vx1+vy1*vy1+vz1*vz1))+Cz;
	thx3 = vx1/vz1*(thz3-Cz)+Cx;
	thy3 = vy1/vz1*(thz3-Cz)+Cy;

	thz4 = -1.0*sqrt(CR*CR*vz1*vz1/(vx1*vx1+vy1*vy1+vz1*vz1))+Cz;
	thx4 = vx1/vz1*(thz4-Cz)+Cx;
	thy4 = vy1/vz1*(thz4-Cz)+Cy;
		

	//////////////////////////////////////////////////////////////////////////
	outx[1] = thx1;
	outy[1] = thy1;
	outz[1] = thz1;

	outx[2] = thx2;
	outy[2] = thy2;
	outz[2] = thz2;

	outx[3] = thx3;
	outy[3] = thy3;
	outz[3] = thz3;

	outx[4] = thx4;
	outy[4] = thy4;
	outz[4] = thz4;
	//////////////////////////////////////////////////////////////////////////

	// 上表面5个点
	double tvd = sqrt(VX*VX+VY*VY+VZ*VZ)*para[7];		// Z - Z'

	outz[5] = tvd*VZ + outz[0];
	outx[5] = tvd*VX + outx[0];
	outy[5] = tvd*VY + outy[0];

	outz[6] = tvd*VZ + outz[1];
	outx[6] = tvd*VX + outx[1];
	outy[6] = tvd*VY + outy[1];

	outz[7] = tvd*VZ + outz[2];
	outx[7] = tvd*VX + outx[2];
	outy[7] = tvd*VY + outy[2];

	outz[8] = tvd*VZ + outz[3];
	outx[8] = tvd*VX + outx[3];
	outy[8] = tvd*VY + outy[3];

	outz[9] = tvd*VZ + outz[4];
	outx[9] = tvd*VX + outx[4];
	outy[9] = tvd*VY + outy[4];
	//////////////////////////////////////////////////////////////////////////
	// 下表面5个点
	outz[10] = -1.0*tvd*VZ + outz[0];
	outx[10] = -1.0*tvd*VX + outx[0];
	outy[10] = -1.0*tvd*VY + outy[0];
		
	outz[11] = -1.0*tvd*VZ + outz[1];
	outx[11] = -1.0*tvd*VX + outx[1];
	outy[11] = -1.0*tvd*VY + outy[1];
		
	outz[12] = -1.0*tvd*VZ + outz[2];
	outx[12] = -1.0*tvd*VX + outx[2];
	outy[12] = -1.0*tvd*VY + outy[2];
		
	outz[13] = -1.0*tvd*VZ + outz[3];
	outx[13] = -1.0*tvd*VX + outx[3];
	outy[13] = -1.0*tvd*VY + outy[3];
		
	outz[14] = -1.0*tvd*VZ + outz[4];
	outx[14] = -1.0*tvd*VX + outx[4];
	outy[14] = -1.0*tvd*VY + outy[4];

	// 底面法向量
	outx[15] = VX;
	outy[15] = VY;
	outz[15] = VZ;

	return 1;
}


// 出口函数
double FitObject(int mode, double *inx, double *iny, double *inz, int innum, double *para, int &paranum, double *outerrinfo, double *pointerror, double *outx, double *outy, double *outz, int &infopnum)
{
	double err = 0.0;
	
	switch (mode)
	{
	case 1:
		err = FitLine(inx,iny,inz,innum,para,paranum,outerrinfo,pointerror);
		GetInfoPointsfromFitLine(inx,iny,inz,innum,para,paranum/2,outx,outy,outz);
		infopnum = 2;
		break;
	case 2:
		err = FitPlane(inx,iny,inz,innum,para,paranum,outerrinfo,pointerror);
		GetInfoPointsfromFitPlane(inx,iny,inz,innum,para,1,outx,outy,outz);
		infopnum = 4;
		break;
	case 3:
		err = FitCircle(inx,iny,inz,innum,para,paranum,outerrinfo,pointerror);
		GetInfoPointsfromFitCircle(inx,iny,inz,innum,para,1,outx,outy,outz);
		infopnum = 1;
		break;
	case 4:
		err = FitElipse(inx,iny,inz,innum,para,paranum,outerrinfo,pointerror);
		GetInfoPointsfromFitElipse(inx,iny,inz,innum,para,1,outx,outy,outz);
		infopnum = 1;
		break;
	case 5:
		err = FitSphere(inx,iny,inz,innum,para,paranum,outerrinfo,pointerror);
		GetInfoPointsfromFitSphere(inx,iny,inz,innum,para,1,outx,outy,outz);
		infopnum = 2;
		break;
	case 6:
		err = FitCylinder(inx,iny,inz,innum,para,paranum,outerrinfo,pointerror);
		GetInfoPointsfromFitCylinder(inx,iny,inz,innum,para,1,outx,outy,outz);
		infopnum = 16;
		break;
	default:
		break;
	}
	
	return err;
}





///////////////////////////////////////////////////////////////////////////////



void cstrq(double *a,int n,double *q,double *b,double *c) //约化实对称矩阵为对称三对角阵 利用Householder变换
//a--双精度实型二维数组，体积n X n，存放n阶实对称矩阵A
//n--整形变量。实对称矩阵A的阶数
//q--双精度实型二维数组，体积n X n。返回Householder变换的乘积矩阵Q
//b--双精度实型一维数组，长度n，返回对称三对角阵的主对角元素
//c--双精度实型一维数组，长度n，其中前n-1个元素返回对称三对角阵中的次对角线元素
// int n;
//double a[],q[],b[],c[];
{ int i,j,k,u;//,v;
double h,f,g,h2;
for (i=0; i<=n-1; i++)
for (j=0; j<=n-1; j++)
{ u=i*n+j; q[u]=a[u];}
for (i=n-1; i>=1; i--)
{ h=0.0;
if (i>1)
for (k=0; k<=i-1; k++)
{ u=i*n+k; h=h+q[u]*q[u];}
if (h+1.0==1.0)
{ c[i]=0.0;
if (i==1) c[i]=q[i*n+i-1];
b[i]=0.0;
}
else
{ c[i]=sqrt(h);
u=i*n+i-1;
if (q[u]>0.0) c[i]=-c[i];
h=h-q[u]*c[i];
q[u]=q[u]-c[i];
f=0.0;
for (j=0; j<=i-1; j++)
{ q[j*n+i]=q[i*n+j]/h;
g=0.0;
for (k=0; k<=j; k++)
g=g+q[j*n+k]*q[i*n+k];
if (j+1<=i-1)
for (k=j+1; k<=i-1; k++)
g=g+q[k*n+j]*q[i*n+k];
c[j]=g/h;
f=f+g*q[j*n+i];
}
h2=f/(h+h);
for (j=0; j<=i-1; j++)
{ f=q[i*n+j];
g=c[j]-h2*f;
c[j]=g;
for (k=0; k<=j; k++)
{ u=j*n+k;
q[u]=q[u]-f*c[k]-g*q[i*n+k];
}
}
b[i]=h;
}
}
for (i=0; i<=n-2; i++) c[i]=c[i+1];
c[n-1]=0.0;
b[0]=0.0;
for (i=0; i<=n-1; i++)
{ if ((b[i]!=0.0)&&(i-1>=0))
for (j=0; j<=i-1; j++)
{ g=0.0;
for (k=0; k<=i-1; k++)
g=g+q[i*n+k]*q[k*n+j];
for (k=0; k<=i-1; k++)
{ u=k*n+j;
q[u]=q[u]-g*q[k*n+i];
}
}
u=i*n+i;
b[i]=q[u]; q[u]=1.0;
if (i-1>=0)
for (j=0; j<=i-1; j++)
{ q[i*n+j]=0.0; q[j*n+i]=0.0;}
}
return;
}

int csstq(int n,double *b,double *c,double *q,double eps,int l) 
//用变形QR方法计算实对称三对角阵的全部特征值与相应的特征向量
//n--整型变量。实对称三对角阵的阶数
//b--双精度实型一维数组，长度n。存放n阶对称三对角阵的主对角线上的各元素；返回时存放全部特征值
//c--双精度一维数组，长度n。前n-1个元素存放n阶对称三对角阵的次对角线上的元素
//q--双精度实型二维数组，体积n X n。返回特征向量组
//eps--双精度实型变量。迭代过程中控制精度要求
//l--整型变量。最大迭代次数
//int n,l;
//double b[],c[],q[],eps;
{ int i,j,k,m,it,u,v;
double d,f,h,g,p,r,e,s;
c[n-1]=0.0; d=0.0; f=0.0;
for (j=0; j<=n-1; j++)
{ it=0;
h=eps*(fabs(b[j])+fabs(c[j]));
if (h>d) d=h;
m=j;
while ((m<=n-1)&&(fabs(c[m])>d)) m=m+1;
if (m!=j)
{ do
{ if (it==l)
{ printf("fail\n");
return(-1);
}
it=it+1;
g=b[j];
p=(b[j+1]-g)/(2.0*c[j]);
r=sqrt(p*p+1.0);
if (p>=0.0) b[j]=c[j]/(p+r);
else b[j]=c[j]/(p-r);
h=g-b[j];
for (i=j+1; i<=n-1; i++)
b[i]=b[i]-h;
f=f+h; p=b[m]; e=1.0; s=0.0;
for (i=m-1; i>=j; i--)
{ g=e*c[i]; h=e*p;
if (fabs(p)>=fabs(c[i]))
{ e=c[i]/p; r=sqrt(e*e+1.0);
c[i+1]=s*p*r; s=e/r; e=1.0/r;
}
else
{ e=p/c[i]; r=sqrt(e*e+1.0);
c[i+1]=s*c[i]*r;
s=1.0/r; e=e/r;
}
p=e*b[i]-s*g;
b[i+1]=h+s*(e*g+s*b[i]);
for (k=0; k<=n-1; k++)
{ u=k*n+i+1; v=u-1;
h=q[u]; q[u]=s*q[v]+e*h;
q[v]=e*q[v]-s*h;
}
}
c[j]=s*p; b[j]=e*p;
}
while (fabs(c[j])>d);
}
b[j]=b[j]+f;
}
for (i=0; i<=n-1; i++)
{ k=i; p=b[i];
if (i+1<=n-1)
{ j=i+1;
while ((j<=n-1)&&(b[j]<=p))
{ k=j; p=b[j]; j=j+1;}
}
if (k!=i)
{ b[k]=b[i]; b[i]=p;
for (j=0; j<=n-1; j++)
{ u=j*n+i; v=j*n+k;
p=q[u]; q[u]=q[v]; q[v]=p;
}
}
}
return(1);
}





int Function_Foury(Point3D *A, Point3D *B, int n, double *R, double *T)
{
	double* a = new double[3*n];					// x 0-i, y 0-i, z 0-i
	double* b = new double[3*n];
	double averA[3]={0.0};
	double averB[3]={0.0};			// average of x or y or z
	int i,j;
	for (i=0;i<n;i++)
	{
		a[i] = A[i].x;
		a[i+n] = A[i].y;
		a[i+2*n] = A[i].z;
		averA[0] += A[i].x/((double)n);
		averA[1] += A[i].y/((double)n);
		averA[2] += A[i].z/((double)n);
		b[i] = B[i].x;
		b[i+n] = B[i].y;
		b[i+2*n] = B[i].z;
		averB[0] += B[i].x/((double)n);
		averB[1] += B[i].y/((double)n);
		averB[2] += B[i].z/((double)n);
	}
	double* ap = new double[3*n];							//3*n
	double* bp = new double[3*n];							//3*n
	for (i=0;i<3*n;i++)
	{
		ap[i] = a[i] - averA[(int)(i/n)];
		bp[i] = b[i] - averB[(int)(i/n)];
	}
//////////////////////////////////////////////////////////////////////////			above is for setting a,b,averA,averB,ap,bp
	double* Trbp = new double[3*n];						// matrix bp transfer to bp t	 as n*3
	//Trans(bp,3,n,Trbp);
	transpose(bp,Trbp,3,n);

	double Aa[9],Ba[9],TrAa[9];							// compute the matrix Aa and Ba in the paper
	//brmul(ap,Trbp,3,n,3,Aa);
	mult(ap,Trbp,Aa,3,n,3);

	for(i=0;i<9;i++)
	{
		Aa[i]=Aa[i]/((double)n);
	}
	//Trans(Aa,3,3,TrAa);
	//Minus(Aa,TrAa,3,3,Ba);

	double TV[3];
	TV[0] = Ba[5];
	TV[1] = Ba[6];
	TV[2] = Ba[1];
//////////////////////////////////////////////////////////////////////////			above is for computing Aa, Ba, delta(TV)
	double Ca[16],Tmp[9];
	double tr = Aa[0]+Aa[4]+Aa[8];								//trace of Aa
	Ca[0] = tr;
	for (i=0;i<3;i++)
	{
		Ca[4*(i+1)] = Ca[i+1] = TV[i];
	}
	for(i=0;i<9;i++)
	{
		Tmp[i] = Aa[i] + TrAa[i];
	}
	for(i=0;i<3;i++)
	{
		Tmp[4*i] = Tmp[4*i] - tr;
	}
	for (i=0;i<9;i++)
	{
		Ca[i+5+(int)(i/3)] = Tmp[i];
	}
//////////////////////////////////////////////////////////////////////////			above is for computing Ca
	int m,l=60;
	double q[4][4],qb[4],qc[4];
	double eps=0.000001;
	cstrq(&Ca[0],4,&q[0][0],qb,qc);
	m=csstq(4,qb,qc,&q[0][0],eps,l);
	double q0,q1,q2,q3;
	if(m>0)
	{
		if(qb[3]>0)
		{
			q0=q[0][3];
			q1=q[1][3];
			q2=q[2][3];
			q3=q[3][3];
		}
		else
		{
			q0=-q[0][3];
			q1=-q[1][3];
			q2=-q[2][3];
			q3=-q[3][3];
		}
	}
	else
	{
		return FALSE;
	}
	R[0]=pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2);
	R[1]=2*(q1*q2-q0*q3);
	R[2]=2*(q1*q3+q0*q2);
	R[3]=2*(q1*q2+q0*q3);
	R[4]=pow(q0,2)+pow(q2,2)-pow(q1,2)-pow(q3,2);
	R[5]=2*(q2*q3-q0*q1);
	R[6]=2*(q1*q3-q0*q2);
	R[7]=2*(q2*q3+q0*q1);
	R[8]=pow(q0,2)+pow(q3,2)-pow(q1,2)-pow(q2,2);

	double* tmpT = new double[3*n];
	double* Ta = new double[3*n];
	//brmul(R,a,3,3,n,tmpT);
	//Minus(b,tmpT,3,n,Ta);
	for (i=0;i<3;i++)
	{
		T[i] = 0.0;
		for(j=0;j<n;j++)
		{
			T[i] += Ta[i*n+j]/((double)n);
		}
	}
//////////////////////////////////////////////////////////////////////////			above is computing q0,q1,q2,q3 for attain R and T
	for (i=0;i<n;i++)
	{
		tmpT[i] = b[i] - T[0];
		tmpT[i+n] = b[i+n] - T[1];
		tmpT[i+2*n] = b[i+2*n] - T[2];
	}
	double RN[9];
	double err[3]={0.0};
	for (i=0;i<9;i++)
	{
		RN[i] = R[i];
	}
	//brinv(RN,3);
	//brmul(RN,tmpT,3,3,n,Ta);
	for (i=0;i<n;i++)
	{
		err[0] += (Ta[i] - a[i])*(Ta[i] - a[i])/((double)n);
		err[1] += (Ta[i+n] - a[i+n])*(Ta[i+n] - a[i+n])/((double)n);
		err[2] += (Ta[i+2*n] - a[i+2*n])*(Ta[i+2*n] - a[i+2*n])/((double)n);
	}
//////////////////////////////////////////////////////////////////////////			above is RN*(b-T) = a
	//brmul(R,a,3,3,n,tmpT);
	for (i=0;i<3*n;i++)
	{
		Ta[i] = tmpT[i]+T[(int)(i/n)];
	}

	for (i=0;i<n;i++)
	{
		err[0] += (Ta[i] - b[i])*(Ta[i] - b[i])/((double)n);
		err[1] += (Ta[i+n] - b[i+n])*(Ta[i+n] - b[i+n])/((double)n);
		err[2] += (Ta[i+2*n] - b[i+2*n])*(Ta[i+2*n] - b[i+2*n])/((double)n);
	}

//////////////////////////////////////////////////////////////////////////			compute error as R*a+T = b
	delete a;
	delete b;
	delete ap;
	delete bp;
	delete Trbp;
	delete Ta;
	delete tmpT;
	return TRUE;
}

double GetFoury_Error(Point3D *inp, Point3D *outp, int num, double *R, double *T)
{
	double terr=0.0;
	int i=0;
	for (i=0;i<num;i++)
	{
		double errin[3] = {0.0};
		errin[0] = inp[i].x*R[0]+inp[i].y*R[1]+inp[i].z*R[2]+T[0]-outp[i].x;
		errin[1] = inp[i].x*R[3]+inp[i].y*R[4]+inp[i].z*R[5]+T[1]-outp[i].y;
		errin[2] = inp[i].x*R[6]+inp[i].y*R[7]+inp[i].z*R[8]+T[2]-outp[i].z;
		terr += sqrt(errin[0]*errin[0]+errin[1]*errin[1]+errin[2]*errin[2]);
	}
	return terr;
}



void Compara(double *invect, double *outvect, double *R, int flatg)
{
	int i=0;
	//
	double vx = invect[0];
	double vy = invect[1];
	double vz = invect[2];
	
	if (flatg == 11)		// X轴或YZ平面
	{
		R[0] = vx;
		R[1] = vy;
		R[2] = vz;
		R[3] = -1.0*vy;
		R[4] = 1.0-vy*vy/(1.0+vx);
		R[5] = -1.0*vy*vz/(1.0+vx);
		R[6] = -1.0*vz;
		R[7] = -1.0*vy*vz/(1.0+vx);
		R[8] = 1.0-vz*vz/(1.0+vx);
	}
	else if (flatg == 12)		// Y轴或XZ平面
	{
		R[0] = 1.0-vx*vx/(1.0+vy);
		R[1] = -1.0*vx;
		R[2] = -1.0*vx*vz/(1.0+vy);
		R[3] = vx;
		R[4] = vy;
		R[5] = vz;
		R[6] = -1.0*vx*vz/(1.0+vy);
		R[7] = -1.0*vz;
		R[8] = 1.0-vz*vz/(1.0+vy);
	}
	else if (flatg == 13)		// Z轴或XY平面
	{
		R[0] = 1.0-vx*vx/(1.0+vz);
		R[1] = -1.0*vx*vy/(1.0+vz);
		R[2] = -1.0*vx;
		R[3] = -1.0*vx*vy/(1.0+vz);
		R[4] = 1.0-vy*vy/(1.0+vz);
		R[5] = -1.0*vy;
		R[6] = vx;
		R[7] = vy;
		R[8] = vz;
	}
	
}





bool ComparaObyO(Point3D inp, Point3D outp, int diver, int num, double *R, double *T)
{
	double tmpR[9]={0.0},tmpT[3]={0.0};
	tmpR[0] = tmpR[4] = tmpR[8] = 1.0;
	
	int i=0,j;
	
	double vc[3] = {inp.x,inp.y,inp.z};
	double vout[3];
	
	if (diver==10)		// 原点
	{
		for (j=0;j<9;j++)
		{
			R[j] = tmpR[j];
		}
		T[0] = outp.x-inp.x;
		T[1] = outp.y-inp.y;
		T[2] = outp.z-inp.z;
	}
	else if (diver==21 || diver==33)
	{
		if (vc[0]<0.0)
		{
			vc[0] = -1.0*vc[0];
			vc[1] = -1.0*vc[1];
			vc[2] = -1.0*vc[2];
		}
		Compara(vc,vout,R,11);
		T[0] = 0.0;
		T[1] = 0.0;
		T[2] = 0.0;
	}
	else if (diver==22 || diver==32)
	{
		if (vc[1]<0.0)
		{
			vc[0] = -1.0*vc[0];
			vc[1] = -1.0*vc[1];
			vc[2] = -1.0*vc[2];
		}
		Compara(vc,vout,R,12);
		T[0] = 0.0;
		T[1] = 0.0;
		T[2] = 0.0;
	}
	else if (diver==23 || diver==31)
	{
		if (vc[2]<0.0)
		{
			vc[0] = -1.0*vc[0];
			vc[1] = -1.0*vc[1];
			vc[2] = -1.0*vc[2];
		}
		Compara(vc,vout,R,13);
		T[0] = 0.0;
		T[1] = 0.0;
		T[2] = 0.0;
	}
	
	
	
	return 1;
}


bool TranRT(Point3D *inp, Point3D *outp, int num, double *R, double *T)
{
	// 先平移，后旋转
	int i=0;
	double pointin[3], pointout[3] ;
	for (i=0;i<num;i++)
	{
		// Pin = R*Pout + T ;
		// Pout = inv(R)*(Pin-T)
		pointin[0] = inp[i].x+T[0] ;
		pointin[1] = inp[i].y+T[1] ;
		pointin[2] = inp[i].z+T[2] ;
		
		
		// Pout = R*Pin + T
		//brmul(R,pointin,3,3,1,pointout);
		mult(R,pointin,pointout,3,3,1);
		
		
		outp[i].x = pointout[0];
		outp[i].y = pointout[1];
		outp[i].z = pointout[2];
	}
	
	return 1 ; 
}


bool TransferRT(double *inx, double *iny, double *inz, double *otx, double *oty, double *otz, int innum, double *R, double *T)
{
	// 先平移，后旋转
	int i=0;
	double pointin[3], pointout[3] ;
	for (i=0;i<innum;i++)
	{
		// Pin = R*Pout + T ;
		// Pout = inv(R)*(Pin-T)
		pointin[0] = inx[i]+T[0] ;
		pointin[1] = iny[i]+T[1] ;
		pointin[2] = inz[i]+T[2] ;
		
		
		// Pout = R*Pin + T
		//brmul(R,pointin,3,3,1,pointout);
		mult(R,pointin,pointout,3,3,1);
		
		
		otx[i] = pointout[0];
		oty[i] = pointout[1];
		otz[i] = pointout[2];
	}
	
	return 1 ; 
}




bool ComputeParameter(double *inx, double *iny, double *inz, int *pairarray, int *objectinfo, int innum, double *para, int &paranum)
{
	// para 前9个是旋转矩阵，后3个是平移矩阵
	double R[9]={0.0},T[3]={0.0};
	R[0] = R[4] = R[8] = 1.0;

	int i=0,j=0;
	memset(para,0,paranum);
	for (i=0;i<6;i++)
	{
		para[3+9*i] = 1.0;
		para[4+3+9*i] = 1.0;
		para[8+3+9*i] = 1.0;
	}


	// 建立索引表, 并对数据进行排序和分类
	Point3D *thispt = new Point3D[innum];
	Point3D *difipt = new Point3D[innum];
	int index[20]={0};			// 起始位置，个数
	int ttindex=0;

	for (i=0;i<innum;i++)			// 原点
	{
		if (pairarray[i]==11)
		{
			int tmpindex = index[1];
			if (tmpindex==0)
			{
				index[0] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[1]++;
			ttindex++;
		}
	}
	for (i=0;i<innum;i++)			// X轴
	{
		if (pairarray[i]==21)
		{
			int tmpindex = index[3];
			if (tmpindex==0)
			{
				index[2] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[3]++;
			ttindex++;
		}
	}

	for (i=0;i<innum;i++)			// Y轴
	{
		if (pairarray[i]==22)
		{
			int tmpindex = index[5];
			if (tmpindex==0)
			{
				index[4] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[5]++;
			ttindex++;
		}
	}
	for (i=0;i<innum;i++)			// Z轴
	{
		if (pairarray[i]==23)
		{
			int tmpindex = index[7];
			if (tmpindex==0)
			{
				index[6] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[7]++;
			ttindex++;
		}
	}
	for (i=0;i<innum;i++)			// XY平面
	{
		if (pairarray[i]==31)
		{
			int tmpindex = index[9];
			if (tmpindex==0)
			{
				index[8] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[9]++;
			ttindex++;
		}
	}
	for (i=0;i<innum;i++)			// XZ平面
	{
		if (pairarray[i]==32)
		{
			int tmpindex = index[11];
			if (tmpindex==0)
			{
				index[10] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[11]++;
			ttindex++;
		}
	}
	for (i=0;i<innum;i++)			// YZ平面
	{
		if (pairarray[i]==33)
		{
			int tmpindex = index[13];
			if (tmpindex==0)
			{
				index[12] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			index[13]++;
			ttindex++;
		}
	}

	for (i=0;i<innum;i++)			// 自定义 点间				注意是从0开始还是从1开始
	{
		if (pairarray[i]>100 
			&& pairarray[i]<101+pairarray[0])
		{
			int tmpindex = index[15];
			if (tmpindex==0)
			{
				index[14] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			int ono = pairarray[i]-100;
			difipt[ttindex].x = inx[ono];
			difipt[ttindex].y = iny[ono];
			difipt[ttindex].z = inz[ono];
			index[15]++;
			ttindex++;
		}
	}

	// 自定义消除
	for (i=0;i<innum;i++)			// 自定义 线间
	{
		if (pairarray[i]>100+objectinfo[0] 
			&& pairarray[i]<101+objectinfo[0]+objectinfo[1])
		{
			int tmpindex = index[17];
			if (tmpindex==0)
			{
				index[16] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			int ono = pairarray[i]-100;
			difipt[ttindex].x = inx[ono];
			difipt[ttindex].y = iny[ono];
			difipt[ttindex].z = inz[ono];
			index[17]++;
			ttindex++;
		}
	}

	for (i=0;i<innum;i++)			// 自定义 面间
	{
		if (pairarray[i]>100+objectinfo[0]+objectinfo[1]
			&& pairarray[i]<101+objectinfo[0]+objectinfo[1]+objectinfo[2])
		{
			int tmpindex = index[19];
			if (tmpindex==0)
			{
				index[18] = ttindex;
			}
			thispt[ttindex].x = inx[i];
			thispt[ttindex].y = iny[i];
			thispt[ttindex].z = inz[i];
			int ono = pairarray[i]-100;
			difipt[ttindex].x = inx[ono];
			difipt[ttindex].y = iny[ono];
			difipt[ttindex].z = inz[ono];
			index[19]++;
			ttindex++;
		}
	}

	for (i=0;i<innum;i++)
	{
		for (j=i;j<innum;j++)
		{
			if ((difipt[i].x == thispt[j].x) 
				&& (difipt[i].y == thispt[j].y)
				&& (difipt[i].z == thispt[j].z)
				&& (difipt[j].x == thispt[i].x)
				&& (difipt[j].y == thispt[i].y)
				&& (difipt[j].z == thispt[i].z)	)
			{
				difipt[j].x = difipt[j].y = difipt[j].z = 0.0;
				thispt[j].x = thispt[j].y = thispt[j].z = 0.0;
			}
		}
	}

	Point3D *pointinarray = new Point3D[innum];
	Point3D *pointotarray = new Point3D[innum];
	int pointnum = 0;
	Point3D *vectinarray = new Point3D[innum];
	Point3D *vectotarray = new Point3D[innum];
	int vectnum = 0;
	Point3D *pointinarrayfree = new Point3D[innum];
	Point3D *pointotarrayfree = new Point3D[innum];
	int pointnumfree = 0;
	Point3D *vectinarrayfree = new Point3D[innum];
	Point3D *vectotarrayfree = new Point3D[innum];
	int vectnumfree = 0;

	// 先求向量的，计算最小残差的RT
	for (i=index[0];i<index[0]+index[1];i++)		// 原点
	{
		pointinarray[pointnum].x = thispt[i].x;
		pointinarray[pointnum].y = thispt[i].y;
		pointinarray[pointnum].z = thispt[i].z;
		pointotarray[pointnum].x = 0.0;
		pointotarray[pointnum].y = 0.0;
		pointotarray[pointnum].z = 0.0;
		pointnum++;
	}
	for (i=index[2];i<index[2]+index[3];i++)		// X轴
	{
		vectinarray[vectnum].x = thispt[i].x;
		vectinarray[vectnum].y = thispt[i].y;
		vectinarray[vectnum].z = thispt[i].z;
		vectotarray[vectnum].x = 1.0;
		vectotarray[vectnum].y = 0.0;
		vectotarray[vectnum].z = 0.0;
		vectnum++;
	}
	for (i=index[4];i<index[4]+index[5];i++)		// Y轴
	{
		vectinarray[vectnum].x = thispt[i].x;
		vectinarray[vectnum].y = thispt[i].y;
		vectinarray[vectnum].z = thispt[i].z;
		vectotarray[vectnum].x = 0.0;
		vectotarray[vectnum].y = 1.0;
		vectotarray[vectnum].z = 0.0;
		vectnum++;
	}
	for (i=index[6];i<index[6]+index[7];i++)		// Z轴
	{
		vectinarray[vectnum].x = thispt[i].x;
		vectinarray[vectnum].y = thispt[i].y;
		vectinarray[vectnum].z = thispt[i].z;
		vectotarray[vectnum].x = 0.0;
		vectotarray[vectnum].y = 0.0;
		vectotarray[vectnum].z = 1.0;
		vectnum++;
	}
	for (i=index[8];i<index[8]+index[9];i++)		// XY
	{
		vectinarray[vectnum].x = thispt[i].x;
		vectinarray[vectnum].y = thispt[i].y;
		vectinarray[vectnum].z = thispt[i].z;
		vectotarray[vectnum].x = 0.0;
		vectotarray[vectnum].y = 0.0;
		vectotarray[vectnum].z = 1.0;
		vectnum++;
	}
	for (i=index[10];i<index[10]+index[11];i++)		// XZ
	{
		vectinarray[vectnum].x = thispt[i].x;
		vectinarray[vectnum].y = thispt[i].y;
		vectinarray[vectnum].z = thispt[i].z;
		vectotarray[vectnum].x = 0.0;
		vectotarray[vectnum].y = 0.0;
		vectotarray[vectnum].z = 1.0;
		vectnum++;
	}
	for (i=index[12];i<index[12]+index[13];i++)		// YZ
	{
		vectinarray[vectnum].x = thispt[i].x;
		vectinarray[vectnum].y = thispt[i].y;
		vectinarray[vectnum].z = thispt[i].z;
		vectotarray[vectnum].x = 0.0;
		vectotarray[vectnum].y = 0.0;
		vectotarray[vectnum].z = 1.0;
		vectnum++;
	}
	for (i=index[14];i<index[14]+index[15];i++)
	{
		if (thispt[i].x!=difipt[i].x 
			|| thispt[i].y!=difipt[i].y 
			|| thispt[i].z!=difipt[i].z)
		{
			pointinarrayfree[pointnumfree].x = thispt[i].x;
			pointinarrayfree[pointnumfree].y = thispt[i].y;
			pointinarrayfree[pointnumfree].z = thispt[i].z;
			pointotarrayfree[pointnumfree].x = difipt[i].x;
			pointotarrayfree[pointnumfree].y = difipt[i].y;
			pointotarrayfree[pointnumfree].z = difipt[i].z;
			pointnumfree++;
		}
	}
	for (i=index[16];i<index[16]+index[17];i++)
	{
		if (thispt[i].x!=difipt[i].x 
			|| thispt[i].y!=difipt[i].y 
			|| thispt[i].z!=difipt[i].z)
		{
			vectinarrayfree[vectnumfree].x = thispt[i].x;
			vectinarrayfree[vectnumfree].y = thispt[i].y;
			vectinarrayfree[vectnumfree].z = thispt[i].z;
			vectotarrayfree[vectnumfree].x = difipt[i].x;
			vectotarrayfree[vectnumfree].y = difipt[i].y;
			vectotarrayfree[vectnumfree].z = difipt[i].z;
			vectnumfree++;
		}
	}
	for (i=index[18];i<index[18]+index[19];i++)
	{
		if (thispt[i].x!=difipt[i].x 
			|| thispt[i].y!=difipt[i].y 
			|| thispt[i].z!=difipt[i].z)
		{
			vectinarrayfree[vectnumfree].x = thispt[i].x;
			vectinarrayfree[vectnumfree].y = thispt[i].y;
			vectinarrayfree[vectnumfree].z = thispt[i].z;
			vectotarrayfree[vectnumfree].x = difipt[i].x;
			vectotarrayfree[vectnumfree].y = difipt[i].y;
			vectotarrayfree[vectnumfree].z = difipt[i].z;
			vectnumfree++;
		}
	}

	double tR[9]={0.0},tT[3]={0.0};
	tR[0] = tR[4] = tR[8] = 1.0;

	// Function_Foury
	// 这里，完全基准物体按照次序求算，自由物体间按四元数方法
	Point3D *tmpinv = new Point3D[vectnumfree];
	Point3D *tmpotv = new Point3D[vectnumfree];
	Point3D *mininv = new Point3D[vectnumfree];
	Point3D *minotv = new Point3D[vectnumfree];

	int classnum = (int)pow(2.0,vectnumfree);
	int *classarray = new int[classnum*vectnumfree];		// 存各级的排列方式,0或1

	if (vectnumfree>2)
	{
		for (i=0;i<classnum;i++)
		{
			for (j=0;j<vectnumfree;j++)
			{
				int tp = (int)(pow(2.0,j+1));
				if (i%tp>=pow(2.0,j))
				{
					classarray[i*vectnumfree+j] = 1;
				}
				else
				{
					classarray[i*vectnumfree+j] = -1;
				}
			}
		}
		
		double tErr = 1000000.0;
		for (i=0;i<classnum;i++)			// 向量有正反，取最小的值
		{
			for (j=0;j<vectnumfree;j++)
			{
				tmpinv[j].x = classarray[i*vectnumfree+j]*vectinarrayfree[j].x;
				tmpinv[j].y = classarray[i*vectnumfree+j]*vectinarrayfree[j].y;
				tmpinv[j].z = classarray[i*vectnumfree+j]*vectinarrayfree[j].z;
				tmpotv[j].x = vectotarrayfree[j].x;
				tmpotv[j].y = vectotarrayfree[j].y;
				tmpotv[j].z = vectotarrayfree[j].z;
			}
			Function_Foury(tmpinv,tmpotv,vectnumfree,tR,tT);
			double thisE = GetFoury_Error(tmpinv,tmpotv,vectnumfree,tR,tT);
			if (tErr > thisE)
			{
				tErr = thisE;
				for (j=0;j<vectnumfree;j++)
				{
					mininv[j].x = tmpinv[j].x;
					mininv[j].y = tmpinv[j].y;
					mininv[j].z = tmpinv[j].z;
					minotv[j].x = tmpotv[j].x;
					minotv[j].y = tmpotv[j].y;
					minotv[j].z = tmpotv[j].z;
				}
			}
		}
		
		Function_Foury(mininv,minotv,vectnumfree,tR,tT);
		double tE = GetFoury_Error(mininv,minotv,vectnumfree,tR,tT);
		

	}


	TranRT(pointinarrayfree,pointinarrayfree,pointnumfree,tR,tT);
	TranRT(pointinarray,pointinarray,pointnum,tR,tT);
	TranRT(vectinarray,vectinarray,vectnum,tR,tT);

	double nR[9] = {0.0}, nT[3] = {0.0};
	nR[0] = nR[4] = nR[8] = 1.0;
	// 点与点之间，不存在向量的方向性
	if (pointnumfree>2)
	{
		Function_Foury(pointinarrayfree,pointotarrayfree,pointnumfree,nR,nT);
		double tpE = GetFoury_Error(pointinarrayfree,pointotarrayfree,pointnumfree,nR,nT);
	}


	TranRT(pointinarray,pointinarray,pointnum,nR,nT);
	TranRT(vectinarray,vectinarray,vectnum,nR,nT);

	double pR[9] = {0.0}, pT[3] = {0.0};
	//brmul(nR,tR,3,3,3,pR);
	//brmul(nR,tT,3,3,1,pT);
	mult(nR,tR,pR,3,3,3);
	mult(nR,tT,pT,3,3,1);
	pT[0] += nT[0];
	pT[1] += nT[1];
	pT[2] += nT[2];

	// 向量相
	for (i=0;i<vectnum;i++)
	{
		if (vectotarray[i].x==1.0 && vectotarray[i].y==0.0 && vectotarray[i].z==0.0)
		{
			ComparaObyO(vectinarray[i],vectotarray[i],21,1,nR,nT);
		}
		else if (vectotarray[i].x==0.0 && vectotarray[i].y==1.0 && vectotarray[i].z==0.0)
		{
			ComparaObyO(vectinarray[i],vectotarray[i],22,1,nR,nT);
		}
		else if (vectotarray[i].x==0.0 && vectotarray[i].y==0.0 && vectotarray[i].z==1.0)
		{
			ComparaObyO(vectinarray[i],vectotarray[i],23,1,nR,nT);
		}
		TranRT(vectinarray,vectinarray,vectnum,nR,nT);
		for (j=0;j<9;j++)
		{
			tR[j] = pR[j];
		}
		for (j=0;j<3;j++)
		{
			tT[j] = pT[j];
		}
		//brmul(nR,tR,3,3,3,pR);
		//brmul(nR,tT,3,3,1,pT);
		mult(nR,tR,pR,3,3,3);
		mult(nR,tT,pT,3,3,1);
		pT[0] += nT[0];
		pT[1] += nT[1];
		pT[2] += nT[2];
	}



	nT[0] = nT[1] = nT[2] = 0.0;
	// 原点相关，只取第一个
	ComparaObyO(pointinarray[0],pointotarray[0],10,1,nR,nT);
	pT[0] += nT[0];
	pT[1] += nT[1];
	pT[2] += nT[2];


	for (i=0;i<9;i++)
	{
		para[i] = pR[i];
	}
	for (i=0;i<3;i++)
	{
		para[i+9] = pT[i];
	}

	Point3D *tinpt = new Point3D[innum];
	Point3D *totpt = new Point3D[innum];
	int ptcnt=0;
	for (i=0;i<vectnumfree;i++)
	{
		tinpt[ptcnt].x = vectinarrayfree[i].x;
		tinpt[ptcnt].y = vectinarrayfree[i].y;
		tinpt[ptcnt].z = vectinarrayfree[i].z;
		totpt[ptcnt].x = vectotarrayfree[i].x;
		totpt[ptcnt].y = vectotarrayfree[i].y;
		totpt[ptcnt].z = vectotarrayfree[i].z;
		ptcnt++;
	}
	for (i=0;i<pointnumfree;i++)
	{
		tinpt[ptcnt].x = pointinarrayfree[i].x;
		tinpt[ptcnt].y = pointinarrayfree[i].y;
		tinpt[ptcnt].z = pointinarrayfree[i].z;
		totpt[ptcnt].x = pointotarrayfree[i].x;
		totpt[ptcnt].y = pointotarrayfree[i].y;
		totpt[ptcnt].z = pointotarrayfree[i].z;
		ptcnt++;
	}
	for (i=0;i<vectnum;i++)
	{
		tinpt[ptcnt].x = vectinarray[i].x;
		tinpt[ptcnt].y = vectinarray[i].y;
		tinpt[ptcnt].z = vectinarray[i].z;
		totpt[ptcnt].x = vectotarray[i].x;
		totpt[ptcnt].y = vectotarray[i].y;
		totpt[ptcnt].z = vectotarray[i].z;
		ptcnt++;
	}
	for (i=0;i<pointnum;i++)
	{
		tinpt[ptcnt].x = pointinarray[i].x;
		tinpt[ptcnt].y = pointinarray[i].y;
		tinpt[ptcnt].z = pointinarray[i].z;
		totpt[ptcnt].x = pointotarray[i].x;
		totpt[ptcnt].y = pointotarray[i].y;
		totpt[ptcnt].z = pointotarray[i].z;
		ptcnt++;
	}

	if (ptcnt>2)
	{
		Function_Foury(pointinarrayfree,pointotarrayfree,pointnumfree,nR,nT);
		double tpE = GetFoury_Error(pointinarrayfree,pointotarrayfree,pointnumfree,nR,nT);
	}
	else
	{
		nR[0] = nR[4] = nR[8] = 1.0;
		nR[1] = nR[2] = nR[3] = nR[5] = nR[6] = nR[7] = 0.0;
		nT[0] = nT[1] = nT[2] = 0.0;
	}

	delete tmpinv;
	delete tmpotv;
	delete mininv;
	delete minotv;
	delete classarray;

	if (pointnum<1 && pointnumfree<3 && vectnum<1 && vectnumfree<3)
	{
		return false;
	}


	delete pointinarray;
	delete pointotarray;
	delete vectinarray;
	delete vectotarray;
	delete difipt;
	delete thispt;



	return true;
}