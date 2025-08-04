#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>

#include "ramses.h"
#include "tree.h"
#include "Memory.h"
#include "params.h"

#define icopysign(a) copysignf(1.,a)

#define den(i,j,k) (den[(i)+mx*(long)((j)+ny*(k))])
#ifdef XYZDBL
#define FABS(a) fabs(a)
#define RINT(a) rint(a)
#else
#define FABS(a) fabsf(a)
#define RINT(a) rintf(a)
#endif

#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif


void findDen(SimpleBasicParticleType *bp,int np, float *densph, double xmin, double ymin, 
		double zmin, double cellsize,
		float *den,
		int nx, int ny, int nz
		){
	long mx = 2*(nx/2+1);
	long i,j,k;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0;i<np;i++){
		double xp = (bp[i].x - xmin)/cellsize;
		double yp = (bp[i].y - ymin)/cellsize;
		double zp = (bp[i].z - zmin)/cellsize;
		int nearx = RINT(xp);
		int neary = RINT(yp);
		int nearz = RINT(zp);
		float xmin = xp - nearx;
		float ymin = yp - neary;
		float zmin = zp - nearz;
		long xsign = icopysign(xmin);
		long ysign = icopysign(ymin);
		long zsign = icopysign(zmin);
		long i1 = nearx;
		long i2 = nearx+xsign;
		long i3 = nearx-xsign;
		long j1 = neary;
		long j2 = neary+ysign;
		long j3 = neary-ysign;
		long k1 = nearz;
		long k2 = nearz+zsign;
		long k3 = nearz-zsign;
		float xd1 = FABS(xmin);
		float yd1 = FABS(ymin);
		float zd1 = FABS(zmin);
		float wx1 = (0.75-xd1*xd1);
		float wy1 = (0.75-yd1*yd1);
		float wz1 = (0.75-zd1*zd1);
		float wx3 = 0.5*(0.25+xd1*(xd1-1.));
		float wx2 = wx3 + xd1;
		float wy3 = 0.5*(0.25+yd1*(yd1-1.));
		float wy2 = wy3 + yd1;
		float wz3 = 0.5*(0.25+zd1*(zd1-1.));
		float wz2 = wz3 + zd1;

		float wx1wy1 = wx1*wy1;
		float wx2wy1 = wx2*wy1; 
		float wx3wy1 = wx3*wy1;
		float wx1wy2 = wx1*wy2;
		float wx2wy2 = wx2*wy2;
		float wx3wy2 = wx3*wy2;
		float wx1wy3 = wx1*wy3;
		float wx2wy3 = wx2*wy3;
		float wx3wy3 = wx3*wy3;

		float density = 0;
		density += den(i1,j1,k1)*wx1wy1*wz1;
		density += den(i2,j1,k1)*wx2wy1*wz1; 
		density += den(i3,j1,k1)* wx3wy1*wz1;
        density += den(i1,j2,k1)* wx1wy2*wz1;
        density += den(i2,j2,k1)* wx2wy2*wz1;
        density += den(i3,j2,k1)* wx3wy2*wz1;
        density += den(i1,j3,k1)* wx1wy3*wz1;
        density += den(i2,j3,k1)* wx2wy3*wz1;
        density += den(i3,j3,k1)* wx3wy3*wz1;

        density += den(i1,j1,k2)* wx1wy1*wz2;
        density += den(i2,j1,k2) *  wx2wy1*wz2;
        density += den(i3,j1,k2) *  wx3wy1*wz2;
        density += den(i1,j2,k2) *  wx1wy2*wz2;
        density += den(i2,j2,k2) *  wx2wy2*wz2;
        density += den(i3,j2,k2) *  wx3wy2*wz2;
        density += den(i1,j3,k2) *  wx1wy3*wz2;
        density += den(i2,j3,k2) *  wx2wy3*wz2;
        density += den(i3,j3,k2) *  wx3wy3*wz2;

        density += den(i1,j1,k3) *  wx1wy1*wz3;
        density += den(i2,j1,k3) *  wx2wy1*wz3;
        density += den(i3,j1,k3) *  wx3wy1*wz3;
        density += den(i1,j2,k3) *  wx1wy2*wz3;
        density += den(i2,j2,k3) *  wx2wy2*wz3;
        density += den(i3,j2,k3) *  wx3wy2*wz3;
        density += den(i1,j3,k3) *  wx1wy3*wz3;
        density += den(i2,j3,k3) *  wx2wy3*wz3;
        density += den(i3,j3,k3) *  wx3wy3*wz3;


		densph[i] = density;

	}
	DEBUGPRINT0("Now exit the density-interpolation Job\n");
}


void assign_density_TSC(SimpleBasicParticleType *bp, int np, float *den, 
		int nx, int ny, int nz,
		double Xmin,
		double Ymin,
		double Zmin,
		double cellsize){
	int nid,myid=0;
	int NZWIDTH=4;
	double pmas0,p05;
//	long long i,j,k;

//	double rngc;
	int usndsize,dgetsize,dsndsize,ugetsize,dest,src;
	int ierror,stag,rtag,tag;
	int local_z_start,local_ny_after_transpose,local_y_start_after_transpose,
		total_local_size;
	long long mx,mgrid,ngrid,slicesize;
	long long *nxperiodic,*nyperiodic,*nzperiodic;
	long long *nxp,*nyp,*nzp;

	double volfact = (0.001/TSC_CELL_SIZE)*(0.001/TSC_CELL_SIZE)*(0.001/TSC_CELL_SIZE);

//	rngc = (double) nx *(double)ny *(double)nz;

//	if(myid==0) printf("Now in the OMP tsc: %d %d %d %g\n", local_nz,nstart_z,np, pmas);

	/*
	nxperiodic = (long long *)Malloc(sizeof(long long)*(nx+10),PPTR(nxperiodic)); nxp = nxperiodic+5;
	nyperiodic = (long long *)Malloc(sizeof(long long)*(ny+10),PPTR(nyperiodic)); nyp = nyperiodic+5;
	nzperiodic = (long long *)Malloc(sizeof(long long)*(nz_per+10),PPTR(nzperiodic)); nzp = nzperiodic+5;
	*/
	{
		long long i;
		/*
		for(i=-5;i<nx+5;i++) nxp[i] = (i+nx)%nx;
		for(i=-5;i<ny+5;i++) nyp[i] = (i+ny)%ny;
		for(i=-5;i<nz_per+5;i++) nzp[i] = (i+nz_per)%nz_per;
		*/


		mx = 2*(nx/2+1);
		slicesize = mx * ny;
		ngrid = slicesize * nz;

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i=0;i<ngrid;i++) den[i] = 0;

	}


	int nxyz = ( (nx>ny) ? nx: ny);
	nxyz =  ( (nxyz > nz) ? nxyz: nz);
	int nthreads, mthreads;
	double hwidth,xwidth;
#ifdef _OPENMP
#pragma omp parallel shared(nthreads)
#endif
	{
#ifdef _OPENMP
#pragma omp master
#endif
		{
			nthreads = omp_get_num_threads();
		}
	}
	xwidth = (double)nxyz/(double)nthreads;
	hwidth = 0.5*xwidth;
	if(hwidth < 5) {
		mthreads = (double)nxyz/(2*5.);
	}
	else {
		mthreads = nthreads;
	}

//	if(myid==0) printf("Now in the OMP tsc with thread= %d\n", mthreads);
	signed char *target = (signed char*)Malloc(sizeof(short)*np, PPTR(target));
#ifdef _OPENMP
	{
#ifdef XYZDBL
		double Xwidth = (double)nxyz/(double)mthreads;
#else
		double Xwidth = (double)nxyz/(double)mthreads;
#endif
		long long i;
#pragma omp parallel for num_threads(mthreads)
		for(i=0;i<np;i++){
			if(bp[i].type != TYPE_STAR) continue;
#ifdef XYZDBL
			double xp;
#else
			double xp;
#endif
			if(nxyz == nx) {
				xp = (bp[i].x-Xmin)/cellsize;
			}
			else if(nxyz == ny) {
				xp = (bp[i].y-Ymin)/cellsize;
			}
			else {
				xp = (bp[i].z-Zmin)/cellsize;
			}

			int j = (int)(xp/Xwidth) + 1;
			int k = (int)(xp/Xwidth+0.5L) + 1;
			if(j >= mthreads+1){
				target[i] = mthreads;
			}
			else if(k ==j) {
				target[i] = -j;
			}
			else {
				target[i] = j;
			}
		}
	}
#endif

#if defined(MULTIMASS)
#ifdef _OPENMP
#pragma omp parallel private(pmas0,p05) num_threads(mthreads)
#endif
#else
#ifdef _OPENMP
#pragma omp parallel num_threads(mthreads)
#endif
#endif
	{
		int idthread = omp_get_thread_num() + 1;
		long long i;
		for(i=0;i<np;i++){
			if(bp[i].type != TYPE_STAR) continue;
			if(target[i] != -idthread) continue;
#ifdef XYZDBL
			double xp,yp,zp;
#else
			float xp,yp,zp;
#endif
			float wx1wy1,wx2wy1,wx3wy1;
			float wx1wy2,wx2wy2,wx3wy2;
			float wx1wy3,wx2wy3,wx3wy3;
			float xmin,ymin,zmin;
			long long xsign,ysign,zsign;
			long long nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
			float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
			float xd1,yd1,zd1;
			pmas0= bp[i].mass;
			p05 = pmas0 * 0.5;
			xp = (bp[i].x-Xmin)/cellsize;
			yp = (bp[i].y-Ymin)/cellsize;
			zp = (bp[i].z-Zmin)/cellsize;
			nearx = RINT(xp);
			neary = RINT(yp);
			nearz = RINT(zp);
			xmin = xp - nearx;
			ymin = yp - neary;
			zmin = zp - nearz;
			xsign = icopysign(xmin);
			ysign = icopysign(ymin);
			zsign = icopysign(zmin);
	
			/*
			i1 = nxp[nearx];
			i2 = nxp[nearx+xsign];
			i3 = nxp[nearx-xsign];
			*/
			i1 = nearx;
			i2 = nearx+xsign;
			i3 = nearx-xsign;
			j1 = neary;
			j2 = neary+ysign;
			j3 = neary-ysign;
			k1 = nearz;
			k2 = nearz+zsign;
			k3 = nearz-zsign;
			xd1 = FABS(xmin);
			yd1 = FABS(ymin);
			zd1 = FABS(zmin);
	
		
			wx1 = (0.75-xd1*xd1)*pmas0;
			wy1 =  0.75-yd1*yd1;
			wz1 =  0.75-zd1*zd1;
			wx3 = p05*(0.25+xd1*(xd1-1.));
			wx2 = wx3 + pmas0*xd1;
			wy3 = 0.5*(0.25+yd1*(yd1-1.));
			wy2 = wy3 + yd1;
			wz3 = 0.5*(0.25+zd1*(zd1-1.));
			wz2 = wz3 + zd1;
		
			wx1wy1 = wx1*wy1;
			wx2wy1 = wx2*wy1;
			wx3wy1 = wx3*wy1;
			wx1wy2 = wx1*wy2;
			wx2wy2 = wx2*wy2;
			wx3wy2 = wx3*wy2;
			wx1wy3 = wx1*wy3;
			wx2wy3 = wx2*wy3;
			wx3wy3 = wx3*wy3; 
			den(i1,j1,k1) += wx1wy1*wz1; 
			den(i2,j1,k1) += wx2wy1*wz1; 
			den(i3,j1,k1) += wx3wy1*wz1; 
			den(i1,j2,k1) += wx1wy2*wz1; 
			den(i2,j2,k1) += wx2wy2*wz1; 
			den(i3,j2,k1) += wx3wy2*wz1; 
			den(i1,j3,k1) += wx1wy3*wz1; 
			den(i2,j3,k1) += wx2wy3*wz1; 
			den(i3,j3,k1) += wx3wy3*wz1;
	
		    den(i1,j1,k2) += wx1wy1*wz2; 
			den(i2,j1,k2) += wx2wy1*wz2; 
			den(i3,j1,k2) += wx3wy1*wz2; 
			den(i1,j2,k2) += wx1wy2*wz2; 
			den(i2,j2,k2) += wx2wy2*wz2; 
			den(i3,j2,k2) += wx3wy2*wz2; 
			den(i1,j3,k2) += wx1wy3*wz2; 
			den(i2,j3,k2) += wx2wy3*wz2; 
			den(i3,j3,k2) += wx3wy3*wz2;
	
		    den(i1,j1,k3) += wx1wy1*wz3; 
			den(i2,j1,k3) += wx2wy1*wz3; 
			den(i3,j1,k3) += wx3wy1*wz3; 
			den(i1,j2,k3) += wx1wy2*wz3; 
			den(i2,j2,k3) += wx2wy2*wz3; 
			den(i3,j2,k3) += wx3wy2*wz3; 
			den(i1,j3,k3) += wx1wy3*wz3; 
			den(i2,j3,k3) += wx2wy3*wz3; 
			den(i3,j3,k3) += wx3wy3*wz3;
	
		}
	}
#if defined(MULTIMASS)
#ifdef _OPENMP
#pragma omp parallel private(pmas0,p05) num_threads(mthreads)
#endif
#else
#ifdef _OPENMP
#pragma omp parallel num_threads(mthreads)
#endif
#endif
	{
		int idthread = omp_get_thread_num()+1;
		long long i;
		for(i=0;i<np;i++){
			if(target[i] != idthread) continue;
#ifdef XYZDBL
			double xp,yp,zp;
#else
			float xp,yp,zp;
#endif
			float wx1wy1,wx2wy1,wx3wy1;
			float wx1wy2,wx2wy2,wx3wy2;
			float wx1wy3,wx2wy3,wx3wy3;
			float xmin,ymin,zmin;
			long long xsign,ysign,zsign;
			long long nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
			float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
			float xd1,yd1,zd1;
			pmas0= bp[i].mass;
			p05 = pmas0 * 0.5;
			xp = (bp[i].x-Xmin)/cellsize;
			yp = (bp[i].y-Ymin)/cellsize;
			zp = (bp[i].z-Zmin)/cellsize;
			nearx = RINT(xp);
			neary = RINT(yp);
			nearz = RINT(zp);
			xmin = xp - nearx;
			ymin = yp - neary;
			zmin = zp - nearz;
			xsign = icopysign(xmin);
			ysign = icopysign(ymin);
			zsign = icopysign(zmin);
	
			i1 = nearx;
			i2 = nearx+xsign;
			i3 = nearx-xsign;
			j1 = neary;
			j2 = neary+ysign;
			j3 = neary-ysign;
			k1 = nearz;
			k2 = nearz+zsign;
			k3 = nearz-zsign;
			xd1 = FABS(xmin);
			yd1 = FABS(ymin);
			zd1 = FABS(zmin);
	
		
			wx1 = (0.75-xd1*xd1)*pmas0;
			wy1 =  0.75-yd1*yd1;
			wz1 =  0.75-zd1*zd1;
			wx3 = p05*(0.25+xd1*(xd1-1.));
			wx2 = wx3 + pmas0*xd1;
			wy3 = 0.5*(0.25+yd1*(yd1-1.));
			wy2 = wy3 + yd1;
			wz3 = 0.5*(0.25+zd1*(zd1-1.));
			wz2 = wz3 + zd1;
		
			wx1wy1 = wx1*wy1;
			wx2wy1 = wx2*wy1;
			wx3wy1 = wx3*wy1;
			wx1wy2 = wx1*wy2;
			wx2wy2 = wx2*wy2;
			wx3wy2 = wx3*wy2;
			wx1wy3 = wx1*wy3;
			wx2wy3 = wx2*wy3;
			wx3wy3 = wx3*wy3; 
			den(i1,j1,k1) += wx1wy1*wz1; 
			den(i2,j1,k1) += wx2wy1*wz1; 
			den(i3,j1,k1) += wx3wy1*wz1; 
			den(i1,j2,k1) += wx1wy2*wz1; 
			den(i2,j2,k1) += wx2wy2*wz1; 
			den(i3,j2,k1) += wx3wy2*wz1; 
			den(i1,j3,k1) += wx1wy3*wz1; 
			den(i2,j3,k1) += wx2wy3*wz1; 
			den(i3,j3,k1) += wx3wy3*wz1;
	
		    den(i1,j1,k2) += wx1wy1*wz2; 
			den(i2,j1,k2) += wx2wy1*wz2; 
			den(i3,j1,k2) += wx3wy1*wz2; 
			den(i1,j2,k2) += wx1wy2*wz2; 
			den(i2,j2,k2) += wx2wy2*wz2; 
			den(i3,j2,k2) += wx3wy2*wz2; 
			den(i1,j3,k2) += wx1wy3*wz2; 
			den(i2,j3,k2) += wx2wy3*wz2; 
			den(i3,j3,k2) += wx3wy3*wz2;
	
		    den(i1,j1,k3) += wx1wy1*wz3; 
			den(i2,j1,k3) += wx2wy1*wz3; 
			den(i3,j1,k3) += wx3wy1*wz3; 
			den(i1,j2,k3) += wx1wy2*wz3; 
			den(i2,j2,k3) += wx2wy2*wz3; 
			den(i3,j2,k3) += wx3wy2*wz3; 
			den(i1,j3,k3) += wx1wy3*wz3; 
			den(i2,j3,k3) += wx2wy3*wz3; 
			den(i3,j3,k3) += wx3wy3*wz3;
	
		}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(long i=0;i<ngrid;i++) den[i] *= volfact; // return to density in cubic kpc

	DEBUGPRINT0("Now exit the main tsc\n");





	Free(target);
}
