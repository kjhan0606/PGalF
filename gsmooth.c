#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<srfftw.h>
#include<sfftw.h>
#include "ramses.h"


void gaussian_Smoothing(float *denGrid,int nx,int ny,int nz, 
		double cellsize, float RG){
	long i,j,k,mx,mxh;
	mx = 2*(nx/2+1);
	mxh = (nx/2+1);

	long ncells =  mx*ny*nz;

	rfftwnd_plan p, pinv;
	double scale;
	double NORM;

	scale = 1.L/((double)nx*(double)ny*(double)nz);

	p = rfftw3d_create_plan(nz,ny,nx,FFTW_REAL_TO_COMPLEX,
			FFTW_ESTIMATE | FFTW_IN_PLACE);
	pinv = rfftw3d_create_plan(nz, ny, nx, FFTW_COMPLEX_TO_REAL,
                                FFTW_ESTIMATE | FFTW_IN_PLACE);

	/*
	fftw_real *a = (fftw_real*)denGrid;
	*/
	fftw_real *a = (fftw_real*)fftw_malloc(sizeof(fftw_real)*ncells);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ncells;i++) a[i] = denGrid[i];
	fftw_complex *A = (fftw_complex*)a;

	rfftwnd_one_real_to_complex(p,a,NULL);

	RG = RG/cellsize;

	double expfac = -(M_PI*RG)*(M_PI*RG);
	int nzh = nz/2;
	int nyh = ny/2;
	int nxh = nx/2;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<nz;k++){
		long k1=k;
		if(k > nzh) k1 = k-nz;
		double zfact = exp(k1*k1*expfac/(double)(nz*nz));
		for(j=0;j<ny;j++){
			long j1 = j;
			if(j > nyh) j1 = j-ny;
			double yfact = exp(j1*j1*expfac/(double)(ny*ny));
			for(i=0;i<=nx/2;i++){
				long i2 = i*2;
				double xfact = exp(i*i*expfac/(double)(nx*nx));
				double factor = xfact*yfact*zfact;
				a[i2+mx*(j+ny*k)] *= factor;
				a[i2+1+mx*(j+ny*k)] *= factor;
			}
		}
	}

	rfftwnd_one_complex_to_real(pinv,A, NULL);


#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ncells;i++) denGrid[i] = a[i]*scale;
	fftw_free(a);
	rfftwnd_destroy_plan(p);
	rfftwnd_destroy_plan(pinv);
	DEBUGPRINT0("Now exit the Gaussian Smoothing Job\n");
}
