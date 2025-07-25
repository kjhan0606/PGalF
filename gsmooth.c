#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<srfftw.h>


void gaussian_Smoothing(float *denGrid,int nx,int ny,int nz, 
		float cellsize, float RG){
	long i,j,k,mx,mmx;
	mx = 2*(nx/2+1);
	mmx = (nx/2+1);

	long ncells =  mx*ny*nz;

	rfftwnd_plan p, pinv;
	fftw_real scale;
	double NORM;

	scale = 1.L/(ncells);

	p = rfftw3d_create_plan(nx,ny,nz,FFTW_REAL_TO_COMPLEX,
			FFTW_ESTIMATE | FFTW_IN_PLACE);
	pinv = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL,
                                FFTW_ESTIMATE);

	fftw_real *a = (fftw_real*)denGrid;
	fftw_complex *A = (fftw_complex*)a;

	rfftwnd_one_real_to_complex(p,a,NULL);

	double expfac = -(M_PI*RG)*(M_PI*RG);
	int nz2 = nz/2;
	int ny2 = ny/2;
	int nx2 = nx/2;
	for(k=0;k<nz;k++){
		long k1=k;
		if(k > nz2) k1 = nz-k;
		for(j=0;j<ny;j++){
			long j1 = j;
			if(j > ny2) j1 = ny-j;
			for(i=0;i<nx/2;i++){
				long i1 = i;
				long kk = i1*i1 + j1*j1 + k1*k1;
				double frac = exp(kk*expfac);
				double factor = frac*scale;
				A[i+mmx*(j+ny*k)].re *= factor;
				A[i+mmx*(j+ny*k)].im *= factor;
			}
		}
	}
	rfftwnd_one_complex_to_real(pinv,A, NULL);


	rfftwnd_destroy_plan(p);
	rfftwnd_destroy_plan(pinv);
}
