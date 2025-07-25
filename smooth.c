/*
 * pgcc -o smooth smooth.c -lm -L/user/kjhan/fftwfinal/lib -I/user/kjhan/fftwfinal/include -lsrfftw -lsfftw -lsfftw
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <srfftw.h>
#include "hst.h"
#define EDGE 64
int N,M,N0,M0;
fftw_real *a,*b,*c;
fftw_complex *A, *B, *C;


long long mx,mmx;

void main(argc,argv)
int argc;
char **argv;
{
     rfftwnd_plan p, pinv;
     fftw_real scale;
     float x;
     double dx, dy, rfil, rfil2, rtrunc=64, rtrunc2, dr2;
     double rfil3;
     double NORM, hstnorm, hstsum;
     long  i, j, i0, j0, offset;
     int nx, ny;
     int gaussflag;
     FILE *hstfile;



     gaussflag = 0;
     if( argc == 2 ) {
         rfil = atof(argv[1]);
         rfil2 = rfil*rfil;
         gaussflag = 1;
         if( rfil < 0 ) gaussflag = 2;
     }
    fprintf(stderr,"gaussflag %d\n",gaussflag);

     fread(&nx,sizeof(int),1,stdin); fread(&ny,sizeof(int),1,stdin);
     N0 = nx;
     M0 = ny;
     N = (N0+2*EDGE);
     M = (M0+2*EDGE);
     mx = 2*(N/2+1);
     mmx = (N/2+1);

     rtrunc2 = rtrunc*rtrunc;
     scale = 1.0 / (M * N);
     p    = rfftw2d_create_plan(M, N, FFTW_REAL_TO_COMPLEX,
                                FFTW_ESTIMATE | FFTW_IN_PLACE);
     pinv = rfftw2d_create_plan(M, N, FFTW_COMPLEX_TO_REAL,
                                FFTW_ESTIMATE);

     a = (fftw_real*)malloc(sizeof(fftw_real)*M*(2*(N/2+1)));
     b = (fftw_real*)malloc(sizeof(fftw_real)*M*(2*(N/2+1)));
     c = (fftw_real*)malloc(sizeof(fftw_real)*M*N);
     C = (fftw_complex*)malloc(sizeof(fftw_complex)*M*(N/2+1));

     /* aliases for accessing complex transform outputs: */
     A = (fftw_complex*) a;
     B = (fftw_complex*) b;

     for(i=0; i<M; i++) {
         for(j=0; j<N; j++) {
            a[i*mx+j] = 0.0;
            b[i*mx+j] = 0.0;
        }
    }

     if( gaussflag == 1 ) {
         nx = (int) (24*rfil); ny = (int) (24*rfil);
        nx = M; ny = N;
     }
     else {
         hstfile = fopen("hstpsf0.dat","r");
         fread(&nx,sizeof(int),1,hstfile); fread(&ny,sizeof(int),1,hstfile);
     }
     NORM = 1.0/(2.0*M_PI*rfil2);
     hstsum = 0.0;
     for(i=0; i<nx; i++) {
         for(j=0; j<ny; j++) {
            float hst;

            offset = j + i*nx;
            i0 = (i - nx/2 + M) % M;
            j0 = (j - ny/2 + N) % N;
            dx = (double) (i - nx/2);
            dy = (double) (j - ny/2);
            dr2 = dx*dx + dy*dy;
            if( gaussflag == 1 ) {
                b[i0*(long)mx+j0] = NORM*exp(-0.5*dr2/rfil2);
            }
            else {
                fread(&hst,sizeof(float),1,hstfile);
                b[i0*(long)mx+j0] = hst;
                if( gaussflag == 2 ) {
                    b[i0*(long)mx+j0] *= NORM*exp(-0.5*dr2/rfil2);

                }
            }
         }
     }
     for (i = 0; i < M0; ++i) {
        int i0;

        i0 = i+EDGE;
         for (j = 0; j < N0; ++j) {
            int j0;

            j0 = j + EDGE;
            fread(&x,sizeof(float),1,stdin);
            if(!isnan(x)) a[i0*(long )mx+j0] = (fftw_real) x;
         }
     }

     rfftwnd_one_real_to_complex(p, a, NULL);
     rfftwnd_one_real_to_complex(p, b, NULL);

     for (i = 0; i < M; ++i)
          for (j = 0; j < N/2+1; ++j) {
               long long ij = i*(long)(N/2+1) + j;
               C[i*mmx+j].re = (A[ij].re * B[ij].re
                             - A[ij].im * B[ij].im) * scale;
               C[i*mmx+j].im = (A[ij].re * B[ij].im
                             + A[ij].im * B[ij].re) * scale;
          }


     /* enhance small scale density contrast: 0.8pixel */
     rfil3 = 1.0*1.0;
     NORM = 1.0/(2.0*M_PI*rfil3);
     hstsum = 0.0;
     for(i=0; i<nx; i++) {
         for(j=0; j<ny; j++) {
            float hst;

            offset = j + i*nx;
            i0 = (i - nx/2 + M) % M;
            j0 = (j - ny/2 + N) % N;
            dx = (double) (i - nx/2);
            dy = (double) (j - ny/2);
            dr2 = dx*dx + dy*dy;
            if( gaussflag == 1 ) {
                b[i0*(long)mx+j0] = NORM*exp(-0.5*dr2/rfil3);
            }
            else {
                fread(&hst,sizeof(float),1,hstfile);
                b[i0*(long)mx+j0] = hst;
                if( gaussflag == 2 ) {
                    b[i0*(long)mx+j0] *= NORM*exp(-0.5*dr2/rfil3);
                }
            }
         }
     }
     rfftwnd_one_real_to_complex(p, b, NULL);
     for (i = 0; i < M; ++i)
          for (j = 0; j < N/2+1; ++j) {
               long long ij = i*(long)(N/2+1) + j;
               C[i*mmx+j].re += (A[ij].re * B[ij].re
                             - A[ij].im * B[ij].im) * scale;
               C[i*mmx+j].im += (A[ij].re * B[ij].im
                             + A[ij].im * B[ij].re) * scale;
          }



     /* inverse transform to get c, the convolution of a and b;
        this has the side effect of overwriting C */
     rfftwnd_one_complex_to_real(pinv, C, c);

     fwrite(&N0,sizeof(int),1,stdout); fwrite(&M0,sizeof(int),1,stdout);
     for (i = EDGE; i < M0+EDGE; ++i) {
         for (j = EDGE; j < N0+EDGE; ++j) {
            x = (float) c[i*N+j];
         fwrite(&x,sizeof(float),1,stdout);
         }
     }
     rfftwnd_destroy_plan(p);
     rfftwnd_destroy_plan(pinv);
}
