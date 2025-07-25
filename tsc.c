#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>

#include"ramses.h"
#include"tree.h"
#include"params.h"

#define icopysign(a) copysignf(1.,a)

#define den(i,j,k) (denGrid[ (i) + nx *( (j)+ ny*(k))]);



void assign_density_TSC(SimpleBasicParticleType *bp, int np, float *denGrid, 
		int nx, int ny, int nz){
	int i,j,k;
	for(i=0;i<np;i++){
		float wx1wy1,wx2wy1,wx3wy1;
        float wx1wy2,wx2wy2,wx3wy2;
        float wx1wy3,wx2wy3,wx3wy3;
		float pmas0 = bp[i].mass;
		float p05 = pmas0 *0.5;
		float xp,yp,zp;
		x = bp[i].x;
		y = bp[i].y;
		z = bp[i].z;
		int nearx,neary,nearz;
		nearx = xp;
		neary = yp;
		nearz = zp;
		float xmin = xp - nearx;
		float ymin = yp - neary;
		float zmin = zp - nearz;
		float xsign = icopysign(xmin);
		float ysign = icopysign(xmin);
		float zsign = icopysign(xmin);
		int i1,i2,i3,j1,j2,j3,k1,k2,k3;
		i1 = nxp[nearx];
        i2 = nxp[nearx+xsign];
        i3 = nxp[nearx-xsign];

        j1 = nyp[neary];
        j2 = nyp[neary+ysign];
        j3 = nyp[neary-ysign];

        k1 = nzp[nearz];
        k2 = nzp[nearz+zsign];
        k3 = nzp[nearz-zsign];

        xd1 = fabsf(xmin);
        yd1 = fabsf(ymin);
        zd1 = fabsf(zmin);
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
