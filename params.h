//*******************
// The number of interation to determine 
// the boundedness of member particles
#define BOUNDITER 4
//*******************
//


//*******************
// fractional error to measure the core density
#define COREDENRESOLUTION (1.e2L) 


//*******************
// lowest stellar-density peak
// in unit of h^2 Msun/ckpc^3
#define PEAKTHRESHOLD 2000.L
//*******************

//--------------------------------
// please tune these six parameters for higher resolution simulations 
#define MINCORENMEM 100 // the minimum number of star/dm particles to make a core 
#define NUMNEIGHBOR 8 // the number of neighbors to build the neighbor network 
#define NSHELLDIVIDE 10 // the number of division of non-core particles 
#define MERGINGPEAKLENGTH 4.e-3  // the separation limit of peaks to merge in cMpc/h 
#define MINCORESTARMASS -1 // the minimum stellar mass for a core 
#define MINSTELLARMASS 1.e7  // the minimun stellar mass of the FoF halo for galaxy finding with stellar density 
//--------------------------------

//*******************
// Maximum number of cores 
#define MAXNUMCORE 1000000
//*******************


//*******************
// (obsolete) The number of nearby stars 
// to measure stellar density
#define NUMNEARDEN 10
//*******************

//*******************
// (obsolete) The number of stellar neighbors
// for the core detection to find density peaks 
#define NUMSTELLARNEIGHBORS 30 
//*******************


//*******************
// (obsolete) minimum smoothing length 
// for stellar density in cMpc/h
#define MIN_CONST_R_SMOOTHING 0.005 
//*******************

//*******************
// The cellsize for TSC of the stellar density 
// in unit of cMpc/h
#define TSC_CELL_SIZE 0.005 
//*******************

//*******************
// Gaussian Smoothing Length 
// for stellar density in unit of cMpc/h. 
// This should be not smaller 
// than 2*TSC_CELL_SIZE.
// Note that the smoothing may 
// lower the stellar density 
// and PEAKTHRESHOLD should be 
// lowered accordingly.
#define Gaussian_Smoothing_Length 0.010 
//*******************
