#define BOUNDITER 4
#define COREDENRESOLUTION (1.e-3L) // fractional error with respect to the core density 
#define PEAKTHRESHOLD 1.E11L   // lowest peak stellar density in unit of Msun/ckpc^3 

/* please tune these five parameters for higher resolution simulations */
#define MINCORENMEM 100 /* the minimum number of star/dm particles to make a core */
#define NUMNEIGHBOR 15 /* the number of neighbors to build the neighbor network */
#define NSHELLDIVIDE 10 /* the number of division of non-core particles */
#define MERGINGPEAKLENGTH 4.e-3  /* the separation limit of peaks to merge in cMpc/h */
#define MINCORESTARMASS -1 /* the minimum stellar mass for a core */
#define MINSTELLARMASS 1.e7  /* the minimun stellar mass of the FoF halo for galaxy finding with stellar density */

#define MAXNUMCORE 1000000
//#define NUMNEARDEN 15
#define NUMNEARDEN 10

#define NUMSTELLARNEIGHBORS 30 // This is for the core detection which help 
							//	to find density peaks only using star particles.
