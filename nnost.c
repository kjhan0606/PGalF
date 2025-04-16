/******************************************************************************/
/******************************************************************************/
/*******                NEAREST NEIGHBOR SEARCHING ALGORITHM         **********/
/******************************************************************************/
/******************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<unistd.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#include "Memory.h"

#include "ramses.h"
#include "tree.h"
#include "params.h"


#define nodesize dist


int nullfct0();
int nullfct1();

float W4(float ,float );

#ifndef _OPENMP
#define Omp_get_thread_num() nullfct0()
#define Omp_get_num_threads() nullfct1()
#else
#define Omp_get_thread_num() omp_get_thread_num()
#define Omp_get_num_threads() omp_get_num_threads()
#endif


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/*
#define YES 1
#define NO 0
*/


#define SubCellDicision(a,b) ((a)>(b)? 1:0)
#define DivideNode(ThisNode,nparticles) ((ThisNode->nodesize > 0.5*MINCELLWIDTH ? YES:NO) && (nparticles>=MIN_CELL_PARTICLE_NUM ? YES:NO))

TStruct *calThreadFreeNodeStart(size_t navail, size_t twork, size_t mys,size_t myf, 
		TStruct *work, TStruct *nextFreeNode){
	TStruct *threadFreeNode, *tmp;
	size_t i, tnp;
	tnp = 0;
	tmp = work;
	for(i=0;i<twork;i++) {
		tnp += tmp->Nparticle;
		tmp ++;
	}
	float ratio = (float)navail/(float)tnp;
	size_t rnp = 0;
	tmp = work;
	for(i=0;i<mys;i++){
		rnp += tmp->Nparticle;
		tmp++;
	}
	threadFreeNode = nextFreeNode + (size_t)(rnp*ratio);
	return threadFreeNode;
}


TStruct *divide_node_Near(TStruct *motherNode, TStruct *nextFreeNode, int recursiveflag)
{
	TStruct *p2tree, tmpnode[8];
	void *lastOffspring;
	TPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
	size_t i, j, k, mnode, mx,my,mz;
	float tmpx,tmpy,tmpz,tmpdist2, distmax;
	float ptlmass;

	nodeparticles = motherNode->daughter;
	motherNode->type = TYPE_TREE;
	motherNode->mass = 0;
	motherNode->monox = motherNode->monoy = motherNode->monoz = 0;

	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		ptlmass = p2ptl->mass;
		motherNode->mass += ptlmass;
		motherNode->monox += ptlmass * p2ptl->x;
		motherNode->monoy += ptlmass * p2ptl->y;
		motherNode->monoz += ptlmass * p2ptl->z;
	}
	motherNode->monox /= motherNode->mass;
	motherNode->monoy /= motherNode->mass;
	motherNode->monoz /= motherNode->mass;
	distmax = -1.E20;
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		tmpx = p2ptl->x - motherNode->monox;
		tmpy = p2ptl->y - motherNode->monoy;
		tmpz = p2ptl->z - motherNode->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax, tmpdist2);
	}
	motherNode->nodesize = sqrt(distmax);

	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
	}
	for(p2ptl=nodeparticles;p2ptl;) {
		mx = SubCellDicision(p2ptl->x,motherNode->monox);
		my = SubCellDicision(p2ptl->y,motherNode->monoy);
		mz = SubCellDicision(p2ptl->z,motherNode->monoz);
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++;
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	for(i=0;i<8;i++) if(tmpnode[i].Nparticle >0) break;

	TStruct *firstDaughter,*newDaughter;
	firstDaughter = newDaughter =  nextFreeNode;
	// Link mother node to the first daughter or particle (first offsprig) 
	if( divideThisNode(motherNode, tmpnode[i].Nparticle) )
		motherNode->daughter = (void*)firstDaughter;
	else 
		motherNode->daughter = (void*) tmpnode[i].daughter;

	// make links between the current offspring
	// lastOffsping: the terminal point of the current offsprings
	lastOffspring = NULL;
	for(i=0;i<8;i++) {
		if( divideThisNode(motherNode, tmpnode[i].Nparticle) ){
			newDaughter->daughter = tmpnode[i].daughter;
			newDaughter->Nparticle = tmpnode[i].Nparticle;
			if(lastOffspring) ((GENERAL_TPtl_POINTER*) lastOffspring)->sibling = newDaughter;
			lastOffspring = (void *)newDaughter; // update lastOffspring
			newDaughter++;
		}
		else if(tmpnode[i].Nparticle >0) {
			tmpptr = tmpnode[i].daughter;
			if(lastOffspring) ((GENERAL_TPtl_POINTER*)lastOffspring)->sibling = tmpptr;
			for(;tmpptr;tmpptr=tmpptr->sibling) lastOffspring = (void*) tmpptr;
		}
	}
	// close by setting the sibling of the last offspring to the mother's sibling.
	((GENERAL_TPtl_POINTER *)lastOffspring)->sibling = motherNode->sibling;
	nextFreeNode = newDaughter;


	if(recursiveflag == RECURSIVE){
		TStruct *NextJobNode = firstDaughter;
		for(i=0;i<8;i++){
			if(divideThisNode(motherNode, tmpnode[i].Nparticle) ){
				nextFreeNode = divide_node_Near(NextJobNode,nextFreeNode, recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return nextFreeNode;
}

void Make_Tree_Near(
		TStruct *TREE_START, 
		size_t nnode, 
		TPtlStruct *ptl, 
		size_t np, 
		int recursiveflag){
	size_t i;
	TStruct *nextFreeNode = TREE_START+1;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + i +1;
		ptl[i].type = TYPE_PTL;
	}
	ptl[np-1].sibling = NULL;
	TREE_START->daughter = &(ptl[0]);
	TREE_START->Nparticle = np;
	if(recursiveflag == RECURSIVE) nextFreeNode = divide_node_Near(TREE_START, nextFreeNode, recursiveflag);
	else if(recursiveflag == SERIALIZED){
		TStruct *work;
		for(work=TREE_START;nextFreeNode-work >0; work++){
			nextFreeNode = divide_node_Near(work, nextFreeNode, SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{
		TStruct *work = TREE_START;
		do{
			nextFreeNode = divide_node_Near(work, nextFreeNode, PTHREAD);
			work ++;
		}
		while( work < nextFreeNode &&  (nextFreeNode-work) < 64);
		size_t twork = (nextFreeNode-work);
		if(twork <=0) return;

#ifdef _OPENMP
#pragma omp parallel 
#endif
		{
			int pid = Omp_get_thread_num();
			int npid = Omp_get_num_threads();
			size_t mys, myf;
			size_t worksize = (twork + npid - 1)/npid;
			size_t j;
			mys = worksize *pid;
			myf = MIN(worksize *(pid+1), twork);
			TStruct *threadnextFreeNode = nextFreeNode + 
				pid * ( (nnode - (nextFreeNode-TREE_START)+npid-1)/ npid);
			size_t navail =  nnode - (nextFreeNode-TREE_START);
			threadnextFreeNode = calThreadFreeNodeStart(navail, twork, 
					mys,myf, work, nextFreeNode);
			for(j=mys;j<myf;j++) {
				threadnextFreeNode = divide_node_Near(work+j,threadnextFreeNode,RECURSIVE);
			}
		}
	}
}



int near_open(particle *point, TStruct *tree, int npneigh, float maxdist , int Num_neighbor){
	dptype tmpx,tmpy,tmpz,dist2, r2, dist, r;
	if(npneigh >= Num_neighbor) {
		tmpx = point->x - tree->monox;
		tmpy = point->y - tree->monoy;
		tmpz = point->z - tree->monoz;
		r = tree->nodesize;
		dist = sqrtf(tmpx*tmpx+ tmpy*tmpy + tmpz*tmpz);
		if(dist-r > maxdist) return NO;
		else return YES;
	}
	else return YES;
}


typedef struct nearestneighbor{
	dptype dist2;
	float mass;
	int indx;
} Neighbor;



#define INSERT(dist2, maxdist, maxdist2, npneigh, neighbor) do{\
	for(i=0;i<npneigh;i++) if(neighbor[i].dist2>dist2) break;\
	for(j=npneigh-1;j>=i;j--) neighbor[j+1] = neighbor[j];\
	neighbor[i].dist2 = dist2;\
	neighbor[i].mass = ( (TPtlStruct*)ptr)->mass;\
	neighbor[i].indx = ( (TPtlStruct*)ptr)-ptl;\
	npneigh ++;\
	npneigh=MIN(Num_neighbor,npneigh);\
	maxdist2 = (neighbor[npneigh-1].dist2);\
	maxdist = sqrtf(maxdist2);\
}while(0)


int Find_Near(
		particle *point, 
		int Num_neighbor, 
		TStruct *tree, 
		TPtlStruct *ptl, 
		float *maxr, 
		int *nearindx,
		dptype *DIST2,
		float *MASS
		){
	int i,j,k;
	dptype dist2, tmpx,tmpy,tmpz;
	dptype maxdist, maxdist2;
	void *ptr;
	int npneigh=0;
	Neighbor neighbor[MAX_NUM_NEAR]={0};
	maxdist = maxdist2 = 1.E23;
	ptr = (void*)tree;
	while(ptr != NULL){
		switch( ((TYPE*)ptr)->type) {
			case TYPE_TREE:
				switch( near_open(point, ptr, npneigh, maxdist, Num_neighbor)){
					case YES:
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default:
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default:
				tmpx = point->x - ((TPtlStruct *)ptr)->x;
				tmpy = point->y - ((TPtlStruct *)ptr)->y;
				tmpz = point->z - ((TPtlStruct *)ptr)->z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				if(npneigh < Num_neighbor || dist2<maxdist2){
					INSERT(dist2, maxdist, maxdist2, npneigh, neighbor);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	*maxr = maxdist;
	for(i=0;i<npneigh;i++) {
		nearindx[i] = neighbor[i].indx;
		DIST2[i] = neighbor[i].dist2;
		MASS[i] = neighbor[i].mass;
	}
	return npneigh;
}
int nearConstOpen(particle *point, TStruct *tree, dptype constR){
	dptype tmpx = point->x - tree->monox;
	dptype tmpy = point->y - tree->monoy; 
	dptype tmpz = point->z - tree->monoz;
	dptype r = tree->nodesize;
    dptype dist = sqrtf(tmpx*tmpx+ tmpy*tmpy + tmpz*tmpz); 
	if(dist-r > constR) return NO; 
	else return YES;
}
dptype getDenConstR (
        particle *point,
        TStruct *tree,
        TPtlStruct *ptl,
		dptype constR
        ){
	dptype tmpx,tmpy,tmpz,dist2,mass;
    dptype den = 0;
    void *ptr = (void*)tree;
    while(ptr != NULL){
        switch( ((TYPE*)ptr)->type) {
            case TYPE_TREE:
                switch( nearConstOpen(point, ptr, constR)){
                    case YES:
                        ptr = (void *)(((TStruct*)ptr)->daughter);
                        break;
                    default:
                        ptr = (void *)(((TStruct*)ptr)->sibling);
                }
                break;
            default:
                tmpx = point->x - ((TPtlStruct *)ptr)->x;
                tmpy = point->y - ((TPtlStruct *)ptr)->y;
                tmpz = point->z - ((TPtlStruct *)ptr)->z;
                dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
                mass = ((TPtlStruct *)ptr)->mass;
                if(dist2<constR){
                    den += W4(dist2,constR/2.) * mass;
                }
                ptr = (void *) (((TPtlStruct*)ptr)->sibling);
        }
    }
    return den;
}

int Find_Near2(
		particle *point, 
		int Num_neighbor, 
		TStruct *tree, 
		TPtlStruct *ptl, 
		size_t nptl, 
		float *maxr, 
		int *nearindx,
		dptype *DIST2,
		float *MASS
		){
	size_t i,j,k;
	size_t iptl = 0;
	dptype dist2, tmpx,tmpy,tmpz;
	dptype maxdist, maxdist2;
	void *ptr;
	int npneigh=0;
	Neighbor neighbor[nptl];
	maxdist = maxdist2 = 1.E23;
	ptr = (void*)tree;
	while(ptr != NULL){
		switch( ((TYPE*)ptr)->type) {
			case TYPE_TREE:
				switch( near_open(point, ptr, npneigh, maxdist, Num_neighbor)){
					case YES:
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default:
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default:
				tmpx = point->x - ((TPtlStruct *)ptr)->x;
				tmpy = point->y - ((TPtlStruct *)ptr)->y;
				tmpz = point->z - ((TPtlStruct *)ptr)->z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				if(iptl < Num_neighbor || dist2<maxdist2){
					neighbor[iptl].dist2 = dist2;
					neighbor[iptl].mass = ptl->mass;
					neighbor[iptl].indx = ptl->indx;
					iptl ++;
					if(iptl>Num_neighbor) maxdist2 = MIN(maxdist2, dist2);
					else maxdist2 = MAX(maxdist2, dist2);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	npneigh = 0;
	for(i=0;i<iptl;i++){
		if(neighbor[i].dist2 <= maxdist2) {
			nearindx[npneigh] = neighbor[i].indx;
			DIST2[npneigh] = neighbor[i].dist2;
			MASS[npneigh] = neighbor[i].mass;
			npneigh ++;
		}
	}
	return npneigh;
}

double f(double);
#define AMPW41 (-0.2387324)
#define AMPW42 (0.07957747)
float W4(float dist,float h){
    float twomr,r;
    r = dist/h;
    if(r<1.){
        /*
        return 1/PI*(1.-1.5*dist*dist+0.75*dist*dist*dist);
        */
        return AMPW41*((2.-r)*r*r-1.333333)/(h*h*h);
    }
    else if(r < 2.){
        twomr = (2.-r)/h;
        /*
        return 1/PI*0.25*twomr*twomr*twomr;
        */
        return AMPW42*twomr*twomr*twomr;
    }
    else {
        return 0.;
    }
}
#undef AMPW41
#undef AMPW42


void findsphdensity(SimpleBasicParticleType *bp,int np,int *nearindex, int Numnear,
        float *densph){
    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    dptype *dist2;
	float *mass;
    float fplmf,ptlmass;
    size_t i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    particle *p;
	p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }
    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    for(i=0;i<np;i++){
        ptl[i].type = TYPE_PTL;
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
        ptl[i].sibling = ptl+(i+1);
        ptl[i].mass = bp[i].mass;
    }
    ptl[np-1].sibling = NULL;
	size_t nnode = MAX(np/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));
	int recursiveflag;
	if(np > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;

	recursiveflag = RECURSIVE;
	Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
    float mindist;
    mindist = 2.e25;
#pragma omp parallel for private(i,j,k) schedule(guided)
    for(i=0;i<np;i++){
        int tmpindx[Numnear];
        dptype dist2[Numnear];
        float mass[Numnear];
        int res;
        float neardist;
        k = i*Numnear;
        res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2, mass);
        densph[i] = 0;
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
            dist2[j] = sqrtf(dist2[j]);
            densph[i] += W4(dist2[j],neardist/2.);
        }
    }
    Free(TREE);Free(ptl);Free(p);
    /*
    printf("ins=", tmpindx);
    fprintf(stdout,"Terminated with success\n");fflush(stdout);
    */
#ifdef DEBUG
    printf("exiting findsphdensity\n");
#endif
}

void starfindsphdensity(SimpleBasicParticleType *bp,int np,int *nearindex, int Numnear,
        float *densph){
    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    dptype *dist2;
	float *mass;
    float fplmf,ptlmass;
    size_t i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    int mp;
    particle *p;


    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }

    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    mp = 0;
    for(i=0;i<np;i++){
        if(bp[i].type == TYPE_STAR){
            ptl[mp].type = TYPE_PTL;
            ptl[mp].x = bp[i].x;
            ptl[mp].y = bp[i].y;
            ptl[mp].z = bp[i].z;
            ptl[mp].sibling = ptl+(mp+1);
            ptl[mp].mass = bp[i].mass;
            mp++;
        }
    }
    ptl[mp-1].sibling = NULL;
	size_t nnode = MAX(mp/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));

	int recursiveflag;
	if(np > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;

	recursiveflag = RECURSIVE;
    Make_Tree_Near(TREE,nnode, ptl,mp,recursiveflag);
    DEBUGPRINT("after Star_Make_Tree_near for nstar = %d\n",mp);
    float mindist;
    mindist = 2.e25;
//    int NumNearDen = NUMNEARDEN;
#pragma omp parallel for private(i,j,k) schedule(dynamic)
    for(i=0;i<np;i++){
        int res;
        int tmpindx[NUMNEARDEN];
        dptype tmpd2[NUMNEARDEN];
        float tmpmass[NUMNEARDEN];
        float neardist;
        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
        densph[i] = 0;
        for(j=0;j<res;j++){
            tmpd2[j] = sqrtf(tmpd2[j]);
            densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
        }
    }
    DEBUGPRINT("End of allocating the stellar density for np = %d\n",np);

	// now find the k-nearest neighbors.
    for(i=0;i<np;i++){
        ptl[i].type = TYPE_PTL;
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
        ptl[i].sibling = ptl+i+1;
        ptl[i].mass = bp[i].mass;
    }
    ptl[np-1].sibling = NULL;
	recursiveflag = RECURSIVE;
	nnode = MAX(np/2, 65*10000);
	TREE = Realloc(TREE, sizeof(TStruct)*nnode);
    Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
    DEBUGPRINT("before Make_Tree_near2 with np= %d\n", np);
#pragma omp parallel for private(i,j,k) schedule(dynamic)
    for(i=0;i<np;i++){
        int res;
        int tmpindx[NUMNEIGHBOR];
        dptype dist2[NUMNEIGHBOR];
        float mass[NUMNEIGHBOR];
        float neardist;
        k = i*Numnear;
        res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2, mass);
        if(res != Numnear){
            printf("error occurred %d %d : %d\n", res, Numnear, np);
            exit(99);
        }
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
        }
    }
//    Free(mass);Free(dist2); Free(h);
    Free(TREE);Free(ptl);Free(p);
    DEBUGPRINT0("end of starfinddensity()\n");
}

void findStellarCore(
		SimpleBasicParticleType *bp,int np,
		int Numnear,
        float *densph,
		Coretype **Core,
		int *NumCore,  // return value for number of core
		int maxnumcore,
		int **Nearindex){
	Coretype *core = *Core;
	int numcore;
//    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    float fplmf,ptlmass;
    size_t i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    int mp;
	int nstar;
    particle *p;


    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }

    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    nstar = 0;
    for(i=0;i<np;i++){
        if(bp[i].type == TYPE_STAR){
            ptl[nstar].type = TYPE_PTL;
            ptl[nstar].x = bp[i].x;
            ptl[nstar].y = bp[i].y;
            ptl[nstar].z = bp[i].z;
            ptl[nstar].sibling = ptl+(nstar+1);
            ptl[nstar].mass = bp[i].mass;
            ptl[nstar].indx = i;
            nstar++;
        }
    }
    ptl[nstar-1].sibling = NULL;
	size_t nnode = MAX(nstar/2,65*10000);
    DEBUGPRINT("nstar = %d, nnode= %d\n",nstar, nnode);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));


	int recursiveflag;
	if(nnode > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;
	recursiveflag = RECURSIVE;
    Make_Tree_Near(TREE,nnode, ptl,nstar,recursiveflag);
    DEBUGPRINT("after Star_Make_Tree_near for nstar = %d\n",nstar);
//  all to star
#pragma omp parallel for private(i,j) schedule(guided)
    for(i=0;i<np;i++){
        int res;
        int tmpindx[NUMNEARDEN];
        dptype tmpd2[NUMNEARDEN];
        float tmpmass[NUMNEARDEN];
        float neardist;
        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
#if defined mADV   || defined ADV
		if(neardist < MIN_CONST_R_SMOOTHING){
			dptype constR = MIN_CONST_R_SMOOTHING;
			densph[i] = getDenConstR(p+i, TREE,ptl, constR);
		}
		else {
       		densph[i] = 0;
	        for(j=0;j<res;j++){ 
				tmpd2[j] = sqrtf(tmpd2[j]); 
				densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
			}
		}
#else
        densph[i] = 0;
        for(j=0;j<res;j++){
            tmpd2[j] = sqrtf(tmpd2[j]);
            densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
        }
#endif
    }




	// now find the k-nearest neighbors.
	mp = 0;
    for(i=0;i<np;i++){
		if(bp[i].type == TYPE_STAR){
			p[mp].x = bp[i].x;
			p[mp].y = bp[i].y;
			p[mp].z = bp[i].z;
			p[mp].indx = i;
			mp++;
		}
    }
	int numSNeigh = MIN(nstar, NUMSTELLARNEIGHBORS);
	int *nearindex = *Nearindex = (int*)Malloc(sizeof(int)*numSNeigh*nstar,PPTR(*Nearindex));
//  star to star
#pragma omp parallel for private(i,j,k) schedule(guided)
    for(i=0;i<nstar;i++){
        int res;
        int tmpindx[NUMSTELLARNEIGHBORS];
        float neardist;
		dptype dist2[NUMSTELLARNEIGHBORS];
		float mass[NUMSTELLARNEIGHBORS];
        k = i*numSNeigh;
        res = Find_Near(p+i,numSNeigh,TREE,ptl,&neardist,tmpindx,
				dist2, mass);
        if(res != numSNeigh){
            DEBUGPRINT("error occurred %d %d : %d\n", res, numSNeigh, np);
            exit(99);
        }
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
        }
    }

//	Free(h);
    Free(TREE);
    DEBUGPRINT0("end of stellarneighbor()\n");
	int istar = 0;
	numcore  = 0;
	for(i=0;i<np;i++){
		if(bp[i].type == TYPE_STAR && densph[i]>PEAKTHRESHOLD){
			k = istar*numSNeigh;
			int iflag = 0;
			for(j=0;j<numSNeigh;j++){
				if(densph[i] < densph[ptl[nearindex[k+j]].indx]) {
					iflag = 1;
					break;
				}
			}
			if(iflag ==0){
				int ibp = p[istar].indx;
				core[numcore].peak = ibp;
				core[numcore].cx = bp[ibp].x;
				core[numcore].cy = bp[ibp].y;
				core[numcore].cz = bp[ibp].z;
				core[numcore].density = densph[i];
				numcore ++;
				if(numcore > maxnumcore-10){
					maxnumcore += MAXNUMCORE;
					*Core = Realloc(*Core, sizeof(Coretype)*maxnumcore);
					core = *Core;
				}
			}
		}
		if(bp[i].type == TYPE_STAR) istar ++;
	}
	Free(ptl); Free(p);
	DEBUGPRINT("The number of cores : %d and before MergingPeak\n", numcore);

	int MergingPeak(SimpleBasicParticleType *, int, Coretype *, int, int);
	if(numcore >10) numcore = MergingPeak(bp,np, core, numcore,0);
	core = *Core = Realloc(*Core, sizeof(Coretype)*numcore);

	DEBUGPRINT("The number of cores : %d and after MergingPeak\n", numcore);

	nearindex = *Nearindex = Realloc(*Nearindex, sizeof(int)*np*Numnear);
    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
#pragma omp parallel for 
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }
    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
#pragma omp parallel for 
    for(i=0;i<np;i++){ 
		ptl[i].type = TYPE_PTL; 
		ptl[i].x = bp[i].x; 
		ptl[i].y = bp[i].y; 
		ptl[i].z = bp[i].z; 
		ptl[i].sibling = ptl+(i+1); 
		ptl[i].mass = bp[i].mass; 
    }
    ptl[np-1].sibling = NULL;
	nnode = MAX(np/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));

	if(nnode > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;
//  all to all
	recursiveflag = RECURSIVE;
    Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
    DEBUGPRINT("after All_Make_Tree_near for np = %d\n",np);
#pragma omp parallel for private(i,j,k) schedule(guided)
    for(i=0;i<np;i++){
//		DEBUGPRINT("p%d  is now being processed: %d\n", i, recursiveflag);
        int res;
        int tmpindx[NUMNEIGHBOR];
        float mass[NUMNEIGHBOR];
        dptype dist2[NUMNEIGHBOR];
        float neardist;
        k = i*Numnear;
        res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,
				dist2, mass);
        if(res != Numnear){
            DEBUGPRINT("error occurred %d %d : %d\n", res, Numnear, np);
            exit(99);
        }
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
        }
    }

	Free(TREE);
	Free(ptl);
	Free(p);


	*NumCore = numcore;

	return;
}

