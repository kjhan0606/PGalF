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
#define OFFSET 0.001
Box findbox(TPtlStruct *ptl, int nend){
    float xmin,ymin,zmin,xmax,ymax,zmax;
    float width;
    Box box;
    int i;
    box.x = box.y = box.z = 0.;
    width=xmax=ymax=zmax = -1.E25;
    xmin=ymin=zmin = 1.E25;
    for(i=0;i<nend;i++){
        xmin = MIN(xmin,ptl[i].x);
        ymin = MIN(ymin,ptl[i].y);
        zmin = MIN(zmin,ptl[i].z);
        xmax = MAX(xmax,ptl[i].x);
        ymax = MAX(ymax,ptl[i].y);
        zmax = MAX(zmax,ptl[i].z);
    }
    width = MAX(width,xmax-xmin);
    width = MAX(width,ymax-ymin);
    width = MAX(width,zmax-zmin);
    box.x = xmin;
    box.y = ymin;
    box.z = zmin;
    box.width = width;
    box.x = box.x -box.width*OFFSET;
    box.y = box.y -box.width*OFFSET;
    box.z = box.z -box.width*OFFSET;
    box.width = box.width + box.width*OFFSET*2.;
    return box;
}



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


TStruct *new_divide_node_Near(TStruct *motherNode, TStruct *nextFreeNode, int recursiveflag)
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
				nextFreeNode = new_divide_node_Near(NextJobNode,nextFreeNode, recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return nextFreeNode;
}

void new_Make_Tree_Near(
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
	if(recursiveflag == RECURSIVE) nextFreeNode = new_divide_node_Near(TREE_START, nextFreeNode, recursiveflag);
	else if(recursiveflag == SERIALIZED){
		TStruct *work;
		for(work=TREE_START;nextFreeNode-work >0; work++){
			nextFreeNode = new_divide_node_Near(work, nextFreeNode, SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{
		TStruct *work = TREE_START;
		do{
			nextFreeNode = new_divide_node_Near(work, nextFreeNode, PTHREAD);
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
				threadnextFreeNode = new_divide_node_Near(work+j,threadnextFreeNode,RECURSIVE);
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

int Find_Near_test(
		particle *point, 
		int Num_neighbor, 
		TStruct *tree, 
		TPtlStruct *ptl, 
		float *maxr, 
		int *nearindx,
		dptype *DIST2,
		float *MASS, 
		float *densph
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
					DEBUGPRINT("xyz= %g %g %g / xyz= %g %g %g with id= %ld\n",
							point->x,point->y,
							point->z,((TPtlStruct *)ptr)->x,((TPtlStruct *)ptr)->y,
							((TPtlStruct *)ptr)->z,
							( (TPtlStruct*)ptr)-ptl);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	*maxr = maxdist;
	for(i=0;i<npneigh;i++) {
		nearindx[i] = neighbor[i].indx;
		DIST2[i] = neighbor[i].dist2;
		MASS[i] = neighbor[i].mass;
		DEBUGPRINT("final detection %d nearindx: %d dist2= %g den= %g\n", i, nearindx[i], DIST2[i], densph[nearindx[i]]);
	}
	return npneigh;
}

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
		dptype constR,
		int *count
        ){
	int icount = 0;
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
					icount ++;
                }
                ptr = (void *) (((TPtlStruct*)ptr)->sibling);
        }
    }
	*count = icount;
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
    long i,j,k;
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
	new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
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
    long i,j,k;
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
    new_Make_Tree_Near(TREE,nnode, ptl,mp,recursiveflag);
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
    new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
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
    long i,j,k;
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
    new_Make_Tree_Near(TREE,nnode, ptl,nstar,recursiveflag);
    DEBUGPRINT("after Star_Make_Tree_near for nstar = %d\n",nstar);
//  all to star
#pragma omp parallel for private(i,j) schedule(guided)
    for(i=0;i<np;i++){
        int res, icount;
        int tmpindx[NUMNEARDEN];
        dptype tmpd2[NUMNEARDEN];
        float tmpmass[NUMNEARDEN];
        float neardist;
#if defined mADV   || defined ADV
		dptype constR = MIN_CONST_R_SMOOTHING;
		densph[i] = getDenConstR(p+i, TREE,ptl, constR, &icount);
		if(icount < NUMNEARDEN){
	        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
       		densph[i] = 0;
	        for(j=0;j<res;j++){ 
				tmpd2[j] = sqrtf(tmpd2[j]); 
				densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
			}
		}
#else
        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
        densph[i] = 0;
        for(j=0;j<res;j++){
            tmpd2[j] = sqrtf(tmpd2[j]);
            densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
        }
#endif
    }
    DEBUGPRINT0("End of getDenConstR()\n");




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
    new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
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

void lagFindStellarCore(
		SimpleBasicParticleType *bp,
		int np,
		int Numnear,
        float *densph,
		Coretype **Core,
		int *NumCore,  // return value for number of core
		int maxnumcore,
		int **Nearindex // index array for neighbor network
		){
	Coretype *core = *Core;
	int numcore;
//    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    float fplmf,ptlmass;
    long i,j,k;
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





	int *nearindex;
	size_t nnode;
	int recursiveflag; 


	if(1){
		float RG = Gaussian_Smoothing_Length;
		int nbuff = NCELLBUFF;
		double cellsize = TSC_CELL_SIZE;
		int nx,ny,nz;
		double xmin,ymin,zmin,xmax,ymax,zmax;
		float *denGrid;
		LinkedListGrid *linkedListGrid;

		xmin = ymin = zmin = 1e20;
		xmax = ymax = zmax = -1e20;
		for(i=0;i<np;i++){
			xmin = MIN(xmin, bp[i].x);
			ymin = MIN(ymin, bp[i].y);
			zmin = MIN(zmin, bp[i].z);
			xmax = MAX(xmax, bp[i].x);
			ymax = MAX(ymax, bp[i].y);
			zmax = MAX(zmax, bp[i].z);
		}
		xmax = xmax + (xmax-xmin)*0.0001; // making a buffer to avoid the numerical exception
		ymax = ymax + (ymax-ymin)*0.0001;
		zmax = zmax + (zmax-zmin)*0.0001;
		nx = (xmax-xmin)/cellsize + nbuff; // buffer for the periodic boundary condition
		ny = (ymax-ymin)/cellsize + nbuff;
		nz = (zmax-zmin)/cellsize + nbuff;
		xmin -= cellsize*nbuff/2.; ymin -= cellsize*nbuff/2.; zmin -= cellsize*nbuff/2.;

		int mx = 2*(nx/2+1);
		long ncells = mx*ny*nz;
		denGrid = (float*)Malloc(sizeof(float)*ncells,PPTR(denGrid));
		{	
			void assign_density_TSC(SimpleBasicParticleType *, int, float *, int, int, int,
				double, double, double, double);
			assign_density_TSC(bp, np, denGrid, nx,ny,nz,xmin,ymin,zmin,cellsize);
		}
		if(0){

			FILE *wp = fopen("denmap.out","w");
			fwrite(&mx, sizeof(int), 1, wp);
			fwrite(&ny, sizeof(int), 1, wp);
			fwrite(&nz, sizeof(int), 1, wp);
			fwrite(denGrid, sizeof(float), mx*ny*nz,wp);
			fclose(wp);
		}
		{
			void gaussian_Smoothing(float *,int ,int ,int , double , float );
			gaussian_Smoothing(denGrid,nx,ny,nz, cellsize, RG);
		}
		if(1){
			FILE *wp = fopen("gS.denmap.out","w");
			fwrite(&mx, sizeof(int), 1, wp);
			fwrite(&ny, sizeof(int), 1, wp);
			fwrite(&nz, sizeof(int), 1, wp);
			fwrite(denGrid, sizeof(float), mx*ny*nz,wp);
			fclose(wp);
		}
		{
			void findDen(SimpleBasicParticleType *, int, float *, double, double, double,
					double, float *, int, int, int);
			findDen(bp,np, densph, xmin,ymin,zmin, cellsize, denGrid, nx,ny,nz);
		}


		linkedListGrid = (LinkedListGrid*)Malloc(sizeof(LinkedListGrid)
				*ncells, PPTR(linkedListGrid));
		for(i=0;i<ncells;i++) {
			linkedListGrid[i].bp = NULL;
			linkedListGrid[i].np = 0;
		}
		int halfnbuff = NCELLBUFF/2;
		for(i=0;i<np;i++){
			if(bp[i].type == TYPE_STAR){
				long ir = rint((bp[i].x-xmin)/cellsize);
				long jr = rint((bp[i].y-ymin)/cellsize);
				long kr = rint((bp[i].z-zmin)/cellsize);
				long ioff = ir+mx*(jr+ny*kr);
				SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
				linkedListGrid[ioff].bp = bp+i;
				linkedListGrid[ioff].np ++;
				bp[i].bp = tmp;
			}
		}
		DEBUGPRINT0("After building LinkedList\n");
		numcore = 0;
		int nthreads;
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
		int num_cores[nthreads];
		for(i=0;i<nthreads;i++) num_cores[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=nbuff/2;k<nz-nbuff/2;k++){
			int thread_id = omp_get_thread_num();
			for(j=nbuff/2;j<ny-nbuff/2;j++){
				for(i=nbuff/2;i<nx-nbuff/2;i++){
					long ioff = i + mx*(long)(j+ny*k);
					if(denGrid[ioff] > PEAKTHRESHOLD){
						int peakflag=1;
						int i1,j1,k1;
						for(k1=-1;k1<2;k1++) for(j1=-1;j1<2;j1++) for(i1=-1;i1<2;i1++){
							long joff = ioff+i1 + mx*(long)(j1 + ny*k1);
							if(joff != ioff && denGrid[ioff]<=denGrid[joff]){
								peakflag = 0;
								break;
							}
						}
						if(peakflag == 1){
							/*
							SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
							float maxden = -1.e20;
							int ibp;
							while(tmp){
								ibp = tmp-bp;
								if(densph[ibp] > maxden){
									maxden = densph[ibp];
									tmp = tmp->bp;
								}
							}
							core[numcore].peak = ibp; 
							core[numcore].cx = bp[ibp].x; 
							core[numcore].cy = bp[ibp].y; 
							core[numcore].cz = bp[ibp].z; 
							core[numcore].density = maxden;
							numcore ++;
							*/
							num_cores[thread_id] ++; 
						}
					}
				}
			}
		}
		numcore = 0;
		for(i=0;i<nthreads;i++) numcore += num_cores[i];
		DEBUGPRINT("The number of cores: %d in initial finding\n", numcore);

		int ioff_cores[nthreads];
		ioff_cores[0] =0;
		for(i=1;i<nthreads;i++) ioff_cores[i] = ioff_cores[i-1] + num_cores[i-1];
		numcore =0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) firstprivate(numcore)
#endif
        for(k=nbuff/2;k<nz-nbuff/2;k++){
            int thread_id = omp_get_thread_num();
			int koff = ioff_cores[thread_id];
            for(j=nbuff/2;j<ny-nbuff/2;j++){
                for(i=nbuff/2;i<nx-nbuff/2;i++){
                    long ioff = i + mx*(long)(j+ny*k);
                    if(denGrid[ioff] > PEAKTHRESHOLD){
                        int peakflag=1;
                        int i1,j1,k1;
                        for(k1=-1;k1<2;k1++) for(j1=-1;j1<2;j1++) for(i1=-1;i1<2;i1++){
                            long joff = ioff+i1 + mx*(long)(j1 + ny*k1);
                            if(joff != ioff && denGrid[ioff]<=denGrid[joff]){
                                peakflag = 0;
                                break;
                            }
                        }
                        if(peakflag == 1){
                            SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
                            float maxden = -1.e20;
                            long ibp;
							long jbp;
                            while(tmp){
                                ibp = tmp-bp;
                                if(densph[ibp] > maxden){
                                    maxden = densph[ibp];
									jbp = ibp;
                                }
                                tmp = tmp->bp;
                            }
                            core[koff+numcore].peak = jbp; 
                            core[koff+numcore].cx = bp[jbp].x; 
                            core[koff+numcore].cy = bp[jbp].y; 
                            core[koff+numcore].cz = bp[jbp].z; 
                            core[koff+numcore].density = maxden;
                            numcore ++;
                        }
                    }
                }
            }
        }
		numcore = 0;
		for(i=0;i<nthreads;i++) numcore += num_cores[i];
		DEBUGPRINT("The number of cores: %d before MergingPeak\n", numcore);

		Free(linkedListGrid);
		Free(denGrid);
	}


	DEBUGPRINT("C70 has peak id= %d  den= %g before merging\n", core[70].peak, core[70].density);

	int MergingPeak(SimpleBasicParticleType *, int, Coretype *, int, int);
	if(numcore >10) numcore = MergingPeak(bp,np, core, numcore,0);
	core = *Core = Realloc(*Core, sizeof(Coretype)*numcore);

	DEBUGPRINT("C70 has peak id= %d  den= %g after merging\n", core[70].peak, core[70].density);



	DEBUGPRINT("The number of cores : %d and after MergingPeak\n", numcore);

	nearindex = *Nearindex = Realloc(*Nearindex, sizeof(int)*(long)np*Numnear);
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
		ptl[i].mass = 1;
    }
    ptl[np-1].sibling = NULL;
	nnode = MAX(np/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));

//  all to all
	if(1){
		if(nnode > 65*10000) recursiveflag = PTHREAD;
		else recursiveflag = RECURSIVE;
		recursiveflag = SERIALIZED;
	    new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
	}
	else { 
		Box box;
		Box findbox(TPtlStruct *, int);
		box = findbox(ptl, np);
		old_Make_Tree_Near(TREE, ptl,np,box);
	}
    DEBUGPRINT("after All_Make_Tree_near for np = %d\n",np);
#pragma omp parallel for private(i,j,k) schedule(guided)
    for(i=0;i<np;i++){
        int res;
        int tmpindx[NUMNEIGHBOR];
        float mass[NUMNEIGHBOR];
        dptype dist2[NUMNEIGHBOR];
        float neardist;
        k = i*Numnear;
		if(i==90999959){
        	res = Find_Near_test(p+i,Numnear,TREE,ptl,&neardist,tmpindx,
				dist2, mass, densph);
		}
		else {
        	res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,
				dist2, mass);
		}
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

