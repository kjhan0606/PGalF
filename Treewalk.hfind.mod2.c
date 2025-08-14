#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "ramses.h"
#include "tree.h"

/*
#include "fof.h"
*/

#include "mathsub.h"
#include "force_spline.h"
extern float diff[NSPLINE][3],slope[NSPLINE][3],ran2nran;
/*
float diff[NSPLINE][3],slope[NSPLINE][3],ran2nran;
*/
enum boolean fof_open(particle p,FoFTStruct *tree, float fof_link){
    float tmpx,tmpy,tmpz;
    float dist2,dist,r,diffr;
    float ratio, link02;
    tmpx = p.x-tree->monox;
    tmpy = p.y-tree->monoy;
    tmpz = p.z-tree->monoz;
    r = tree->dist;
    dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
    dist = sqrt(dist2);
    diffr = dist - r;
	link02 = 0.5*(p.link02 + tree->maxlink02);
	float minlink02 =0.5*(p.link02 + tree->minlink02);
	if(minlink02 > (dist+r)) return ENCLOSE;
	else if(diffr <= link02) return YES;
    else return NO;
}

enum boolean treeopen(particle *p,TStruct *tree, float theta2){
	dptype tmpx,tmpy,tmpz; 
	dptype dist2; 
	dptype ratio; 
	tmpx = p->x-tree->monox;
	tmpy = p->y-tree->monoy;
	tmpz = p->z-tree->monoz;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	ratio = tree->dist2/theta2; 
	if(dist2 < ratio) return YES;
	else return NO;
}
#define CELLFORCE(p,ptr,force,ran2nran) \
{\
	int ntmp;   \
	dptype tmpx,tmpy,tmpz,dist2,dist;   \
	dptype xx,yy,zz,xy,xz,yz,tmpxx; \
	dptype qxx1,qxx2,qxx3,qxy,tmpdist;   \
	dptype fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   \
	TStruct *pointer;   \
	pointer = (TStruct *) ptr;   \
	tmpx = pointer->monox - p->x;   \
	tmpy = pointer->monoy - p->y; \
	tmpz = pointer->monoz - p->z;   \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	xx = pointer->quad[0];   \
	yy = pointer->quad[1];   \
	zz = pointer->quad[2];   \
	xy = pointer->quad[3];   \
	xz = pointer->quad[4];   \
	yz = pointer->quad[5];   \
	tmpxx = xx + yy + zz;   \
	qxx1 = xx*tmpx+xy*tmpy+xz*tmpz;   \
	qxx2 = xy*tmpx+yy*tmpy+yz*tmpz;   \
	qxx3 = xz*tmpx+yz*tmpy+zz*tmpz;   \
	qxy = xx*tmpx*tmpx+yy*tmpy*tmpy+zz*tmpz*tmpz+   \
		2.*(xy*tmpx*tmpy+xz*tmpx*tmpz+yz*tmpy*tmpz); \
	dist = sqrt(dist2);\
	ntmp = dist/ran2nran;   \
	tmpdist = (dist-ntmp*ran2nran);   \
	fplmf1 = diff[ntmp][0] + slope[ntmp][0]*tmpdist;   \
	fplmf2 = diff[ntmp][1] + slope[ntmp][1]*tmpdist;   \
	fplmf3 = diff[ntmp][2] + slope[ntmp][2]*tmpdist;   \
	tmptmp = pointer->mass * fplmf1 + tmpxx*fplmf2 +   \
		qxy*fplmf3;   \
	force->x += tmpx*tmptmp+2.*qxx1*fplmf2; \
	force->y += tmpy*tmptmp+2.*qxx2*fplmf2; \
	force->z += tmpz*tmptmp+2.*qxx3*fplmf2; \
}

#define PARTICLEFORCE(p,ptr,force,ran2nran) \
{ \
	int ntmp;  \
	dptype tmpx,tmpy,tmpz,dist2,dist;  \
	dptype fplmf;  \
	dptype ptlmass; \
	TPtlStruct *pointer; \
	pointer = (TPtlStruct *) ptr; \
	ptlmass = pointer->mass; \
	tmpx = pointer->x - p->x;  \
	tmpy = pointer->y - p->y;  \
	tmpz = pointer->z - p->z;  \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	dist = sqrt(dist2);\
	ntmp = dist/ran2nran;  \
	fplmf = diff[ntmp][0] + slope[ntmp][0] * (dist-ntmp*ran2nran);  \
	fplmf *= ptlmass; \
	force->x +=  tmpx * fplmf;  \
	force->y +=  tmpy * fplmf;  \
	force->z +=  tmpz * fplmf;  \
}

#define CELLPLUMPOTENT(p,ptr,potent,epsilon2) \
{\
	int ntmp;   \
	dptype tmpx,tmpy,tmpz,dist2,dist;   \
	dptype xx,yy,zz,xy,xz,yz,tmpxx; \
	dptype qxx1,qxx2,qxx3,qxy,tmpdist;   \
	dptype fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   \
	dptype idist2,sqrtdist2,isqrtdist2; \
	TStruct *pointer;   \
	pointer = (TStruct *) ptr;   \
	tmpx = p->x - pointer->monox;   \
	tmpy = p->y - pointer->monoy;   \
	tmpz = p->z - pointer->monoz;   \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	xx = pointer->quad[0];   \
	yy = pointer->quad[1];   \
	zz = pointer->quad[2];   \
	xy = pointer->quad[3];   \
	xz = pointer->quad[4];   \
	yz = pointer->quad[5];   \
	tmpxx = xx + yy + zz;   \
	qxy = xx*tmpx*tmpx+yy*tmpy*tmpy+zz*tmpz*tmpz+   \
		2.*(xy*tmpx*tmpy+xz*tmpx*tmpz+yz*tmpy*tmpz); \
	dist = sqrt(dist2);\
	ntmp = dist/ran2nran;   \
	ntmp = MIN(ntmp,NSPLINE-1);\
	tmpdist = (dist-ntmp*ran2nran);   \
	fplmf1 = diff[ntmp][0] + slope[ntmp][0]*tmpdist;   \
	fplmf2 = diff[ntmp][1] + slope[ntmp][1]*tmpdist;   \
	fplmf3 = diff[ntmp][2] + slope[ntmp][2]*tmpdist;   \
	potent += pointer->mass * fplmf1 + tmpxx*fplmf2 + \
		qxy*fplmf3; \
}

#define PARTICLEPLUMPOTENT(p,ptr,potent,epsilon2) \
{\
	float tmpx,tmpy,tmpz;   \
	TPtlStruct *pointer; \
	float fplmf,ptlmass,dist2,dist; \
	int ntmp; \
	pointer = (TPtlStruct *) ptr; \
	ptlmass = pointer->mass; \
	tmpx = p->x -  pointer->x; \
	tmpy = p->y -  pointer->y; \
	tmpz = p->z -  pointer->z; \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	if(dist2 >0) fplmf = -ptlmass/sqrt(dist2+epsilon2);\
	else fplmf = 0;\
	potent += fplmf; \
}


/* 
 * *p is the position at which you want to calculate force using tree
 * theta2 is the opening angle for tree walk
 * *tree is the tree structure.
 */

void treeforce(particle *p,float theta2,TStruct *tree, TPtlStruct *ptl,
		pforce *force){
	void *ptr,*optr;
	ptr = (void *)tree;
	force->x = force->y = force->z = 0.;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(treeopen(p,ptr,theta2)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLFORCE(p,ptr,force,ran2nran);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEFORCE(p,ptr,force,ran2nran);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return;
}
float treeplumpotential(particle *p,float theta2,TStruct *tree,
		TPtlStruct *ptl){
	void *ptr;
	float potent=0.;
	float epsilon2;
	epsilon2 = EPSILON * EPSILON;
	ptr = (void *) tree;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(treeopen(p,ptr,theta2)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLPLUMPOTENT(p,ptr,potent,epsilon2);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEPLUMPOTENT(p,ptr,potent,epsilon2);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return potent;
}

int ScoopUpParticles(void *optr, particle *linked, int oncount){
	int ncount=oncount;
	void *ptr = optr;
	void *terminal = ptr->sibling;
	while(ptr != terminal){
		switch(((TYPE*)ptr)->type){
			case TYPE_TREE:
				ptr = ((FoFTStruct*)ptr)->daughter;
				break;
			default : 
				FoFTPtlStruct *tmp = (FoFTPtlStruct*)ptr;
				if(((FoFTPtlStruct*)ptr)->included == NO){
					linked[ncount].x = tmp->x;
					linked[ncount].y = tmp->y;
					linked[ncount].z = tmp->z;
					linked[ncount].link02 = tmp->link02; 
					tmp->included = YES; 
					ncount ++;
				}
				ptr = (void*)(tmp->sibling);
		}
	}
	((FoFTStruct *)ptr)->sibling = ((FoFTStruct *)ptr)->daughter;
	return ncount;
}
int ompFoFGroup(particle *p,POSTYPE fof_link,FoFTStruct *tree,
        FoFTPtlStruct *ptl,int nptl, particle *linked){
    int ncount, now;
    void *ptr,*optr,*nptr;
    POSTYPE tmpx,tmpy,tmpz;
    POSTYPE fof_link2,dist2;
    particle point;
    ncount = now = 0;
	point = *p;
    do {
        ptr = (void*) tree;
        while(ptr != NULL){
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
                    switch(fof_open(point,ptr,fof_link)){
                        case NO:
                            ptr = (void *)(((FoFTStruct*)ptr)->sibling);
                            break;
                        default :
                            ptr = (void *)(((FoFTStruct*)ptr)->daughter);
                    }
                    break;
                default :
                    if(((FoFTPtlStruct*)ptr)->included == NO) {
                        tmpx = point.x - ((FoFTPtlStruct*)ptr)->x;
                        tmpy = point.y - ((FoFTPtlStruct*)ptr)->y;
                        tmpz = point.z - ((FoFTPtlStruct*)ptr)->z;
                        dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						dist2 = sqrt(dist2);
                        if(dist2 <= 0.5*(point.link02+((FoFTPtlStruct*)ptr)->link02))
						{
                            linked[ncount].x = ((FoFTPtlStruct*)ptr)->x;
                            linked[ncount].y = ((FoFTPtlStruct*)ptr)->y;
                            linked[ncount].z = ((FoFTPtlStruct*)ptr)->z;
                            linked[ncount].link02 = ((FoFTPtlStruct*)ptr)->link02;
#pragma omp atomic write
                            ((FoFTPtlStruct*)ptr)->included = YES;
                            ncount ++;
                        }
                    }
                    ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
            }
        }
        point = linked[now];
        now ++;
    } while( now <= ncount);
    return ncount;
}

typedef struct pthreadStruct{
	POSTYPE fof_link;
	FoFTStruct *tree;
	FoFTPlStruct *ptl;
	int nptl, id, nid, haloid;
} pthreadStruct;

int pthreadFoFGroup(void *args){
	pthreadStruct *p = (pthreadStruct *)args;
	POSTYPE fof_link = p->fof_link;
	FoFTSTruct *tree = p->tree;
	FoFTPtlStruct *ptl = p->ptl;
	int nptl = p->nptl;
	int pid = p->id;
	int nid = p->nid;
	int haloid = p->haloid;
	void *ptr;
	{
		ptr = (void *) tree;
		while(ptr){
			switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
                    switch(fof_open(point,ptr,fof_link)){
						case NO:
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
							break;
						default:
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
					}
					break;
				default:
					FoFTPtlStruct *tmp = (FoFTPtlStruct*)ptr;
					int ioff = tmp - ptl;
					if(ioff%nid == pid && tmp->included == NO) {
						POSTYPE tmpx = point.x - tmp->x;
                        POSTYPE tmpy = point.y - tmp->y;
                        POSTYPE tmpz = point.z - tmp->z;
                        POSTYPE dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
                        dist2 = sqrt(dist2);
                        if(dist2 <= 0.5*(point.link02+tmp->link02)){
							tmp->included =NEW;
							tmp->haloindx =haloid;
						}
					}
					ptr = (void*)(tmp->sibling);
			}

		}
	}
}

int new_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree,
        FoFTPtlStruct *ptl,particle *linked){
    int ncount, now;
    void *ptr,*optr,*nptr;
    POSTYPE tmpx,tmpy,tmpz;
    POSTYPE fof_link2,dist2;
    particle point;
    ncount = now = 0;
	point = *p;
    do {
        ptr = (void*) tree;
        while(ptr != NULL){
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
                    switch(fof_open(point,ptr,fof_link)){
						case ENCLOSE:
							ncount = ScoopUpParticles(ptr,linked,ncount);
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling); 
                        case YES:
                            ptr = (void *)(((FoFTStruct*)ptr)->daughter);
                            break;
                        default :
                            ptr = (void *)(((FoFTStruct*)ptr)->sibling);
                    }
                    break;
                default :
                    if(((FoFTPtlStruct*)ptr)->included == NO) {
                        tmpx = point.x - ((FoFTPtlStruct*)ptr)->x;
                        tmpy = point.y - ((FoFTPtlStruct*)ptr)->y;
                        tmpz = point.z - ((FoFTPtlStruct*)ptr)->z;
                        dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						dist2 = sqrt(dist2);
                        if(dist2 <= 0.5*(point.link02+((FoFTPtlStruct*)ptr)->link02))
						{
                            linked[ncount].x = ((FoFTPtlStruct*)ptr)->x;
                            linked[ncount].y = ((FoFTPtlStruct*)ptr)->y;
                            linked[ncount].z = ((FoFTPtlStruct*)ptr)->z;
                            linked[ncount].link02 = ((FoFTPtlStruct*)ptr)->link02;
                            ((FoFTPtlStruct*)ptr)->included = YES;
                            ncount ++;
                        }
                    }
                    ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
            }
        }
        point = linked[now];
        now ++;
    } while( now <= ncount);
    return ncount;
}

int destroy_new_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree,
        FoFTPtlStruct *ptl,particle *linked){
    int ncount, now;
    void *ptr,*optr,*nptr;
    POSTYPE tmpx,tmpy,tmpz;
    POSTYPE fof_link2,dist2;
    particle point;
    ncount = now = 0;
    point = *p;
    do {
        ptr = (void*) tree;
        optr = (void*) tree;
        while(ptr != NULL){
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr =  ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(fof_open(point,ptr,fof_link)){ 
						case ENCLOSE:
							ncount = ScoopUpParticles(ptr,linked,ncount);
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling); 
							break;
						case YES: 
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter); 
							break; 
						default : 
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling); 
					}
                    break;
                default :
                    if(((FoFTPtlStruct*)ptr)->included == YES) {
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else{
                        tmpx = point.x - ((FoFTPtlStruct*)ptr)->x;
                        tmpy = point.y - ((FoFTPtlStruct*)ptr)->y;
                        tmpz = point.z - ((FoFTPtlStruct*)ptr)->z;
                        dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
                        dist2 = sqrt(dist2);
                        if(dist2 <= 0.5*(point.link02+((FoFTPtlStruct*)ptr)->link02))
                        {
                            linked[ncount].x = ((FoFTPtlStruct*)ptr)->x;
                            linked[ncount].y = ((FoFTPtlStruct*)ptr)->y;
                            linked[ncount].z = ((FoFTPtlStruct*)ptr)->z;
                            linked[ncount].link02 = ((FoFTPtlStruct*)ptr)->link02;
                            ((FoFTPtlStruct*)ptr)->included = YES;

                            ncount ++;
							nptr = ((FoFTPtlStruct *)ptr)->sibling;
							EraseFromTree(optr,ptr,nptr);
                        }
						else optr = ptr;
                    	ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
                    }
            }
        }
        point = linked[now];
        now ++;
    } while( now <= ncount);
    return ncount;
}


void FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,size_t np,Box box){
	FoFBeginEndTree beginend;
	FoFTStruct *NewTree;
	/*
	TPtlStruct *ptr;
	ptr = ptl;
	while(ptr!= NULL){
		ptr=ptr->sibling;
	}
	*/
	TREE_START->sibling = NULL;
	beginend = FoF_divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
}

void Make_Tree(TStruct *TREE_START,TPtlStruct *ptl,int np,Box box){
	BeginEndTree beginend;
	TStruct *NewTree;
	/*
	TPtlStruct *ptr;
	ptr = ptl;
	while(ptr!= NULL){
		ptr=ptr->sibling;
	}
	*/
	TREE_START->sibling = NULL;
	beginend = divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
}

FoFBeginEndTree FoF_divide_node(FoFTStruct *TREE_START,FoFTStruct *NewTree, 
		FoFTPtlStruct *ptl, Box box,FoFTStruct *ThisTree){ 
	FoFBeginEndTree beginend;
	FoFTStruct *p2tree,tmpnode[8];
	FoFTStruct *NowCal;
	void *from_sibling;
	FoFTPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	float x0,y0,z0,inv_halfw,halfw;
	float tmpx,tmpy,tmpz,tmpdist2,distmax;
	float ptlmass;
	int count;

	/*   */
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->x0 = box.x;
	ThisTree->y0 = box.y;
	ThisTree->z0 = box.z;
	ThisTree->maxlink02 = ThisTree->minlink02 = 0;
	ThisTree->Nparticle = 0;
	ThisTree->monox= ThisTree->monoy= ThisTree->monoz= 0.;

	p2ptl = ptl;
	while(p2ptl != NULL){
		ThisTree->Nparticle ++;
		ThisTree->monox += p2ptl->x;
		ThisTree->monoy += p2ptl->y;
		ThisTree->monoz += p2ptl->z;
		p2ptl = p2ptl->sibling;
	}
	ThisTree->monox /= ThisTree->Nparticle;
	ThisTree->monoy /= ThisTree->Nparticle;
	ThisTree->monoz /= ThisTree->Nparticle;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		tmpx = p2ptl->x - ThisTree->monox;
		tmpy = p2ptl->y - ThisTree->monoy;
		tmpz = p2ptl->z - ThisTree->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		ThisTree->maxlink02 = MAX(ThisTree->maxlink02,p2ptl->link02);
		ThisTree->minlink02 = MIN(ThisTree->minlink02,p2ptl->link02);
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	ThisTree->dist = sqrt(distmax);
	p2ptl = ptl;
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	/* initialize temporary tree array */
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		tmpbox[i].x = x0+IMOD(i,2)*halfw;
		tmpbox[i].y = y0+(IMOD(i,4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}
	while(p2ptl != NULL){
		mx = ( p2ptl->x >= ThisTree->monox ? 1:0);
		my = ( p2ptl->y >= ThisTree->monoy ? 1:0);
		mz = ( p2ptl->z >= ThisTree->monoz ? 1:0);
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	/* Making link from Mother Node */
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle > 0) break;
	}
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			NewTree->daughter = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
			from_sibling = NewTree;
			NewTree++;
			count ++;
		}
		else if(tmpnode[i].Nparticle > 0 ){
			tmpptr = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = tmpptr;
			while(tmpptr != NULL){
				from_sibling = tmpptr;
				tmpptr = tmpptr->sibling;
			}
		}
	}
	((GENERAL_TPtl_POINTER*)from_sibling)->sibling = ThisTree->sibling;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			beginend = FoF_divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
BeginEndTree divide_node(TStruct *TREE_START,TStruct *NewTree, 
		TPtlStruct *ptl, Box box,TStruct *ThisTree){ 
	BeginEndTree beginend;
	TStruct *p2tree,tmpnode[8];
	TStruct *NowCal;
	void *from_sibling;
	TPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	float x0,y0,z0,inv_halfw,halfw;
	float tmpx,tmpy,tmpz,tmpdist2,distmax;
	float ptlmass;
	int count;

	/*   */
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->x0 = box.x;
	ThisTree->y0 = box.y;
	ThisTree->z0 = box.z;
	ThisTree->mass = 0;
	ThisTree->Nparticle = 0;
	ThisTree->monox= ThisTree->monoy= ThisTree->monoz= 0.;
	ThisTree->quad[0] = ThisTree->quad[1] = ThisTree->quad[2] = 0.;
	ThisTree->quad[3] = ThisTree->quad[4] = ThisTree->quad[5] = 0.;
	ThisTree->mrr = 0.;

	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		ThisTree->mass += ptlmass;
		ThisTree->Nparticle ++;
		ThisTree->monox += ptlmass * p2ptl->x;
		ThisTree->monoy += ptlmass * p2ptl->y;
		ThisTree->monoz += ptlmass * p2ptl->z;
		p2ptl = p2ptl->sibling;
	}
	ThisTree->monox /= ThisTree->mass;
	ThisTree->monoy /= ThisTree->mass;
	ThisTree->monoz /= ThisTree->mass;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		tmpx = p2ptl->x - ThisTree->monox;
		tmpy = p2ptl->y - ThisTree->monoy;
		tmpz = p2ptl->z - ThisTree->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		ThisTree->mrr += ptlmass * tmpdist2;
		ThisTree->quad[0] += ptlmass * tmpx * tmpx;
		ThisTree->quad[1] += ptlmass * tmpy * tmpy;
		ThisTree->quad[2] += ptlmass * tmpz * tmpz;
		ThisTree->quad[3] += ptlmass * tmpx * tmpy;
		ThisTree->quad[4] += ptlmass * tmpx * tmpz;
		ThisTree->quad[5] += ptlmass * tmpy * tmpz;
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	p2ptl = ptl;
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	/* initialize temporary tree array */
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		tmpbox[i].x = x0+IMOD(i,2)*halfw;
		tmpbox[i].y = y0+(IMOD(i,4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}
	while(p2ptl != NULL){
		mx = (int)((p2ptl->x - x0)*inv_halfw);
		my = (int)((p2ptl->y - y0)*inv_halfw);
		mz = (int)((p2ptl->z - z0)*inv_halfw);
		mx = MIN(1,MAX(0,mx));
		my = MIN(1,MAX(0,my));
		mz = MIN(1,MAX(0,mz));
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	/* Making link from Mother Node */
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle > 0) break;
	}
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			NewTree->daughter = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
			from_sibling = NewTree;
			NewTree++;
			count ++;
		}
		else if(tmpnode[i].Nparticle > 0 ){
			tmpptr = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = tmpptr;
			while(tmpptr != NULL){
				from_sibling = tmpptr;
				tmpptr = tmpptr->sibling;
			}
		}
	}
	((GENERAL_TPtl_POINTER*)from_sibling)->sibling = ThisTree->sibling;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			beginend = divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
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
Box foffindbox(FoFTPtlStruct *ptl, int nend){
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
#undef OFFSET

