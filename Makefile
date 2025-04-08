#--------------------------------------------------
#--- Makefile for MPI
#--------------------------------------------------

#--- Objects to build 
#OBJS		= gfind.o Memory2.o force_spline.mod3.o Treewalk.hfind.mod2.o \
				subhaloden.mod6.o nrutil.o b2l.o spline.mod2.o \
				subnear.o find_near.o  Treewalk.near.o  utils.o mkRtidal.o

OBJS		= gfind.o Memory2.o force_spline.mod3.o ost.o \
				subhaloden.mod6.o nrutil.o b2l.o spline.mod2.o \
				nnost.o  utils.o mkRtidal.o






#				hfind.sub.mod3.o tree_fortran.o fof.mod8.o   indexx_integer.o  \

PROGS =  gfind.exe
#RANLIB = ranlib
###########################################################################
# GNU CC
###########################################################################
#FC = mpif77
#CC = mpicc
#INCLUDES =  -I/u/kjh2000/include
#DFLAGS = -g #--fast-math #-g -q32  -qtune=pwr3 -qarch=pwr3 -bmaxdata:0x70000000 -bmaxstack:0x70000000 
#FFLAGS       = $(DFLAGS) #-qextname  #-WF,-DPMFORCE
#CFLAGS       = $(DFLAGS) #-DNMEG=500  
#LDFLAGS       = $(DFLAGS) #-qextname  
###########################################################################
# IBM SP3
###########################################################################
#FC = mpxlf
#CC = mpcc
#INCLUDES =  -I/u/kjh2000/include
#DFLAGS = -O3  -q32  -qtune=pwr3 -qarch=pwr3 -bmaxdata:0x70000000 -bmaxstack:0x70000000 
#FFLAGS       = $(DFLAGS) -qextname  #-WF,-DPMFORCE
#CFLAGS       = $(DFLAGS) #-DNMEG=500  
#LDFLAGS       = $(DFLAGS) -qextname  


#FFLAGS       = $(DFLAGS) -fast -inline  all -DPMFORCE 
#FFLAGS1       = $(DFLAGS) -fast  -inline all -nofor_main  -DPMFORCE 
#CFLAGS       = $(DFLAGS) -fast -inline  all  -DPMFORCE  -DCHECK_INFORM 

#FFLAGS       = $(DFLAGS) -g  -DPMFORCE
#FFLAGS1       = $(DFLAGS) -g  -nofor_main -DPMFORCE
#CFLAGS       = $(DFLAGS) -g  -DCHECK_INFORM -DPMFORCE

#MPI_LINKS    = 
##       -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm
##LIBS = -lm  -lmpi
#LIBS = -L/u/kjh2000/lib -lm # -lhdf5 -lz 
##--- C Compiler information
##  Leave the rest untouched
#
#
##--- Suffix-based compilation rules
#.SUFFIXES:
#.SUFFIXES: .exe .o .c .f .F
#
##rules to build binary from source
#
#.c.o :
#	$(CC) $(CFLAGS) $(INCLUDES) -c $<
#.f.o :
#	$(FC) $(FFLAGS) $(INCLUDES) -c $<
#.F.o :
#	$(FC) $(FFLAGS) $(INCLUDES) -c $<
#
#
#
##--- Targets
#SUBDIRS = Tree

all: 
	\rm -f $(PROGS)
	$(MAKE) $(PROGS) 
new: clean 
	$(MAKE) all
clean:
	rm -f $(OBJS)
	rm -f $(PROGS)
$(PROGS): $(OBJS) 
include Rules.make
