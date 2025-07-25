#--------------------------------------------------
#--- Makefile for MPI
#--------------------------------------------------

#--- Objects to build 
OBJS		= gfind.o Memory2.o force_spline.mod3.o ost.o \
				subhaloden.mod6.o nrutil.o b2l.o spline.mod2.o \
				nnost.o  utils.o mkRtidal.o tsc_omp2.o gsmooth.o

PROGS =  gfind.exe
all: 
	\rm -f $(PROGS)
	$(MAKE) $(PROGS) 
new: clean 
	$(MAKE) all
clean:
	rm -f $(OBJS)
	rm -f $(PROGS)
$(PROGS): $(OBJS) 

AR = ar  rcv
RANLIB = ranlib
FFTW = /home/kjhan/local/
OPT = -g  -qopenmp -DINDEX -DVarPM   -DXYZDBL -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 \
		-DNDUST=4  -DDEBUG -DADV 

INCLUDES = -I$(FFTW)/include
FC = mpifort
CC = mpicc

DFLAGS = -DNMEG=400000L #-D_LARGE_FILES #-DCHECK_INFORM
FFLAGS = $(FDFLAGS) $(OPT)  -mcmodel=medium
CFLAGS = $(DFLAGS)  $(OPT) $(INCLUDES)  #-DDEBUG


LIBS = -L$(FFTW)/lib -lsrfftw -lsfftw -lm

#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c .f .F

#rules to build binary from source
.c.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.f.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.for.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.c.exe :
	$(CC) $(OPT) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.f.exe :
	$(FC) $(OPT) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.for.exe :
	$(FC) $(OPT) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.F.exe :
	$(FC) $(OPT) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

#--- Targets
