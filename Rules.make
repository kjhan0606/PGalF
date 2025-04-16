#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#


AR = ar  rcv
RANLIB = ranlib
#OPT = -O3  -qopenmp -DDEBUG -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL #-DDEBUG
#OPT = -g  -qopenmp -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL #-DDEBUG
#OPT = -O3  -qopenmp -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL #-DDEBUG
OPT = -g  -qopenmp -DINDEX -DVarPM   -DXYZDBL -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 -DNDUST=4  -DDEBUG -DADV



#BINDIR = /applic/compilers/intel/11.1/mpi/openmpi/1.2.6/bin/

#############################################
# List of compilation directives
#############################################
# 
# -DNTREEMIN=nnn  - minimum number of particles in tree - default NTREEMIN=4
# -DGROUPNUM=nnn  - tree grouping number - default GROUPNUM=10
# -DEPSSQ=fff     - softening length squared for PM force correction 
#                 - default=4.0e-8
# -DNMEG=nnn      - number of megabytes of storage per processor - default=40
# -DINCLUDE_TREE_FORCE - define if you want tree code corrections
# -DTREE               - define if you want tree code corrections
# -DGIGA               - run in single mass mode for maximum memory usage
# -DTREEFIX            - correct treebuilding code for closely spaced particles
# -DMOVIE              - include image rendering code
# -DACCOUT             - test mode - read in file 'rvtmp' and dump accelerations
# -DDEBUG              - output debugging information
# -DVERBOSE            - output even more debugging information
# -DBIT64              - use on machines using 64 bit addressing 
#                       - i.e. Compaq alpha
# -DAIX                - use on IBM SP3 machines
# -DDYNAMIC_LOAD_BALANCE  - use on systems with heterogeneous processors
# -DOBSERVATION_MODE   - save halo data at every ObservationStride'th step
# -DSAVESLICE          - save slice density file at every time step
######################################################
# CURRENTLY NOT CHECKED - April 23, 2001 - jjd
#############################################################
# -DENOUGH_SPACE   : if there is enough space for memory
#    -DFDA_FORCE_ARRAY_SWP  : force array + swap space
#    -DTSC_SWAP : tsc swap
#############################################################
#
#####################################
# System specific definitions
##################################
# GNU CC
##################################
#FC = mpif77
#CC = mpicc
#INCLUDES = -I/packages/include
#FDFLAGS =
#DFLAGS = -DNMEG=1350 -D_FILE_OFFSET_BITS=64  -D_LARGEFILE_SOURCE
#CFLAGS = -O3 -ffast-math $(DFLAGS)
#LDFLAGS = -O3 -ffast-math
#LIBS = -L/packages/lib -L/usr/local/lib\
#	-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm 
##################################
# PG compiler
##################################
FC = mpiifort
CC = mpiicc

#DFLAGS = -DNMEG=3300L -D_LARGE_FILES -DENDIANCHANGE  -DCHECK_INFORM
#DFLAGS = -DNMEG=7000L -D_LARGE_FILES -DENDIANCHANGE   #-DCHECK_INFORM
DFLAGS = -DNMEG=400000L -D_LARGE_FILES #-DCHECK_INFORM
#
#FFLAGS = $(FDFLAGS) -O3 -qextname -q32 -qtune=pwr3 -qarch=pwr3 -bmaxdata:2000000000 -bmaxstack:2000000000
#CFLAGS = -O3 $(DFLAGS)  -q32 -qtune=pwr3 -qarch=pwr3 -bmaxdata:2000000000 -bmaxstack:2000000000
#LDFLAGS = -O3 -qextname -q32 -qtune=pwr3 -qarch=pwr3 -bmaxdata:2000000000 -bmaxstack:2000000000
FFLAGS = $(FDFLAGS) $(OPT)  -mcmodel=medium
CFLAGS = $(DFLAGS)  $(OPT)  #-DDEBUG
LDFLAGS = $(OPT) -mcmodel=medium


LIBS = -lm
##################################
# IBM SP3
##################################
#FC = mpxlf_r
#CC = mpcc_r
#FDFLAGS = -WF,-DINCLUDE_TREE_FORCE -WF,-DAIX -WF,-DCHECK_INFORM
#DFLAGS = -DNMEG=5900 -D_LARGE_FILES  -DCHECK_INFORM
#
#FFLAGS = $(FDFLAGS) -O3 -qextname -q64 -qtune=auto -qarch=auto 
#CFLAGS = -O3 $(DFLAGS)  -q64 -qtune=auto -qarch=auto 
#LDFLAGS = -O3 -qextname -q64 -qtune=auto -qarch=auto

#LIBS = -L/u/kjh2000/lib -L./PREWORK -L./WORKNEW -lmpi -lm
#################################
# Compaq Alpha
#################################
#FC = f77
#CC = cc 
#INCLUDES = -I/home/kjhan/dolphin/fftw/include
#FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DTREE
#DFLAGS = -DBIT64 -DNMEG=256 -DINCLUDE_TREE_FORCE -DTREE -DGIGA -DTREEFIX
#FFLAGS = $(FDFLAGS)  -fast -nofor_main
#CFLAGS = $(DFLAGS)  -fast
#LDFLAGS =  -fast -nofor_main
#LIBS = -L/home/kjhan/dolphin/fftw/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm -lz -lmpi -lm
##################################

#--- C Compiler information
#  Leave the rest untouched


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
	$(CC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.f.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.for.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.F.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

#--- Targets
