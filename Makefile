#
# Unix/Linux GCC Makefile for pose estimation
#

CC=gcc
LEVMAR_PATH=/home/lourakis/levmar-src/levmar-2.6/ # CHANGE THIS TO POINT TO YOUR COMPILED COPY OF LEVMAR
INCLUDES=-I$(LEVMAR_PATH) -Imlsl
CFLAGS=$(INCLUDES) -O3 -funroll-loops -Wall #-g #-pg
LAPACKLIBS_PATH=/usr/local/lib # WHEN USING LAPACK, CHANGE THIS TO WHERE YOUR COMPILED LIBS ARE!
LDFLAGS=-L. -L$(LEVMAR_PATH) -L$(LAPACKLIBS_PATH)

LIBOBJS=buckets.o lqs.o ransac.o rngs.o prosac.o deal.o p3p.o p4pf.o planep4p.o svd3.o polysolve.o\
poseproj.o posest.o align.o lhm.o sam.o mlsl/mlsl.o mlsl/mt19937ar.o mlsl/redblack.o mlsl/sobolseq.o

LIBSRCS=buckets.c lqs.c ransac.c rngs.c prosac.c deal.c p3p.c p4pf.c planep4p.c svd3.c polysolve.c\
poseproj.c posest.c align.c lhm.c sam.c mlsl/mlsl.c mlsl/mt19937ar.c mlsl/redblack.c mlsl/sobolseq.c

DEMOOBJS=posest_demo.o
DEMOSRCS=posest_demo.c
AR=ar
RANLIB=ranlib
RM=rm -f
MAKE=make
MAPLE=maple

# On systems with a FORTRAN (not f2c'ed) version of LAPACK, -lf2c is
# not necessary; on others, -lf2c is equivalent to -lF77 -lI77
LAPACKLIBS=-llapack -lblas -lf2c

all: libposest.a posest_demo

libposest.a: $(LIBOBJS)
	$(AR) crv libposest.a $(LIBOBJS)
	$(RANLIB) libposest.a

poseproj.o: poseproj.h
posest.o: util.h posest.h poseproj.h lqs.h
lqs.o: lqs.h compiler.h
p3p.o: compiler.h p3p.h polysolve.h
p4pf.o: compiler.h p4pf.h
polysolver.o: polysolve.h
    
posest_demo: $(DEMOOBJS) libposest.a
	$(CC) $(LDFLAGS) $(DEMOOBJS) -o posest_demo -lposest -llevmar $(LAPACKLIBS) -lm 

posest_demo.o: posest.h

#poseproj.c: clc-poseproj.mpl
#	$(MAPLE) <  $<

clean:
	-$(RM) $(LIBOBJS) $(DEMOOBJS) gmon.out
	-cd matlab; make clean

cleanall: clean
	-$(RM) libposest.a

depend:
	makedepend $(INCLUDES) -f Makefile $(LIBSRCS) $(DEMOSRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

