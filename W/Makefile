# Makefile for Vito
CC= icc
#FC= gfortran
FC= ifort
LDFLAGS= -lm -L. -liotk -g -mkl
FFLAGS= -Ofast -g -traceback -fpe0 
SRC = Sloss.f PointR.f paths.f gjel.f 
OBJECTS = Sloss.o PointR.o paths.o gjel.o 


.PHONY: clean

all: Loss

Loss: $(OBJECTS)
	${FC} -I. $(OBJECTS) ${LDFLAGS} -o Loss

Sloss.o: Sloss.f
	${FC} ${FFLAGS} -c Sloss.f90

PointR.o: PointR.f
	${FC} ${FFLAGS} -c PointR_mkl.f90 -o PointR.o

paths.o: paths.f
	${FC} ${FFLAGS} -c paths.f90

gjel.o: gjel.f
	${FC} ${FFLAGS} -c gjel.f

clean:
	rm  *.o Loss
