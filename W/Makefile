# Makefile for Vito
CC= icc
#FC= gfortran
FC= mpiifort
LDFLAGS= -lrt -lm -L. -liotk -g
FFLAGS= -Ofast -g
SRC = Sloss.f PointR.f paths.f gjel.f 
OBJECTS = Sloss.o PointR.o paths.o gjel.o 


.PHONY: clean

all: Loss

Loss: $(OBJECTS)
	${FC} -I. $(OBJECTS) ${LDFLAGS} -o Loss

Sloss.o: Sloss.f
	${FC} ${FFLAGS} -c Sloss.f90

PointR.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f90

paths.o: paths.f
	${FC} ${FFLAGS} -c paths.f90

gjel.o: gjel.f
	${FC} ${FFLAGS} -c gjel.f

clean:
	rm  *.o Loss
