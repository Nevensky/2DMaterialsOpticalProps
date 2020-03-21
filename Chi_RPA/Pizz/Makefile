# Makefile for Vito
CC= gcc
#FC= gfortran
FC=mpif90
LDFLAGS= -lrt -lm -L. -liotk -g
FFLAGS=-O3 -g
SRC = Optabs.f PointR.f paths.f gjel.f 
OBJECTS = Optabs.o  PointR.o paths.o gjel.o 


.PHONY: clean

all: Loss

Loss: $(OBJECTS)
	${FC} $(OBJECTS) ${LDFLAGS} -o Loss

Optabs.o: Optabs.f
	${FC} ${FFFLAGS} -c  Optabs.f

PointR.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f

paths.o: paths.f
	${FC} ${FFLAGS} -c paths.f

gjel.o: gjel.f
	${FC} ${FFLAGS} -c gjel.f

clean:
	rm  *.o Loss
