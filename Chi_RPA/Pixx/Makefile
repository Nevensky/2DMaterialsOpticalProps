# Makefile for Vito
CC= gcc
#FC= gfortran
#FC=mpif90
FC=ifort
# LDFLAGS= -lm -L. -liotk -g -traceback #-lrt 
LDFLAGS = -L. -liotk -mkl 
LDFLAGS_PARA = -qopenmp -L. -liotk -mkl
DEBUG = -g -traceback -fpe0 -fp-model precise -fp-stack-check -heap-arrays # -debug all -debug-parameters all -DD -C -warn all -check all
FFLAGS = -Ofast ${DEBUG}
SRC = PointR.f Optabs.f paths.f gjel.f 
OBJECTS = PointR.o Optabs.o paths.o gjel.o 
OBJECTS_LEGACY = PointR_legacy.o Optabs_legacy.o paths.o gjel.o 
OBJECTS_PARA = PointR.o Optabs_para.o paths.o gjel.o 
OBJECTS_QE6 = PointR.o Optabs_para_qe6.o paths.o gjel.o 


.PHONY: clean

all: Loss

qe6: Loss_qe6

serial: Loss_serial

legacy: Loss_legacy

Loss: $(OBJECTS)
	${FC} $(OBJECTS_PARA) ${LDFLAGS_PARA} -o Loss

Loss_serial: $(OBJECTS)
	${FC} $(OBJECTS) ${LDFLAGS} -o Loss

Loss_qe6: $(OBJECTS)
	${FC} $(OBJECTS_QE6) ${LDFLAGS_PARA} -o Loss

Loss_legacy: $(OBJECTS)
	${FC} $(OBJECTS_LEGACY) ${LDFLAGS} -o Loss

Optabs.o: Optabs.f
	${FC} ${FFLAGS} -c  Optabs.f90

Optabs_para.o: Optabs.f
	${FC} ${FFLAGS} ${PARA} -c  Optabs_para.f90 -o Optabs_para.o

Optabs_para_qe6.o: Optabs.f
	${FC} ${FFLAGS} ${PARA} -c  Optabs_qe6.f90	-o Optabs_para_qe6.f90

Optabs_legacy.o: Optabs.f
	${FC} ${FFLAGS} -c  Optabs.f -o Optabs_legacy.o

PointR.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f90 -o PointR.o
	${FC} ${FFLAGS} -c PointR.f90

PointR_legacy.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f -o PointR_legacy.o

paths.o: paths.f
	${FC} ${FFLAGS} -c paths.f90

gjel.o: gjel.f
	${FC} ${FFLAGS} -c gjel.f

clean:
	rm  *.o Loss fbz_check.dat modpointr.mod
