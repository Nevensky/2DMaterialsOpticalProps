# Makefile for Vito
CC= gcc
#FC= gfortran
#FC=mpif90
FC=ifort
# LDFLAGS= -lm -L. -liotk -g -traceback #-lrt 
LIBIOTK=-L. -liotk
LDFLAGS = -qmkl 
LDFLAGS_PARA = -qopenmp -qmkl
OPTREPORT=-qopt-report-phase=vec
DEBUG = -g -traceback # -warn all -check all -debug all#-fpe0 -fp-model precise -fp-stack-check -heap-arrays # -debug all -debug-parameters all -DD -C -warn all -check all
FFLAGS = -Ofast -fpp  ${DEBUG} # ${OPTREPORT}
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

Loss: $(OBJECTS_PARA)
	${FC} $(OBJECTS_PARA) ${LDFLAGS_PARA} ${LIBIOTK} -o Loss

Loss_serial: $(OBJECTS)
	${FC} $(OBJECTS) ${LDFLAGS} ${LIBIOTK} -o Loss

Loss_qe6: $(OBJECTS_QE6)
	${FC} $(OBJECTS_QE6) ${LDFLAGS_PARA} -o Loss

Loss_legacy: $(OBJECTS_LEGACY)
	${FC} $(OBJECTS_LEGACY) ${LDFLAGS} ${LIBIOTK} -o Loss

Optabs.o: Optabs.f
	${FC} ${FFLAGS} -c Optabs.f90

Optabs_para.o: Optabs.f
	${FC} ${FFLAGS} ${LDFLAGS_PARA} ${LIBIOTK} -c Optabs.f90 -o Optabs_para.o

Optabs_para_qe6.o: Optabs.f
	${FC} ${FFLAGS} ${LDFLAGS_PARA} -c Optabs_qe6.f90 -o Optabs_para_qe6.o

Optabs_legacy.o: Optabs.f
	${FC} ${FFLAGS} -c Optabs.f -o Optabs_legacy.o

PointR.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f90 -o PointR.o
	${FC} ${FFLAGS} -c PointR.f90

PointR_legacy.o: PointR.f
	${FC} ${FFLAGS} -c PointR_old.f90 -o PointR_legacy.o

paths.o: paths.f
	${FC} ${FFLAGS} -c paths.f90

gjel.o: gjel.f
	${FC} ${FFLAGS} -c gjel.f

clean:
	rm  *.o Loss fbz_check.dat modpointr.mod
