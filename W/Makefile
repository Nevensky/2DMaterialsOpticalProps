# Makefile for Vito
CC= icc
#FC= gfortran
# FC= mpiifort
FC= ifort
# LDFLAGS= -lm -L. -liotk -mkl
LDFLAGS_PARA = -qopenmp -L. -liotk -mkl # -Wl,-stack_size,0x100000000 
LDFLAGS= -L. -liotk -mkl 
PARA = -qopenmp
DEBUG = -g -traceback -fpe0 -fp-model precise -fp-stack-check -heap-arrays #-debug all -debug-parameters all -DD -C -warn all -check all
# PARA2= -parallel -opt-report-phase=par -opt-report:5
FFLAGS = -Ofast ${DEBUG}
#OMPFLAGS = -qopenmp -debug extended -heap-array
SRC = PointR.f Sloss.f paths.f gjel.f 
OBJECTS = PointR.o Sloss.o paths.o gjel.o 
OBJECTS_LEGACY = Sloss_legacy.o PointR_legacy.o paths.o gjel.o 
OBJECTS_PARA = PointR.o Sloss_para.o paths.o gjel.o  # PointR_mkl.o
OBJECTS_QE6 = PointR.o Sloss_para_qe6.o paths.o gjel.o  # PointR_mkl.o


.PHONY: clean

all: Loss

qe6: Loss_qe6

serial: Loss_serial

legacy: Loss_legacy

Loss: $(OBJECTS_PARA)
	${FC} -I. $(OBJECTS_PARA) ${LDFLAGS_PARA} -o Loss

Loss_serial: $(OBJECTS)
	${FC} -I. $(OBJECTS) ${LDFLAGS} -o Loss

Loss_qe6: $(OBJECTS_QE6)
	${FC} -I. $(OBJECTS_QE6) ${LDFLAGS_PARA} -o Loss

Loss_legacy: $(OBJECTS_LEGACY)
	${FC} -I. $(OBJECTS_LEGACY) ${LDFLAGS} -o Loss

Sloss.o: Sloss.f
	${FC} ${FFLAGS}  -c Sloss.f90

Sloss_para.o: Sloss.f
	${FC} ${FFLAGS} ${PARA} -c Sloss.f90 -o Sloss_para.o

Sloss_para_qe6.o: Sloss.f
	${FC} ${FFLAGS} ${PARA} -c Sloss_qe6.f90 -o Sloss_para_qe6.o

Sloss_legacy.o: Sloss.f
	${FC} ${FFLAGS} -c Sloss.f -o Sloss_legacy.o

PointR.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f90 -o PointR.o
	${FC} ${FFLAGS} -c PointR.f90

PointR_mkl.o: PointR.f
	${FC} ${FFLAGS} -c PointR_mkl.f90 -o PointR_mkl.o
	${FC} ${FFLAGS} -c PointR_mkl.f90

PointR_legacy.o: PointR.f
	${FC} ${FFLAGS} -c PointR.f -o PointR_legacy.o

paths.o: paths.f
	${FC} ${FFLAGS} -c paths.f90

gjel.o: gjel.f
	${FC} ${FFLAGS} -c gjel.f

clean:
	rm  *.o Loss fbz_check.dat modpointr.mod
