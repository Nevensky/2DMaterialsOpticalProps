# Makefile for Vito
DEBUG=-warn all -check all -debug all -g
FC=ifort
#FFLAGS=-Ofast
FFLAGS=-qmkl -Ofast -fpp #${DEBUG}
SRC = constants.f90 io_qe.f90 PointR_mkl.f90 matrix_inverse.f90 brillouin_zone.f90 Gamma.f90 RPA.f90 W.f90 Ladder.f90
OBJECTS = constants io_qe PointR brillouin_zone Gamma RPA W Ladder


.PHONY: clean

all: Modules

Modules: $(OBJECTS)

constants : constants.f90
	${FC} ${FFLAGS} -c constants.f90

io_qe : io_qe.f90
	${FC} ${FFLAGS} -c io_qe.f90

PointR : PointR_mkl.f90
	${FC} ${FFLAGS} -c PointR_mkl.f90

matrix_inverse : matrix_inverse.f90
	${FC} ${FFLAGS} -c matrix_inverse.f90

brillouin_zone : brillouin_zone.f90
	${FC} ${FFLAGS} -c  brillouin_zone.f90

Gamma : Gamma.f90
	${FC} ${FFLAGS} -c Gamma.f90

RPA : RPA.f90
	${FC} ${FFLAGS} -c RPA.f90

W : W.f90
	${FC} ${FFLAGS} -c W.f90

Ladder : Ladder.f90
	${FC} ${FFLAGS} -c Ladder.f90

clean:
	rm  *.o *__genmod.f90 *__genmod.mod