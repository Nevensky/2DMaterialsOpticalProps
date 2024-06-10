# Makefile for Vito
DEBUG=-warn all -check all -debug all
FC=ifort
#FFLAGS=-Ofast
OMP=-qopenmp
FFLAGS=-qmkl -O2 -fpp -g -traceback -check bounds ${OMP} #${DEBUG}
SRC = constants.f90 statistics.f90 utility.f90 matrix_inverse.f90 io_qe.f90 io_xml.f90 io_wannier90.f90  local_field.f90 PointR_mkl.f90 brillouin_zone.f90 qe_gksort.f90 propagators.f90 vertices.f90 Ladder.f90
OBJECTS = constants statistics utility matrix_inverse io_qe io_xml io_wannier90 local_field PointR brillouin_zone qe_gksort propagators vertices Ladder


.PHONY: clean

all: Modules

Modules: $(OBJECTS)

constants : constants.f90
	${FC} ${FFLAGS} -c constants.f90

statistics : statistics.f90
	${FC} ${FFLAGS} -c statistics.f90

io_qe : io_qe.f90
	${FC} ${FFLAGS} -c io_qe.f90

io_xml : io_qe.f90
	${FC} ${FFLAGS} -c io_xml.f90

io_wannier90 : io_wannier90.f90
	${FC} ${FFLAGS} -c io_wannier90.f90

utility : utility.f90
	${FC} ${FFLAGS} -c utility.f90

local_field : local_field.f90
	${FC} ${FFLAGS} -c local_field.f90

PointR : PointR_mkl.f90
	${FC} ${FFLAGS} -c PointR_mkl.f90

matrix_inverse : matrix_inverse.f90
	${FC} -qmkl -Ofast -fpp -c matrix_inverse.f90

brillouin_zone : brillouin_zone.f90
	${FC} ${FFLAGS} -c brillouin_zone.f90

qe_gksort : qe_gksort.f90
	${FC} ${FFLAGS} -c qe_gksort.f90

propagators : propagators.f90
	${FC} ${FFLAGS} -c propagators.f90

vertices : vertices.f90
	${FC} ${FFLAGS} -c vertices.f90

Ladder : Ladder.f90
	${FC} ${FFLAGS} -c Ladder.f90

clean:
	rm  *.o *__genmod.f90 *__genmod.mod *.mod
