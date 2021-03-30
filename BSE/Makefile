# Makefile for Vito
# FC= gfortran
FC = ifort
DEBUG = -g -traceback -fpe0 -check bounds
FFLAGS = -mkl -Ofast ${DEBUG}
OBJECTS = Photonew.f90 gjel.f
# OBJECTS_LEGACY = Photonew.f gjel.f


all: $(OBJECTS)
	${FC} ${FFLAGS} $(OBJECTS) -o Photon


clean:
	rm  *.o Photon