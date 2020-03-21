# Makefile for Vito
CC= gcc
FC= gfortran
CFLAGS= -lrt -lm -O3
FFLAGS=	${CFLAGS}
OBJECTS = Photonew.f gjel.f
compile: $(OBJECTS)
	${FC} ${FFLAGS} $(OBJECTS) -o Foton
