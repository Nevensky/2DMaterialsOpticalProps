FC=ifort
# FFLAGS=-Ofast
OMP=-qopenmp
FFLAGS_BASE=-qmkl -O2 -fpp -g -traceback -check bounds ${OMP}
FFLAGS_DEBUG=-warn all -check all -debug all

# Base names of the source files of Modules
BASE_NAMES = constants statistics utility matrix_inverse io_qe io_xml io_wannier90 local_field PointR_mkl brillouin_zone qe_gksort propagators vertices Ladder

# Source files and object files
SRC = $(addprefix Modules/, $(addsuffix .f90, $(BASE_NAMES)))
OBJECTS = $(addprefix Modules/, $(addsuffix .o, $(BASE_NAMES)))

# Ensures make will run the commands associated with these targets even if there are files named clean.f90 or all.f90 in the directory.
.PHONY: clean all

all: FFLAGS=${FFLAGS_BASE} 
all: $(OBJECTS)

debug: FFLAGS=${FFLAGS_BASE} ${FFLAGS_DEBUG} 
debug: $(OBJECTS)

Modules/%.o: Modules/%.f90
	$(FC) $(FFLAGS) -c $< -o $@ -module Modules

clean:
	rm -f Modules/*.o Modules/*.mod Modules/*.genmod