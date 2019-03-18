MICROPHYSICS = WEAKLIB

include ../../../Build/Makefile_Build

THORNADO_DIR ?= $(HOME)/dgHydro/thornado
WEAKLIB_DIR ?= $(HOME)/dgHydro/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: ComputeEigensystem_NuclearEOS

ComputeEigensystem_NuclearEOS: \
	$(thornado) \
	$(weaklib) \
	UtilitiesModule_NuclearEOS.o \
	ComputeEigensystem_NuclearEOS.o 
	$(FLINKER) $(FLAGS) -o ComputeEigensystem_NuclearEOS_$(MACHINE) \
	$(thornado) \
	$(weaklib) \
	UtilitiesModule_NuclearEOS.o \
	ComputeEigensystem_NuclearEOS.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ComputeEigensystem_NuclearEOS_$(MACHINE)
	rm -rf *.dSYM
