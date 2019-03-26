MICROPHYSICS   = WEAKLIB
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

all: Relaxation

Relaxation: \
	$(thornado) \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_beta.o \
	InitializationModule.o \
	Relaxation.o
	$(FLINKER) $(FLAGS) -o Relaxation_$(MACHINE) \
	$(thornado) \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_beta.o \
	InitializationModule.o \
	Relaxation.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Relaxation_$(MACHINE)
