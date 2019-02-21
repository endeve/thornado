MICROPHYSICS   = WEAKLIB
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

all: Relaxation

Relaxation: \
	$(thornado) \
	InitializationModule.o \
	Relaxation.o
	$(FLINKER) $(FLAGS) -o Relaxation_$(MACHINE) \
	$(thornado) \
	InitializationModule.o \
	Relaxation.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Relaxation_$(MACHINE)
