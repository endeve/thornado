MICROPHYSICS   = WEAKLIB
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: Relaxation

Relaxation: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	Relaxation.o
	$(FLINKER) $(FLAGS) -o Relaxation_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	Relaxation.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Relaxation_$(MACHINE)
