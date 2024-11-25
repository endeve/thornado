MICROPHYSICS = WEAKLIB
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_Flavor

.DEFAULT_GOAL := all

all: Decoherence

Decoherence: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	Decoherence.o
	$(FLINKER) $(FLAGS) -o Decoherence_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	Decoherence.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Decoherence_$(MACHINE)

Decoherence.o: \
  Decoherence.f90
