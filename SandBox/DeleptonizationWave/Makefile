MICROPHYSICS = WEAKLIB
#MOMENT_CLOSURE = MINERBO
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: DeleptonizationWave

DeleptonizationWave: \
	$(thornado) \
	$(weaklib) \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave_$(MACHINE) \
	$(thornado) \
	$(weaklib) \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f DeleptonizationWave_$(MACHINE)
