MICROPHYSICS = WEAKLIB
#MOMENT_CLOSURE = MINERBO
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

include ../../../Build/Makefile_Build

all: DeleptonizationWave

DeleptonizationWave: \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave_$(MACHINE) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f DeleptonizationWave_$(MACHINE)
