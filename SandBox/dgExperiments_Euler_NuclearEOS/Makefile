MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: ApplicationDriver

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ComputeEigensystem_NuclearEOS: \
	$(weaklib) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeEigensystem_NuclearEOS.o
	$(FLINKER) $(FLAGS) -o ComputeEigensystem_NuclearEOS_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeEigensystem_NuclearEOS.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f ComputeEigensystem_NuclearEOS_$(MACHINE)
	rm -rf *.dSYM
