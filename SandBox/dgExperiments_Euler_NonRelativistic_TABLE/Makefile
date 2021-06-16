MICROPHYSICS   = WEAKLIB
#GRAVITY_SOLVER = POSEIDON_NEWTON

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

#DEFINES += $(PP)-DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
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

EigensystemTest: \
	$(weaklib) \
	$(thornado) \
        EigensystemTest.o
	$(FLINKER) $(FLAGS) -o EigensystemTest_$(MACHINE) \
        $(weaklib) \
        $(thornado) \
        EigensystemTest.o \
        $(LIBRARIES)

ComputeLeftState_NuclearEOS: \
	$(weaklib) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeLeftState_NuclearEOS.o
	$(FLINKER) $(FLAGS) -o ComputeLeftState_NuclearEOS_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeLeftState_NuclearEOS.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f ComputeEigensystem_NuclearEOS_$(MACHINE)
	rm -f EigensystemTest_$(MACHINE)
	rm -f ComputeLeftState_NuclearEOS_$(MACHINE)
	rm -rf *.dSYM
