MICROPHYSICS   = WEAKLIB
#GRAVITY_SOLVER = POSEIDON_NEWTON

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

ifeq ($(USE_GPU),FALSE)
  POSEIDON_DIR ?= $(HOME)/poseidon
  include $(POSEIDON_DIR)/Build/Make.Poseidon_Native
endif

#DEFINES += $(PP)-DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ProgenitorModule.o \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ProgenitorModule.o \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ComputeEigensystem_NuclearEOS: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeEigensystem_NuclearEOS.o
	$(FLINKER) $(FLAGS) -o ComputeEigensystem_NuclearEOS_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeEigensystem_NuclearEOS.o \
	$(LIBRARIES)

EigensystemTest: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
        EigensystemTest.o
	$(FLINKER) $(FLAGS) -o EigensystemTest_$(MACHINE) \
        $(weaklib) \
	$(POSEIDON_o) \
        $(thornado) \
        EigensystemTest.o \
        $(LIBRARIES)

ComputeLeftState_NuclearEOS: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	UtilitiesModule_NuclearEOS.o \
	ComputeLeftState_NuclearEOS.o
	$(FLINKER) $(FLAGS) -o ComputeLeftState_NuclearEOS_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
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
