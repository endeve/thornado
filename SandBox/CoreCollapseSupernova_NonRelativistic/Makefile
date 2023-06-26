MICROPHYSICS           = WEAKLIB
GRAVITY_SOLVER         = POSEIDON_NEWTON
MOMENT_CLOSURE         = MINERBO
NEUTRINO_MATTER_SOLVER = EMAB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

POSEIDON_DIR ?= $(HOME)/poseidon
include $(POSEIDON_DIR)/Build/Make.Poseidon_Native

all: ApplicationDriver

ApplicationDriver: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	TwoMoment_PositivityLimiterModule_Old.o \
	ProgenitorModule.o \
	InitializationModule.o \
	TimeSteppingModule_CCSN.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	TwoMoment_PositivityLimiterModule_Old.o \
	ProgenitorModule.o \
	InitializationModule.o \
	TimeSteppingModule_CCSN.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -rf *.dSYM
