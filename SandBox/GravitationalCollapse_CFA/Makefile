MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= $(HOME)/thornado
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: ApplicationDriver TestPressure

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

TestPressure: \
	$(weaklib) \
	$(thornado) \
	TestPressure.o
	$(FLINKER) $(FLAGS) -o TestPressure_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TestPressure.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM
	rm -f  TestPressure_$(MACHINE)
	rm -rf TestPressure_$(MACHINE).dSYM
