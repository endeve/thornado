include ../../../Build/Makefile_Build

all: ApplicationDriver ComputePrimitiveTest

ApplicationDriver: \
	$(thornado) \
	Euler_GR_PositivityLimiterModule.o \
	Euler_GR_UtilitiesModule.o \
	Euler_GR_CharacteristicDecompositionModule.o \
	Euler_GR_BoundaryConditionsModule.o \
	Euler_GR_SlopeLimiterModule.o \
	InitializationModule_GR.o \
	Euler_GR_dgDiscretizationModule.o \
	TimeSteppingModule_SSPRK.o \
	Euler_GR_TallyModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	Euler_GR_PositivityLimiterModule.o \
	Euler_GR_UtilitiesModule.o \
	Euler_GR_CharacteristicDecompositionModule.o \
	Euler_GR_BoundaryConditionsModule.o \
	Euler_GR_SlopeLimiterModule.o \
	InitializationModule_GR.o \
	Euler_GR_dgDiscretizationModule.o \
	TimeSteppingModule_SSPRK.o \
	Euler_GR_TallyModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ComputePrimitiveTest: \
	$(thornado) \
	Euler_GR_UtilitiesModule.o \
	ComputePrimitiveTest.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
	$(thornado) \
	Euler_GR_UtilitiesModule.o \
	ComputePrimitiveTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f ComputePrimitiveTest_$(MACHINE)
	rm -f *.dat
	rm -rf *.dSYM
