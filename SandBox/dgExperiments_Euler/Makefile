include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	Euler_UtilitiesModule.o \
	Euler_CharacteristicDecompositionModule.o \
	InitializationModule.o \
	Euler_BoundaryConditionsModule.o \
	Euler_SlopeLimiterModule.o \
	Euler_PositivityLimiterModule.o \
	dgDiscretizationModule_Euler.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	Euler_UtilitiesModule.o \
	Euler_CharacteristicDecompositionModule.o \
	InitializationModule.o \
	Euler_BoundaryConditionsModule.o \
	Euler_SlopeLimiterModule.o \
	Euler_PositivityLimiterModule.o \
	dgDiscretizationModule_Euler.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
