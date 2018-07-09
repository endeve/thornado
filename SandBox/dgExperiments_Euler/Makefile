include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	EulerEquationsUtilitiesModule_Beta.o \
	CharacteristicDecompositionModule.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler.o \
	PositivityLimiterModule_Euler.o \
	dgDiscretizationModule_Euler.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	EulerEquationsUtilitiesModule_Beta.o \
	CharacteristicDecompositionModule.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler.o \
	PositivityLimiterModule_Euler.o \
	dgDiscretizationModule_Euler.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
