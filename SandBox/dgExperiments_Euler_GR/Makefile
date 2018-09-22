include ../../../Build/Makefile_Build

all: ApplicationDriver RiemannProblem StandingAccretionShock ComputePrimitiveTest

ApplicationDriver: \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

RiemannProblem: \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o
	$(FLINKER) $(FLAGS) -o RiemannProblem_$(MACHINE) \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o \
	$(LIBRARIES)

StandingAccretionShock: \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	StandingAccretionShock.o
	$(FLINKER) $(FLAGS) -o StandingAccretionShock_$(MACHINE) \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	StandingAccretionShock.o \
	$(LIBRARIES)

ComputePrimitiveTest: \
	$(thornado) \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	ComputePrimitiveTest.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
	$(thornado) \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	ComputePrimitiveTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f RiemannProblem_$(MACHINE)
	rm -f StandingAccretionShock_$(MACHINE)
	rm -f ComputePrimitiveTest_$(MACHINE)
	rm -f *.dat
	rm -rf *.dSYM
