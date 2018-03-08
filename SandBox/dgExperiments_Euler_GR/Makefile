include ../../../Build/Makefile_Build

all: RiemannProblem StandingAccretionShock

RiemannProblem: \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o
	$(FLINKER) $(FLAGS) -o RiemannProblem_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o \
	$(LIBRARIES)

StandingAccretionShock: \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	StandingAccretionShock.o
	$(FLINKER) $(FLAGS) -o StandingAccretionShock_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	StandingAccretionShock.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f RiemannProblem_$(MACHINE)
	rm -f StandingAccretionShock_$(MACHINE)
	rm -f *.dat
	rm -rf *.dSYM
