include ../../../Build/Makefile_Build

all: ApplicationDriver RiemannProblem

ApplicationDriver: \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	ApplicationDriver.o \
	$(LIBRARIES)

RiemannProblem: \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	PositivityLimiterModule.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
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
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f RiemannProblem_$(MACHINE)
	rm -f *.dat
	rm -rf *.dSYM
