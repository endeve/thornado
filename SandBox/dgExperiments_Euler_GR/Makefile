include ../../../Build/Makefile_Build

all: ApplicationDriver RiemannProblem

ApplicationDriver: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule_Euler_GR.o \
	ApplicationDriver.o \
	$(LIBRARIES)

RiemannProblem: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	TimeSteppingModule_SSPRK.o \
	dgDiscretizationModule_Euler_GR.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o
	$(FLINKER) $(FLAGS) -o RiemannProblem_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	SlopeLimiterModule_Euler_GR.o \
	InitializationModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	TimeSteppingModule_SSPRK.o \
	dgDiscretizationModule_Euler_GR.o \
	RiemannProblemInitializer.o \
	RiemannProblem.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f RiemannProblem_$(MACHINE)
