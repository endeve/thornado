MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
#MOMENT_CLOSURE = KERSHAW_BL

include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleE_Lagrange.o \
	GeometryComputationModuleE_Beta.o \
	ReferenceElementModule_Lagrange.o \
	PositivityLimiterModule.o \
	TimeSteppingModule_IMEX_RK.o \
	ClosureModule_M1.o \
	MomentEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule.o \
	dgDiscretizationModule_Collisions.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleE_Lagrange.o \
	GeometryComputationModuleE_Beta.o \
	ReferenceElementModule_Lagrange.o \
	PositivityLimiterModule.o \
	TimeSteppingModule_IMEX_RK.o \
	ClosureModule_M1.o \
	MomentEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule.o \
	dgDiscretizationModule_Collisions.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
