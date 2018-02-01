MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

include ../../../Build/Makefile_Build

all: StreamingSineWave

StreamingSineWave: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleE_Lagrange.o \
	GeometryComputationModuleE_Beta.o \
	ReferenceElementModule_Lagrange.o \
	TimeSteppingModule_IMEX_RK.o \
	ClosureModule_M1.o \
	MomentEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule.o \
	dgDiscretizationModule_Collisions.o \
	StreamingSineWave.o
	$(FLINKER) $(FLAGS) -o StreamingSineWave_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleE_Lagrange.o \
	GeometryComputationModuleE_Beta.o \
	ReferenceElementModule_Lagrange.o \
	TimeSteppingModule_IMEX_RK.o \
	ClosureModule_M1.o \
	MomentEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule.o \
	dgDiscretizationModule_Collisions.o \
	StreamingSineWave.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f StreamingSineWave_$(MACHINE)
