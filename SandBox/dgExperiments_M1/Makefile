include ../../../Build/Makefile_Build

all: rhsTest

rhsTest: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleE_Lagrange.o \
	GeometryComputationModuleE_Beta.o \
	ReferenceElementModule_Lagrange.o \
	MomentEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule.o \
	rhsTest.o
	$(FLINKER) $(FLAGS) -o rhsTest_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleE_Lagrange.o \
	GeometryComputationModuleE_Beta.o \
	ReferenceElementModule_Lagrange.o \
	MomentEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	dgDiscretizationModule.o \
	rhsTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f rhsTest_$(MACHINE)
