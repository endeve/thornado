include ../../../Build/Makefile_Build

all: rhsTest

rhsTest: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	EulerEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	dgDiscretizationModule_Euler.o \
	rhsTest.o
	$(FLINKER) $(FLAGS) -o rhsTest_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	EulerEquationsUtilitiesModule_Beta.o \
	InitializationModule.o \
	dgDiscretizationModule_Euler.o \
	rhsTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f rhsTest_$(MACHINE)
