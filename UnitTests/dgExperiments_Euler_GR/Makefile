include ../../../Build/Makefile_Build

all: rhsTest_GR

rhsTest_GR: \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	rhsTest_GR.o
	$(FLINKER) $(FLAGS) -o rhsTest_GR_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX.o \
	ReferenceElementModuleX_Lagrange.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	rhsTest_GR.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f rhsTest_GR_$(MACHINE)
