include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta.o \
        InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	ReferenceElementModuleX_Lagrange.o \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta.o \
        InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
