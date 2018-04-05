include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta.o \
        InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	GeometryComputationModule_Beta.o \
	EulerEquationsUtilitiesModule_Beta.o \
        InitializationModule.o \
	BoundaryConditionsModule_Beta.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
