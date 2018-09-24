include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	PositivityLimiterModule_Euler_GR.o \
	EulerEquationsUtilitiesModule_Beta_GR.o \
	CharacteristicDecompositionModule_GR.o \
	BoundaryConditionsModule_Beta.o \
	SlopeLimiterModule_Euler_GR.o \
	DataFileReaderModule.o \
	InitializationModule_GR.o \
	dgDiscretizationModule_Euler_GR.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f *.dat
	rm -rf *.dSYM
