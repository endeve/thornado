include ../../../Build/Makefile_Build

all: ApplicationDriver SlopeLimiterDiagnostics #ComputePrimitiveTest

ApplicationDriver: \
	$(thornado) \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	AccretionShockDiagnosticsModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	AccretionShockDiagnosticsModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

#ComputePrimitiveTest: \
#	$(thornado) \
#	ComputePrimitiveTest.o
#	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
#	$(thornado) \
#	ComputePrimitiveTest.o \
#	$(LIBRARIES)

SlopeLimiterDiagnostics: \
	$(thornado) \
	InitializationModule_Relativistic.o \
	SlopeLimiterDiagnostics.o
	$(FLINKER) $(FLAGS) -o SlopeLimiterDiagnostics_$(MACHINE) \
	$(thornado) \
	InitializationModule_Relativistic.o \
	SlopeLimiterDiagnostics.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM
#	rm -f  ComputePrimitiveTest_$(MACHINE)
#	rm -rf ComputePrimitiveTest_$(MACHINE).dSYM
	rm -f  SlopeLimiterDiagnostics_$(MACHINE)
	rm -rf SlopeLimiterDiagnostics_$(MACHINE).dSYM
