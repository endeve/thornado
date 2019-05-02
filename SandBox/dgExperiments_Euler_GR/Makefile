include ../../../Build/Makefile_Build

all: ApplicationDriver ComputePrimitiveTest

ApplicationDriver: \
	$(thornado) \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ComputePrimitiveTest: \
	$(thornado) \
	ComputePrimitiveTest.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
	$(thornado) \
	ComputePrimitiveTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f ComputePrimitiveTest_$(MACHINE)
	rm -f *.dat
	rm -rf *.dSYM
