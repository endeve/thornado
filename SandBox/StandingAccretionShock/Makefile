include ../../../Build/Makefile_Build

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	GravitySolutionModule_Newtonian_PointMass_Beta.o \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	GravitySolutionModule_Newtonian_PointMass_Beta.o \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	Euler_TallyModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
