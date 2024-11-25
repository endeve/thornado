include ../../../Build/Makefile_Build

POSEIDON_DIR ?= $(HOME)/poseidon
include $(POSEIDON_DIR)/Build/Make.Poseidon_Native

all: ApplicationDriver

ApplicationDriver: \
	$(POSEIDON_o) \
	$(thornado) \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(POSEIDON_o) \
	$(thornado) \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

Poseidon_UtilitiesModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  EquationOfStateModule_IDEAL.o \
  EquationOfStateModule_TABLE.o \
  Poseidon_UtilitiesModule.F90

InitializationModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GravitySolutionModule_XCFC_Poseidon.o \
  GeometryComputationModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  UnitsModule.o \
  UtilitiesModule.o \
  Poseidon_UtilitiesModule.o \
  InitializationModule_Relativistic.f90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  TimersModule_Euler.o \
  Poseidon_UtilitiesModule.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.f90

ApplicationDriver.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  EquationOfStateModule.o \
  ProgramHeaderModule.o \
  InitializationModule_Relativistic.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  InputOutputModuleHDF.o \
  FluidFieldsModule.o \
  GeometryFieldsModule.o \
  GravitySolutionModule_XCFC_Poseidon.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  UnitsModule.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  Poseidon_UtilitiesModule.o \
  ApplicationDriver.f90

