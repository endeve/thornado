MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= $(HOME)/thornado
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: ApplicationDriver TestPressure

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

TestPressure: \
	$(weaklib) \
	$(thornado) \
	TestPressure.o
	$(FLINKER) $(FLAGS) -o TestPressure_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TestPressure.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM
	rm -f  TestPressure_$(MACHINE)
	rm -rf TestPressure_$(MACHINE).dSYM

ApplicationDriver.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  EquationOfStateModule.o \
  EquationOfStateModule_TABLE.o \
  ProgramHeaderModule.o \
  InitializationModule_Relativistic.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
  Euler_UtilitiesModule_Relativistic.o \
  InputOutputModuleHDF.o \
  FluidFieldsModule.o \
  GeometryFieldsModule.o \
  GravitySolutionModule_CFA_Poseidon.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  UnitsModule.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  Poseidon_UtilitiesModule.o \
  ApplicationDriver.f90

ProgenitorModule.o: \
  KindModule.o \
  UnitsModule.o \
  UtilitiesModule.o \
  ProgenitorModule.f90

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
  MeshModule.o \
  UtilitiesModule.o \
  GravitySolutionModule_CFA_Poseidon.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  Poseidon_UtilitiesModule.o \
  EquationOfStateModule.o \
  ProgenitorModule.o \
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
