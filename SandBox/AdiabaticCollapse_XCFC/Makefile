MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= $(HOME)/thornado
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

POSEIDON_DIR ?= $(HOME)/poseidon
include $(POSEIDON_DIR)/Build/Make.Poseidon_Native

.DEFAULT_GOAL := all

all: ApplicationDriver# ComputePrimitiveTest TestPressure

ApplicationDriver: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ProgenitorModule.o \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ProgenitorModule.o \
	Poseidon_UtilitiesModule.o \
	InitializationModule_Relativistic.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ComputePrimitiveTest: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ComputePrimitiveTest.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ComputePrimitiveTest.o \
	$(LIBRARIES)

TestPressure: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	TestPressure.o
	$(FLINKER) $(FLAGS) -o TestPressure_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
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
	rm -f  ComputePrimitiveTest_$(MACHINE)
	rm -rf ComputePrimitiveTest_$(MACHINE).dSYM

ProgenitorModule.o: \
  KindModule.o \
  UnitsModule.o \
  UtilitiesModule.o \
  ProgenitorModule.f90

Poseidon_UtilitiesModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  UnitsModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  EquationOfStateModule.o \
  Euler_ErrorModule.o \
  TimersModule_Euler.o \
  Poseidon_UtilitiesModule.F90

InitializationModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  MeshModule.o \
  UtilitiesModule.o \
  GravitySolutionModule_XCFC_Poseidon.o \
  Poseidon_UtilitiesModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
  Euler_UtilitiesModule_Relativistic.o \
  EquationOfStateModule.o \
  ProgenitorModule.o \
  InitializationModule_Relativistic.f90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GravitySolutionModule_XCFC_Poseidon.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  Poseidon_UtilitiesModule.o \
  TimersModule_Euler.o \
  Euler_dgDiscretizationModule.o \
  Euler_TallyModule_Relativistic.o \
  TimeSteppingModule_SSPRK.f90

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
  GravitySolutionModule_XCFC_Poseidon.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  UnitsModule.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  Poseidon_UtilitiesModule.o \
  ApplicationDriver.f90

ComputePrimitiveTest.o: \
  KindModule.o \
  UnitsModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule.o \
  EquationOfStateModule_TABLE.o \
  Euler_UtilitiesModule_Relativistic.o \
  ComputePrimitiveTest.f90

TestPressure.o: \
  KindModule.o \
  EquationOfStateModule.o \
  UnitsModule.o \
  TestPressure.f90
