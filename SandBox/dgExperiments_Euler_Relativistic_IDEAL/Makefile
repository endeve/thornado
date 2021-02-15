THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

DEFINES += -DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver #ComputePrimitiveTest SlopeLimiterDiagnostics

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
	rm -f  ComputePrimitiveTest_$(MACHINE)
	rm -f  SlopeLimiterDiagnostics_$(MACHINE)
	rm -rf *.dSYM

ApplicationDriver.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  EquationOfStateModule.o \
  ProgramHeaderModule.o \
  GeometryComputationModule.o \
  InitializationModule_Relativistic.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  InputOutputModuleHDF.o \
  FluidFieldsModule.o \
  GeometryFieldsModule.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  UnitsModule.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  ApplicationDriver.F90

InitializationModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  UnitsModule.o \
  QuadratureModule.o \
  PolynomialBasisModule_Lagrange.o \
  InitializationModule_Relativistic.f90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  TimersModule_Euler.o \
  TimeSteppingModule_SSPRK.F90
