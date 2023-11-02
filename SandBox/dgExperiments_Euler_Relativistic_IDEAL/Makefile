THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

#DEFINES += -DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver \
	ComputePrimitiveTest \
	ComputePrimitiveTest_Single \
	SlopeLimiterDiagnostics

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

ComputePrimitiveTest_Single: \
	$(thornado) \
	ComputePrimitiveTest_Single.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_Single_$(MACHINE) \
	$(thornado) \
	ComputePrimitiveTest_Single.o \
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
	rm -f  ComputePrimitiveTest_Single_$(MACHINE)
	rm -f  SlopeLimiterDiagnostics_$(MACHINE)
	rm -rf *.dSYM

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
  InitializationModule_Relativistic.F90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  TimersModule_Euler.o \
  Euler_dgDiscretizationModule.o \
  Euler_TallyModule_Relativistic.o \
  TimeSteppingModule_SSPRK.F90

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

ComputePrimitiveTest.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  GeometryComputationModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  Euler_ErrorModule.o \
  TimersModule_Euler.o \
  ComputePrimitiveTest.F90

ComputePrimitiveTest_Single.o: \
  KindModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  Euler_ErrorModule.o \
  UnitsModule.o \
  ComputePrimitiveTest_Single.F90

SlopeLimiterDiagnostics.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  EquationOfStateModule.o \
  ProgramHeaderModule.o \
  GeometryComputationModule.o \
  InitializationModule_Relativistic.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  InputOutputModuleHDF.o \
  FluidFieldsModule.o \
  GeometryFieldsModule.o \
  SlopeLimiterDiagnostics.f90
