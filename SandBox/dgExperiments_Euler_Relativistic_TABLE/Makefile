MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= ../../../weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

#DEFINES += -DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver \
	ComputePrimitiveTest \
	ComputePrimitiveTest_Single

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ComputePrimitiveTest: \
	$(weaklib) \
	$(thornado) \
	ComputePrimitiveTest.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ComputePrimitiveTest.o \
	$(LIBRARIES)

ComputePrimitiveTest_Single: \
	$(weaklib) \
	$(thornado) \
	ComputePrimitiveTest_Single.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_Single_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ComputePrimitiveTest_Single.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM
	rm -f  ComputePrimitiveTest_$(MACHINE)
	rm -rf ComputePrimitiveTest_$(MACHINE).dSYM

InitializationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  UnitsModule.o \
  InitializationModule.F90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
  TimersModule_Euler.o \
  Euler_dgDiscretizationModule.o \
  Euler_TallyModule_Relativistic.o \
  TimeSteppingModule_SSPRK.F90

ApplicationDriver.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  EquationOfStateModule.o \
  EquationOfStateModule_TABLE.o \
  GeometryFieldsModule.o \
  GeometryComputationModule.o \
  FluidFieldsModule.o \
  InitializationModule.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
  Euler_UtilitiesModule_Relativistic.o \
  InputOutputModuleHDF.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  UnitsModule.o \
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
  ProgramInitializationModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  EquationOfStateModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  Euler_ErrorModule.o \
  UnitsModule.o \
  ComputePrimitiveTest_Single.F90
