MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= ../../../weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

DEFINES += -DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver

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

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

ApplicationDriver.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  EquationOfStateModule.o \
  EquationOfStateModule_TABLE.o \
  ProgramHeaderModule.o \
  GeometryComputationModule.o \
  InitializationModule.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
  Euler_UtilitiesModule_Relativistic.o \
  InputOutputModuleHDF.o \
  FluidFieldsModule.o \
  GeometryFieldsModule.o \
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  UnitsModule.o \
  ApplicationDriver.F90

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
  InitializationModule.f90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
  TimersModule_Euler.o \
  Euler_dgDiscretizationModule.o \
  Euler_TallyModule_Relativistic.o \
  TimeSteppingModule_SSPRK.f90
