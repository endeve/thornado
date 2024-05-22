THORNAOD_DIR ?= ../../../

include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

.DEFAULT_GOAL := all

all: ApplicationDriver

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

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

InitializationModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryComputationModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  UnitsModule.o \
  UtilitiesModule.o \
  InitializationModule_Relativistic.f90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  TimersModule_Euler.o \
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
  Euler_dgDiscretizationModule.o \
  TimeSteppingModule_SSPRK.o \
  UnitsModule.o \
  Euler_TallyModule_Relativistic.o \
  TimersModule_Euler.o \
  ApplicationDriver.f90

