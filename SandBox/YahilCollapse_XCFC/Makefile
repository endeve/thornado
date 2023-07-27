THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

POSEIDON_DIR ?= $(HOME)/poseidon
include $(POSEIDON_DIR)/Build/Make.Poseidon_Native

#DEFINES += -DTHORNADO_DEBUG

.DEFAULT_GOAL := all

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
  ReferenceElementModuleX.o \
  MeshModule.o \
  GravitySolutionModule_XCFC_Poseidon.o \
  Poseidon_UtilitiesModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  EquationOfStateModule_IDEAL.o \
  Euler_UtilitiesModule_Relativistic.o \
  UnitsModule.o \
  UtilitiesModule.o \
  Euler_SlopeLimiterModule_Relativistic_IDEAL.o \
  Euler_PositivityLimiterModule_Relativistic_IDEAL.o \
  InitializationModule_Relativistic.f90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GravitySolutionModule_XCFC_Poseidon.o \
  FluidFieldsModule.o \
  Euler_SlopeLimiterModule_Relativistic_TABLE.o \
  Euler_PositivityLimiterModule_Relativistic_TABLE.o \
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
  Poseidon_UtilitiesModule.o \
  TimersModule_Euler.o \
  ApplicationDriver.f90

