MOMENT_CLOSURE = MINERBO

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_GR

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	Euler_BoundaryConditionsModule_Relativistic.o \
	TwoMoment_UtilitiesModule_Relativistic.o \
	TwoMoment_DiscretizationModule_Streaming_Relativistic.o \
        TwoMoment_PositivityLimiterModule_Relativistic.o \
	TwoMoment_TimeSteppingModule_Relativistic.o \
	InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	Euler_BoundaryConditionsModule_Relativistic.o \
	TwoMoment_UtilitiesModule_Relativistic.o \
	TwoMoment_DiscretizationModule_Streaming_Relativistic.o \
        TwoMoment_PositivityLimiterModule_Relativistic.o \
	TwoMoment_TimeSteppingModule_Relativistic.o \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

ApplicationDriver.o: \
  KindModule.o \
  ApplicationDriver.f90

InitializationModule.o: \
  KindModule.o \
  InitializationModule.F90

TwoMoment_UtilitiesModule_Relativistic.o: \
  KindModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule_Relativistic.F90

TwoMoment_DiscretizationModule_Streaming_Relativistic.o: \
  KindModule.o \
  TwoMoment_DiscretizationModule_Streaming_Relativistic.F90

TwoMoment_TimeSteppingModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_PositivityLimiterModule_Relativistic.o \
  TwoMoment_TimeSteppingModule_Relativistic.F90

TwoMoment_PositivityLimiterModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule_Lagrange.o \
  ReferenceElementModuleX_Lagrange.o \
  LinearAlgebraModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  Euler_UtilitiesModule_Relativistic.o \
  TwoMoment_PositivityLimiterModule_Relativistic.F90


