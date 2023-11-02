MOMENT_CLOSURE  = MINERBO

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_FMC

DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	TwoMoment_FieldsModule_FMC.o \
	TwoMoment_InputOutputModule_FMC.o \
	TwoMoment_UtilitiesModule_FMC.o \
  TwoMoment_PositivityLimiterModule_FMC.o \
  TwoMoment_DiscretizationModule_Streaming_FMC.o \
  TwoMoment_DiscretizationModule_Collisions_FMC.o \
  TwoMoment_OpacityModule_FMC.o \
  TwoMoment_TimeSteppingModule_FMC.o \
  InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	TwoMoment_FieldsModule_FMC.o \
	TwoMoment_InputOutputModule_FMC.o \
	TwoMoment_UtilitiesModule_FMC.o \
  TwoMoment_PositivityLimiterModule_FMC.o \
  TwoMoment_DiscretizationModule_Streaming_FMC.o \
  TwoMoment_DiscretizationModule_Collisions_FMC.o \
  TwoMoment_OpacityModule_FMC.o \
  TwoMoment_TimeSteppingModule_FMC.o \
  InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld *.lst *.acc.s

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

ApplicationDriver.o: \
  ApplicationDriver.F90

InitializationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleZ.o \
  MeshModule.o \
  TwoMoment_FieldsModule_FMC.o \
  InitializationModule.F90

TwoMoment_FieldsModule_FMC.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TwoMoment_FieldsModule_FMC.F90

TwoMoment_InputOutputModule_FMC.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule.o \
  MeshModule.o \
  InputOutputUtilitiesModule.o \
  TwoMoment_FieldsModule_FMC.o \
  TwoMoment_InputOutputModule_FMC.F90

TwoMoment_UtilitiesModule_FMC.o: \
  KindModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule_FMC.F90

TwoMoment_TimeSteppingModule_FMC.o: \
  TwoMoment_PositivityLimiterModule_FMC.o \
  TwoMoment_DiscretizationModule_Streaming_FMC.o \
  TwoMoment_TimeSteppingModule_FMC.F90

TwoMoment_DiscretizationModule_Streaming_FMC.o: \
  TwoMoment_DiscretizationModule_Streaming_FMC.F90

TwoMoment_DiscretizationModule_Collisions_FMC.o: \
  TwoMoment_DiscretizationModule_Collisions_FMC.F90

TwoMoment_PositivityLimiterModule_FMC.o: \
  TwoMoment_PositivityLimiterModule_FMC.F90

TwoMoment_OpacityModule_FMC.o: \
  TwoMoment_OpacityModule_FMC.F90
