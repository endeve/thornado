MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

.DEFAULT_GOAL := all

all: TwoMomentClosure PrimitiveConserved TimeStepper ApplicationDriver

TwoMomentClosure: \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMomentClosure.o
	$(FLINKER) $(FLAGS) -o TwoMomentClosure_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMomentClosure.o \
	$(LIBRARIES)

PrimitiveConserved: \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	PrimitiveConserved.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	PrimitiveConserved.o \
	$(LIBRARIES)

TimeStepper: \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_OpacityModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	TimeStepper.o
	$(FLINKER) $(FLAGS) -o TimeStepper_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_OpacityModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	TimeStepper.o \
	$(LIBRARIES)

ApplicationDriver: \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_OpacityModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_OpacityModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  TwoMomentClosure_$(MACHINE)
	rm -rf TwoMomentClosure_$(MACHINE).dSYM
	rm -f  PrimitiveConserved_$(MACHINE)
	rm -rf PrimitiveConserved_$(MACHINE).dSYM
	rm -f  TimeStepper_$(MACHINE)
	rm -rf TimeStepper_$(MACHINE).dSYM
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

TwoMomentClosure.o: \
  KindModule.o \
  TwoMomentClosure.f90

PrimitiveConserved.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  PrimitiveConserved.f90

TimeStepper.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  GeometryComputationModule.o \
  GeometryComputationModuleE.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_TimeSteppingModule_OrderV.o \
  TimeStepper.f90

ApplicationDriver.o: \
  KindModule.o \
  ApplicationDriver.f90

InitializationModule.o: \
  KindModule.o \
  InitializationModule.F90

TwoMoment_TimersModule_OrderV.o: \
  KindModule.o \
  TwoMoment_TimersModule_OrderV.F90

TwoMoment_UtilitiesModule.o: \
  KindModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule.F90

TwoMoment_TroubledCellIndicatorModule.o: \
  KindModule.o \
  TwoMoment_TroubledCellIndicatorModule.F90

TwoMoment_SlopeLimiterModule_OrderV.o: \
  KindModule.o \
  TwoMoment_SlopeLimiterModule_OrderV.F90

TwoMoment_PositivityLimiterModule_OrderV.o: \
  KindModule.o \
  TwoMoment_PositivityLimiterModule_OrderV.F90

TwoMoment_OpacityModule_OrderV.o: \
  KindModule.o \
  TwoMoment_OpacityModule_OrderV.F90

TwoMoment_DiscretizationModule_Streaming_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.F90

TwoMoment_DiscretizationModule_Collisions_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.F90

TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.F90

TwoMoment_TimeSteppingModule_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.o \
  TwoMoment_TimeSteppingModule_OrderV.F90
