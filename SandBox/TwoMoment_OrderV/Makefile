MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

.DEFAULT_GOAL := all

all: PrimitiveConserved TimeStepper

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
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	TimeStepper.o
	$(FLINKER) $(FLAGS) -o TimeStepper_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	TimeStepper.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f PrimitiveConserved_$(MACHINE)
	rm -f TimeStepper_$(MACHINE)

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

TwoMoment_UtilitiesModule.o: \
  KindModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule.F90

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
