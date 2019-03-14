MICROPHYSICS = WEAKLIB
#MOMENT_CLOSURE = MINERBO
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: DeleptonizationWave

DeleptonizationWave: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f DeleptonizationWave_$(MACHINE)

include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies
DeleptonizationWave.o: \
   KindModule.o \
   ProgramHeaderModule.o \
   UnitsModule.o \
   ProgramInitializationModule.o \
   ReferenceElementModuleX.o \
   ReferenceElementModuleX_Lagrange.o \
   ReferenceElementModuleE.o \
   ReferenceElementModuleE_Lagrange.o \
   ReferenceElementModule.o \
   ReferenceElementModule_Lagrange.o \
   GeometryFieldsModule.o \
   GeometryComputationModule.o \
   GeometryFieldsModuleE.o \
   GeometryComputationModuleE.o \
   FluidFieldsModule.o \
   RadiationFieldsModule.o \
   EquationOfStateModule_TABLE.o \
   OpacityModule_TABLE.o \
   NeutrinoOpacitiesModule.o \
   NeutrinoOpacitiesComputationModule.o \
   TimeSteppingModule_Castro.o \
   InitializationModule.o \
   InputOutputModuleHDF.o \
   TwoMoment_ClosureModule.o \
   TwoMoment_PositivityLimiterModule.o \
   DeleptonizationWave.f90

InitializationModule.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  TwoMoment_UtilitiesModule.o \
  InitializationModule.f90

TimeSteppingModule_Castro.o: \
   KindModule.o \
   ProgramHeaderModule.o \
   FluidFieldsModule.o \
   RadiationFieldsModule.o \
   TwoMoment_DiscretizationModule_Streaming.o \
   TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
   TwoMoment_PositivityLimiterModule.o \
   PhysicalConstantsModule.o \
   GeometryFieldsModuleE.o \
   GeometryFieldsModule.o \
   TimeSteppingModule_Castro.f90
