MICROPHYSICS           = WEAKLIB
MOMENT_CLOSURE         = MAXIMUM_ENTROPY_CB
NEUTRINO_MATTER_SOLVER = FIXED_POINT_COUPLED
TWOMOMENT_ORDER        = ORDER_1

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/DeleptonizationWave
VPATH += $(THORNADO_DIR)/SandBox/Interface_FLASH

#DEFINES += $(PP)-DTHORNADO_DEBUG
DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: DeleptonizationWave DeleptonizationWave1D

DeleptonizationWave: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	DeleptonizationWave.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	DeleptonizationWave.o \
	$(LIBRARIES)

DeleptonizationWave1D: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	DeleptonizationWave1D.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave1D_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	DeleptonizationWave1D.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld
	rm -f *.cub *.ptx

clobber: clean
	rm -f DeleptonizationWave_$(MACHINE)
	rm -f DeleptonizationWave1D_$(MACHINE)

DeleptonizationWave.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o \
  ProgramInitializationModule.o \
  TimersModule.o \
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
  TimeSteppingModule_Flash.o \
  InitializationModule.o \
  InputOutputModuleHDF.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o \
  DeleptonizationWave.F90

DeleptonizationWave1D.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o \
  ProgramInitializationModule.o \
  TimersModule.o \
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
  TimeSteppingModule_Flash.o \
  InitializationModule.o \
  InputOutputModuleHDF.o \
  Euler_UtilitiesModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o \
  DeleptonizationWave1D.F90

InitializationModule.o: \
  KindModule.o \
  UnitsModule.o \
  UtilitiesModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  NeutrinoOpacitiesComputationModule.o \
  TwoMoment_UtilitiesModule.o \
  InitializationModule.F90

TimeSteppingModule_Flash.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TimersModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
  TwoMoment_PositivityLimiterModule.o \
  PhysicalConstantsModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  TimeSteppingModule_Flash.F90
