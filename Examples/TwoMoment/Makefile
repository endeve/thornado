MICROPHYSICS           = WEAKLIB
MOMENT_CLOSURE         = MINERBO
NEUTRINO_MATTER_SOLVER = FIXED_POINT_NESTED_AA
TWOMOMENT_ORDER        = ORDER_1

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/Examples/TwoMoment
VPATH += $(THORNADO_DIR)/SandBox/Interface_FLASH

#DEFINES += $(PP)-DTHORNADO_DEBUG
DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_Flash.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld *.lst *.acc.s

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -rf INLINE_$(MACHINE)

ApplicationDriver.o: \
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
  NeutrinoOpacitiesComputationModule.o \
  TimeSteppingModule_Flash.o \
  InitializationModule.o \
  InputOutputModuleHDF.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o \
  ApplicationDriver.F90

InitializationModule.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModule.o \
  MeshModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  TwoMoment_UtilitiesModule.o \
  EquationOfStateModule_TABLE.o \
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
