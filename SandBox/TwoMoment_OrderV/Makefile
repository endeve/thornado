MOMENT_CLOSURE  = MINERBO
MICROPHYSICS   ?=
GRAVITY_SOLVER ?=
TWOMOMENT_ORDER = ORDER_V

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
ifeq ($(MICROPHYSICS),WEAKLIB)
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies
endif

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: TwoMomentClosure PrimitiveConserved \
     ApplicationDriver ApplicationDriver_Neutrinos

TwoMomentClosure: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMomentClosure.o
	$(FLINKER) $(FLAGS) -o TwoMomentClosure_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMomentClosure.o \
	$(LIBRARIES)

PrimitiveConserved: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	PrimitiveConserved.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	PrimitiveConserved.o \
	$(LIBRARIES)

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_OpacityModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_OpacityModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ApplicationDriver_Neutrinos: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule_Neutrinos.o \
	ApplicationDriver_Neutrinos.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_Neutrinos_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule_Neutrinos.o \
	ApplicationDriver_Neutrinos.o \
	$(LIBRARIES)

ApplicationDriver_CCSN: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	ProgenitorModule.o \
	InitializationModule_CCSN.o \
	ApplicationDriver_CCSN.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_CCSN_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	ProgenitorModule.o \
	InitializationModule_CCSN.o \
	ApplicationDriver_CCSN.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  TwoMomentClosure_$(MACHINE)
	rm -rf TwoMomentClosure_$(MACHINE).dSYM
	rm -f  PrimitiveConserved_$(MACHINE)
	rm -rf PrimitiveConserved_$(MACHINE).dSYM
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM
	rm -f  ApplicationDriver_Neutrinos_$(MACHINE)
	rm -rf ApplicationDriver_Neutrinos_$(MACHINE).dSYM
	rm -f  ApplicationDriver_CCSN_$(MACHINE)
	rm -rf ApplicationDriver_CCSN_$(MACHINE).dSYM

TwoMomentClosure.o: \
  KindModule.o \
  UtilitiesModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMomentClosure.f90

PrimitiveConserved.o: \
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
  RadiationFieldsModule.o \
  InputOutputModuleHDF.o \
  UtilitiesModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  PrimitiveConserved.f90

ApplicationDriver.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  InputOutputModuleHDF.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.o \
  TwoMoment_OpacityModule_OrderV.o \
  TwoMoment_TimeSteppingModule_OrderV.o \
  InitializationModule.o \
  TwoMoment_TimersModule_OrderV.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  GeometryComputationModule.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  GeometryComputationModuleE.o \
  ReferenceElementModuleZ.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  EquationOfStateModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  ApplicationDriver.F90

ApplicationDriver_Neutrinos.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  Euler_SlopeLimiterModule_NonRelativistic_TABLE.o \
  Euler_PositivityLimiterModule_NonRelativistic_TABLE.o \
  RadiationFieldsModule.o \
  InputOutputModuleHDF.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
  TwoMoment_TimeSteppingModule_OrderV.o \
  InitializationModule_Neutrinos.o \
  TwoMoment_TallyModule_OrderV.o \
  EquationOfStateModule_TABLE.o \
  OpacityModule_TABLE.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModuleZ.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  ProgramInitializationModule.o \
  TwoMoment_TimersModule_OrderV.o \
  GeometryComputationModule.o \
  GeometryComputationModuleE.o \
  TwoMoment_ClosureModule.o \
  ApplicationDriver_Neutrinos.F90

ApplicationDriver_CCSN.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  Euler_SlopeLimiterModule_NonRelativistic_TABLE.o \
  Euler_PositivityLimiterModule_NonRelativistic_TABLE.o \
  RadiationFieldsModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
  GravitySolutionModule_Newtonian_Poseidon.o \
  TwoMoment_TimeSteppingModule_OrderV.o \
  InputOutputModuleHDF.o \
  ApplicationDriver_CCSN.f90

InitializationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleZ.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  RadiationFieldsModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_OpacityModule_OrderV.o \
  InitializationModule.F90

InitializationModule_Neutrinos.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  RadiationFieldsModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  EquationOfStateModule_TABLE.o \
  InitializationModule_Neutrinos.F90

TwoMoment_TimersModule_OrderV.o: \
  KindModule.o \
  TwoMoment_TimersModule_OrderV.F90

TwoMoment_UtilitiesModule_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  RadiationFieldsModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule_OrderV.F90

TwoMoment_TroubledCellIndicatorModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TwoMoment_TimersModule_OrderV.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  PolynomialBasisModuleX_Lagrange.o \
  RadiationFieldsModule.o \
  TwoMoment_BoundaryConditionsModule.o \
  TwoMoment_TroubledCellIndicatorModule.F90

TwoMoment_SlopeLimiterModule_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TwoMoment_TimersModule_OrderV.o \
  LinearAlgebraModule.o \
  UtilitiesModule.o \
  ReferenceElementModuleX.o \
  PolynomialBasisMappingModule.o \
  MeshModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_BoundaryConditionsModule.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  TwoMoment_SlopeLimiterModule_OrderV.F90

TwoMoment_PositivityLimiterModule_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TwoMoment_TimersModule_OrderV.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  LinearAlgebraModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_PositivityLimiterModule_OrderV.F90

TwoMoment_OpacityModule_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  MeshModule.o \
  RadiationFieldsModule.o \
  TwoMoment_OpacityModule_OrderV.F90

TwoMoment_DiscretizationModule_Streaming_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TwoMoment_TimersModule_OrderV.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_BoundaryConditionsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  RadiationFieldsModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_BoundaryConditionsModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.F90

TwoMoment_DiscretizationModule_Collisions_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  RadiationFieldsModule.o \
  TwoMoment_TimersModule_OrderV.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_OpacityModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.F90

TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  EquationOfStateModule_TABLE.o \
  TwoMoment_NeutrinoMatterSolverModule_OrderV.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.F90

TwoMoment_NeutrinoMatterSolverModule_OrderV.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  TimersModule.o \
  ArrayUtilitiesModule.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleE.o \
  MeshModule.o \
  RadiationFieldsModule.o \
  EquationOfStateModule_TABLE.o \
  NeutrinoOpacitiesComputationModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_NeutrinoMatterSolverModule_OrderV.F90

TwoMoment_TallyModule_OrderV.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModule.o \
  MeshModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  RadiationFieldsModule.o \
  TwoMoment_TallyModule_OrderV.F90

TwoMoment_TimeSteppingModule_OrderV.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_TimersModule_OrderV.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  TwoMoment_SlopeLimiterModule_OrderV.o \
  TwoMoment_PositivityLimiterModule_OrderV.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.o \
  TwoMoment_TimeSteppingModule_OrderV.F90
