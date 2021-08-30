MOMENT_CLOSURE  = MINERBO
MICROPHYSICS   ?=
GRAVITY_SOLVER ?=
TWOMOMENT_ORDER = ORDER_V

THORNADO_DIR ?= ../../..
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
ifeq ($(MICROPHYSICS),WEAKLIB)
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies
endif

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

Tasmanian_link = -L/home/plc/TASMANIAN_cuda/lib
#Tasmanian_libs = -ltasmaniandream -ltasmaniansparsegrid  /usr/lib/libblas.so /usr/lib/gcc/x86_64-linux-gnu/5/libgomp.so /usr/lib/x86_64-linux-gnu/libpthread.so
Tasmanian_libs = -ltasmaniandream -ltasmaniansparsegrid
Tarmanian_fortran_libs = -ltasmanianfortran90 -ltasmanianfortran $(Tasmanian_libs) -lstdc++

INCLUDE_TASMANIAN_atlasmath = -I/home/plc/TASMANIAN_cuda/include
LIBRARIES_TASMANIAN_atlasmath = $(Tasmanian_link) $(Tarmanian_fortran_libs)

INCLUDE_TASMANIAN = $(INCLUDE_TASMANIAN_$(MACHINE))
LIBRARIES_TASMANIAN = $(LIBRARIES_TASMANIAN_$(MACHINE))

INCLUDE += $(INCLUDE_TASMANIAN)
LIBRARIES += $(LIBRARIES_TASMANIAN)

DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: ApplicationDriver_Neutrinos ApplicationDriver_Neutrinos_SG


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

ApplicationDriver_Neutrinos_SG: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV_SG.o \
	TwoMoment_NeutrinoMatterSolverModule_OrderV_SG.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule_Neutrinos.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	ApplicationDriver_Neutrinos_SG.o 
	$(FLINKER) $(FLAGS) -o ApplicationDriver_Neutrinos_SG_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_TimersModule_OrderV.o \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_TroubledCellIndicatorModule.o \
	TwoMoment_SlopeLimiterModule_OrderV.o \
	TwoMoment_PositivityLimiterModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV_SG.o \
	TwoMoment_NeutrinoMatterSolverModule_OrderV_SG.o \
	TwoMoment_TallyModule_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule_Neutrinos.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	ApplicationDriver_Neutrinos_SG.o \
	$(LIBRARIES)

NeutrinoOpacities: \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacitiesSparseComputationModule.o \
	NeutrinoOpacities.o
	$(FLINKER) $(FLAGS) -o NeutrinoOpacities_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacitiesSparseComputationModule.o \
	NeutrinoOpacities.o \
	$(LIBRARIES)

TasmanianTest: \
	$(weaklib) \
	$(thornado) \
	example_sparse_grids.o
	$(FLINKER) $(FLAGS) -o TasmanianTest_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	example_sparse_grids.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f NeutrinoOpacities_$(MACHINE)
	rm -f TasmanianTest_$(MACHINE)
	rm -f Relaxation_$(NEUTRINO_MATTER_SOLVER)_$(MACHINE)
	rm -f  ApplicationDriver_Neutrinos_$(MACHINE)
	rm -rf ApplicationDriver_Neutrinos_$(MACHINE).dSYM
	rm -f  ApplicationDriver_Neutrinos_SG_$(MACHINE)
	rm -rf ApplicationDriver_Neutrinos_SG_$(MACHINE).dSYM

NeutrinoOpacities.o: \
	KindModule.o \
	UnitsModule.o \
	ProgramInitializationModule.o \
	MeshModule.o \
	UtilitiesModule.o \
	EquationOfStateModule_TABLE.o \
	OpacityModule_TABLE.o \
	RadiationFieldsModule.o \
	NeutrinoOpacitiesComputationModule.o \
	DeviceModule.o \
	NeutrinoOpacitiesSparseComputationModule.o \
  NeutrinoOpacities.f90

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

ApplicationDriver_Neutrinos_SG.o: \
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
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV_SG.o \
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
	NeutrinoOpacitiesSparseComputationModule.o \
  ApplicationDriver_Neutrinos_SG.F90

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

TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV_SG.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  EquationOfStateModule_TABLE.o \
  TwoMoment_NeutrinoMatterSolverModule_OrderV_SG.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV_SG.F90

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

TwoMoment_NeutrinoMatterSolverModule_OrderV_SG.o: \
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
	NeutrinoOpacitiesSparseComputationModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_NeutrinoMatterSolverModule_OrderV_SG.F90

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


NeutrinoOpacitiesSparseComputationModule.o: \
	KindModule.o \
	UnitsModule.o \
	ProgramHeaderModule.o \
	DeviceModule.o \
	LinearAlgebraModule.o \
	ReferenceElementModuleX.o \
	ReferenceElementModule.o \
	MeshModule.o \
	EquationOfStateModule_TABLE.o \
	OpacityModule_TABLE.o \
	RadiationFieldsModule.o \
	NeutrinoOpacitiesModule.o \
	NeutrinoOpacitiesSparseComputationModule.F90


example_sparse_grids.o: \
  example_sparse_grids.f90
