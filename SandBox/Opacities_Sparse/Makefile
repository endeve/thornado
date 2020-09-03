MICROPHYSICS           = WEAKLIB
MOMENT_CLOSURE         = MAXIMUM_ENTROPY_CB
NEUTRINO_MATTER_SOLVER = FIXED_POINT_NESTED_AA

THORNADO_DIR ?= ../../..
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/Interface_CASTRO

Tasmanian_link = -L/home/plc/.tasmanian_dev/lib/
Tasmanian_libs = -ltasmaniandream -ltasmaniansparsegrid  /usr/lib/libblas.so /usr/lib/gcc/x86_64-linux-gnu/5/libgomp.so /usr/lib/x86_64-linux-gnu/libpthread.so
Tarmanian_fortran_libs = -ltasmanianfortran90 $(Tasmanian_libs) -lstdc++

INCLUDE_TASMANIAN_paullaiu_gnu = -I/home/plc/.tasmanian_dev/include/
LIBRARIES_TASMANIAN_paullaiu_gnu = $(Tasmanian_link) $(Tarmanian_fortran_libs)

INCLUDE_TASMANIAN = $(INCLUDE_TASMANIAN_$(MACHINE))
LIBRARIES_TASMANIAN = $(LIBRARIES_TASMANIAN_$(MACHINE))

INCLUDE += $(INCLUDE_TASMANIAN)
LIBRARIES += $(LIBRARIES_TASMANIAN)

#DEFINES += $(PP)-DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: DeleptonizationWave1D

DeleptonizationWave1D: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave1D_SG.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave1D_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	InitializationModule.o \
	TimeSteppingModule_Castro.o \
	DeleptonizationWave1D_SG.o \
	$(LIBRARIES)

Relaxation: \
	$(weaklib) \
	$(thornado) \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	InitializationModule.o \
	Relaxation_SG.o
	$(FLINKER) $(FLAGS) -o \
	Relaxation_$(NEUTRINO_MATTER_SOLVER)_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	InitializationModule.o \
	Relaxation_SG.o \
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


DeleptonizationWave1D_SG.o: \
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
  TimeSteppingModule_Castro.o \
  InitializationModule.o \
  InputOutputModuleHDF.o \
  Euler_UtilitiesModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o \
  DeleptonizationWave1D_SG.F90


Relaxation_SG.o: \
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
	InitializationModule.o \
	InputOutputModuleHDF.o \
	TwoMoment_ClosureModule.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o \
	Relaxation_SG.f90

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
	NeutrinoOpacitiesSparseComputationModule.o \
	InitializationModule.f90

TimeSteppingModule_Castro.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  TimersModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o \
  TwoMoment_PositivityLimiterModule.o \
  PhysicalConstantsModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  TimeSteppingModule_Castro.F90

TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.o: \
	KindModule.o \
	UnitsModule.o \
	UtilitiesModule.o \
	ProgramHeaderModule.o \
	TimersModule.o \
	ReferenceElementModule.o \
	ReferenceElementModuleE.o \
	ReferenceElementModuleX.o \
	MeshModule.o \
	GeometryFieldsModule.o \
	GeometryFieldsModuleE.o \
	FluidFieldsModule.o \
	RadiationFieldsModule.o \
	EquationOfStateModule_TABLE.o \
	NeutrinoOpacitiesModule.o \
	NeutrinoOpacitiesComputationModule.o \
	NeutrinoOpacitiesSparseComputationModule.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_SG.F90

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
