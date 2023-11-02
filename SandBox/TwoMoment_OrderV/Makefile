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

POSEIDON_DIR ?= $(HOME)/Poseidon
ifeq ($(GRAVITY_SOLVER),POSEIDON_NEWTON)
include $(POSEIDON_DIR)/Build/Make.Poseidon_Native
endif

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: TwoMomentClosure PrimitiveConserved \
     ApplicationDriver ApplicationDriver_Neutrinos

TwoMomentClosure: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	TwoMomentClosure.o
	$(FLINKER) $(FLAGS) -o TwoMomentClosure_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	TwoMomentClosure.o \
	$(LIBRARIES)

PrimitiveConserved: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	PrimitiveConserved.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	PrimitiveConserved.o \
	$(LIBRARIES)

ApplicationDriver: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ApplicationDriver_Neutrinos: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule_Neutrinos.o \
	ApplicationDriver_Neutrinos.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_Neutrinos_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule_Neutrinos.o \
	ApplicationDriver_Neutrinos.o \
	$(LIBRARIES)

ApplicationDriver_CCSN: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ProgenitorModule.o \
	InitializationModule_CCSN.o \
	ApplicationDriver_CCSN.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_CCSN_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	ProgenitorModule.o \
	InitializationModule_CCSN.o \
	ApplicationDriver_CCSN.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld *.lst *.i *.s

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
  TwoMoment_UtilitiesModule.o \
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
  TwoMoment_UtilitiesModule.o \
  PrimitiveConserved.f90

ApplicationDriver.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  InputOutputModuleHDF.o \
  TwoMoment_UtilitiesModule.o \
  TwoMoment_SlopeLimiterModule.o \
  TwoMoment_PositivityLimiterModule.o \
  TwoMoment_DiscretizationModule_Collisions.o \
  TwoMoment_OpacityModule.o \
  TwoMoment_TimeSteppingModule.o \
  InitializationModule.o \
  TwoMoment_TimersModule.o \
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
  TwoMoment_UtilitiesModule.o \
  TwoMoment_SlopeLimiterModule.o \
  TwoMoment_PositivityLimiterModule.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
  TwoMoment_TimeSteppingModule.o \
  InitializationModule_Neutrinos.o \
  TwoMoment_TallyModule.o \
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
  TwoMoment_TimersModule.o \
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
  TwoMoment_UtilitiesModule.o \
  TwoMoment_SlopeLimiterModule.o \
  TwoMoment_PositivityLimiterModule.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
  GravitySolutionModule_Newtonian_Poseidon.o \
  TwoMoment_TimeSteppingModule.o \
  InputOutputModuleHDF.o \
  ApplicationDriver_CCSN.F90

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
  TwoMoment_UtilitiesModule.o \
  TwoMoment_OpacityModule.o \
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
  TwoMoment_UtilitiesModule.o \
  EquationOfStateModule_TABLE.o \
  InitializationModule_Neutrinos.F90
