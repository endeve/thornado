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

VPATH += $(THORNADO_DIR)/SandBox/MuonCollapse: \
VPATH += $(WEAKLIB_DIR)/TableCreator/EOSSource/EquationsOfState/BaryonsPlusLeptons/CommonModules: \
..

include $(WEAKLIB_DIR)/TableCreator/EOSSource/EquationsOfState/BaryonsPlusLeptons/CommonModules/Makefile_CommonModules_Dependencies
include $(WEAKLIB_DIR)/TableCreator/EOSSource/EquationsOfState/BaryonsPlusLeptons/CommonModules/Makefile_CommonModules_ObjectFiles


#DEFINES += $(PP)THORNADO_EULER_NOGPU
#DEFINES += -DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: ApplicationDriver ApplicationDriverComponents \
	ApplicationDriver_Neutrinos ApplicationDriverComponents_Neutrinos

ApplicationDriver: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) $(WL_VERSION_FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ApplicationDriverComponents: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule.o \
	ApplicationDriverComponents.o
	$(FLINKER) $(FLAGS) $(WL_VERSION_FLAGS) -o ApplicationDriverComponents_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule.o \
	ApplicationDriverComponents.o \
	$(LIBRARIES)

ApplicationDriver_Neutrinos: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule_Neutrinos.o \
	ApplicationDriver_Neutrinos.o
	$(FLINKER) $(FLAGS) $(WL_VERSION_FLAGS) -o ApplicationDriver_Neutrinos_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule_Neutrinos.o \
	ApplicationDriver_Neutrinos.o \
	$(LIBRARIES)
	
ApplicationDriverComponents_Neutrinos: \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule_Neutrinos.o \
	ApplicationDriverComponents_Neutrinos.o
	$(FLINKER) $(FLAGS) $(WL_VERSION_FLAGS) -o ApplicationDriverComponents_Neutrinos_$(MACHINE) \
	$(weaklib) \
	$(POSEIDON_o) \
	$(thornado) \
	InitializationModule_Neutrinos.o \
	ApplicationDriverComponents_Neutrinos.o \
	$(LIBRARIES)
clean:
	rm -f *.o *.mod *.ld *.lst *.i *.s

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM
	rm -f  ApplicationDriverComponents_$(MACHINE)
	rm -rf ApplicationDriverComponents_$(MACHINE).dSYM
	rm -f  ApplicationDriver_Neutrinos_$(MACHINE)
	rm -rf ApplicationDriver_Neutrinos_$(MACHINE).dSYM
	rm -f  ApplicationDriverComponents_Neutrinos_$(MACHINE)
	rm -rf ApplicationDriverComponents_Neutrinos_$(MACHINE).dSYM
	
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
  EquationOfStateModule_TABLE.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  ApplicationDriver.F90

ApplicationDriverComponents.o: \
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
  EquationOfStateComponentsModule_TABLE.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_TroubledCellIndicatorModule.o \
  ApplicationDriverComponents.F90

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

ApplicationDriverComponents_Neutrinos.o: \
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
  InitializationModuleComponents_Neutrinos.o \
  TwoMoment_TallyModule.o \
  EquationOfStateComponentsModule_TABLE.o \
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
  ApplicationDriverComponents_Neutrinos.F90

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

InitializationModuleComponents_Neutrinos.o: \
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
  EquationOfStateComponentsModule_TABLE.o \
  InitializationModule_Neutrinos.F90
  
ifeq ($(PE_ENV),CRAY)
TwoMoment_PositivityLimiterModule.o : TwoMoment_PositivityLimiterModule.F90
	$(FORTRAN) -c $(DEFINES) $(FLAGS) $(filter-out -qrealsize=8,$(F90FLAGS)) -h ipa1 $< -o $@
Euler_XCFC_UtilitiesModule.o : Euler_XCFC_UtilitiesModule.f90
	$(FORTRAN) -c $(DEFINES) $(FLAGS) $(filter-out -qrealsize=8,$(F90FLAGS)) -h ipa1 $< -o $@
endif
