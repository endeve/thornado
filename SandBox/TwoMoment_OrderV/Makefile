MOMENT_CLOSURE = MINERBO
#MICROPHYSICS   = WEAKLIB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

.DEFAULT_GOAL := all

all: TwoMomentClosure PrimitiveConserved TimeStepper \
     ApplicationDriver ApplicationDriverNeutrinos

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
	TwoMoment_TimeSteppingModule_OrderV.o \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

ApplicationDriverNeutrinos: \
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
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	ApplicationDriverNeutrinos.o
	$(FLINKER) $(FLAGS) -o ApplicationDriverNeutrinos_$(MACHINE) \
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
	TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.o \
	TwoMoment_TimeSteppingModule_OrderV.o \
	ApplicationDriverNeutrinos.o \
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
	rm -f  ApplicationDriverNeutrinos_$(MACHINE)
	rm -rf ApplicationDriverNeutrinos_$(MACHINE).dSYM

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

TimeStepper.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ProgramInitializationModule.o \
  TimersModule.o \
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
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule_OrderV.o \
  TwoMoment_TimeSteppingModule_OrderV.o \
  TimeStepper.f90

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
  ApplicationDriver.f90

ApplicationDriverNeutrinos.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  GeometryFieldsModuleE.o \
  ApplicationDriverNeutrinos.f90

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
  ProgramHeaderModule.o \
  GeometryFieldsModuleE.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  Euler_UtilitiesModule_NonRelativistic.o \
  EquationOfStateModule_TABLE.o \
  TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV.F90

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
