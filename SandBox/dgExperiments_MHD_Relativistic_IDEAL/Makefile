THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SANDBOX/DGEXPERIMENTS_MHD_RELATIVISTIC_IDEAL

DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
        MHD_UtilitiesModule_Relativistic.o \
        MHD_UtilitiesModule.o \
        MHD_PerturbationModule.o \
        MHD_BoundaryConditionsModule.o \
        MHD_TallyModule_Relativistic.o \
        MHD_TallyModule.o \
        MHD_CharacteristicDecompositionModule_Relativistic_IDEAL.o \
        MHD_CharacteristicDecompositionModule.o \
        MHD_SlopeLimiterModule_Relativistic_IDEAL.o \
        MHD_SlopeLimiterModule.o \
        MHD_DiscretizationModule_Relativistic.o \
        TimeSteppingModule_SSPRK.o \
        InitializationModule.o \
        ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
        MHD_UtilitiesModule_Relativistic.o \
        MHD_UtilitiesModule.o \
        MHD_PerturbationModule.o \
        MHD_BoundaryConditionsModule.o \
        MHD_TallyModule_Relativistic.o \
        MHD_TallyModule.o \
        MHD_CharacteristicDecompositionModule_Relativistic_IDEAL.o \
        MHD_CharacteristicDecompositionModule.o \
        MHD_SlopeLimiterModule_Relativistic_IDEAL.o \
        MHD_SlopeLimiterModule.o \
        MHD_DiscretizationModule_Relativistic.o \
        TimeSteppingModule_SSPRK.o \
        InitializationModule.o \
        ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf *.dSYM

ApplicationDriver.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  MHD_PerturbationModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  GeometryComputationModule.o \
  InitializationModule.o \
  UnitsModule.o \
  MHD_TallyModule_Relativistic.o \
  TimeSteppingModule_SSPRK.o \
  ApplicationDriver.F90

InitializationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  EquationOfStateModule_IDEAL.o \
  MHD_UtilitiesModule_Relativistic.o \
  UnitsModule.o \
  QuadratureModule.o \
  PolynomialBasisModule_Lagrange.o \
  UtilitiesModule.o \
  InitializationModule.F90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  MagnetofluidFieldsModule.o \
  MHD_DiscretizationModule_Relativistic.o \
  MHD_TallyModule_Relativistic.o \
  TimeSteppingModule_SSPRK.F90

MHD_BoundaryConditionsModule.o: \
  KindModule.o \
  MeshModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  InputOutputUtilitiesModule.o \
  UtilitiesModule.o \
  MagnetofluidFieldsModule.o \
  MHD_BoundaryConditionsModule.F90

MHD_DiscontinuityDetectionModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  LinearAlgebraModule.o \
  PolynomialBasisModuleX_Lagrange.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule.o \
  EquationOfStateModule.o \
  UtilitiesModule.o \
  TimersModule_Euler.o \
  MHD_ErrorModule.o \
  MHD_DiscontinuityDetectionModule.F90

MHD_DiscretizationModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_BoundaryConditionsModule.o \
  MHD_UtilitiesModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  EquationOfStateModule.o \
  MHD_DiscretizationModule_Relativistic.F90

MHD_TallyModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  UnitsModule.o \
  MHD_TallyModule_Relativistic.F90

MHD_TallyModule.o: \
  KindModule.o \
  MagnetofluidFieldsModule.o \
  MHD_TallyModule_Relativistic.o \
  MHD_TallyModule.F90

MHD_SlopeLimiterModule.o: \
  KindModule.o \
  MHD_BoundaryConditionsModule.o \
  MHD_SlopeLimiterModule_Relativistic_IDEAL.o \
  MHD_SlopeLimiterModule.F90

MHD_SlopeLimiterModule_Relativistic_IDEAL.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  UtilitiesModule.o \
  PolynomialBasisModuleX_Lagrange.o \
  PolynomialBasisModuleX_Legendre.o \
  PolynomialBasisMappingModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_BoundaryConditionsModule.o \
  MHD_CharacteristicDecompositionModule_Relativistic_IDEAL.o \
  MHD_DiscontinuityDetectionModule.o \
  MHD_SlopeLimiterModule_Relativistic_IDEAL.F90

MHD_CharacteristicDecompositionModule_Relativistic_IDEAL.o: \
  KindModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  EquationOfStateModule.o \
  MHD_CharacteristicDecompositionModule_Relativistic_IDEAL.f90

MHD_CharacteristicDecompositionModule.o: \
  KindModule.o \
  MagnetofluidFieldsModule.o \
  MHD_CharacteristicDecompositionModule_Relativistic_IDEAL.o \
  MHD_CharacteristicDecompositionModule.F90

MHD_UtilitiesModule.o: \
  KindModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  MHD_UtilitiesModule.F90

MHD_UtilitiesModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_BoundaryConditionsModule.o \
  EquationOfStateModule.o \
  EquationOfStateModule_TABLE.o \
  UnitsModule.o \
  MHD_UtilitiesModule_Relativistic.F90

MHD_ErrorModule.o: \
  KindModule.o \
  UnitsModule.o \
  UtilitiesModule.o \
  GeometryFieldsModule.o \
  MHD_ErrorModule.f90

MHD_PerturbationModule.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  EquationOfStateModule.o \
  MHD_PerturbationModule.f90
