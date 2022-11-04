THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SANDBOX/DGEXPERIMENTS_MHD_RELATIVISTIC_IDEAL

DEFINES += $(PP)THORNADO_EULER_NOGPU

.DEFAULT_GOAL := all

all: ApplicationDriver ComputePrimitiveTest

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
        MHD_UtilitiesModule_Relativistic.o \
        MHD_UtilitiesModule.o \
        MHD_BoundaryConditionsModule.o \
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
        MHD_BoundaryConditionsModule.o \
        MHD_SlopeLimiterModule_Relativistic_IDEAL.o \
        MHD_SlopeLimiterModule.o \
        MHD_DiscretizationModule_Relativistic.o \
        TimeSteppingModule_SSPRK.o \
        InitializationModule.o \
        ApplicationDriver.o \
	$(LIBRARIES)

ComputePrimitiveTest: \
	$(thornado) \
        MHD_UtilitiesModule_Relativistic.o \
      	ComputePrimitiveTest.o
	$(FLINKER) $(FLAGS) -o ComputePrimitiveTest_$(MACHINE) \
        $(thornado) \
        MHD_UtilitiesModule_Relativistic.o \
        ComputePrimitiveTest.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -f  ComputePrimitiveTest_$(MACHINE)
	rm -rf *.dSYM

ApplicationDriver.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  GeometryComputationModule.o \
  InitializationModule.o \
  UnitsModule.o \
  ApplicationDriver.F90

ComputePrimitiveTest.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  GeometryComputationModule.o \
  MagnetofluidFieldsModule.o \
  EquationOfStateModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  ComputePrimitiveTest.F90

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
  InitializationModule.F90

TimeSteppingModule_SSPRK.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  MagnetofluidFieldsModule.o \
  MHD_DiscretizationModule_Relativistic.o \
  TimeSteppingModule_SSPRK.F90

MHD_BoundaryConditionsModule.o: \
  KindModule.o \
  MeshModule.o \
  ProgramHeaderModule.o \
  MHD_UtilitiesModule.o \
  MagnetofluidFieldsModule.o \
  MHD_BoundaryConditionsModule.F90

MHD_dgDiscretizationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  MHD_BoundaryConditionsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  EquationOfStateModule.o \
  MHD_DiscretizationModule_Relativistic.F90

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
  MHD_SlopeLimiterModule_Relativistic_IDEAL.F90

MHD_UtilitiesModule.o: \
  KindModule.o \
  MagnetofluidFieldsModule.o \
  MHD_UtilitiesModule_Relativistic.o \
  MHD_UtilitiesModule.F90

MHD_UtilitiesModule_Relativistic.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  MagnetofluidFieldsModule.o \
  EquationOfStateModule.o \
  EquationOfStateModule_TABLE.o \
  UnitsModule.o \
  MHD_UtilitiesModule_Relativistic.F90
