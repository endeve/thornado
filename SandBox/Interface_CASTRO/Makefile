MACHINE = $(THORNADO_MACHINE)

OPT_LEVEL = DEBUG
FLAGS     = $(FLAGS_$(OPT_LEVEL))

#MOMENT_CLOSURE = MINERBO
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
#MOMENT_CLOSURE = KERSHAW_BL

MICROPHYSICS = WEAKLIB

include ./Makefile_WeakLib
include ../../Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

VPATH += \
$(THORNADO_DIR)/Modules/Numbers: \
$(THORNADO_DIR)/Modules/Units: \
$(THORNADO_DIR)/Modules/ProgramHeader: \
$(THORNADO_DIR)/Modules/Library: \
$(THORNADO_DIR)/Modules/Mesh: \
$(THORNADO_DIR)/Modules/Geometry: \
$(THORNADO_DIR)/Modules/Fields: \
$(THORNADO_DIR)/Modules/EquationOfState: \
$(THORNADO_DIR)/Modules/Opacities: \
$(THORNADO_DIR)/Modules/TwoMoment: \
.

all: Driver

Driver: $(weaklib) \
	KindModule.o \
	PhysicalConstantsModule.o \
	UnitsModule.o \
	ProgramHeaderModule.o \
	TimersModule.o \
	DeviceModule.o \
	LinearAlgebraModule.o \
	QuadratureModule.o \
	ReferenceElementModuleX.o \
	ReferenceElementModuleE.o \
	ReferenceElementModule.o \
	LagrangePolynomialsModule.o \
	PolynomialBasisModuleX_Lagrange.o \
	PolynomialBasisModule_Lagrange.o \
	ReferenceElementModuleX_Lagrange.o \
	ReferenceElementModuleE_Lagrange.o \
	ReferenceElementModule_Lagrange.o \
	SubcellReconstructionModule.o \
	MeshModule.o \
	GeometryFieldsModule.o \
	GeometryComputationModule.o \
	FluidFieldsModule.o \
	GeometryFieldsModuleE.o \
	GeometryComputationModuleE.o \
	RadiationFieldsModule.o \
	EquationOfStateModule_TABLE.o \
	OpacityModule_TABLE.o \
	NeutrinoOpacitiesModule.o \
	NeutrinoOpacitiesComputationModule.o \
	TwoMoment_MeshRefinementModule.o \
	TwoMoment_ClosureModule.o \
	TwoMoment_BoundaryConditionsModule.o \
	TwoMoment_UtilitiesModule.o \
	TwoMoment_PositivityLimiterModule.o \
	TwoMoment_DiscretizationModule_Streaming.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
	ThornadoInitializationModule.o \
	DiscretizationModule_Dummy.o \
	TimeSteppingModule_Castro.o \
	Driver.o
	$(FLINKER) $(FLAGS) -o Driver_$(MACHINE) \
	$(weaklib) \
	KindModule.o \
	PhysicalConstantsModule.o \
	UnitsModule.o \
	ProgramHeaderModule.o \
	TimersModule.o \
	DeviceModule.o \
	LinearAlgebraModule.o \
	QuadratureModule.o \
	ReferenceElementModuleX.o \
	ReferenceElementModuleE.o \
	ReferenceElementModule.o \
	LagrangePolynomialsModule.o \
	PolynomialBasisModuleX_Lagrange.o \
	PolynomialBasisModule_Lagrange.o \
	ReferenceElementModuleX_Lagrange.o \
	ReferenceElementModuleE_Lagrange.o \
	ReferenceElementModule_Lagrange.o \
	SubcellReconstructionModule.o \
	MeshModule.o \
	GeometryFieldsModule.o \
	GeometryComputationModule.o \
	FluidFieldsModule.o \
	GeometryFieldsModuleE.o \
	GeometryComputationModuleE.o \
	RadiationFieldsModule.o \
	EquationOfStateModule_TABLE.o \
	OpacityModule_TABLE.o \
	NeutrinoOpacitiesModule.o \
	NeutrinoOpacitiesComputationModule.o \
	TwoMoment_ClosureModule.o \
	TwoMoment_MeshRefinementModule.o \
	TwoMoment_BoundaryConditionsModule.o \
	TwoMoment_UtilitiesModule.o \
	TwoMoment_PositivityLimiterModule.o \
	TwoMoment_DiscretizationModule_Streaming.o \
	TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
	ThornadoInitializationModule.o \
	DiscretizationModule_Dummy.o \
	TimeSteppingModule_Castro.o \
	Driver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Driver_$(MACHINE)
