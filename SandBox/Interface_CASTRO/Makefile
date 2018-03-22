MACHINE = $(THORNADO_MACHINE)

OPT_LEVEL = OPTIMIZE
FLAGS     = $(FLAGS_$(OPT_LEVEL))

#MOMENT_CLOSURE = MINERBO
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
#MOMENT_CLOSURE = KERSHAW_BL

include ../../Build/Makefile_Build

VPATH = \
$(THORNADO_DIR)/Modules/Numbers: \
$(THORNADO_DIR)/Modules/ProgramHeader: \
$(THORNADO_DIR)/Modules/Library: \
$(THORNADO_DIR)/Modules/Mesh: \
$(THORNADO_DIR)/Modules/Geometry: \
$(THORNADO_DIR)/Modules/Fields: \
.

all: Driver

Driver: \
	KindModule.o \
	ProgramHeaderModule.o \
	QuadratureModule.o \
	ReferenceElementModuleX.o \
	ReferenceElementModuleE.o \
	ReferenceElementModule.o \
	PolynomialBasisModule_Lagrange.o \
	ReferenceElementModuleX_Lagrange.o \
	ReferenceElementModuleE_Lagrange.o \
	ReferenceElementModule_Lagrange.o \
	MeshModule.o \
	GeometryFieldsModule.o \
	FluidFieldsModule.o \
	GeometryFieldsModuleE.o \
	RadiationFieldsModule.o \
	ThornadoInitializationModule.o \
	Driver.o
	$(FLINKER) $(FLAGS) -o Driver_$(MACHINE) \
	KindModule.o \
	ProgramHeaderModule.o \
	QuadratureModule.o \
	ReferenceElementModuleX.o \
	ReferenceElementModuleE.o \
	ReferenceElementModule.o \
	PolynomialBasisModule_Lagrange.o \
	ReferenceElementModuleX_Lagrange.o \
	ReferenceElementModuleE_Lagrange.o \
	ReferenceElementModule_Lagrange.o \
	MeshModule.o \
	GeometryFieldsModule.o \
	FluidFieldsModule.o \
	GeometryFieldsModuleE.o \
	RadiationFieldsModule.o \
	ThornadoInitializationModule.o \
	Driver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Driver_$(MACHINE)
