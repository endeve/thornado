MICROPHYSICS = WEAKLIB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies


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

all: NeutrinoOpacities

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
