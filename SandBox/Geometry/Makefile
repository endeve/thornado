MOMENT_CLOSURE = MINERBO
#MICROPHYSICS   = WEAKLIB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/Geometry

.DEFAULT_GOAL := all

all: ApplicationDriverGeometry



ApplicationDriverGeometry: \
	$(weaklib) \
	$(thornado) \
        WeakDerivatives.o \
	ApplicationDriverGeometry.o
	$(FLINKER) $(FLAGS) -o ApplicationDriverGeometry_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
        WeakDerivatives.o \
	ApplicationDriverGeometry.o \
	$(LIBRARIES)




clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriverGeometry_$(MACHINE)
	rm -rf ApplicationDriverGeometry_$(MACHINE).dSYM

ApplicationDriverGeometry.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramHeaderModule.o \
  GeometryFieldsModule.o \
  WeakDerivatives.o \
  ApplicationDriverGeometry.f90

WeakDerivatives.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  LinearAlgebraModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  ReferenceElementModuleX.o \
  Euler_UtilitiesModule_Relativistic.o \
  WeakDerivatives.F90
