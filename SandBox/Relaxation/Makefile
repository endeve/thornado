MICROPHYSICS           = WEAKLIB
MOMENT_CLOSURE         = MAXIMUM_ENTROPY_CB
NEUTRINO_MATTER_SOLVER = NEWTON

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

#DEFINES += $(PP)-DTHORNADO_DEBUG

.DEFAULT_GOAL := all

all: Relaxation

Relaxation: \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	Relaxation.o
	$(FLINKER) $(FLAGS) -o \
	Relaxation_$(NEUTRINO_MATTER_SOLVER)_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InitializationModule.o \
	Relaxation.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Relaxation_$(NEUTRINO_MATTER_SOLVER)_$(MACHINE)

Relaxation.o: \
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
   TwoMoment_DiscretizationModule_Collisions_Neutrinos.o \
   Relaxation.f90

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
   InitializationModule.f90
