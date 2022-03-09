MICROPHYSICS    = WEAKLIB
TWOMOMENT_ORDER = ORDER_1

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies

DEFINES += $(PP)THORNADO_EULER_NOGPU
#DEFINES += $(PP)THORNADO_DEBUG

.DEFAULT_GOAL := all

all: NeutrinoOpacities NeutrinoOpacities_Profile

NeutrinoOpacities: \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities.o
	$(FLINKER) $(FLAGS) -o NeutrinoOpacities_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities.o \
	$(LIBRARIES)

NeutrinoOpacities_Profile: \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities_Profile.o
	$(FLINKER) $(FLAGS) -o NeutrinoOpacities_Profile_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities_Profile.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f NeutrinoOpacities_$(MACHINE)

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
  NeutrinoOpacities.f90
