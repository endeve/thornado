MICROPHYSICS    = WEAKLIB
TWOMOMENT_ORDER = ORDER_V

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

all: NeutrinoOpacities #NeutrinoOpacities_Profile

NeutrinoOpacities: \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities.o
	$(FLINKER) $(FLAGS) -o NeutrinoOpacities_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities.o \
	$(LIBRARIES)

NeutrinoOpacities_Muons: \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities_Muons.o
	$(FLINKER) $(FLAGS) -o NeutrinoOpacities_Muons_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	NeutrinoOpacities_Muons.o \
	$(LIBRARIES)

#NeutrinoOpacities_Profile: \
#	$(weaklib) \
#	$(thornado) \
#	NeutrinoOpacities_Profile.o
#	$(FLINKER) $(FLAGS) -o NeutrinoOpacities_Profile_$(MACHINE) \
#	$(weaklib) \
#	$(thornado) \
#	NeutrinoOpacities_Profile.o \
#	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld *.s

clobber: clean
	rm -f NeutrinoOpacities_$(MACHINE)

NeutrinoOpacities.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramInitializationModule.o \
  MeshModule.o \
  UtilitiesModule.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  EquationOfStateModule_TABLE.o \
  OpacityModule_TABLE.o \
  RadiationFieldsModule.o \
  NeutrinoOpacitiesComputationModule.o \
  DeviceModule.o \
  TimersModule.o \
  NeutrinoOpacities.F90

NeutrinoOpacities_Muons.o: \
  KindModule.o \
  UnitsModule.o \
  ProgramInitializationModule.o \
  MeshModule.o \
  UtilitiesModule.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  EquationOfStateModule_TABLE.o \
  OpacityModule_TABLE.o \
  RadiationFieldsModule.o \
  NeutrinoOpacitiesComputationModule.o \
  TimersModule.o \
  NeutrinoOpacities_Muons.F90

ifeq ($(PE_ENV),CRAY)
TwoMoment_PositivityLimiterModule.o : TwoMoment_PositivityLimiterModule.F90
	$(FORTRAN) -c $(DEFINES) $(FLAGS) $(filter-out -qrealsize=8,$(F90FLAGS)) -h ipa1 $< -o $@
endif
