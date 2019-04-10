MICROPHYSICS = WEAKLIB
#MOMENT_CLOSURE = MINERBO
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
#MOMENT_CLOSURE = KERSHAW_BL

THORNADO_DIR ?= $(HOME)/thornado
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: ApplicationDriver DeleptonizationWave

ApplicationDriver: \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	TwoMoment_TallyModule.o \
	TimeSteppingModule_IMEX_RK.o \
	InitializationModule.o \
	dgDiscretizationModule_Collisions.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	TwoMoment_TallyModule.o \
	TimeSteppingModule_IMEX_RK.o \
	InitializationModule.o \
	dgDiscretizationModule_Collisions.o \
	ApplicationDriver.o \
	$(LIBRARIES)

DeleptonizationWave: \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	GeometryComputationModule_Beta.o \
	GeometryComputationModuleE_Beta.o \
	MomentEquationsUtilitiesModule_Beta.o \
	ClosureModule_M1.o \
	PositivityLimiterModule.o \
	TimeSteppingModule_IMEX_RK.o \
	InitializationModule.o \
	dgDiscretizationModule_Collisions_Neutrinos.o \
	DeleptonizationWave.o
	$(FLINKER) $(FLAGS) -o DeleptonizationWave_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	ProgenitorModule.o \
	GeometryComputationModule_Beta.o \
	GeometryComputationModuleE_Beta.o \
	MomentEquationsUtilitiesModule_Beta.o \
	ClosureModule_M1.o \
	PositivityLimiterModule.o \
	TimeSteppingModule_IMEX_RK.o \
	InitializationModule.o \
	dgDiscretizationModule_Collisions_Neutrinos.o \
	DeleptonizationWave.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
	rm -f DeleptonizationWave_$(MACHINE)
