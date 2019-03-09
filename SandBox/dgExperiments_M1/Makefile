#MICROPHYSICS = WEAKLIB
MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
#MOMENT_CLOSURE = KERSHAW_BL

include ../../../Build/Makefile_Build

all: ApplicationDriver DeleptonizationWave

ApplicationDriver: \
	$(thornado) \
	ProgenitorModule.o \
	TwoMoment_TallyModule.o \
	TimeSteppingModule_IMEX_RK.o \
	InitializationModule.o \
	dgDiscretizationModule_Collisions.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	ProgenitorModule.o \
	TwoMoment_TallyModule.o \
	TimeSteppingModule_IMEX_RK.o \
	InitializationModule.o \
	dgDiscretizationModule_Collisions.o \
	ApplicationDriver.o \
	$(LIBRARIES)

DeleptonizationWave: \
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
