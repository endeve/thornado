MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = LEVERMORE
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
#MOMENT_CLOSURE = KERSHAW_BL

THORNADO_DIR ?= $(HOME)/thornado
include $(THORNADO_DIR)/Build/Makefile_Build

all: ApplicationDriver

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

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f ApplicationDriver_$(MACHINE)
