MOMENT_CLOSURE = MINERBO

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_GR

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	Euler_BoundaryConditionsModule_Relativistic.o \
	InitializationModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	Euler_BoundaryConditionsModule_Relativistic.o \
	InitializationModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  ApplicationDriver_$(MACHINE)
	rm -rf ApplicationDriver_$(MACHINE).dSYM

ApplicationDriver.o: \
  KindModule.o \
  ApplicationDriver.f90

InitializationModule.o: \
  KindModule.o \
  InitializationModule.F90
