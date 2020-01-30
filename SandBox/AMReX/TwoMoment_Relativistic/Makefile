MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/AMReX/TwoMoment_Relativistic

.DEFAULT_GOAL := all

all: PrimitiveConserved

PrimitiveConserved: \
	$(thornado) \
	TwoMoment_UtilitiesModule_Relativistic.o \
	PrimitiveConserved.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_Relativistic.o \
	PrimitiveConserved.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  PrimitiveConserved_$(MACHINE)
	rm -rf PrimitiveConserved_$(MACHINE).dSYM

PrimitiveConserved.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  PrimitiveConserved.f90
