MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

DEFINES += -DTHORNADO_USE_AMREX

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/AMReX/TwoMoment_Relativistic

.DEFAULT_GOAL := all

all: PrimitiveConserved_Richardson

PrimitiveConserved_Richardson: \
	$(thornado) \
	TwoMoment_UtilitiesModule_Relativistic.o \
        TwoMoment_TimersModule_Relativistic.o    \
	PrimitiveConserved_Richardson.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_Richardson_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_Relativistic.o \
        TwoMoment_TimersModule_Relativistic.o    \
	PrimitiveConserved_Richardson.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f  PrimitiveConserved_Richardson_$(MACHINE)
	rm -rf PrimitiveConserved_Richardson_$(MACHINE).dSYM

PrimitiveConserved_Richardson.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  PrimitiveConserved_Richardson.F90
