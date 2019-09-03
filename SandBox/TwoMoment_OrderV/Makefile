MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

.DEFAULT_GOAL := all

all: PrimitiveConserved

PrimitiveConserved: \
	$(thornado) \
	PrimitiveConserved.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_$(MACHINE) \
	$(thornado) \
	PrimitiveConserved.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f PrimitiveConserved_$(MACHINE)

PrimitiveConserved.o: \
  KindModule.o \
  ProgramInitializationModule.o \
  PrimitiveConserved.f90
