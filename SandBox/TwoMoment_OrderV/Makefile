MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

VPATH += $(THORNADO_DIR)/SandBox/TwoMoment_OrderV

.DEFAULT_GOAL := all

all: PrimitiveConserved

PrimitiveConserved: \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
	PrimitiveConserved.o
	$(FLINKER) $(FLAGS) -o PrimitiveConserved_$(MACHINE) \
	$(thornado) \
	TwoMoment_UtilitiesModule_OrderV.o \
	TwoMoment_DiscretizationModule_Streaming_OrderV.o \
	TwoMoment_DiscretizationModule_Collisions_OrderV.o \
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

TwoMoment_UtilitiesModule.o: \
  KindModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_UtilitiesModule.F90

TwoMoment_DiscretizationModule_Streaming_OrderV.o: \
  KindModule.o \
  TwoMoment_DiscretizationModule_Streaming_OrderV.F90

TwoMoment_DiscretizationModule_Collisions_OrderV.o: \
  KindModule.o \
  TwoMoment_DiscretizationModule_Collisions_OrderV.F90
