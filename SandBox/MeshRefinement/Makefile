TWOMOMENT_ORDER = ORDER_V

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

DEFINES += $(PP)THORNADO_EULER_NOGPU
#DEFINES += $(PP)THORNADO_DEBUG

.DEFAULT_GOAL := all

all: MeshRefinementTest_TwoMoment

MeshRefinementTest_TwoMoment: \
	$(thornado) \
	MeshRefinementTest_TwoMoment.o
	$(FLINKER) $(FLAGS) -o MeshRefinementTest_TwoMoment_$(MACHINE) \
	$(thornado) \
	MeshRefinementTest_TwoMoment.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld *.s

clobber: clean
	rm -f MeshRefinementTest_TwoMoment_$(MACHINE)

MeshRefinementTest_TwoMoment.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ProgramInitializationModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  MeshModule.o \
  UtilitiesModule.o \
  TwoMoment_MeshRefinementModule.o \
  DeviceModule.o \
  TimersModule.o \
  MeshRefinementTest_TwoMoment.F90
