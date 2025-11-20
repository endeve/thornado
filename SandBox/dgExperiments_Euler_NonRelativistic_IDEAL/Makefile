THORNADO_DIR ?= ../../..
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

.DEFAULT_GOAL := all

all: ApplicationDriver

ApplicationDriver: \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	CellMergingModule.o \
	ApplicationDriver.o
	$(FLINKER) $(FLAGS) -o ApplicationDriver_$(MACHINE) \
	$(thornado) \
	InitializationModule.o \
	TimeSteppingModule_SSPRK.o \
	CellMergingModule.o \
	ApplicationDriver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f   ApplicationDriver_$(MACHINE)
	rm -rf  ApplicationDriver_$(MACHINE).dSYM

CellMergingModule.o: \
	Euler_PositivityLimiterModule_NonRelativistic_IDEAL.o\
	CellMergingModule.F90
