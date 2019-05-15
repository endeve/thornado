MICROPHYSICS = WEAKLIB
#MOMENT_CLOSURE = MINERBO
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB
#MOMENT_CLOSURE = MAXIMUM_ENTROPY_BL
MOMENT_CLOSURE = KERSHAW_BL

THORNADO_DIR ?= $(HOME)/thornado
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: Driver

Driver: \
	$(weaklib) \
	$(thornado) \
	Driver.o
	$(FLINKER) $(FLAGS) -o Driver_$(MACHINE)_$(MOMENT_CLOSURE)  \
	$(weaklib) \
	$(thornado) \
	Driver.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Driver_$(MACHINE)_$(MOMENT_CLOSURE)
