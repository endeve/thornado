# AMREX_HOME defines the directory in which we will find all the AMReX code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME    ?= $(AMREX_DIR)
THORNADO_HOME ?= $(THORNADO_DIR)
WEAKLIB_HOME  ?= $(WEAKLIB_DIR)
HDF5_HOME     ?= $(HDF5_DIR)

HDF5_INC       = $(HDF5_HOME)/serial/include
HDF5_LIB       = $(HDF5_HOME)/serial/lib

DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_MPI3  = FALSE

USE_F_INTERFACES = TRUE

MEM_PROFILE  = FALSE
TINY_PROFILE = FALSE

COMP      = gnu

DIM       = 1

EBASE     = main

WEAK ?=TABLE

ifeq ($(WEAK),TABLE)
#  USE_HDF5 = TRUE
  INCLUDE_LOCATIONS += $(HDF5_INC)
  LIBRARY_LOCATIONS += $(HDF5_LIB)
  LIBRARIES += -lhdf5_fortran -lhdf5
  include $(THORNADO_HOME)/SandBox/AMReX/TwoMoment/Make.weaklib
endif


INCLUDE_LOCATIONS += /usr/include
LIBRARY_LOCATIONS += /usr/lib/x86_64-linux-gnu/lapack
LIBRARIES += -llapack -lblas

include $(THORNADO_HOME)/SandBox/AMReX/TwoMoment/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package

### Checkpoint ###
#include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
#include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
