# AMREX_HOME defines the directory in which we will find all the AMReX code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME    ?= $(AMREX_DIR)
THORNADO_HOME ?= $(THORNADO_DIR)
WEAKLIB_HOME  ?= $(WEAKLIB_DIR)
HDF5_HOME     ?= $(HDF5_DIR)
LAPACK_HOME   ?= $(LAPACK_DIR)

HDF5_INC       = $(HDF5_HOME)/include
HDF5_LIB       = $(HDF5_HOME)/lib

LAPACK_INC     = $(LAPACK_HOME)/include
LAPACK_LIB     = $(LAPACK_HOME)/lib

EOS            ?= IDEAL
REFINE_MESH    ?= TRUE

DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE

#COMP      = gnu

#MEM_PROFILE  = FALSE
#TINY_PROFILE = FALSE

DIM       = 1

EBASE     = ApplicationDriver

INCLUDE_LOCATIONS += $(HDF5_INC)
LIBRARY_LOCATIONS += $(HDF5_LIB)
LIBRARIES += -lhdf5_fortran -lhdf5

ifeq ($(EOS),TABLE)
  include $(THORNADO_HOME)/SandBox/AMReX/Euler_MR/Make.weaklib
endif

INCLUDE_LOCATIONS += $(LAPACK_INC)
LIBRARY_LOCATIONS += $(LAPACK_LIB)
LIBRARIES += -llapack -lblas

# Set pre-processor directives
# This is also where, e.g., USE_HDF5 and USE_PETSC are used
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Include basic AMReX framework
include $(AMREX_HOME)/Src/Base/Make.package

ifeq ($(REFINE_MESH),TRUE)
  # Requires amrex fork from github.com:dunhamsj/amrex
  # Branch: MeshRefinement_DG
  include $(AMREX_HOME)/Src/DG/Make.package
endif

# Include directories and source files needed for amrex_parallel_module
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/AmrCore/Make.package

# Include directories and source files needed from thornado
include $(THORNADO_HOME)/SandBox/AMReX/Euler_MR/Make.thornado

# Include rules for compilation
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
