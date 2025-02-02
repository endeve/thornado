# XXX_HOME defines the directory in which all the XXX source code is located
# If you set XXX_HOME as an environment variable, this line will be ignored
HDF5_HOME     ?= $(HDF5_DIR)
LAPACK_HOME   ?= $(LAPACK_DIR)
WEAKLIB_HOME  ?= $(WEAKLIB_DIR)
AMREX_HOME    ?= $(AMREX_DIR)
THORNADO_HOME ?= $(THORNADO_DIR)

DEBUG = FALSE

USE_MPI = TRUE
USE_OMP = FALSE

DIM = 1

DEFINES += -DTHORNADO_USE_AMREX
#DEFINES += -DTHORNADO_OMP
DEFINES += -DHYDRO_RELATIVISTIC
DEFINES += -DHYDRO_RIEMANN_SOLVER_HLL
DEFINES += -DMICROPHYSICS_WEAKLIB
DEFINES += -DTHORNADO_NOTRANSPORT

include Make.package

# Include directories and source files needed from thornado, etc.
include $(THORNADO_HOME)/SandBox/AMReX/Build/Make.weaklib
include $(THORNADO_HOME)/SandBox/AMReX/Build/Make.Euler

include $(THORNADO_HOME)/SandBox/AMReX/Make.thornado
