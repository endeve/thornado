#!/usr/bin/env python3

import os

### Begin user input

CXX        = "mpic++"
FC         = "mpif90"
FLINKER    = "mpif90"
HDF5_INC   = "HDF5_DIR/include"
HDF5_LIB   = "HDF5_DIR/lib"
LAPACK_INC = "LAPACK_DIR/include"
LAPACK_LIB = "LAPACK_DIR/lib"

### End user input

env = os.environ

try:
    THORNADO_DIR = env['THORNADO_DIR']
except:
    print( '>>> environment variable THORNADO_DIR not set. Exiting...' )
    exit()

try:
    THORNADO_MACHINE = env['THORNADO_MACHINE']
except:
    print( '>>> environment variable THORNADO_MACHINE not set. Exiting...' )
    exit()

st = \
f"""## Makefile definitions for {THORNADO_MACHINE}

C       = {CXX} --c++11
FORTRAN = {FC} -cpp
FLINKER = {FLINKER}

DEBUG    = -g -O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized -Wunused \\
                  -ffpe-trap=invalid,zero,overflow,underflow  \\
                  -ffpe-summary=invalid,zero,overflow,underflow \\
                  -fallow-argument-mismatch -finit-real=snan -ftrapv
OPTIMIZE = -g -O2 -fallow-argument-mismatch

SUFFIX_f90 =

MDEFS =
PP    = -D

ifeq ($(USE_OMP),TRUE)
  OPENMP = -fopenmp
else ifeq ($(USE_OMP_OL),TRUE)
  OPENMP = -fopenmp
endif

ifeq ($(USE_OACC),TRUE)
  OPENACC = -fopenacc
endif

INCLUDE_CUDA     =
INCLUDE_HDF5     = -I{HDF5_INC}
INCLUDE_LAPACK   = -I{LAPACK_INC}
INCLUDE_MAGMA    =
INCLUDE_PETSC    =
INCLUDE_ROCM     =

LIBRARIES_CUDA     =
LIBRARIES_HDF5     = -L{HDF5_LIB} -lhdf5_fortran -lhdf5
LIBRARIES_LAPACK   = -L{LAPACK_LIB} -llapack -lblas
LIBRARIES_MAGMA    =
LIBRARIES_PETSC    =
LIBRARIES_ROCM     =
"""

fileName \
  = THORNADO_DIR + "/Build/Machines/Makefile_{:}".format( THORNADO_MACHINE )

ow = True
if ( os.path.isfile( fileName ) ):
    yn = input( '{:} exists. Overwrite? [y/N]: '.format( fileName ) )
    ow = False
    if ( yn == 'y' ): ow = True

if ( ow ):
    with open( fileName, 'w' ) as f:
        f.write(st)
    print( '\n  Created {:}'.format( fileName ) )
else:
    print( '\n  Not overwriting.' )
