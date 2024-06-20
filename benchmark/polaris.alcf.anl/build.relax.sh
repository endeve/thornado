#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function load_set_common(){

  module purge

   module load nvhpc/23.9
   #module load nvhpc/23.3
   module restore
   #module load nvhpc/21.9
   #ml  nvhpc/23.1

   mr=${1:-1}
   date0=`date +%m-%d-%Y`
   hname=`hostname`

   export OP_LEVEL=O3
   export LOG_FILE=relax_a100.${hname}_mpi${mr}_${date0}.log
   rm -f $LOG_FILE
   export APP_NAME=ApplicationDriver_Neutrinos
#   export EXASTAR_HOME=/home/mthavappiragasam/lcf/ExaStar
   export EXASTAR_HOME=/grand/projects/Performance/mathit/ExaStar
#   export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.12.0-nvhpc
   export HDF5_ROOT=/home/mthavappiragasam/lcf/ExaStar/hdf5/1.12.0-nvhpc
   #export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.10.7-a100
   export HDF5_LIB=${HDF5_ROOT}/lib
   export HDF5_INC=${HDF5_ROOT}/include
   export THORNADO_DIR=${EXASTAR_HOME}/thornado
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export CUDA_ROOT=${NVIDIA_PATH}/cuda
   #export CUDA_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/cuda
   #export CUDA_ROOT=/soft/compilers/nvhpc/Linux_x86_64/22.1/compilers
   export NETLIB_LAPACK_ROOT=${NVIDIA_PATH}/math_libs/lib64
   #export NETLIB_LAPACK_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib
   export MAGMA_ROOT=/soft/libraries/magma-2.7.0
   export THORNADO_MACHINE=polaris_nvhpc

   }

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   rm ./${APP_NAME}_${THORNADO_MACHINE}

   module list   |& tee -a $LOG_FILE

   make clean -f ${THORNADO_DIR}/Makefile
   ( time make -j8 -f ${THORNADO_DIR}/Makefile $APP_NAME MICROPHYSICS=WEAKLIB USE_OACC=TRUE  USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE ) |& tee -a $LOG_FILE
   #( time make -j8 $APP_NAME MICROPHYSICS=WEAKLIB USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE ) |& tee -a $LOG_FILE

}

load_set_common
buildApp
