#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function load_set_common(){

   #module purge
   
   #module load oneapi/eng-compiler/2022.01.30.005
   #module load oneapi/eng-compiler/2022.06.30.002
   #module load oneapi/eng-compiler/2022.10.15.006
   #module restore
   #module load oneapi/eng-compiler/2023.05.15.007
   #module load oneapi/release/2023.12.15.001 #Seems losing perfomance
   #module load oneapi/eng-compiler/2023.10.15.002
   #module load mpich/51.2/icc-all-pmix-gpu

   export OP_LEVEL=O3
   export LOG_FILE=sineWave_iprof_.${OP_LEVEL}.10.15.006
   rm $LOG_FILE
   export APP_NAME=ApplicationDriver
   export EXASTAR_HOME=/home/mthavappiragasam/lcf/ExaStar
   #export HDF5_ROOT=${EXASTAR_HOME}/hdf5-1.10.7
   export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.13.3-intel
   export HDF5_INC=${HDF5_ROOT}/include
   export HDF5_LIB=${HDF5_ROOT}/lib
   export THORNADO_DIR=${EXASTAR_HOME}/thornado-bench
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export THORNADO_MACHINE=sunspot_intel

   export IGC_OverrideOCLMaxParamSize=4096
   #export OMP_NUM_THREADS=1
}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   rm ./${APP_NAME}_${THORNADO_MACHINE}

   echo $MKLROOT |& tee -a $LOG_FILE 
   module list   |& tee -a $LOG_FILE

   make clean -f ${THORNADO_DIR}/Makefile
   ( time make -f ${THORNADO_DIR}/Makefile $APP_NAME USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a $LOG_FILE
}

load_set_common
buildApp

