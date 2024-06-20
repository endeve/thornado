#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function load_set_common(){


   #module purge
  
   #module load nvhpc/21.9
   module load nvhpc/23.3
   #module restore

   mr=${1:-1}
   date0=`date +%m-%d-%Y`
   hname=`hostname`

   export OP_LEVEL=O3
   export LOG_FILE=sw_a100.${hname}_mpi${mr}_${date0}.log
   rm -f $LOG_FILE
   export APP_NAME=ApplicationDriver
   export EXASTAR_HOME=/home/mthavappiragasam/lcf/ExaStar
   export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.12.0-nvhpc
   #export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.10.7-a100
   export HDF5_LIB=${HDF5_ROOT}/lib 
   export HDF5_INC=${HDF5_ROOT}/include
   #export THORNADO_DIR=${EXASTAR_HOME}/thornado-new
   export THORNADO_DIR=/grand/projects/Performance/mathit/ExaStar/thornado
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export CUDA_ROOT=${NVIDIA_PATH}/cuda
   #export CUDA_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/cuda
   #export CUDA_ROOT=/soft/compilers/nvhpc/Linux_x86_64/22.1/compilers
   export NETLIB_LAPACK_ROOT=${NVIDIA_PATH}/math_libs/lib64
   #export NETLIB_LAPACK_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib
   export MAGMA_ROOT=/soft/libraries/magma-2.7.0
   export THORNADO_MACHINE=polaris_nvhpc
   #export OMP_NUM_THREADS=1
}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   rm ./${APP_NAME}_${THORNADO_MACHINE}
   
   module list   |& tee -a $LOG_FILE

   make clean
   ( time make -j8 $APP_NAME USE_OMP_OL=TRUE  USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE USE_ONEMKL=FALSE) |& tee -a $LOG_FILE
   #( time make -j8 $APP_NAME USE_OACC=TRUE  USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE USE_ONEMKL=FALSE) |& tee -a $LOG_FILE
   #( time make -j8 $APP_NAME MICROPHYSICS=WEAKLIB USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE ) |& tee -a $LOG_FILE

}

###########################################################################################
###  Run Script
###########################################################################################

function runApp(){
   mr=${1:-1}

   #nsys0="nsys profile --verbose --stats=true -t mpi,openacc,cudnn,cuda,cublas,osrt,nvtx -o ${logfile}.nsys --force-overwrite true  "
   nsys0=""
   #nsys0="ncu -f -o $THORNADO_DIR/ssw_small_10 --set full"
   #nsys0="ncu -f -o $THORNADO_DIR/ssw_16 --set full --target-processes all"

   #export PGI_ACC_DEBUG=1
   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   ##export LIBOMPTARGET_PLUGIN=OPENCL
   #export LIBOMPTARGET_DEBUG=0
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   #export OMP_TARGET_OFFLOAD=DISABLED
   #export OMP_TARGET_OFFLOAD=MANDATORY
   #export OMP_NUM_THREADS=1

   module list |& tee -a $LOG_FILE

   #export LIBOPMTARGET_INFO=2
    #( time $nsys0 mpiexec -n ${mr} --bind-to core ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}
    #( time $nsys0 mpirun -n ${mr} --bind-to core ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}
   # ( time mpiexec -n 1  --cpu-bind depth  ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}
    #( time mpiexec -n ${mr} --bind-to core --mca btl openib,vader,self  -n 1 ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}
    #( time  ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}_nvprof
    ( time  $nsys0  ./${APP_NAME}_${THORNADO_MACHINE} 8 8 8 ) |& tee -a ${LOG_FILE}

}

###########################################################################################
###  Application Compile and Run
###########################################################################################


load_set_common
if [[ "$1" == -[rR]* ]]; then
   if [ -f "${APP_NAME}_${THORNADO_MACHINE}" ];then
      runApp
   else
      echo "The executable does not exist", ${APP_NAME}_${THORNADO_MACHINE}
   fi
elif [[ "$1" == -[bB]* ]]; then   
   buildApp
else
   buildApp
   if [ -f "${APP_NAME}_${THORNADO_MACHINE}" ];then
      runApp
   else
      echo "The executable does not exist", ${APP_NAME}_${THORNADO_MACHINE}
   fi
fi

