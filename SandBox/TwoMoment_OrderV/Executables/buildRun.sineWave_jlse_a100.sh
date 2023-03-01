#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function load_set_common(){

   module use /soft/modulefiles
   module use -p  /soft/restricted/CNDA/modulefiles

   #module purge
   #module use /soft/restricted/CNDA/modulefiles

   module load nvhpc/nvhpc/22.7

   mr=${1:-1}
   date0=`date +%m-%d-%Y`
   hname=`hostname`

   export OP_LEVEL=O3
   export LOG_FILE=sw_a100.${hname}_mpi${mr}_${date0}.log
   rm -f $LOG_FILE
   export APP_NAME=ApplicationDriver
   export EXASTAR_HOME=/home/mthavappiragasam/lcf/ExaStar/
   export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.12.0-nvhpc
   #export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.10.7-a100
   export HDF5_LIB=${HDF5_ROOT}/lib 
   export HDF5_INC=${HDF5_ROOT}/include
   export THORNADO_DIR=${EXASTAR_HOME}/thornado
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export CUDA_ROOT=/soft/compilers/nvhpc/Linux_x86_64/22.7/compilers
   export NETLIB_LAPACK_ROOT=/soft/compilers/nvhpc/Linux_x86_64/22.7/compilers/lib
   export THORNADO_MACHINE=jlse_A100_nvhpc
   #export OMP_NUM_THREADS=1
}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   rm ./${APP_NAME}_${THORNADO_MACHINE}
   
   module list   |& tee -a $LOG_FILE

   make clean
   ( time make -j8 $APP_NAME USE_OMP_OL=TRUE  USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE ) |& tee -a $LOG_FILE
   #( time make -j8 $APP_NAME USE_OACC=TRUE  USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE ) |& tee -a $LOG_FILE
   #( time make -j8 $APP_NAME MICROPHYSICS=WEAKLIB USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE ) |& tee -a $LOG_FILE

}

###########################################################################################
###  Run Script
###########################################################################################

function runApp(){

   export PGI_ACC_DEBUG=0
   #export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   ##export LIBOMPTARGET_PLUGIN=OPENCL
   #export LIBOMPTARGET_DEBUG=0
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   #export OMP_TARGET_OFFLOAD=DISABLED
   #export OMP_TARGET_OFFLOAD=MANDATORY
   #export OMP_NUM_THREADS=1

   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   rm -rf ../Output
   mkdir  ../Output
   ulimit -s unlimited

   module list |& tee -a $LOG_FILE

   nsys0="nsys profile --verbose --stats=true -t mpi,openacc,cudnn,cuda,cublas,osrt,nvtx -o ${logfile}.nsys --force-overwrite true  "
   #nsys0="iprof"
   #nsys0=""

   #export LIBOPMTARGET_INFO=2
    #( time mpirun -n ${mr} --bind-to core --mca btl openib,vader,self  -n 1 ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}
    #( time  nvprof ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}_nvprof
    ( time $nsys0 mpirun -n ${mr} --bind-to core ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a ${LOG_FILE}

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

