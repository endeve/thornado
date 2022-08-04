#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function load_set_common(){

   module purge
   
   #module load oneapi/eng-compiler/2022.01.30.005
   #module load oneapi/eng-compiler/2022.01.30.006
   #module load oneapi/eng-compiler/2022.01.30.007
   #module load oneapi/eng-compiler/2022.06.30.002
   #module load nightly-compiler/2022.07.05
   #module load nightly-compiler/2022.07.26
   #module load nightly-compiler/2022.07.27
   #module load nightly-compiler/2022.08.02
   module load nightly-compiler/2022.08.03
   #module switch -f intel_compute_runtime/release/agama-prerelease-475 neo/agama-prerelease/531-22.28.023672.11691-main
   #module switch -f intel_compute_runtime/release/agama-prerelease-475 neo/agama-prerelease/534-22.28.023672.11688-main
   #module switch -f intel_compute_runtime/release/agama-prerelease-475 neo/agama-prerelease/535-22.28.023711.11742-main
   #module switch -f intel_compute_runtime/release/agama-prerelease-475 neo/agama-prerelease/546-22.30.023803.11807-main
   #module switch -f intel_compute_runtime/release/agama-prerelease-475 neo/agama-prerelease/552-22.30.023828.11834-main
   module switch -f intel_compute_runtime/release/agama-prerelease-475 neo/agama-prerelease/554-22.31.023852.11876-main


   export OP_LEVEL=O1
   export LOG_FILE=sineWave.${OP_LEVEL}.0803.umd554
   #export LOG_FILE=sineWave.${OP_LEVEL}.eng06.30.002
   rm $LOG_FILE
   export APP_NAME=ApplicationDriver
   export EXASTAR_HOME=/localdisk/quanshao
   export HDF5_INC=${EXASTAR_HOME}/ExaStar/hdf57/include
   export HDF5_LIB=${EXASTAR_HOME}/ExaStar/hdf57/lib64
   export THORNADO_DIR=${EXASTAR_HOME}/ExaStar/thornado
   export WEAKLIB_DIR=${EXASTAR_HOME}/ExaStar/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/ExaStar/weaklib-tables
   export THORNADO_MACHINE=beacon_intel
   export IGC_OverrideOCLMaxParamSize=4096
   export MPIR_CVAR_ENABLE_GPU=0
##   export LIBOMPTARGET_LEVEL_ZERO_COMMAND_BATCH=copy,8
   #export OMP_NUM_THREADS=1
}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   rm ./${APP_NAME}_${THORNADO_MACHINE}

   echo $MKLROOT |& tee -a $LOG_FILE 
   module list   |& tee -a $LOG_FILE

   make clean
   ( time make $APP_NAME USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a $LOG_FILE
}

###########################################################################################
###  Run Script
###########################################################################################

function runApp(){

   module load iprof

   export LTTNG_HOME=$EXASTAR_HOME
   mkdir -p $LTTNG_HOME
   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   export LIBOMPTARGET_PLUGIN=LEVEL0
   export SYCL_DEVICE_FILTER=LEVEL_ZERO
   ##export LIBOMPTARGET_PLUGIN=OPENCL
   export LIBOMPTARGET_DEBUG=0
   export EnableImplicitScaling=1
   export ZE_AFFINITY_MASK=0
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   #export OMP_TARGET_OFFLOAD=DISABLED
   export OMP_TARGET_OFFLOAD=MANDATORY
   #export OMP_TARGET_OFFLOAD=DISABLED
   #unset OMP_TARGET_OFFLOAD
   export OMP_NUM_THREADS=1
   ulimit -s unlimited
   #ulimit -n 20480
   ## The following seems working well for the SineWaveStream app.
   export LIBOMPTARGET_LEVEL0_MEMORY_POOL=device,16,32

   module list |& tee -a $LOG_FILE

# For vtune
   source /sharedjf/mjh/tools/Intel_VTune_Profiler_2022.3.0_nda/env/vars.sh
   VT_OUTPUT=vtune07June2022
   rm -rf $VT_OUTPUT

# echo some env variables to $LOG_FILE   

   echo "ZE_AFFINITY_MASK="${ZE_AFFINITY_MASK}                                |& tee -a $LOG_FILE
   echo "EnableImplicitScaling="${EnableImplicitScaling}                      |& tee -a $LOG_FILE
   echo "LIBOMPTARGET_LEVEL0_MEMORY_POOL="${LIBOMPTARGET_LEVEL0_MEMORY_POOL}  |& tee -a $LOG_FILE


   ( time ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE
   #( time iprof ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE

   #( time  gdb-oneapi ./${APP_NAME}_${THORNADO_MACHINE}) |& tee $LOG_FILE
   #valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-fds=yes ./${APP_NAME}_${THORNADO_MACHINE}|& tee -a $LOG_FILE
   #vtune -collect gpu-hotspots -knob characterization-mode=global-local-accesses -r vtune-06032022-02 ./${APP_NAME}_${THORNADO_MACHINE} |& tee -a $LOG_FILE
   #(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -ring-buffer 10 -r $VT_OUTPUT ./${APP_NAME}_${THORNADO_MACHINE}) |& tee -a $LOG_FILE
    #(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -r $VT_OUTPUT ./${APP_NAME}_${THORNADO_MACHINE}) |& tee -a $LOG_FILE

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

