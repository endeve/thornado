#!/bin/bash
###
# A script for testing/evaluating the results of Thornado on sunspot.alcf.anl apps streaming sinewave(ssw) and relaxation
# Argument-options: 
# 0. no arguments tests the both apps
# 1. -s* or -S* to run the app ssw only
# 2. -r* or -R* to run the app relaxation only
# 3. -o* or O* to collect numbers for already obtained results  
###
               
function set_common(){

   #module purge

   #module load oneapi/eng-compiler/2022.01.30.005
   #module load oneapi/eng-compiler/2022.06.30.002
   #module load oneapi/eng-compiler/2022.10.15.006
   #module restore
   #module load oneapi/eng-compiler/2023.05.15.007
   #module load oneapi/release/2023.12.15.001 #Seems losing perfomance
   module load oneapi/eng-compiler/2023.10.15.002
   #module load mpich/51.2/icc-all-pmix-gpu

   export EXASTAR_HOME=/home/mthavappiragasam/lcf/ExaStar
   #export HDF5_ROOT=${EXASTAR_HOME}/hdf5-1.10.7
   export HDF5_ROOT=${EXASTAR_HOME}/hdf5/1.13.3-intel
   export HDF5_INC=${HDF5_ROOT}/include
   export HDF5_LIB=${HDF5_ROOT}/lib
   export THORNADO_DIR=${EXASTAR_HOME}/thornado-bench
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export THORNADO_MACHINE=sunspot_intel
   #export THORNADO_MACHINE=sunspot.alcf.anl

   export IGC_OverrideOCLMaxParamSize=4096
}

function runApp(){
   module list |& tee -a $LOG_FILE

   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   export LIBOMPTARGET_PLUGIN=LEVEL0
   export LIBOMPTARGET_DEBUG=0
   export EnableImplicitScaling=0
   export ZE_AFFINITY_MASK=0.0
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   #export OMP_TARGET_OFFLOAD=DISABLED
   export OMP_TARGET_OFFLOAD=MANDATORY
   export OMP_NUM_THREADS=1
   ulimit -s unlimited
   #ulimit -n 20480
   export LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL=device,256,128,32768

   # echo some env variables to $LOG_FILE
   echo "ZE_AFFINITY_MASK="${ZE_AFFINITY_MASK}                                |& tee -a $LOG_FILE
   echo "EnableImplicitScaling="${EnableImplicitScaling}                      |& tee -a $LOG_FILE
   echo "LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL="${LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL}  |& tee -a $LOG_FILE
   echo "IGC_EnableZEBinary="${IGC_EnableZEBinary}
   ( time  ./${APP_NAME}_${THORNADO_MACHINE} $nX ) |& tee -a ${LOG_FILE}

}
                   
grids=("16 16 16") 
gridNames=("16")
appNames=(ApplicationDriver ApplicationDriver_Neutrinos)
theoreticalNames=(streamingSineWave relaxation)
error_order_relax=(-13 -13 -9 -9 -9 -9)
error_order_ssw_16=(-16 -16 -16 -16 -16 -16)
#error_order_ssw_8=(-2 -16)
   
function run_app(){
   
   export APP_NAME=${appNames[jj]}
   echo " Running ${logFiles[jj]} "
   
   if [ -f "${APP_NAME}_${THORNADO_MACHINE}" ];then
                echo " ${appNames[jj]}"
                #runApp
            else
                echo "The executable does not exist", ${APP_NAME}_${THORNADO_MACHINE}
   fi
}


function test_ssw(){
        jj=0
        run_app
        ri=0
        inf_error=`grep 'Sp ' sineWave.16 -A 6 | tail -6 | grep -o '...$'`
        #while read line; do echo "LINE: '${line}'"; done <<< "${inf_error}"
                while read line;
                do
                        if [ ${error_order_ssw_16[$ri]} -lt ${line} ]; then
                                echo "Result seems wrong for the species '$ri', order of the infinity-norm error is ${line} but it should be '${error_order_ssw_16[$ri]}'"
                                break
                        fi
                        ri=`expr $ri + 1`
                        #echo "LINE: '${line}'"
                done <<< "${inf_error}"

                if [ $ri -ge 6 ]; then
                   echo " Test is successful for the app, '${theoreticalNames[jj]}'"
                   else
                   echo " Something wrong with the test app, '${theoreticalNames[jj]}'"
                fi
}

function test_relax(){
        jj=1
        run_app
        ri=0
        inf_error=`grep 'Inf Error' relax.16 | grep -o '...$'`
        #while read line; do echo "LINE: '${line}'"; done <<< "${inf_error}"
                while read line;
                do
                        if [ ${error_order_relax[$ri]} -lt ${line} ]; then
                                echo "Result seems wrong for the species '$ri', order of the infinity-norm error is ${line} but it should be '${error_order_relax[$ri]}'"
                                break
                        fi
                        ri=`expr $ri + 1`
                        #echo "LINE: '${line}'"
                done <<< "${inf_error}"

                if [ $ri -ge 6 ]; then
                   echo " Test is successful for the app, '${theoreticalNames[jj]}'"
                   else
                   echo " Something wrong with the test app, '${theoreticalNames[jj]}'"
                fi
}

set_common

if [[ "$1" == -[sS]* ]]; then
   test_ssw
elif [[ "$1" == -[rR]* ]]; then
   test_relax
elif [[ "$1" == -[aA]* ]]; then
   test_ssw
   test_relax
else
   echo " Pass a valued argument '-[sS]*' to test ssw only OR '-[rR]*' to test relaxation only  OR '-[aA]*' to tes the both apps "
fi
