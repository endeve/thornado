#!/bin/bash

function set_common(){
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

function runApp(){
   module list |& tee -a $LOG_FILE
   # echo some env variables to $LOG_FILE   

   echo "ZE_AFFINITY_MASK="${ZE_AFFINITY_MASK}                                |& tee -a $LOG_FILE
   echo "EnableImplicitScaling="${EnableImplicitScaling}                      |& tee -a $LOG_FILE
   echo "LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL="${LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL}  |& tee -a $LOG_FILE
   echo "IGC_EnableZEBinary="${IGC_EnableZEBinary}
   ( time  ./${APP_NAME}_${THORNADO_MACHINE} $nX ) |& tee -a ${LOG_FILE}

}

set_common

grids=("8 8 8" "16 16 16")
gridNames=("8", "16")
appNames=(ApplicationDriver ApplicationDriver_Neutrinos)
logFiles=(sineWave relax)

   timeFOMLog="timeFOM_${THORNADO_MACHINE}.txt"
   rm -rf $timeFOMLog

   echo "-----------------------------------------------------------">>$timeFOMLog
   echo "|      AppName  |   Grid    |   Time(sec)  |    FOM       |">>$timeFOMLog
   echo "-----------------------------------------------------------">>$timeFOMLog

for ((jj=0; jj<${#appNames[@]}; jj++));
do
   export APP_NAME=${appNames[jj]}
   for ((ii=0; ii<${#grids[@]}; ii++));
   do
        nX=${grids[ii]}
        export LOG_FILE=${logFiles[jj]}_${gridNames[ii]}.${THORNADO_MACHINE}
        if [[ "$1" == -[rR]* ]]; then
            echo ""
            echo "Gathering performance data from $LOG_FILE"
            echo ""
        else
            rm $LOG_FILE
            echo " Running ${logFiles[jj]} "

            if [ -f "${APP_NAME}_${THORNADO_MACHINE}" ];then
                echo "${grids[ii]} ${appNames[jj]}"
                runApp
            else
                echo "The executable does not exist", ${APP_NAME}_${THORNADO_MACHINE}
            fi
        fi
        currTime=`grep IMEX_Time $LOG_FILE |cut -d':' -f2`
        ## currTime=`echo $currTime |cut -d ' ' -f1`
        currTime=`printf "%.6f" $currTime`
        currTime=`printf "%12.4e" $currTime`

        currFOM=`grep "FOM" $LOG_FILE |cut -d':' -f2`
        currFOM=`printf "%.6f" $currFOM`
        currFOM=`printf "%12.4e" $currFOM`

        caseName=`printf "%-10s" ${logFiles[jj]}`
        gg=`printf "%-3s" ${grids[ii]}`
        echo "|     $caseName | $gg | $currTime | $currFOM |" >>$timeFOMLog
        echo "-----------------------------------------------------------">>$timeFOMLog

   done
done
echo "cat $timeFOMLog"
cat $timeFOMLog
