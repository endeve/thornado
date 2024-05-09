#!/bin/bash

function set_common(){

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
