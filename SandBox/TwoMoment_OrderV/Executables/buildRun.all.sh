#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function set_common(){


   export EXASTAR_HOME=/localdisk/quanshao/ExaStar
   export HDF5_INC=${EXASTAR_HOME}/hdf57/include
   export HDF5_LIB=${EXASTAR_HOME}/hdf57/lib64
   export THORNADO_DIR=${EXASTAR_HOME}/thornado-mergingMaster
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib-merging
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export THORNADO_MACHINE=beacon_intel
   export IGC_OverrideOCLMaxParamSize=4096
   export MPIR_CVAR_ENABLE_GPU=0

## for running

   #export MKL_VERBOSE=2
   export LTTNG_HOME=$EXASTAR_HOME
   mkdir -p $LTTNG_HOME
   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   #export LIBOMPTARGET_PLUGIN=LEVEL0
   #export ONEAPI_DEVICE_FILTER=level_zero:gpu
   ##export LIBOMPTARGET_PLUGIN=OPENCL
#   export LIBOMPTARGET_DEBUG=1
   #export EnableImplicitScaling=1
   export ZE_AFFINITY_MASK=0.0
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   #export OMP_TARGET_OFFLOAD=DISABLED
   export OMP_TARGET_OFFLOAD=MANDATORY
   #export OMP_TARGET_OFFLOAD=DISABLED
   #unset OMP_TARGET_OFFLOAD
   export OMP_NUM_THREADS=1
   ulimit -s unlimited
   #ulimit -n 20480
   ## The following seems working well for the SineWaveStream app.
   export LIBOMPTARGET_LEVEL0_MEMORY_POOL=device,128,64,16384
   #export LIBOMPTARGET_LEVEL_ZERO_COMMAND_BATCH=copy,8
   #export OMP_NUM_THREADS=1
}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   echo $MKLROOT |& tee -a $LOG_FILE 
   module list   |& tee -a $LOG_FILE

   make clean
   ( time make -j 16 $APP_NAME ${USER_OPTION} USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a $LOG_FILE
}

###########################################################################################
###  Run Script
###########################################################################################

function runApp(){

   module list |& tee -a $LOG_FILE
# For vtune
##   source /sharedjf/mjh/tools/Intel_VTune_Profiler_2022.3.0_nda/env/vars.sh

# echo some env variables to $LOG_FILE   

   echo "ZE_AFFINITY_MASK="${ZE_AFFINITY_MASK}                                |& tee -a $LOG_FILE
   echo "EnableImplicitScaling="${EnableImplicitScaling}                      |& tee -a $LOG_FILE
   echo "LIBOMPTARGET_LEVEL0_MEMORY_POOL="${LIBOMPTARGET_LEVEL0_MEMORY_POOL}  |& tee -a $LOG_FILE

   if [[ "$ACTION" == "iprof" ]]; then
      module load iprof
      ( time iprof ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE
   elif [[ "$ACTION" == "onetrace" ]]; then
      module use /nfs/pdx/home/roymoore/modules
      module load onetrace
      (time onetrace -h -d  ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE
      #(time onetrace -h -d -v ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE
   elif [[ "$ACTION" == "vtune" ]]; then
      (time (vtune -collect gpu-hotspots -knob characterization-mode=global-local-accesses -data-limit=0 -r ${VT_OUTPUT} ./${APP_NAME}_${THORNADO_MACHINE})) |& tee -a $OUTPUT_LOG
   elif [[ "$ACTION" == "advisor" ]]; then
      module use /nfs/pdx/home/mheckel/modules/modulefiles_nightly
      module load nightly-advisor/23.1.0.613762  ## VERY slow and require old binary. 
      #module load nightly-advisor/23.1.0.613901
      time(advisor --collect=roofline --data-limit=0 --profile-gpu --project-dir=/localdisk/quanshao/ExaStar/thornado-dev/roofline04-03-762 -- ./${APP_NAME}_${THORNADO_MACHINE}) |& tee -a $OUTPUT_LOG
   else
      ( time ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE
   fi

   #(time /nfs/pdx/home/mheckel/pti-gpu/tools/bin/onetrace -h -d -v ./${APP_NAME}_${THORNADO_MACHINE} ) |& tee -a $LOG_FILE
   #( time  gdb-oneapi ./${APP_NAME}_${THORNADO_MACHINE}) |& tee $OUTPUT_LOG
   #valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-fds=yes ./${APP_NAME}_${THORNADO_MACHINE}|& tee -a $OUTPUT_LOG
   #(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -ring-buffer 10 -r $VT_OUTPUT ./${APP_NAME}_${THORNADO_MACHINE}) |& tee -a $OUTPUT_LOG
    #(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -r $VT_OUTPUT ./${APP_NAME}_${THORNADO_MACHINE}) |& tee -a $OUTPUT_LOG
##   vtune-backend --allow-remote-access --enable-server-profiling --reset-passphrase --web-port 8080 --data-directory=${VT_OUTPUT}

   echo "Log file:" $LOG_FILE "writting finished"
}


###########################################################################################
###  Main 
###########################################################################################

module purge

#export A21_SDK_MKLROOT_OVERRIDE=/exaperf/nightly/mkl-cev/2022.11.02 ## Latest nightly, i.e. 10.06, uses this mkl
export A21_SDK_MKLROOT_OVERRIDE=/exaperf/nightly/mkl-cev/2023.04.19 ## Latest nightly, i.e. 2023-05-01 use this mkl 

#export IGC_EnableZEBinary=0
#export IGC_ForceOCLSIMDWidth=16
#export LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=0
#export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=0

ACTION="iprof"
#ACTION="perf"
#ACTION="onetrace"
#ACTION="advisor"
#ACTION="vtune"
ACTION=""
if [[ -n $ACTION ]];then
   faction="-$ACTION"
fi
#export BASE_DATE="2023.03.10"
export BASE_DATE="2023.04.01"
#export COMPILER_DATE="2023.05.03"
#export COMPILER_DATE="2023.05.08"
export COMPILER_DATE="2023.05.15"
#export COMPILER_DATE="2023.06.22"

#export COMPILER_DATE="2023.03.30"
#export COMPILER_DATE="2023.06.11"
#export COMPILER_DATE="2023.06.07"
#export COMPILER_DATE="2023.06.05"
#export COMPILER_DATE="2023.06.06"
#export COMPILER_DATE="2023.06.04"
#export COMPILER_DATE="2023.05.21"
#export COMPILER_DATE="2023.05.22"
#export COMPILER_DATE="2023.05.23"
#export COMPILER_DATE="2023.05.24"
#export COMPILER_DATE="2023.05.29"
#export COMPILER_DATE="2023.06.01"
#export COMPILER_DATE="2023.03.10"
#export COMPILER_DATE="2022.12.30.002"
export AADEBUG=""

#export ONEAPI_MODULE_OVERRIDE=oneapi/eng-compiler/2022.12.30.003
module load nightly-compiler/${COMPILER_DATE}

#UMD=""
#UMD="neo/agama-devel-sp3/611-23.09.25812.14-609"
#UMD="neo/agama-devel-sp3/609-23.09.25812.14-609"
#UMD="neo/agama-devel-sp3/610-23.09.25812.14-609"
#UMD="neo/agama-devel-sp3/608-23.05.25593.18-i606"
#UMD="neo/agama-devel-sp3/619-23.09.25812.15-619"
#UMD="neo/agama-devel-sp3/625-23.13.26032.7-624"
#UMD="neo/agama-devel-sp3/622-23.09.25812.15-622"
#UMD="neo/agama-devel-sp3/628-23.13.26032.8-626"
#UMD="neo/agama-devel-sp3/627-23.13.26032.8-626"
#UMD="neo/agama-devel-sp3/636-23.13.26032.22-631"
#UMD="neo/agama-devel-sp3/637-23.13.26032.26-637"
#UMD="neo/agama-devel-sp3/639-23.13.26032.26-637"
#UMD="neo/agama-devel-sp3/644-23.13.26032.30-644"
#UMD="neo/agama-devel-sp3/646-23.13.26032.30-646"
#UMD="neo/agama-devel-sp3/648-23.17.26241.13-648"
#UMD="neo/agama-devel-sp3/650-23.17.26241.13-649"
#UMD="neo/agama-devel-sp3/651-23.17.26241.14-651"
#UMD="neo/agama-devel-sp3/652-23.17.26241.14-651"
#UMD="neo/agama-devel-sp3/653-23.17.26241.15-653"
#UMD="intel_compute_runtime/release/agama-devel-602"
#UMD="neo/agama-devel-sp3/657-23.17.26241.18-657"
#UMD="neo/agama-devel-sp3/658-23.17.26241.19-658"
#UMD="neo/agama-devel-sp3/661-23.17.26241.21-661"
#UMD="neo/agama-devel-sp3/663-23.17.26241.21-662"
#UMD="neo/agama-devel-sp3/664-23.17.26241.21-664"
#UMD="neo/agama-devel-sp3/666-23.17.26241.22-665"

#UMD="neo/agama-devel-sp3/627-23.13.26032.8-626"
#UMD="neo/agama-devel-sp3/670-23.17.26241.24-670"
#UMD="neo/agama-devel-sp3/671-23.22.26516.8-671"
#UMD="neo/agama-devel-sp3/672-23.22.26516.8-672"
#UMD="neo/agama-devel-sp3/673-23.22.26516.8-673"
UMD="neo/agama-devel-sp3/674-23.22.26516.8-673"
#UMD="neo/agama-devel-sp3/673-23.22.26516.8-673"
umdf=""
export useAGRF="FALSE"
if [[ -n $UMD ]]; then
   module switch -f intel_compute_runtime/release/agama-devel-551 $UMD
   if [[ $UMD == intel_compute* ]]; then
      umdf="-devel602"
      export useAGRF="FALSE"
   else
      umdf=`echo $UMD |cut -d '/' -f3`
      umdf=`echo $umdf |cut -d '-' -f1`
      umdf="-umd$umdf"
#      umdf="-Tdebug1IMM1-umd$umdf"
      export useAGRF="TRUE"
   fi
fi

#export COMPILER_DATE="2023.05.03"
#module load nightly-compiler/${COMPILER_DATE}
#export COMPILER_DATE=".2023.05.15.003-rc11"
#module load oneapi/eng-compiler/${COMPILER_DATE}
#export FI_PROVIDER=sockets
#export useAGRF="TRUE"
#module switch -f intel_compute_runtime/release/agama-devel-627 neo/agama-devel-sp3/627-23.13.26032.8-626

#module switch -f intel_compute_runtime/release/agama-devel-551 neo/agama-devel-sp3/627-23.13.26032.8-626 
#module switch -f intel_compute_runtime/release/agama-devel-551 intel_compute_runtime/release/agama-devel-627 
#umdf="umd627NoIMMmkl2"

#if action is empty, performance comparison will be done. otherwise there is no performance comparison and just run the app using such as onetrace, vtune etc. so action can be "", "onetrace", "iprof", "vtune", 

#opLevels=(O3)
#grids=("[8,8,8]")
#gridNames=("")
#appNames=(ApplicationDriver)
#logFiles=(sineWave)
#CaseNames=(SineWaveStreaming)
#userOptions=("")
#gridLines=(85)

#opLevels=(O3)
#grids=("[8,8,8]")
#gridNames=("")
#appNames=(ApplicationDriver_Neutrinos)
#logFiles=(relax)
#CaseNames=(Relaxation)
#userOptions=("MICROPHYSICS=WEAKLIB")
#gridLines=(127)

#opLevels=(O0 O1 O2 O3)
opLevels=(O3)
grids=("[8,8,8]" "[16,16,16]")
gridNames=("" "-xN16")
appNames=(ApplicationDriver ApplicationDriver_Neutrinos)
logFiles=(sineWave relax)
CaseNames=(SineWaveStreaming Relaxation)
userOptions=("" "MICROPHYSICS=WEAKLIB")
gridLines=(85 127)

#grids=("[8,8,8]")
#gridNames=("")
#appNames=(ApplicationDriver)
#logFiles=(sineWave)
#CaseNames=(SineWaveStreaming)
#userOptions=("")
#gridLines=(85)
#appNames=(ApplicationDriver_Neutrinos)
#logFiles=(relax)
#CaseNames=(Relaxation)
#userOptions=("MICROPHYSICS=WEAKLIB")
#gridLines=(127)

if [[ "$ACTION" == "vtune" ]]; then
   opLevels=(O3)
   grids=("[8,8,8]")
   gridNames=("")
   appNames=(ApplicationDriver)
   logFiles=(sineWave)
   CaseNames=(SineWaveStreaming)
   userOptions=("")
   gridLines=(85)

fi

set_common

timeFOMLog="timeFOM_${COMPILER_DATE}.txt${umdf}$AADEBUG"
if [[ -z $ACTION ]];then
   rm -rf $timeFOMLog
   echo "                                             Time(seconds)                         |              Figure of Merit (FOM)">>$timeFOMLog
   echo "AppName     Grid      OpLevel :  ${COMPILER_DATE}   ${BASE_DATE}    TimeDiff   Percentage  |    ${COMPILER_DATE}   ${BASE_DATE}    FOM-Diff   Percentage">>$timeFOMLog
   echo "-----------------------------    ------------------------------------------------       ------------------------------------------------">>$timeFOMLog
fi

for ((jj=0; jj<${#appNames[@]}; jj++));
do
   export APP_NAME=${appNames[jj]}
   for ((ii=0; ii<${#grids[@]}; ii++));
   do

      sed -i "${gridLines[jj]}s/.*/      nX  =${grids[ii]}/" ../${appNames[jj]}.F90
      for op in "${opLevels[@]}";
      do 
         if [[ "$op" == "O0" && "${grids[ii]}" == "[16,16,16]" ]];then
            continue
         fi

         export OP_LEVEL=$op
         export LOG_FILE=${logFiles[jj]}.${OP_LEVEL}.${COMPILER_DATE}.ms69${umdf}${gridNames[ii]}${faction}$AADEBUG
#         export LOG_BASE=${logFiles[jj]}.${OP_LEVEL}.${BASE_DATE}.ms69-umd602${gridNames[ii]}$AADEBUG
         export LOG_BASE=${logFiles[jj]}.${OP_LEVEL}.${BASE_DATE}.ms69-devel602${gridNames[ii]}$AADEBUG
         #export LOG_BASE=${logFiles[jj]}.${OP_LEVEL}.${BASE_DATE}.ms69${umdf}${gridNames[ii]}$AADEBUG
         export USER_OPTION=${userOptions[jj]}

         if [[ "$ACTION" == "vtune" ]]; then
            export VT_OUTPUT=vtune_${logFiles[jj]}.${COMPILER_DATE}${umdf}
            rm -rf $VT_OUTPUT
         fi

         if [[ "$1" == "FOM" ]];then
            echo ""
            echo "Gathering performance data from $LOG_FILE and $LOG_BASE"
            echo ""
         else 
            #echo $USER_OPTION
            echo "Building and/or running" ${logFiles[jj]} "using Op-level "${OP_LEVEL} 
            rm $LOG_FILE

            if [[ "$1" == -[rR]* ]]; then
               echo ""
               echo "Runing ${CaseNames[jj]} ..."
               echo ""
               if [ -f "${APP_NAME}_${THORNADO_MACHINE}" ];then
                  echo "$op ${grids[ii]} ${appNames[jj]}"
                  runApp
               else
                  echo "The executable does not exist", ${APP_NAME}_${THORNADO_MACHINE}
               fi
            elif [[ "$1" == -[bB]* ]]; then
               echo ""
               echo "Compiling ${CaseNames[jj]} ..."
               echo ""
               rm ${APP_NAME}_${THORNADO_MACHINE}
               buildApp
            else
               rm ${APP_NAME}_${THORNADO_MACHINE}
               echo ""
               echo "Compiling and Running ${CaseNames[jj]} ..."
               echo ""
               buildApp
               if [ -f "${APP_NAME}_${THORNADO_MACHINE}" ];then
                  echo "$op ${grids[ii]} ${appNames[jj]}"
                  runApp
               else
                  echo "The executable does not exist", ${APP_NAME}_${THORNADO_MACHINE}
               fi
            fi   
         fi

         ## compare IMEX_TIME to the BASE_DATE
         if [[ -z $ACTION ]];then
            baseTime=`grep Timer_IMEX $LOG_BASE |cut -d':' -f2`
            baseTime=`echo $baseTime |cut -d ' ' -f1`
            baseTime=`printf "%.6f" $baseTime`
            currTime=`grep Timer_IMEX $LOG_FILE |cut -d':' -f2`
            currTime=`echo $currTime |cut -d ' ' -f1`
            currTime=`printf "%.6f" $currTime`
            diffTime=`echo ${currTime}-${baseTime}|bc -l`
            percentage=`echo 100*${diffTime}/${baseTime}|bc -l`

            percentage=`printf "%8.2f" $percentage`
            currTime=`printf "%12.4e" $currTime`
            baseTime=`printf "%12.4e" $baseTime`
            diffTime=`printf "%12.4e" $diffTime`

            baseFOM=`grep "FOM is" $LOG_BASE |cut -d':' -f2`
            baseFOM=`printf "%.6f" $baseFOM`
            currFOM=`grep "FOM is" $LOG_FILE |cut -d':' -f2`
            currFOM=`printf "%.6f" $currFOM`
            diffFOM=`echo ${currFOM}-${baseFOM}|bc -l`
            percentFOM=`echo 100*${diffFOM}/${baseFOM}|bc -l`

            percentFOM=`printf "%8.2f" $percentFOM`
            baseFOM=`printf "%12.4e" $baseFOM`
            currFOM=`printf "%12.4e" $currFOM`
            diffFOM=`printf "%12.4e" $diffFOM`

            caseName=`printf "%-10s" ${logFiles[jj]}`
            gg=`printf "%-10s" ${grids[ii]}`
            echo "$caseName $gg   $OP_LEVEL    :$currTime $baseTime $diffTime $percentage%       $currFOM $baseFOM $diffFOM $percentFOM%" >>$timeFOMLog
         fi
   done
done
done

if [[ -z $ACTION ]];then
   echo
   echo " Performance Comparison between compiler $COMPILER_DATE and $BASE_DATE"
   echo 

   echo "cat $timeFOMLog"
   cat $timeFOMLog
fi
