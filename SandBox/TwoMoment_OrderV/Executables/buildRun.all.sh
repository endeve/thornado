#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function set_common(){


   export EXASTAR_HOME=/localdisk/quanshao/ExaStar
   export HDF5_INC=${EXASTAR_HOME}/hdf57/include
   export HDF5_LIB=${EXASTAR_HOME}/hdf57/lib64
   export THORNADO_DIR=${EXASTAR_HOME}/thornado-nre
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/toAurora/weaklib-tables
   export THORNADO_MACHINE=beacon_intel

   export IGC_OverrideOCLMaxParamSize=4096
   export numTeamCap="-mllvm -vpo-paropt-atomic-free-red-global-buf-size=4096"

## for running

#   export MPIR_CVAR_ENABLE_GPU=0
   #export MKL_VERBOSE=2
   export LTTNG_HOME=$EXASTAR_HOME
   mkdir -p $LTTNG_HOME
   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   #export LIBOMPTARGET_PLUGIN=LEVEL0
   #export ONEAPI_DEVICE_FILTER=level_zero:gpu
   #export LIBOMPTARGET_DEBUG=4
   #export EnableImplicitScaling=1
   export ZE_AFFINITY_MASK=0.0
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   export OMP_TARGET_OFFLOAD=MANDATORY
   export OMP_NUM_THREADS=1
   ulimit -s unlimited
   #ulimit -n 20480
   export LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL=device,128,64,16384
   export FI_PROVIDER=sockets

}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   echo $MKLROOT |& tee -a $LOG_FILE 
   module list   |& tee -a $LOG_FILE

   make clean
   ( time make -j 16 $APP_NAME ${USER_OPTION} USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a $LOG_FILE
   mv ${APP_NAME}_${THORNADO_MACHINE} $exeName
   echo "${exeName} has been build" 
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
   echo "LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL="${LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL}  |& tee -a $LOG_FILE
   echo "IGC_EnableZEBinary="${IGC_EnableZEBinary}

   if [[ "$ACTION" == "iprof" ]]; then
      module load iprof/0.11.2
      ( time iprof ./${exeName} ) |& tee -a $LOG_FILE
   elif [[ "$ACTION" == "onetrace" ]]; then
      module use /nfs/pdx/home/roymoore/modules
      module load onetrace
      (time onetrace -h -d  ./${exeName} ) |& tee -a $LOG_FILE
      #(time onetrace -h -d -v ./${exeName} ) |& tee -a $LOG_FILE
   elif [[ "$ACTION" == "vtune" ]]; then
      set -x
       vtune -version
      #(time (vtune -collect gpu-hotspots -knob characterization-mode=global-local-accesses -data-limit=0 -r ${VT_OUTPUT} ./${exeName})) |& tee -a $OUTPUT_LOG
      (time (vtune -collect gpu-hotspots -data-limit=0 -r ${VT_OUTPUT} ./${exeName})) |& tee -a $OUTPUT_LOG
      set +x
   elif [[ "$ACTION" == "advisor" ]]; then
      module use /nfs/pdx/home/mheckel/modules/modulefiles_nightly
      module load nightly-advisor/23.2.0.614354
      #module load nightly-advisor/23.1.0.613762  ## VERY slow and require old binary. 
      #module load nightly-advisor/23.1.0.613901
      time(advisor --collect=roofline --data-limit=0 --profile-gpu --project-dir=/localdisk/quanshao/ExaStar/thornado-ms69/SandBox/TwoMoment_OrderV/Executables/advisor-sineWave-umd692 -- ./${exeName}) |& tee -a $OUTPUT_LOG
   else
#      echo "FI_PROVIDER="$FI_PROVIDER
#      gdb-oneapi ./ApplicationDriver_beacon_intel
      ( time ./${exeName} ) |& tee -a $LOG_FILE
   fi

   #(time /nfs/pdx/home/mheckel/pti-gpu/tools/bin/onetrace -h -d -v ./${exeName} ) |& tee -a $LOG_FILE
   #( time  gdb-oneapi ./${exeName}) |& tee $OUTPUT_LOG
   #valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-fds=yes ./${exeName}|& tee -a $OUTPUT_LOG
   #(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -ring-buffer 10 -r $VT_OUTPUT ./${exeName}) |& tee -a $OUTPUT_LOG
    #(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -r $VT_OUTPUT ./${exeName}) |& tee -a $OUTPUT_LOG
##   vtune-backend --allow-remote-access --enable-server-profiling --reset-passphrase --web-port 8080 --data-directory=${VT_OUTPUT}

   echo "Log file:" $LOG_FILE "writting finished"
}


###########################################################################################
###  Main 
###########################################################################################

#ACTION="iprof"
#ACTION="perf"
#ACTION="onetrace"
#ACTION="advisor"
#ACTION="vtune"
#export ONEAPI_MODULE_OVERRIDE=oneapi/eng-compiler/2023.05.15.007
#ACTION=""

if [[ -n $ACTION ]];then
   faction="-$ACTION"
fi

BASE_CMP_DATE="eng-23/05.15.007"
#BASE_UMD="-dev682.20"
export AADEBUG=""
export useAGRF="TRUE"
module purge
module load oneapi/eng-compiler/2023.10.15.002
COMPILER_DATE="eng-23.10.15.002" 
umdf=""

#if action is empty, performance comparison will be done. otherwise there is no performance comparison and just run the app using such as onetrace, vtune etc. so action can be "", "onetrace", "iprof", "vtune", 
#opLevels=(O0 O1 O2 O3)
opLevels=(O3)
#grids=("[64, 1, 1]" "[8, 8, 8]" "[16, 16, 16]")
#gridNames=("64-1-1" "" "-xN16")
grids=( "[8, 8, 8]" "[16, 16, 16]")
gridNames=("-xN8" "-xN16")
#grids=( "[8, 8, 8]")
#gridNames=("")
appNames=(ApplicationDriver ApplicationDriver_Neutrinos)
logFiles=(sineWaveNRE relaxNRE)
CaseNames=(SineWaveStreaming Relaxation)
userOptions=("" "MICROPHYSICS=WEAKLIB")
gridLines=(84 129)

#grids=("[8,8,8]")
#gridNames=("")
#appNames=(ApplicationDriver)
#logFiles=(sineWaveJIRA)
#CaseNames=(SineWaveStreaming)
#userOptions=("")
#gridLines=(84)
#appNames=(ApplicationDriver_Neutrinos)
#logFiles=(relaxBugReduce)
#CaseNames=(Relaxation)
#userOptions=("MICROPHYSICS=WEAKLIB")
#gridLines=(129)

set_common

timeFOMLog="timeFOM_${COMPILER_DATE}.txt${umdf}$AADEBUG"
if [[ -z $ACTION || "$1" == "FOM" ]];then
   rm -rf $timeFOMLog
   echo "                                                        Time(seconds)                             |                      Figure of Merit (FOM)">>$timeFOMLog
   echo "AppName     Grid      OpLevel :  ${COMPILER_DATE}$umdf   ${BASE_CMP_DATE}${BASE_UMD}    TimeDiff   Percentage   |   ${COMPILER_DATE}${umdf}   ${BASE_CMP_DATE}${BASE_UMD}    FOM-Diff   Percentage">>$timeFOMLog
   echo "-----------------------------    --------------------------------------------------------------       --------------------------------------------------------------">>$timeFOMLog
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
         export LOG_BASE=${logFiles[jj]}.${OP_LEVEL}.${BASE_CMP_DATE}${BASE_UMD}${gridNames[ii]}$AADEBUG
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
            export exeName=${APP_NAME}_${THORNADO_MACHINE}${gridNames[ii]}

            if [[ "$1" == -[rR]* ]]; then
               LOG_FILE=${logFiles[jj]}.${OP_LEVEL}.${COMPILER_DATE}${umdf}${gridNames[ii]}${faction}$AADEBUG
               export LOG_FILE
               rm $LOG_FILE
               echo ""
               echo "Runing ${CaseNames[jj]} with a grid of ${grids[ii]} ..."
               echo ""
               if [ -f "${exeName}" ];then
                  echo "./${exeName}}"
                  runApp
               else
                  echo "The executable does not exist", ${exeName}
               fi
            elif [[ "$1" == -[bB]* ]]; then
               LOG_FILE=${logFiles[jj]}-BLD.${OP_LEVEL}.${COMPILER_DATE}${umdf}${gridNames[ii]}${faction}$AADEBUG
               export LOG_FILE
               rm $LOG_FILE
               rm ${exeName}
               echo ""
               echo "Compiling ${CaseNames[jj]} with a grid of ${grids[ii]} ..."
               echo ""
               buildApp
            else
               LOG_FILE=${logFiles[jj]}.${OP_LEVEL}.${COMPILER_DATE}${umdf}${gridNames[ii]}${faction}$AADEBUG
               export LOG_FILE
               rm $LOG_FILE
               rm ${exeName}
               echo ""
               echo "Compiling and Running ${CaseNames[jj]} with a grid of ${grids[ii]} ..."
               echo ""
               buildApp
               if [ -f "${exeName}" ];then
                  echo "./${exeName}"
                  runApp
               else
                  echo "The executable does not exist", ${exeName}
               fi
            fi   
         fi

         ## compare IMEX_TIME to the BASE_DATE
         if [[ -z $ACTION && "$1" != -[bB]* ]]; then
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
            echo "$caseName $gg   $OP_LEVEL    :   $currTime        $baseTime     $diffTime $percentage%          $currFOM        $baseFOM      $diffFOM $percentFOM%" >>$timeFOMLog
         fi
   done
done
done

if [[ -z $ACTION && "$1" != -[bB]* ]];then
   echo
   echo " Performance Comparison between compiler $COMPILER_DATE and $BASE_DATE"
   echo 

   echo "cat $timeFOMLog"
   cat $timeFOMLog
fi
