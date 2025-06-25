#!/bin/bash

#set -x

###########################################################################################
###  Common Modules and Enviroment Variables.
###########################################################################################

function set_common(){


   export EXASTAR_HOME=/localdisk/quanshao/ExaStar
   export HDF5_INC=${EXASTAR_HOME}/hdf57/include
   export HDF5_LIB=${EXASTAR_HOME}/hdf57/lib64
   export THORNADO_DIR=${EXASTAR_HOME}/thornado-anl
   export WEAKLIB_DIR=${EXASTAR_HOME}/weaklib
   export WEAKLIB_TABLES_DIR=${EXASTAR_HOME}/weaklib-tables
   export THORNADO_MACHINE=sunspot_intel

   export IGC_OverrideOCLMaxParamSize=4096
   #export numTeamCap="-mllvm -vpo-paropt-atomic-free-red-global-buf-size=4096"
   export numTeamCap="-mllvm -vpo-paropt-atomic-free-reduction-slm=true" 

## for running

   export FI_PROVIDER=sockets
   export ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE
   export MPIR_CVAR_ENABLE_GPU=0
   export ZE_AFFINITY_MASK=0.0
   export LTTNG_HOME=$EXASTAR_HOME
   mkdir -p $LTTNG_HOME
   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH
   export OMP_TARGET_OFFLOAD=MANDATORY
   export OMP_NUM_THREADS=1
   ulimit -s unlimited
   export LIBOMPTARGET_LEVEL_ZERO_MEMORY_POOL=device,128,64,16384
#   export PrintVerboseGenericControlFlowLog=1
#   export IGC_StackOverflowDetection=1

#   export IGC_ShaderDumpEnable=1
#   export IGC_ShowFullVectorsInShaderDumps=1
   #export LIBOMPTARGET_LEVEL_ZERO_COMMAND_BATCH=copy,8
   #export IGC_EnableZEBinary=0
   #export IGC_ForceOCLSIMDWidth=16
   #export LIBOMPTARGET_LEVEL_ZERO_USE_IMMEDIATE_COMMAND_LIST=0
   #export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=0
   #export MKL_VERBOSE=2
   #export LIBOMPTARGET_PLUGIN=LEVEL0
   #export ONEAPI_DEVICE_FILTER=level_zero:gpu
   ##export LIBOMPTARGET_PLUGIN=OPENCL
#   export LIBOMPTARGET_DEBUG=1
#   export LIBOMPTARGET_INFO=4
   #export EnableImplicitScaling=1
   #export LIBOMPTARGET_PLUGIN_PROFILE=T
   #ulimit -n 20480
}

###########################################################################################
###  Make Script
###########################################################################################

function buildApp(){

   echo $MKLROOT |& tee -a $LOG_FILE 
   module list   |& tee -a $LOG_FILE

   make clean
   ( time make $APP_NAME ${USER_OPTION} USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a $LOG_FILE
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

   if [[ "$ACTION" == "iprof" ]]; then
      module load iprof/0.11.2
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

#export A21_SDK_MKLROOT_OVERRIDE=/exaperf/nightly/mkl-cev/2022.11.02 ## Latest nightly, i.e. 2022-10-06, uses this mkl

#ACTION="iprof"
#ACTION="perf"
#ACTION="onetrace"
#ACTION="advisor"
#ACTION="vtune"
#ACTION=""

if [[ -n $ACTION ]];then
   faction="-$ACTION"
fi
#BASE_DATE="2023.10.15.002"
BASE_DATE="2025.04.27"
BASE_MKL="mkl2025.04.24"
BASE_UMD="umd1130"
#MKL_BASE_DATE="" ## A underline is need before the date string for clarity
#export AADEBUG="-g"
export useAGRF="TRUE"


#export A21_SDK_MKLROOT_OVERRIDE=/exaperf/nightly/mkl-cev_nightly/2024.06.02
#COMPILER_DATE=2024.06.29
if [[ "$1" == "22" ]];then
    COMPILER_DATE=2024.04.22
elif [[ "$1" == "23" ]];then
    COMPILER_DATE=2024.04.23
fi    
#module load nightly-compiler/${COMPILER_DATE}
MKL_DATE=2025.05.29
#COMPILER_DATE=2024.12.10
ml use /opt/software/oneapi-lkg/modulefiles/oneapi-lkg
module load 2025.2.0
module load neo/agama-hotfix-devel-1099-sp4/12-25.05.32567.18-1099
module load mpich/52.2-256/icc-sockets-gpu
ml list
#module load nightly-mkl-cev_nightly/$MKL_DATE
#module load  nightly-mkl-develop/$MKL_DATE
#module load nightly-mkl-cev_rls/$MKL_DATE
#COMPILER_DATE=2025.02.06
#LD_LIBRARY_PATH=/exaperf/nightly/compiler/2024.03.06/linux/lib/x86_64-unknown-linux-gnu:$LD_LIBRARY_PATH
#
#module switch -f nightly-compiler/2025.05.07 nightly-compiler/2025.05.06

#if action is empty, performance comparison will be done. otherwise there is no performance comparison and just run the app using such as onetrace, vtune etc. so action can be "", "onetrace", "iprof", "vtune", 

#opLevels=(O3)
#grids=("[8,8,8]")
#grids=("[16,16,16]")
#gridNames=("")
#gridNames=("-xN16")
#appNames=(ApplicationDriver)
#logFiles=(sineWaveShaderDump)
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
#gridLines=(130)

#opLevels=(O0 O1 O2 O3)
opLevels=(O3)
grids=("[8,8,8]" "[16,16,16]")
gridNames=("" "-xN16")
appNames=(ApplicationDriver ApplicationDriver_Neutrinos)
logFiles=(SineWave Relax)
CaseNames=(SineWaveStreaming Relaxation)
userOptions=("" "MICROPHYSICS=WEAKLIB")
gridLines=(84 202)

#grids=("[16,16,16]")
#gridNames=("-xN16")
#grids=("[8,8,8]")
#gridNames=("")
#appNames=(ApplicationDriver)
#logFiles=(sineWave)
#CaseNames=(SineWaveStreaming)
#userOptions=("")
#gridLines=(84)
#appNames=(ApplicationDriver_Neutrinos)
#logFiles=(relaxHang)
#CaseNames=(Relaxation)
#userOptions=("MICROPHYSICS=WEAKLIB")
#gridLines=(202)

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

timeFOMLog="TimeFOM-${COMPILER_DATE}-mkl${MKL_DATE}.txt${umdf}$AADEBUG"

if [[ -z $ACTION || "$1" == "FOM" ]];then
   rm -rf $timeFOMLog
   echo "                                                        Time(seconds)                              |                      Figure of Merit (FOM)">>$timeFOMLog
   echo "AppName     Grid      OpLevel :  ${COMPILER_DATE}$umdf   ${BASE_DATE}${BASE_UMD}    TimeDiff   Percentage   |   ${COMPILER_DATE}${umdf}   ${BASE_DATE}${BASE_UMD}    FOM-Diff   Percentage">>$timeFOMLog
   echo "                     MKL Date :  $MKL_DATE               $MKL_BASE_DATE">>$timeFOMLog

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
         export LOG_FILE=${logFiles[jj]}.${OP_LEVEL}.${COMPILER_DATE}-mkl${MKL_DATE}${umdf}${gridNames[ii]}${faction}$AADEBUG
         export LOG_BASE=${logFiles[jj]}.${OP_LEVEL}.${BASE_DATE}-${BASE_MKL}-${BASE_UMD}${gridNames[ii]}$AADEBUG
	 echo $LOG_FILE    $LOG_BASE
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
            echo "$caseName $gg   $OP_LEVEL    :   $currTime        $baseTime     $diffTime $percentage%          $currFOM        $baseFOM      $diffFOM $percentFOM%" >>$timeFOMLog
         fi
   done
done
done
echo "  " |& tee -a $timeFOMLog
ml list |& tee -a $timeFOMLog
tail -n 2 $timeFOMLog |wc -c|xargs -I {} truncate $timeFOMLog -s -{}
ifx -what bb.f90 |& tee -a $timeFOMLog
rm a.out
if [[ -z $ACTION ]];then
   echo
   echo "Performance Comparison between compiler $COMPILER_DATE and $BASE_DATE"
   echo 

   echo "cat $timeFOMLog"
   cat $timeFOMLog
fi
