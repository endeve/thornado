#!/bin/bash

module purge

#module load oneapi/eng-compiler/2022.01.30.005
#module load oneapi/eng-compiler/2022.01.30.006
#module load oneapi/eng-compiler/2022.01.30.007
module load nightly-compiler/2022.07.01

#module switch -f intel_compute_runtime/release/agama-prerelease-402 neo/agama-prerelease/491-22.22.023368.11404-main
module load iprof

export LTTNG_HOME=/localdisk/quanshao
mkdir -p $LTTNG_HOME

export LD_LIBRARY_PATH=/localdisk/quanshao/ExaStar/hdf57/lib64:$LD_LIBRARY_PATH

export THORNADO_DIR=/localdisk/quanshao/ExaStar/thornado
export WEAKLIB_DIR=/localdisk/quanshao/ExaStar/weaklib
export WEAKLIB_TABLES_DIR=/localdisk/quanshao/ExaStar/weaklib-tables
export THORNADO_MACHINE=beacon_intel
#export OMP_NUM_THREADS=1

export LIBOMPTARGET_PLUGIN=LEVEL0
export SYCL_DEVICE_FILTER=LEVEL_ZERO
##export LIBOMPTARGET_PLUGIN=OPENCL
export LIBOMPTARGET_DEBUG=0
unset EnableWalkerPartition
#export EnableWalkerPartition=1
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
export LIBOMPTARGET_LEVEL0_MEMORY_POOL=device,16,32
export IGC_OverrideOCLMaxParamSize=4096

VT_OUTPUT=vtune07June2022
rm -rf $VT_OUTPUT

rm output.log
module list |& tee output.log
#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-fds=yes ./ApplicationDriver_beacon_intel  |& tee -a output.val.log
#vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -r vtune-06032022 ./ApplicationDriver_beacon_intel 
source /sharedjf/mjh/tools/Intel_VTune_Profiler_2022.3.0_nda/env/vars.sh
#vtune -collect gpu-hotspots -knob characterization-mode=global-local-accesses -r vtune-06032022-02 ./ApplicationDriver_beacon_intel 
#(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -ring-buffer 10 -r $VT_OUTPUT ./ApplicationDriver_beacon_intel ) |& tee -a output.log
#(vtune -collect gpu-hotspots -knob target-gpu=0:154:0.0 -r $VT_OUTPUT ./ApplicationDriver_beacon_intel ) |& tee -a output.log

( time iprof ./ApplicationDriver_beacon_intel ) |& tee -a output.log
#( time  gdb ./ApplicationDriver_beacon_intel ) |& tee output.log
#( time  gdb-oneapi ./ApplicationDriver_beacon_intel ) |& tee output.log
