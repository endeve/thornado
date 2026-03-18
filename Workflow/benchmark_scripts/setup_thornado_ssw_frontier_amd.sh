#!/usr/bin/env bash

#############################################
# Thornado streaming sine wave benchmark setup script
#
# Usage examples:
#   INSTALL_ROOT=/lustre/orion/stf006/scratch/$USER PROJ=stf006 TIME=00:10:00 ./setup_thornado_ssw.sh
#############################################

# ---- helpers ----
die() { echo "ERROR: $*" >&2; exit 1; }

# ---- User-tweakable knobs (can be overridden via env) ----
CPE_VERSION="${CPE_VERSION:-25.09}"
ROCM_VERSION="${ROCM_VERSION:-7.12.0}"
PROJ="${PROJ:-stf006}"
PRGENV="${PRGENV:-amd}"
HOST_SHORT="${LMOD_SYSTEM_NAME:-unknownhost}"
AFAR_VERSION="${AFAR_VERSION:-23.1.0-7.12.0}"

# Where to put repos + large runtime tables
## (weaklib-tables is >20 GB and read at run-time, so don't use $HOME)
INSTALL_ROOT="${INSTALL_ROOT:-$PWD/thornado_ssw_test}"
[[ -n "$INSTALL_ROOT" ]] || die "INSTALL_ROOT is empty. Set INSTALL_ROOT to a writable directory"
mkdir -p "$INSTALL_ROOT" || die "Cannot create INSTALL_ROOT=$INSTALL_ROOT"

# Slurm run settings
PARTITION="batch"
QOS="debug"
TIME="00:10:00"
NODES="1"
NTASKS="1"
CPUS_PER_TASK="1"
NTASKS_PER_GPU="1"

## Git coordinates
THORNADO_REPO="git@github.com:endeve/thornado.git"
THORNADO_SHA="dbb167440333903fce1c84b711f6b687825fc3cb" # HEAD of olcf_test_suite on 03/18/2026

# ---- Setup the build/run environment ----

## modules
module reset
source /opt/cray/pe/cpe/${CPE_VERSION}/restore_lmod_system_defaults.sh
module use /sw/crusher/ums/compilers/modulefiles/afar-prgenv-modules/modulefiles
module load PrgEnv-${PRGENV}
module load cray-mpich
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load afar-prgenv/${AFAR_VERSION}
export OLCF_HDF5_ROOT=/lustre/orion/world-shared/stf006/jaharris/sw/frontier/amd/${AFAR_DROP_VERSION}/hdf5-parallel-1.14.3

## can load my hipfort if system-built unavailable
## NOTE: hipfort module not needed with AFAR
#if ! module avail hipfort/${ROCM_VERSION} &>/dev/null; then
#  module use /ccs/home/jaharris/modulefiles/${HOST_SHORT}
#fi
#module load hipfort/${ROCM_VERSION}
if [[ -z ${OLCF_HIPFORT_ROOT+x} ]]; then
  if [[ -d "${ROCM_PATH}/include/hipfort" ]]; then
    export OLCF_HIPFORT_ROOT="${ROCM_PATH}"
  elif [[ -d "${ROCM_PATH}/llvm/include/hipfort" ]]; then
    export OLCF_HIPFORT_ROOT="${ROCM_PATH}/llvm"
  else
    die "Cannot find hipfort installation"
  fi
fi

## Need this if not using default ROCm version
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

COMPILER="amd"

export THORNADO_MACHINE="${HOST_SHORT}_${COMPILER}"

THORNADO_APP="ApplicationDriver"
RUN_OUTPUT_FILE="ssw.out"

export THORNADO_DIR="${INSTALL_ROOT}/thornado"

## Set the build/run directory
BUILDDIR="${THORNADO_DIR}/SandBox/TwoMoment_OrderV/Executables"
RUNDIR="${THORNADO_DIR}/SandBox/TwoMoment_OrderV/Executables"

# ---- Get the code (at a specific SHA)  ----

git_clone_checkout() {
  local repo="$1" repo_dir="$2" sha="$3"
  if [[ ! -d "$repo_dir" ]]; then
    git clone -v "$repo" "$repo_dir"
  fi
  git -C "$repo_dir" fetch --all --tags --prune
  git -C "$repo_dir" checkout --detach "$sha"
}

git_clone_checkout "$THORNADO_REPO" "$THORNADO_DIR" "$THORNADO_SHA"

# ---- Setup runtime tables + symlinks ----

mkdir -p "${RUNDIR}/../Output"
 
# ---- Build ----

make -C ${BUILDDIR} -j 16 ${THORNADO_APP} \
  OPT_LEVEL=OPTIMIZE \
  USE_GPU=TRUE USE_OMP_OL=TRUE USE_ROCM=TRUE USE_HIP=TRUE

EXEC="${THORNADO_APP}_${THORNADO_MACHINE}"
[[ -x "${BUILDDIR}/${EXEC}" ]] || die "Build succeeded but executable not found: ${BUILDDIR}/${EXEC}"
cp -v "${BUILDDIR}/${EXEC}" "${RUNDIR}/."
 
# ---- Run ----

cd "$RUNDIR"

SRUNCMD=(
  srun -u
  -A "$PROJ"
  -p "$PARTITION"
  -q "$QOS"
  -t "$TIME"
  -N "$NODES"
  -n "$NTASKS"
  -c "$CPUS_PER_TASK"
  --ntasks-per-gpu="$NTASKS_PER_GPU"
  --gpu-bind=closest
)
echo "Running: ${SRUNCMD[*]} ./${EXEC}"
"${SRUNCMD[@]}" "./${EXEC}" > >(tee "${RUN_OUTPUT_FILE}") 2>&1

cd -

echo "Checking results"
python3 thornado_check_ssw.py "$RUNDIR/$RUN_OUTPUT_FILE"
echo "exit=$?"

