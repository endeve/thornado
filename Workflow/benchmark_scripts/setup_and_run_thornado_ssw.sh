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
ROCM_VERSION="${ROCM_VERSION:-6.4.2}"
PROJ="${PROJ:-stf006}"
PRGENV="${PRGENV:-cray}"
HOST_SHORT="${LMOD_SYSTEM_NAME:-unknownhost}"

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
THORNADO_SHA="61764141124259a1d394ff29b66198d5d8ee71de" # HEAD of olcf_test_suite on 01/16/2026

# ---- Setup the build/run environment ----

## modules
module load cpe/${CPE_VERSION}
module load PrgEnv-$PRGENV
module load cray-hdf5-parallel
module load craype-accel-amd-gfx90a
module load rocm/${ROCM_VERSION}

## can load my hipfort if system-built unavailable
if ! module avail hipfort/${ROCM_VERSION} &>/dev/null; then
  module use /ccs/home/jaharris/modulefiles/${HOST_SHORT}
fi
module load hipfort/${ROCM_VERSION}

## Need this if not using default ROCm version
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

COMPILER="${LMOD_FAMILY_COMPILER}" # cce, amd, afar-amd, etc

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
echo "Running: ${SRUN[*]} ./${EXEC}"
"${SRUN[@]}" "./${EXEC}" > >(tee "${RUN_OUTPUT_FILE}") 2>&1

cd -

echo "Checking results"
python3 thornado_check_ssw.py "$RUNDIR/$RUN_OUTPUT_FILE"
echo "exit=$?"

