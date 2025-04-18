#!/bin/bash

export THORNADO_MACHINE=$1

if [[ $THORNADO_MACHINE == mac* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == titan* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

  module unload fftw cray-hdf5 cray-petsc silo subversion
  module unload pgi gcc cce pathscale
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale PrgEnv-intel

elif [[ $THORNADO_MACHINE == summit* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

  module unload xl spectrum-mpi hsi xalt lsf-tools darshan-runtime
  module unload DefApps

elif [[ $THORNADO_MACHINE == ascent* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

  module unload xl spectrum-mpi hsi xalt lsf-tools darshan-runtime
  module unload DefApps

elif [[ $THORNADO_MACHINE == frontier* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

elif [[ $THORNADO_MACHINE == perlmutter* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

elif [[ $THORNADO_MACHINE == darter* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload fftw cray-hdf5 cray-petsc
  module unload pgi gcc cce intel
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-intel

elif [[ $THORNADO_MACHINE == beacon* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload pgi gcc cce intel
  module unload PE-gnu PE-intel

elif [[ $THORNADO_MACHINE == acf* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload PE-intel PE-gnu

elif [[ $THORNADO_MACHINE == acf* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module purge
  module load GNU/8.2.0-2.31.1 OpenMPI/3.1.3 HDF5/1.10.4

elif [[ $THORNADO_MACHINE == isaac* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module unload PE-intel
  module unload intel-mpi
  module unload intel-compilers

elif [[ $THORNADO_MACHINE == sn1987b* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == bbarker* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == juliana* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == kristopher* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == sjdunham* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == kkadoogan* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == jrober* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE
  export POSEIDON_DIR=~/Poseidon

elif [[ $THORNADO_MACHINE == rmurph* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == accre* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  module purge

elif [[ $THORNADO_MACHINE == ranchu* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == zelledge* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == ranchuair* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == mcarpe21* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == dpochik* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == lucalin* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == paullaiu_gnu* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

fi

if [[ $THORNADO_MACHINE == mac_gnu ]]; then

  echo
  export WEAKLIB_DIR=
  export POSEIDON_DIR=

elif [[ $THORNADO_MACHINE == titan_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5
  module load subversion

elif [[ $THORNADO_MACHINE == titan_cray ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5
  module load subversion

elif [[ $THORNADO_MACHINE == summit_pgi ]]; then

  echo

  module load pgi/20.4
  module load spectrum-mpi
  module load hdf5
  module load netlib-lapack
  module load essl

elif [[ $THORNADO_MACHINE == summit_nvhpc ]]; then

  echo

  module load nvhpc/23.9
  module load spectrum-mpi
  module load hdf5
  module load netlib-lapack
  module load essl

elif [[ $THORNADO_MACHINE == summit_xl ]]; then

  echo

  module load xl
  module load spectrum-mpi
  module load hdf5/1.10.7
  module load netlib-lapack/3.8.0
  module load essl

elif [[ $THORNADO_MACHINE == summit_gcc ]]; then

  echo

  module load gcc/12.1.0
  module load spectrum-mpi
  module load hdf5
  module load netlib-lapack
  module load essl
  module load cuda

elif [[ $THORNADO_MACHINE == ascent_pgi ]]; then

  echo

  module load pgi
  module load spectrum-mpi
  module load hdf5
  module load netlib-lapack
  module load essl
  module load cuda

elif [[ $THORNADO_MACHINE == frontier_cce ]]; then

  echo

  module restore

  module use /ccs/home/jaharris/modulefiles/frontier

  module load cpe/24.11
  module load PrgEnv-cray cray-hdf5-parallel craype-accel-amd-gfx90a 
  module load rocm/6.2.4 hipfort/6.2.4 
  module unload darshan-runtime

  export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

  module -t list
  
elif [[ $THORNADO_MACHINE == frontier_gcc ]]; then

  echo
  
  module load cpe/23.09
  module load PrgEnv-gnu
  module load cray-hdf5-parallel
  module unload darshan-runtime

  module -t list
  
elif [[ $THORNADO_MACHINE == perlmutter_nvhpc ]]; then

  echo

  module load PrgEnv-nvidia
  module load cray-hdf5-parallel
  module unload cpu
  module load gpu
  
elif [[ $THORNADO_MACHINE == perlmutter_cce ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5-parallel
  module unload cpu
  module load gpu

elif [[ $THORNADO_MACHINE == perlmutter_gcc ]]; then

  echo

  module load PrgEnv-gnu
  module load cudatoolkit
  module load cray-hdf5-parallel
  module unload cpu
  module load gpu

elif [[ $THORNADO_MACHINE == darter_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5

elif [[ $THORNADO_MACHINE == darter_cray ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5

elif [[ $THORNADO_MACHINE == beacon_intel ]]; then

  echo

  module load PE-intel
  module load hdf5/1.8.14

elif [[ $THORNADO_MACHINE == acf_gnu ]]; then

  echo

  module load PE-gnu
  module swap intel-mpi/2018.1.163 openmpi/3.0.0-gcc6.3.0
  module load hdf5
  module load lapack

elif [[ $THORNADO_MACHINE == isaac_intel ]]; then

  echo

  module load intel-compilers
  module load intel-mpi
  module load PE-intel
  module load hdf5/1.10.8-intel
  module load netlib-lapack/3.8.0

  export LAPACK_DIR=/spack/spack-0.16.3/apps/linux-rhel8-cascadelake/gcc-10.2.0/netlib-lapack-3.8.0-xh7uzdm4yijch4yghw3fo2tvv3wrq553

elif [[ $THORNADO_MACHINE == sn1987b ]]; then

  echo

elif [[ $THORNADO_MACHINE == sn1987b ]]; then

  echo

elif [[ $THORNADO_MACHINE == bbarker ]]; then

  echo

  module purge
  module load GNU/8.2.0-2.31.1
  module load OpenMPI/3.1.3
  module load HDF5/1.10.4

elif [[ $THORNADO_MACHINE == juliana ]]; then

  echo
  export WEAKLIB_DIR=~/Desktop/thornado_files/weaklib
  export THORNADO_DIR=~/Desktop/thornado_files/thornado_new

elif [[ $THORNADO_MACHINE == kristopher ]]; then

  echo
  export THORNADO_DIR=~/Desktop/Work/ORNL/thornado
  export POSEIDON_DIR=~/Desktop/Work/ORNL/poseidon

elif [[ $THORNADO_MACHINE == sjdunham ]]; then

  echo

elif [[ $THORNADO_MACHINE == kkadoogan_gnu ]]; then

  echo

elif [[ $THORNADO_MACHINE == jrober ]]; then

  echo

elif [[ $THORNADO_MACHINE == rmurph ]]; then

  echo

elif [[ $THORNADO_MACHINE == accre_gnu ]]; then

  echo

  module load GCC/10.2.0
  module load OpenMPI/4.0.5
  module load ScaLAPACK/2.1.0
  module load HDF5/1.10.7

  export LAPACK_DIR=$EBROOTSCALAPACK

elif [[ $THORNADO_MACHINE == ranchu ]]; then

  echo

elif [[ $THORNADO_MACHINE == zelledge ]]; then

  echo
  export WEAKLIB_DIR=/mnt/c/Users/Zack/weaklib/
  export THORNADO_DIR=/mnt/c/Users/Zack/thornado/

elif [[ $THORNADO_MACHINE == ranchuair ]]; then

  echo

elif [[ $THORNADO_MACHINE == mcarpe21 ]]; then

  echo
  export THORNADO_DIR=~/thornado
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/hdf5/lib/

elif [[ $THORNADO_MACHINE == dpochik ]]; then

  echo
  export WEAKLIB_DIR=/home/dpochik/weaklib
  export THORNADO_DIR=/home/dpochik/thornado_neos

elif [[ $THORNADO_MACHINE == lucalin ]]; then

  echo
  export WEAKLIB_DIR=/home/luca/EOS_tables/weaklib
  export THORNADO_DIR=/home/luca/thornado

elif [[ $THORNADO_MACHINE == paullaiu_gnu ]]; then

  echo
  export THORNADO_DIR=/home/plc/Documents/Git/thornado
  export WEAKLIB_DIR=/home/plc/Documents/Git/Neutrino/weaklib

else

  echo "  WARNING: Unknown machine " $THORNADO_MACHINE

fi
