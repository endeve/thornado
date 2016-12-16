#!/bin/bash

export THORNADO_MACHINE=$1

if [[ $THORNADO_MACHINE == titan* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

  module unload fftw cray-hdf5 cray-petsc silo subversion
  module unload pgi gcc cce pathscale
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale PrgEnv-intel

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

elif [[ $THORNADO_MACHINE == sn1987b* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == sjdunham* ]]; then

echo
echo "INFO: Setting environment for" $THORNADO_MACHINE

fi


if [[ $THORNADO_MACHINE == titan_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5
  module load subversion

elif [[ $THORNADO_MACHINE == titan_cray ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5
  module load subversion

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

elif [[ $THORNADO_MACHINE == sn1987b ]]; then

  echo

elif [[ $THORNADO_MACHINE == sjdunham ]]; then

echo

else

  echo "  WARNING: Unknown machine " $THORNADO_MACHINE

fi
