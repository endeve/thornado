#!/bin/bash

  ml cuda nvhpc essl hdf5 magma netlib-lapack

  export APP_NAME=ApplicationDriver_Neutrinos
  export THORNADO_MACHINE=summit_nvhpc


   rm ./${APP_NAME}_${THORNADO_MACHINE}
   module list   
   make clean

( time make -j8 $APP_NAME MICROPHYSICS=WEAKLIB USE_OACC=TRUE  USE_GPU=TRUE USE_CUDA=TRUE USE_CUBLAS=TRUE )


   export PGI_ACC_DEBUG=0
   export LD_LIBRARY_PATH=${HDF5_LIB}:$LD_LIBRARY_PATH

 ( time mpirun -n 1 --bind-to core ./${APP_NAME}_${THORNADO_MACHINE} )
