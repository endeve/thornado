#!/usr/bin/env bash

for nn in 2 3; do
  for nx in 8 16 32 64 128 256 512 1024; do

     formatted_nn=$(printf "%02d" $nn)
     formatted_nx=$(printf "%04d" $nx)

     input_file="PoissonSolverTest3_Newtonian_nN${formatted_nn}_nX${formatted_nx}.inputs"
     output_file="PoissonSolverTest3_Newtonian_nN${formatted_nn}_nX${formatted_nx}.out"

     mpiexec -n 1 .././main1d.gnu.DEBUG.MPI.ex $input_file > $output_file &
  done
done
