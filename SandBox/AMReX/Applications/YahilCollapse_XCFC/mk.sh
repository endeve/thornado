#!/usr/bin/env bash

rm -rf amrex* *plt* *chk* *.dat* *.old.* && \
make DEBUG=FALSE -j12 && \
mpiexec -quiet -n 1 ./main1d.gnu.MPI.ex YahilCollapse_XCFC.inputs_multiLevel
rm -rf *.old.* Yahi*dat
