#!/usr/bin/env bash

File_t=t_AMReX_7lvl_512.dat
File_dt=dt_AMReX_7lvl_512.dat
echo "Time" > ${File_t}
echo "dt" > ${File_dt}
for i in Ad*.chk*; do
  awk 'NR==4{ print; }' ${i}/Header | cut -d' ' -f1-1 >> ${File_dt}
  awk 'NR==5{ print; }' ${i}/Header | cut -d' ' -f1-1 >> ${File_t}
done
