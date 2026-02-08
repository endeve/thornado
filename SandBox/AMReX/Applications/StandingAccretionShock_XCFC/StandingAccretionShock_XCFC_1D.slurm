#!/bin/bash

#SBATCH --mail-user=samueljdunham+astro@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --exclude=cn1609
#SBATCH --time=0-00:05:00
#SBATCH --mem=128G
#SBATCH --output=StandingAccretionShock_XCFC_1D.%j.out
#SBATCH --job-name=SAS1D

# Adapted from: https://github.com/AMReX-Astro/workflow/blob/main/
#               job_scripts/cori-haswell/cori_haswell.MPI.slurm

ROOT=StandingAccretionShock_XCFC_1D
INFILE=${ROOT}.inputs
OUTFILE=${ROOT}.${SLURM_JOB_ID}.out
SLURMFILE=${ROOT}.slurm
CHKFILEROOT=${ROOT}.chk
EXEC=./main1d.gnu.MPI.ex

echo "##### StandingAccretionShock_XCFC_1D.slurm #####" >> ${OUTFILE}
echo "" >> ${OUTFILE}
cat ${SLURMFILE} >> ${OUTFILE}

echo "" >> ${OUTFILE}
cat GNUmakefile >> ${OUTFILE}

echo "" >> ${OUTFILE}
cat Make.package >> ${OUTFILE}

echo "" >> ${OUTFILE}
cat ${INFILE} >> ${OUTFILE}

echo "" >> ${OUTFILE}
cat ../TaggingModule.f90 >> ${OUTFILE}
echo "" >> ${OUTFILE}

function findChkFile {

  # find_chk_file takes a single argument -- the wildcard pattern
  # for checkpoint files to look through
  chk=$1

  # find the latest 2 restart files.  This way if the latest didn't
  # complete we fall back to the previous one.
  tempFiles=$(find . -maxdepth 1 -name "${chk}" -print | sort | tail -2)
  restartFile=""
  for f in ${tempFiles}
  do
    echo ${f}
    # the Header is the last thing written -- if it's there,
    # update the restart file
    if [ -f ${f}/Header ]; then
      restartFile="${f: -8}"
    fi
  done
}

# look for 8-digit chk files
findChkFile "${CHKFILEROOT}????????"

# restartString will be empty if no chk files are found -- i.e. new run
if [ "${restartFile}" = "" ]; then
  restartString=""
else
  restartString="thornado.iRestart=${restartFile}"
fi

module restore thornado
srun --mpi=pmi2 ${EXEC} ${INFILE} ${restartString}

