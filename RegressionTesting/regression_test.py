# Code based on regression_test.py written by Brandon Barker for
# https://github.com/AstroBarker/sitar/tree/main.

#!/usr/bin/env python3

import os
from subprocess import call, STDOUT
import shutil
import sys
import warnings

import glob
import numpy as np
import h5py as h5

# ============
# Utility
# ============


def soft_equiv(val, ref, tol=1.0e-5):
  """
  Soft equivalence of value and reference.
  """
  tiny = 1.0e-12
  numerator = np.fabs(val - ref)
  denominator = max(np.fabs(ref), tiny)

  if numerator / denominator > tol:
    return False
  else:
    return True


# End soft_equiv


class Regression:
  """
  regression testing class
  """

  def __init__(
    self='./',
    src_dir='../',
    build_dir="../SandBox/dgExperiments_Euler_Relativistic_IDEAL/Executables",
    executable="./ApplicationDriver_",
  ):
    self.src_dir = src_dir
    self.build_dir = os.path.relpath(build_dir)
    self.executable = executable
    self.data_dir = os.path.join(src_dir, "GoldData")

  # End __init__

  def __str__(self):
    print("========== REGRESSION TESTING ==========")
    print(f"Source Dir  : {self.src_dir}")
    print(f"Build Dir   : {self.build_dir}")
    print(f"Executable  : {self.executable}")
    print("========================================")
    return "\n"

  # End __str__

  def build_code(self):
    os.chdir(self.build_dir)

    # Base configure options
    configure_options = []

    make_call = []
    make_call.append("make")

    # Configure
    call(make_call)

  # End build_code

  def run_code(self):
    self.executable = glob.glob("./" + self.executable + '*')[0]
    if not os.path.isfile(self.executable):
      print(f"Executable '{self.executable}' does not exist!")
      sys.exit(os.EX_SOFTWARE)
    run_cmd = []  # empty now, can accomodate mpi runs
    outfile = open("out.dat", "w")
    call(self.executable, stdout=outfile, stderr=STDOUT)

  # End run_code

  def load_output(self, fn):
    """
    load simulation output
    """
    os.chdir("../Output")

    H5_H = h5.File( fn, 'r' )

    PF_D = H5_H[ 'Fluid Fields/Primitive/Comoving Baryon Density' ]

    return PF_D

  # End load_output

  def compare_gold(self, fn, fn_gold):
    """
    compare to gold data
    """
    # load sim data
    PF_D = self.load_output(fn)

    # load gold
    PF_D_Gold = self.load_output(fn_gold)

    # comparison
    PF_D_Status = np.array_equal(PF_D, PF_D_Gold)

    success = all([PF_D_Status])

    if success:
      print("Test Passed! :)")
      return os.EX_OK
    else:
      print("Some test failed! :(")
      print(f"Comoving Baryon Density : {PF_D_Status}")
      return os.EX_SOFTWARE


# End Regression


def main():
  """
  run regression tests for CaI, CaII, CaIII
  """
  reg = Regression()
  reg.build_code()
  reg.run_code()

  fn = "Advection_FluidFields_000101.h5"
  fn_gold = "Advection_FluidFields_000101.h5"

  return reg.compare_gold(fn, fn_gold)


if __name__ == "__main__":
  main()
