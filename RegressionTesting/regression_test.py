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

    PF_D  = H5_H[ 'Fluid Fields/Primitive/Comoving Baryon Density' ][()]
    PF_V1 = H5_H[ 'Fluid Fields/Primitive/Three-Velocity (1)'      ][()]
    PF_V2 = H5_H[ 'Fluid Fields/Primitive/Three-Velocity (2)'      ][()]
    PF_V3 = H5_H[ 'Fluid Fields/Primitive/Three-Velocity (3)'      ][()]
    PF_E  = H5_H[ 'Fluid Fields/Primitive/Internal Energy Density' ][()]

    return PF_D, PF_V1, PF_V2, PF_V3, PF_E

  # End load_output

  def compare_gold(self, fn, fn_gold):
    """
    compare to gold data
    """
    # load sim data
    PF_D, PF_V1, PF_V2, PF_V3, PF_E  = self.load_output(fn)

    # load gold
    PF_D_Gold, PF_V1_Gold, PF_V2_Gold, PF_V3_Gold, \
    PF_E_Gold  = self.load_output(fn_gold)

    # comparison
    PF_D_Status_Eq  = np.allclose(PF_D,  PF_D_Gold,  atol = 0.0, rtol = 1e-13)
    PF_V1_Status_Eq = np.allclose(PF_V1, PF_V1_Gold, atol = 0.0, rtol = 1e-13)
    PF_V2_Status_Eq = np.allclose(PF_V2, PF_V2_Gold, atol = 0.0, rtol = 1e-13)
    PF_V3_Status_Eq = np.allclose(PF_V3, PF_V3_Gold, atol = 0.0, rtol = 1e-13)
    PF_E_Status_Eq  = np.allclose(PF_E,  PF_E_Gold,  atol = 0.0, rtol = 1e-13)

    equal_success = all([PF_D_Status_Eq,  PF_V1_Status_Eq, \
                         PF_V2_Status_Eq, PF_V3_Status_Eq,
                         PF_E_Status_Eq])

    if equal_success:
      print("Soft Equality Tests Passed!")
      return os.EX_OK
    else:
      print("Failure!")
      print(f"Comoving Baryon Density : {PF_D_Status_Eq }")
      print(f"Three-Velocity (1)      : {PF_V1_Status_Eq}")
      print(f"Three-Velocity (2)      : {PF_V2_Status_Eq}")
      print(f"Three-Velocity (3)      : {PF_V3_Status_Eq}")
      print(f"Internal Energy Density : {PF_E_Status_Eq }")
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
