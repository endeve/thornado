# Thornado oneMKL -- Summary
This is created based oneMKL branch of https://github.com/endeve/thornado.git   
commit: `b19c9a4e1ac4204b4e5cd4526beddd851605f021`  
Application used: SandBox/TwoMoment_OrderV   
The code was first cloned on Feb 22 2022, and the suggested extra compiling flags are:  
  `USE_OMP_OL=TRUE USE_GPU=TRUE USE_ONEMKL=TRUE`
The code is located at `$HOME/ExaStar`  and we mainly work on `$HOME/ExaStar/thornado/SandBox/TwoMoment_OrderV`  
   (Previous working code: https://gitlab.devtools.intel.com/A21-NRE-EPA/NRE-Codes/thornado)  
   

# Thornado oneMKL -- How to Build and Run
1. Two extra external packages are needed   
   a)  HDF5   
   b)  weaklib-table from ORNL. source: git clone https://code.ornl.gov/astro/weaklib-tables.git   
More information on the external packages, please visit: https://gitlab.devtools.intel.com/williamr/thornado-builder  
2. Enviromental variable for weaklib-table needs to be set: export WEAKLIB_TABLES_DIR=$HOME/ExaStar/weaklib-tables
3. Run `softlink.sh` inside `$HOME/ExaStar/thornado/SandBox/TwoMoment_OrderV/Executables` to create links to the external libraries. Be sure to edit the content of this script as your path to the libraries may be different.
4. Copy `/nfs/site/home/quanshao/ExaStar/mpifort` to your `$HOME/ExaStar` and modify `$HOME/ExaStar/thornado/Build/Machines/Makefile_beacon_intel` to reflect the correct path of `mpifort`
5. To compile,  run `build.anl.sh` inside `$HOME/ExaStar/thornado/SandBox/TwoMoment_OrderV/Executables`. Be sure to make modifications to reflect the path of your source code and change the module to your favorite one
6. To run apps, run `run.anl.active.sh` inside `$HOME/ExaStar/thornado/SandBox/TwoMoment_OrderV/Executables`. Be sure to make modifications to reflect the path of your source code and change the module to your favorite one
