#!/usr/bin/env python3
import subprocess

nN = [2, 3]
nX = [8, 16, 32, 64, 128, 256, 512, 1024]

for nn in nN:
    for nx in nX:
        snn = str(nn).zfill(2)
        snx = str(nx).zfill(4)
        fem_degree = nn - 1

        fname = 'PoissonSolverTest3_Newtonian_nN{}_nX{}.inputs'.format(snn, snx)

        filename = \
('''##### PoissonSolverTest_Newtonian.inputs_singleLevel #####

# For all LOGICAL types, use 0 for .FALSE. and 1 for .TRUE.
# For all REAL types, use "e" and not "d", i.e. 1.0e3

debug.DEBUG = 1

# ###############

inputs.rho_c       = 1.50e2 # g/cc
inputs.SolarRadius = 7.00e5 # kilometers
inputs.R_star      = 1.0 # solar radii
inputs.r_c         = 0.2 # solar radii

##################

thornado.ProgramName = "PoissonSolverTest_Newtonian"

thornado.iRestart = -1

thornado.UsePhysicalUnits = 1

thornado.t_end   = 1.5e+2 # For this test t_end is a dummy

thornado.PlotFileNameRoot = "PoissonSolverTest3_Newtonian_nN{snn}_nX{snx}.plt"
thornado.iCycleW = 1

thornado.CheckpointFileNameRoot = "PoissonSolverTest3_Newtonian_nN{snn}_nX{snx}.chk"
thornado.iCycleChk = 1

thornado.nNodes = {nn}

thornado.bcX         = 30 00 00 # 30 is 3 is inner and 0 is outer and 3 means reflecting 0 is leave it fixed
geometry.is_periodic = 0  1  1

geometry.coord_sys = 2
geometry.prob_lo   = 000.0 0.0               0.0
geometry.prob_hi   = 14.0e5 3.1415926535897931 6.2831853071795862

amr.n_cell                  = {nx} 01 01
thornado.swX                = 01  00 00
#amr.max_grid_size_x         = 512
#amr.blocking_factor_x       = 512
amr.max_level               = 0
amr.UseFluxCorrection_Euler = 0
amr.UseAMR                  = 0
amr.TagCriteria             = 0.0
amr.n_error_buf             = 0
amr.ref_ratio               = 2

# Poseidon parameters
GS.EvolveGravity    = 1
poseidon.fem_degree = {fem_degree}
poseidon.l_limit    = 0

# Slope limiter parameters
SL.UseSlopeLimiter_Euler           =  0 # was 1
SL.SlopeLimiterMethod_Euler        = "TVD"
SL.BetaTVD_Euler                   = 1.75e+0
SL.BetaTVB_Euler                   = 0.00e+0
SL.SlopeTolerance_Euler            = 1.00e-6
SL.UseCharacteristicLimiting_Euler = 0
SL.UseTroubledCellIndicator_Euler  = 0
SL.LimiterThresholdParameter_Euler = 5.00e-3
SL.UseConservativeCorrection_Euler = 1

# Positivity limiter parameters
PL.UsePositivityLimiter_Euler = 0
PL.Min_1_Euler                = 1.0e-13
PL.Min_2_Euler                = 1.0e-13

# Equation of state parameters
EoS.EquationOfState = "IDEAL"
EoS.Gamma_IDEAL     = 1.30
''').format(snn=snn, snx=snx, nn=nn, nx=nx, fem_degree=fem_degree)

        with open(fname, 'w') as fl:
            fl.write(filename)
print('Input files generated sucessfully')

print('*************')
# Ask if the user wants to run the simulation
run_sims = input("Do you want to run the simulations (yes/no or y/n)? : ").strip().lower()
if run_sims in ['yes', 'y']:
    print("Running run_sims.bash ...")
    subprocess.run(["bash", "run_sims.bash"], check=True)
    print("run_sims.bash has been executed successfully.")
    #print("The nodal data is stored as ")
else:
    print("Skipping running the sims ")
