PROGRAM RiemannProblem1D_SphericalSymmetry

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
            = 'RiemannProblem1D_SphericalSymmetry', &
           nX_Option &
             = [ 1024, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = 1, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           EquationOfState_Option &
             = 'IDEAL', &
           Gamma_IDEAL_Option &
             = 1.4_DP, &
           FluidSolver_Option &
             = 'Euler_DG', &
           FluidRiemannSolver_Option &
             = 'HLL', &
           EvolveFluid_Option &
             = .TRUE., &
           nStages_SSP_RK_Option &
             = 1 )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D_SphericalSymmetry
