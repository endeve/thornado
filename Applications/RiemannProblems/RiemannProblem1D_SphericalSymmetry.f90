PROGRAM RiemannProblem1D_SphericalSymmetry

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeRiemannProblem1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
            = 'RiemannProblem1D_SphericalSymmetry', &
           nX_Option &
             = [ 1000, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 3, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 2.0_DP, Pi, 2 * Pi ], &
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
             = 'HLLC', &
           EvolveFluid_Option &
             = .TRUE., &
           ApplySlopeLimiter_Option &
             = .TRUE., &
           BetaTVB_Option &
             = 0.0d0, &
           BetaTVD_Option &
             = 2.0d0, &
           ApplyPositivityLimiter_Option &
             = .FALSE., &
           nStages_SSP_RK_Option &
             = 1 )

  CALL InitializeRiemannProblem1D &
         ( D_L = 1.000_DP, V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ], P_L = 1.0_DP, &
           D_R = 0.125_DP, V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ], P_R = 0.1_DP, &
           X_D_Option = 1.0_DP )

  CALL EvolveFields &
         ( t_begin  = 0.0d-0,  &
           t_end    = 5.0d-1, &
           dt_write = 2.5d-2, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D_SphericalSymmetry
