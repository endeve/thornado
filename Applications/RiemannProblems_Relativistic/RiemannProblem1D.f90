PROGRAM RiemannProblem1D

  USE KindModule, ONLY: &
    DP
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
             = 'RiemannProblem1D_Relativistic', &
           nX_Option &
             = [ 2000, 1, 1 ], &
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
           EquationOfState_Option &
             = 'IDEAL', &
           Gamma_IDEAL_Option &
             = 5.0_DP / 3.0_DP, &
           FluidSolver_Option &
             = 'Euler_DG_GR', &
           FluidRiemannSolver_Option &
             = 'LLF', &
           EvolveFluid_Option &
             = .TRUE., &
           ApplySlopeLimiter_Option &
             = .FALSE., &
           BetaTVB_Option &
             = 0.0_DP, &
           BetaTVD_Option &
             = 1.0_DP, &
           ApplyPositivityLimiter_Option &
             = .FALSE., &
           nStages_SSP_RK_Option &
             = 1 )

  CALL InitializeRiemannProblem1D &
         ( D_L = 1.0_DP, V_L = [ 0.0_DP, 0.0_DP, 0.0_DP], P_L = 1.0d+3, &
           D_R = 1.0_DP, V_R = [ 0.0_DP, 0.0_DP, 0.0_DP], P_R = 1.0d-2, &
           X_D_Option = 0.5_DP )

  CALL EvolveFields &
         ( t_begin  = 0.0_DP, &
           t_end    = 3.5d-1, &
           dt_write = 3.5d-3, &
           UpdateFields = SSP_RK, &
           dt_fixed_Option = 1.0d-4 )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D
