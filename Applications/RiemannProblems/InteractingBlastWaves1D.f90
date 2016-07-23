PROGRAM InteractingBlastWaves1D

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeInteractingBlastWaves1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'InteractingBlastWaves1D', &
           nX_Option &
             = [ 400, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 3, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = 1, &
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
           nStages_SSP_RK_Option &
             = 1 )

  CALL InitializeInteractingBlastWaves1D

  CALL EvolveFields &
         ( t_begin = 0.0_DP, t_end = 3.8d-2, dt_write = 1.0d-3, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM InteractingBlastWaves1D
