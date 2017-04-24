PROGRAM LinearWaves1D

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE SmoothProblemsInitializationModule, ONLY: &
    InitializeLinearWaves1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'LinearWaves1D', &
           nX_Option &
             = [ 16, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 1, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = 3, &
           EquationOfState_Option &
             = 'IDEAL', &
           Gamma_IDEAL_Option &
             = 5.0_DP / 3.0_DP, &
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
             = 1.0d0, &
           ApplyPositivityLimiter_Option &
             = .FALSE., &
           nStages_SSP_RK_Option &
             = 3 )

  CALL InitializeLinearWaves1D &
         ( 'AcousticWaves', Amplitude_Option = 1.0d-5 )

  CALL EvolveFields &
         ( t_begin  = 0.0d-0, &
           t_end    = 1.0d+1, &
           dt_write = 1.0d-1, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM LinearWaves1D
