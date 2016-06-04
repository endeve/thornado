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
         ( ProgramName_Option = 'LinearWaves1D', &
           nX_Option &
             = [ 64, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option = [ 1, 0, 0 ], &
           xL_Option = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option = 3, &
           EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = 5.0_DP / 3.0_DP, &
           FluidSolver_Option = 'Euler_DG', &
           FluidRiemannSolver_Option = 'HLLC', &
           EvolveFluid_Option = .TRUE., &
           nStagesSSPRK_Option = 3 )

  CALL InitializeLinearWaves1D( 'EntropyWaves' )

  CALL EvolveFields &
         ( t_begin = 0.0_DP, t_end = 2.0d-0, dt_write = 1.0d-1, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM LinearWaves1D
