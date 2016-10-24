PROGRAM ShockEntropyWaveInteraction1D

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeShockEntropyWaveInteraction1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'ShockEntropyWaveInteraction1D', &
           nX_Option &
             = [ 2000, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ - 5.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ + 5.0_DP, 1.0_DP, 1.0_DP ], &
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
           BetaTVB_Option &
             = 5.0d1, &
           nStages_SSP_RK_Option &
             = 1 )

  CALL InitializeShockEntropyWaveInteraction1D

  CALL EvolveFields &
         ( t_begin = 0.0_DP, t_end = 1.8d-0, dt_write = 5.0d-2, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM ShockEntropyWaveInteraction1D
