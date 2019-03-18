PROGRAM ShockEntropyWaveInteraction1D_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeShockEntropyWaveInteraction1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  REAL(DP) :: StartTime

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'ShockEntropyWaveInteraction1D_NuclearEOS', &
           nX_Option &
             = [ 100, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ - 10.0_DP, 0.0_DP, 0.0_DP ] * Kilometer, &
           xR_Option &
             = [ + 10.0_DP, 1.0_DP, 1.0_DP ] * Kilometer, &
           nNodes_Option &
             = 3, &
           ActivateUnits_Option &
             = .TRUE., &
           EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5', &
           FluidSolver_Option &
             = 'Euler_DG', &
           FluidRiemannSolver_Option &
             = 'HLLC', &
           EvolveFluid_Option &
             = .TRUE., &
           ApplySlopeLimiter_Option &
             = .TRUE., &
           BetaTVB_Option &
             = 0.0d2, &
           BetaTVD_Option &
             = 1.8d0, &
           ApplyPositivityLimiter_Option &
             = .TRUE., &
           nStages_SSP_RK_Option &
             = 3 )

  CALL InitializeShockEntropyWaveInteraction1D &
         ( StartTime = StartTime )

  CALL EvolveFields &
         ( t_begin  = StartTime, &
           t_end    = 1.25d-1 * Millisecond, &
           dt_write = 5.00d-3 * Millisecond, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM ShockEntropyWaveInteraction1D_NuclearEOS
