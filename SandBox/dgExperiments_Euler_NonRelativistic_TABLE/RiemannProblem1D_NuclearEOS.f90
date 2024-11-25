PROGRAM RiemannProblem1D_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Gram, &
    Second, &
    Millisecond, &
    Microsecond, &
    Erg
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeRiemannProblem1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  REAL(DP) :: StartTime

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RiemannProblem1D_NuclearEOS', &
           nX_Option &
             = [ 128, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer, &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer, &
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

  CALL InitializeRiemannProblem1D &
         ( D_L  = 1.00d12 * Gram / Centimeter**3, &
           V_L  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer / Second, &
           P_L  = 1.0d32 * Erg / Centimeter**3, &
           Ye_L = 0.3_DP, &
           D_R  = 1.25d11 * Gram / Centimeter**3, &
           V_R  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer / Second, &
           P_R  = 1.0d31 * Erg / Centimeter**3, &
           Ye_R = 0.3_DP, &
           StartTime = StartTime )

  CALL EvolveFields &
         ( t_begin  = StartTime, &
           t_end    = 2.5d-0 * Microsecond, &
           dt_write = 5.0d-2 * Microsecond, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D_NuclearEOS
