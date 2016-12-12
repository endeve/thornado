PROGRAM RiemannProblem1D_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    Gram, &
    Second, &
    Microsecond, &
    Erg
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeRiemannProblem1D_NuclearEOS
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK
  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RiemannProblem1D_NuclearEOS', &
           nX_Option &
             = [ 200, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer, &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer, &
           nNodes_Option &
             = 2, &
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
             = 2.0d0, &
           ApplyPositivityLimiter_Option &
             = .FALSE., &
           nStages_SSP_RK_Option &
             = 2 )

  CALL InitializeRiemannProblem1D_NuclearEOS &
         ( D_L  = 1.00d12 * Gram / Centimeter**3, &
           V_L  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer / Second, &
           P_L  = 1.0d32 * Erg / Centimeter**3, &
           Ye_L = 0.4_DP, &
           D_R  = 1.25d11 * Gram / Centimeter**3, &
           V_R  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer / Second, &
           P_R  = 1.0d31 * Erg / Centimeter**3, &
           Ye_R = 0.4_DP )

  CALL EvolveFields &
         ( t_begin  = 0.0d-0 * Microsecond, &
           t_end    = 2.0d-0 * Microsecond, &
           dt_write = 1.0d-1 * Microsecond, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D_NuclearEOS
