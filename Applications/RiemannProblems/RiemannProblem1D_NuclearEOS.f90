PROGRAM RiemannProblem1D_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE UnitsModule, ONLY: &
    Kilometer

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RiemannProblem1D_NuclearEOS', &
           nX_Option &
             = [ 1024, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer, &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer, &
           nNodes_Option &
             = 1, &
           ActivateUnits_Option &
             = .TRUE., &
           EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5', &
           FluidSolver_Option &
             = 'Euler_DG', &
           FluidRiemannSolver_Option &
             = 'LLF', &
           EvolveFluid_Option &
             = .TRUE., &
           nStagesSSPRK_Option &
             = 1 )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D_NuclearEOS
