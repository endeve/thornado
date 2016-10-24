PROGRAM RiemannProblem1D_CylindricalSymmetry

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
            = 'RiemannProblem1D_CylindricalSymmetry', &
           nX_Option &
             = [ 2048, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.0_DP, 2 * Pi, 1.0_DP ], &
           nNodes_Option &
             = 1, &
           CoordinateSystem_Option &
             = 'CYLINDRICAL', &
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
  
  CALL InitializeRiemannProblem1D &
         ( D_L = 0.00_DP, V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ], P_L = 0.00_DP, &
           D_R = 1.00_DP, V_R = [ -1.0_DP, 0.0_DP, 0.0_DP ], P_R = 1.0D-6, X_D_Option = 0.0_DP )

  CALL EvolveFields &
         ( t_begin = 0.0_DP, t_end = 1.0_DP, dt_write = 0.01_DP, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM RiemannProblem1D_CylindricalSymmetry
