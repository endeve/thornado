PROGRAM Sedov_SphericalSymmetry

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE RiemannProblemInitializationModule, ONLY: &
    InitializeSedov
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
            = 'Sedov_SphericalSymmetry', &
           nX_Option &
             = [ 1024, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 2, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 1.2_DP, Pi,     4.0_DP ], &
           nNodes_Option &
             = 1, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
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
           ApplySlopeLimiter_Option &
             = .TRUE., &
           BetaTVB_Option &
             = 1.0d2, &
           BetaTVD_Option &
             = 2.0d0, &
           ApplyPositivityLimiter_Option &
             = .TRUE., &
           nStages_SSP_RK_Option &
             = 1 )
  
  CALL InitializeSedov( E_0 = 1.0_DP, nDetonationCells_Option = 3 )

  CALL EvolveFields &
         ( t_begin  = 0.00_DP, &
           t_end    = 1.00_DP, &
           dt_write = 0.05_DP, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM Sedov_SphericalSymmetry
