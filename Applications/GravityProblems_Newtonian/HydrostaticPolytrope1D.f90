PROGRAM HydrostaticPolytrope1D

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE GravityProblemsInitializationModule, ONLY: &
    InitializeHydrostaticPolytrope
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK
  USE GravitySolutionModule, ONLY: &
    SolveGravity

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
            = 'HydrostaticPolytrope1D', &
           nX_Option &
             = [ 32, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 0.8_DP, Pi,     4.0_DP ], &
           nNodes_Option &
             = 3, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           EquationOfState_Option &
             = 'IDEAL', &
           Gamma_IDEAL_Option &
             = 2.0_DP, &
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
             = 1.8d0, &
           ApplyPositivityLimiter_Option &
             = .TRUE., &
           EvolveGravity_Option &
             = .TRUE., &
           GravitySolver_Option &
             = 'Newtonian_Poseidon', &
           nStages_SSP_RK_Option &
             = 3 )

  CALL InitializeHydrostaticPolytrope

  CALL SolveGravity

  CALL EvolveFields &
         ( t_begin  = 0.0_DP, &
           t_end    = 1.0d+1, &
           dt_write = 1.0d-1, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM HydrostaticPolytrope1D
