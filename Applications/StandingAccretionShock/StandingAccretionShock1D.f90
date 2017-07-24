PROGRAM StandingAccretionShock1D

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE InitializationModule, ONLY: &
    InitializeStandingAccretionShock
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK
  USE GravitySolutionModule, ONLY: &
    SolveGravity

  IMPLICIT NONE

  ! --- Problem Parameters ---

  REAL(DP), PARAMETER :: AccretionRate     = 4.0_DP * pi
  REAL(DP), PARAMETER :: GravitationalMass = 0.5_DP
  REAL(DP), PARAMETER :: ShockRadius       = 1.0_DP
  REAL(DP), PARAMETER :: PolytropicGamma   = 4.0_DP / 3.0_DP
  REAL(DP), PARAMETER :: MachNumber        = 3.0d1

  CALL InitializeProgram &
         ( ProgramName_Option &
            = 'StandingAccretionShock1D', &
           nX_Option &
             = [ 128, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.2_DP, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 2.0_DP, Pi,     4.0_DP ], &
           zoomX_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           EquationOfState_Option &
             = 'IDEAL', &
           Gamma_IDEAL_Option &
             = PolytropicGamma, &
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
             = 'Newtonian_PointMass', &
           PointMass_Option &
             = GravitationalMass, &
           nStages_SSP_RK_Option &
             = 2 )

  CALL InitializeStandingAccretionShock &
         ( AccretionRate, &
           GravitationalMass, &
           ShockRadius, &
           PolytropicGamma, &
           MachNumber )

  CALL SolveGravity

  CALL EvolveFields &
         ( t_begin  = 0.00_DP, &
           t_end    = 5.00_DP, &
           dt_write = 0.50_DP, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM StandingAccretionShock1D
