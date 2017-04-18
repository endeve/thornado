PROGRAM GravitationalCollapse1D

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    Kilometer, &
    MeV, &
    Millisecond
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE InitializationModule, ONLY: &
    InitializeGravitationalCollapse
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SSP_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
            = 'GravitationalCollapse1D', &
           nX_Option &
             = [ 256, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.0d0 * Kilometer, 0.0_DP, 0.0_DP ], &
           xR_Option &
             = [ 8.0d3 * Kilometer, Pi,     4.0_DP ], &
           zoomX_Option &
             = [ 1.020059256924853_DP, 1.0_DP, 1.0_DP ], &
           nE_Option &
             = 1, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 1.0d2 * MeV, &
           ZoomE_Option &
             = 1.0_DP, &
           nNodes_Option &
             = 2, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
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
           RadiationSolver_Option &
             = 'M1_DG', &
           RadiationRiemannSolver_Option &
             = 'HLL', &
           EvolveRadiation_Option &
             = .FALSE., &
           ApplySlopeLimiter_Option &
             = .TRUE., &
           BetaTVB_Option &
             = 0.0d0, &
           BetaTVD_Option &
             = 1.0d0, &
           ApplyPositivityLimiter_Option &
             = .FALSE., &
           EvolveGravity_Option &
             = .TRUE., &
           GravitySolver_Option &
             = 'Newtonian_Poseidon', &
           nStages_SSP_RK_Option &
             = 2 )

  CALL InitializeGravitationalCollapse &
         ( ProgenitorFile = 'WH07_15M_Sun.h5' )

  CALL EvolveFields &
         ( t_begin  = 0.0d+0 * Millisecond, &
           t_end    = 4.0d+2 * Millisecond, &
           dt_write = 5.0d-1 * Millisecond, &
           UpdateFields = SSP_RK )

  CALL FinalizeProgram

END PROGRAM GravitationalCollapse1D
