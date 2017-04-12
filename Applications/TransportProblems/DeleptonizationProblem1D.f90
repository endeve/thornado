PROGRAM DeleptonizationProblem1D

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    Kilometer, &
    MeV, &
    Kelvin, &
    Microsecond, &
    Millisecond
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TransportProblemsInitializationModule, ONLY: &
    InitializeDeleptonizationProblem1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SI_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'DeleptonizationProblem1D', &
           nX_Option &
             = [ 100, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.0d0 * Kilometer, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 2.0d2 * Kilometer, Pi,    4.0d0 ], &
           nE_Option &
             = 12, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 5.0d1 * MeV, &
           ZoomE_Option &
             = 1.15_DP, &
           nNodes_Option &
             = 3, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5', &
           Opacity_Option &
             = 'TABLE', &
           OpacityTableName_Option &
             = 'OpacityTable.h5', &
           FluidRadiationCoupling_Option &
             = 'EmissionAbsorption', &
           EvolveFluid_Option &
             = .TRUE., &
           RadiationSolver_Option &
             = 'M1_DG', &
           RadiationRiemannSolver_Option &
             = 'HLL', &
           EvolveRadiation_Option &
             = .TRUE., &
           ApplySlopeLimiter_Option &
             = .TRUE., &
           BetaTVB_Option &
             = 0.0d0, &
           BetaTVD_Option &
             = 1.8d0, &
           ApplyPositivityLimiter_Option &
             = .TRUE., &
           nStages_SI_RK_Option &
             = 2 )

  CALL InitializeDeleptonizationProblem1D &
         ( Temperature_Option &
             = 1.0d11 * Kelvin, &
           ElectronFraction_Option &
             = 0.3_DP )

  CALL EvolveFields &
         ( t_begin  = 0.0d+0 * Millisecond, &
           t_end    = 2.5d-1 * Millisecond, &
           dt_write = 5.0d-2 * Millisecond, &
           UpdateFields = SI_RK )

  CALL FinalizeProgram

END PROGRAM DeleptonizationProblem1D
