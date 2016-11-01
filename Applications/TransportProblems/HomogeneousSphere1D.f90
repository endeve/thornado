PROGRAM HomogeneousSphere1D

  USE KindModule, ONLY: &
    DP, Pi, TwoPi
  USE UnitsModule, ONLY: &
    Kilometer, &
    MeV, &
    Microsecond
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TransportProblemsInitializationModule, ONLY: &
    InitializeHomogeneousSphere1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SI_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'HomogeneousSphere1D', &
           nX_Option &
             = [ 100, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.0d0 * Kilometer, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 2.0d2 * Kilometer, Pi,    TwoPi ], &
           nE_Option &
             = 10, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 1.0d2 * MeV, &
           ZoomE_Option &
             = 1.25_DP, &
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
           Opacity_Option &
             = 'TABLE', &
           OpacityTableName_Option &
             = 'OpacityTable.h5', &
           FluidRadiationCoupling_Option &
             = 'EmissionAbsorption', &
           EvolveFluid_Option &
             = .FALSE., &
           RadiationSolver_Option &
             = 'M1_DG', &
           EvolveRadiation_Option &
             = .TRUE., &
           nStages_SI_RK_Option &
             = 2 )

  CALL InitializeHomogeneousSphere1D

  CALL EvolveFields &
         ( t_begin  = 0.0d+0 * Microsecond, &
           t_end    = 1.0d+0 * Microsecond, &
           dt_write = 1.0d+1 * Microsecond, &
           UpdateFields = SI_RK )

  CALL FinalizeProgram

END PROGRAM HomogeneousSphere1D
