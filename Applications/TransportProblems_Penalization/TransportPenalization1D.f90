PROGRAM TransportPenalization1D

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
  USE InitializationModule, ONLY: &
    InitializeTransportProblem1D
  USE InputOutputModule, ONLY: &
    WriteFields1D

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'TransportPenalization1D', &
           nX_Option &
             = [ 64, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.0d0 * Kilometer, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 1.0d2 * Kilometer, Pi,    4.0d0 ], &
           zoomX_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ], &
           nE_Option &
             = 20, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 3.0d2 * MeV, &
           ZoomE_Option &
             = 1.239020750754000_DP, &
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
           ApplyPositivityLimiter_Option &
             = .TRUE. )

  CALL InitializeTransportProblem1D

  CALL WriteFields1D &
         ( Time = 0.0_DP, WriteFluidFields_Option = .TRUE., &
           WriteRadiationFields_Option = .TRUE. )

  CALL FinalizeProgram

END PROGRAM TransportPenalization1D
