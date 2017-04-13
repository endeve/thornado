PROGRAM RelaxationNES

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Kilometer, &
    MeV, &
    Gram, &
    Centimeter, &
    Millisecond
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE FluidRadiationCouplingInitializationModule, ONLY: &
    InitializeRelaxationNES
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    BackwardEuler

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'RelaxationNES', &
           nX_Option &
             = [ 100, 1, 1 ], &
           swX_Option &
             = [ 0, 0, 0 ], &
           bcX_Option &
             = [ 0, 0, 0 ], &
           xL_Option &
             = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer, &
           xR_Option &
             = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer, &
           nE_Option &
             = 20, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 3.0d2 * MeV, &
           ZoomE_Option &
             = 1.2_DP, &
           nNodes_Option &
             = 3, &
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
             = 'NES', &
           EvolveFluid_Option &
             = .TRUE., &
           EvolveRadiation_Option &
             = .TRUE. )

  CALL InitializeRelaxationNES &
         ( Density &
             = 1.0d14 * Gram / Centimeter**3, &
           Temperature &
             = 2.1d1 * MeV, &
           ElectronFraction &
             = 0.25_DP )

  CALL EvolveFields &
         ( t_begin         = 0.0d+0 * Millisecond, &
           t_end           = 2.0d+0 * Millisecond, &
           dt_write        = 1.0d-2 * Millisecond, &
           UpdateFields    = BackwardEuler, &
           dt_fixed_Option = 1.0d-3 * Millisecond )

  CALL FinalizeProgram

END PROGRAM RelaxationNES
