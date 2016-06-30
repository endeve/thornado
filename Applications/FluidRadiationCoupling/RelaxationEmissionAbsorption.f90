PROGRAM RelaxationEmissionAbsorption

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Kilometer, &
    MeV, &
    Gram, &
    Centimeter, &
    Kelvin, &
    Microsecond, &
    PlanckConstant, &
    SpeedOfLight
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram
  USE FluidRadiationCouplingInitializationModule, ONLY: &
    InitializeRelaxation
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    BackwardEuler

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option = 'RelaxationEmissionAbsorption', &
           nX_Option = [ 1, 1, 1 ], swX_Option = [ 0, 0, 0 ], &
           bcX_Option = [ 0, 0, 0 ], &
           xL_Option = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer, &
           xR_Option = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer, &
           nE_Option = 30, &
           eL_Option = 0.0d0 * MeV, &
           eR_Option = 3.0d2 * MeV, &
           ZoomE_Option = 1.1_DP, &
           nNodes_Option = 2, &
           ActivateUnits_Option = .TRUE., &
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
           EvolveFluid_Option = .TRUE., &
           EvolveRadiation_Option = .TRUE. )

  CALL InitializeRelaxation &
         ( D   = 1.00d14 * Gram / Centimeter**3, &
           T   = 1.00d11 * Kelvin, &
           Ye  = 0.30_DP, &
           E_0 = 5.0d1 * MeV )

  CALL EvolveFields &
         ( t_begin  = 0.0d+0 * Microsecond, &
           t_end    = 1.0d+1 * Microsecond, &
           dt_write = 1.0d-1 * Microsecond, &
           UpdateFields = BackwardEuler, &
           dt_fixed_Option = 1.0d-1 * Microsecond )

END PROGRAM RelaxationEmissionAbsorption
