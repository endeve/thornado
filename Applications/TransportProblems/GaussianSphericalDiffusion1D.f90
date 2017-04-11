PROGRAM GaussianSphericalDiffusion1D

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Kilometer, &
    MeV, &
    Microsecond, &
    Millisecond
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TransportProblemsInitializationModule, ONLY: &
    InitializeGaussianSphericalDiffusion1D
  USE TimeSteppingModule, ONLY: &
    EvolveFields, &
    SI_RK

  IMPLICIT NONE

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'GaussianSphericalDiffusion1D', &
           nX_Option &
             = [ 32, 1, 1 ], &
           swX_Option &
             = [ 1, 0, 0 ], &
           bcX_Option &
             = [ 10, 0, 0 ], &
           xL_Option &
             = [ 0.0d0 * Kilometer, 0.0d0, 0.0d0 ], &
           xR_Option &
             = [ 5.0d1 * Kilometer, Pi,    4.0d0 ], &
           nE_Option &
             = 1, &
           eL_Option &
             = 0.0d0 * MeV, &
           eR_Option &
             = 2.0d2 * MeV, &
           ZoomE_Option &
             = 1.0_DP, &
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
             = 'ElasticScattering', &
           EvolveFluid_Option &
             = .FALSE., &
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

  CALL InitializeGaussianSphericalDiffusion1D &
         ( t_0 = 3.0d-0 * Millisecond, &
           BackgroundConditions_Option = '01' )

  CALL EvolveFields &
         ( t_begin  = 0.0d+0 * Millisecond, &
           t_end    = 1.0d+1 * Millisecond, &
           dt_write = 1.0d-0 * Millisecond, &
           UpdateFields = SI_RK )

  CALL FinalizeProgram

END PROGRAM GaussianSphericalDiffusion1D
