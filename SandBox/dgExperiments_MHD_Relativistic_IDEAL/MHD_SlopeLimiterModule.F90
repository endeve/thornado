MODULE MHD_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP
  USE MHD_BoundaryConditionsModule, ONLY: &
    iApplyBC_MHD_Both
  USE MHD_SlopeLimiterModule_Relativistic_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_MHD
  PUBLIC :: FinalizeSlopeLimiter_MHD
  PUBLIC :: ApplySlopeLimiter_MHD


CONTAINS


  SUBROUTINE InitializeSlopeLimiter_MHD &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, SlopeLimiterMethod_Option, &
      LimiterThresholdParameter_Option, &
      UseConservativeCorrection_Option, Verbose_Option, &
      EvolveOnlyMagnetic_Option )

    REAL(DP),         INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP),         INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,          INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, &
      UseConservativeCorrection_Option, &
      Verbose_Option, &
      EvolveOnlyMagnetic_Option
    REAL(DP),         INTENT(in), OPTIONAL :: &
      LimiterThresholdParameter_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      SlopeLimiterMethod_Option

    CALL InitializeSlopeLimiter_MHD_Relativistic_IDEAL &
           ( UseSlopeLimiter_Option,           &
             SlopeLimiterMethod_Option,        &
             BetaTVD_Option,                   &
             BetaTVB_Option,                   &
             SlopeTolerance_Option,            &
             UseConservativeCorrection_Option, &
             UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option,  &
             LimiterThresholdParameter_Option, &
             Verbose_Option,                   &
             EvolveOnlyMagnetic_Option )

  END SUBROUTINE InitializeSlopeLimiter_MHD


  SUBROUTINE FinalizeSlopeLimiter_MHD

    CALL FinalizeSlopeLimiter_MHD_Relativistic_IDEAL

  END SUBROUTINE FinalizeSlopeLimiter_MHD


  SUBROUTINE ApplySlopeLimiter_MHD &
    ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    REAL(DP), INTENT(in)           :: t
    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    LOGICAL :: SuppressBC
    INTEGER :: iApplyBC(3)

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    iApplyBC = iApplyBC_MHD_Both
    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

    CALL ApplySlopeLimiter_MHD_Relativistic_IDEAL &
           ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
             SuppressBC_Option = SuppressBC, iApplyBC_Option = iApplyBC )

  END SUBROUTINE ApplySlopeLimiter_MHD


END MODULE MHD_SlopeLimiterModule
