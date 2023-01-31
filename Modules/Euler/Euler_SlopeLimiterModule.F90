MODULE Euler_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP
  USE Euler_BoundaryConditionsModule, ONLY: &
    iApplyBC_Euler_Both

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

  USE Euler_SlopeLimiterModule_Relativistic_TABLE

#else

  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE

#endif

#else

#ifdef HYDRO_RELATIVISTIC

  USE Euler_SlopeLimiterModule_Relativistic_IDEAL

#else

  USE Euler_SlopeLimiterModule_NonRelativistic_IDEAL

#endif

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler
  PUBLIC :: FinalizeSlopeLimiter_Euler
  PUBLIC :: ApplySlopeLimiter_Euler


CONTAINS


  SUBROUTINE InitializeSlopeLimiter_Euler &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, SlopeLimiterMethod_Option, &
      LimiterThresholdParameter_Option, &
      UseConservativeCorrection_Option, Verbose_Option )

    REAL(DP),         INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP),         INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,          INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, &
      UseConservativeCorrection_Option, &
      Verbose_Option
    REAL(DP),         INTENT(in), OPTIONAL :: &
      LimiterThresholdParameter_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      SlopeLimiterMethod_Option

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

    CALL InitializeSlopeLimiter_Euler_Relativistic_TABLE &
           ( UseSlopeLimiter_Option,           &
             SlopeLimiterMethod_Option,        &
             BetaTVD_Option,                   &
             BetaTVB_Option,                   &
             SlopeTolerance_Option,            &
             UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option,  &
             LimiterThresholdParameter_Option, &
             UseConservativeCorrection_Option, &
             Verbose_Option )

#else

    CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
           ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
             UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option, &
             LimiterThresholdParameter_Option, &
             UseConservativeCorrection_Option, &
             Verbose_Option )

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
           ( UseSlopeLimiter_Option,           &
             SlopeLimiterMethod_Option,        &
             BetaTVD_Option,                   &
             BetaTVB_Option,                   &
             SlopeTolerance_Option,            &
             UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option,  &
             LimiterThresholdParameter_Option, &
             UseConservativeCorrection_Option, &
             Verbose_Option )

#else

    CALL InitializeSlopeLimiter_Euler_NonRelativistic_IDEAL &
           ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
             UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option, &
             LimiterThresholdParameter_Option, &
             UseConservativeCorrection_Option, &
             Verbose_Option )

#endif

#endif

  END SUBROUTINE InitializeSlopeLimiter_Euler


  SUBROUTINE FinalizeSlopeLimiter_Euler

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

    CALL FinalizeSlopeLimiter_Euler_Relativistic_TABLE

#else

    CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

#else

    CALL FinalizeSlopeLimiter_Euler_NonRelativistic_IDEAL

#endif

#endif

  END SUBROUTINE FinalizeSlopeLimiter_Euler


  SUBROUTINE ApplySlopeLimiter_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

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

    iApplyBC = iApplyBC_Euler_Both
    IF( PRESENT( iApplyBC_Option ) ) &
      iApplyBC = iApplyBC_Option

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

    CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
             SuppressBC_Option = SuppressBC, iApplyBC_Option = iApplyBC )

#else

    CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
             SuppressBC_Option = SuppressBC, iApplyBC_Option = iApplyBC )

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
             SuppressBC_Option = SuppressBC, iApplyBC_Option = iApplyBC )

#else

    CALL ApplySlopeLimiter_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
             SuppressBC_Option = SuppressBC, iApplyBC_Option = iApplyBC )

#endif

#endif

  END SUBROUTINE ApplySlopeLimiter_Euler


END MODULE Euler_SlopeLimiterModule
