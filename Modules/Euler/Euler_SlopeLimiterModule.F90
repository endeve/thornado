MODULE Euler_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_SlopeLimiter

#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_SlopeLimiterModule_NonRelativistic

#elif HYDRO_RELATIVISTIC

  USE Euler_SlopeLimiterModule_Relativistic

#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializeSlopeLimiter
  PUBLIC :: Euler_FinalizeSlopeLimiter
  PUBLIC :: Euler_ApplySlopeLimiter


CONTAINS


  SUBROUTINE Euler_InitializeSlopeLimiter &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, LimiterThresholdParameter_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option, &
      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      LimiterThresholdParameter_Option

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_InitializeSlopeLimiter_NonRelativistic &
           ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
             UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option, &
             LimiterThresholdParameter_Option, Verbose_Option )

#elif HYDRO_RELATIVISTIC

    CALL Euler_InitializeSlopeLimiter_Relativistic &
           ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
             UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option, &
             LimiterThresholdParameter_Option, Verbose_Option )

#endif

  END SUBROUTINE Euler_InitializeSlopeLimiter


  SUBROUTINE Euler_FinalizeSlopeLimiter

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_FinalizeSlopeLimiter_NonRelativistic

#elif HYDRO_RELATIVISTIC

    CALL Euler_FinalizeSlopeLimiter_Relativistic

#endif

  END SUBROUTINE Euler_FinalizeSlopeLimiter


  SUBROUTINE Euler_ApplySlopeLimiter &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL :: SuppressBC

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ApplySlopeLimiter_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressBC )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ApplySlopeLimiter_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressBC )

#endif

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE Euler_ApplySlopeLimiter


END MODULE Euler_SlopeLimiterModule
