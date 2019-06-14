MODULE Euler_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP







  USE Euler_SlopeLimiterModule_Relativistic




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



    CALL Euler_InitializeSlopeLimiter_Relativistic &
           ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
             UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
             UseTroubledCellIndicator_Option, &
             LimiterThresholdParameter_Option, Verbose_Option )



  END SUBROUTINE Euler_InitializeSlopeLimiter


  SUBROUTINE Euler_FinalizeSlopeLimiter







    CALL Euler_FinalizeSlopeLimiter_Relativistic



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








    CALL Euler_ApplySlopeLimiter_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressBC )



  END SUBROUTINE Euler_ApplySlopeLimiter


END MODULE Euler_SlopeLimiterModule

