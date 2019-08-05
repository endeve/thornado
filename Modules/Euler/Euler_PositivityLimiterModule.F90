MODULE Euler_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_PositivityLimiter


#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_PositivityLimiterModule_NonRelativistic

#elif HYDRO_RELATIVISTIC

  USE Euler_PositivityLimiterModule_Relativistic

#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializePositivityLimiter
  PUBLIC :: Euler_FinalizePositivityLimiter
  PUBLIC :: Euler_ApplyPositivityLimiter


CONTAINS


  SUBROUTINE Euler_InitializePositivityLimiter &
    ( Min_1_Option, Min_2_Option, UsePositivityLimiter_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_InitializePositivityLimiter_NonRelativistic &
           ( Min_1_Option, Min_2_Option, &
             UsePositivityLimiter_Option, Verbose_Option )

#elif HYDRO_RELATIVISTIC

    CALL Euler_InitializePositivityLimiter_Relativistic &
           ( Min_1_Option, Min_2_Option, &
             UsePositivityLimiter_Option, Verbose_Option )

#endif

  END SUBROUTINE Euler_InitializePositivityLimiter


  SUBROUTINE Euler_FinalizePositivityLimiter

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_FinalizePositivityLimiter_NonRelativistic

#elif HYDRO_RELATIVISTIC

    CALL Euler_FinalizePositivityLimiter_Relativistic

#endif

  END SUBROUTINE Euler_FinalizePositivityLimiter


  SUBROUTINE Euler_ApplyPositivityLimiter &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ApplyPositivityLimiter_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#elif HYDRO_RELATIVISTIC

    CALL Euler_ApplyPositivityLimiter_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#endif

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE Euler_ApplyPositivityLimiter


END MODULE Euler_PositivityLimiterModule
