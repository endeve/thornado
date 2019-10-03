MODULE Euler_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_PositivityLimiter

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE

#elif defined HYDRO_NONRELATIVISTIC

  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL

#elif defined HYDRO_RELATIVISTIC

  USE Euler_PositivityLimiterModule_Relativistic_IDEAL

#else

  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL

#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializePositivityLimiter
  PUBLIC :: Euler_FinalizePositivityLimiter
  PUBLIC :: Euler_ApplyPositivityLimiter


CONTAINS


  SUBROUTINE Euler_InitializePositivityLimiter &
    ( UsePositivityLimiter_Option, Verbose_Option, &
      Min_1_Option, Min_2_Option, Min_3_Option, &
      Max_1_Option, Max_2_Option, Max_3_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option, Max_2_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_3_Option, Max_3_Option

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL Euler_InitializePositivityLimiter_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option, Min_3_Option, &
             Max_1_Option, Max_2_Option, Max_3_Option )

#elif defined HYDRO_NONRELATIVISTIC

    CALL Euler_InitializePositivityLimiter_NonRelativistic &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_InitializePositivityLimiter_Relativistic &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option )

#else

    CALL Euler_InitializePositivityLimiter_NonRelativistic &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option )

#endif

  END SUBROUTINE Euler_InitializePositivityLimiter


  SUBROUTINE Euler_FinalizePositivityLimiter

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL Euler_FinalizePositivityLimiter_NonRelativistic_TABLE

#elif defined HYDRO_NONRELATIVISTIC

    CALL Euler_FinalizePositivityLimiter_NonRelativistic

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_FinalizePositivityLimiter_Relativistic

#else

    CALL Euler_FinalizePositivityLimiter_NonRelativistic

#endif

  END SUBROUTINE Euler_FinalizePositivityLimiter


  SUBROUTINE Euler_ApplyPositivityLimiter &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, iErr_Option )

    INTEGER,  INTENT(in)             :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)             :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)          :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    INTEGER, INTENT(inout), OPTIONAL :: &
      iErr_Option

    INTEGER :: iErr = 0
    IF( PRESENT( iErr_Option ) ) iErr = iErr_Option

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL Euler_ApplyPositivityLimiter_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#elif defined HYDRO_NONRELATIVISTIC

    CALL Euler_ApplyPositivityLimiter_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ApplyPositivityLimiter_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, iErr_Option = iErr )

#else

    CALL Euler_ApplyPositivityLimiter_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#endif

    IF( PRESENT( iErr_Option ) ) iErr_Option = iErr

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE Euler_ApplyPositivityLimiter


END MODULE Euler_PositivityLimiterModule
