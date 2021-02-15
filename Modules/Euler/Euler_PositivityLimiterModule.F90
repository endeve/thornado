MODULE Euler_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

  USE Euler_PositivityLimiterModule_Relativistic_TABLE

#else

  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE

#endif

#else

#ifdef HYDRO_RELATIVISTIC

  USE Euler_PositivityLimiterModule_Relativistic_IDEAL

#else

  USE Euler_PositivityLimiterModule_NonRelativistic_IDEAL

#endif

#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler
  PUBLIC :: FinalizePositivityLimiter_Euler
  PUBLIC :: ApplyPositivityLimiter_Euler


CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler &
    ( UsePositivityLimiter_Option, Verbose_Option, &
      Min_1_Option, Min_2_Option, Min_3_Option, &
      Max_1_Option, Max_2_Option, Max_3_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option, Max_2_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_3_Option, Max_3_Option

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

    CALL InitializePositivityLimiter_Euler_Relativistic_TABLE &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option, Min_3_Option, &
             Max_1_Option, Max_2_Option, Max_3_Option )

#else

    CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option, Min_3_Option, &
             Max_1_Option, Max_2_Option, Max_3_Option )

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL InitializePositivityLimiter_Euler_Relativistic_IDEAL &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option )

#else

    CALL InitializePositivityLimiter_Euler_NonRelativistic_IDEAL &
           ( UsePositivityLimiter_Option, Verbose_Option, &
             Min_1_Option, Min_2_Option )

#endif

#endif

  END SUBROUTINE InitializePositivityLimiter_Euler


  SUBROUTINE FinalizePositivityLimiter_Euler

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

    CALL FinalizePositivityLimiter_Euler_Relativistic_TABLE

#else

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL FinalizePositivityLimiter_Euler_Relativistic_IDEAL

#else

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL

#endif

#endif

  END SUBROUTINE FinalizePositivityLimiter_Euler


  SUBROUTINE ApplyPositivityLimiter_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)              :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)              :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)           :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout), OPTIONAL :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#else

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#else

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

#endif

#endif

  END SUBROUTINE ApplyPositivityLimiter_Euler


END MODULE Euler_PositivityLimiterModule
