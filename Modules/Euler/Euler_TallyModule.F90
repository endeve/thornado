MODULE Euler_TallyModule

  USE KindModule, ONLY: &
    DP

  USE Euler_TallyModule_NonRelativistic_IDEAL
  USE Euler_TallyModule_NonRelativistic_TABLE
  USE Euler_TallyModule_Relativistic_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializeTally
  PUBLIC :: Euler_FinalizeTally
  PUBLIC :: Euler_ComputeTally


CONTAINS


  SUBROUTINE Euler_InitializeTally( iX_B0, iX_E0, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    REAL(DP), INTENT(in) :: &
      U(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL Euler_InitializeTally_NonRelativistic_TABLE( iX_B0, iX_E0, G, U )

#elif defined HYDRO_NONRELATIVISTIC

    CALL Euler_InitializeTally_NonRelativistic( iX_B0, iX_E0, G, U )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_InitializeTally_Relativistic( iX_B0, iX_E0, G, U )

#else

    CALL Euler_InitializeTally_NonRelativistic( iX_B0, iX_E0, G, U )

#endif

  END SUBROUTINE Euler_InitializeTally


  SUBROUTINE Euler_FinalizeTally

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL Euler_FinalizeTally_NonRelativistic_TABLE

#elif defined HYDRO_NONRELATIVISTIC

    CALL Euler_FinalizeTally_NonRelativistic

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_FinalizeTally_Relativistic

#else

    CALL Euler_FinalizeTally_NonRelativistic

#endif

  END SUBROUTINE Euler_FinalizeTally


  SUBROUTINE Euler_ComputeTally &
    ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    REAL(DP), INTENT(in) :: &
      U(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    REAL(DP), INTENT(in) :: &
      Time
    INTEGER,  INTENT(in), OPTIONAL :: &
      iState_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      DisplayTally_Option

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL Euler_ComputeTally_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )


#elif defined HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeTally_NonRelativistic &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ComputeTally_Relativistic &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

#else

    CALL Euler_ComputeTally_NonRelativistic &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

#endif

  END SUBROUTINE Euler_ComputeTally


END MODULE Euler_TallyModule
