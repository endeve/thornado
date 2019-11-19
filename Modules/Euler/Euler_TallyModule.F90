MODULE Euler_TallyModule

  USE KindModule, ONLY: &
    DP

  USE Euler_TallyModule_NonRelativistic_IDEAL
  USE Euler_TallyModule_NonRelativistic_TABLE
  USE Euler_TallyModule_Relativistic_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler
  PUBLIC :: FinalizeTally_Euler
  PUBLIC :: ComputeTally_Euler


CONTAINS


  SUBROUTINE InitializeTally_Euler( iX_B0, iX_E0, G, U, SuppressTally_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)           :: &
      G(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    REAL(DP), INTENT(in)           :: &
      U(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressTally_Option

    LOGICAL :: SuppressTally

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL InitializeTally_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, G, U, SuppressTally )

#elif defined HYDRO_NONRELATIVISTIC

    CALL InitializeTally_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, G, U, SuppressTally )

#elif defined HYDRO_RELATIVISTIC

    CALL InitializeTally_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, G, U, SuppressTally )

#else

    CALL InitializeTally_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, G, U, SuppressTally )

#endif

  END SUBROUTINE InitializeTally_Euler


  SUBROUTINE FinalizeTally_Euler

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL FinalizeTally_Euler_NonRelativistic_TABLE

#elif defined HYDRO_NONRELATIVISTIC

    CALL FinalizeTally_Euler_NonRelativistic_IDEAL

#elif defined HYDRO_RELATIVISTIC

    CALL FinalizeTally_Euler_Relativistic_IDEAL

#else

    CALL FinalizeTally_Euler_NonRelativistic_IDEAL

#endif

  END SUBROUTINE FinalizeTally_Euler


  SUBROUTINE ComputeTally_Euler &
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

    CALL ComputeTally_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )


#elif defined HYDRO_NONRELATIVISTIC

    CALL ComputeTally_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeTally_Euler_Relativistic_IDEAL &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

#else

    CALL ComputeTally_Euler_NonRelativistic_IDEAL &
           ( iX_B0, iX_E0, G, U, Time, iState_Option, DisplayTally_Option )

#endif

  END SUBROUTINE ComputeTally_Euler


END MODULE Euler_TallyModule
