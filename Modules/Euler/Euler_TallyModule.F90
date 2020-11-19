MODULE Euler_TallyModule

  USE KindModule, ONLY: &
    DP

#ifdef HYDRO_RELATIVISTIC

  USE Euler_TallyModule_Relativistic

#else

  USE Euler_TallyModule_NonRelativistic

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler
  PUBLIC :: FinalizeTally_Euler
  PUBLIC :: ComputeTally_Euler


CONTAINS


  SUBROUTINE InitializeTally_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressTally_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(in)           :: &
      U(:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressTally_Option

    LOGICAL :: SuppressTally

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

#ifdef HYDRO_RELATIVISTIC

    CALL InitializeTally_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressTally )

#else

    CALL InitializeTally_Euler_NonRelativistic &
           ( iX_B0, iX_E0,             &
             G(:,iX_B0(1):iX_E0(1),    &
                 iX_B0(2):iX_E0(2),    &
                 iX_B0(3):iX_E0(3),:), &
             U(:,iX_B0(1):iX_E0(1),    &
                 iX_B0(2):iX_E0(2),    &
                 iX_B0(3):iX_E0(3),:), &
             SuppressTally )

#endif

  END SUBROUTINE InitializeTally_Euler


  SUBROUTINE FinalizeTally_Euler

#ifdef HYDRO_RELATIVISTIC

    CALL FinalizeTally_Euler_Relativistic

#else

    CALL FinalizeTally_Euler_NonRelativistic

#endif

  END SUBROUTINE FinalizeTally_Euler


  SUBROUTINE ComputeTally_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time, &
      iState_Option, DisplayTally_Option )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      Time
    INTEGER,  INTENT(in), OPTIONAL :: &
      iState_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      DisplayTally_Option

#ifdef HYDRO_RELATIVISTIC

    CALL ComputeTally_Euler_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time )

#else

    CALL ComputeTally_Euler_NonRelativistic &
           ( iX_B0, iX_E0,             &
             G(:,iX_B0(1):iX_E0(1),    &
                 iX_B0(2):iX_E0(2),    &
                 iX_B0(3):iX_E0(3),:), &
             U(:,iX_B0(1):iX_E0(1),    &
                 iX_B0(2):iX_E0(2),    &
                 iX_B0(3):iX_E0(3),:), &
              Time, iState_Option, DisplayTally_Option )

#endif

  END SUBROUTINE ComputeTally_Euler


END MODULE Euler_TallyModule
