MODULE TwoMoment_TallyModule

  USE KindModule, ONLY: &
    DP, Zero
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE RadiationFieldsModule, ONLY: &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nSpecies

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_TwoMoment
  PUBLIC :: FinalizeTally_TwoMoment
  PUBLIC :: ComputeTally_TwoMoment

  CHARACTER(256)        :: TallyFileName
  REAL(DP), ALLOCATABLE :: TwoMomentTally(:,:,:)

CONTAINS


  SUBROUTINE InitializeTally_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)           :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)           :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)           :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)           :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: FileUnit

    ALLOCATE( TwoMomentTally(1:nCR,1:nSpecies,0:1) )

    TallyFileName = '../Output/TwoMomentTally.dat'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

    WRITE( FileUnit, '(5(A14,x))' ) 'Time', 'N', 'G1', 'G2', 'G3'

    CLOSE( FileUnit )

    CALL ComputeTally_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, Time = Zero, &
             iState_Option = 0, DisplayTally_Option = .FALSE. )

    CALL ComputeTally_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, Time = Zero, &
             iState_Option = 1, DisplayTally_Option = .TRUE. )

  END SUBROUTINE InitializeTally_TwoMoment


  SUBROUTINE FinalizeTally_TwoMoment

    DEALLOCATE( TwoMomentTally )

  END SUBROUTINE FinalizeTally_TwoMoment


  SUBROUTINE ComputeTally_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, Time, &
      iState_Option, DisplayTally_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)           :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)           :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)           :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)           :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(in)           :: &
      Time
    INTEGER,  INTENT(in), OPTIONAL :: &
      iState_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      DisplayTally_Option

    LOGICAL :: DisplayTally
    INTEGER :: iState
    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iS

    IF( PRESENT( iState_Option ) )THEN
      iState = iState_Option
    ELSE
      iState = 1
    END IF

    IF( PRESENT( DisplayTally_Option ) )THEN
      DisplayTally = DisplayTally_Option
    ELSE
      DisplayTally = .FALSE.
    END IF

    TwoMomentTally(:,:,iState) = Zero

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)

                TwoMomentTally(iCR,iS,iState) &
                  = TwoMomentTally(iCR,iS,iState) &
                      + DOT_PRODUCT &
                          ( Weights_q, U(:,iZ1,iZ2,iZ3,iZ4,iCR,iS) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    IF( DisplayTally )THEN

      CALL DisplayTally_TwoMoment
      CALL WriteTally_TwoMoment( Time )

    END IF

  END SUBROUTINE ComputeTally_TwoMoment


  SUBROUTINE DisplayTally_TwoMoment

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'INFO: Two-Moment Tally'
    WRITE(*,*)
    WRITE(*,'(A8,A18,ES18.10E3)') &
      '', 'Total Number = ', &
      TwoMomentTally(iCR_N,1,1)
    WRITE(*,'(A8,A18,ES18.10E3)') &
      '', 'Relative Change = ', &
      ( TwoMomentTally(iCR_N,1,1) - TwoMomentTally(iCR_N,1,0) ) &
        / TwoMomentTally(iCR_N,1,0)
    WRITE(*,*)

  END SUBROUTINE DisplayTally_TwoMoment


  SUBROUTINE WriteTally_TwoMoment( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    OPEN( NEWUNIT=FileUnit, FILE=TRIM( TallyFileName ), ACCESS='APPEND' )

    WRITE( FileUnit, '(5(ES14.5,x))' ) &
      Time, &
      TwoMomentTally(iCR_N, 1,1), &
      TwoMomentTally(iCR_G1,1,1), &
      TwoMomentTally(iCR_G2,1,1), &
      TwoMomentTally(iCR_G3,1,1)

    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_TwoMoment


END MODULE TwoMoment_TallyModule
