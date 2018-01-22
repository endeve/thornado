MODULE BoundaryConditionsModule_Beta

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    bcZ, swZ
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Radiation

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Radiation &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    CALL ApplyBC_Radiation_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    CALL ApplyBC_Radiation_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    CALL ApplyBC_Radiation_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

  END SUBROUTINE ApplyBoundaryConditions_Radiation


  SUBROUTINE ApplyBC_Radiation_X1( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iS, iCR, iZ1, iZ2, iZ3, iZ4

    SELECT CASE ( bcZ(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary ---

                  U(:,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ_E0(2)-(iZ2-1),iZ3,iZ4,iCR,iS)

                  ! --- Outer Boundary ---

                  U(:,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ_B0(2)+(iZ2-1),iZ3,iZ4,iCR,iS)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = iZ_B0(3), iZ_E0(3)
              DO iZ2 = 1, swZ(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary ---

                  U(:,iZ1,iZ_B0(2)-iZ2,iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ_B0(2),iZ3,iZ4,iCR,iS)

                  ! --- Outer Boundary ---

                  U(:,iZ1,iZ_E0(2)+iZ2,iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ_E0(2),iZ3,iZ4,iCR,iS)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Radiation X1: ', bcZ(2)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_Radiation_X1


  SUBROUTINE ApplyBC_Radiation_X2( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iS, iCR, iZ1, iZ2, iZ3, iZ4

    SELECT CASE ( bcZ(3) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = 1, swZ(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary ---

                  U(:,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ2,iZ_E0(3)-(iZ3-1),iZ4,iCR,iS)

                  ! --- Outer Boundary ---

                  U(:,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ2,iZ_B0(3)+(iZ3-1),iZ4,iCR,iS)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iZ4 = iZ_B0(4), iZ_E0(4)
            DO iZ3 = 1, swZ(3)
              DO iZ2 = iZ_B0(2), iZ_E0(2)
                DO iZ1 = iZ_B0(1), iZ_E0(1)

                  ! --- Inner Boundary ---

                  U(:,iZ1,iZ2,iZ_B0(3)-iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ2,iZ_B0(3),iZ4,iCR,iS)

                  ! --- Outer Boundary ---

                  U(:,iZ1,iZ2,iZ_E0(3)+iZ3,iZ4,iCR,iS) &
                    = U(:,iZ1,iZ2,iZ_E0(3),iZ4,iCR,iS)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Radiation X2: ', bcZ(3)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_Radiation_X2


  SUBROUTINE ApplyBC_Radiation_X3( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(inout) :: &
      U(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

  END SUBROUTINE ApplyBC_Radiation_X3


END MODULE BoundaryConditionsModule_Beta
