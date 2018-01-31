MODULE BoundaryConditionsModule_Beta

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    bcX, swX
  USE FluidFieldsModule, ONLY: &
    nCF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Fluid

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Fluid &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ApplyBC_Fluid_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL ApplyBC_Fluid_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL ApplyBC_Fluid_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

  END SUBROUTINE ApplyBoundaryConditions_Fluid


  SUBROUTINE ApplyBC_Fluid_X1( iX_B0, iX_E0, iX_B1, iX_E1, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3

    SELECT CASE ( bcX(1) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

      DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
            DO iX1 = 1, swX(1)

              ! --- Inner Boundary ---

              U(:,iX_B0(1)-iX1,iX2,iX3,iCF) &
                = U(:,iX_E0(1)-(iX1-1),iX2,iX3,iCF)

              ! --- Outer Boundary ---

              U(:,iX_E0(1)+iX1,iX2,iX3,iCF) &
                = U(:,iX_B0(1)+(iX1-1),iX2,iX3,iCF)

            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

      DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
            DO iX1 = 1, swX(1)

              ! --- Inner Boundary ---

              U(:,iX_B0(1)-iX1,iX2,iX3,iCF) &
                = U(:,iX_B0(1),iX2,iX3,iCF)

              ! --- Outer Boundary ---

              U(:,iX_E0(1)+iX1,iX2,iX3,iCF) &
                = U(:,iX_E0(1),iX2,iX3,iCF)

            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Fluid X1: ', bcX(1)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_Fluid_X1


  SUBROUTINE ApplyBC_Fluid_X2( iX_B0, iX_E0, iX_B1, iX_E1, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3

    SELECT CASE ( bcX(2) )

    CASE ( 0 ) ! No Boundary Condition

    CASE ( 1 ) ! Periodic

      DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = 1, swX(2)
            DO iX1 = iX_B0(1), iX_E0(1)

              ! --- Inner Boundary ---

              U(:,iX1,iX_B0(2)-iX2,iX3,iCF) &
                = U(:,iX1,iX_E0(2)-(iX2-1),iX3,iCF)

              ! --- Outer Boundary ---

              U(:,iX1,iX_E0(2)+iX2,iX3,iCF) &
                = U(:,iX1,iX_B0(2)+(iX2-1),iX3,iCF)

            END DO
          END DO
        END DO
      END DO

    CASE ( 2 ) ! Homogeneous

      DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = 1, swX(2)
            DO iX1 = iX_B0(1), iX_E0(1)

              ! --- Inner Boundary ---

              U(:,iX1,iX_B0(2)-iX2,iX3,iCF) &
                = U(:,iX1,iX_B0(2),iX3,iCF)

              ! --- Outer Boundary ---

              U(:,iX1,iX_E0(2)+iX2,iX3,iCF) &
                = U(:,iX1,iX_E0(2),iX3,iCF)

            END DO
          END DO
        END DO
      END DO

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A5,A45,I2.2)') &
        '', 'Invalid Boundary Condition for Fluid X2: ', bcX(2)
      STOP

    END SELECT

  END SUBROUTINE ApplyBC_Fluid_X2


  SUBROUTINE ApplyBC_Fluid_X3( iX_B0, iX_E0, iX_B1, iX_E1, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

  END SUBROUTINE ApplyBC_Fluid_X3


END MODULE BoundaryConditionsModule_Beta
