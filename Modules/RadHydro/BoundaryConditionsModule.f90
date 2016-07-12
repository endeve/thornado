MODULE BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, bcX, nE
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, uPF, nPF, uAF, nAF
  USE RadiationFieldsModule, ONLY: &
    uCR, nCR, nSpecies
  USE ApplicationBoundaryConditionsModule, ONLY: &
    ApplyApplicationBoundaryConditions_Radiation_X1

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Fluid
  PUBLIC :: ApplyBoundaryConditions_Radiation

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Fluid

    CALL ApplyBoundaryConditions_Fluid_X1

  END SUBROUTINE ApplyBoundaryConditions_Fluid


  SUBROUTINE ApplyBoundaryConditions_Fluid_X1

    INTEGER :: iX2, iX3, iCF, iPF, iAF

    SELECT CASE ( bcX(1) )

      CASE ( 0 ) ! No Boundary Condition

      CASE ( 1 ) ! Periodic

        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)

            DO iCF = 1, nCF

              uCF(:,0,iX2,iX3,iCF) &
                = uCF(:,nX(1),iX2,iX3,iCF)

              uCF(:,nX(1)+1,iX2,iX3,iCF) &
                = uCF(:,1,iX2,iX3,iCF)

            END DO

            DO iPF = 1, nPF

              uPF(:,0,iX2,iX3,iPF) &
                = uPF(:,nX(1),iX2,iX3,iPF)

              uPF(:,nX(1)+1,iX2,iX3,iPF) &
                = uPF(:,1,iX2,iX3,iPF)

            END DO

            DO iAF = 1, nAF

              uAF(:,0,iX2,iX3,iAF) &
                = uAF(:,nX(1),iX2,iX3,iAF)

              uAF(:,nX(1)+1,iX2,iX3,iAF) &
                = uAF(:,1,iX2,iX3,iAF)

            END DO

          END DO
        END DO

      CASE ( 2 ) ! Homogeneous

        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)

            DO iCF = 1, nCF

              uCF(:,0,iX2,iX3,iCF) &
                = uCF(:,1,iX2,iX3,iCF)

              uCF(:,nX(1)+1,iX2,iX3,iCF) &
                = uCF(:,nX(1),iX2,iX3,iCF)

            END DO

            DO iPF = 1, nPF

              uPF(:,0,iX2,iX3,iPF) &
                = uPF(:,1,iX2,iX3,iPF)

              uPF(:,nX(1)+1,iX2,iX3,iPF) &
                = uPF(:,nX(1),iX2,iX3,iPF)

            END DO

            DO iAF = 1, nAF

              uAF(:,0,iX2,iX3,iAF) &
                = uAF(:,1,iX2,iX3,iAF)

              uAF(:,nX(1)+1,iX2,iX3,iAF) &
                = uAF(:,nX(1),iX2,iX3,iAF)

            END DO

          END DO
        END DO

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A41,I2.2)') &
          '', 'Invalid Boundary Condition for Fluid X1: ', bcX(1)
        STOP

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_Fluid_X1


  SUBROUTINE ApplyBoundaryConditions_Radiation( Time )

    REAL(DP), INTENT(in) :: Time

    IF( bcX(1) == 10 )THEN

      CALL ApplyApplicationBoundaryConditions_Radiation_X1( Time )

    ELSE

      CALL ApplyBoundaryConditions_Radiation_X1

    END IF

  END SUBROUTINE ApplyBoundaryConditions_Radiation


  SUBROUTINE ApplyBoundaryConditions_Radiation_X1

    INTEGER :: iS, iX2, iX3, iE, iCR

    SELECT CASE ( bcX(1) )

      CASE ( 0 ) ! No Boundary Condition

      CASE ( 1 ) ! Periodic

        DO iS = 1, nSpecies
          DO iX3 = 1, nX(3)
            DO iX2 = 1, nX(2)
              DO iE = 1, nE

                DO iCR = 1, nCR

                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)

                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,1,iX2,iX3,iCR,iS)

                END DO

              END DO
            END DO
          END DO
        END DO

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A45,I2.2)') &
          '', 'Invalid Boundary Condition for Radiation X1: ', bcX(1)
        STOP

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_Radiation_X1


END MODULE BoundaryConditionsModule
