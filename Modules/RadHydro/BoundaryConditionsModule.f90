MODULE BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, bcX, nE
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, iCF_S1, &
    uPF, nPF, iPF_V1, &
    uAF, nAF
  USE RadiationFieldsModule, ONLY: &
    uCR, nCR, nSpecies
  USE ApplicationBoundaryConditionsModule, ONLY: &
    ApplyApplicationBoundaryConditions_Fluid_X1, &
    ApplyApplicationBoundaryConditions_Radiation_X1

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Fluid
  PUBLIC :: ApplyBoundaryConditions_Radiation

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Fluid

    IF( bcX(1) == 10 )THEN

      CALL ApplyApplicationBoundaryConditions_Fluid_X1

    ELSE

      CALL ApplyBoundaryConditions_Fluid_X1

    END IF

  END SUBROUTINE ApplyBoundaryConditions_Fluid


  SUBROUTINE ApplyBoundaryConditions_Fluid_X1

    INTEGER :: iX2, iX3, iCF, iPF, iAF
    INTEGER :: iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNodeX, jNodeX

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

      CASE ( 3 ) ! Reflecting

        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
                  jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

                  ! -- Conserved --

                  DO iCF = 1, nCF

                    uCF(iNodeX,0,      iX2,iX3,iCF) &
                      = uCF(jNodeX,1,    iX2,iX3,iCF)
                    uCF(iNodeX,nX(1)+1,iX2,iX3,iCF) &
                      = uCF(jNodeX,nX(1),iX2,iX3,iCF)

                  END DO

                  uCF(iNodeX,0,      iX2,iX3,iCF_S1) &
                    = - uCF(jNodeX,1,    iX2,iX3,iCF_S1)
                  uCF(iNodeX,nX(1)+1,iX2,iX3,iCF_S1) &
                    = - uCF(jNodeX,nX(1),iX2,iX3,iCF_S1)

                  ! -- Primitive --

                  DO iPF = 1, nPF

                    uPF(iNodeX,0,      iX2,iX3,iPF) &
                      = uPF(jNodeX,1,    iX2,iX3,iPF)
                    uPF(iNodeX,nX(1)+1,iX2,iX3,iPF) &
                      = uPF(jNodeX,nX(1),iX2,iX3,iPF)

                  END DO

                  uPF(iNodeX,0,      iX2,iX3,iPF_V1) &
                    = - uPF(jNodeX,1,    iX2,iX3,iPF_V1)
                  uPF(iNodeX,nX(1)+1,iX2,iX3,iPF_V1) &
                    = - uPF(jNodeX,nX(1),iX2,iX3,iPF_V1)

                  ! -- Auxiliary --

                  DO iAF = 1, nAF

                    uAF(iNodeX,0,      iX2,iX3,iAF) &
                      = uAF(jNodeX,1,    iX2,iX3,iAF)
                    uAF(iNodeX,nX(1)+1,iX2,iX3,iAF) &
                      = uAF(jNodeX,nX(1),iX2,iX3,iAF)

                  END DO

                END DO
              END DO
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


  SUBROUTINE ApplyBoundaryConditions_Radiation( Time, LimiterBC_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: LimiterBC_Option

    IF( bcX(1) == 10 )THEN

      CALL ApplyApplicationBoundaryConditions_Radiation_X1 &
             ( Time, LimiterBC_Option )

    ELSE

      CALL ApplyBoundaryConditions_Radiation_X1

      CALL ApplyBoundaryConditions_Radiation_X2

      CALL ApplyBoundaryConditions_Radiation_X3

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

      CASE ( 2 ) ! Homogeneous

        DO iS = 1, nSpecies
          DO iX3 = 1, nX(3)
            DO iX2 = 1, nX(2)
              DO iE = 1, nE

                DO iCR = 1, nCR

                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,1,iX2,iX3,iCR,iS)

                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)

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


  SUBROUTINE ApplyBoundaryConditions_Radiation_X2

    INTEGER :: iS, iX2, iX3, iE, iCR

    SELECT CASE ( bcX(2) )

      CASE ( 0 ) ! No Boundary Condition

      CASE ( 1 ) ! Periodic

!!$        DO iS = 1, nSpecies
!!$          DO iX3 = 1, nX(3)
!!$            DO iX2 = 1, nX(2)
!!$              DO iE = 1, nE
!!$
!!$                DO iCR = 1, nCR
!!$
!!$                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)
!!$
!!$                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,1,iX2,iX3,iCR,iS)
!!$
!!$                END DO
!!$
!!$              END DO
!!$            END DO
!!$          END DO
!!$        END DO

      CASE ( 2 ) ! Homogeneous

!!$        DO iS = 1, nSpecies
!!$          DO iX3 = 1, nX(3)
!!$            DO iX2 = 1, nX(2)
!!$              DO iE = 1, nE
!!$
!!$                DO iCR = 1, nCR
!!$
!!$                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,1,iX2,iX3,iCR,iS)
!!$
!!$                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)
!!$
!!$                END DO
!!$
!!$              END DO
!!$            END DO
!!$          END DO
!!$        END DO

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A45,I2.2)') &
          '', 'Invalid Boundary Condition for Radiation X2: ', bcX(2)
        STOP

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_Radiation_X2


  SUBROUTINE ApplyBoundaryConditions_Radiation_X3

    INTEGER :: iS, iX2, iX3, iE, iCR

    SELECT CASE ( bcX(3) )

      CASE ( 0 ) ! No Boundary Condition

      CASE ( 1 ) ! Periodic

!!$        DO iS = 1, nSpecies
!!$          DO iX3 = 1, nX(3)
!!$            DO iX2 = 1, nX(2)
!!$              DO iE = 1, nE
!!$
!!$                DO iCR = 1, nCR
!!$
!!$                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)
!!$
!!$                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,1,iX2,iX3,iCR,iS)
!!$
!!$                END DO
!!$
!!$              END DO
!!$            END DO
!!$          END DO
!!$        END DO

      CASE ( 2 ) ! Homogeneous

!!$        DO iS = 1, nSpecies
!!$          DO iX3 = 1, nX(3)
!!$            DO iX2 = 1, nX(2)
!!$              DO iE = 1, nE
!!$
!!$                DO iCR = 1, nCR
!!$
!!$                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,1,iX2,iX3,iCR,iS)
!!$
!!$                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
!!$                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)
!!$
!!$                END DO
!!$
!!$              END DO
!!$            END DO
!!$          END DO
!!$        END DO

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A45,I2.2)') &
          '', 'Invalid Boundary Condition for Radiation X3: ', bcX(3)
        STOP

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_Radiation_X3


END MODULE BoundaryConditionsModule
