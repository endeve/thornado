MODULE BoundaryConditionsModule

  USE ProgramHeaderModule, ONLY: &
    nX, bcX
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, uPF, nPF, uAF, nAF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Fluid

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

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_Fluid_X1


END MODULE BoundaryConditionsModule
