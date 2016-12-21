MODULE GravitySolutionUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    VolJacX
  USE FluidFieldsModule, ONLY: &
    WeightsF, &
    uPF, iPF_D

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTotalBaryonMass

CONTAINS


  SUBROUTINE ComputeTotalBaryonMass( Mass )

    REAL(DP), INTENT(out) :: Mass

    INTEGER :: iX1, iX2, iX3

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    Mass = 0.0_DP
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          PRINT*
          PRINT*, "dX = ", dX1(iX1), dX2(iX2), dX3(iX3)
          PRINT*, "VJ = ", VolJacX(:,iX1,iX2,iX3)
          PRINT*, "D  = ", uPF(:,iX1,iX2,iX3,iPF_D)
          PRINT*, "W  = ", WeightsF, SUM( WeightsF )

          Mass &
            = Mass &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsF(:) * uPF(:,iX1,iX2,iX3,iPF_D) &
                             * VolJacX(:,iX1,iX2,iX3) )

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTotalBaryonMass


END MODULE GravitySolutionUtilitiesModule
