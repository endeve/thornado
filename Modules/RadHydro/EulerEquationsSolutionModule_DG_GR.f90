MODULE EulerEquationsSolutionModule_DG_GR

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, &
    uAF, nAF
  USE EulerEquationsUtilitiesModule_GR, ONLY: &
    ComputePrimitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeRHS_Euler_DG_GR

CONTAINS


  SUBROUTINE ComputeRHS_Euler_DG_GR( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3

    PRINT*, "ComputeRHS_Euler_DG_GR Begin"

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uGF(:,iX1,iX2,iX3,1:nGF), &
                   uPF(:,iX1,iX2,iX3,1:nPF), uAF(:,iX1,iX2,iX3,1:nAF) )

        END DO
      END DO
    END DO

    PRINT*, "ComputeRHS_Euler_DG_GR End"

  END SUBROUTINE ComputeRHS_Euler_DG_GR


END MODULE EulerEquationsSolutionModule_DG_GR
