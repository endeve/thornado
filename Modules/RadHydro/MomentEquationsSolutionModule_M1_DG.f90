MODULE MomentEquationsSolutionModule_M1_DG

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeRHS_M1_DG

CONTAINS


  SUBROUTINE ComputeRHS_M1_DG( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

  END SUBROUTINE ComputeRHS_M1_DG


END MODULE MomentEquationsSolutionModule_M1_DG
