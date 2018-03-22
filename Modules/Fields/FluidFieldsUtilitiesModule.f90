MODULE FluidFieldsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  IMPLICIT NONE
  PRIVATE

CONTAINS


  PURE REAL(DP) FUNCTION CellAverage( uF, VJ )

    REAL(DP), DIMENSION(1:nDOFX), INTENT(in) :: uF, VJ

!!$    CellAverage &
!!$      = DOT_PRODUCT( WeightsF, uF * VJ ) &
!!$        / DOT_PRODUCT( WeightsF, VJ )

    RETURN
  END FUNCTION CellAverage


END MODULE FluidFieldsUtilitiesModule
