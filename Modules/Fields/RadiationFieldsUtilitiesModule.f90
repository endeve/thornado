MODULE RadiationFieldsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CellAverage

CONTAINS


  PURE REAL(DP) FUNCTION CellAverage( uR, VJ )

    REAL(DP), DIMENSION(1:nDOF), INTENT(in) :: uR, VJ

!!$    CellAverage &
!!$      = DOT_PRODUCT( WeightsR, uR * VJ ) &
!!$        / DOT_PRODUCT( WeightsR, VJ )

    RETURN
  END FUNCTION CellAverage


END MODULE RadiationFieldsUtilitiesModule
