MODULE RadiationFieldsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOF

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Weights

  PUBLIC :: InitializeRadiationFieldsUtilities
  PUBLIC :: FinalizeRadiationFieldsUtilities
  PUBLIC :: CellAverage

CONTAINS


  SUBROUTINE InitializeRadiationFieldsUtilities( w_E, w_X1, w_X2, w_X3 )

    REAL(DP), DIMENSION(:), INTENT(in) :: w_E, w_X1, w_X2, w_X3

    ALLOCATE( Weights(nDOF) )

  END SUBROUTINE InitializeRadiationFieldsUtilities


  SUBROUTINE FinalizeRadiationFieldsUtilities

    DEALLOCATE( Weights )

  END SUBROUTINE FinalizeRadiationFieldsUtilities


  PURE REAL(DP) FUNCTION CellAverage( uR, VJ )

    REAL(DP), DIMENSION(1:nDOF), INTENT(in) :: uR, VJ

    CellAverage &
      = DOT_PRODUCT( Weights, uR * VJ ) &
        / DOT_PRODUCT( Weights, VJ )

  END FUNCTION CellAverage


END MODULE RadiationFieldsUtilitiesModule
