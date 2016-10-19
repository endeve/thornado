MODULE RadiationFieldsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nNodesE, &
    nDOF
  USE UtilitiesModule, ONLY: &
    NodeNumber

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

    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNode

    ALLOCATE( Weights(nDOF) )

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE  = 1, nNodesE

            iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

            Weights(iNode) &
              = w_E(iNodeE) * w_X1(iNodeX1) * w_X2(iNodeX2) * w_X3(iNodeX3)

          END DO
        END DO
      END DO
    END DO

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
