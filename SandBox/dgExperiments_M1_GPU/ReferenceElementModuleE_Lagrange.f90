MODULE ReferenceElementModuleE_Lagrange

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
    nDOFE
  USE PolynomialBasisModule_Lagrange, ONLY: &
    LagrangeP
  USE ReferenceElementModuleE, ONLY: &
    NodesE, &
    NodesE_L

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LE_L2G

  PUBLIC :: InitializeReferenceElementE_Lagrange
  PUBLIC :: FinalizeReferenceElementE_Lagrange

CONTAINS


  SUBROUTINE InitializeReferenceElementE_Lagrange

    INTEGER :: iNodeE, jNodeE

    ALLOCATE( LE_L2G(nDOFE,nDOFE) )

    DO jNodeE = 1, nDOFE

      DO iNodeE = 1, nDOFE

        LE_L2G(iNodeE,jNodeE) &
          = LagrangeP( NodesE(iNodeE), jNodeE, NodesE_L, nNodesE )

      END DO

    END DO

  END SUBROUTINE InitializeReferenceElementE_Lagrange


  SUBROUTINE FinalizeReferenceElementE_Lagrange

    DEALLOCATE( LE_L2G )

  END SUBROUTINE FinalizeReferenceElementE_Lagrange


END MODULE ReferenceElementModuleE_Lagrange
