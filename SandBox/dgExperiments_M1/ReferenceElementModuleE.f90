MODULE ReferenceElementModuleE

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nNodesE

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: NodesE,   WeightsE
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: NodesE_L, WeightsE_L

  PUBLIC :: InitializeReferenceElementE
  PUBLIC :: FinalizeReferenceElementE

CONTAINS


  SUBROUTINE InitializeReferenceElementE

    ! --- Gaussian Quadrature Points and Weights ---

    ALLOCATE( NodesE(nNodesE), WeightsE(nNodesE) )

    CALL GetQuadrature( nNodesE, NodesE, WeightsE )

    ! --- Lobatto Quadrature Points and Weights ---

    ALLOCATE( NodesE_L(nNodesE), WeightsE_L(nNodesE) )

    CALL GetQuadrature( nNodesE, NodesE_L, WeightsE_L, 'Lobatto' )

  END SUBROUTINE InitializeReferenceElementE


  SUBROUTINE FinalizeReferenceElementE

    DEALLOCATE( NodesE,   WeightsE )
    DEALLOCATE( NodesE_L, WeightsE_L )

  END SUBROUTINE FinalizeReferenceElementE


END MODULE ReferenceElementModuleE
