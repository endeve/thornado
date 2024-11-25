MODULE ReferenceElementModule

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nNodesE, &
    nDOF

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nDOF_X1
  INTEGER :: nDOF_X2
  INTEGER :: nDOF_X3
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: NodeNumberTable1D3D
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: NodeNumberTable1D3D_X1
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: NodeNumberTable1D3D_X2
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: NodeNumberTable1D3D_X3
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesE,  WeightsE
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesX1, WeightsX1
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesX2, WeightsX2
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesX3, WeightsX3
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Weights_q
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Weights_X1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Nodes_q

  PUBLIC :: InitializeReferenceElement
  PUBLIC :: FinalizeReferenceElement

CONTAINS


  SUBROUTINE InitializeReferenceElement

    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode

    nDOF_X1 = nNodesX(2) * nNodesX(3) * nNodesE
    nDOF_X2 = nNodesX(1) * nNodesX(3) * nNodesE
    nDOF_X3 = nNodesX(1) * nNodesX(2) * nNodesE

    ALLOCATE( NodeNumberTable1D3D(4,nDOF) )

    iNode = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            NodeNumberTable1D3D(1:4,iNode) &
              = [ iNodeE, iNodeX1, iNodeX2, iNodeX3 ]

            iNode = iNode + 1

          END DO
        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTable1D3D_X1(3,nDOF_X1) )

    iNode = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeE = 1, nNodesE

          NodeNumberTable1D3D_X1(1:3,iNode) &
            = [ iNodeE, iNodeX2, iNodeX3 ]

          iNode = iNode + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTable1D3D_X2(3,nDOF_X2) )

    iNode = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX1 = 1, nNodesX(1)
        DO iNodeE = 1, nNodesE

          NodeNumberTable1D3D_X2(1:3,iNode) &
            = [ iNodeE, iNodeX1, iNodeX3 ]

          iNode = iNode + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTable1D3D_X3(3,nDOF_X3) )

    iNode = 1
    DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)
        DO iNodeE = 1, nNodesE

          NodeNumberTable1D3D_X3(1:3,iNode) &
            = [ iNodeE, iNodeX1, iNodeX2 ]

          iNode = iNode + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodesE (nNodesE),    WeightsE (nNodesE) )
    ALLOCATE( NodesX1(nNodesX(1)), WeightsX1(nNodesX(1)) )
    ALLOCATE( NodesX2(nNodesX(2)), WeightsX2(nNodesX(2)) )
    ALLOCATE( NodesX3(nNodesX(3)), WeightsX3(nNodesX(3)) )

    CALL GetQuadrature( nNodesE,    NodesE,  WeightsE )
    CALL GetQuadrature( nNodesX(1), NodesX1, WeightsX1 )
    CALL GetQuadrature( nNodesX(2), NodesX2, WeightsX2 )
    CALL GetQuadrature( nNodesX(3), NodesX3, WeightsX3 )

    ALLOCATE( Weights_q(nDOF) )
    ALLOCATE( Nodes_q(4,nDOF) )

    DO iNode = 1, nDOF

      iNodeE  = NodeNumberTable1D3D(1,iNode)
      iNodeX1 = NodeNumberTable1D3D(2,iNode)
      iNodeX2 = NodeNumberTable1D3D(3,iNode)
      iNodeX3 = NodeNumberTable1D3D(4,iNode)

      Weights_q(iNode) &
        = WeightsE(iNodeE) * WeightsX1(iNodeX1) &
            * WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

      Nodes_q(1:4,iNode) &
        = [ NodesE (iNodeE),  NodesX1(iNodeX1), &
            NodesX2(iNodeX2), NodesX3(iNodeX3) ]

    END DO

    ALLOCATE( Weights_X1(nDOF_X1) )

    DO iNode = 1, nDOF_X1

      iNodeE  = NodeNumberTable1D3D_X1(1,iNode)
      iNodeX2 = NodeNumberTable1D3D_X1(2,iNode)
      iNodeX3 = NodeNumberTable1D3D_X1(3,iNode)

      Weights_x1(iNode) &
        = WeightsE(iNodeE) * WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

    END DO

  END SUBROUTINE InitializeReferenceElement


  SUBROUTINE FinalizeReferenceElement

    DEALLOCATE( NodeNumberTable1D3D )
    DEALLOCATE( NodeNumberTable1D3D_X1 )
    DEALLOCATE( NodeNumberTable1D3D_X2 )
    DEALLOCATE( NodeNumberTable1D3D_X3 )
    DEALLOCATE( NodesE,  WeightsE )
    DEALLOCATE( NodesX1, WeightsX1 )
    DEALLOCATE( NodesX2, WeightsX2 )
    DEALLOCATE( NodesX3, WeightsX3 )
    DEALLOCATE( Nodes_q, Weights_q )
    DEALLOCATE( Weights_X1 )

  END SUBROUTINE FinalizeReferenceElement


END MODULE ReferenceElementModule
