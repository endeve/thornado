MODULE ReferenceElementModuleX

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX

  IMPLICIT NONE
  PRIVATE

  INTEGER,                               PUBLIC :: nDOFX_X1
  INTEGER,                               PUBLIC :: nDOFX_X2
  INTEGER,                               PUBLIC :: nDOFX_X3
  INTEGER,  DIMENSION(:,:), ALLOCATABLE, PUBLIC :: NodeNumberTableX
  INTEGER,  DIMENSION(:,:), ALLOCATABLE, PUBLIC :: NodeNumberTableX_X1
  INTEGER,  DIMENSION(:,:), ALLOCATABLE, PUBLIC :: NodeNumberTableX_X2
  INTEGER,  DIMENSION(:,:), ALLOCATABLE, PUBLIC :: NodeNumberTableX_X3
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: NodesX1, WeightsX1
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: NodesX2, WeightsX2
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: NodesX3, WeightsX3
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: WeightsX_X1
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: WeightsX_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: NodesX_q

  PUBLIC :: InitializeReferenceElementX
  PUBLIC :: FinalizeReferenceElementX

CONTAINS


  SUBROUTINE InitializeReferenceElementX

    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX

    nDOFX_X1 = nNodesX(2) * nNodesX(3)
    nDOFX_X2 = nNodesX(1) * nNodesX(3)
    nDOFX_X3 = nNodesX(1) * nNodesX(2)

    ALLOCATE( NodeNumberTableX(3,nDOFX) )

    iNodeX = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          NodeNumberTableX(1:3,iNodeX) &
            = [ iNodeX1, iNodeX2, iNodeX3 ]

          iNodeX = iNodeX + 1

        END DO
      END DO
    END DO

    ALLOCATE( NodeNumberTableX_X1(2,nDOFX_X1) )

    iNodeX = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)

        NodeNumberTableX_X1(1:2,iNodeX) &
          = [ iNodeX2, iNodeX3 ]

        iNodeX = iNodeX + 1

      END DO
    END DO

    ALLOCATE( NodeNumberTableX_X2(2,nDOFX_X2) )

    iNodeX = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX1 = 1, nNodesX(1)

        NodeNumberTableX_X2(1:2,iNodeX) &
          = [ iNodeX1, iNodeX3 ]

        iNodeX = iNodeX + 1

      END DO
    END DO

    ALLOCATE( NodeNumberTableX_X3(2,nDOFX_X3) )

    iNodeX = 1
    DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        NodeNumberTableX_X3(1:2,iNodeX) &
          = [ iNodeX1, iNodeX2 ]

        iNodeX = iNodeX + 1

      END DO
    END DO

    ALLOCATE( NodesX1(nNodesX(1)), WeightsX1(nNodesX(1)) )
    ALLOCATE( NodesX2(nNodesX(2)), WeightsX2(nNodesX(2)) )
    ALLOCATE( NodesX3(nNodesX(3)), WeightsX3(nNodesX(3)) )

    CALL GetQuadrature( nNodesX(1), NodesX1, WeightsX1 )
    CALL GetQuadrature( nNodesX(2), NodesX2, WeightsX2 )
    CALL GetQuadrature( nNodesX(3), NodesX3, WeightsX3 )

    ALLOCATE( WeightsX_X1(nDOFX_X1) )

    DO iNodeX = 1, nDOFX_X1

      iNodeX2 = NodeNumberTableX_X1(1,iNodeX)
      iNodeX3 = NodeNumberTableX_X1(2,iNodeX)

      WeightsX_X1(iNodeX) = WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

    END DO

    ALLOCATE( WeightsX_q(nDOFX) )
    ALLOCATE( NodesX_q(3,nDOFX) )

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      WeightsX_q(iNodeX) &
        = WeightsX1(iNodeX1) * WeightsX2(iNodeX2) * WeightsX3(iNodeX3)

      NodesX_q(1:3,iNodeX) &
        = [ NodesX1(iNodeX1), NodesX2(iNodeX2), NodesX3(iNodeX3) ]

    END DO

  END SUBROUTINE InitializeReferenceElementX


  SUBROUTINE FinalizeReferenceElementX

    DEALLOCATE( NodeNumberTableX )
    DEALLOCATE( NodeNumberTableX_X1 )
    DEALLOCATE( NodeNumberTableX_X2 )
    DEALLOCATE( NodeNumberTableX_X3 )
    DEALLOCATE( NodesX1, WeightsX1 )
    DEALLOCATE( NodesX2, WeightsX2 )
    DEALLOCATE( NodesX3, WeightsX3 )
    DEALLOCATE( WeightsX_X1 )
    DEALLOCATE( WeightsX_q )
    DEALLOCATE( NodesX_q )

  END SUBROUTINE FinalizeReferenceElementX


END MODULE ReferenceElementModuleX
