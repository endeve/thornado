MODULE ReferenceElementModuleZ

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    GetQuadrature
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOFZ

  IMPLICIT NONE
  PRIVATE

  INTEGER,               PUBLIC :: nDOFZ_Z1
  INTEGER,               PUBLIC :: nDOFZ_Z2
  INTEGER,               PUBLIC :: nDOFZ_Z3
  INTEGER,               PUBLIC :: nDOFZ_Z4
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableZ(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableZ_Z1(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableZ_Z2(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableZ_Z3(:,:)
  INTEGER,  ALLOCATABLE, PUBLIC :: NodeNumberTableZ_Z4(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesZ1(:), WeightsZ1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesZ2(:), WeightsZ2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesZ3(:), WeightsZ3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: NodesZ4(:), WeightsZ4(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsZ_q (:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsZ_Z1(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsZ_Z2(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsZ_Z3(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: WeightsZ_Z4(:)

  PUBLIC :: InitializeReferenceElementZ
  PUBLIC :: FinalizeReferenceElementZ

CONTAINS


  SUBROUTINE InitializeReferenceElementZ

    INTEGER :: iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4, iNodeZ

    nDOFZ_Z1 =              nNodesZ(2) * nNodesZ(3) * nNodesZ(4)
    nDOFZ_Z2 = nNodesZ(1)              * nNodesZ(3) * nNodesZ(4)
    nDOFZ_Z3 = nNodesZ(1) * nNodesZ(2)              * nNodesZ(4)
    nDOFZ_Z4 = nNodesZ(1) * nNodesZ(2) * nNodesZ(3)

    ALLOCATE( NodeNumberTableZ(4,nDOFZ) )

    iNodeZ = 0
    DO iNodeZ4 = 1, nNodesZ(4)
    DO iNodeZ3 = 1, nNodesZ(3)
    DO iNodeZ2 = 1, nNodesZ(2)
    DO iNodeZ1 = 1, nNodesZ(1)

      iNodeZ = iNodeZ + 1

      NodeNumberTableZ(1:4,iNodeZ) &
        = [ iNodeZ1, iNodeZ2, iNodeZ3, iNodeZ4 ]

    END DO
    END DO
    END DO
    END DO

    ALLOCATE( NodeNumberTableZ_Z1(3,nDOFZ_Z1) )

    iNodeZ = 0
    DO iNodeZ4 = 1, nNodesZ(4)
    DO iNodeZ3 = 1, nNodesZ(3)
    DO iNodeZ2 = 1, nNodesZ(2)

      iNodeZ = iNodeZ + 1

      NodeNumberTableZ_Z1(1:3,iNodeZ) &
        = [ iNodeZ2, iNodeZ3, iNodeZ4 ]

    END DO
    END DO
    END DO

    ALLOCATE( NodeNumberTableZ_Z2(3,nDOFZ_Z2) )

    iNodeZ = 0
    DO iNodeZ4 = 1, nNodesZ(4)
    DO iNodeZ3 = 1, nNodesZ(3)
    DO iNodeZ1 = 1, nNodesZ(1)

      iNodeZ = iNodeZ + 1

      NodeNumberTableZ_Z2(1:3,iNodeZ) &
        = [ iNodeZ1, iNodeZ3, iNodeZ4 ]

    END DO
    END DO
    END DO

    ALLOCATE( NodeNumberTableZ_Z3(3,nDOFZ_Z3) )

    iNodeZ = 0
    DO iNodeZ4 = 1, nNodesZ(4)
    DO iNodeZ2 = 1, nNodesZ(2)
    DO iNodeZ1 = 1, nNodesZ(1)

      iNodeZ = iNodeZ + 1

      NodeNumberTableZ_Z3(1:3,iNodeZ) &
        = [ iNodeZ1, iNodeZ2, iNodeZ4 ]

    END DO
    END DO
    END DO

    ALLOCATE( NodeNumberTableZ_Z4(3,nDOFZ_Z4) )

    iNodeZ = 0
    DO iNodeZ3 = 1, nNodesZ(3)
    DO iNodeZ2 = 1, nNodesZ(2)
    DO iNodeZ1 = 1, nNodesZ(1)

      iNodeZ = iNodeZ + 1

      NodeNumberTableZ_Z4(1:3,iNodeZ) &
        = [ iNodeZ1, iNodeZ2, iNodeZ3 ]

    END DO
    END DO
    END DO

    ALLOCATE( NodesZ1(nNodesZ(1)), WeightsZ1(nNodesZ(1)) )
    ALLOCATE( NodesZ2(nNodesZ(2)), WeightsZ2(nNodesZ(2)) )
    ALLOCATE( NodesZ3(nNodesZ(3)), WeightsZ3(nNodesZ(3)) )
    ALLOCATE( NodesZ4(nNodesZ(4)), WeightsZ4(nNodesZ(4)) )

    CALL GetQuadrature( nNodesZ(1), NodesZ1, WeightsZ1 )
    CALL GetQuadrature( nNodesZ(2), NodesZ2, WeightsZ2 )
    CALL GetQuadrature( nNodesZ(3), NodesZ3, WeightsZ3 )
    CALL GetQuadrature( nNodesZ(4), NodesZ4, WeightsZ4 )

    ALLOCATE( WeightsZ_q(nDOFZ) )

    DO iNodeZ = 1, nDOFZ

      iNodeZ1 = NodeNumberTableZ(1,iNodeZ)
      iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
      iNodeZ3 = NodeNumberTableZ(3,iNodeZ)
      iNodeZ4 = NodeNumberTableZ(4,iNodeZ)

      WeightsZ_q(iNodeZ) &
        = WeightsZ1(iNodeZ1) * WeightsZ2(iNodeZ2) &
            * WeightsZ3(iNodeZ3) * WeightsZ4(iNodeZ4)

    END DO

    ALLOCATE( WeightsZ_Z1(nDOFZ_Z1) )

    DO iNodeZ = 1, nDOFZ_Z1

      iNodeZ2 = NodeNumberTableZ_Z1(1,iNodeZ)
      iNodeZ3 = NodeNumberTableZ_Z1(2,iNodeZ)
      iNodeZ4 = NodeNumberTableZ_Z1(3,iNodeZ)

      WeightsZ_Z1(iNodeZ) &
        = WeightsZ2(iNodeZ2) * WeightsZ3(iNodeZ3) * WeightsZ4(iNodeZ4)

    END DO

    ALLOCATE( WeightsZ_Z2(nDOFZ_Z2) )

    DO iNodeZ = 1, nDOFZ_Z2

      iNodeZ1 = NodeNumberTableZ_Z2(1,iNodeZ)
      iNodeZ3 = NodeNumberTableZ_Z2(2,iNodeZ)
      iNodeZ4 = NodeNumberTableZ_Z2(3,iNodeZ)

      WeightsZ_Z2(iNodeZ) &
        = WeightsZ1(iNodeZ1) * WeightsZ3(iNodeZ3) * WeightsZ4(iNodeZ4)

    END DO

    ALLOCATE( WeightsZ_Z3(nDOFZ_Z3) )

    DO iNodeZ = 1, nDOFZ_Z3

      iNodeZ1 = NodeNumberTableZ_Z3(1,iNodeZ)
      iNodeZ2 = NodeNumberTableZ_Z3(2,iNodeZ)
      iNodeZ4 = NodeNumberTableZ_Z3(3,iNodeZ)

      WeightsZ_Z3(iNodeZ) &
        = WeightsZ1(iNodeZ1) * WeightsZ2(iNodeZ2) * WeightsZ4(iNodeZ4)

    END DO

    ALLOCATE( WeightsZ_Z4(nDOFZ_Z4) )

    DO iNodeZ = 1, nDOFZ_Z4

      iNodeZ1 = NodeNumberTableZ_Z4(1,iNodeZ)
      iNodeZ2 = NodeNumberTableZ_Z4(2,iNodeZ)
      iNodeZ3 = NodeNumberTableZ_Z4(3,iNodeZ)

      WeightsZ_Z4(iNodeZ) &
        = WeightsZ1(iNodeZ1) * WeightsZ2(iNodeZ2) * WeightsZ3(iNodeZ3)

    END DO

  END SUBROUTINE InitializeReferenceElementZ


  SUBROUTINE FinalizeReferenceElementZ

    DEALLOCATE( NodeNumberTableZ )
    DEALLOCATE( NodeNumberTableZ_Z1 )
    DEALLOCATE( NodeNumberTableZ_Z2 )
    DEALLOCATE( NodeNumberTableZ_Z3 )
    DEALLOCATE( NodeNumberTableZ_Z4 )
    DEALLOCATE( NodesZ1, WeightsZ1 )
    DEALLOCATE( NodesZ2, WeightsZ2 )
    DEALLOCATE( NodesZ3, WeightsZ3 )
    DEALLOCATE( NodesZ4, WeightsZ4 )
    DEALLOCATE( WeightsZ_q )
    DEALLOCATE( WeightsZ_Z1 )
    DEALLOCATE( WeightsZ_Z2 )
    DEALLOCATE( WeightsZ_Z3 )
    DEALLOCATE( WeightsZ_Z4 )

  END SUBROUTINE FinalizeReferenceElementZ


END MODULE ReferenceElementModuleZ
