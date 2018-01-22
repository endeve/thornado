MODULE ReferenceElementModule_Lagrange

  USE KindModule, ONLY: &
    DP, Half
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
    nNodesX, &
    nDOF
  USE UtilitiesModule, ONLY: &
    WriteMatrix
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_E,  L_X1,  L_X2,  L_X3, &
         dL_X1, dL_X2, dL_X3
  USE ReferenceElementModule_Beta, ONLY: &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    NodeNumberTable, &
    NodeNumberTable_X1, &
    NodeNumberTable_X2, &
    NodeNumberTable_X3, &
    NodesE, &
    NodesX1, &
    NodesX2, &
    NodesX3

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  REAL(DP) :: wTime
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC :: i_L_X1
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC :: j_L_X1
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: L_X1_Dn_CRS
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X1_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X1_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X2_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X2_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X3_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X3_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdX1_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdX2_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdX3_q

  PUBLIC :: InitializeReferenceElement_Lagrange
  PUBLIC :: FinalizeReferenceElement_Lagrange

CONTAINS


  SUBROUTINE InitializeReferenceElement_Lagrange

    INTEGER :: k
    INTEGER :: iNode, iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER :: jNode, jNodeE, jNodeX1, jNodeX2, jNodeX3
    INTEGER :: iNode_X1, iNodeE_X1, iNodeX2_X1, iNodeX3_X1
    INTEGER :: iNode_X2, iNodeE_X2, iNodeX1_X2, iNodeX3_X2
    INTEGER :: iNode_X3, iNodeE_X3, iNodeX1_X3, iNodeX2_X3

    wTime = MPI_WTIME( )

    ALLOCATE( i_L_X1(nDOF) )
    ALLOCATE( j_L_X1(nDOF) )
    ALLOCATE( L_X1_Dn_CRS(nDOF) )

    ALLOCATE( L_X1_Dn(nDOF_X1,nDOF) )
    ALLOCATE( L_X1_Up(nDOF_X1,nDOF) )

    ALLOCATE( L_X2_Dn(nDOF_X2,nDOF) )
    ALLOCATE( L_X2_Up(nDOF_X2,nDOF) )

    ALLOCATE( L_X3_Dn(nDOF_X3,nDOF) )
    ALLOCATE( L_X3_Up(nDOF_X3,nDOF) )

    k = 1
    !$OMP PARALLEL DO PRIVATE &
    !$OMP&              ( iNode,iNodeE,iNodeX1,iNodeX2,iNodeX3, &
    !$OMP&                iNode_X1,iNodeE_X1,iNodeX2_X1,iNodeX3_X1 )
    DO iNode_X1 = 1, nDOF_X1

      iNodeE_X1  = NodeNumberTable_X1(1,iNode_X1)
      iNodeX2_X1 = NodeNumberTable_X1(2,iNode_X1)
      iNodeX3_X1 = NodeNumberTable_X1(3,iNode_X1)

      DO iNode = 1, nDOF

        iNodeE  = NodeNumberTable(1,iNode)
        iNodeX1 = NodeNumberTable(2,iNode)
        iNodeX2 = NodeNumberTable(3,iNode)
        iNodeX3 = NodeNumberTable(4,iNode)

        L_X1_Dn(iNode_X1,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X1) ) &
            * L_X1(iNodeX1) % P( - Half ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X1) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X1) )

        IF( iNodeE  == iNodeE_X1  .AND. &
            iNodeX2 == iNodeX2_X1 .AND. &
            iNodeX3 == iNodeX3_X1 ) &
        THEN

          L_X1_Dn_CRS(k) &
            = L_E(iNodeE) % P( NodesE(iNodeE_X1) ) &
              * L_X1(iNodeX1) % P( - Half ) &
              * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X1) ) &
              * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X1) )

          i_L_X1(k) = iNode_X1
          j_L_X1(k) = iNode

          k = k + 1

        END IF

        L_X1_Up(iNode_X1,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X1) ) &
            * L_X1(iNodeX1) % P( + Half ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X1) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X1) )

      END DO

    END DO
    !$OMP END PARALLEL DO

!    PRINT*,"MAX/MIN L_X1_CRS = ", &
!      MAXVAL( ABS( L_X1_Dn_CRS ) ), MINVAL( ABS( L_X1_Dn_CRS ) )

    k = 1
    !$OMP PARALLEL DO PRIVATE &
    !$OMP&              ( iNode,iNodeE,iNodeX1,iNodeX2,iNodeX3,  &
    !$OMP&                iNode_X2,iNodeE_X2,iNodeX1_X2,iNodeX3_X2, &
    !$OMP&                iNode_X3,iNodeE_X3,iNodeX1_X3,iNodeX2_X3 )
    DO iNode = 1, nDOF

      iNodeE  = NodeNumberTable(1,iNode)
      iNodeX1 = NodeNumberTable(2,iNode)
      iNodeX2 = NodeNumberTable(3,iNode)
      iNodeX3 = NodeNumberTable(4,iNode)

      DO iNode_X2 = 1, nDOF_X2

        iNodeE_X2  = NodeNumberTable_X2(1,iNode_X2)
        iNodeX1_X2 = NodeNumberTable_X2(2,iNode_X2)
        iNodeX3_X2 = NodeNumberTable_X2(3,iNode_X2)

        L_X2_Dn(iNode_X2,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X2) ) &
            * L_X1(iNodeX1) % P( NodesX1(iNodeX1_X2) ) &
            * L_X2(iNodeX2) % P( - Half ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X2) )

        L_X2_Up(iNode_X2,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X2) ) &
            * L_X1(iNodeX1) % P( NodesX1(iNodeX1_X2) ) &
            * L_X2(iNodeX2) % P( + Half ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X2) )

      END DO

      DO iNode_X3 = 1, nDOF_X3

        iNodeE_X3  = NodeNumberTable_X3(1,iNode_X3)
        iNodeX1_X3 = NodeNumberTable_X3(2,iNode_X3)
        iNodeX2_X3 = NodeNumberTable_X3(3,iNode_X3)

        L_X3_Dn(iNode_X3,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X3) ) &
            * L_X1(iNodeX1) % P( NodesX1(iNodeX1_X3) ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X3) ) &
            * L_X3(iNodeX3) % P( - Half )

        L_X3_Up(iNode_X3,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X3) ) &
            * L_X1(iNodeX1) % P( NodesX1(iNodeX1_X3) ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X3) ) &
            * L_X3(iNodeX3) % P( + Half )

      END DO

    END DO
    !$OMP END PARALLEL DO

    ALLOCATE( dLdX1_q(nDOF,nDOF) )
    ALLOCATE( dLdX2_q(nDOF,nDOF) )
    ALLOCATE( dLdX3_q(nDOF,nDOF) )

    !$OMP PARALLEL DO PRIVATE( iNode, iNodeE, iNodeX1, iNodeX2, iNodeX3, &
    !$OMP&                     jNode, jNodeE, jNodeX1, jNodeX2, jNodeX3 )
    DO jNode = 1, nDOF

      jNodeE  = NodeNumberTable(1,jNode)
      jNodeX1 = NodeNumberTable(2,jNode)
      jNodeX2 = NodeNumberTable(3,jNode)
      jNodeX3 = NodeNumberTable(4,jNode)

      DO iNode = 1, nDOF

        iNodeE  = NodeNumberTable(1,iNode)
        iNodeX1 = NodeNumberTable(2,iNode)
        iNodeX2 = NodeNumberTable(3,iNode)
        iNodeX3 = NodeNumberTable(4,iNode)

        dLdX1_q(iNode,jNode) &
          = L_E(jNodeE) % P( NodesE(iNodeE) ) &
            * dL_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            *  L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            *  L_X3(jNodeX3) % P( NodesX3(iNodeX3) )

        dLdX2_q(iNode,jNode) &
          = L_E(jNodeE) % P( NodesE(iNodeE) ) &
            *  L_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            * dL_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            *  L_X3(jNodeX3) % P( NodesX3(iNodeX3) )

        dLdX3_q(iNode,jNode) &
          = L_E(jNodeE) % P( NodesE(iNodeE) ) &
            *  L_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            *  L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            * dL_X3(jNodeX3) % P( NodesX3(iNodeX3) )

      END DO

    END DO
    !$OMP END PARALLEL DO

    wTime = MPI_WTIME( ) - wTime

    WRITE(*,*)
    WRITE(*,'(A4,A,ES10.4E2)') &
      '', 'InitializeReferenceElement_Lagrange: ', wTime
    WRITE(*,*)

    CALL WriteMatrix( nDOF_X1, nDOF, L_X1_Dn, 'L_X1_Dn.dat' )
    CALL WriteMatrix( nDOF_X1, nDOF, L_X1_Up, 'L_X1_Up.dat' )
    CALL WriteMatrix( nDOF_X2, nDOF, L_X2_Dn, 'L_X2_Dn.dat' )
    CALL WriteMatrix( nDOF_X2, nDOF, L_X2_Up, 'L_X2_Up.dat' )
    CALL WriteMatrix( nDOF_X3, nDOF, L_X3_Dn, 'L_X3_Dn.dat' )
    CALL WriteMatrix( nDOF_X3, nDOF, L_X3_Up, 'L_X3_Up.dat' )

    CALL WriteMatrix( nDOF, nDOF, dLdX1_q, 'dLdX1_q.dat' )
    CALL WriteMatrix( nDOF, nDOF, dLdX2_q, 'dLdX2_q.dat' )
    CALL WriteMatrix( nDOF, nDOF, dLdX3_q, 'dLdX3_q.dat' )

  END SUBROUTINE InitializeReferenceElement_Lagrange


  SUBROUTINE FinalizeReferenceElement_Lagrange

    DEALLOCATE( i_L_X1 )
    DEALLOCATE( j_L_X1 )
    DEALLOCATE( L_X1_Dn_CRS )

    DEALLOCATE( L_X1_Dn )
    DEALLOCATE( L_X1_Up )
    DEALLOCATE( L_X2_Dn )
    DEALLOCATE( L_X2_Up )
    DEALLOCATE( L_X3_Dn )
    DEALLOCATE( L_X3_Up )

    DEALLOCATE( dLdX1_q )
    DEALLOCATE( dLdX2_q )
    DEALLOCATE( dLdX3_q )

  END SUBROUTINE FinalizeReferenceElement_Lagrange


END MODULE ReferenceElementModule_Lagrange
