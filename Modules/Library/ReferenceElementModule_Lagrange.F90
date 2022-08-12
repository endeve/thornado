MODULE ReferenceElementModule_Lagrange

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesE, nNodesX, nDOF
  USE ReferenceElementModule, ONLY: &
    nDOF_E, nDOF_X1, nDOF_X2, nDOF_X3, &
    NodeNumberTable, &
    NodeNumberTable_E, &
    NodeNumberTable_X1, &
    NodeNumberTable_X2, &
    NodeNumberTable_X3, &
    NodesE, NodesX1, NodesX2, NodesX3
  USE PolynomialBasisModule_Lagrange, ONLY: &
     L_E,  L_X1,  L_X2,  L_X3, &
    dL_E, dL_X1, dL_X2, dL_X3

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_E_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_E_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X1_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X1_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X2_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X2_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X3_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L_X3_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdE_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdX1_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdX2_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLdX3_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: InterpMat

  PUBLIC :: InitializeReferenceElement_Lagrange
  PUBLIC :: FinalizeReferenceElement_Lagrange

CONTAINS


  SUBROUTINE InitializeReferenceElement_Lagrange

    INTEGER :: iNode, iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER :: jNode, jNodeE, jNodeX1, jNodeX2, jNodeX3
    INTEGER :: iNode_E ,            iNodeX1_E , iNodeX2_E , iNodeX3_E
    INTEGER :: iNode_X1, iNodeE_X1            , iNodeX2_X1, iNodeX3_X1
    INTEGER :: iNode_X2, iNodeE_X2, iNodeX1_X2            , iNodeX3_X2
    INTEGER :: iNode_X3, iNodeE_X3, iNodeX1_X3, iNodeX2_X3
    INTEGER :: i, iOS, nDOF_K

    ALLOCATE( L_E_Dn(nDOF_E,nDOF) )
    ALLOCATE( L_E_Up(nDOF_E,nDOF) )

    ALLOCATE( L_X1_Dn(nDOF_X1,nDOF) )
    ALLOCATE( L_X1_Up(nDOF_X1,nDOF) )

    ALLOCATE( L_X2_Dn(nDOF_X2,nDOF) )
    ALLOCATE( L_X2_Up(nDOF_X2,nDOF) )

    ALLOCATE( L_X3_Dn(nDOF_X3,nDOF) )
    ALLOCATE( L_X3_Up(nDOF_X3,nDOF) )


    DO iNode_E = 1, nDOF_E

      iNodeX1_E = NodeNumberTable_E(1,iNode_E)
      iNodeX2_E = NodeNumberTable_E(2,iNode_E)
      iNodeX3_E = NodeNumberTable_E(3,iNode_E)

      DO iNode = 1, nDOF

        iNodeE  = NodeNumberTable(1,iNode)
        iNodeX1 = NodeNumberTable(2,iNode)
        iNodeX2 = NodeNumberTable(3,iNode)
        iNodeX3 = NodeNumberTable(4,iNode)

        L_E_Dn(iNode_E,iNode) &
          = L_E(iNodeE) % P( - Half ) &
            * L_X1(iNodeX1) % P( NodesX1(iNodeX1_E) ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_E) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_E) )

        L_E_Up(iNode_E,iNode) &
          = L_E(iNodeE) % P( + Half ) &
            * L_X1(iNodeX1) % P( NodesX1(iNodeX1_E) ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_E) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_E) )

      END DO

    END DO


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

        L_X1_Up(iNode_X1,iNode) &
          = L_E(iNodeE) % P( NodesE(iNodeE_X1) ) &
            * L_X1(iNodeX1) % P( + Half ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X1) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X1) )

      END DO

    END DO


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

    ALLOCATE( dLdE_q (nDOF,nDOF) )
    ALLOCATE( dLdX1_q(nDOF,nDOF) )
    ALLOCATE( dLdX2_q(nDOF,nDOF) )
    ALLOCATE( dLdX3_q(nDOF,nDOF) )

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

        dLdE_q(iNode,jNode) &
          = dL_E(jNodeE) % P( NodesE(iNodeE) ) &
            * L_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            * L_X3(jNodeX3) % P( NodesX3(iNodeX3) )

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

    ! --- Total Number of Points on Elements (Interior and Faces) ---

    nDOF_K = nNodesE * PRODUCT( nNodesX )

    IF( nNodesE > 1 )THEN

      nDOF_K = nDOF_K + 2 * PRODUCT( nNodesX )

    END IF

    DO i = 1, 3

      IF( nNodesX(i) > 1 )THEN

        nDOF_K &
          = nDOF_K + 2 * nNodesE * PRODUCT( nNodesX, MASK = [1,2,3] .NE. i )

      END IF

    END DO

    ALLOCATE( InterpMat(nDOF_K,nDOF) )

    InterpMat = Zero

    DO jNode = 1, nDOF

      InterpMat(jNode,jNode) = One

      iOS = nDOF

      IF( nNodesE > 1 )THEN

        InterpMat(iOS+1:iOS+nDOF_E,jNode) = L_E_Dn(1:nDOF_E,jNode)

        iOS = iOS + nDOF_E

        InterpMat(iOS+1:iOS+nDOF_E,jNode) = L_E_Up(1:nDOF_E,jNode)

        iOS = iOS + nDOF_E

      END IF

      IF( nNodesX(1) > 1 )THEN

         InterpMat(iOS+1:iOS+nDOF_X1,jNode) = L_X1_Dn(1:nDOF_X1,jNode)

         iOS = iOS + nDOF_X1

         InterpMat(iOS+1:iOS+nDOF_X1,jNode) = L_X1_Up(1:nDOF_X1,jNode)

         iOS = iOS + nDOF_X1

      END IF

      IF( nNodesX(2) > 1 )THEN

         InterpMat(iOS+1:iOS+nDOF_X2,jNode) = L_X2_Dn(1:nDOF_X2,jNode)

         iOS = iOS + nDOF_X2

         InterpMat(iOS+1:iOS+nDOF_X2,jNode) = L_X2_Up(1:nDOF_X2,jNode)

         iOS = iOS + nDOF_X2

      END IF

      IF( nNodesX(3) > 1 )THEN

         InterpMat(iOS+1:iOS+nDOF_X3,jNode) = L_X3_Dn(1:nDOF_X3,jNode)

         iOS = iOS + nDOF_X3

         InterpMat(iOS+1:iOS+nDOF_X3,jNode) = L_X3_Up(1:nDOF_X3,jNode)

      END IF

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: L_X1_Dn, L_X2_Dn, L_X3_Dn, L_E_Dn, &
    !$OMP          L_X1_Up, L_X2_Up, L_X3_Up, L_E_Up, &
    !$OMP          dLdX1_q, dLdX2_q, dLdX3_q, dLdE_q, &
    !$OMP          InterpMat )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( L_X1_Dn, L_X2_Dn, L_X3_Dn, L_E_Dn, &
    !$ACC         L_X1_Up, L_X2_Up, L_X3_Up, L_E_Up, &
    !$ACC         dLdX1_q, dLdX2_q, dLdX3_q, dLdE_q, &
    !$ACC         InterpMat )
#endif

  END SUBROUTINE InitializeReferenceElement_Lagrange


  SUBROUTINE FinalizeReferenceElement_Lagrange

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: L_X1_Dn, L_X2_Dn, L_X3_Dn, L_E_Dn, &
    !$OMP               L_X1_Up, L_X2_Up, L_X3_Up, L_E_Up, &
    !$OMP               dLdX1_q, dLdX2_q, dLdX3_q, dLdE_q, &
    !$OMP               InterpMat )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( L_X1_Dn, L_X2_Dn, L_X3_Dn, L_E_Dn, &
    !$ACC         L_X1_Up, L_X2_Up, L_X3_Up, L_E_Up, &
    !$ACC         dLdX1_q, dLdX2_q, dLdX3_q, dLdE_q, &
    !$ACC         InterpMat )
#endif

    DEALLOCATE( L_E_Dn )
    DEALLOCATE( L_E_Up )
    DEALLOCATE( L_X1_Dn )
    DEALLOCATE( L_X1_Up )
    DEALLOCATE( L_X2_Dn )
    DEALLOCATE( L_X2_Up )
    DEALLOCATE( L_X3_Dn )
    DEALLOCATE( L_X3_Up )

    DEALLOCATE( dLdE_q )
    DEALLOCATE( dLdX1_q )
    DEALLOCATE( dLdX2_q )
    DEALLOCATE( dLdX3_q )

    DEALLOCATE( InterpMat )

  END SUBROUTINE FinalizeReferenceElement_Lagrange


END MODULE ReferenceElementModule_Lagrange
