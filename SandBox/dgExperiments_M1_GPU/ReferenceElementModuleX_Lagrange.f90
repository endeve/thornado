MODULE ReferenceElementModuleX_Lagrange

  USE KindModule, ONLY: &
    DP, Half
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE UtilitiesModule, ONLY: &
    WriteMatrix
  USE PolynomialBasisModule_Lagrange, ONLY: &
     L_X1,  L_X2,  L_X3, &
    dL_X1, dL_X2, dL_X3, &
    LagrangeP
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    NodeNumberTableX, &
    NodeNumberTableX_X1, &
    NodeNumberTableX_X2, &
    NodeNumberTableX_X3, &
    NodesX1,  NodesX2,  NodesX3, &
    NodesLX1, NodesLX2, NodesLX3

  IMPLICIT NONE
  PRIVATE

  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_X1_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_X1_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_X2_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_X2_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_X3_Dn
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_X3_Up
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLXdX1_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLXdX2_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dLXdX3_q
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LX_L2G

  PUBLIC :: InitializeReferenceElementX_Lagrange
  PUBLIC :: FinalizeReferenceElementX_Lagrange

CONTAINS


  SUBROUTINE InitializeReferenceElementX_Lagrange

    INTEGER :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    INTEGER :: jNodeX, jNodeX1, jNodeX2, jNodeX3
    INTEGER :: iNodeX_X1, iNodeX2_X1, iNodeX3_X1
    INTEGER :: iNodeX_X2, iNodeX1_X2, iNodeX3_X2
    INTEGER :: iNodeX_X3, iNodeX1_X3, iNodeX2_X3

    ALLOCATE( LX_X1_Dn(nDOFX_X1,nDOFX) )
    ALLOCATE( LX_X1_Up(nDOFX_X1,nDOFX) )

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      DO iNodeX_X1 = 1, nDOFX_X1

        iNodeX2_X1 = NodeNumberTableX_X1(1,iNodeX_X1)
        iNodeX3_X1 = NodeNumberTableX_X1(2,iNodeX_X1)

        LX_X1_Dn(iNodeX_X1,iNodeX) &
          = L_X1(iNodeX1) % P( - Half ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X1) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X1) )

        LX_X1_Up(iNodeX_X1,iNodeX) &
          = L_X1(iNodeX1) % P( + Half ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X1) ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X1) )

      END DO

    END DO

    ALLOCATE( LX_X2_Dn(nDOFX_X2,nDOFX) )
    ALLOCATE( LX_X2_Up(nDOFX_X2,nDOFX) )

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      DO iNodeX_X2 = 1, nDOFX_X2

        iNodeX1_X2 = NodeNumberTableX_X2(1,iNodeX_X2)
        iNodeX3_X2 = NodeNumberTableX_X2(2,iNodeX_X2)

        LX_X2_Dn(iNodeX_X2,iNodeX) &
          = L_X1(iNodeX1) % P( NodesX1(iNodeX1_X2) ) &
            * L_X2(iNodeX2) % P( - Half ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X2) )

        LX_X2_Up(iNodeX_X2,iNodeX) &
          = L_X1(iNodeX1) % P( NodesX1(iNodeX1_X2) ) &
            * L_X2(iNodeX2) % P( + Half ) &
            * L_X3(iNodeX3) % P( NodesX3(iNodeX3_X2) )

      END DO

    END DO

    ALLOCATE( LX_X3_Dn(nDOFX_X3,nDOFX) )
    ALLOCATE( LX_X3_Up(nDOFX_X3,nDOFX) )

    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      iNodeX2 = NodeNumberTableX(2,iNodeX)
      iNodeX3 = NodeNumberTableX(3,iNodeX)

      DO iNodeX_X3 = 1, nDOFX_X3

        iNodeX1_X3 = NodeNumberTableX_X3(1,iNodeX_X3)
        iNodeX2_X3 = NodeNumberTableX_X3(2,iNodeX_X3)

        LX_X3_Dn(iNodeX_X3,iNodeX) &
          = L_X1(iNodeX1) % P( NodesX1(iNodeX1_X3) ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X3) ) &
            * L_X3(iNodeX3) % P( - Half )

        LX_X3_Up(iNodeX_X3,iNodeX) &
          = L_X1(iNodeX1) % P( NodesX1(iNodeX1_X3) ) &
            * L_X2(iNodeX2) % P( NodesX2(iNodeX2_X3) ) &
            * L_X3(iNodeX3) % P( + Half )

      END DO

    END DO

    ALLOCATE( dLXdX1_q(nDOFX,nDOFX) )
    ALLOCATE( dLXdX2_q(nDOFX,nDOFX) )
    ALLOCATE( dLXdX3_q(nDOFX,nDOFX) )

    DO jNodeX = 1, nDOFX

      jNodeX1 = NodeNumberTableX(1,jNodeX)
      jNodeX2 = NodeNumberTableX(2,jNodeX)
      jNodeX3 = NodeNumberTableX(3,jNodeX)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        dLXdX1_q(iNodeX,jNodeX) &
          =  dL_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            * L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            * L_X3(jNodeX3) % P( NodesX3(iNodeX3) )

        dLXdX2_q(iNodeX,jNodeX) &
          =    L_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            * dL_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            *  L_X3(jNodeX3) % P( NodesX3(iNodeX3) )

        dLXdX3_q(iNodeX,jNodeX) &
          =    L_X1(jNodeX1) % P( NodesX1(iNodeX1) ) &
            *  L_X2(jNodeX2) % P( NodesX2(iNodeX2) ) &
            * dL_X3(jNodeX3) % P( NodesX3(iNodeX3) )

      END DO

    END DO

    ALLOCATE( LX_L2G(nDOFX,nDOFX) )

    DO jNodeX = 1, nDOFX

      jNodeX1 = NodeNumberTableX(1,jNodeX)
      jNodeX2 = NodeNumberTableX(2,jNodeX)
      jNodeX3 = NodeNumberTableX(3,jNodeX)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        LX_L2G(iNodeX,jNodeX) &
          = LagrangeP  ( NodesX1(iNodeX1), jNodeX1, NodesLX1, nNodesX(1) ) &
            * LagrangeP( NodesX2(iNodeX2), jNodeX2, NodesLX2, nNodesX(2) ) &
            * LagrangeP( NodesX3(iNodeX3), jNodeX3, NodesLX3, nNodesX(3) )

      END DO

    END DO

    CALL WriteMatrix( nDOFX_X1, nDOFX, LX_X1_Dn, 'LX_X1_Dn.dat' )
    CALL WriteMatrix( nDOFX_X1, nDOFX, LX_X1_Up, 'LX_X1_Up.dat' )
    CALL WriteMatrix( nDOFX_X2, nDOFX, LX_X2_Dn, 'LX_X2_Dn.dat' )
    CALL WriteMatrix( nDOFX_X2, nDOFX, LX_X2_Up, 'LX_X2_Up.dat' )
    CALL WriteMatrix( nDOFX_X3, nDOFX, LX_X3_Dn, 'LX_X3_Dn.dat' )
    CALL WriteMatrix( nDOFX_X3, nDOFX, LX_X3_Up, 'LX_X3_Up.dat' )

    CALL WriteMatrix( nDOFX, nDOFX, dLXdX1_q, 'dLXdX1_q.dat' )
    CALL WriteMatrix( nDOFX, nDOFX, dLXdX2_q, 'dLXdX2_q.dat' )
    CALL WriteMatrix( nDOFX, nDOFX, dLXdX3_q, 'dLXdX3_q.dat' )

    CALL WriteMatrix( nDOFX, nDOFX, LX_L2G, 'LX_L2G.dat' )

  END SUBROUTINE InitializeReferenceElementX_Lagrange


  SUBROUTINE FinalizeReferenceElementX_Lagrange

    DEALLOCATE( LX_X1_Dn )
    DEALLOCATE( LX_X1_Up )
    DEALLOCATE( LX_X2_Dn )
    DEALLOCATE( LX_X2_Up )
    DEALLOCATE( LX_X3_Dn )
    DEALLOCATE( LX_X3_Up )
    DEALLOCATE( dLXdX1_q )
    DEALLOCATE( dLXdX2_q )
    DEALLOCATE( dLXdX3_q )
    DEALLOCATE( LX_L2G )

  END SUBROUTINE FinalizeReferenceElementX_Lagrange


END MODULE ReferenceElementModuleX_Lagrange
