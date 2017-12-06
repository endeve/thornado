PROGRAM dgMatrixTest

  USE KindModule, ONLY: &
    DP, Pi
  USE QuadratureModule, ONLY: &
    InitializeQuadratures, &
    GetQuadrature
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE PolynomialBasisModule_Lagrange, ONLY: &
    LagrangeP, &
    dLagrangeP

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: i1, i2, i3, iNode
  INTEGER  :: j1, j2, j3, jNode
  INTEGER  :: k1, k2, k3, kNode
  INTEGER  :: q1, q2, q3
  INTEGER  :: mpierr
  INTEGER  :: NodeNumberTable2D_G4(2,16)
  INTEGER  :: NodeNumberTable2D_G5(2,25)
  INTEGER  :: NodeNumberTable3D(3,64)
  REAL(DP) :: wTime
  REAL(DP) :: xG4(4), wG4(4), xG5(5), wG5(5)
  REAL(DP) :: xL4(4), wL4(4)
  REAL(DP) :: xL(3), xR(3), dX(3)
  REAL(DP) :: xN_L(3,64), xN(3,64)
  REAL(DP) :: J_L(64), J  (64)
  REAL(DP) :: g11_L(64), g22_L(64), g33_L(64)
  REAL(DP) :: g11  (64), g22  (64), g33  (64)
  REAL(DP) :: sg11_L(64), sg22_L(64), sg33_L(64)
  REAL(DP) :: sg11  (64), sg22  (64), sg33  (64)
  REAL(DP) :: w(64), w1(16), w2(16), w3(16)
  REAL(DP) :: L_L2G(64,64)
  REAL(DP) :: L_G2L(64,64)
  REAL(DP) :: L1L_4(16,64), L1H_4(16,64)
  REAL(DP) :: L1L_5(25,64), L1H_5(25,64)
  REAL(DP) :: L2L_4(16,64), L2H_4(16,64)
  REAL(DP) :: L2L_5(25,64), L2H_5(25,64)
  REAL(DP) :: L3L_4(16,64), L3H_4(16,64)
  REAL(DP) :: L3L_5(25,64), L3H_5(25,64)
  REAL(DP) :: dLdX1(64,64), dLdX2(64,64), dLdX3(64,64)
  REAL(DP) :: M0ij_4(64,64), M0ij_5(64,64)
  REAL(DP) :: M1ijL_4(64,16), M1ijH_4(64,16)
  REAL(DP) :: M2ijL_4(64,16), M2ijH_4(64,16)
  REAL(DP) :: M3ijL_4(64,16), M3ijH_4(64,16)
  REAL(DP) :: Mij_4 (64,64), Mij_5 (64,64)
  REAL(DP) :: S1ij_4(64,64), S1ij_5(64,64)
  REAL(DP) :: S2ij_4(64,64), S2ij_5(64,64)
  REAL(DP) :: S3ij_4(64,64), S3ij_5(64,64)

  CALL MPI_INIT( mpierr )

  ! --- Create Node Number Tables ---

  iNode = 1
  DO i3 = 1, 4
    DO i2 = 1, 4
      DO i1 = 1, 4

        NodeNumberTable3D(1:3,iNode) = [ i1, i2, i3 ]

        iNode = iNode + 1

      END DO
    END DO
  END DO

  iNode = 1
  DO i3 = 1, 4
    DO i2 = 1, 4

      NodeNumberTable2D_G4(1:2,iNode) = [ i2, i3 ]

      iNode = iNode + 1

    END DO
  END DO

  iNode = 1
  DO i3 = 1, 5
    DO i2 = 1, 5

      NodeNumberTable2D_G5(1:2,iNode) = [ i2, i3 ]

      iNode = iNode + 1

    END DO
  END DO

  CALL InitializeQuadratures

  ! --- Gaussian Quadrature ---

  CALL GetQuadrature( 4, xG4, wG4, 'Gaussian' )
  CALL GetQuadrature( 5, xG5, wG5, 'Gaussian' )

  ! --- Lobatto Quadrature ---

  CALL GetQuadrature( 4, xL4, wL4, 'Lobatto' )

  ! --- Define Element ---

  xL = [ 0.0_DP, 0.0_DP, 1.0_DP * Pi ]
  xR = [ 1.0_DP, Pi/4.0_DP, 2.0_DP * Pi ]
  dX = xR - xL

  CALL WriteVector( 3, xL, 'xL.dat' )
  CALL WriteVector( 3, xR, 'xH.dat' )
  CALL WriteVector( 3, dX, 'dX.dat' )

  ! --- Compute Geometry Fields in Lobatto Points ---

  DO iNode = 1, 64

    i1 = NodeNumberTable3D(1,iNode)
    i2 = NodeNumberTable3D(2,iNode)
    i3 = NodeNumberTable3D(3,iNode)

    xN_L(:,iNode) = xL + dX * ( 0.5_DP + [ xL4(i1), xL4(i2), xL4(i3) ] )
    xN  (:,iNode) = xL + dX * ( 0.5_DP + [ xG4(i1), xG4(i2), xG4(i3) ] )

    g11_L(iNode) = 1.0_DP
    g22_L(iNode) = xN_L(1,iNode)**2
    g33_L(iNode) = xN_L(1,iNode)**2 * SIN( xN_L(2,iNode) )**2

    J_L(iNode) = SQRT( ( xN_L(1,iNode)**2 * SIN( xN_L(2,iNode) ) )**2 )

    sg11_L(iNode) = 1.0_DP
    sg22_L(iNode) = xN_L(1,iNode)
    sg33_L(iNode) = xN_L(1,iNode) * SIN( xN_L(2,iNode) )

  END DO

  ! --- Write Global Node Coordinates ---

  CALL WriteVector( 64, xN(1,:), 'x1.dat' )
  CALL WriteVector( 64, xN(2,:), 'x2.dat' )
  CALL WriteVector( 64, xN(3,:), 'x3.dat' )

  ! --- Transformation Matrix Lobatto -> Gaussian: L_L2G ---
  ! --- Transformation Matrix Gaussian -> Lobatto: L_G2L ---

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      L_L2G(iNode,jNode) &
        = LagrangeP  (xG4(i1),j1,xL4,4) &
          * LagrangeP(xG4(i2),j2,xL4,4) &
          * LagrangeP(xG4(i3),j3,xL4,4)

      L_G2L(iNode,jNode) &
        = LagrangeP  (xL4(i1),j1,xG4,4) &
          * LagrangeP(xL4(i2),j2,xG4,4) &
          * LagrangeP(xL4(i3),j3,xG4,4)

    END DO

  END DO

  CALL WriteMatrix( 64, 64, L_L2G, 'L_L2G.dat' )

  CALL WriteMatrix( 64, 64, L_G2L, 'L_G2L.dat' )

  ! --- Compute Geometry Fields in Gaussian Points ---

  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, g11_L, 1, 0.0_DP, g11, 1 )
  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, g22_L, 1, 0.0_DP, g22, 1 )
  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, g33_L, 1, 0.0_DP, g33, 1 )
  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, J_L,   1, 0.0_DP, J,   1 )

  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, sg11_L, 1, 0.0_DP, sg11, 1 )
  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, sg22_L, 1, 0.0_DP, sg22, 1 )
  CALL DGEMV( 'N', 64, 64, 1.0_DP, L_L2G, 64, sg33_L, 1, 0.0_DP, sg33, 1 )

  CALL WriteVector( 64, g11, 'g11.dat' )
  CALL WriteVector( 64, g22, 'g22.dat' )
  CALL WriteVector( 64, g33, 'g33.dat' )
  CALL WriteVector( 64, J_L, 'J_L.dat' )
  CALL WriteVector( 64, J,   'J.dat'   )
  CALL WriteVector( 64, sg11, 'sg11.dat' )
  CALL WriteVector( 64, sg22, 'sg22.dat' )
  CALL WriteVector( 64, sg33, 'sg33.dat' )

  ! --- Compute Volume Integral Weights ---

  DO iNode = 1, 64

    i1 = NodeNumberTable3D(1,iNode)
    i2 = NodeNumberTable3D(2,iNode)
    i3 = NodeNumberTable3D(3,iNode)

    w(iNode) = 0.0_DP

    DO q3 = 1, 4
      DO q2 = 1, 4
        DO q1 = 1, 4

          w(iNode) &
            = w(iNode) &
                + wG4(q1) * wG4(q2) * wG4(q3) &
                  * LagrangeP(XG4(q1),i1,xG4,4) &
                  * LagrangeP(XG4(q2),i2,xG4,4) &
                  * LagrangeP(XG4(q3),i3,xG4,4)

        END DO
      END DO
    END DO

  END DO

  CALL WriteVector( 64, w, 'w.dat' )

  ! --- Compute Surface Integral Weights ---

  ! --- Dimension 1 ---

  DO iNode = 1, 16

    i2 = NodeNumberTable2D_G4(1,iNode)
    i3 = NodeNumberTable2D_G4(2,iNode)

    w1(iNode) = 0.0_DP

    DO q3 = 1, 4
      DO q2 = 1, 4

        w1(iNode) &
          = w1(iNode) &
              + wG4(q2) * wG4(q3) &
                * LagrangeP(xG4(q2),i2,XG4,4) &
                * LagrangeP(xG4(q3),i3,XG4,4)

      END DO
    END DO

  END DO

  ! --- Dimension 2 ---

  DO iNode = 1, 16

    i1 = NodeNumberTable2D_G4(1,iNode)
    i3 = NodeNumberTable2D_G4(2,iNode)

    w2(iNode) = 0.0_DP

    DO q3 = 1, 4
      DO q1 = 1, 4

        w2(iNode) &
          = w2(iNode) &
              + wG4(q1) * wG4(q3) &
                * LagrangeP(xG4(q1),i1,XG4,4) &
                * LagrangeP(xG4(q3),i3,XG4,4)

      END DO
    END DO

  END DO

  ! --- Dimension 3 ---

  DO iNode = 1, 16

    i1 = NodeNumberTable2D_G4(1,iNode)
    i2 = NodeNumberTable2D_G4(2,iNode)

    w3(iNode) = 0.0_DP

    DO q2 = 1, 4
      DO q1 = 1, 4

        w3(iNode) &
          = w3(iNode) &
              + wG4(q1) * wG4(q2) &
                * LagrangeP(xG4(q1),i1,XG4,4) &
                * LagrangeP(xG4(q2),i2,XG4,4)

      END DO
    END DO

  END DO

  CALL WriteVector( 16, w1, 'w1.dat' )
  CALL WriteVector( 16, w2, 'w2.dat' )
  CALL WriteVector( 16, w3, 'w3.dat' )

  ! --- Interpolation Matrices ---

  ! --- Dimension 1 ---

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 16

      i2 = NodeNumberTable2D_G4(1,iNode)
      i3 = NodeNumberTable2D_G4(2,iNode)

      L1L_4(iNode,jNode) &
        = LagrangeP( - 0.5_DP,j1,xG4,4) &
          * LagrangeP(xG4(i2),j2,xG4,4) &
          * LagrangeP(xG4(i3),j3,xG4,4)

      L1H_4(iNode,jNode) &
        = LagrangeP( + 0.5_DP,j1,xG4,4) &
          * LagrangeP(xG4(i2),j2,xG4,4) &
          * LagrangeP(xG4(i3),j3,xG4,4)

    END DO

    DO iNode = 1, 25

      i2 = NodeNumberTable2D_G5(1,iNode)
      i3 = NodeNumberTable2D_G5(2,iNode)

      L1L_5(iNode,jNode) &
        = LagrangeP( - 0.5_DP,j1,xG4,4) &
          * LagrangeP(xG5(i2),j2,xG4,4) &
          * LagrangeP(xG5(i3),j3,xG4,4)

      L1H_5(iNode,jNode) &
        = LagrangeP( + 0.5_DP,j1,xG4,4) &
          * LagrangeP(xG5(i2),j2,xG4,4) &
          * LagrangeP(xG5(i3),j3,xG4,4)

    END DO

  END DO

  ! --- Dimension 2 ---

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 16

      i1 = NodeNumberTable2D_G4(1,iNode)
      i3 = NodeNumberTable2D_G4(2,iNode)

      L2L_4(iNode,jNode) &
        = LagrangeP  (xG4(i1),j1,xG4,4) &
          * LagrangeP(-0.5_DP,j2,xG4,4) &
          * LagrangeP(xG4(i3),j3,xG4,4)

      L2H_4(iNode,jNode) &
        = LagrangeP  (xG4(i1),j1,xG4,4) &
          * LagrangeP(+0.5_DP,j2,xG4,4) &
          * LagrangeP(xG4(i3),j3,xG4,4)

    END DO

    DO iNode = 1, 25

      i1 = NodeNumberTable2D_G5(1,iNode)
      i3 = NodeNumberTable2D_G5(2,iNode)

      L2L_5(iNode,jNode) &
        = LagrangeP  (xG5(i1),j1,xG4,4) &
          * LagrangeP(-0.5_DP,j2,xG4,4) &
          * LagrangeP(xG5(i3),j3,xG4,4)

      L2H_5(iNode,jNode) &
        = LagrangeP  (xG5(i1),j1,xG4,4) &
          * LagrangeP(+0.5_DP,j2,xG4,4) &
          * LagrangeP(xG5(i3),j3,xG4,4)

    END DO

  END DO

  ! --- Dimension 3 ---

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 16

      i1 = NodeNumberTable2D_G4(1,iNode)
      i2 = NodeNumberTable2D_G4(2,iNode)

      L3L_4(iNode,jNode) &
        = LagrangeP  (xG4(i1),j1,xG4,4) &
          * LagrangeP(xG4(i2),j2,xG4,4) &
          * LagrangeP(-0.5_DP,j3,xG4,4)

      L3H_4(iNode,jNode) &
        = LagrangeP  (xG4(i1),j1,xG4,4) &
          * LagrangeP(xG4(i2),j2,xG4,4) &
          * LagrangeP(+0.5_DP,j3,xG4,4)

    END DO

    DO iNode = 1, 25

      i1 = NodeNumberTable2D_G5(1,iNode)
      i2 = NodeNumberTable2D_G5(2,iNode)

      L3L_5(iNode,jNode) &
        = LagrangeP  (xG5(i1),j1,xG4,4) &
          * LagrangeP(xG5(i2),j2,xG4,4) &
          * LagrangeP(-0.5_DP,j3,xG4,4)

      L3H_5(iNode,jNode) &
        = LagrangeP  (xG5(i1),j1,xG4,4) &
          * LagrangeP(xG5(i2),j2,xG4,4) &
          * LagrangeP(+0.5_DP,j3,xG4,4)

    END DO

  END DO

  CALL WriteMatrix( 16, 64, L1L_4, 'L1L_4.dat' )
  CALL WriteMatrix( 16, 64, L1H_4, 'L1H_4.dat' )

  CALL WriteMatrix( 25, 64, L1L_5, 'L1L_5.dat' )
  CALL WriteMatrix( 25, 64, L1H_5, 'L1H_5.dat' )

  CALL WriteMatrix( 16, 64, L2L_4, 'L2L_4.dat' )
  CALL WriteMatrix( 16, 64, L2H_4, 'L2H_4.dat' )

  CALL WriteMatrix( 25, 64, L2L_5, 'L2L_5.dat' )
  CALL WriteMatrix( 25, 64, L2H_5, 'L2H_5.dat' )

  CALL WriteMatrix( 16, 64, L3L_4, 'L3L_4.dat' )
  CALL WriteMatrix( 16, 64, L3H_4, 'L3H_4.dat' )

  CALL WriteMatrix( 25, 64, L3L_5, 'L3L_5.dat' )
  CALL WriteMatrix( 25, 64, L3H_5, 'L3H_5.dat' )

  ! --- Derivative Matrices ---

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      dLdX1(iNode,jNode) &
        = dLagrangeP (xG4(i1),j1,xG4,4) &
          * LagrangeP(xG4(i2),j2,xG4,4) &
          * LagrangeP(xG4(i3),j3,xG4,4)

      dLdX2(iNode,jNode) &
        = LagrangeP   (xG4(i1),j1,xG4,4) &
          * dLagrangeP(xG4(i2),j2,xG4,4) &
          * LagrangeP (xG4(i3),j3,xG4,4)

      dLdX3(iNode,jNode) &
        = LagrangeP   (xG4(i1),j1,xG4,4) &
          * LagrangeP (xG4(i2),j2,xG4,4) &
          * dLagrangeP(xG4(i3),j3,xG4,4)

    END DO

  END DO

  CALL WriteMatrix( 64, 64, dLdX1, 'dLdX1.dat' )
  CALL WriteMatrix( 64, 64, dLdX2, 'dLdX2.dat' )
  CALL WriteMatrix( 64, 64, dLdX3, 'dLdX3.dat' )

  ! --- Flux Integration Matrices ---

  ! --- Dimension 1 ---

  DO jNode = 1, 16

    j2 = NodeNumberTable2D_G4(1,jNode)
    j3 = NodeNumberTable2D_G4(2,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      M1ijL_4(iNode,jNode) = 0.0_DP
      M1ijH_4(iNode,jNode) = 0.0_DP

      DO q3 = 1, 4
        DO q2 = 1, 4

          M1ijL_4(iNode,jNode) &
            = M1ijL_4(iNode,jNode) &
                + wG4(q2) * wG4(q3) &
                  ! ---
                  * LagrangeP(-0.5_DP,i1,xG4,4) &
                  * LagrangeP(xG4(q2),i2,xG4,4) &
                  * LagrangeP(xG4(q3),i3,xG4,4) &
                  ! ---
                  * LagrangeP(xG4(q2),j2,xG4,4) &
                  * LagrangeP(xG4(q3),j3,xG4,4)

          M1ijH_4(iNode,jNode) &
            = M1ijH_4(iNode,jNode) &
                + wG4(q2) * wG4(q3) &
                  ! ---
                  * LagrangeP(+0.5_DP,i1,xG4,4) &
                  * LagrangeP(xG4(q2),i2,xG4,4) &
                  * LagrangeP(xG4(q3),i3,xG4,4) &
                  ! ---
                  * LagrangeP(xG4(q2),j2,xG4,4) &
                  * LagrangeP(xG4(q3),j3,xG4,4)

        END DO
      END DO

    END DO

  END DO

  ! --- Dimension 2 ---

  DO jNode = 1, 16

    j1 = NodeNumberTable2D_G4(1,jNode)
    j3 = NodeNumberTable2D_G4(2,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      M2ijL_4(iNode,jNode) = 0.0_DP
      M2ijH_4(iNode,jNode) = 0.0_DP

      DO q3 = 1, 4
        DO q1 = 1, 4

          M2ijL_4(iNode,jNode) &
            = M2ijL_4(iNode,jNode) &
                + wG4(q1) * wG4(q3) &
                  ! ---
                  * LagrangeP(xG4(q1),i1,xG4,4) &
                  * LagrangeP(-0.5_DP,i2,xG4,4) &
                  * LagrangeP(xG4(q3),i3,xG4,4) &
                  ! ---
                  * LagrangeP(xG4(q1),j1,xG4,4) &
                  * LagrangeP(xG4(q3),j3,xG4,4)

          M2ijH_4(iNode,jNode) &
            = M2ijH_4(iNode,jNode) &
                + wG4(q1) * wG4(q3) &
                  ! ---
                  * LagrangeP(xG4(q1),i1,xG4,4) &
                  * LagrangeP(+0.5_DP,i2,xG4,4) &
                  * LagrangeP(xG4(q3),i3,xG4,4) &
                  ! ---
                  * LagrangeP(xG4(q1),j1,xG4,4) &
                  * LagrangeP(xG4(q3),j3,xG4,4)

        END DO
      END DO

    END DO

  END DO

  ! --- Dimension 3 ---

  DO jNode = 1, 16

    j1 = NodeNumberTable2D_G4(1,jNode)
    j2 = NodeNumberTable2D_G4(2,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      M3ijL_4(iNode,jNode) = 0.0_DP
      M3ijH_4(iNode,jNode) = 0.0_DP

      DO q2 = 1, 4
        DO q1 = 1, 4

          M3ijL_4(iNode,jNode) &
            = M3ijL_4(iNode,jNode) &
                + wG4(q1) * wG4(q2) &
                  ! ---
                  * LagrangeP(xG4(q1),i1,xG4,4) &
                  * LagrangeP(xG4(q2),i2,xG4,4) &
                  * LagrangeP(-0.5_DP,i3,xG4,4) &
                  ! ---
                  * LagrangeP(xG4(q1),j1,xG4,4) &
                  * LagrangeP(xG4(q2),j2,xG4,4)

          M3ijH_4(iNode,jNode) &
            = M3ijH_4(iNode,jNode) &
                + wG4(q1) * wG4(q2) &
                  ! ---
                  * LagrangeP(xG4(q1),i1,xG4,4) &
                  * LagrangeP(xG4(q2),i2,xG4,4) &
                  * LagrangeP(+0.5_DP,i3,xG4,4) &
                  ! ---
                  * LagrangeP(xG4(q1),j1,xG4,4) &
                  * LagrangeP(xG4(q2),j2,xG4,4)

        END DO
      END DO

    END DO

  END DO

  CALL WriteMatrix( 64, 16, M1ijL_4, 'M1ijL_4.dat' )
  CALL WriteMatrix( 64, 16, M1ijH_4, 'M1ijH_4.dat' )

  CALL WriteMatrix( 64, 16, M2ijL_4, 'M2ijL_4.dat' )
  CALL WriteMatrix( 64, 16, M2ijH_4, 'M2ijH_4.dat' )

  CALL WriteMatrix( 64, 16, M3ijL_4, 'M3ijL_4.dat' )
  CALL WriteMatrix( 64, 16, M3ijH_4, 'M3ijH_4.dat' )

  ! --- Compute Triple Index Mass Matrices ---

  wTime = MPI_WTIME( )

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      M0ij_4(iNode,jNode) = 0.0_DP
      M0ij_5(iNode,jNode) = 0.0_DP

      DO kNode = 1, 64

        k1 = NodeNumberTable3D(1,kNode)
        k2 = NodeNumberTable3D(2,kNode)
        k3 = NodeNumberTable3D(3,kNode)

        ! --- 4-Point Gaussian Quadrature

        DO q3 = 1, 4
          DO q2 = 1, 4
            DO q1 = 1, 4

              M0ij_4(iNode,jNode) &
                = M0ij_4(iNode,jNode) &
                  + wG4(q1)*wG4(q2)*wG4(q3) &
                    ! --- 
                    * LagrangeP(xG4(q1),k1,xG4,4) &
                    * LagrangeP(xG4(q2),k2,xG4,4) &
                    * LagrangeP(xG4(q3),k3,xG4,4) &
                    ! ---
                    * LagrangeP(xG4(q1),i1,xG4,4) &
                    * LagrangeP(xG4(q2),i2,xG4,4) &
                    * LagrangeP(xG4(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG4(q1),j1,xG4,4) &
                    * LagrangeP(xG4(q2),j2,xG4,4) &
                    * LagrangeP(xG4(q3),j3,xG4,4) &
                    * J(kNode)

            END DO
          END DO
        END DO

        ! --- 5-Point Gaussian Quadrature

        DO q3 = 1, 5
          DO q2 = 1, 5
            DO q1 = 1, 5

              M0ij_5(iNode,jNode) &
                = M0ij_5(iNode,jNode) &
                  + wG5(q1)*wG5(q2)*wG5(q3) &
                    ! --- 
                    * LagrangeP(xG5(q1),k1,xG4,4) &
                    * LagrangeP(xG5(q2),k2,xG4,4) &
                    * LagrangeP(xG5(q3),k3,xG4,4) &
                    ! ---
                    * LagrangeP(xG5(q1),i1,xG4,4) &
                    * LagrangeP(xG5(q2),i2,xG4,4) &
                    * LagrangeP(xG5(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG5(q1),j1,xG4,4) &
                    * LagrangeP(xG5(q2),j2,xG4,4) &
                    * LagrangeP(xG5(q3),j3,xG4,4) &
                    * J(kNode)

            END DO
          END DO
        END DO

      END DO
    END DO
  END DO

  wTime = MPI_WTIME( ) - wTime

  PRINT*
  PRINT*, "  M(0): ", wTime
  PRINT*

  ! --- Remove Numerical Noise ---

  DO jNode = 1, 64
    DO iNode = 1, 64

      IF( ABS( M0ij_4(iNode,jNode) ) < 1.0d-12 ) &
        M0ij_4(iNode,jNode) = 0.0_DP

      IF( ABS( M0ij_5(iNode,jNode) ) < 1.0d-12 ) &
        M0ij_5(iNode,jNode) = 0.0_DP

    END DO
  END DO

  CALL WriteMatrix( 64, 64, M0ij_4, 'M0ij_4.dat' )

  CALL WriteMatrix( 64, 64, M0ij_5, 'M0ij_5.dat' )

  ! --- Compute Mass and Stiffness Matrices ---

  wTime = MPI_WTIME( )

  wTime = MPI_WTIME( )

  DO jNode = 1, 64

    j1 = NodeNumberTable3D(1,jNode)
    j2 = NodeNumberTable3D(2,jNode)
    j3 = NodeNumberTable3D(3,jNode)

    DO iNode = 1, 64

      i1 = NodeNumberTable3D(1,iNode)
      i2 = NodeNumberTable3D(2,iNode)
      i3 = NodeNumberTable3D(3,iNode)

      ! --- 4-Point Gaussian Quadrature

      Mij_4 (iNode,jNode) = 0.0_DP
      S1ij_4(iNode,jNode) = 0.0_DP
      S2ij_4(iNode,jNode) = 0.0_DP
      S3ij_4(iNode,jNode) = 0.0_DP

      DO q3 = 1, 4
        DO q2 = 1, 4
          DO q1 = 1, 4

            Mij_4(iNode,jNode) &
              = Mij_4(iNode,jNode) &
                  + wG4(q1)*wG4(q2)*wG4(q3) &
                    ! ---
                    * LagrangeP(xG4(q1),i1,xG4,4) &
                    * LagrangeP(xG4(q2),i2,xG4,4) &
                    * LagrangeP(xG4(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG4(q1),j1,xG4,4) &
                    * LagrangeP(xG4(q2),j2,xG4,4) &
                    * LagrangeP(xG4(q3),j3,xG4,4)

            S1ij_4(iNode,jNode) &
              = S1ij_4(iNode,jNode) &
                  + wG4(q1)*wG4(q2)*wG4(q3) &
                    ! ---
                    * dLagrangeP(xG4(q1),i1,xG4,4) &
                    *  LagrangeP(xG4(q2),i2,xG4,4) &
                    *  LagrangeP(xG4(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG4(q1),j1,xG4,4) &
                    * LagrangeP(xG4(q2),j2,xG4,4) &
                    * LagrangeP(xG4(q3),j3,xG4,4)

            S2ij_4(iNode,jNode) &
              = S2ij_4(iNode,jNode) &
                  + wG4(q1)*wG4(q2)*wG4(q3) &
                    ! ---
                    *  LagrangeP(xG4(q1),i1,xG4,4) &
                    * dLagrangeP(xG4(q2),i2,xG4,4) &
                    *  LagrangeP(xG4(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG4(q1),j1,xG4,4) &
                    * LagrangeP(xG4(q2),j2,xG4,4) &
                    * LagrangeP(xG4(q3),j3,xG4,4)

            S3ij_4(iNode,jNode) &
              = S3ij_4(iNode,jNode) &
                  + wG4(q1)*wG4(q2)*wG4(q3) &
                    ! ---
                    *  LagrangeP(xG4(q1),i1,xG4,4) &
                    *  LagrangeP(xG4(q2),i2,xG4,4) &
                    * dLagrangeP(xG4(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG4(q1),j1,xG4,4) &
                    * LagrangeP(xG4(q2),j2,xG4,4) &
                    * LagrangeP(xG4(q3),j3,xG4,4)

          END DO
        END DO
      END DO

      ! --- 5-Point Gaussian Quadrature

      Mij_5 (iNode,jNode) = 0.0_DP
      S1ij_5(iNode,jNode) = 0.0_DP
      S2ij_5(iNode,jNode) = 0.0_DP
      S3ij_5(iNode,jNode) = 0.0_DP

      DO q3 = 1, 5
        DO q2 = 1, 5
          DO q1 = 1, 5

            Mij_5(iNode,jNode) &
              = Mij_5(iNode,jNode) &
                  + wG5(q1)*wG5(q2)*wG5(q3) &
                    ! ---
                    * LagrangeP(xG5(q1),i1,xG4,4) &
                    * LagrangeP(xG5(q2),i2,xG4,4) &
                    * LagrangeP(xG5(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG5(q1),j1,xG4,4) &
                    * LagrangeP(xG5(q2),j2,xG4,4) &
                    * LagrangeP(xG5(q3),j3,xG4,4)

            S1ij_5(iNode,jNode) &
              = S1ij_5(iNode,jNode) &
                  + wG5(q1)*wG5(q2)*wG5(q3) &
                    ! ---
                    * dLagrangeP(xG5(q1),i1,xG4,4) &
                    *  LagrangeP(xG5(q2),i2,xG4,4) &
                    *  LagrangeP(xG5(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG5(q1),j1,xG4,4) &
                    * LagrangeP(xG5(q2),j2,xG4,4) &
                    * LagrangeP(xG5(q3),j3,xG4,4)

            S2ij_5(iNode,jNode) &
              = S2ij_5(iNode,jNode) &
                  + wG5(q1)*wG5(q2)*wG5(q3) &
                    ! ---
                    *  LagrangeP(xG5(q1),i1,xG4,4) &
                    * dLagrangeP(xG5(q2),i2,xG4,4) &
                    *  LagrangeP(xG5(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG5(q1),j1,xG4,4) &
                    * LagrangeP(xG5(q2),j2,xG4,4) &
                    * LagrangeP(xG5(q3),j3,xG4,4)

            S3ij_5(iNode,jNode) &
              = S3ij_5(iNode,jNode) &
                  + wG5(q1)*wG5(q2)*wG5(q3) &
                    ! ---
                    *  LagrangeP(xG5(q1),i1,xG4,4) &
                    *  LagrangeP(xG5(q2),i2,xG4,4) &
                    * dLagrangeP(xG5(q3),i3,xG4,4) &
                    ! ---
                    * LagrangeP(xG5(q1),j1,xG4,4) &
                    * LagrangeP(xG5(q2),j2,xG4,4) &
                    * LagrangeP(xG5(q3),j3,xG4,4)

          END DO
        END DO
      END DO

    END DO

  END DO

  wTime = MPI_WTIME( ) - wTime

  PRINT*
  PRINT*, "  M and Si: ", wTime
  PRINT*

  CALL WriteMatrix( 64, 64, Mij_4,  'Mij_4.dat'  )

  CALL WriteMatrix( 64, 64, Mij_5,  'Mij_5.dat'  )

  CALL WriteMatrix( 64, 64, S1ij_4, 'S1ij_4.dat' )

  CALL WriteMatrix( 64, 64, S1ij_5, 'S1ij_5.dat' )

  CALL WriteMatrix( 64, 64, S1ij_4, 'S1ij_4.dat' )

  CALL WriteMatrix( 64, 64, S1ij_5, 'S1ij_5.dat' )

  CALL WriteMatrix( 64, 64, S2ij_4, 'S2ij_4.dat' )

  CALL WriteMatrix( 64, 64, S2ij_5, 'S2ij_5.dat' )

  CALL WriteMatrix( 64, 64, S3ij_4, 'S3ij_4.dat' )

  CALL WriteMatrix( 64, 64, S3ij_5, 'S3ij_5.dat' )

  CALL MPI_FINALIZE( mpierr )

END PROGRAM dgMatrixTest
