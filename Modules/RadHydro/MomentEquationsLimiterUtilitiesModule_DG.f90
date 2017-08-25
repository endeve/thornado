MODULE MomentEquationsLimiterUtilitiesModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOF, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_E, L_X1, L_X2, L_X3
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Radiation, &
    MapModalToNodal_Radiation
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    WeightsR, &
    uCR, iCR_N, nCR, &
    Discontinuity

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeDiscontinuityDetector
  PUBLIC :: DetectDiscontinuities
  PUBLIC :: PackSpatialNodes
  PUBLIC :: UnpackSpatialNodes
  PUBLIC :: SolveTheta_Bisection

  LOGICAL, PARAMETER :: UseDiscontinuityDetector = .TRUE.
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Lagrange_X1_P
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Lagrange_X1_M

CONTAINS


  SUBROUTINE InitializeDiscontinuityDetector

    INTEGER  :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode
    INTEGER  :: jNodeE, jNodeX1, jNodeX2, jNodeX3, jNode
    REAL(DP) :: eta_E, eta_X1, eta_X2, eta_X3

    ALLOCATE( Lagrange_X1_P(nDOF,nDOF) )
    ALLOCATE( Lagrange_X1_M(nDOF,nDOF) )

    DO jNodeX3 = 1, nNodesX(3)
      DO jNodeX2 = 1, nNodesX(2)
        DO jNodeX1 = 1, nNodesX(1)
          DO jNodeE = 1, nNodesE

            jNode = NodeNumber( jNodeE, jNodeX1, jNodeX2, jNodeX3 )

            eta_E  = MeshE    % Nodes(jNodeE)
            eta_X1 = MeshX(1) % Nodes(jNodeX1)
            eta_X2 = MeshX(2) % Nodes(jNodeX2)
            eta_X3 = MeshX(3) % Nodes(jNodeX3)

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    Lagrange_X1_P(iNode,jNode) &
                      = L_E(iNodeE) % P( eta_E ) &
                        * L_X1(iNodeX1) % P( eta_X1 + 1.0_DP ) &
                          * L_X2(iNodeX2) % P( eta_X2 ) &
                            * L_X3(iNodeX3) % P( eta_X3 )

                    Lagrange_X1_M(iNode,jNode) &
                      = L_E(iNodeE) % P( eta_E ) &
                        * L_X1(iNodeX1) % P( eta_X1 - 1.0_DP ) &
                          * L_X2(iNodeX2) % P( eta_X2 ) &
                            * L_X3(iNodeX3) % P( eta_X3 )

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeDiscontinuityDetector


  SUBROUTINE DetectDiscontinuities

    INTEGER :: iS, iE, iX1, iX2, iX3, iCR, k
    REAL(DP) :: F_M, F_P
    REAL(DP), DIMENSION(nCR) :: uCR_A
    REAL(DP), DIMENSION(nDOF,nCR) :: uCR_M
    REAL(DP), DIMENSION(nDOF,nCR) :: uCR_M_P_X1
    REAL(DP), DIMENSION(nDOF,nCR) :: uCR_M_N_X1
    REAL(DP), PARAMETER :: alpha = 1.5_DP

    Discontinuity(:,:,:,:) = 0.0_DP

    IF( .NOT. UseDiscontinuityDetector )THEN

      Discontinuity(:,:,:,:) = 1.0_DP
      RETURN

    END IF

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          !$OMP PARALLEL DO PRIVATE &
          !$OMP&              ( iX1, iE, iCR, k, F_M, F_P, uCR_A, &
          !$OMP&                uCR_M_P_X1, uCR_M, uCR_M_N_X1 )
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              ! --- Map To Modal Representation ---

              DO iCR = 1, nCR

                CALL MapNodalToModal_Radiation &
                       ( uCR(:,iE,iX1-1,iX2,iX3,iCR,iS), uCR_M_P_X1(:,iCR) )
                CALL MapNodalToModal_Radiation &
                       ( uCR(:,iE,iX1,  iX2,iX3,iCR,iS), uCR_M     (:,iCR) )
                CALL MapNodalToModal_Radiation &
                       ( uCR(:,iE,iX1+1,iX2,iX3,iCR,iS), uCR_M_N_X1(:,iCR) )

              END DO

              ! --- Cell-Averaged Quantities ---

              uCR_A(1:nCR) = uCR_M(1,1:nCR)

              ! --- Detect Discontinuity with Harten's Sub-Cell Method ---

              F_M = 0.0_DP
              F_P = 0.0_DP
              DO k = 1, nDOF

                F_M &
                  = F_M &
                    + WeightsR(k) &
                      * SUM( Lagrange_X1_M(:,k) &
                             * uCR(:,iE,iX1+1,iX2,iX3,iCR_N,iS) )

                F_P &
                  = F_P &
                    + WeightsR(k) &
                      * SUM( Lagrange_X1_P(:,k) &
                             * uCR(:,iE,iX1-1,iX2,iX3,iCR_N,iS) )

              END DO
              F_M = F_M - uCR_A(iCR_N)
              F_P = F_P - uCR_A(iCR_N)

              IF( F_M * F_P <= 0.0_DP )THEN
                LOOP: DO k = 1, nDOF

                  IF( ALL( [ ABS( uCR_M     (k,iCR_N) ), &
                             ABS( uCR_M_P_X1(k,iCR_N) ), &
                             ABS( uCR_M_N_X1(k,iCR_N) ) ] * 1.0d12 &
                           < 1.0_DP ) ) CYCLE LOOP

                  IF( ABS( uCR_M(k,iCR_N) ) &
                      > alpha * ABS( uCR_M_P_X1(k,iCR_N) ) .AND. &
                      ABS( uCR_M(k,iCR_N) ) &
                      > alpha * ABS( uCR_M_N_X1(k,iCR_N) ) ) &
                  THEN
                    Discontinuity(iE,iX1,iX2,iX3) = 1.0_DP
                    EXIT LOOP
                  END IF

                END DO LOOP
              END IF

            END DO ! iE
          END DO ! iX1
          !$OMP END PARALLEL DO
        END DO ! iX2
      END DO ! iX3

    END DO ! iS

  END SUBROUTINE DetectDiscontinuities


  SUBROUTINE PackSpatialNodes( iNodeE, u, uX )

    INTEGER,                    INTENT(in)  :: iNodeE
    REAL(DP), DIMENSION(nDOF),  INTENT(in)  :: u
    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uX

    INTEGER :: iNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNodeX, iNode

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iNodeX &
            = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
          iNode &
            = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

          uX(iNodeX) = u(iNode)

        END DO
      END DO
    END DO

  END SUBROUTINE PackSpatialNodes


  SUBROUTINE UnpackSpatialNodes( iNodeE, u, uX )

    INTEGER,                    INTENT(in)  :: iNodeE
    REAL(DP), DIMENSION(nDOF),  INTENT(out) :: u
    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uX

    INTEGER :: iNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNodeX, iNode

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iNodeX &
            = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
          iNode &
            = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

          u(iNode) = uX(iNodeX)

        END DO
      END DO
    END DO

  END SUBROUTINE UnpackSpatialNodes


  SUBROUTINE SolveTheta_Bisection &
    ( N_K, G1_K, G2_K, G3_K, dN, dG1, dG2, dG3, Tol_G, Theta )

    REAL(DP), INTENT(in)  :: N_K, G1_K, G2_K, G3_K, dN, dG1, dG2, dG3, Tol_G
    REAL(DP), INTENT(out) :: Theta

    LOGICAL  :: Converged
    INTEGER  :: Iter
    REAL(DP) :: a, b, c, ab, Phi_a, Phi_b, Phi_c, Phi_0
    REAL(DP) :: Phi_a0, Phi_b0
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol_ab  = 1.0d-8
    REAL(DP), PARAMETER :: Tol_Phi = 1.0d-8

    a = 0.0_DP
    Phi_a &
      = ThetaFun &
          ( N_K, G1_K, G2_K, G3_K, dN, dG1, dG2, dG3, Tol_G, a )

    Phi_a0 &
      = Phi_a

    b = 1.0_DP
    Phi_b &
      = ThetaFun &
          ( N_K, G1_K, G2_K, G3_K, dN, dG1, dG2, dG3, Tol_G, b )

    Phi_b0 &
      = Phi_b

    Phi_0 = Phi_a
    ab    = b - a

    Converged = .FALSE.
    Iter      = 0
    DO WHILE ( .NOT. Converged )

      Iter = Iter + 1

      ab = 0.5_DP * ab
      c  = a + ab

      Phi_c &
      = ThetaFun &
          ( N_K, G1_K, G2_K, G3_K, dN, dG1, dG2, dG3, Tol_G, c )

      IF( Phi_a * Phi_c < 0.0_DP )THEN

        b     = c
        Phi_b = Phi_c

      ELSE

        a     = c
        Phi_a = Phi_c

      END IF

      IF( ab < Tol_ab .AND. Phi_a / Phi_0 < Tol_Phi ) &
        Converged = .TRUE.

      IF( Iter > MaxIter .AND. .NOT. Converged )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'SolveTheta_Bisection (Moment Equations)'
        WRITE(*,'(A6,A21,I4.4,A11)') &
          '', 'No Convergence After ', Iter, ' Iterations'
        WRITE(*,*)
        WRITE(*,*) '  a, b, ab = ', a, b, ab
        WRITE(*,*) '  Phi_a0, Phi_b0 = ', Phi_a0, Phi_b0
        WRITE(*,*) '  Phi_ai, Phi_bi = ', Phi_a,  Phi_b
        WRITE(*,*) '  N_K, G_K = ', N_K, G1_K, G2_K, G3_K
        WRITE(*,'(A8,ES20.6e3)') 'Tol_G = ', Tol_G
        STOP
      END IF

    END DO

    Theta = a

  END SUBROUTINE SolveTheta_Bisection


  PURE REAL(DP) FUNCTION ThetaFun &
    ( N_K, G1_K, G2_K, G3_K, dN, dG1, dG2, dG3, Tol_G, Theta )

    REAL(DP), INTENT(in) :: N_K, G1_K, G2_K, G3_K
    REAL(DP), INTENT(in) :: dN, dG1, dG2, dG3, Tol_G, Theta

    ThetaFun &
      = N_K+Theta*dN &
        - SQRT( DOT_PRODUCT &
                  ( [G1_K+Theta*dG1,G2_K+Theta*dG2,G3_K+Theta*dG3], &
                    [G1_K+Theta*dG1,G2_K+Theta*dG2,G3_K+Theta*dG3] ) ) &
        - Tol_G

    RETURN
  END FUNCTION ThetaFun


END MODULE MomentEquationsLimiterUtilitiesModule_DG
