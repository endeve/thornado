MODULE EulerEquationsLimiterUtilitiesModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_S1, iCF_S2, iCF_S3, nCF, &
    Shock

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeDiscontinuityDetector
  PUBLIC :: DetectDiscontinuities
  PUBLIC :: SolveTheta_Bisection

  LOGICAL, PARAMETER :: UseDiscontinuityDetector = .TRUE.
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Lagrange_X1_P
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Lagrange_X1_M

CONTAINS


  SUBROUTINE InitializeDiscontinuityDetector

    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    INTEGER  :: jNodeX1, jNodeX2, jNodeX3, jNodeX
    REAL(DP) :: eta_X1, eta_X2, eta_X3

    ALLOCATE( Lagrange_X1_P(nDOFX,nDOFX) )
    ALLOCATE( Lagrange_X1_M(nDOFX,nDOFX) )

    DO jNodeX3 = 1, nNodesX(3)
      DO jNodeX2 = 1, nNodesX(2)
        DO jNodeX1 = 1, nNodesX(1)

          jNodeX = NodeNumberX( jNodeX1, jNodeX2, jNodeX3 )

          eta_X1 = MeshX(1) % Nodes(jNodeX1)
          eta_X2 = MeshX(2) % Nodes(jNodeX2)
          eta_X3 = MeshX(3) % Nodes(jNodeX3)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                Lagrange_X1_P(iNodeX,jNodeX) &
                  = L_X1(iNodeX1) % P( eta_X1 + 1.0_DP ) &
                    * L_X2(iNodeX2) % P( eta_X2 ) &
                      * L_X3(iNodeX3) % P( eta_X3 )

                Lagrange_X1_M(iNodeX,jNodeX) &
                  = L_X1(iNodeX1) % P( eta_X1 - 1.0_DP ) &
                    * L_X2(iNodeX2) % P( eta_X2 ) &
                      * L_X3(iNodeX3) % P( eta_X3 )

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeDiscontinuityDetector


  SUBROUTINE DetectDiscontinuities

    INTEGER  :: iX1, iX2, iX3, iCF, k
    REAL(DP) :: F_M, F_P
    REAL(DP), PARAMETER :: alpha = 1.5_DP
    REAL(DP), DIMENSION(nCF)       :: uCF_A
    REAL(DP), DIMENSION(nDOFX,nCF) :: uCF_M
    REAL(DP), DIMENSION(nDOFX,nCF) :: uCF_M_P_X1
    REAL(DP), DIMENSION(nDOFX,nCF) :: uCF_M_N_X1

    Shock(:,:,:) = 0.0_DP

    IF( .NOT. UseDiscontinuityDetector )THEN

      Shock(:,:,:) = 1.0_DP
      RETURN

    END IF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        !$OMP PARALLEL DO PRIVATE &
        !$OMP&              ( iX1, iCF, k, F_M, F_P, uCF_A, &
        !$OMP&                uCF_M_P_X1, uCF_M, uCF_M_N_X1 )
        DO iX1 = 1, nX(1)

          ! --- Map To Modal Representation ---

          DO iCF = 1, nCF

            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1-1,iX2,iX3,iCF), uCF_M_P_X1(:,iCF) )
            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1,  iX2,iX3,iCF), uCF_M     (:,iCF) )
            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1+1,iX2,iX3,iCF), uCF_M_N_X1(:,iCF) )

          END DO

          ! --- Cell-Averaged Quantities ---

          uCF_A(1:nCF) = uCF_M(1,1:nCF)

          ! --- Detect Discontinuity with Harten's Sub-Cell Method ---

          LOOP: DO iCF = 1, nCF

            IF( ANY( iCF == [ iCF_S1, iCF_S2, iCF_S3 ] ) ) CYCLE LOOP

            F_M = 0.0_DP
            F_P = 0.0_DP
            DO k = 1, nDOFX

              F_M &
                = F_M &
                  + WeightsX_q(k) &
                    * SUM( Lagrange_X1_M(:,k) * uCF(:,iX1+1,iX2,iX3,iCF) )

              F_P &
                = F_P &
                  + WeightsX_q(k) &
                    * SUM( Lagrange_X1_P(:,k) * uCF(:,iX1-1,iX2,iX3,iCF) )

            END DO
            F_M = F_M - uCF_A(iCF)
            F_P = F_P - uCF_A(iCF)

            IF( F_M * F_P <= 0.0_DP )THEN
              DO k = 1, MIN( 2, nDOFX )
                IF( ABS( uCF_M(k,iCF) ) &
                    > alpha * ABS( uCF_M_P_X1(k,iCF) ) .AND. &
                    ABS( uCF_M(k,iCF) ) &
                    > alpha * ABS( uCF_M_N_X1(k,iCF) ) ) &
                THEN
                  Shock(iX1,iX2,iX3) = 1.0_DP
                  EXIT LOOP
                END IF
              END DO
            END IF

          END DO LOOP

        END DO
        !$OMP END PARALLEL DO
      END DO
    END DO

  END SUBROUTINE DetectDiscontinuities


  SUBROUTINE SolveTheta_Bisection &
    ( D_K, S1_K, S2_K, S3_K, E_K, dD, dS1, dS2, dS3, dE, Tol_E, Theta )

    REAL(DP), INTENT(in)  :: D_K, S1_K, S2_K, S3_K, E_K
    REAL(DP), INTENT(in)  :: dD,  dS1,  dS2,  dS3,  dE, Tol_E
    REAL(DP), INTENT(out) :: Theta

    LOGICAL  :: Converged
    INTEGER  :: Iter
    REAL(DP) :: a, b, c, ab, Phi_a, Phi_b, Phi_c, Phi_0
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol_ab  = 1.0d-8
    REAL(DP), PARAMETER :: Tol_Phi = 1.0d-8

    a = 0.0_DP
    Phi_a &
      = ThetaFun &
          ( D_K, S1_K, S2_K, S3_K, E_K, dD, dS1, dS2, dS3, dE, Tol_E, a )

    b = 1.0_DP
    Phi_b &
      = ThetaFun &
          ( D_K, S1_K, S2_K, S3_K, E_K, dD, dS1, dS2, dS3, dE, Tol_E, b )

    Phi_0 = Phi_a
    ab    = b - a

    Converged = .FALSE.
    Iter      = 0
    DO WHILE ( .NOT. Converged )

      Iter = Iter + 1

      ab = 0.5_DP * ab
      c  = a + ab

      Phi_c  &
        = ThetaFun &
            ( D_K, S1_K, S2_K, S3_K, E_K, dD, dS1, dS2, dS3, dE, Tol_E, c )

      IF( Phi_a * Phi_c < 0.0_DP )THEN

        b     = c
        Phi_b = Phi_c

      ELSE

        a     = c
        Phi_a = Phi_c

      END IF

      IF( ab < Tol_ab .AND. ABS( Phi_a ) / Phi_0 < Tol_Phi ) &
        Converged = .TRUE.

      IF( Iter > MaxIter .AND. .NOT. Converged )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'SolveTheta_Bisection (Euler Equations)'
        WRITE(*,'(A6,A21,I4.4,A11)') &
          '', 'No Convergence After ', Iter, ' Iterations'
        WRITE(*,*)
      END IF

    END DO

    Theta = a

  END SUBROUTINE SolveTheta_Bisection


  PURE REAL(DP) FUNCTION ThetaFun &
    ( D_K, S1_K, S2_K, S3_K, E_K, dD, dS1, dS2, dS3, dE, Tol_E, Theta )

    REAL(DP), INTENT(in) :: D_K, S1_K, S2_K, S3_K, E_K
    REAL(DP), INTENT(in) :: dD, dS1, dS2, dS3, dE, Tol_E, Theta

    ThetaFun &
      = E_K+Theta*dE &
        - 0.5_DP*DOT_PRODUCT &
                   ( [S1_K+Theta*dS1,S2_K+Theta*dS2,S3_K+Theta*dS3], &
                     [S1_K+Theta*dS1,S2_K+Theta*dS2,S3_K+Theta*dS3] ) &
          / ( D_K+Theta*dD ) - Tol_E

    RETURN
  END FUNCTION ThetaFun


END MODULE EulerEquationsLimiterUtilitiesModule_DG
