MODULE MomentEquationsSlopeLimiterModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOF, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    MinModB
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Radiation_X, &
    MapModalToNodal_Radiation_X
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    WeightsGX, &
    VolJacX, VolX
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR, &
    Discontinuity
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeEigenvectors_L, &
    ComputeEigenvectors_R
  USE MomentEquationsLimiterUtilitiesModule_DG, ONLY: &
    DetectDiscontinuities, &
    PackSpatialNodes, &
    UnpackSpatialNodes
  USE MomentEquationsLimiterModule_DG, ONLY: &
    ApplySlopeLimiter, &
    BetaTVD, &
    BetaTVB, &
    LegendreX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplySlopeLimiter_M1_DG

  LOGICAL,  PARAMETER :: Debug         = .FALSE.
  LOGICAL,  PARAMETER :: Componentwise = .FALSE.
  REAL(DP), PARAMETER :: Tol_TVD       = 1.0d-2

CONTAINS


  SUBROUTINE ApplySlopeLimiter_M1_DG

    LOGICAL :: LimitPolynomial
    INTEGER :: iS, iE, iX1, iX2, iX3, iCR
    INTEGER :: iNodeE
    INTEGER :: LWORK, INFO
    REAL(DP), DIMENSION(1)   :: d
    REAL(DP), DIMENSION(2)   :: c
    REAL(DP), DIMENSION(5)   :: WORK
    REAL(DP), DIMENSION(2,2) :: A0, A
    REAL(DP), DIMENSION(1,2) :: B0, B
    REAL(DP), DIMENSION(nCR) :: uCR_A, uCR_A_P_X1, uCR_A_N_X1
    REAL(DP), DIMENSION(nCR) :: SlopeDifference
    REAL(DP), DIMENSION(nCR) :: uCR_K_0, uCR_K_1
    REAL(DP), DIMENSION(nCR,nCR) :: L1, R1
    REAL(DP), DIMENSION(nDOFX,nCR) :: uCR_X, uCR_X_P_X1, uCR_X_N_X1
    REAL(DP), DIMENSION(nDOFX,nCR) :: uCR_M, uCR_M_P_X1, uCR_M_N_X1
    REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: uCR_X1, uCR_X1_T

    IF( nDOF == 1 ) RETURN

    IF( .NOT. ApplySlopeLimiter ) RETURN

    CALL DetectDiscontinuities

    ALLOCATE &
      ( uCR_X1  (1:nCR,1:nNodesE,1:nE,1:nX(1),1:nX(2),1:nX(3)), &
        uCR_X1_T(1:nCR,1:nNodesE,1:nE,1:nX(1),1:nX(2),1:nX(3)) )

    ASSOCIATE( dX1 => MeshX(1) % Width(1:nX(1)) )

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              ! --- Allow limiting if element is tagged as discontinuity

              IF( .NOT. ( Discontinuity(iE,iX1,iX2,iX3) == 1.0_DP ) ) &
                CYCLE

              IF( Debug )THEN
                WRITE(*,*)
                WRITE(*,'(A6,A)') '', 'ApplySlopeLimiter_M1_DG'
                WRITE(*,*)
              END IF

              DO iNodeE = 1, nNodesE

                DO iCR = 1, nCR

                  CALL PackSpatialNodes &
                         ( iNodeE, uCR(:,iE,iX1-1,iX2,iX3,iCR,iS), &
                           uCR_X_P_X1(:,iCR) )
                  CALL PackSpatialNodes &
                         ( iNodeE, uCR(:,iE,iX1,  iX2,iX3,iCR,iS), &
                           uCR_X     (:,iCR) )
                  CALL PackSpatialNodes &
                         ( iNodeE, uCR(:,iE,iX1+1,iX2,iX3,iCR,iS), &
                           uCR_X_N_X1(:,iCR) )

                END DO

                ! --- Limiting Using Modal Spatial Representation ---

                DO iCR = 1, nCR

                  ! --- Map To Modal Representation ---

                  CALL MapNodalToModal_Radiation_X &
                         ( uCR_X_P_X1(:,iCR), uCR_M_P_X1(:,iCR) )
                  CALL MapNodalToModal_Radiation_X &
                         ( uCR_X     (:,iCR), uCR_M     (:,iCR) )
                  CALL MapNodalToModal_Radiation_X &
                         ( uCR_X_N_X1(:,iCR), uCR_M_N_X1(:,iCR) )

                END DO

                ! --- Cell-Averaged Quantities ---

                uCR_A     (1:nCR) = uCR_M     (1,1:nCR)
                uCR_A_P_X1(1:nCR) = uCR_M_P_X1(1,1:nCR)
                uCR_A_N_X1(1:nCR) = uCR_M_N_X1(1,1:nCR)

                ! --- Slope From Modal Representation ---

                CALL ComputeEigenvectors_L &
                       ( uCR_A(iCR_N),  uCR_A(iCR_G1), &
                         uCR_A(iCR_G2), uCR_A(iCR_G3), &
                         L1, Componentwise )

                uCR_X1(1:nCR,iNodeE,iE,iX1,iX2,iX3) &
                  = MATMUL( L1, uCR_M(2,1:nCR) ) ! X1-Dimension

                ! --- Compute Limited Slopes ---

                uCR_X1_T(1:nCR,iNodeE,iE,iX1,iX2,iX3) &
                  = MinModB &
                      ( uCR_X1(1:nCR,iNodeE,iE,iX1,iX2,iX3), &
                        BetaTVD * MATMUL( L1, uCR_A(1:nCR) - uCR_A_P_X1 ), &
                        BetaTVD * MATMUL( L1, uCR_A_N_X1 - uCR_A(1:nCR) ), &
                        dX1(iX1), BetaTVB )

                IF( Debug )THEN
                  PRINT*
                  PRINT*, "iE, iX1, iX2, iX3 = ", iE, iX1, iX2, iX3
                  PRINT*, "  iNodeE   = ", iNodeE
                  PRINT*, "  uCR_X1   = ", &
                    uCR_X1(1:nCR,iNodeE,iE,iX1,iX2,iX3)
                  PRINT*, "  uCR_X1_L = ", &
                    BetaTVD * MATMUL( L1, uCR_A(1:nCR) - uCR_A_P_X1 )
                  PRINT*, "  uCR_X1_R = ", &
                    BetaTVD * MATMUL( L1, uCR_A_N_X1 - uCR_A(1:nCR) )
                  PRINT*, "  uCR_X1_T = ", &
                    uCR_X1_T(1:nCR,iNodeE,iE,iX1,iX2,iX2)
                  PRINT*
                END IF

              END DO

            END DO
          END DO
        END DO
      END DO

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              ! --- Allow limiting if element is tagged as discontinuity

              IF( .NOT. ( Discontinuity(iE,iX1,iX2,iX3) == 1.0_DP ) ) &
                CYCLE

              DO iNodeE = 1, nNodesE

                DO iCR = 1, nCR

                  SlopeDifference(iCR) &
                    = ABS( uCR_X1(iCR,iNodeE,iE,iX1,iX2,iX3) &
                           - uCR_X1_T(iCR,iNodeE,iE,iX1,iX2,iX3) ) &
                      / MAX( ABS( uCR_X1  (iCR,iNodeE,iE,iX1,iX2,iX3) ), &
                             ABS( uCR_X1_T(iCR,iNodeE,iE,iX1,iX2,iX3) ), &
                             TINY(1.0_DP) )

                END DO

                LimitPolynomial = ANY( SlopeDifference(1:nCR) > Tol_TVD )

                IF( LimitPolynomial )THEN

                  IF( Debug )THEN

                    WRITE(*,*)
                    WRITE(*,'(A4,A,1I2.2)') &
                      '', 'Limiting Radiation Field: '
                    WRITE(*,'(A6,A,4I5.4)') '', &
                      'iE, iX1, iX2, iX3 = ', iE, iX1, iX2, iX3
                    WRITE(*,'(A6,A,I5.4)')  '', &
                      'iNodeE = ', iNodeE
                    WRITE(*,*)
                    PRINT*, "  uCR_X1    = ", &
                      uCR_X1  (1:nCR,iNodeE,iE,iX1,iX2,iX3)
                    PRINT*, "  uCR_X1_T  = ", &
                      uCR_X1_T(1:nCR,iNodeE,iE,iX1,iX2,iX3)
                    PRINT*, "  slopeDiff = ", SlopeDifference

                  END IF


                  DO iCR = 1, nCR

                    CALL PackSpatialNodes &
                           ( iNodeE, uCR(:,iE,iX1,iX2,iX3,iCR,iS), &
                             uCR_X(:,iCR) )

                    uCR_K_0(iCR) &
                      = SUM( WeightsGX(:) * uCR_X(:,iCR) &
                               * VolJacX(:,iX1,iX2,iX3) ) / VolX(iX1,iX2,iX3)

                    CALL MapNodalToModal_Radiation_X &
                           ( uCR_X(:,iCR), uCR_M(:,iCR) )

                  END DO

                  ! --- Cell-Averaged Quantities ---

                  uCR_A(1:nCR) = uCR_M(1,1:nCR)

                  ! --- Back to Conserved Variables ---

                  CALL ComputeEigenvectors_R &
                         ( uCR_A(iCR_N),  uCR_A(iCR_G1), &
                           uCR_A(iCR_G2), uCR_A(iCR_G3), &
                           R1, Componentwise )

                  uCR_M(:,1:nCR) = 0.0_DP
                  uCR_M(1,1:nCR) & ! -- Cell-Average
                    = uCR_A(1:nCR)
                  uCR_M(2,1:nCR) & ! -- Slope X1-Direction
                    = MATMUL( R1, uCR_X1_T(1:nCR,iNodeE,iE,iX1,iX2,iX3) )

                  ! --- Correct Coefficients with Constrained Least Squares ---
                  ! --- to Presevre Conserved Quantities (i.e., N,G1,G2,G3) ---

                  A0(1:2,1) = [ 1.0_DP, 0.0_DP ]
                  A0(1:2,2) = [ 0.0_DP, 1.0_DP ]
                  B0(1,1) &
                    = SUM( WeightsGX(:) * LegendreX(:,1) &
                             * VolJacX(:,iX1,iX2,iX3) ) / VolX(iX1,iX2,iX3)
                  B0(1,2) &
                    = SUM( WeightsGX(:) * LegendreX(:,2) &
                             * VolJacX(:,iX1,iX2,iX3) ) / VolX(iX1,iX2,iX3)

                  DO iCR = 1, nCR

                    A = A0
                    B = B0
                    c = uCR_M(1:2,iCR)
                    d = uCR_K_0  (iCR)

                    LWORK = SIZE( WORK )
                    CALL DGGLSE( SIZE( A, 1 ), SIZE( A, 2 ), SIZE( B, 1 ), &
                                 A, SIZE( A, 1 ), B, SIZE( B, 1 ), c, d, &
                                 uCR_M(1:2,iCR), WORK, LWORK, INFO )

                  END DO

                  ! --- Back to Nodal Representation ---

                  DO iCR = 1, nCR

                    CALL MapModalToNodal_Radiation_X &
                           ( uCR_X(:,iCR), uCR_M(:,iCR) )

                    uCR_K_1(iCR) &
                      = SUM( WeightsGX(:) * uCR_X(:,iCR) &
                               * VolJacX(:,iX1,iX2,iX3) ) / VolX(iX1,iX2,iX3)

                    CALL UnpackSpatialNodes &
                           ( iNodeE, uCR(:,iE,iX1,iX2,iX3,iCR,iS), &
                             uCR_X(:,iCR) )

                  END DO

                  IF( Debug )THEN
                    PRINT*
                    PRINT*, "  |duCR| = ", ABS( uCR_K_1(:) - uCR_K_0(:) )
                    PRINT*
                  END IF

                END IF

              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

    END ASSOCIATE ! dX1, etc.

    DEALLOCATE( uCR_X1, uCR_X1_T )

  END SUBROUTINE ApplySlopeLimiter_M1_DG


END MODULE MomentEquationsSlopeLimiterModule_DG
