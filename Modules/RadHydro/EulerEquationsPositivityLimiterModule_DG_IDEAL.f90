MODULE EulerEquationsPositivityLimiterModule_DG_IDEAL

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE PolynomialBasisModule_Legendre, ONLY: &
    evalPX
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, nCF, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, &
    nAF
  USE EulerEquationsUtilitiesModule, ONLY: &
    Primitive
  USE EulerEquationsLimiterUtilitiesModule_DG, ONLY: &
    SolveTheta_Bisection
  USE EulerEquationsLimiterModule_DG, ONLY: &
    ApplyPositivityLimiter, &
    nPositivePoints, &
    Points_X1, Points_X2, Points_X3, &
    Lagrange

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PARAMETER :: Debug = .FALSE.
  REAL(DP), PARAMETER :: Tol_D = 1.0d-12
  REAL(DP), PARAMETER :: Tol_E = 1.0d-12

  PUBLIC :: ApplyPositivityLimiter_Euler_DG_IDEAL

CONTAINS


  SUBROUTINE ApplyPositivityLimiter_Euler_DG_IDEAL

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iPoint, iCF
    REAL(DP) :: Theta, Theta_P, dD, dS1, dS2, dS3, dE
    REAL(DP), DIMENSION(1:nCF)               :: uCF_K, uCF_K_0, uCF_K_1
    REAL(DP), DIMENSION(1:nPF)               :: uPF_A
    REAL(DP), DIMENSION(nDOFX,nCF)           :: uCF_M
    REAL(DP), DIMENSION(nPositivePoints,nCF) :: uCF_P
    REAL(DP), DIMENSION(nPositivePoints,nPF) :: uPF_P
    REAL(DP), DIMENSION(nPositivePoints,nAF) :: uAF_P

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. ApplyPositivityLimiter ) RETURN

    IF( Debug )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'ApplyPositivityLimiter_Euler_DG_IDEAL'
      WRITE(*,*)
    END IF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iPoint = 1, nPositivePoints
            DO iCF = 1, nCF

              uCF_P(iPoint,iCF) &
                = DOT_PRODUCT &
                    ( uCF(:,iX1,iX2,iX3,iCF), Lagrange(:,iPoint) )

            END DO

            uPF_P(iPoint,:) = Primitive( uCF_P(iPoint,:) )

          END DO

          IF( ANY( uPF_P(:,iPF_D) < Tol_D ) &
                .OR. ANY( uPF_P(:,iPF_E) < Tol_E ) )THEN

            ! --- Limiting Using Modal Representation ---

            DO iCF = 1, nCF

              CALL MapNodalToModal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

              uCF_K(iCF) &
                = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF) &
                         * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                  / SUM( WeightsX_q(:) * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

              uCF_K_0(iCF) &
                = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF) &
                         * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

            END DO

            ! --- Ensure Positive Mass Density ---

            IF( ANY( uPF_P(:,iPF_D) < Tol_D ) )THEN

              Theta = 1.0_DP
              DO iPoint = 1, nPositivePoints
                IF( uCF_P(iPoint,iCF_D) < Tol_D ) &
                  Theta &
                    = MIN( Theta, &
                           ( Tol_D - uCF_K(iCF_D) ) &
                           / ( uCF_P(iPoint,iCF_D) - uCF_K(iCF_D) ) )
              END DO

              IF( Theta < 1.0_DP )THEN

                IF( Debug )THEN
                  PRINT*
                  PRINT*, "iX1,iX2,iX3 = ", iX1, iX2, iX3
                  PRINT*, "Theta_1 = ", Theta
                  PRINT*, " D_K = ", uCF_K(iCF_D)
                  PRINT*, " D_p = ", uCF_P(:,iCF_D)
                  PRINT*
                END IF

                uCF_M(1,iCF_D) &
                  = Theta * uCF_M(1,iCF_D) &
                    + (1.0_DP - Theta) * uCF_K(iCF_D)

                uCF_M(2:nDOFX,iCF_D) &
                  = Theta * uCF_M(2:nDOFX,iCF_D)

                DO iPoint = 1, nPositivePoints
                  DO iCF = 1, nCF
                    uCF_P(iPoint,iCF) &
                      = evalPX( uCF_M(:,iCF), Points_X1(iPoint), &
                                Points_X2(iPoint), Points_X3(iPoint) )
                  END DO

                  uPF_P(iPoint,:) = Primitive( uCF_P(iPoint,:) )

                END DO

              END IF

            END IF

            ! --- Ensure Positive Energy Density ---

            IF( ANY( uPF_P(:,iPF_E) < Tol_E ) )THEN

              Theta = 1.0_DP
              DO iPoint = 1, nPositivePoints

                IF( uPF_P(iPoint,iPF_E) >= Tol_E ) CYCLE

                dD  = uCF_P(iPoint,iCF_D)  - uCF_K(iCF_D)
                dS1 = uCF_P(iPoint,iCF_S1) - uCF_K(iCF_S1)
                dS2 = uCF_P(iPoint,iCF_S2) - uCF_K(iCF_S2)
                dS3 = uCF_P(iPoint,iCF_S3) - uCF_K(iCF_S3)
                dE  = uCF_P(iPoint,iCF_E)  - uCF_K(iCF_E)

                CALL SolveTheta_Bisection &
                       ( uCF_K(iCF_D),  uCF_K(iCF_S1), uCF_K(iCF_S2), &
                         uCF_K(iCF_S3), uCF_K(iCF_E), dD, dS1, dS2, dS3, dE, &
                         Tol_E, Theta_P )

                Theta = MIN( Theta, Theta_P )

              END DO

              IF( Theta < 1.0_DP )THEN

                IF( Debug )THEN
                  PRINT*
                  PRINT*, "iX1,iX2,iX3 = ", iX1, iX2, iX3
                  PRINT*, "Theta_2 = ", Theta
                  PRINT*, "Tol_E   = ", Tol_E
                  PRINT*, " E_K = ", uCF_K(iCF_E)
                  PRINT*, " E_p = ", uCF_P(:,iCF_E)
                  PRINT*
                END IF

                DO iCF = 1, nCF

                  uCF_M(1,iCF) &
                    = Theta * uCF_M(1,iCF) &
                        + (1.0_DP - Theta) * uCF_K(iCF)

                  uCF_M(2:nDOFX,iCF) &
                    = Theta * uCF_M(2:nDOFX,iCF)

                END DO

              END IF

            END IF

            ! --- Back to Nodal Representation ---

            DO iCF = 1, nCF

              CALL MapModalToNodal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

              uCF_K_1(iCF) &
                = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF) &
                         * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

            END DO

            ! --- Check Conservation ---

            IF( Debug )THEN
              PRINT*
              PRINT*, '|duCF| = ', ABS( uCF_K_1 - uCF_K_0 )
              PRINT*
              IF( ANY( ABS( uCF_K_1 - uCF_K_0 ) > 1.0d-10 ) ) &
                STOP
            END IF

            ! --- Check Positivity ---

            uPF_A(1:nPF) = Primitive( uCF_M(1,1:nCF) )

            DO iPoint = 1, nPositivePoints
              DO iCF = 1, nCF

                uCF_P(iPoint,iCF) &
                  = DOT_PRODUCT &
                      ( uCF(:,iX1,iX2,iX3,iCF), Lagrange(:,iPoint) )

              END DO

              uPF_P(iPoint,:) = Primitive( uCF_P(iPoint,:) )

            END DO

            IF( ANY( uPF_P(:,iPF_D) <= 0.0_DP ) &
                  .OR. ANY( uPF_P(:,iPF_E) <= 0.0_DP ) )THEN

              PRINT*
              PRINT*, "Problem with Positivity Limiter!"
              PRINT*, "  iX1, iX2, iX3 = ", iX1, iX2, iX3
              PRINT*, "  Theta = ", Theta
              PRINT*
              PRINT*, "  Conserved Fields (Nodal):"
              PRINT*, "  D_N  = ", uCF_P(:,iCF_D)
              PRINT*, "  S1_N = ", uCF_P(:,iCF_S1)
              PRINT*, "  S2_N = ", uCF_P(:,iCF_S2)
              PRINT*, "  S3_N = ", uCF_P(:,iCF_S3)
              PRINT*, "  E_N  = ", uCF_P(:,iCF_E)
              PRINT*
              PRINT*, "  Primitive Fields (Nodal):"
              PRINT*, "  D_N  = ", uPF_P(:,iPF_D)
              PRINT*, "  V1_N = ", uPF_P(:,iPF_V1)
              PRINT*, "  V2_N = ", uPF_P(:,iPF_V2)
              PRINT*, "  V3_N = ", uPF_P(:,iPF_V3)
              PRINT*, "  E_N  = ", uPF_P(:,iPF_E)
              PRINT*
              PRINT*, "  Cell-Averages (Conserved):"
              PRINT*, "  D_A  = ", uCF_M(1,iCF_D)
              PRINT*, "  S1_A = ", uCF_M(1,iCF_S1)
              PRINT*, "  S2_A = ", uCF_M(1,iCF_S2)
              PRINT*, "  S3_A = ", uCF_M(1,iCF_S3)
              PRINT*, "  E_A  = ", uCF_M(1,iCF_E)
              PRINT*
              PRINT*, "  Cell-Averages (Primitive):"
              PRINT*, "  D_A  = ", uPF_A(iPF_D)
              PRINT*, "  V1_A = ", uPF_A(iPF_V1)
              PRINT*, "  V2_A = ", uPF_A(iPF_V2)
              PRINT*, "  V3_A = ", uPF_A(iPF_V3)
              PRINT*, "  E_A  = ", uPF_A(iPF_E)
              PRINT*

              STOP

            END IF

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE ApplyPositivityLimiter_Euler_DG_IDEAL


END MODULE EulerEquationsPositivityLimiterModule_DG_IDEAL
