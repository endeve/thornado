MODULE EulerEquationsPositivityLimiterModule_DG_TABLE

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    Centimeter, &
    Gram, &
    Erg, &
    Kelvin
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
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, nCF, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_T, iAF_E, iAF_Ye
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, Min_T, Min_Y, &
    ComputeSpecificInternalEnergy_TABLE
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

  LOGICAL, PARAMETER :: Debug = .FALSE.

  PUBLIC :: ApplyPositivityLimiter_Euler_DG_TABLE

CONTAINS


  SUBROUTINE ApplyPositivityLimiter_Euler_DG_TABLE

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iPoint, iCF
    REAL(DP) :: Theta, Theta_P
    REAL(DP) :: dD, dS1, dS2, dS3, dE, dNe
    REAL(DP) :: Tol_E
    REAL(DP), DIMENSION(1:nCF)               :: uCF_K
    REAL(DP), DIMENSION(1:nPF)               :: uPF_K, uPF_A
    REAL(DP), DIMENSION(nDOFX,nCF)           :: uCF_M
    REAL(DP), DIMENSION(nPositivePoints)     :: MinNe
    REAL(DP), DIMENSION(nPositivePoints)     :: MinPF_E
    REAL(DP), DIMENSION(nPositivePoints)     :: dMinPF_E
    REAL(DP), DIMENSION(nPositivePoints)     :: MinAF_E
    REAL(DP), DIMENSION(nPositivePoints)     :: dMinAF_EdD
    REAL(DP), DIMENSION(nPositivePoints)     :: dMinAF_EdT
    REAL(DP), DIMENSION(nPositivePoints)     :: dMinAF_EdY
    REAL(DP), DIMENSION(nPositivePoints,nCF) :: uCF_P
    REAL(DP), DIMENSION(nPositivePoints,nPF) :: uPF_P
    REAL(DP), DIMENSION(nPositivePoints,nAF) :: uAF_P

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. ApplyPositivityLimiter ) RETURN

    IF( Debug )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'ApplyPositivityLimiter_Euler_DG_TABLE'
      WRITE(*,*)
    END IF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          ! --- Ensure Positive Mass Density ---

          DO iPoint = 1, nPositivePoints

            uCF_P(iPoint,iCF_D) &
              = DOT_PRODUCT( uCF(:,iX1,iX2,iX3,iCF_D ), Lagrange(:,iPoint) )

          END DO
          uPF_P(:,iPF_D) = uCF_P(:,iCF_D)

          IF( ANY( uPF_P(:,iPF_D) < Min_D ) )THEN

            ! --- Limiting Using Modal Representation ---

            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1,iX2,iX3,iCF_D), uCF_M(:,iCF_D) )

            uCF_K(iCF_D) &
              = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF_D) &
                       * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                  / SUM( WeightsX_q(:) * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

            Theta = 1.0_DP
            DO iPoint = 1, nPositivePoints

              IF( uCF_P(iPoint,iCF_D) < Min_D ) &
                Theta &
                  = MIN( Theta, &
                         ( Min_D - uCF_K(iCF_D) ) &
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

              CALL MapModalToNodal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF_D ), uCF_M(:,iCF_D ) )

            END IF

          END IF

          ! --- Ensure Positive Electron Density ---

          DO iPoint = 1, nPositivePoints

            uCF_P(iPoint,iCF_D) &
              = DOT_PRODUCT( uCF(:,iX1,iX2,iX3,iCF_D ), Lagrange(:,iPoint) )

            uCF_P(iPoint,iCF_Ne) &
              = DOT_PRODUCT( uCF(:,iX1,iX2,iX3,iCF_Ne), Lagrange(:,iPoint) )

          END DO
          uPF_P(:,iPF_D)  = uCF_P(:,iCF_D)
          uPF_P(:,iPF_Ne) = uCF_P(:,iCF_Ne)

          MinNe = 1.000001_DP * Min_Y * uPF_P(:,iPF_D) / AtomicMassUnit

          IF( ANY( uPF_P(:,iPF_Ne) < MinNe(:) ) )THEN

            ! --- Limiting Using Modal Representation ---

            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1,iX2,iX3,iCF_Ne), uCF_M(:,iCF_Ne) )

            uCF_K(iCF_Ne) &
              = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF_Ne) &
                       * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                / SUM( WeightsX_q(:) * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

            Theta = 1.0_DP
            DO iPoint = 1, nPositivePoints

              IF( uCF_P(iPoint,iCF_Ne) < MinNe(iPoint) ) &
                Theta &
                  = MIN( Theta, &
                         ( MinNe(iPoint) - uCF_K(iCF_Ne) ) &
                           / ( uCF_P(iPoint,iCF_Ne) - uCF_K(iCF_Ne) ) )

            END DO

            IF( Theta < 1.0_DP )THEN

              IF( Debug )THEN
                PRINT*
                PRINT*, "iX1,iX2,iX3 = ", iX1, iX2, iX3
                PRINT*, "Theta_2 = ", Theta
                PRINT*, " N_K = ", uCF_K(iCF_Ne)
                PRINT*, " N_p = ", uCF_P(:,iCF_Ne)
                PRINT*
              END IF

              uCF_M(1,iCF_Ne) &
                = Theta * uCF_M(1,iCF_Ne) &
                    + (1.0_DP - Theta) * uCF_K(iCF_Ne)

              uCF_M(2:nDOFX,iCF_Ne) &
                = Theta * uCF_M(2:nDOFX,iCF_Ne)

              CALL MapModalToNodal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF_Ne), uCF_M(:,iCF_Ne) )

            END IF

          END IF

          ! --- Ensure Positive Energy Density ---

          DO iPoint = 1, nPositivePoints

            DO iCF = 1, nCF

              uCF_P(iPoint,iCF) &
                = DOT_PRODUCT( uCF(:,iX1,iX2,iX3,iCF), Lagrange(:,iPoint) )

            END DO
            uPF_P(iPoint,:) = Primitive( uCF_P(iPoint,:) )

          END DO

          uAF_P(:,iAF_T ) = 1.000001_DP * Min_T
          uAF_P(:,iAF_Ye) = AtomicMassUnit * uPF_P(:,iPF_Ne) &
                              / uPF_P(:,iPF_D)

          CALL ComputeSpecificInternalEnergy_TABLE &
                 ( uPF_P(:,iPF_D), uAF_P(:,iAF_T), uAF_P(:,iAF_Ye), &
                   MinAF_E(:), dMinAF_EdD(:), dMinAF_EdT(:), dMinAF_EdY(:) )

          MinPF_E(:) = MinAF_E(:) * uPF_P(:,iPF_D)

          IF( ANY( uPF_P(:,iPF_E) < MinPF_E(:) ) )THEN

            ! --- Limiting Using Modal Representation ---

            DO iCF = 1, nCF

              CALL MapNodalToModal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

              uCF_K(iCF) &
                = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF) &
                         * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                  / SUM( WeightsX_q(:) * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

            END DO

            Theta = 1.0_DP
            DO iPoint = 1, nPositivePoints

              dMinPF_E(iPoint) = 0.0_DP

              IF( uPF_P(iPoint,iPF_E) > MinPF_E(iPoint) ) CYCLE

              dD  = uCF_P(iPoint,iCF_D)  - uCF_K(iCF_D)
              dS1 = uCF_P(iPoint,iCF_S1) - uCF_K(iCF_S1)
              dS2 = uCF_P(iPoint,iCF_S2) - uCF_K(iCF_S2)
              dS3 = uCF_P(iPoint,iCF_S3) - uCF_K(iCF_S3)
              dE  = uCF_P(iPoint,iCF_E)  - uCF_K(iCF_E)
              dNe = uCF_P(iPoint,iCF_Ne) - uCF_K(iCF_Ne)

              uPF_K = Primitive( uCF_K )

              dMinPF_E(iPoint) &
                = - ( MinAF_E(iPoint) &
                      + uPF_P(iPoint,iPF_D ) * dMinAF_EdD(iPoint) &
                      - uAF_P(iPoint,iAF_Ye) * dMinAF_EdY(iPoint) ) * dD &
                  - AtomicMassUnit * dMinAF_EdY(iPoint) * dNe

              Tol_E = MIN( 0.999_DP * uPF_K(iPF_E), &
                           MAX( MinPF_E(iPoint), &
                                MinPF_E(iPoint)+dMinPF_E(iPoint) ) )

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
                PRINT*, "  Theta_3 = ", Theta
                PRINT*, "      E_K = ", uCF_K  (iCF_E) / (Erg/Centimeter**3)
                PRINT*, "      E_p = ", uCF_P(:,iCF_E) / (Erg/Centimeter**3)
                PRINT*, "      e_K = ", uPF_K  (iPF_E) / (Erg/Centimeter**3)
                PRINT*, "      e_p = ", uPF_P(:,iPF_E) / (Erg/Centimeter**3)
                PRINT*, "  MinPF_E = ", MinPF_E (:)    / (Erg/Centimeter**3)
                PRINT*, " dMinPF_E = ", dMinPF_E(:)    / (Erg/Centimeter**3)
                PRINT*, "    Tol_E = ", &
                  MIN( 0.999_DP * uPF_K(iPF_E), &
                       MAX( MinPF_E(:), MinPF_E(:)+dMinPF_E(:) ) ) &
                    / (Erg/Centimeter**3)
                PRINT*
              END IF

              DO iCF = 1, nCF

                uCF_M(1,iCF) &
                  = Theta * uCF_M(1,iCF) &
                      + (1.0_DP - Theta) * uCF_K(iCF)

                uCF_M(2:nDOFX,iCF) &
                  = Theta * uCF_M(2:nDOFX,iCF)

                CALL MapModalToNodal_Fluid &
                       ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

              END DO

            END IF

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

          uAF_P(:,iAF_T ) = 1.000001_DP * Min_T
          uAF_P(:,iAF_Ye) = AtomicMassUnit * uPF_P(:,iPF_Ne) &
                              / uPF_P(:,iPF_D)

          CALL ComputeSpecificInternalEnergy_TABLE &
                 ( uPF_P(:,iPF_D), uAF_P(:,iAF_T), uAF_P(:,iAF_Ye), &
                   MinAF_E(:) )

          MinPF_E(:) = MinAF_E(:) * uPF_P(:,iPF_D)

          IF( ANY( uPF_P(:,iPF_D) < Min_D ) &
                .OR. ANY( uPF_P(:,iPF_E) < MinPF_E(:) ) )THEN

            PRINT*
            PRINT*, "Problem with Positivity Limiter!"
            PRINT*, "  iX1, iX2, iX3 = ", iX1, iX2, iX3
            PRINT*, "  Theta = ", Theta
            PRINT*, "  Min_D  = ", Min_D / (Gram/Centimeter**3)
            PRINT*, "  MinPF_E  = ", MinPF_E(:) / (Erg/Centimeter**3)
            PRINT*
            PRINT*, "  Conserved Fields (Nodal):"
            PRINT*, "  D_N  = ", uCF_P(:,iCF_D) / (Gram/Centimeter**3)
            PRINT*, "  S1_N = ", uCF_P(:,iCF_S1)
            PRINT*, "  S2_N = ", uCF_P(:,iCF_S2)
            PRINT*, "  S3_N = ", uCF_P(:,iCF_S3)
            PRINT*, "  E_N  = ", uCF_P(:,iCF_E) / (Erg/Centimeter**3)
            PRINT*, "  Ne_N = ", uCF_P(:,iCF_Ne)
            PRINT*
            PRINT*, "  Primitive Fields (Nodal):"
            PRINT*, "  D_N  = ", uPF_P(:,iPF_D) / (Gram/Centimeter**3)
            PRINT*, "  V1_N = ", uPF_P(:,iPF_V1)
            PRINT*, "  V2_N = ", uPF_P(:,iPF_V2)
            PRINT*, "  V3_N = ", uPF_P(:,iPF_V3)
            PRINT*, "  E_N  = ", uPF_P(:,iPF_E) / (Erg/Centimeter**3)
            PRINT*, "  Ne_N = ", uPF_P(:,iPF_Ne)
            PRINT*
            PRINT*, "  Cell-Averages (Conserved):"
            PRINT*, "  D_A  = ", uCF_M(1,iCF_D) / (Gram/Centimeter**3)
            PRINT*, "  S1_A = ", uCF_M(1,iCF_S1)
            PRINT*, "  S2_A = ", uCF_M(1,iCF_S2)
            PRINT*, "  S3_A = ", uCF_M(1,iCF_S3)
            PRINT*, "  E_A  = ", uCF_M(1,iCF_E) / (Erg/Centimeter**3)
            PRINT*, "  Ne_A = ", uCF_M(1,iCF_Ne)
            PRINT*
            PRINT*, "  Cell-Averages (Primitive):"
            PRINT*, "  D_A  = ", uPF_A(iPF_D) / (Gram/Centimeter**3)
            PRINT*, "  V1_A = ", uPF_A(iPF_V1)
            PRINT*, "  V2_A = ", uPF_A(iPF_V2)
            PRINT*, "  V3_A = ", uPF_A(iPF_V3)
            PRINT*, "  E_A  = ", uPF_A(iPF_E) / (Erg/Centimeter**3)
            PRINT*, "  Ne_A = ", uPF_A(iPF_Ne)
            PRINT*

            STOP

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE ApplyPositivityLimiter_Euler_DG_TABLE


END MODULE EulerEquationsPositivityLimiterModule_DG_TABLE
