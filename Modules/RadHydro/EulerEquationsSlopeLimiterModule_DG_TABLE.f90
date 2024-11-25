MODULE EulerEquationsSlopeLimiterModule_DG_TABLE

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE UtilitiesModule, ONLY: &
    MinModB
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_SqrtGm
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, iAF_Gm, iAF_Cs, &
    Shock
  USE EquationOfStateModule, ONLY: &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    Primitive
  USE EulerEquationsUtilitiesModule_TABLE, ONLY: &
    ComputeEigenvectors_L, &
    ComputeEigenvectors_R
  USE EulerEquationsLimiterUtilitiesModule_DG, ONLY: &
    DetectDiscontinuities
  USE EulerEquationsLimiterModule_DG, ONLY: &
    ApplySlopeLimiter, &
    BetaTVD, &
    BetaTVB, &
    Legendre

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PARAMETER :: Debug         = .FALSE.
  LOGICAL,  PARAMETER :: Componentwise = .FALSE.
  REAL(DP), PARAMETER :: Tol_TVD       = 1.0d-2

  PUBLIC :: ApplySlopeLimiter_Euler_DG_TABLE

CONTAINS


  SUBROUTINE ApplySlopeLimiter_Euler_DG_TABLE

    LOGICAL :: LimitPolynomial
    INTEGER :: iX1, iX2, iX3, iCF
    INTEGER :: LWORK, INFO
    REAL(DP), DIMENSION(1)   :: d
    REAL(DP), DIMENSION(2)   :: c
    REAL(DP), DIMENSION(5)   :: WORK
    REAL(DP), DIMENSION(2,2) :: A0, A
    REAL(DP), DIMENSION(1,2) :: B0, B
    REAL(DP), DIMENSION(nCF) :: uCF_A, uCF_A_P_X1, uCF_A_N_X1
    REAL(DP), DIMENSION(nCF) :: SlopeDifference
    REAL(DP), DIMENSION(nCF) :: Lambda
    REAL(DP), DIMENSION(nCF) :: uCF_K_0, uCF_K_1
    REAL(DP), DIMENSION(nPF) :: uPF_A
    REAL(DP), DIMENSION(nAF) :: uAF_A
    REAL(DP), DIMENSION(nCF,nCF) :: L1, dFdU, R1
    REAL(DP), DIMENSION(nDOFX,nCF) :: uCF_M, uCF_M_P_X1, uCF_M_N_X1
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: uCF_X1, uCF_X1_T

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. ApplySlopeLimiter ) RETURN

    CALL DetectDiscontinuities

    IF( Debug )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'ApplySlopeLimiter_Euler_DG_TABLE'
      WRITE(*,*)
    END IF

    ALLOCATE &
      ( uCF_X1  (1:nCF,1:nX(1),1:nX(2),1:nX(3)), &
        uCF_X1_T(1:nCF,1:nX(1),1:nX(2),1:nX(3)) )

    ASSOCIATE( dX1 => MeshX(1) % Width(1:nX(1)) )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          ! --- Allow limiting if element marked as shock ---

          IF( .NOT. Shock(iX1,iX2,iX3) == 1.0 ) CYCLE

          ! --- Limiting Using Modal Representation ---

          DO iCF = 1, nCF

            ! --- Map To Modal Representation ---

            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1-1,iX2,iX3,iCF), uCF_M_P_X1(:,iCF) )
            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1,  iX2,iX3,iCF), uCF_M     (:,iCF) )
            CALL MapNodalToModal_Fluid &
                   ( uCF(:,iX1+1,iX2,iX3,iCF), uCF_M_N_X1(:,iCF) )

          END DO

          ! --- Cell-Averaged Quantities ---

          uCF_A     (1:nCF) = uCF_M     (1,1:nCF)
          uCF_A_P_X1(1:nCF) = uCF_M_P_X1(1,1:nCF)
          uCF_A_N_X1(1:nCF) = uCF_M_N_X1(1,1:nCF)

          uPF_A(1:nPF) = Primitive( uCF_A(1:nCF) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF_A(iPF_D ), uPF_A(iPF_E), uPF_A(iPF_Ne), &
                   uAF_A(iAF_P ), uAF_A(iAF_T), uAF_A(iAF_Ye), &
                   uAF_A(iAF_S ), uAF_A(iAF_E), uAF_A(iAF_Gm), &
                   uAF_A(iAF_Cs) )

          ! --- Slope From Modal Representation ---

          CALL ComputeEigenvectors_L &
                 ( [ uPF_A(iPF_D ) ], [ uAF_A(iAF_T ) ], [ uAF_A(iAF_Ye) ], &
                   [ uPF_A(iPF_V1) ], [ uPF_A(iPF_V2) ], [ uPF_A(iPF_V3) ], &
                   Lambda, L1, dFdU, Componentwise )

          uCF_X1(1:nCF,iX1,iX2,iX3) &
            = MATMUL( L1, uCF_M(2,1:nCF) ) ! X1-Dimension

          uCF_X1_T(1:nCF,iX1,iX2,iX2) &
            = MinModB &
                ( uCF_X1(1:nCF,iX1,iX2,iX3), &
                  BetaTVD * MATMUL( L1, uCF_A(1:nCF) - uCF_A_P_X1 ), &
                  BetaTVD * MATMUL( L1, uCF_A_N_X1 - uCF_A(1:nCF) ), &
                  dX1(iX1), BetaTVB )

          IF( Debug )THEN
            PRINT*
            PRINT*, "iX1, iX2, iX3 = ", iX1, iX2, iX3
            PRINT*, "  uCF_X1   = ", uCF_X1(1:nCF,iX1,iX2,iX3)
            PRINT*, "  uCF_X1_L = ", &
              BetaTVD * MATMUL( L1, uCF_A(1:nCF) - uCF_A_P_X1 )
            PRINT*, "  uCF_X1_R = ", &
              BetaTVD * MATMUL( L1, uCF_A_N_X1 - uCF_A(1:nCF) )
            PRINT*, "  uCF_X1_T = ", uCF_X1_T(1:nCF,iX1,iX2,iX2)
            PRINT*
          END IF

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc. 

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          ! --- Allow limiting if element marked as shock ---

          IF( .NOT. Shock(iX1,iX2,iX3) == 1.0 ) CYCLE

          DO iCF = 1, nCF

            SlopeDifference(iCF) &
              = ABS( uCF_X1(iCF,iX1,iX2,iX3) - uCF_X1_T(iCF,iX1,iX2,iX3) ) &
                / MAX( ABS( uCF_X1  (iCF,iX1,iX2,iX3) ), &
                       ABS( uCF_X1_T(iCF,iX1,iX2,iX3) ), TINY(1.0_DP) )

          END DO

          LimitPolynomial = ANY( SlopeDifference(1:nCF) > Tol_TVD )

          IF( LimitPolynomial )THEN

            IF( Debug )THEN

              PRINT*
              PRINT*, "iX1, iX2, iX3 = ", iX1, iX2, iX3
              PRINT*, "  uCF_X1    = ", uCF_X1  (1:nCF,iX1,iX2,iX3)
              PRINT*, "  uCF_X1_T  = ", uCF_X1_T(1:nCF,iX1,iX2,iX3)
              PRINT*, "  slopeDiff = ", SlopeDifference
              PRINT*, "  Legendre(:,1) = ", Legendre(:,1)
              PRINT*, "  Legendre(:,2) = ", Legendre(:,2)
              PRINT*

            END IF

            DO iCF = 1, nCF

              uCF_K_0(iCF) &
                = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF) &
                         * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

              CALL MapNodalToModal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

            END DO

            ! --- Cell-Averaged Quantities ---

            uCF_A(1:nCF) = uCF_M(1,1:nCF)
            uPF_A(1:nPF) = Primitive( uCF_A(1:nCF) )

            CALL ComputeAuxiliary_Fluid &
                   ( uPF_A(iPF_D ), uPF_A(iPF_E), uPF_A(iPF_Ne), &
                     uAF_A(iAF_P ), uAF_A(iAF_T), uAF_A(iAF_Ye), &
                     uAF_A(iAF_S ), uAF_A(iAF_E), uAF_A(iAF_Gm), &
                     uAF_A(iAF_Cs) )

            ! --- Back to Conserved Variables ---

            CALL ComputeEigenvectors_R &
                 ( [ uPF_A(iPF_D ) ], [ uAF_A(iAF_T ) ], [ uAF_A(iAF_Ye) ], &
                   [ uPF_A(iPF_V1) ], [ uPF_A(iPF_V2) ], [ uPF_A(iPF_V3) ], &
                   Lambda, R1, dFdU, Componentwise )

            uCF_M(:,1:nCF) = 0.0_DP
            uCF_M(1,1:nCF) & ! -- Cell-Average
              = uCF_A(1:nCF)
            uCF_M(2,1:nCF) & ! -- Slope X1-Direction
              = MATMUL( R1, uCF_X1_T(1:nCF,iX1,iX2,iX3) )

            ! --- Correct Modal Coefficients with Constrained Least Squares ---
            ! --- to Presevre Conserved Quantities (i.e., D,S1,S2,S3,E,Ne)  ---

            A0(1:2,1) = [ 1.0_DP, 0.0_DP ]
            A0(1:2,2) = [ 0.0_DP, 1.0_DP ]
            B0(1,1) = SUM( WeightsX_q(:) * Legendre(:,1) &
                             * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )
            B0(1,2) = SUM( WeightsX_q(:) * Legendre(:,2) &
                             * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )
            DO iCF = 1, nCF

              A = A0
              B = B0
              c = uCF_M(1:2,iCF)
              d = uCF_K_0  (iCF)

              LWORK = SIZE( WORK )
              CALL DGGLSE( SIZE( A, 1 ), SIZE( A, 2 ), SIZE( B, 1 ), &
                           A, SIZE( A, 1 ), B, SIZE( B, 1 ), c, d, &
                           uCF_M(1:2,iCF), WORK, LWORK, INFO )

            END DO

            ! --- Back to Nodal Representation ---

            DO iCF = 1, nCF

              CALL MapModalToNodal_Fluid &
                     ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

              uCF_K_1(iCF) &
                = SUM( WeightsX_q(:) * uCF(:,iX1,iX2,iX3,iCF) &
                         * uGF(:,iX1,iX2,iX3,iGF_SqrtGm) )

            END DO

            IF( Debug )THEN
              PRINT*
              PRINT*, "  |duCF| = ", ABS( uCF_K_1(:) - uCF_K_0(:) )
              PRINT*
            END IF

          END IF

        END DO
      END DO
    END DO

    DEALLOCATE( uCF_X1, uCF_X1_T )

    IF( Debug )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'ApplySlopeLimiter_Euler_DG Done'
      WRITE(*,*)
    END IF

  END SUBROUTINE ApplySlopeLimiter_Euler_DG_TABLE


END MODULE EulerEquationsSlopeLimiterModule_DG_TABLE
