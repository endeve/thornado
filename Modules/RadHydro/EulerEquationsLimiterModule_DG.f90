MODULE EulerEquationsLimiterModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX
  USE UtilitiesModule, ONLY: &
    GetRoots_Quadratic
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalLX
  USE PolynomialBasisModule_Legendre, ONLY: &
    evalPX
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal, &
    MapModalToNodal
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, nCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, nPF, &
    iAF_P, iAF_Cs, nAF
  USE EquationOfStateModule, ONLY: &
    Auxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    Primitive, &
    ComputeEigenvectors_L, &
    ComputeEigenvectors_R

  IMPLICIT NONE
  PRIVATE

  INTEGER                               :: nPoints
  REAL(DP), PARAMETER                   :: BetaTVD = 1.00_DP
  REAL(DP), PARAMETER                   :: BetaTVB = 50.0_DP
  REAL(DP), PARAMETER                   :: Tol_TVD = 2.0d-2
  REAL(DP), PARAMETER                   :: Tol_D = 1.0d-13
  REAL(DP), PARAMETER                   :: Tol_E = 1.0d-13
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_X1
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_X2
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_X3
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: uCF_P, uPF_P, uAF_P
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: uCF_M, uCF_N

  PUBLIC :: InitializeLimiters_Euler_DG
  PUBLIC :: ApplySlopeLimiter_Euler_DG
  PUBLIC :: ApplyPositivityLimiter_Euler_DG

CONTAINS


  SUBROUTINE InitializeLimiters_Euler_DG

    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iPoint
    REAL(DP), DIMENSION(:), ALLOCATABLE :: NodesX1

    ALLOCATE( NodesX1(nNodesX(1)+2) )

    NodesX1 = [ - 0.5_DP, MeshX(1) % Nodes, + 0.5_DP ]

    nPoints = ( nNodesX(1) + 2 ) * nNodesX(2) * nNodesX(3)

    IF( nNodesX(2) > 1 ) &
      nPoints = nPoints + 2 * nNodesX(1) * nNodesX(3)

    IF( nNodesX(3) > 1 ) &
      nPoints = nPoints + 2 * nNodesX(1) * nNodesX(2)

    ALLOCATE( uCF_P(nPoints,nCF), uPF_P(nPoints,nPF), uAF_P(nPoints,nAF) )

    ALLOCATE( uCF_M(nDOFX,nCF), uCF_N(nDOFX,nCF) )

    ! --- Coordinates of Points Where Positivity is Required:

    ALLOCATE( Points_X1(nPoints), Points_X2(nPoints), Points_X3(nPoints) )

    iPoint = 0

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, SIZE( NodesX1 )

          iPoint = iPoint + 1

          Points_X1(iPoint) = NodesX1(iNodeX1)
          Points_X2(iPoint) = MeshX(2) % Nodes(iNodeX2)
          Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

        END DO
      END DO
    END DO

    IF( nNodesX(2) > 1 )THEN

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)

          iPoint = iPoint + 1

          Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
          Points_X2(iPoint) = - 0.5_DP
          Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

        END DO
      END DO

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)

          iPoint = iPoint + 1

          Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
          Points_X2(iPoint) = + 0.5_DP
          Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

        END DO
      END DO

    END IF

    IF( nNodesX(3) > 1 )THEN

      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iPoint = iPoint + 1

          Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
          Points_X2(iPoint) = MeshX(2) % Nodes(iNodeX2)
          Points_X3(iPoint) = - 0.5_DP

        END DO
      END DO

      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iPoint = iPoint + 1

          Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
          Points_X2(iPoint) = MeshX(1) % Nodes(iNodeX1)
          Points_X3(iPoint) = + 0.5_DP

        END DO
      END DO

    END IF

    DEALLOCATE( NodesX1 )

  END SUBROUTINE InitializeLimiters_Euler_DG


  SUBROUTINE ApplySlopeLimiter_Euler_DG

    LOGICAL                  :: LimitPolynomial
    INTEGER                  :: iX1, iX2, iX3, iCF, i
    REAL(DP), DIMENSION(nCF) :: uCF_A, uCF_A_P, uCF_A_N
    REAL(DP), DIMENSION(nPF) :: uPF_A
    REAL(DP), DIMENSION(nAF) :: uAF_A
    REAL(DP), DIMENSION(nCF,nCF) :: L1, R1
    REAL(DP), DIMENSION(nDOFX,nCF) :: uCF_M_P, uCF_M_N
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      uCF_1, uCF_1_T

    IF( nDOFX == 1 ) RETURN

    ALLOCATE( uCF_1  (nX(1),nX(2),nX(3),nCF), &
              uCF_1_T(nX(1),nX(2),nX(3),nCF) )

    ASSOCIATE( dX1 => MeshX(1) % Width(1:nX(1)) )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          ! --- Limiting Using Modal Representation ---

          DO iCF = 1, nCF

            CALL MapNodalToModal &
                   ( uCF(:,iX1-1,iX2,iX3,iCF), uCF_M_P(:,iCF) )
            CALL MapNodalToModal &
                   ( uCF(:,iX1,  iX2,iX3,iCF), uCF_M  (:,iCF) )
            CALL MapNodalToModal &
                   ( uCF(:,iX1+1,iX2,iX3,iCF), uCF_M_N(:,iCF) )

          END DO

          ! --- Cell-Averaged Quantities ---

          uCF_A = uCF_M(1,:)
          uPF_A = Primitive( uCF_A )
          uAF_A = Auxiliary_Fluid( uPF_A )

          uCF_A_P = uCF_M_P(1,:)
          uCF_A_N = uCF_M_N(1,:)

          ! --- Slope From Modal Representation ---

          CALL ComputeEigenvectors_L &
                 ( uPF_A(iPF_V1), uPF_A(iPF_V2), uPF_A(iPF_V3), &
                   uPF_A(iPF_E), uAF_A(iAF_P), uAF_A(iAF_Cs), L1 )

          uCF_1(iX1,iX2,iX3,1:nCF) &
            = MATMUL( L1, uCF_M(2,1:nCF) ) ! X1-Dimension

          uCF_1_T(iX1,iX2,iX2,1:nCF) &
            = MinModB &
                ( uCF_1(iX1,iX2,iX3,1:nCF), &
                  BetaTVD * MATMUL( L1, ( uCF_A - uCF_A_P ) ), &
                  BetaTVD * MATMUL( L1, ( uCF_A_N - uCF_A ) ), dX1(iX1) )

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc. 

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          LimitPolynomial = .FALSE.
          LimitPolynomial &
            = ANY( ABS( uCF_1(iX1,iX2,iX3,:) &
                          - uCF_1_T(iX1,iX2,iX3,:) ) > Tol_TVD )

          IF( LimitPolynomial )THEN

            DO iCF = 1, nCF

              CALL MapNodalToModal &
                     ( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

            END DO

            ! --- Cell-Averaged Quantities ---

            uCF_A = uCF_M(1,:)
            uPF_A = Primitive( uCF_A )
            uAF_A = Auxiliary_Fluid( uPF_A )

            ! --- Back to Conserved Variables ---

            CALL ComputeEigenvectors_R &
                   ( uPF_A(iPF_V1), uPF_A(iPF_V2), uPF_A(iPF_V3), &
                     uPF_A(iPF_E), uAF_A(iAF_P), uAF_A(iAF_Cs), R1 )

            uCF_M(:,1:nCF) = 0.0_DP
            uCF_M(1,1:nCF) & ! -- Cell-Average
              = uCF_A
            uCF_M(2,1:nCF) & ! -- Slope X1-Direction
              = MATMUL( R1, uCF_1_T(iX1,iX2,iX3,1:nCF) )

            ! --- Back to Nodal Representation ---

            DO iCF = 1, nCF

              CALL MapModalToNodal( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

            END DO

          END IF

        END DO
      END DO
    END DO

    DEALLOCATE( uCF_1, uCF_1_T )

  END SUBROUTINE ApplySlopeLimiter_Euler_DG


  REAL(DP) PURE ELEMENTAL FUNCTION MinMod2( a, b )

    REAL(DP), INTENT(in) :: a, b

    IF( a * b > 0.0_DP )THEN
      IF( ABS( a ) < ABS( b ) )THEN
        MinMod2 = a
      ELSE
        MinMod2 = b
      END IF
    ELSE
      MinMod2 = 0.0_DP
    END IF

    RETURN
  END FUNCTION MinMod2


  REAL(DP) PURE ELEMENTAL FUNCTION MinMod( a, b, c )

    REAL(DP), INTENT(in) :: a, b, c

    MinMod = MinMod2( a, MinMod2( b, c ) )

    RETURN
  END FUNCTION MinMod


  REAL(DP) PURE ELEMENTAL FUNCTION MinModB( a, b, c, dx )

    REAL(DP), INTENT(in) :: a, b, c, dx

    IF( ABS( a ) < BetaTVB * dx**2 )THEN

      MinModB = a

    ELSE

      MinModB = MinMod( a, b, c )

    END IF

  END FUNCTION MinModB


  SUBROUTINE ApplyPositivityLimiter_Euler_DG

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iPoint, iCF, iNode
    REAL(DP) :: Theta, dD, dS1, dS2, dS3, dE
    REAL(DP) :: a, b, c, r1, r2
    REAL(DP), DIMENSION(1:nPF) :: uPF_A

    IF( nDOFX == 1 ) RETURN

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iPoint = 1, nPoints
            DO iCF = 1, nCF

              uCF_P(iPoint,iCF) &
                = evalLX( uCF(:,iX1,iX2,iX3,iCF), Points_X1(iPoint), &
                          Points_X2(iPoint), Points_X3(iPoint) )

            END DO

            uPF_P(iPoint,:) = Primitive( uCF_P(iPoint,:) )

          END DO

          IF( ANY( uPF_P(:,iPF_D) < Tol_D ) &
                .OR. ANY( uPF_P(:,iPF_E) < Tol_E ) )THEN

            ! --- Limiting Using Modal Representation ---

            DO iCF = 1, nCF

              CALL MapNodalToModal( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

            END DO

            ! --- Ensure Positive Mass Density ---

            Theta = 1.0_DP
            DO iPoint = 1, nPoints
              IF( uCF_P(iPoint,iCF_D) < Tol_D ) &
                Theta &
                  = MIN( Theta, &
                         ( Tol_D - uCF_M(1,iCF_D) ) &
                           / ( uCF_P(iPoint,iCF_D) - uCF_M(1,iCF_D) ) )
            END DO

            IF( Theta < 1.0_DP )THEN

              uCF_M(2:nDOFX,iCF_D) &
                = Theta * uCF_M(2:nDOFX,iCF_D)

              DO iPoint = 1, nPoints
                DO iCF = 1, nCF
                  uCF_P(iPoint,iCF) &
                    = evalPX( uCF_M(:,iCF), Points_X1(iPoint), &
                              Points_X2(iPoint), Points_X3(iPoint) )
                END DO

                uPF_P(iPoint,:) = Primitive( uCF_P(iPoint,:) )

              END DO

            END IF

            ! --- Ensure Positive Energy Density ---

            IF( ANY( uPF_P(:,iPF_E) < Tol_E ) )THEN

              Theta = 1.0_DP
              DO iPoint = 1, nPoints

                dD  = uCF_P(iPoint,iCF_D)  - uCF_M(1,iCF_D)
                dS1 = uCF_P(iPoint,iCF_S1) - uCF_M(1,iCF_S1)
                dS2 = uCF_P(iPoint,iCF_S2) - uCF_M(1,iCF_S2)
                dS3 = uCF_P(iPoint,iCF_S3) - uCF_M(1,iCF_S3)
                dE  = uCF_P(iPoint,iCF_E)  - uCF_M(1,iCF_E)

                a = dD * dE - 0.5_DP * ( dS1**2 + dS2**2 + dS3**3 )
                b = dE * uCF_M(1,iCF_D) + dD * uCF_M(1,iCF_E) &
                      - ( dS1 * uCF_M(1,iCF_S1) + dS2 * uCF_M(1,iCF_S2) &
                            + dS3 * uCF_M(1,iCF_S3) )
                c = uCF_M(1,iCF_D) * uCF_M(1,iCF_E) &
                      - 0.5_DP * ( uCF_M(1,iCF_S1)**2 + uCF_M(1,iCF_S2)**2 &
                                     + uCF_M(1,iCF_S3) + Tol_E )

                CALL GetRoots_Quadratic( a, b, c, r1, r2 )

                IF( r1 < 0.0_DP .AND. r2 < 0.0_DP )THEN
                  Theta = 0.0_DP
                ELSE
                  IF( r1 > 0.0_DP ) Theta = MIN( r1, Theta )
                  IF( r2 > 0.0_DP ) Theta = MIN( r2, Theta )
                END IF

              END DO

              IF( Theta < 1.0_DP )THEN

                DO iCF = 1, nCF
                  uCF_M(2:nDOFX,iCF) &
                    = Theta * uCF_M(2:nDOFX,iCF)
                END DO

              END IF

            END IF

            ! --- Back to Nodal Representation ---

            DO iCF = 1, nCF

              CALL MapModalToNodal( uCF(:,iX1,iX2,iX3,iCF), uCF_M(:,iCF) )

            END DO

            uPF_A(1:nPF) = Primitive( uCF_M(1,1:nCF) )

            DO iPoint = 1, nPoints
              DO iCF = 1, nCF

                uCF_P(iPoint,iCF) &
                  = evalLX( uCF(:,iX1,iX2,iX3,iCF), Points_X1(iPoint), &
                            Points_X2(iPoint), Points_X3(iPoint) )

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

  END SUBROUTINE ApplyPositivityLimiter_Euler_DG


END MODULE EulerEquationsLimiterModule_DG
