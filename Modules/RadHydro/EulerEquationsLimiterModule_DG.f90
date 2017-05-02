MODULE EulerEquationsLimiterModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE PolynomialBasisModule_Legendre, ONLY: &
    P_X1, P_X2, P_X3, &
    IndPX_Q
  USE MeshModule, ONLY: &
    MeshX
  USE EulerEquationsLimiterUtilitiesModule_DG, ONLY: &
    InitializeDiscontinuityDetector

  IMPLICIT NONE
  PRIVATE

  LOGICAL,                               PUBLIC :: ApplySlopeLimiter
  REAL(DP),                              PUBLIC :: BetaTVD
  REAL(DP),                              PUBLIC :: BetaTVB
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Legendre

  LOGICAL,                               PUBLIC :: ApplyPositivityLimiter
  INTEGER,                               PUBLIC :: nPositivePoints
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_X1
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_X2
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_X3  
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Lagrange

  PUBLIC :: InitializeLimiters_Euler_DG

CONTAINS


  SUBROUTINE InitializeLimiters_Euler_DG &
               ( ApplySlopeLimiter_Option, BetaTVB_Option, BetaTVD_Option, &
                 ApplyPositivityLimiter_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP), INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,  INTENT(in), OPTIONAL :: ApplyPositivityLimiter_Option

    ! --- Limiter Parameters ---

    ApplySlopeLimiter = .TRUE.
    IF( PRESENT( ApplySlopeLimiter_Option ) ) &
      ApplySlopeLimiter = ApplySlopeLimiter_Option

    BetaTVB = 0.0_DP
    IF( PRESENT( BetaTVB_Option ) ) &
      BetaTVB = BetaTVB_Option

    BetaTVD = 1.8_DP
    IF( PRESENT( BetaTVD_Option ) ) &
      BetaTVD = BetaTVD_Option

    ApplyPositivityLimiter = .TRUE.
    IF( PRESENT( ApplyPositivityLimiter_Option ) ) &
      ApplyPositivityLimiter = ApplyPositivityLimiter_Option

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'InitializeLimiters_Euler_DG'
    WRITE(*,*)
    WRITE(*,'(A7,A20,L1)') &
      '', 'ApplySlopeLimiter = ', ApplySlopeLimiter
    WRITE(*,'(A7,A10,ES8.2E2)') '', 'BetaTVB = ', BetaTVB
    WRITE(*,'(A7,A10,ES8.2E2)') '', 'BetaTVD = ', BetaTVD
    WRITE(*,'(A7,A25,L1)') &
      '', 'ApplyPositivityLimiter = ', ApplyPositivityLimiter
    WRITE(*,*)

    CALL InitializeDiscontinuityDetector

    CALL InitializeSlopeLimiter_Euler_DG

    CALL InitializePositivityLimiter_Euler_DG

  END SUBROUTINE InitializeLimiters_Euler_DG


  SUBROUTINE InitializeSlopeLimiter_Euler_DG

    INTEGER :: iPol
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX

    ! --- Legendre Polynomials in Gaussian Quadrature Points ---

    ALLOCATE( Legendre(nDOFX,nDOFX) )

    DO iPol = 1, nDOFX

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

            Legendre(iNodeX,iPol) &
                = P_X1  (IndPX_Q(1,iPol)) % P( MeshX(1) % Nodes(iNodeX1) ) &
                  * P_X2(IndPX_Q(2,iPol)) % P( MeshX(2) % Nodes(iNodeX2) ) &
                  * P_X3(IndPX_Q(3,iPol)) % P( MeshX(3) % Nodes(iNodeX3) )

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeSlopeLimiter_Euler_DG


  SUBROUTINE InitializePositivityLimiter_Euler_DG

    INTEGER :: iPoint
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP), DIMENSION(:), ALLOCATABLE :: NodesX1

    ! --- Number of Points Where Positivity is Required ---

    ALLOCATE( NodesX1(nNodesX(1)+2) )

    NodesX1 = [ - 0.5_DP, MeshX(1) % Nodes, + 0.5_DP ]

    nPositivePoints = ( nNodesX(1) + 2 ) * nNodesX(2) * nNodesX(3)

    IF( nNodesX(2) > 1 ) &
      nPositivePoints = nPositivePoints + 2 * nNodesX(1) * nNodesX(3)

    IF( nNodesX(3) > 1 ) &
      nPositivePoints = nPositivePoints + 2 * nNodesX(1) * nNodesX(2)

    ! --- Coordinates of Points Where Positivity is Required ---

    ALLOCATE &
      ( Points_X1(nPositivePoints), &
        Points_X2(nPositivePoints), &
        Points_X3(nPositivePoints) )

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

    ! --- Lagrange Polynomials Evaluated in Positive Points ---

    ALLOCATE( Lagrange(nDOFX,nPositivePoints) )

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

          DO iPoint = 1, nPositivePoints

            Lagrange(iNodeX,iPoint) &
              = L_X1(iNodeX1) % P( Points_X1(iPoint) ) &
                  * L_X2(iNodeX2) % P( Points_X2(iPoint) ) &
                      * L_X3(iNodeX3) % P( Points_X3(iPoint) )

          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE InitializePositivityLimiter_Euler_DG


END MODULE EulerEquationsLimiterModule_DG
