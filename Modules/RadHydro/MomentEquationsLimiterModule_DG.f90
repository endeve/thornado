MODULE MomentEquationsLimiterModule_DG

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
  USE PolynomialBasisModule_Legendre, ONLY: &
    P_X1, P_X2, P_X3, IndPX_Q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE MomentEquationsLimiterUtilitiesModule_DG, ONLY: &
    InitializeDiscontinuityDetector

  IMPLICIT NONE
  PRIVATE

  LOGICAL,                               PUBLIC :: ApplySlopeLimiter
  REAL(DP),                              PUBLIC :: BetaTVD
  REAL(DP),                              PUBLIC :: BetaTVB
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: LegendreX

  LOGICAL,                               PUBLIC :: ApplyPositivityLimiter
  INTEGER,                               PUBLIC :: nPositivePoints
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_E
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_X1
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_X2
  REAL(DP), DIMENSION(:),   ALLOCATABLE, PUBLIC :: Points_X3
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Lagrange

  PUBLIC :: InitializeLimiters_M1_DG

CONTAINS


  SUBROUTINE InitializeLimiters_M1_DG &
               ( ApplySlopeLimiter_Option, BetaTVB_Option, BetaTVD_Option, &
                 ApplyPositivityLimiter_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP), INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,  INTENT(in), OPTIONAL :: ApplyPositivityLimiter_Option

    
    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode, iNodeX
    INTEGER :: jNodeE, jNodeX1, jNodeX2, jNodeX3, jNode
    REAL(DP) :: eta_E, eta_X1, eta_X2, eta_X3
    

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
    WRITE(*,'(A5,A)') '', 'InitializeLimiters_M1_DG'
    WRITE(*,*)
    WRITE(*,'(A7,A20,L1)') &
      '', 'ApplySlopeLimiter = ', ApplySlopeLimiter
    WRITE(*,'(A7,A10,ES8.2E2)') '', 'BetaTVB = ', BetaTVB
    WRITE(*,'(A7,A10,ES8.2E2)') '', 'BetaTVD = ', BetaTVD
    WRITE(*,'(A7,A25,L1)') &
      '', 'ApplyPositivityLimiter = ', ApplyPositivityLimiter
    WRITE(*,*)

    CALL InitializeDiscontinuityDetector

    CALL InitializeSlopeLimiter_M1_DG

    CALL InitializePositivityLimiter_M1_DG

  END SUBROUTINE InitializeLimiters_M1_DG


  SUBROUTINE InitializeSlopeLimiter_M1_DG

    INTEGER :: iPol
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNode, iNodeX

    ! --- Legendre Polynomials in Gaussian Spatial Quadrature Points ---

    ALLOCATE( LegendreX(nDOFX,nDOFX) )

    DO iPol = 1, nDOFX

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)

            iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

            LegendreX(iNodeX,iPol) &
                = P_X1  (IndPX_Q(1,iPol)) % P( MeshX(1) % Nodes(iNodeX1) ) &
                  * P_X2(IndPX_Q(2,iPol)) % P( MeshX(2) % Nodes(iNodeX2) ) &
                  * P_X3(IndPX_Q(3,iPol)) % P( MeshX(3) % Nodes(iNodeX3) )

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeSlopeLimiter_M1_DG


  SUBROUTINE InitializePositivityLimiter_M1_DG

    INTEGER :: iPoint
    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP), DIMENSION(:), ALLOCATABLE :: NodesX1

    ! --- Number of Points Where Positivity is Required ---

    ALLOCATE( NodesX1(nNodesX(1)+2) )

    NodesX1 = [ - 0.5_DP, MeshX(1) % Nodes, + 0.5_DP ]

    nPositivePoints &
      = nNodesE * ( nNodesX(1) + 2 ) * nNodesX(2) * nNodesX(3)

    IF( nNodesX(2) > 1 ) &
      nPositivePoints &
        = nPositivePoints + 2 * nNodesE * nNodesX(1) * nNodesX(3)

    IF( nNodesX(3) > 1 ) &
      nPositivePoints &
        = nPositivePoints + 2 * nNodesE * nNodesX(2) * nNodesX(3)

    ! --- Coordinates of Points Where Positivity is Required:

    ALLOCATE &
      ( Points_E (nPositivePoints), Points_X1(nPositivePoints), &
        Points_X2(nPositivePoints), Points_X3(nPositivePoints) )

    iPoint = 0

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, SIZE( NodesX1 )
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = NodesX1(iNodeX1)
            Points_X2(iPoint) = MeshX(2) % Nodes(iNodeX2)
            Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

          END DO
        END DO
      END DO
    END DO

    IF( nNodesX(2) > 1 )THEN

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = - 0.5_DP
            Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

          END DO
        END DO
      END DO

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = + 0.5_DP
            Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

          END DO
        END DO
      END DO

    END IF

    IF( nNodesX(3) > 1 )THEN

      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = MeshX(2) % Nodes(iNodeX2)
            Points_X3(iPoint) = - 0.5_DP

          END DO
        END DO
      END DO

      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X3(iPoint) = + 0.5_DP

          END DO
        END DO
      END DO

    END IF

    DEALLOCATE( NodesX1 )

    ! --- Lagrange Polynomials Evaluated in Positive Points ---

    ALLOCATE( Lagrange(nDOF,nPositivePoints) )

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

            DO iPoint = 1, nPositivePoints

              Lagrange(iNode,iPoint) &
                = L_E(iNodeE) % P( Points_E(iPoint) ) &
                    * L_X1(iNodeX1) % P( Points_X1(iPoint) ) &
                        * L_X2(iNodeX2) % P( Points_X2(iPoint) ) &
                            * L_X3(iNodeX3) % P( Points_X3(iPoint) )

            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializePositivityLimiter_M1_DG


END MODULE MomentEquationsLimiterModule_DG
