PROGRAM PolynomialBasisMappingTest

  USE KindModule, ONLY: &
    DP, TwoPi
  USE UtilitiesModule, ONLY: &
    NodeNumber, &
    WriteVector, &
    WriteMatrix
  USE ProgramHeaderModule, ONLY: &
    nNodes, nNodesE, nNodesX, nDOF, nDOFX, &
    nDims, nDimsX
  USE QuadratureModule, ONLY: &
    InitializeQuadratures, &
    GetQuadrature
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange, &
    EvalL
  USE PolynomialBasisModule_Legendre, ONLY: &
    InitializePolynomialBasis_Legendre, &
    EvalP
  USE PolynomialBasisMappingModule, ONLY: &
    InitializePolynomialBasisMapping, &
    MapNodalToModal_Radiation, &
    MapModalToNodal_Radiation

  IMPLICIT NONE

  INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode
  INTEGER :: iPointE, iPointX1
  INTEGER, PARAMETER :: nPointsE = 64, nPointsX1 = 64
  REAL(DP), PARAMETER :: eL = - 0.5_DP, eR = + 0.5_DP
  REAL(DP), PARAMETER :: xL = - 0.5_DP, xR = + 0.5_DP
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: U, U_2, C, C_2
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesE,  WeightsE
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesX1, WeightsX1
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesX2, WeightsX2
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: NodesX3, WeightsX3
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: PointsE
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: PointsX1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: U_2D_L, U_2_2D_L
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: U_2D_P, U_2_2D_P

  nNodes     = 4
  nNodesE    = nNodes
  nNodesX(1) = nNodes
  nNodesX(2) = 1
  nNodesX(3) = 1
  nDOFX      = PRODUCT( nNodesX )
  nDOF       = nDOFX * nNodesE
  nDimsX     = 1
  nDims      = 2

  CALL InitializeQuadratures

  CALL InitializePolynomialBasis_Lagrange

  CALL InitializePolynomialBasis_Legendre

  ALLOCATE( U(nDOF), U_2(nDOF), C(nDOF), C_2(nDOF) )
  ALLOCATE( NodesE(nNodesE),     WeightsE(nNodesE) )
  ALLOCATE( NodesX1(nNodesX(1)), WeightsX1(nNodesX(1)) )
  ALLOCATE( NodesX2(nNodesX(2)), WeightsX2(nNodesX(2)) )
  ALLOCATE( NodesX3(nNodesX(3)), WeightsX3(nNodesX(3)) )
  ALLOCATE( U_2D_L(nPointsE,nPointsX1), U_2_2D_L(nPointsE,nPointsX1) )
  ALLOCATE( U_2D_P(nPointsE,nPointsX1), U_2_2D_P(nPointsE,nPointsX1) )

  CALL GetQuadrature( nNodesE,    NodesE,  WeightsE )
  CALL GetQuadrature( nNodesX(1), NodesX1, WeightsX1 )
  CALL GetQuadrature( nNodesX(2), NodesX2, WeightsX2 )
  CALL GetQuadrature( nNodesX(3), NodesX3, WeightsX3 )

  CALL InitializePolynomialBasisMapping &
         ( NodesE, NodesX1, NodesX2, NodesX3 )

  DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)
        DO iNodeE = 1, nNodesE

          iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

          U(iNode) = COS( TwoPi * ( eL - NodesE(iNodeE) ) ) &
                     * SIN( TwoPi * ( xL - NodesX1(iNodeX1) ) )

        END DO
      END DO
    END DO
  END DO

  CALL MapNodalToModal_Radiation( U, C )

  ! --- Damp Higher Modal Coefficients ---

  C_2(1)      = C(1)
  C_2(2:nDOF) = 0.9_DP * C(2:nDOF)

  CALL MapModalToNodal_Radiation( U_2, C_2 )

  ALLOCATE( PointsE(nPointsE) )
  DO iPointE = 1, nPointsE
    PointsE(iPointE) &
      = eL + REAL(iPointE-1)/REAL(nPointsE-1)
  END DO

  ALLOCATE( PointsX1(nPointsX1) )
  DO iPointX1 = 1, nPointsX1
    PointsX1(iPointX1) &
      = xL + REAL(iPointX1-1)/REAL(nPointsX1-1)
  END DO

  DO iPointX1 = 1, nPointsX1
    DO iPointE = 1, nPointsE

      U_2D_L(iPointE,iPointX1) &
        = EvalL( U, PointsE(iPointE), PointsX1(iPointX1), 0.0_DP, 0.0_DP )

      U_2D_P(iPointE,iPointX1) &
        = EvalP( C, PointsE(iPointE), PointsX1(iPointX1), 0.0_DP, 0.0_DP )

      U_2_2D_L(iPointE,iPointX1) &
        = EvalL( U_2, PointsE(iPointE), PointsX1(iPointX1), 0.0_DP, 0.0_DP )

      U_2_2D_P(iPointE,iPointX1) &
        = EvalP( C_2, PointsE(iPointE), PointsX1(iPointX1), 0.0_DP, 0.0_DP )

    END DO
  END DO

  CALL WriteVector( nPointsE,  PointsE,  'E.dat' )
  CALL WriteVector( nPointsX1, PointsX1, 'X.dat' )
  CALL WriteVector( nDOF, U, 'U.dat' )
  CALL WriteVector( nDOF, C, 'C.dat' )
  CALL WriteMatrix( nPointsE, nPointsX1, U_2D_L,   'U_2D_L.dat' )
  CALL WriteMatrix( nPointsE, nPointsX1, U_2D_P,   'U_2D_P.dat' )
  CALL WriteMatrix( nPointsE, nPointsX1, U_2_2D_L, 'U_2_2D_L.dat' )
  CALL WriteMatrix( nPointsE, nPointsX1, U_2_2D_P, 'U_2_2D_P.dat' )

  DEALLOCATE( U, U_2, C, C_2 )
  DEALLOCATE( NodesE,  WeightsE )
  DEALLOCATE( NodesX1, WeightsX1 )
  DEALLOCATE( NodesX2, WeightsX2 )
  DEALLOCATE( NodesX3, WeightsX3 )
  DEALLOCATE( U_2D_L, U_2_2D_L )
  DEALLOCATE( U_2D_P, U_2_2D_P )
  DEALLOCATE( PointsE )
  DEALLOCATE( PointsX1 )

END PROGRAM PolynomialBasisMappingTest
