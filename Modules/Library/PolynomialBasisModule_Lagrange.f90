MODULE PolynomialBasisModule_Lagrange

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOF, &
    nNodesE, nNodesX, nNodes
  USE QuadratureModule, ONLY: &
    xG1, xG2, xG3, xG4
  USE UtilitiesModule, ONLY: &
    NodeNumber, &
    NodeNumberX, &
    WriteMatrix

  IMPLICIT NONE
  PRIVATE

  INTERFACE
    PURE REAL(DP) FUNCTION Basis( x0 )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(IN) :: x0
    END FUNCTION Basis
  END INTERFACE

  TYPE :: PolynomialBasisType
    PROCEDURE (Basis), POINTER, NOPASS :: P
  END TYPE PolynomialBasisType

  REAL(DP),                  DIMENSION(:,:), ALLOCATABLE         :: nodes
  INTEGER,                   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: IndL_Q
  INTEGER,                   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: IndLx_Q
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: L_E,  dL_E
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: L_X1, dL_X1
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: L_X2, dL_X2
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: L_X3, dL_X3

  PUBLIC :: InitializePolynomialBasis_Lagrange
  PUBLIC :: evalL
  PUBLIC :: evalL_X1
  PUBLIC :: evalLX
  PUBLIC :: evalLX_X1

CONTAINS


  SUBROUTINE InitializePolynomialBasis_Lagrange

    ALLOCATE( nodes(nNodes,nNodes) )
    CALL InitializeNodes

    ALLOCATE( L_E(nNodesE), dL_E(nNodesE) )
    CALL InitializeBasis( L_E, dL_E )

    ALLOCATE( L_X1(nNodesX(1)), dL_X1(nNodesX(1)) )
    CALL InitializeBasis( L_X1, dL_X1 )

    ALLOCATE( L_X2(nNodesX(2)), dL_X2(nNodesX(2)) )
    CALL InitializeBasis( L_X2, dL_X2 )

    ALLOCATE( L_X3(nNodesX(3)), dL_X3(nNodesX(3)) )
    CALL InitializeBasis( L_X3, dL_X3 )

    ALLOCATE( IndL_Q (0:3,nDOF) )
    ALLOCATE( IndLX_Q(1:3,nDOF) )
    CALL InitializeIndices_TensorProductBasis

    CALL ComputeMassMatrix

  END SUBROUTINE InitializePolynomialBasis_Lagrange


  SUBROUTINE InitializeNodes

    SELECT CASE ( nNodes )

      CASE ( 1 )

        CALL SetNodes( xG1 )

      CASE ( 2 )

        CALL SetNodes( xG2 )

      CASE ( 3 )

        CALL SetNodes( xG3 )

      CASE ( 4 )

        CALL SetNodes( xG4 )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,*) "  Invalid nNodes  "
        STOP

    END SELECT

  END SUBROUTINE InitializeNodes


  SUBROUTINE SetNodes( x )

    REAL(DP), DIMENSION(nNodes), INTENT(inout) :: x

    INTEGER :: i, j, k

    nodes(nNodes,:) = x

    DO j = 1, nNodes
      k = 1
      DO i = 1, nNodes
        IF( i == j ) CYCLE
        nodes(k,j) = x(i)
        k = k + 1
      END DO
    END DO

  END SUBROUTINE SetNodes


  SUBROUTINE InitializeBasis( P, dP )

    TYPE(PolynomialBasisType), DIMENSION(:) :: P, dP

    SELECT CASE ( SIZE( P ) )

      CASE ( 1 )

        P (1) % P => Unity
        dP(1) % P => Empty

      CASE ( 2 )

        P (1) % P =>  Lagrange_1
        P (2) % P =>  Lagrange_2
        dP(1) % P => dLagrange_1
        dP(2) % P => dLagrange_2

      CASE ( 3 )

        P (1) % P =>  Lagrange_1
        P (2) % P =>  Lagrange_2
        P (3) % P =>  Lagrange_3
        dP(1) % P => dLagrange_1
        dP(2) % P => dLagrange_2
        dP(3) % P => dLagrange_3

      CASE ( 4 )

        P (1) % P =>  Lagrange_1
        P (2) % P =>  Lagrange_2
        P (3) % P =>  Lagrange_3
        P (4) % P =>  Lagrange_4
        dP(1) % P => dLagrange_1
        dP(2) % P => dLagrange_2
        dP(3) % P => dLagrange_3
        dP(4) % P => dLagrange_4

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,*) "  Invalid Basis Size  "
        STOP

    END SELECT

  END SUBROUTINE InitializeBasis


  SUBROUTINE InitializeIndices_TensorProductBasis

    INTEGER :: i, iNodeE, iNodeX1, iNodeX2, iNodeX3

    i = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)

          IndLX_Q(1:3,i) = [ iNodeX1, iNodeX2, iNodeX3 ]
          i = i + 1

        END DO
      END DO
    END DO

    i = 1
    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            IndL_Q(0:3,i) = [ iNodeE, iNodeX1, iNodeX2, iNodeX3 ]
            i = i + 1

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeIndices_TensorProductBasis


  SUBROUTINE ComputeMassMatrix

    REAL(DP), DIMENSION(1:nNodesE,   1:nNodesE)    :: M_E
    REAL(DP), DIMENSION(1:nNodesX(1),1:nNodesX(1)) :: M_X1
    REAL(DP), DIMENSION(1:nNodesX(2),1:nNodesX(2)) :: M_X2
    REAL(DP), DIMENSION(1:nNodesX(3),1:nNodesX(3)) :: M_X3

  END SUBROUTINE ComputeMassMatrix


  PURE REAL(DP) FUNCTION evalL( u, E, X1, X2, X3 )

    REAL(DP), INTENT(in) :: u(1:nDOF), E, X1, X2, X3

    INTEGER  :: i, iE, iX1, iX2, iX3
    REAL(DP) :: TMP_X1, TMP_X2, TMP_X3

    i = 1
    evalL = 0.0_DP
    DO iX3 = 1, nNodesX(3)
      TMP_X3 = L_X3(iX3) % P( X3 )
      DO iX2 = 1, nNodesX(2)
        TMP_X2 = L_X2(iX2) % P( X2 )
        DO iX1 = 1, nNodesX(1)
          TMP_X1 = L_X1(iX1) % P( X1 )
          DO iE = 1, nNodesE

            evalL &
              = evalL &
                + L_E(iE) % P( E ) * u(i) &
                    * TMP_X1 * TMP_X2 * TMP_X3

            i = i + 1

          END DO
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION evalL


  PURE REAL(DP) FUNCTION evalL_X1( u, X1, iE, iX2, iX3 )

    REAL(DP), INTENT(in) :: u(1:nDOF), X1
    INTEGER,  INTENT(in) :: iE, iX2, iX3

    INTEGER :: iX1, iNode

    evalL_X1 = 0.0_DP
    DO iX1 = 1, nNodesX(1)

      iNode = NodeNumber( iE, iX1, iX2, iX3 )

      evalL_X1 &
        = evalL_X1 &
            + L_X1(iX1) % P( X1 ) * u(iNode)

    END DO

    RETURN
  END FUNCTION evalL_X1


  PURE REAL(DP) FUNCTION evalLX( u, X1, X2, X3 )

    REAL(DP), INTENT(in) :: u(1:nDOFX), X1, X2, X3

    INTEGER  :: i, iX1, iX2, iX3
    REAL(DP) :: TMP_X2, TMP_X3

    i = 1
    evalLX = 0.0_DP
    DO iX3 = 1, nNodesX(3)
      TMP_X3 = L_X3(iX3) % P( X3 )
      DO iX2 = 1, nNodesX(2)
        TMP_X2 = L_X2(iX2) % P( X2 )
        DO iX1 = 1, nNodesX(1)

          evalLX &
            = evalLX &
              + L_X1(iX1) % P( X1 ) * u(i) &
                  * TMP_X2 * TMP_X3

          i = i + 1

        END DO
      END DO
    END DO

    RETURN
  END FUNCTION evalLX


  PURE REAL(DP) FUNCTION evalLX_X1( u, X1, iX2, iX3 )

    REAL(DP), INTENT(in) :: u(1:nDOFX), X1
    INTEGER,  INTENT(in) :: iX2, iX3

    INTEGER :: iX1, iNode

    evalLX_X1 = 0.0_DP
    DO iX1 = 1, nNodesX(1)

      iNode = NodeNumberX( iX1, iX2, iX3 )

      evalLX_X1 &
        = evalLX_X1 &
            + L_X1(iX1) % P( X1 ) * u(iNode)

    END DO

    RETURN
  END FUNCTION evalLX_X1


  !***************************************************************************
  !  Elements of Lagrange Polynomial Basis
  !***************************************************************************


  PURE REAL(DP) FUNCTION Unity( x )
    REAL(DP), INTENT(IN) :: x
    Unity = 1.0_DP
    RETURN
  END FUNCTION Unity


  PURE REAL(DP) FUNCTION Lagrange_1( x )
    REAL(DP), INTENT(IN) :: x
    Lagrange_1 = Lagrange( x, nodes(:,1) )
    RETURN
  END FUNCTION Lagrange_1


  PURE REAL(DP) FUNCTION Lagrange_2( x )
    REAL(DP), INTENT(IN) :: x
    Lagrange_2 = Lagrange( x, nodes(:,2) )
    RETURN
  END FUNCTION Lagrange_2


  PURE REAL(DP) FUNCTION Lagrange_3( x )
    REAL(DP), INTENT(IN) :: x
    Lagrange_3 = Lagrange( x, nodes(:,3) )
    RETURN
  END FUNCTION Lagrange_3


  PURE REAL(DP) FUNCTION Lagrange_4( x )
    REAL(DP), INTENT(IN) :: x
    Lagrange_4 = Lagrange( x, nodes(:,4) )
    RETURN
  END FUNCTION Lagrange_4


  PURE REAL(DP) FUNCTION Lagrange( x, xx )

    REAL(DP),                    INTENT(in) :: x
    REAL(DP), DIMENSION(nNodes), INTENT(in) :: xx

    INTEGER :: i

    Lagrange = 1.0_DP
    DO i = 1, nNodes - 1
      Lagrange = Lagrange * ( x - xx(i) ) / ( xx(nNodes) - xx(i) )
    END DO

  END FUNCTION Lagrange


  !***************************************************************************
  !  Derivative of Lagrange Polynomial Basis
  !***************************************************************************


  PURE REAL(DP) FUNCTION Empty( x )

    REAL(DP), INTENT(IN) :: x

    Empty = 0.0_DP

  END FUNCTION Empty


  PURE REAL(DP) FUNCTION dLagrange_1( x )

    REAL(DP), INTENT(IN) :: x

    dLagrange_1 = dLagrange( x, nodes(:,1) )

  END FUNCTION dLagrange_1


  PURE REAL(DP) FUNCTION dLagrange_2( x )

    REAL(DP), INTENT(IN) :: x

    dLagrange_2 = dLagrange( x, nodes(:,2) )

  END FUNCTION dLagrange_2


  PURE REAL(DP) FUNCTION dLagrange_3( x )

    REAL(DP), INTENT(IN) :: x

    dLagrange_3 = dLagrange( x, nodes(:,3) )

  END FUNCTION dLagrange_3


  PURE REAL(DP) FUNCTION dLagrange_4( x )

    REAL(DP), INTENT(IN) :: x

    dLagrange_4 = dLagrange( x, nodes(:,4) )

  END FUNCTION dLagrange_4


  PURE REAL(DP) FUNCTION dLagrange( x, xx )

    REAL(DP),                    INTENT(IN) :: x
    REAL(DP), DIMENSION(nNodes), INTENT(IN) :: xx

    INTEGER  :: i, j
    REAL(DP) :: Denominator, Numerator

    Denominator = 1.0_DP
    DO i = 1, nNodes - 1
      Denominator = Denominator * ( xx(nNodes) - xx(i) )
    END DO

    dLagrange = 0.0_DP
    DO i = 1, nNodes - 1
      Numerator = 1.0_DP
      DO j = 1, nNodes - 1
        IF( j == i ) CYCLE
        Numerator = Numerator * ( x - xx(j) )
      END DO
      dLagrange = dLagrange + Numerator / Denominator
    END DO

  END FUNCTION dLagrange


END MODULE PolynomialBasisModule_Lagrange
