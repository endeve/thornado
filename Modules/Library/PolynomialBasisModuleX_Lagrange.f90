MODULE PolynomialBasisModuleX_Lagrange

  USE KindModule, ONLY: &
    DP, Zero, One
  USE QuadratureModule, ONLY: &
    InitializeQuadratures, &
    xG1, xG2, xG3, xG4
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nNodesX, nNodes

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

  REAL(DP),                  ALLOCATABLE         :: nodes(:,:)
  INTEGER,                   ALLOCATABLE, PUBLIC :: IndLX_Q(:,:)
  TYPE(PolynomialBasisType), ALLOCATABLE, PUBLIC :: L_X1(:), dL_X1(:)
  TYPE(PolynomialBasisType), ALLOCATABLE, PUBLIC :: L_X2(:), dL_X2(:)
  TYPE(PolynomialBasisType), ALLOCATABLE, PUBLIC :: L_X3(:), dL_X3(:)

  PUBLIC :: InitializePolynomialBasisX_Lagrange

CONTAINS


  SUBROUTINE InitializePolynomialBasisX_Lagrange

    CALL InitializeQuadratures

    ALLOCATE( nodes(nNodes,nNodes) )
    CALL InitializeNodes

    ALLOCATE( L_X1(nNodesX(1)), dL_X1(nNodesX(1)) )
    CALL InitializeBasis( L_X1, dL_X1 )

    ALLOCATE( L_X2(nNodesX(2)), dL_X2(nNodesX(2)) )
    CALL InitializeBasis( L_X2, dL_X2 )

    ALLOCATE( L_X3(nNodesX(3)), dL_X3(nNodesX(3)) )
    CALL InitializeBasis( L_X3, dL_X3 )

    ALLOCATE( IndLX_Q(1:3,nDOFX) )
    CALL InitializeIndices_TensorProductBasis

  END SUBROUTINE InitializePolynomialBasisX_Lagrange


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

    INTEGER :: i, iNodeX1, iNodeX2, iNodeX3

    i = 0
    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)
    DO iNodeX1 = 1, nNodesX(1)

      i = i + 1

      IndLX_Q(1:3,i) = [ iNodeX1, iNodeX2, iNodeX3 ]

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeIndices_TensorProductBasis


  !***************************************************************************
  !  Elements of Lagrange Polynomial Basis
  !***************************************************************************


  PURE REAL(DP) FUNCTION Unity( x )
    REAL(DP), INTENT(IN) :: x
    Unity = One
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

    Lagrange = One
    DO i = 1, nNodes - 1
      Lagrange = Lagrange * ( x - xx(i) ) / ( xx(nNodes) - xx(i) )
    END DO

  END FUNCTION Lagrange


  !***************************************************************************
  !  Derivative of Lagrange Polynomial Basis
  !***************************************************************************


  PURE REAL(DP) FUNCTION Empty( x )

    REAL(DP), INTENT(IN) :: x

    Empty = Zero

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

    Denominator = One
    DO i = 1, nNodes - 1
      Denominator = Denominator * ( xx(nNodes) - xx(i) )
    END DO

    dLagrange = Zero
    DO i = 1, nNodes - 1
      Numerator = One
      DO j = 1, nNodes - 1
        IF( j == i ) CYCLE
        Numerator = Numerator * ( x - xx(j) )
      END DO
      dLagrange = dLagrange + Numerator / Denominator
    END DO

  END FUNCTION dLagrange


END MODULE PolynomialBasisModuleX_Lagrange
