MODULE PolynomialBasisModuleX_Legendre

  USE KindModule, ONLY: &
    DP
  USE QuadratureModule, ONLY: &
    InitializeQuadratures, &
    xG5, wG5
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nNodesX, nDimsX, nNodes

  IMPLICIT NONE
  PRIVATE

  INTERFACE
    PURE REAL(DP) FUNCTION Basis( X )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION Basis
  END INTERFACE

  TYPE :: PolynomialBasisType
    PROCEDURE (Basis), POINTER, NOPASS :: P
  END TYPE PolynomialBasisType

  INTEGER,                   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: IndPX_Q
  REAL(DP),                  DIMENSION(:),   ALLOCATABLE, PUBLIC :: MassPX
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: P_X1
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: P_X2
  TYPE(PolynomialBasisType), DIMENSION(:),   ALLOCATABLE, PUBLIC :: P_X3

  PUBLIC :: InitializePolynomialBasisX_Legendre
  PUBLIC :: evalPX

CONTAINS


  SUBROUTINE InitializePolynomialBasisX_Legendre

    CALL InitializeQuadratures

    ALLOCATE( P_X1(nNodesX(1)) )
    CALL InitializeBasis( P_X1 )

    ALLOCATE( P_X2(nNodesX(2)) )
    CALL InitializeBasis( P_X2 )

    ALLOCATE( P_X3(nNodesX(3)) )
    CALL InitializeBasis( P_X3 )

    ALLOCATE( IndPX_Q(1:3,nDOFX) )
    CALL InitializeIndices_TensorProductBasis

    ALLOCATE( MassPX(nDOFX) )
    CALL ComputeMassMatrix

  END SUBROUTINE InitializePolynomialBasisX_Legendre


  SUBROUTINE InitializeBasis( P )

    TYPE(PolynomialBasisType), DIMENSION(:) :: P

    SELECT CASE ( SIZE( P ) )

      CASE ( 1 )

        P(1) % P => Legendre_0

      CASE ( 2 )

        P(1) % P => Legendre_0
        P(2) % P => Legendre_1

      CASE ( 3 )

        P(1) % P => Legendre_0
        P(2) % P => Legendre_1
        P(3) % P => Legendre_2

      CASE ( 4 )

        P(1) % P => Legendre_0
        P(2) % P => Legendre_1
        P(3) % P => Legendre_2
        P(4) % P => Legendre_3

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,*) "  Invalid Basis Size  "
        STOP

    END SELECT

  END SUBROUTINE InitializeBasis


  SUBROUTINE InitializeIndices_TensorProductBasis

    INTEGER :: i, ik, kMax, ik_X1, ik_X2, ik_X3

    kMax = nDimsX * nNodes

    i = 0
    DO ik    = 0, kMax
    DO ik_X3 = 0, nNodesX(3) - 1
    DO ik_X2 = 0, nNodesX(2) - 1
    DO ik_X1 = 0, nNodesX(1) - 1

      IF( SUM( [ ik_X1, ik_X2, ik_X3 ] ) == ik )THEN

        i = i + 1

        IndPX_Q(1:3,i) = [ ik_X1+1, ik_X2+1, ik_X3+1 ]

      END IF

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeIndices_TensorProductBasis


  SUBROUTINE ComputeMassMatrix

    INTEGER :: i, j
    INTEGER :: qX1, qX2, qX3

    MassPX = 0.0_DP
    DO j = 1, nDOFX
    DO i = 1, nDOFX

      IF( i /= j ) cycle ! --- Diagonal Mass Matrix

      DO qX3 = 1, SIZE( xG5 )
      DO qX2 = 1, SIZE( xG5 )
      DO qX1 = 1, SIZE( xG5 )

        MassPX(i) &
          = MassPX(i) &
              + wG5(qX1) * wG5(qX2) * wG5(qX3) &
                  * P_X1(IndPX_Q(1,i)) % P( xG5(qX1) ) &
                  * P_X2(IndPX_Q(2,i)) % P( xG5(qX2) ) &
                  * P_X3(IndPX_Q(3,i)) % P( xG5(qX3) ) &
                  * P_X1(IndPX_Q(1,j)) % P( xG5(qX1) ) &
                  * P_X2(IndPX_Q(2,j)) % P( xG5(qX2) ) &
                  * P_X3(IndPX_Q(3,j)) % P( xG5(qX3) )

      END DO
      END DO
      END DO

    END DO
    END DO

    MassPX = 1.0_DP / MassPX ! --- Store as Inverse

  END SUBROUTINE ComputeMassMatrix


  PURE REAL(DP) FUNCTION evalPX( u, X1, X2, X3 )

    REAL(DP), INTENT(in) :: u(1:nDOFX), X1, X2, X3

    INTEGER :: i

    evalPX = 0.0_DP
    DO i = 1, nDOFX

      evalPX &
        = evalPX &
            + P_X1(IndPX_Q(1,i)) % P( X1 ) &
            * P_X2(IndPX_Q(2,i)) % P( X2 ) &
            * P_X3(IndPX_Q(3,i)) % P( X3 ) &
            * u(i)

    END DO

    RETURN
  END FUNCTION evalPX


  !***************************************************************************
  !  Elements of Legendre Polynomial Basis
  !***************************************************************************


  PURE REAL(DP) FUNCTION Legendre_0( x )
    REAL(DP), INTENT(IN) :: x
    Legendre_0 = 1.0_DP
    RETURN
  END FUNCTION Legendre_0


  PURE REAL(DP) FUNCTION Legendre_1( x )
    REAL(DP), INTENT(IN) :: x
    Legendre_1 = x
    RETURN
  END FUNCTION Legendre_1


  PURE REAL(DP) FUNCTION Legendre_2( x )
    REAL(DP), INTENT(IN) :: x
    Legendre_2 = x**2 - 1.0_DP / 12.0_DP
    RETURN
  END FUNCTION Legendre_2


  PURE REAL(DP) FUNCTION Legendre_3( x )
    REAL(DP), INTENT(IN) :: x
    Legendre_3 = x * ( x**2 - 3.0_DP / 20.0_DP )
    RETURN
  END FUNCTION Legendre_3


END MODULE PolynomialBasisModuleX_Legendre
