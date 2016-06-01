MODULE PolynomialBasisModule_Legendre

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nNodesX
  USE QuadratureModule, ONLY: &
    xG5, wG5
  USE UtilitiesModule, ONLY: &
    WriteMatrix

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

  REAL(DP),                  DIMENSION(:), ALLOCATABLE, PUBLIC :: MassPX
  TYPE(PolynomialBasisType), DIMENSION(:), ALLOCATABLE, PUBLIC :: P_X1
  TYPE(PolynomialBasisType), DIMENSION(:), ALLOCATABLE, PUBLIC :: P_X2
  TYPE(PolynomialBasisType), DIMENSION(:), ALLOCATABLE, PUBLIC :: P_X3

  PUBLIC :: InitializePolynomialBasis_Legendre
  PUBLIC :: evalPX

CONTAINS


  SUBROUTINE InitializePolynomialBasis_Legendre

    ALLOCATE( P_X1(nNodesX(1)) )
    CALL InitializeBasis( P_X1 )

    ALLOCATE( P_X2(nNodesX(2)) )
    CALL InitializeBasis( P_X2 )

    ALLOCATE( P_X3(nNodesX(3)) )
    CALL InitializeBasis( P_X3 )

    ALLOCATE( MassPX(nDOFX) )
    CALL ComputeMassMatrices

  END SUBROUTINE InitializePolynomialBasis_Legendre


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


  SUBROUTINE ComputeMassMatrices

    INTEGER :: iX1, iX2, iX3
    INTEGER :: jX1, jX2, jX3
    INTEGER :: qX1, qX2, qX3
    INTEGER :: i, j
    REAL(DP), DIMENSION(:), ALLOCATABLE :: M_X

    ALLOCATE( M_X(1:nDOFX) )

    i = 0
    DO iX1 = 1, nNodesX(1)
      DO iX2 = 1, nNodesX(2)
        DO iX3 = 1, nNodesX(3)

          i = i + 1

          j = 0
          DO jX1 = 1, nNodesX(1)
            DO jX2 = 1, nNodesX(2)
              DO jX3 = 1, nNodesX(3)

                j = j + 1

                IF( j /= i ) cycle

                M_X(i) = 0.0_DP
                DO qX3 = 1, SIZE( xG5 )
                  DO qX2 = 1, SIZE( xG5 )
                    DO qX1 = 1, SIZE( xG5 )

                      M_X(i) &
                        = M_X(i) &
                            + wG5(qX1) * wG5(qX2) * wG5(qX3) &
                                * P_X1(iX1) % P( xG5(qX1) ) &
                                    * P_X2(iX2) % P( xG5(qX2) ) &
                                        * P_X3(iX3) % P( xG5(qX3) ) &
                                * P_X1(jX1) % P( xG5(qX1) ) &
                                    * P_X2(jX2) % P( xG5(qX2) ) &
                                        * P_X3(jX3) % P( xG5(qX3) )

                    END DO
                  END DO
                END DO

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

    MassPX = 1.0_DP / M_X ! Store as Inverse

    DEALLOCATE( M_X )

  END SUBROUTINE ComputeMassMatrices


  PURE REAL(DP) FUNCTION evalPX( u, X1, X2, X3 )

    REAL(DP), INTENT(in) :: u(1:nDOFX), X1, X2, X3

    INTEGER  :: i, iX1, iX2, iX3
    REAL(DP) :: TMP_X2, TMP_X3

    i = 1
    evalPX = 0.0_DP
    DO iX3 = 1, nNodesX(3)
      TMP_X3 = P_X3(iX3) % P( X3 )
      DO iX2 = 1, nNodesX(2)
        TMP_X2 = P_X2(iX2) % P( X2 )
        DO iX1 = 1, nNodesX(1)

          evalPX &
            = evalPX &
              + P_X1(iX1) % P( X1 ) * u(i) &
                  * TMP_X2 * TMP_X3

          i = i + 1

        END DO
      END DO
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


END MODULE PolynomialBasisModule_Legendre
