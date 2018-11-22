MODULE LagrangePolynomialsModule

  USE KindModule, ONLY: &
    DP, Zero, One

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LagrangeP
  PUBLIC :: dLagrangeP

CONTAINS


  PURE REAL(DP) FUNCTION LagrangeP( x, i, xx, nn )

    INTEGER,  INTENT(in) :: i, nn
    REAL(DP), INTENT(in) :: x, xx(nn)

    INTEGER :: j

    LagrangeP = One
    DO j = 1, nn
      IF( j == i ) CYCLE
      LagrangeP = LagrangeP * ( x - xx(j) ) / ( xx(i) - xx(j) )
    END DO

    RETURN
  END FUNCTION LagrangeP


  PURE REAL(DP) FUNCTION dLagrangeP( x, i, xx, nn )

    INTEGER,  INTENT(in) :: i, nn
    REAL(DP), INTENT(in) :: x, xx(nn)

    INTEGER  :: j, k
    REAL(DP) :: Denominator, Numerator

    Denominator = One
    DO j = 1, nn
      IF( j == i ) CYCLE
      Denominator = Denominator * ( xx(i) - xx(j) )
    END DO

    dLagrangeP = Zero
    DO j = 1, nn
      IF( j == i ) CYCLE
      Numerator = One
      DO k = 1, nn
        IF( k == i .OR. k == j ) CYCLE
        Numerator = Numerator * ( x - xx(k) )
      END DO
      dLagrangeP = dLagrangeP + Numerator / Denominator
    END DO

    RETURN
  END FUNCTION dLagrangeP


END MODULE LagrangePolynomialsModule
