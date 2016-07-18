MODULE KindModule

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: &
    DP = KIND( 1.d0 )
  REAL(DP), PUBLIC, PARAMETER :: &
    Pi     = ACOS( - 1.0_DP ), &
    TwoPi  = 2.0_DP * Pi, &
    FourPi = 4.0_DP * Pi

END MODULE KindModule
