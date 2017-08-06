MODULE KindModule

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: &
    DP = KIND( 1.d0 )
  REAL(DP), PUBLIC, PARAMETER :: &
    Zero  = 0.0_DP, &
    One   = 1.0_DP, &
    Two   = 2.0_DP, &
    Three = 3.0_DP, &
    Four  = 4.0_DP, &
    Five  = 5.0_DP
  REAL(DP), PUBLIC, PARAMETER :: &
    Fifth = 1.0_DP / 5.0_DP, &
    Third = 1.0_DP / 3.0_DP, &
    Half  = 0.5_DP
  REAL(DP), PUBLIC, PARAMETER :: &
    Pi     = ACOS( - 1.0_DP ), &
    TwoPi  = 2.0_DP * Pi, &
    FourPi = 4.0_DP * Pi
  REAL(DP), PUBLIC, PARAMETER :: &
    SqrtTiny = SQRT( TINY( 1.0_DP ) )

END MODULE KindModule
