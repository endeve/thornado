MODULE MF_KindModule

  USE amrex_fort_module, ONLY: &
    amrex_real

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: DP = KIND( 1.e0_amrex_real )

  REAL(DP), PUBLIC, PARAMETER :: Zero     = 0.0_DP
  REAL(DP), PUBLIC, PARAMETER :: Half     = 0.5_DP
  REAL(DP), PUBLIC, PARAMETER :: One      = 1.0_DP
  REAL(DP), PUBLIC, PARAMETER :: Two      = 2.0_DP
  REAL(DP), PUBLIC, PARAMETER :: Three    = 3.0_DP
  REAL(DP), PUBLIC, PARAMETER :: Four     = 4.0_DP
  REAL(DP), PUBLIC, PARAMETER :: Five     = 5.0_DP
  REAL(DP), PUBLIC, PARAMETER :: Pi       = ACOS( -One )
  REAL(DP), PUBLIC, PARAMETER :: TwoPi    = Two * Pi
  REAL(DP), PUBLIC, PARAMETER :: FourPi   = Four * Pi
  REAL(DP), PUBLIC, PARAMETER :: SqrtTiny = SQRT( TINY( One ) )

END MODULE MF_KindModule
