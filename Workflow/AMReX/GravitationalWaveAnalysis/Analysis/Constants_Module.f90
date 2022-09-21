MODULE Constants_Module

     ! Constants and RadHyd definitions

     INTEGER, PARAMETER            :: single = 4
     INTEGER, PARAMETER            :: double = 8
     REAL(KIND=double), PARAMETER  :: zero   = 0.0d0
     REAL(KIND=double), PARAMETER  :: half   = 0.5d0
     REAL(KIND=double), PARAMETER  :: one    = 1.0d0
     REAL(KIND=double), PARAMETER  :: two    = 2.0d0
     REAL(KIND=double), PARAMETER  :: four   = 4.0d0
     REAL(KIND=double), PARAMETER  :: eight  = 8.0d0
     REAL(KIND=double), PARAMETER  :: pi     = 3.1415926535897932385d0
     REAL(KIND=double), PARAMETER  :: pi2    = two * pi
     REAL(KIND=double), PARAMETER  :: pi4    = four * pi
     REAL(KIND=double), PARAMETER  :: coef = sqrt(15.0d0/pi)/eight
                                        ! coefficient due to spherical harmonics decomposition

     ! Natural constants in cgs

     REAL(KIND=double), PARAMETER  :: cvel = 2.99792458d10
     REAL(KIND=double), PARAMETER  :: Gc   = 6.67408d-8
     REAL(KIND=double), PARAMETER  :: epsilon = 1.0d-100
     REAL(KIND=double), PARAMETER  :: MassSol = 1.989d33
     REAL(KIND=double), PARAMETER  :: Lum_units = (cvel**3/ Gc / 32.d0 / pi)


     ! Numerical coefficients for Quadrupole formula

     INTEGER(KIND=single), PARAMETER :: LL   = 4    ! number of harmonics
     INTEGER(KIND=single), PARAMETER :: mth  = 32   !number of Gauss-Legendre nodal point in Y-direction
     INTEGER(KIND=single), PARAMETER :: kphi = 100  !number of Gauss-Legendre nodal point in Z-direction
     INTEGER(KIND=single), DIMENSION(2), PARAMETER :: iopt = 2

     REAL(KIND=double), PARAMETER    :: Nc1 = Gc * four * pi**half  / (cvel**4*sqrt(1.0d1))
     REAL(KIND=double), PARAMETER    :: Nc2 = Gc * eight * pi**half / (cvel**4*sqrt(1.0d1))
     REAL(KIND=double), PARAMETER    :: Nc3 = Gc * eight * pi**half / (cvel**4*sqrt(1.5d1))
     REAL(KIND=double), DIMENSION(5), PARAMETER   :: Nc = [Nc1, Nc2, Nc3, Nc2, Nc1]
     REAL(KIND=double), DIMENSION(2), PARAMETER   :: df = zero

     COMPLEX(KIND=double), PARAMETER :: ui = (0.0d0,1.0d0) ! Imaginary unit


END MODULE Constants_Module
