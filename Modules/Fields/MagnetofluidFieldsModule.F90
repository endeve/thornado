MODULE MagnetofluidFieldsModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  ! --- Conserved Magnetofluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iCM_D   = 1  ! Conserved Baryon Density
  INTEGER, PUBLIC, PARAMETER :: iCM_S1  = 2  ! Conserved Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iCM_S2  = 3  ! Conserved Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iCM_S3  = 4  ! Conserved Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: iCM_E   = 5  ! Conserved Energy Density
  INTEGER, PUBLIC, PARAMETER :: iCM_Ne  = 6  ! Conserved Electron Density
  INTEGER, PUBLIC, PARAMETER :: iCM_B1  = 7  ! Conserved Magnetic Field 1
  INTEGER, PUBLIC, PARAMETER :: iCM_B2  = 8  ! Conserved Magnetic Field 2
  INTEGER, PUBLIC, PARAMETER :: iCM_B3  = 9  ! Conserved Magnetic Field 3
  INTEGER, PUBLIC, PARAMETER :: iCM_Chi = 10 ! Divergence Violation Field
  INTEGER, PUBLIC, PARAMETER :: nCM     = 10 ! n Conserved Magnetofluid Fields

  CHARACTER(32), DIMENSION(nCM), PUBLIC, PARAMETER :: &
    namesCM = [ 'Conserved Baryon Density        ', &
                'Conserved Momentum Density (1)  ', &
                'Conserved Momentum Density (2)  ', &
                'Conserved Momentum Density (3)  ', &
                'Conserved Energy Density        ', &
                'Conserved Electron Density      ', &
                'Conserved Magnetic Field (1)    ', &
                'Conserved Magnetic Field (2)    ', &
                'Conserved Magnetic Field (3)    ', &
                'Divergence Violation Field      ' ]

  CHARACTER(10),  DIMENSION(nCM), PUBLIC, PARAMETER :: &
    ShortNamesCM = [ 'CM_D      ', &
                     'CM_S1     ', &
                     'CM_S2     ', &
                     'CM_S3     ', &
                     'CM_E      ', &
                     'CM_Ne     ', &
                     'CM_B1     ', &
                     'CM_B2     ', &
                     'CM_B3     ', &
                     'CM_Chi    ' ]

  REAL(DP), DIMENSION(nCM), PUBLIC :: unitsCM

  REAL(DP), ALLOCATABLE, PUBLIC :: uCM  (:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rhsCM(:,:,:,:,:)

  ! --- Primitive Magnetofluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iPM_D   = 1  ! Comoving Baryon Density
  INTEGER, PUBLIC, PARAMETER :: iPM_V1  = 2  ! Three-Velocity 1
  INTEGER, PUBLIC, PARAMETER :: iPM_V2  = 3  ! Three-Velocity 2
  INTEGER, PUBLIC, PARAMETER :: iPM_V3  = 4  ! Three-Velocity 3
  INTEGER, PUBLIC, PARAMETER :: iPM_E   = 5  ! Internal Energy Density
  INTEGER, PUBLIC, PARAMETER :: iPM_Ne  = 6  ! Comoving Electron Density
  INTEGER, PUBLIC, PARAMETER :: iPM_B1  = 7  ! Comoving Magnetic Field 1
  INTEGER, PUBLIC, PARAMETER :: iPM_B2  = 8  ! Comoving Magnetic Field 2
  INTEGER, PUBLIC, PARAMETER :: iPM_B3  = 9  ! Comoving Magnetic Field 3
  INTEGER, PUBLIC, PARAMETER :: iPM_Chi = 10 ! Divergence Violation Field
  INTEGER, PUBLIC, PARAMETER :: nPM     = 10 ! n Primitive Magnetofluid Fields

  CHARACTER(32), DIMENSION(nPM), PUBLIC, PARAMETER :: &
    namesPM = [ 'Comoving Baryon Density         ', &
                'Three-Velocity (1)              ', &
                'Three-Velocity (2)              ', &
                'Three-Velocity (3)              ', &
                'Internal Energy Density         ', &
                'Comoving Electron Density       ', &
                'Comoving Magnetic Field (1)     ', &
                'Comoving Magnetic Field (2)     ', &
                'Comoving Magnetic Field (3)     ', &
                'Divergence Violation Field      ' ]

  CHARACTER(10),  DIMENSION(nPM), PUBLIC, PARAMETER :: &
    ShortNamesPM = [ 'PM_D      ', &
                     'PM_V1     ', &
                     'PM_V2     ', &
                     'PM_V3     ', &
                     'PM_E      ', &
                     'PM_Ne     ', &  
                     'PM_B1     ', &
                     'PM_B2     ', &
                     'PM_B3     ', &
                     'PM_Chi    ' ]

  REAL(DP), DIMENSION(nPM), PUBLIC :: unitsPM

  REAL(DP), ALLOCATABLE, PUBLIC :: uPM(:,:,:,:,:)

  ! --- Auxiliary Magnetofluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iAM_P       = 01 ! Pressure
  INTEGER, PUBLIC, PARAMETER :: iAM_Pb      = 02 ! Magnetic Pressure
  INTEGER, PUBLIC, PARAMETER :: iAM_T       = 03 ! Temperature
  INTEGER, PUBLIC, PARAMETER :: iAM_Ye      = 04 ! Electron Fraction
  INTEGER, PUBLIC, PARAMETER :: iAM_S       = 05 ! Entropy Per Baryon
  INTEGER, PUBLIC, PARAMETER :: iAM_E       = 06 ! Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: iAM_h       = 07 ! Specific Enthalpy
  INTEGER, PUBLIC, PARAMETER :: iAM_hb      = 08 ! Specific Magnetic Enthalpy
  INTEGER, PUBLIC, PARAMETER :: iAM_Me      = 09 ! Electron Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAM_Mp      = 10 ! Proton Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAM_Mn      = 11 ! Neutron Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAM_Xp      = 12 ! Proton Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAM_Xn      = 13 ! Neutron Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAM_Xa      = 14 ! Alpha Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAM_Xh      = 15 ! Heavy Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAM_Gm      = 16 ! Ratio of Specific Heats
  INTEGER, PUBLIC, PARAMETER :: iAM_Cs      = 17 ! Sound Speed
  INTEGER, PUBLIC, PARAMETER :: iAM_Ca      = 18 ! Alfven Speed
  INTEGER, PUBLIC, PARAMETER :: iAM_Cf_ub_p = 19 ! Fast Upper Bound (Plus)
  INTEGER, PUBLIC, PARAMETER :: iAM_Cf_ub_m = 20 ! Fast Upper Bound (Minus)
  INTEGER, PUBLIC, PARAMETER :: iAM_EF1     = 21 ! Eulerian Electric Field (1)
  INTEGER, PUBLIC, PARAMETER :: iAM_EF2     = 22 ! Eulerian Electric Field (2)
  INTEGER, PUBLIC, PARAMETER :: iAM_EF3     = 23 ! Eulerian Electric Field (3)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem00   = 24 ! EM Stress Tensor (00)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem11   = 25 ! EM Stress Tensor (11)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem22   = 26 ! EM Stress Tensor (22)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem33   = 27 ! EM Stress Tensor (33)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem12   = 28 ! EM Stress Tensor (12)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem13   = 29 ! EM Stress Tensor (13)
  INTEGER, PUBLIC, PARAMETER :: iAM_Tem23   = 30 ! EM Stress Tensor (23)
  INTEGER, PUBLIC, PARAMETER :: nAM         = 30 ! n AuxiliaryMagneto Fields

  CHARACTER(32), DIMENSION(nAM), PUBLIC, PARAMETER :: &
    namesAM = [ 'Pressure                        ', &
                'Magnetic Pressure               ', &
                'Temperature                     ', &
                'Electron Fraction               ', &
                'Entropy Per Baryon              ', &
                'Specific Internal Energy        ', &
                'Specific Enthalpy               ', &
                'Specific Magnetic Enthalpy      ', &
                'Electron Chemical Potential     ', &
                'Proton Chemical Potential       ', &
                'Neutron Chemical Potential      ', &
                'Proton Mass Fraction            ', &
                'Neutron Mass Fraction           ', &
                'Alpha Mass Fraction             ', &
                'Heavy Mass Fraction             ', &
                'Ratio of Specific Heats (Gamma) ', &
                'Sound Speed                     ', &
                'Alfven Speed                    ', &
                'Fast Upper Bound (Plus)         ', &
                'Fast Upper Bound (Minus)        ', &
                'Eulerian Electric Field (1)     ', &
                'Eulerian Electric Field (2)     ', &
                'Eulerian Electric Field (3)     ', &
                'EM Stress Tensor (00)           ', &
                'EM Stress Tensor (11)           ', &
                'EM Stress Tensor (22)           ', &
                'EM Stress Tensor (33)           ', &
                'EM Stress Tensor (12)           ', &
                'EM Stress Tensor (13)           ', &
                'EM Stress Tensor (23)           ' ]

  CHARACTER(10),  DIMENSION(nAM), PUBLIC, PARAMETER :: &
    ShortNamesAM = [ 'AM_P      ', &
                     'AM_Pb     ', &
                     'AM_T      ', &
                     'AM_Ye     ', &
                     'AM_S      ', &
                     'AM_E      ', &
                     'AM_h      ', &
                     'AM_hb     ', &
                     'AM_Me     ', &
                     'AM_Mp     ', &
                     'AM_Mn     ', &
                     'AM_Xp     ', &
                     'AM_Xn     ', &
                     'AM_Xa     ', &
                     'AM_Xh     ', &
                     'AM_Gm     ', &
                     'AM_Cs     ', &
                     'AM_Ca     ', &
                     'AM_Cf_ub_p', &
                     'AM_Cf_ub_m', &
                     'AM_EF1    ', &
                     'AM_EF2    ', &
                     'AM_EF3    ', &
                     'AM_Tem00  ', &
                     'AM_Tem11  ', &
                     'AM_Tem22  ', &
                     'AM_Tem33  ', &
                     'AM_Tem12  ', &
                     'AM_Tem13  ', &
                     'AM_Tem23  ' ]

  REAL(DP), DIMENSION(nAM), PUBLIC :: unitsAM

  REAL(DP), ALLOCATABLE, PUBLIC :: uAM(:,:,:,:,:)

  ! --- Diagnostic Variables ---

  INTEGER, PUBLIC, PARAMETER :: iDM_TCI   = 01 ! Troubled-Cell Indicator
  INTEGER, PUBLIC, PARAMETER :: iDM_Sh_X1 = 02 ! Shock Detector (X1)
  INTEGER, PUBLIC, PARAMETER :: iDM_Sh_X2 = 03 ! Shock Detector (X2)
  INTEGER, PUBLIC, PARAMETER :: iDM_Sh_X3 = 04 ! Shock Detector (X3)
  INTEGER, PUBLIC, PARAMETER :: iDM_T1    = 05 ! Theta 1
  INTEGER, PUBLIC, PARAMETER :: iDM_T2    = 06 ! Theta 2
  INTEGER, PUBLIC, PARAMETER :: iDM_T3    = 07 ! Theta 3
  INTEGER, PUBLIC, PARAMETER :: iDM_MinE  = 08 ! Minimum Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: iDM_MaxE  = 09 ! Maximum Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: iDM_Div   = 10 ! Divergence of Eulerian Magnetic Field
  INTEGER, PUBLIC, PARAMETER :: nDM       = 10 ! n Diagnostic Magnetofluid Fields

  CHARACTER(32), DIMENSION(nDM), PUBLIC, PARAMETER :: &
    namesDM = [ 'TCI                             ', &
                'Shock (X1)                      ', &
                'Shock (X2)                      ', &
                'Shock (X3)                      ', &
                'Theta 1                         ', &
                'Theta 2                         ', &
                'Theta 3                         ', &
                'Min E                           ', &
                'Max E                           ', &
                'Div                             ' ]

  CHARACTER(10), DIMENSION(nDM), PUBLIC, PARAMETER :: &
    ShortNamesDM = [ 'DM_TCI    ', &
                     'DM_Sh_X1  ', &
                     'DM_Sh_X2  ', &
                     'DM_Sh_X3  ', &
                     'DM_T1     ', &
                     'DM_T2     ', &
                     'DM_T3     ', &
                     'DM_MinE   ', &
                     'DM_MaxE   ', &
                     'DM_Div    ' ]

  REAL(DP), DIMENSION(nDM), PUBLIC :: unitsDM

  REAL(DP), ALLOCATABLE, PUBLIC :: uDM(:,:,:,:,:)

  PUBLIC :: CreateMagnetofluidFields
  PUBLIC :: DestroyMagnetofluidFields
  PUBLIC :: ResetFields_Diagnostic
  PUBLIC :: SetUnitsFields
  PUBLIC :: DescribeFields_Conserved
  PUBLIC :: DescribeFields_Primitive
  PUBLIC :: DescribeFields_Auxiliary
  PUBLIC :: DescribeFields_Diagnostic

CONTAINS


  SUBROUTINE CreateMagnetofluidFields &
    ( nX, swX, CoordinateSystem_Option, Verbose_Option )

    INTEGER,          INTENT(in)           :: nX(3), swX(3)
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: CoordinateSystem_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    CHARACTER(LEN=16) :: CoordinateSystem

    CoordinateSystem = 'CARTESIAN'
    IF( PRESENT( CoordinateSystem_Option ) ) &
      CoordinateSystem = TRIM( CoordinateSystem_Option )

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL CreateFields_Conserved ( nX, swX )
    CALL CreateFields_Primitive ( nX, swX )
    CALL CreateFields_Auxiliary ( nX, swX )
    CALL CreateFields_Diagnostic( nX, swX )

    CALL ResetFields_Diagnostic &
           ( nX, swX, uDM )

    ALLOCATE( rhsCM(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCM) )

    CALL SetUnitsFields( TRIM( CoordinateSystem ), &
                              Verbose_Option = Verbose )

  END SUBROUTINE CreateMagnetofluidFields


  SUBROUTINE CreateFields_Conserved( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFields_Conserved( Verbose )

    ALLOCATE( uCM &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCM) )

    uCM = Zero

  END SUBROUTINE CreateFields_Conserved


  SUBROUTINE DescribeFields_Conserved( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iCM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', ' Fields (Conserved)'
      WRITE(*,*)
      DO iCM = 1, nCM
        WRITE(*,'(A5,A32)') '', TRIM( namesCM(iCM) )
      END DO

    END IF

  END SUBROUTINE DescribeFields_Conserved


  SUBROUTINE CreateFields_Primitive( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFields_Primitive( Verbose )

    ALLOCATE( uPM(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nPM) )

    uPM = Zero

  END SUBROUTINE CreateFields_Primitive


  SUBROUTINE DescribeFields_Primitive( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iPM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', ' Fields (Primitive)'
      WRITE(*,*)
      DO iPM = 1, nPM
        WRITE(*,'(A5,A32)') '', TRIM( namesPM(iPM) )
      END DO

    END IF

  END SUBROUTINE DescribeFields_Primitive


  SUBROUTINE CreateFields_Auxiliary( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFields_Auxiliary( Verbose )

    ALLOCATE( uAM(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nAM) )

    uAM = Zero

  END SUBROUTINE CreateFields_Auxiliary


  SUBROUTINE DescribeFields_Auxiliary( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iAM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A25)') '', ' Fields (Auxiliary)'
      WRITE(*,*)
      DO iAM = 1, nAM
        WRITE(*,'(A5,A32)') '', TRIM( namesAM(iAM) )
      END DO

    END IF

  END SUBROUTINE DescribeFields_Auxiliary


  SUBROUTINE CreateFields_Diagnostic( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFields_Diagnostic( Verbose )

    ALLOCATE( uDM(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nDM) )

    uDM = Zero

  END SUBROUTINE CreateFields_Diagnostic


  SUBROUTINE DescribeFields_Diagnostic( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iDM

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A25)') '', ' Fields (Diagnostic)'
      WRITE(*,*)
      DO iDM = 1, nDM
        WRITE(*,'(A5,A32)') '', TRIM( namesDM(iDM) )
      END DO

    END IF

  END SUBROUTINE DescribeFields_Diagnostic


  SUBROUTINE DestroyMagnetofluidFields

    DEALLOCATE( uCM, rhsCM, uPM, uAM, uDM )

  END SUBROUTINE DestroyMagnetofluidFields


  SUBROUTINE ResetFields_Diagnostic( nX, swX, uDM )

     INTEGER, INTENT(in) :: nX(3), swX(3)

     REAL(DP), INTENT(out) :: uDM( 1:nDOFX, &
                                   1-swX(1):nX(1)+swX(1), &
                                   1-swX(2):nX(2)+swX(2), &
                                   1-swX(3):nX(3)+swX(3), &
                                   1:nDM)

    uDM(:,:,:,:,iDM_TCI)   = Zero
    uDM(:,:,:,:,iDM_Sh_X1) = Zero
    uDM(:,:,:,:,iDM_Sh_X2) = Zero
    uDM(:,:,:,:,iDM_Sh_X3) = Zero
    uDM(:,:,:,:,iDM_T1)    = One
    uDM(:,:,:,:,iDM_T2)    = One
    uDM(:,:,:,:,iDM_T3)    = One
    uDM(:,:,:,:,iDM_Div)   = Zero

  END SUBROUTINE ResetFields_Diagnostic


  SUBROUTINE SetUnitsFields( CoordinateSystem, Verbose_Option )

    USE UnitsModule, ONLY: &
      UnitsActive, &
      Gram, &
      Centimeter, &
      Kilometer, &
      Second, &
      Kelvin, &
      MeV, &
      Erg, &
      BoltzmannConstant, &
      Gauss

    CHARACTER(LEN=*), INTENT(in)           :: CoordinateSystem
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( UnitsActive )THEN

      ! --- Conserved ---

      unitsCM(iCM_D)   = Gram / Centimeter**3

      SELECT CASE( TRIM( CoordinateSystem ) )

        CASE( 'CARTESIAN' )

          unitsCM(iCM_S1) = Gram / Centimeter**2 / Second
          unitsCM(iCM_S2) = Gram / Centimeter**2 / Second
          unitsCM(iCM_S3) = Gram / Centimeter**2 / Second

          unitsCM(iCM_B1)  = Gauss
          unitsCM(iCM_B2)  = Gauss
          unitsCM(iCM_B3)  = Gauss

        CASE( 'CYLINDRICAL' )

          unitsCM(iCM_S1) = Gram / Centimeter**2 / Second
          unitsCM(iCM_S2) = Gram / Centimeter**2 / Second
          unitsCM(iCM_S3) = Gram / Centimeter / Second

          unitsCM(iCM_B1)  = Gauss
          unitsCM(iCM_B2)  = Gauss
          unitsCM(iCM_B3)  = Gauss / Centimeter

        CASE( 'SPHERICAL' )

          unitsCM(iCM_S1) = Gram / Centimeter**2 / Second
          unitsCM(iCM_S2) = Gram / Centimeter / Second
          unitsCM(iCM_S3) = Gram / Centimeter / Second

          unitsCM(iCM_B1)  = Gauss
          unitsCM(iCM_B2)  = Gauss / Centimeter
          unitsCM(iCM_B3)  = Gauss / Centimeter

        CASE DEFAULT

          WRITE(*,*) 'Invalid choice of coordinate system: ', CoordinateSystem
          WRITE(*,*) 'Stopping...'
          STOP

      END SELECT

      unitsCM(iCM_E)   = Erg / Centimeter**3
      unitsCM(iCM_Ne)  = One / Centimeter**3
      unitsCM(iCM_Chi) = Gauss * Kilometer / Second

      ! --- Primitive ---

      unitsPM(iPM_D)  = Gram / Centimeter**3

      SELECT CASE( TRIM( CoordinateSystem ) )

        CASE( 'CARTESIAN' )

          unitsPM(iPM_V1) = Kilometer / Second
          unitsPM(iPM_V2) = Kilometer / Second
          unitsPM(iPM_V3) = Kilometer / Second

          unitsPM(iPM_B1)  = Gauss
          unitsPM(iPM_B2)  = Gauss
          unitsPM(iPM_B3)  = Gauss

        CASE( 'CYLINDRICAL' )

          unitsPM(iPM_V1) = Kilometer / Second
          unitsPM(iPM_V2) = Kilometer / Second
          unitsPM(iPM_V3) = One / Second

          unitsPM(iPM_B1)  = Gauss
          unitsPM(iPM_B2)  = Gauss
          unitsPM(iPM_B3)  = Gauss / Centimeter

        CASE( 'SPHERICAL' )

          unitsPM(iPM_V1) = Kilometer / Second
          unitsPM(iPM_V2) = One / Second
          unitsPM(iPM_V3) = One / Second

          unitsPM(iPM_B1)  = Gauss
          unitsPM(iPM_B2)  = Gauss / Centimeter
          unitsPM(iPM_B3)  = Gauss / Centimeter

        CASE DEFAULT

          WRITE(*,*) 'Invalid choice of coordinate system: ', CoordinateSystem
          WRITE(*,*) 'Stopping...'
          STOP

      END SELECT

      unitsPM(iPM_E)  = Erg / Centimeter**3
      unitsPM(iPM_Ne) = One / Centimeter**3

      unitsPM(iPM_Chi) = Gauss * Kilometer / Second

      ! --- Auxiliary ---

      unitsAM(iAM_P)  = Erg / Centimeter**3
      unitsAM(iAM_Pb) = Erg / Centimeter**3
      unitsAM(iAM_T)  = Kelvin
      unitsAM(iAM_Ye) = One
      unitsAM(iAM_S)  = BoltzmannConstant
      unitsAM(iAM_E)  = Erg / Gram
      unitsAM(iAM_h)  = Erg / Gram
      unitsAM(iAM_hb) = Erg / Gram
      unitsAM(iAM_Me) = MeV
      unitsAM(iAM_Mp) = MeV
      unitsAM(iAM_Mn) = MeV
      unitsAM(iAM_Xp) = One
      unitsAM(iAM_Xn) = One
      unitsAM(iAM_Xa) = One
      unitsAM(iAM_Xh) = One
      unitsAM(iAM_Gm) = One
      unitsAM(iAM_Cs) = Kilometer / Second
      unitsAM(iAM_Ca) = Kilometer / Second
      unitsAM(iAM_Cf_ub_p) = Kilometer / Second
      unitsAM(iAM_Cf_ub_m) = Kilometer / Second
      unitsAM(iAM_Tem00) = Gauss**2 ! Fix later.
      unitsAM(iAM_Tem11) = unitsPM(iPM_B1)**2
      unitsAM(iAM_Tem22) = unitsPM(iPM_B2)**2
      unitsAM(iAM_Tem33) = unitsPM(iPM_B3)**2
      unitsAM(iAM_Tem12) = unitsPM(iPM_B1) * unitsPM(iPM_B2)
      unitsAM(iAM_Tem13) = unitsPM(iPM_B1) * unitsPM(iPM_B3)
      unitsAM(iAM_Tem23) = unitsPM(iPM_B2) * unitsPM(iPM_B3)

      SELECT CASE( TRIM( CoordinateSystem ) )

        CASE( 'CARTESIAN' )

          unitsAM(iAM_EF1) = Gauss ! Need to fix?
          unitsAM(iAM_EF2) = Gauss
          unitsAM(iAM_EF3) = Gauss

        CASE( 'CYLINDRICAL' )

          unitsAM(iAM_EF1) = Gauss
          unitsAM(iAM_EF2) = Gauss
          unitsAM(iAM_EF3) = Gauss / Centimeter

        CASE( 'SPHERICAL' )

          unitsAM(iAM_EF1)  = Gauss
          unitsAM(iAM_EF2)  = Gauss / Centimeter
          unitsAM(iAM_EF3)  = Gauss / Centimeter

        CASE DEFAULT

          WRITE(*,*) 'Invalid choice of coordinate system: ', CoordinateSystem
          WRITE(*,*) 'Stopping...'
          STOP

      END SELECT

      ! --- Diagnostic ---

      unitsDM(iDM_TCI)   = One
      unitsDM(iDM_Sh_X1) = One
      unitsDM(iDM_Sh_X2) = One
      unitsDM(iDM_Sh_X3) = One
      unitsDM(iDM_T1)    = One
      unitsDM(iDM_T2)    = One
      unitsDM(iDM_T3)    = One
      unitsDM(iDM_MinE)  = Erg / Gram
      unitsDM(iDM_MaxE)  = Erg / Gram
      unitsDM(iDM_Div)   = Gauss / Kilometer

    ELSE

      unitsCM = One
      unitsPM = One
      unitsAM = One
      unitsDM = One

    END IF

  END SUBROUTINE SetUnitsFields


END MODULE MagnetofluidFieldsModule
