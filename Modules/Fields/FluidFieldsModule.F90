MODULE FluidFieldsModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  ! --- Conserved Fluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iCF_D  = 1  ! Conserved Baryon Density
  INTEGER, PUBLIC, PARAMETER :: iCF_S1 = 2  ! Conserved Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iCF_S2 = 3  ! Conserved Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iCF_S3 = 4  ! Conserved Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: iCF_E  = 5  ! Conserved Energy Density
  INTEGER, PUBLIC, PARAMETER :: iCF_Ne = 6  ! Conserved Electron Density
  INTEGER, PUBLIC, PARAMETER :: nCF    = 6  ! n Conserved Fluid Fields

  CHARACTER(32), DIMENSION(nCF), PUBLIC, PARAMETER :: &
    namesCF = [ 'Conserved Baryon Density        ', &
                'Conserved Momentum Density (1)  ', &
                'Conserved Momentum Density (2)  ', &
                'Conserved Momentum Density (3)  ', &
                'Conserved Energy Density        ', &
                'Conserved Electron Density      ' ]

  CHARACTER(10),  DIMENSION(nCF), PUBLIC, PARAMETER :: &
    ShortNamesCF = [ 'CF_D      ', &
                     'CF_S1     ', &
                     'CF_S2     ', &
                     'CF_S3     ', &
                     'CF_E      ', &
                     'CF_Ne     ' ]

  REAL(DP), DIMENSION(nCF), PUBLIC :: unitsCF

  REAL(DP), ALLOCATABLE, PUBLIC :: uCF  (:,:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: rhsCF(:,:,:,:,:)

  ! --- Primitive Fluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iPF_D  = 1  ! Comoving Baryon Density
  INTEGER, PUBLIC, PARAMETER :: iPF_V1 = 2  ! Three-Velocity 1
  INTEGER, PUBLIC, PARAMETER :: iPF_V2 = 3  ! Three-Velocity 2
  INTEGER, PUBLIC, PARAMETER :: iPF_V3 = 4  ! Three-Velocity 3
  INTEGER, PUBLIC, PARAMETER :: iPF_E  = 5  ! Internal Energy Density
  INTEGER, PUBLIC, PARAMETER :: iPF_Ne = 6  ! Comoving Electron Density
  INTEGER, PUBLIC, PARAMETER :: nPF    = 6  ! n Primitive Fluid Fields

  CHARACTER(32), DIMENSION(nPF), PUBLIC, PARAMETER :: &
    namesPF = [ 'Comoving Baryon Density         ', &
                'Three-Velocity (1)              ', &
                'Three-Velocity (2)              ', &
                'Three-Velocity (3)              ', &
                'Internal Energy Density         ', &
                'Comoving Electron Density       ' ]

  CHARACTER(10),  DIMENSION(nPF), PUBLIC, PARAMETER :: &
    ShortNamesPF = [ 'PF_D      ', &
                     'PF_V1     ', &
                     'PF_V2     ', &
                     'PF_V3     ', &
                     'PF_E      ', &
                     'PF_Ne     ' ]

  REAL(DP), DIMENSION(nPF), PUBLIC :: unitsPF

  REAL(DP), ALLOCATABLE, PUBLIC :: uPF(:,:,:,:,:)

  ! --- Auxiliary Fluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iAF_P  = 01 ! Pressure
  INTEGER, PUBLIC, PARAMETER :: iAF_T  = 02 ! Temperature
  INTEGER, PUBLIC, PARAMETER :: iAF_Ye = 03 ! Electron Fraction
  INTEGER, PUBLIC, PARAMETER :: iAF_S  = 04 ! Entropy Per Baryon
  INTEGER, PUBLIC, PARAMETER :: iAF_E  = 05 ! Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: iAF_Me = 06 ! Electron Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAF_Mp = 07 ! Proton Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAF_Mn = 08 ! Neutron Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAF_Xp = 09 ! Proton Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAF_Xn = 10 ! Neutron Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAF_Xa = 11 ! Alpha Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAF_Xh = 12 ! Heavy Mass Fraction
  INTEGER, PUBLIC, PARAMETER :: iAF_Gm = 13 ! Ratio of Specific Heats
  INTEGER, PUBLIC, PARAMETER :: iAF_Cs = 14 ! Sound Speed
  INTEGER, PUBLIC, PARAMETER :: nAF    = 14 ! n Auxiliary Fluid Fields

  CHARACTER(32), DIMENSION(nAF), PUBLIC, PARAMETER :: &
    namesAF = [ 'Pressure                        ', &
                'Temperature                     ', &
                'Electron Fraction               ', &
                'Entropy Per Baryon              ', &
                'Specific Internal Energy        ', &
                'Electron Chemical Potential     ', &
                'Proton Chemical Potential       ', &
                'Neutron Chemical Potential      ', &
                'Proton Mass Fraction            ', &
                'Neutron Mass Fraction           ', &
                'Alpha Mass Fraction             ', &
                'Heavy Mass Fraction             ', &
                'Ratio of Specific Heats (Gamma) ', &
                'Sound Speed                     ' ]

  CHARACTER(10),  DIMENSION(nAF), PUBLIC, PARAMETER :: &
    ShortNamesAF = [ 'AF_P      ', &
                     'AF_T      ', &
                     'AF_Ye     ', &
                     'AF_S      ', &
                     'AF_E      ', &
                     'AF_Me     ', &
                     'AF_Mp     ', &
                     'AF_Mn     ', &
                     'AF_Xp     ', &
                     'AF_Xn     ', &
                     'AF_Xa     ', &
                     'AF_Xh     ', &
                     'AF_Gm     ', &
                     'AF_Cs     ' ]

  REAL(DP), DIMENSION(nAF), PUBLIC :: unitsAF

  REAL(DP), ALLOCATABLE, PUBLIC :: uAF(:,:,:,:,:)

  ! --- Diagnostic Variables ---

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: Shock
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: Theta1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: Theta2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: Theta3

  INTEGER, PUBLIC, PARAMETER :: iDF_TCI   = 01 ! Troubled-Cell Indicator
  INTEGER, PUBLIC, PARAMETER :: iDF_Sh_X1 = 02 ! Shock Detector (X1)
  INTEGER, PUBLIC, PARAMETER :: iDF_Sh_X2 = 03 ! Shock Detector (X2)
  INTEGER, PUBLIC, PARAMETER :: iDF_Sh_X3 = 04 ! Shock Detector (X3)
  INTEGER, PUBLIC, PARAMETER :: iDF_T1    = 05 ! Theta 1
  INTEGER, PUBLIC, PARAMETER :: iDF_T2    = 06 ! Theta 2
  INTEGER, PUBLIC, PARAMETER :: iDF_T3    = 07 ! Theta 3
  INTEGER, PUBLIC, PARAMETER :: iDF_MinE  = 08 ! Minimum Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: iDF_MaxE  = 09 ! Maximum Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: nDF       = 09 ! n Diagnostic Fluid Fields

  CHARACTER(32), DIMENSION(nDF), PUBLIC, PARAMETER :: &
    namesDF = [ 'TCI                             ', &
                'Shock (X1)                      ', &
                'Shock (X2)                      ', &
                'Shock (X3)                      ', &
                'Theta 1                         ', &
                'Theta 2                         ', &
                'Theta 3                         ', &
                'Min E                           ', &
                'Max E                           ' ]

  CHARACTER(10), DIMENSION(nDF), PUBLIC, PARAMETER :: &
    ShortNamesDF = [ 'DF_TCI    ', &
                     'DF_Sh_X1  ', &
                     'DF_Sh_X2  ', &
                     'DF_Sh_X3  ', &
                     'DF_T1     ', &
                     'DF_T2     ', &
                     'DF_T3     ', &
                     'DF_MinE   ', &
                     'DF_MaxE   ' ]

  REAL(DP), DIMENSION(nDF), PUBLIC :: unitsDF

  REAL(DP), ALLOCATABLE, PUBLIC :: uDF(:,:,:,:,:)

  PUBLIC :: CreateFluidFields
  PUBLIC :: DestroyFluidFields
  PUBLIC :: ResetFluidFields_Diagnostic
  PUBLIC :: SetUnitsFluidFields
  PUBLIC :: DescribeFluidFields_Conserved
  PUBLIC :: DescribeFluidFields_Primitive
  PUBLIC :: DescribeFluidFields_Auxiliary
  PUBLIC :: DescribeFluidFields_Diagnostic

CONTAINS


  SUBROUTINE CreateFluidFields &
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

    CALL CreateFluidFields_Conserved ( nX, swX )
    CALL CreateFluidFields_Primitive ( nX, swX )
    CALL CreateFluidFields_Auxiliary ( nX, swX )
    CALL CreateFluidFields_Diagnostic( nX, swX )

    CALL ResetFluidFields_Diagnostic &
           ( nX, swX, uDF )

    ALLOCATE( rhsCF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) )

    CALL SetUnitsFluidFields( TRIM( CoordinateSystem ), &
                              Verbose_Option = Verbose )

  END SUBROUTINE CreateFluidFields


  SUBROUTINE CreateFluidFields_Conserved( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFluidFields_Conserved( Verbose )

    ALLOCATE( uCF &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCF) )

    uCF = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uCF )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE(     uCF )
#endif

  END SUBROUTINE CreateFluidFields_Conserved


  SUBROUTINE DescribeFluidFields_Conserved( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iCF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Fluid Fields (Conserved)'
      WRITE(*,*)
      DO iCF = 1, nCF
        WRITE(*,'(A5,A32)') '', TRIM( namesCF(iCF) )
      END DO

    END IF

  END SUBROUTINE DescribeFluidFields_Conserved


  SUBROUTINE CreateFluidFields_Primitive( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFluidFields_Primitive( Verbose )

    ALLOCATE( uPF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nPF) )

    uPF = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uPF )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE(     uPF )
#endif

  END SUBROUTINE CreateFluidFields_Primitive


  SUBROUTINE DescribeFluidFields_Primitive( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iPF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Fluid Fields (Primitive)'
      WRITE(*,*)
      DO iPF = 1, nPF
        WRITE(*,'(A5,A32)') '', TRIM( namesPF(iPF) )
      END DO

    END IF

  END SUBROUTINE DescribeFluidFields_Primitive


  SUBROUTINE CreateFluidFields_Auxiliary( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFluidFields_Auxiliary( Verbose )

    ALLOCATE( uAF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nAF) )

    uAF = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uAF )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE(     uAF )
#endif

  END SUBROUTINE CreateFluidFields_Auxiliary


  SUBROUTINE DescribeFluidFields_Auxiliary( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iAF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Fluid Fields (Auxiliary)'
      WRITE(*,*)
      DO iAF = 1, nAF
        WRITE(*,'(A5,A32)') '', TRIM( namesAF(iAF) )
      END DO

    END IF

  END SUBROUTINE DescribeFluidFields_Auxiliary


  SUBROUTINE CreateFluidFields_Diagnostic( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    CALL DescribeFluidFields_Diagnostic( Verbose )

    ALLOCATE( uDF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nDF) )

    uDF = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uDF )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE(     uDF )
#endif

  END SUBROUTINE CreateFluidFields_Diagnostic


  SUBROUTINE DescribeFluidFields_Diagnostic( Verbose )

    LOGICAL, INTENT(in) :: Verbose

    INTEGER :: iDF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A5,A25)') '', 'Fluid Fields (Diagnostic)'
      WRITE(*,*)
      DO iDF = 1, nDF
        WRITE(*,'(A5,A32)') '', TRIM( namesDF(iDF) )
      END DO

    END IF

  END SUBROUTINE DescribeFluidFields_Diagnostic


  SUBROUTINE DestroyFluidFields

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: uCF, uPF, uAF, uDF, &
    !$OMP               unitsCF, unitsPF, unitsAF, unitsDF )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE(       uCF, uPF, uAF, uDF, &
    !$ACC               unitsCF, unitsPF, unitsAF, unitsDF )
#endif

    DEALLOCATE( uCF, rhsCF, uPF, uAF, uDF )

  END SUBROUTINE DestroyFluidFields


  SUBROUTINE ResetFluidFields_Diagnostic( nX, swX, uDF )

     INTEGER, INTENT(in) :: nX(3), swX(3)

     REAL(DP), INTENT(out) :: uDF( 1:nDOFX, &
                                   1-swX(1):nX(1)+swX(1), &
                                   1-swX(2):nX(2)+swX(2), &
                                   1-swX(3):nX(3)+swX(3), &
                                   1:nDF)

    uDF(:,:,:,:,iDF_TCI)   = Zero
    uDF(:,:,:,:,iDF_Sh_X1) = Zero
    uDF(:,:,:,:,iDF_Sh_X2) = Zero
    uDF(:,:,:,:,iDF_Sh_X3) = Zero
    uDF(:,:,:,:,iDF_T1)    = One
    uDF(:,:,:,:,iDF_T2)    = One
    uDF(:,:,:,:,iDF_T3)    = One

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( uDF )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( uDF )
#endif

  END SUBROUTINE ResetFluidFields_Diagnostic


  SUBROUTINE SetUnitsFluidFields( CoordinateSystem, Verbose_Option )

    USE UnitsModule, ONLY: &
      UnitsActive, &
      Gram, &
      Centimeter, &
      Kilometer, &
      Second, &
      Kelvin, &
      MeV, &
      Erg, &
      BoltzmannConstant

    CHARACTER(LEN=*), INTENT(in)           :: CoordinateSystem
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN

      WRITE(*,*)
#if   defined HYDRO_RIEMANN_SOLVER_HLL
      WRITE(*,'(5x,A)') 'Fluid Riemann Solver: HLL'
#elif defined HYDRO_RIEMANN_SOLVER_HLLC
      WRITE(*,'(5x,A)') 'Fluid Riemann Solver: HLLC'
#elif defined HYDRO_RIEMANN_SOLVER_HYBRID
      WRITE(*,'(5x,A)') 'Fluid Riemann Solver: HYBRID'
#endif

  END IF

    IF( UnitsActive )THEN

      ! --- Conserved ---

      unitsCF(iCF_D)  = Gram / Centimeter**3

      SELECT CASE( TRIM( CoordinateSystem ) )

        CASE( 'CARTESIAN' )

          unitsCF(iCF_S1) = Gram / Centimeter**2 / Second
          unitsCF(iCF_S2) = Gram / Centimeter**2 / Second
          unitsCF(iCF_S3) = Gram / Centimeter**2 / Second

        CASE( 'CYLINDRICAL' )

          unitsCF(iCF_S1) = Gram / Centimeter**2 / Second
          unitsCF(iCF_S2) = Gram / Centimeter**2 / Second
          unitsCF(iCF_S3) = Gram / Centimeter / Second

        CASE( 'SPHERICAL' )

          unitsCF(iCF_S1) = Gram / Centimeter**2 / Second
          unitsCF(iCF_S2) = Gram / Centimeter / Second
          unitsCF(iCF_S3) = Gram / Centimeter / Second

        CASE DEFAULT

          WRITE(*,*) 'Invalid choice of coordinate system: ', CoordinateSystem
          WRITE(*,*) 'Stopping...'
          STOP

      END SELECT

      unitsCF(iCF_E)  = Erg / Centimeter**3
      unitsCF(iCF_Ne) = One / Centimeter**3

      ! --- Primitive ---

      unitsPF(iPF_D)  = Gram / Centimeter**3

      SELECT CASE( TRIM( CoordinateSystem ) )

        CASE( 'CARTESIAN' )

          unitsPF(iPF_V1) = Kilometer / Second
          unitsPF(iPF_V2) = Kilometer / Second
          unitsPF(iPF_V3) = Kilometer / Second

        CASE( 'CYLINDRICAL' )

          unitsPF(iPF_V1) = Kilometer / Second
          unitsPF(iPF_V2) = Kilometer / Second
          unitsPF(iPF_V3) = One / Second

        CASE( 'SPHERICAL' )

          unitsPF(iPF_V1) = Kilometer / Second
          unitsPF(iPF_V2) = One / Second
          unitsPF(iPF_V3) = One / Second

        CASE DEFAULT

          WRITE(*,*) 'Invalid choice of coordinate system: ', CoordinateSystem
          WRITE(*,*) 'Stopping...'
          STOP

      END SELECT

      unitsPF(iPF_E)  = Erg / Centimeter**3
      unitsPF(iPF_Ne) = One / Centimeter**3

      ! --- Auxiliary ---

      unitsAF(iAF_P)  = Erg / Centimeter**3
      unitsAF(iAF_T)  = Kelvin
      unitsAF(iAF_Ye) = One
      unitsAF(iAF_S)  = BoltzmannConstant
      unitsAF(iAF_E)  = Erg / Gram
      unitsAF(iAF_Me) = MeV
      unitsAF(iAF_Mp) = MeV
      unitsAF(iAF_Mn) = MeV
      unitsAF(iAF_Xp) = One
      unitsAF(iAF_Xn) = One
      unitsAF(iAF_Xa) = One
      unitsAF(iAF_Xh) = One
      unitsAF(iAF_Gm) = One
      unitsAF(iAF_Cs) = Kilometer / Second

      ! --- Diagnostic ---
      unitsDF(iDF_TCI)   = One
      unitsDF(iDF_Sh_X1) = One
      unitsDF(iDF_Sh_X2) = One
      unitsDF(iDF_Sh_X3) = One
      unitsDF(iDF_T1)    = One
      unitsDF(iDF_T2)    = One
      unitsDF(iDF_T3)    = One
      unitsDF(iDF_MinE)  = Erg / Gram
      unitsDF(iDF_MaxE)  = Erg / Gram

    ELSE

      unitsCF = One
      unitsPF = One
      unitsAF = One
      unitsDF = One

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: unitsCF, unitsPF, unitsAF, unitsDF)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( unitsCF, unitsPF, unitsAF, unitsDF )
#endif

  END SUBROUTINE SetUnitsFluidFields


END MODULE FluidFieldsModule
