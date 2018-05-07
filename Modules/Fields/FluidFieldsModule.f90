MODULE FluidFieldsModule

  USE KindModule, ONLY: &
    DP
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

  REAL(DP), DIMENSION(nAF), PUBLIC :: unitsAF

  REAL(DP), ALLOCATABLE, PUBLIC :: uAF(:,:,:,:,:)

  ! --- Diagnostic Variables ---

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: Shock

  PUBLIC :: CreateFluidFields
  PUBLIC :: DestroyFluidFields

CONTAINS


  SUBROUTINE CreateFluidFields( nX, swX, Verbose_Option )

    INTEGER, INTENT(in)           :: nX(3), swX(3)
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL CreateFluidFields_Conserved( nX, swX )
    CALL CreateFluidFields_Primitive( nX, swX )
    CALL CreateFluidFields_Auxiliary( nX, swX )

    ALLOCATE( rhsCF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) )

    ALLOCATE( Shock(1:nX(1),1:nX(2),1:nX(3)) )
    Shock = 0.0_DP

    CALL SetUnitsFluidFields

  END SUBROUTINE CreateFluidFields


  SUBROUTINE CreateFluidFields_Conserved( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    INTEGER :: iCF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Fluid Fields (Conserved)'
      WRITE(*,*)
      DO iCF = 1, nCF
        WRITE(*,'(A5,A32)') '', TRIM( namesCF(iCF) )
      END DO
    END IF

    ALLOCATE( uCF &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCF) )

  END SUBROUTINE CreateFluidFields_Conserved


  SUBROUTINE CreateFluidFields_Primitive( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    INTEGER :: iPF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Fluid Fields (Primitive)'
      WRITE(*,*)
      DO iPF = 1, nPF
        WRITE(*,'(A5,A32)') '', TRIM( namesPF(iPF) )
      END DO
    END IF

    ALLOCATE( uPF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nPF) )

  END SUBROUTINE CreateFluidFields_Primitive


  SUBROUTINE CreateFluidFields_Auxiliary( nX, swX )

    INTEGER, INTENT(in) :: nX(3), swX(3)

    INTEGER :: iAF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A24)') '', 'Fluid Fields (Auxiliary)'
      WRITE(*,*)
      DO iAF = 1, nAF
        WRITE(*,'(A5,A32)') '', TRIM( namesAF(iAF) )
      END DO
    END IF

    ALLOCATE( uAF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nAF) )

  END SUBROUTINE CreateFluidFields_Auxiliary


  SUBROUTINE DestroyFluidFields

    DEALLOCATE( uCF, rhsCF, uPF, uAF )
    DEALLOCATE( Shock )

  END SUBROUTINE DestroyFluidFields


  SUBROUTINE SetUnitsFluidFields

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

    IF( UnitsActive )THEN

      ! --- Conserved ---

      unitsCF(iCF_D)  = Gram / Centimeter**3
      unitsCF(iCF_S1) = Gram / Centimeter**2 / Second
      unitsCF(iCF_S2) = Gram / Centimeter**2 / Second
      unitsCF(iCF_S3) = Gram / Centimeter**2 / Second
      unitsCF(iCF_E)  = Erg / Centimeter**3
      unitsCF(iCF_Ne) = 1.0_DP / Centimeter**3

      ! --- Primitive ---

      unitsPF(iPF_D)  = Gram / Centimeter**3
      unitsPF(iPF_V1) = Kilometer / Second
      unitsPF(iPF_V2) = Kilometer / Second
      unitsPF(iPF_V3) = Kilometer / Second
      unitsPF(iPF_E)  = Erg / Centimeter**3
      unitsPF(iPF_Ne) = 1.0_DP / Centimeter**3

      ! --- Auxiliary ---

      unitsAF(iAF_P)  = Erg / Centimeter**3
      unitsAF(iAF_T)  = Kelvin
      unitsAF(iAF_Ye) = 1.0_DP
      unitsAF(iAF_S)  = BoltzmannConstant
      unitsAF(iAF_E)  = Erg / Gram
      unitsAF(iAF_Me) = MeV
      unitsAF(iAF_Mp) = MeV
      unitsAF(iAF_Mn) = MeV
      unitsAF(iAF_Xp) = 1.0_DP
      unitsAF(iAF_Xn) = 1.0_DP
      unitsAF(iAF_Xa) = 1.0_DP
      unitsAF(iAF_Xh) = 1.0_DP
      unitsAF(iAF_Gm) = 1.0_DP
      unitsAF(iAF_Cs) = Kilometer / Second

    ELSE

      unitsCF = 1.0_DP
      unitsPF = 1.0_DP
      unitsAF = 1.0_DP

    END IF

  END SUBROUTINE SetUnitsFluidFields


END MODULE FluidFieldsModule
