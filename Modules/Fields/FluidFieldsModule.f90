MODULE FluidFieldsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nNodesX

  IMPLICIT NONE
  PRIVATE

  ! --- Weights to Integrate Fluid Fields ---

  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: WeightsF
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: WeightsF_X1
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: WeightsF_X2
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: WeightsF_X3

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

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uCF, rhsCF

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

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uPF

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

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uAF

  ! --- Diagnostic Variables ---

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: Shock

  PUBLIC :: CreateFluidFields
  PUBLIC :: DestroyFluidFields

CONTAINS


  SUBROUTINE CreateFluidFields( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    ALLOCATE( WeightsF(1:nDOFX) )
    ALLOCATE( WeightsF_X1(1:nNodesX(2)*nNodesX(3)) )
    ALLOCATE( WeightsF_X2(1:nNodesX(1)*nNodesX(3)) )
    ALLOCATE( WeightsF_X3(1:nNodesX(1)*nNodesX(2)) )

    CALL CreateFluidFields_Conserved( nX, swX )
    CALL CreateFluidFields_Primitive( nX, swX )
    CALL CreateFluidFields_Auxiliary( nX, swX )

    ALLOCATE( Shock(1:nX(1),1:nX(2),1:nX(3)) )
    Shock = 0.0_DP

  END SUBROUTINE CreateFluidFields


  SUBROUTINE CreateFluidFields_Conserved( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iCF

    WRITE(*,*)
    WRITE(*,'(A5,A24)') '', 'Fluid Fields (Conserved)'
    WRITE(*,*)
    DO iCF = 1, nCF
      WRITE(*,'(A5,A32)') '', TRIM( namesCF(iCF) )
    END DO

    ALLOCATE( uCF &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCF) )

    ALLOCATE( rhsCF &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCF) )

  END SUBROUTINE CreateFluidFields_Conserved


  SUBROUTINE CreateFluidFields_Primitive( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iPF

    WRITE(*,*)
    WRITE(*,'(A5,A24)') '', 'Fluid Fields (Primitive)'
    WRITE(*,*)
    DO iPF = 1, nPF
      WRITE(*,'(A5,A32)') '', TRIM( namesPF(iPF) )
    END DO

    ALLOCATE( uPF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nPF) )

  END SUBROUTINE CreateFluidFields_Primitive


  SUBROUTINE CreateFluidFields_Auxiliary( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iAF

    WRITE(*,*)
    WRITE(*,'(A5,A24)') '', 'Fluid Fields (Auxiliary)'
    WRITE(*,*)
    DO iAF = 1, nAF
      WRITE(*,'(A5,A32)') '', TRIM( namesAF(iAF) )
    END DO

    ALLOCATE( uAF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nAF) )

  END SUBROUTINE CreateFluidFields_Auxiliary


  SUBROUTINE DestroyFluidFields

    DEALLOCATE( WeightsF, WeightsF_X1, WeightsF_X2, WeightsF_X3 )
    DEALLOCATE( uCF, rhsCF, uPF, uAF )
    DEALLOCATE( Shock )

  END SUBROUTINE DestroyFluidFields


END MODULE FluidFieldsModule
