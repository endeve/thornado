MODULE NeutrinoOpacitiesModule

  USE KindModule, ONLY: &
    DP, One
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  ! --- Equilibrium Distributions ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesEQ = 'Equilibrium Distribution'
  REAL(DP),                 PUBLIC :: unitsEQ
  REAL(DP),    ALLOCATABLE, PUBLIC :: f_EQ(:,:,:)

  ! --- Electron Capture Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesEC = 'Electron Capture Opacities'
  REAL(DP),                 PUBLIC :: unitsEC
  REAL(DP),    ALLOCATABLE, PUBLIC :: opEC(:,:,:)

  ! --- Elastic Scattering Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesES = 'Elastic Scattering Opacities'
  REAL(DP),                 PUBLIC :: unitsES
  REAL(DP),    ALLOCATABLE, PUBLIC :: opES(:,:,:)

  ! --- Inelastic Scattering Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesIS = 'Inelastic Scattering Opacities'
  REAL(DP),                 PUBLIC :: unitsIS
  REAL(DP),    ALLOCATABLE, PUBLIC :: opIS(:,:,:,:)

  ! --- Pair Processes Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesPP = 'Pair Process Opacities'
  REAL(DP),                 PUBLIC :: unitsPP
  REAL(DP),    ALLOCATABLE, PUBLIC :: opPP(:,:,:,:)

  ! --- Nucleon-nucleon Bremsstrahlung Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesBrem = 'Bremsstrahlung Opacities'
  REAL(DP),                 PUBLIC :: unitsBrem
  REAL(DP),    ALLOCATABLE, PUBLIC :: opBrem(:,:,:,:)

  PUBLIC :: CreateNeutrinoOpacities
  PUBLIC :: DestroyNeutrinoOpacities

CONTAINS


  SUBROUTINE CreateNeutrinoOpacities( nZ, nNodesZ, nSpecies, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER, INTENT(in) :: nZ(4), nNodesZ(4), nSpecies
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'Neutrino Opacities:'
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', TRIM( namesEQ )
      WRITE(*,'(A6,A)') '', TRIM( namesEC )
      WRITE(*,'(A6,A)') '', TRIM( namesES )
    END IF

    ALLOCATE( f_EQ(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )
    ALLOCATE( opEC(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )
    ALLOCATE( opES(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )

    CALL SetUnitsNeutrinoOpacities

  END SUBROUTINE CreateNeutrinoOpacities


  SUBROUTINE DestroyNeutrinoOpacities

    DEALLOCATE( f_EQ, opEC, opES )

  END SUBROUTINE DestroyNeutrinoOpacities


  SUBROUTINE SetUnitsNeutrinoOpacities

    USE UnitsModule, ONLY: &
      UnitsActive, &
      Centimeter

    IF( UnitsActive )THEN

      unitsEQ   = One
      unitsEC   = One / Centimeter
      unitsES   = One / Centimeter
      unitsIS   = One ! --- Not Set Yet
      unitsPP   = One / ( Centimeter * MeV**3 )
      unitsBrem = One / ( Centimeter * MeV**3 )

    ELSE

      unitsEQ   = One
      unitsEC   = One
      unitsES   = One
      unitsIS   = One
      unitsPP   = One
      unitsBrem = One

    END IF

  END SUBROUTINE SetUnitsNeutrinoOpacities


END MODULE NeutrinoOpacitiesModule
