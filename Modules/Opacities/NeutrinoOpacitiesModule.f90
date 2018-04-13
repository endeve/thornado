MODULE NeutrinoOpacitiesModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  ! --- Electron Capture Opacities ---

  CHARACTER(32), PUBLIC :: namesEC = 'Electron Capture Opacities'

  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE, PUBLIC :: opEC

  ! --- Elastic Scattering Opacities ---

  CHARACTER(32), PUBLIC :: namesES = 'Elastic Scattering Opacities'

  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE, PUBLIC :: opES

  ! --- Inelastic Scattering Opacities ---

  CHARACTER(32), PUBLIC :: namesIS = 'Inelastic Scattering Opacities'

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: opIS

  ! --- Pair Processes Opacities ---

  CHARACTER(32), PUBLIC :: namesPP = 'Pair Process Opacities'

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: opPP


  PUBLIC :: CreateNeutrinoOpacities
  PUBLIC :: DestroyNeutrinoOpacities

CONTAINS


  SUBROUTINE CreateNeutrinoOpacities( nZ, nNodesZ, nSpecies, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER, INTENT(in) :: nZ(4), nNodesZ(4), nSpecies
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'Neutrino Opacities:'
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', TRIM( namesEC )
    END IF

    ALLOCATE( opEC(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )

  END SUBROUTINE CreateNeutrinoOpacities


  SUBROUTINE DestroyNeutrinoOpacities

    DEALLOCATE( opEC )

  END SUBROUTINE DestroyNeutrinoOpacities


END MODULE NeutrinoOpacitiesModule
