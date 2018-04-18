MODULE NeutrinoOpacitiesModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: Verbose

  ! --- Equilibrium Distributions ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesEQ = 'Equilibrium Distribution'
  REAL(DP),    ALLOCATABLE, PUBLIC :: f_EQ(:,:,:)

  ! --- Electron Capture Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesEC = 'Electron Capture Opacities'
  REAL(DP),    ALLOCATABLE, PUBLIC :: opEC(:,:,:)

  ! --- Elastic Scattering Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesES = 'Elastic Scattering Opacities'
  REAL(DP),    ALLOCATABLE, PUBLIC :: opES(:,:,:)

  ! --- Inelastic Scattering Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesIS = 'Inelastic Scattering Opacities'
  REAL(DP),    ALLOCATABLE, PUBLIC :: opIS(:,:,:,:)

  ! --- Pair Processes Opacities ---

  CHARACTER(32), PARAMETER, PUBLIC :: namesPP = 'Pair Process Opacities'
  REAL(DP),    ALLOCATABLE, PUBLIC :: opPP(:,:,:,:)

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
      WRITE(*,'(A6,A)') '', TRIM( namesEQ )
      WRITE(*,'(A6,A)') '', TRIM( namesEC )
      WRITE(*,'(A6,A)') '', TRIM( namesES )
    END IF

    ALLOCATE( f_EQ(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )
    ALLOCATE( opEC(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )
    ALLOCATE( opES(nZ(1)*nNodesZ(1),nSpecies,PRODUCT(nZ(2:4)*nNodesZ(2:4))) )

  END SUBROUTINE CreateNeutrinoOpacities


  SUBROUTINE DestroyNeutrinoOpacities

    DEALLOCATE( f_EQ, opEC, opES )

  END SUBROUTINE DestroyNeutrinoOpacities


END MODULE NeutrinoOpacitiesModule
