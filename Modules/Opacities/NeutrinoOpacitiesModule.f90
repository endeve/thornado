MODULE NeutrinoOpacitiesModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  ! --- Electron Capture Opacities ---

  CHARACTER(32) :: namesEC = 'Electron Capture Opacities'

  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE, PUBLIC :: opEC

  ! --- Elastic Scattering Opacities ---

  CHARACTER(32) :: namesES = 'Elastic Scattering Opacities'

  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE, PUBLIC :: opES

  ! --- Inelastic Scattering Opacities ---

  CHARACTER(32) :: namesIS = 'Inelastic Scattering Opacities'

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: opIS

  ! --- Pair Processes Opacities ---

  CHARACTER(32) :: namesPP = 'Pair Process Opacities'

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: opPP


  PUBLIC :: CreateNeutrinoOpacities
  PUBLIC :: DestroyNeutrinoOpacities

CONTAINS


  SUBROUTINE CreateNeutrinoOpacities( nZ, nNodesZ, nSpecies )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER, INTENT(in) :: nZ(4), nNodesZ(4), nSpecies

    ALLOCATE( opEC(nZ(1)*nNodesZ(1),PRODUCT(nZ(2:4)*nNodesZ(2:4)),nSpecies) )

    PRINT*,"SHAPE(opEC) = ", SHAPE( opEC )

  END SUBROUTINE CreateNeutrinoOpacities


  SUBROUTINE DestroyNeutrinoOpacities

    DEALLOCATE( opEC )

  END SUBROUTINE DestroyNeutrinoOpacities


END MODULE NeutrinoOpacitiesModule
