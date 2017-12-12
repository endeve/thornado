MODULE GeometryFieldsModuleE

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFE, nNodesE

  IMPLICIT NONE
  PRIVATE

  ! --- 1D Momentum Space (Energy) Geometry Fields ---

  INTEGER, PUBLIC, PARAMETER :: iGE_Esq = 1 ! E^2
  INTEGER, PUBLIC, PARAMETER :: iGE_Ecb = 2 ! E^3
  INTEGER, PUBLIC, PARAMETER :: nGE     = 2

  CHARACTER(32), DIMENSION(nGE), PUBLIC, PARAMETER :: &
    namesGE = [ 'Energy Squared                              ', &
                'Energy Cubed                                ' ]

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: uGE

  PUBLIC :: CreateGeometryFieldsE
  PUBLIC :: DestroyGeometryFieldsE

CONTAINS


  SUBROUTINE CreateGeometryFieldsE( nE, swE )

    INTEGER, INTENT(in) :: nE, swE

    INTEGER :: iGE

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'Geometry Fields (Energy)'
    WRITE(*,*)
    DO iGE = 1, nGE
      WRITE(*,'(A5,A32)') '', TRIM( namesGE(iGE) )
    END DO

    ALLOCATE( uGE(1:nDOFE,1-swE:nE+swE,1:nGE) )

    uGE(:,:,iGE_Esq) = 0.0_DP
    uGE(:,:,iGE_Ecb) = 0.0_DP

  END SUBROUTINE CreateGeometryFieldsE


  SUBROUTINE DestroyGeometryFieldsE

    DEALLOCATE( uGE )

  END SUBROUTINE DestroyGeometryFieldsE


END MODULE GeometryFieldsModuleE
