MODULE MF_MeshModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodesX
  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    nX, &
    swX, &
    xL, &
    xR, &
    iOS_CPP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CreateMesh_MF
  PUBLIC :: DestroyMesh_MF

CONTAINS


  SUBROUTINE CreateMesh_MF( iLevel, MeshX )

    INTEGER,        INTENT(in)  :: iLevel
    TYPE(MeshType), INTENT(out) :: MeshX(3)

    INTEGER :: iDim, nXX(3)

    nXX = nX

    nXX(1) = 2**( iLevel ) * nX(1)
    IF( amrex_spacedim .GT. 1 ) nXX(2) = 2**( iLevel ) * nX(2)
    IF( amrex_spacedim .GT. 2 ) nXX(3) = 2**( iLevel ) * nX(3)

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nXX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim), iOS_Option = iOS_CPP(iDim) )

    END DO

  END SUBROUTINE CreateMesh_MF


  SUBROUTINE DestroyMesh_MF( MeshX )

    TYPE(MeshType), INTENT(inout) :: MeshX(3)

    INTEGER :: iDim

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE DestroyMesh_MF


END MODULE MF_MeshModule
