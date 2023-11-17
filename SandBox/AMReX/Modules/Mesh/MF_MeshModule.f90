MODULE MF_MeshModule

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nX, &
    swX, &
    xL, &
    xR, &
    nDimsX
  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    iOS_CPP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CreateMesh_MF
  PUBLIC :: DestroyMesh_MF

CONTAINS


  SUBROUTINE CreateMesh_MF( iLevel, MeshX )

    INTEGER,        INTENT(in)  :: iLevel
    TYPE(MeshType), INTENT(out) :: MeshX(3)

    INTEGER :: iDimX, nXX(3)

    nXX = nX

    nXX(1) = 2**( iLevel ) * nX(1)
    IF( nDimsX .GT. 1 ) nXX(2) = 2**( iLevel ) * nX(2)
    IF( nDimsX .GT. 2 ) nXX(3) = 2**( iLevel ) * nX(3)

    DO iDimX = 1, 3

      CALL CreateMesh &
             ( MeshX(iDimX), nXX(iDimX), nNodesX(iDimX), swX(iDimX), &
               xL(iDimX), xR(iDimX), iOS_Option = iOS_CPP(iDimX) )

    END DO

  END SUBROUTINE CreateMesh_MF


  SUBROUTINE DestroyMesh_MF( MeshX )

    TYPE(MeshType), INTENT(inout) :: MeshX(3)

    INTEGER :: iDimX

    DO iDimX = 1, 3

      CALL DestroyMesh( MeshX(iDimX) )

    END DO

  END SUBROUTINE DestroyMesh_MF


END MODULE MF_MeshModule
