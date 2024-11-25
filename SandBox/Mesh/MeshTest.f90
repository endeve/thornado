PROGRAM meshTest

  USE KindModule, ONLY: &
    DP
  USE MeshModule, ONLY: &
    MeshX, MeshE, &
    CreateMesh, &
    DestroyMesh

  IMPLICIT NONE

  INTEGER :: &
    iDim, iNode
  INTEGER, PARAMETER :: &
    nNodes = 4, &
    nE = 20
  INTEGER, DIMENSION(3), PARAMETER :: &
    nX = [ 10, 10, 10 ]
  REAL(DP), PARAMETER :: &
    eL = 0.0d0, &
    eR = 1.0d2
  REAL(DP), DIMENSION(3), PARAMETER :: &
    xL = [ 0.0_DP, 0.0_DP, 0.0_DP ], &
    xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

  ! -- Spatial Grid:

  DO iDim = 1, 3

    CALL CreateMesh( MeshX(iDim), nX(iDim), nNodes, xL(iDim), xR(iDim) )

    WRITE(*,*)
    WRITE(*,'(A4,A7,I1,A8,ES8.2E2)') &
      '', 'Length(', iDim, ')     = ', MeshX(iDim) % Length
    WRITE(*,'(A4,A11,I1,A4,ES8.2E2,A3,ES8.2E2)') &
      '', 'MIN/MAX dx(', iDim, ') = ', &
      MINVAL( MeshX(iDim) % Width ), ' / ', MAXVAL( MeshX(iDim) % Width )
    WRITE(*,'(A4,A11,I1,A4,ES8.2E2,A3,ES8.2E2)') &
      '', 'MIN/MAX xC(', iDim, ') = ', &
      MINVAL( MeshX(iDim) % Center ), ' / ', MAXVAL( MeshX(iDim) % Center )
    WRITE(*,*)
    DO iNode = 1, nNodes
      WRITE(*,'(A6,A6,ES12.4E2,A2,A6,ES12.4E2)') &
        '', 'x_n = ', MeshX(iDim) % Nodes(iNode), &
        '', 'w_n = ', MeshX(iDim) % Weights(iNode)
    END DO

  END DO

  DO iDim = 1, 3

    CALL DestroyMesh( MeshX(iDim) )

  END DO

  ! -- Spectral Grid:

  CALL CreateMesh( MeshE, nE, nNodes, eL, eR, ZoomOption = 1.2385_DP )

  CALL DestroyMesh( MeshE )

END PROGRAM meshTest
