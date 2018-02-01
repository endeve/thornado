MODULE InputOutputUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE MeshModule, ONLY: &
    MeshType, &
    NodeCoordinate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: NodeCoordinates
  PUBLIC :: Field4D

CONTAINS


  FUNCTION NodeCoordinates( Mesh, nX, nN )

    REAL(DP) :: NodeCoordinates(nX*nN)
    TYPE(MeshType), INTENT(in) :: Mesh
    INTEGER,        INTENT(in) :: nX
    INTEGER,        INTENT(in) :: nN

    INTEGER :: i, j, iNode

    iNode = 0
    DO j = 1, nX
      DO i = 1, nN
        iNode = iNode + 1
        NodeCoordinates(iNode) &
          = NodeCoordinate( Mesh, j, i )
      END DO
    END DO

    RETURN
  END FUNCTION NodeCoordinates


  FUNCTION Field4D( F, nX, nN, nDOF, Tab )

    INTEGER,  INTENT(in) :: &
      nX(4), nN(4), nDOF, Tab(4,nDOF)
    REAL(DP), INTENT(in) :: &
      F(nDOF,nX(1),nX(2),nX(3),nX(4))
    REAL(DP) :: &
      Field4D(nX(1)*nN(1),nX(2)*nN(2),nX(3)*nN(3),nX(4)*nN(4))

    INTEGER :: iX1, iX2, iX3, iX4
    INTEGER :: iN1, iN2, iN3, iN4, iNode

    DO iX4 = 1, nX(4)
      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)

            DO iNode = 1, nDOF

              iN1 = Tab(1,iNode)
              iN2 = Tab(2,iNode)
              iN3 = Tab(3,iNode)
              iN4 = Tab(4,iNode)

              Field4D &
                ( (iX1-1)*nN(1)+iN1, (iX2-1)*nN(2)+iN2, &
                  (iX3-1)*nN(3)+iN3, (iX4-1)*nN(4)+iN4 ) &
                = F(iNode,iX1,iX2,iX3,iX4)

            END DO
          END DO
        END DO
      END DO
    END DO

    RETURN
  END FUNCTION Field4D


END MODULE InputOutputUtilitiesModule
