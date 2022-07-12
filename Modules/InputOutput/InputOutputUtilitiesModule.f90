MODULE InputOutputUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE MeshModule, ONLY: &
    MeshType, &
    NodeCoordinate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: NodeCoordinates
  PUBLIC :: Field3D
  PUBLIC :: Field3D_INT
  PUBLIC :: FromField3D
  PUBLIC :: Field4D
  PUBLIC :: FromField4D
  PUBLIC :: Opacity4D

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


  FUNCTION Field3D( F, nX, nN, nDOF, Tab )

    INTEGER,  INTENT(in) :: &
      nX(3), nN(3), nDOF, Tab(3,nDOF)
    REAL(DP), INTENT(in) :: &
      F(nDOF,nX(1),nX(2),nX(3))
    REAL(DP) :: &
      Field3D(nX(1)*nN(1),nX(2)*nN(2),nX(3)*nN(3))

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iN1, iN2, iN3, iNode

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNode = 1, nDOF

        iN1 = Tab(1,iNode)
        iN2 = Tab(2,iNode)
        iN3 = Tab(3,iNode)

        Field3D &
          ( (iX1-1)*nN(1)+iN1, (iX2-1)*nN(2)+iN2, (iX3-1)*nN(3)+iN3 ) &
          = F(iNode,iX1,iX2,iX3)

      END DO

    END DO
    END DO
    END DO

    RETURN
  END FUNCTION Field3D


  FUNCTION Field3D_INT( F, nX, nN, nDOF, Tab )

    INTEGER, INTENT(in) :: &
      nX(3), nN(3), nDOF, Tab(3,nDOF)
    INTEGER, INTENT(in) :: &
      F(nDOF,nX(1),nX(2),nX(3))
    INTEGER             :: &
      Field3D_INT(nX(1)*nN(1),nX(2)*nN(2),nX(3)*nN(3))

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iN1, iN2, iN3, iNode

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNode = 1, nDOF

        iN1 = Tab(1,iNode)
        iN2 = Tab(2,iNode)
        iN3 = Tab(3,iNode)

        Field3D_INT &
          ( (iX1-1)*nN(1)+iN1, (iX2-1)*nN(2)+iN2, (iX3-1)*nN(3)+iN3 ) &
          = F(iNode,iX1,iX2,iX3)

      END DO

    END DO
    END DO
    END DO

    RETURN
  END FUNCTION Field3D_INT


  FUNCTION FromField3D( F, nX, nN, nDOF, Tab )

    INTEGER,  INTENT(in) :: &
      nX(3), nN(3), nDOF, Tab(3,nDOF)
    REAL(DP), INTENT(in) :: &
      F(nX(1)*nN(1),nX(2)*nN(2),nX(3)*nN(3))
    REAL(DP) :: &
      FromField3D(nDOF,nX(1),nX(2),nX(3))

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iN1, iN2, iN3, iNode

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNode = 1, nDOF

        iN1 = Tab(1,iNode)
        iN2 = Tab(2,iNode)
        iN3 = Tab(3,iNode)

        FromField3D(iNode,iX1,iX2,iX3) &
          = F( (iX1-1)*nN(1)+iN1, (iX2-1)*nN(2)+iN2, (iX3-1)*nN(3)+iN3 )

      END DO

    END DO
    END DO
    END DO

    RETURN
  END FUNCTION FromField3D


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


  FUNCTION FromField4D( F, nX, nN, nDOF, Tab )

    INTEGER,  INTENT(in) :: &
      nX(4), nN(4), nDOF, Tab(4,nDOF)
    REAL(DP), INTENT(in) :: &
      F(nX(1)*nN(1),nX(2)*nN(2),nX(3)*nN(3),nX(4)*nN(4))
    REAL(DP) :: &
      FromField4D(nDOF,nX(1),nX(2),nX(3),nX(4))

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

        FromField4D(iNode,iX1,iX2,iX3,iX4) &
          = F( (iX1-1)*nN(1)+iN1, (iX2-1)*nN(2)+iN2, &
               (iX3-1)*nN(3)+iN3, (iX4-1)*nN(4)+iN4 )

      END DO

    END DO
    END DO
    END DO
    END DO

    RETURN
  END FUNCTION FromField4D


  FUNCTION Opacity4D( O, nE, nNE, nDOFE, nX, nNX, nDOFX, TabX )

    INTEGER,  INTENT(in) :: &
      nE,    nNE,    nDOFE
    INTEGER,  INTENT(in) :: &
      nX(3), nNX(3), nDOFX, TabX(3,nDOFX)
    REAL(DP), INTENT(in) :: &
      O(nE*nNE,PRODUCT( nX*nNX ))
    REAL(DP)             :: &
      Opacity4D(nE*nNE,nX(1)*nNX(1),nX(2)*nNX(2),nX(3)*nNX(3))

    INTEGER :: iE, iX1, iX2, iX3, iOS_X
    INTEGER :: iNX1, iNX2, iNX3, iNodeX

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)

      DO iNodeX = 1, nDOFX

        iNX1 = TabX(1,iNodeX)
        iNX2 = TabX(2,iNodeX)
        iNX3 = TabX(3,iNodeX)

        iOS_X = ( (iX3-1)*nX(2)*nX(1) + (iX2-1)*nX(1) + (iX1-1) ) * nDOFX

        Opacity4D &
          (:,(iX1-1)*nNX(1)+iNX1, (iX2-1)*nNX(2)+iNX2, (iX3-1)*nNX(3)+iNX3) &
          = O(:,iOS_X+iNodeX)

      END DO

    END DO
    END DO
    END DO

    RETURN
  END FUNCTION Opacity4D


END MODULE InputOutputUtilitiesModule
