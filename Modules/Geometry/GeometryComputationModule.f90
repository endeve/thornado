MODULE GeometryComputationModule

  USE KindModule, ONLY: &
    DP
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    WeightsGX, VolX, VolJacX, &
    WeightsG,  Vol,  VolJac, VolJacE, d
  USE GeometryBoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_GeometryX, &
    ApplyBoundaryConditions_Geometry

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX
  PUBLIC :: ComputeGeometry

CONTAINS


  SUBROUTINE ComputeGeometryX( nX, nNodesX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, nNodesX, swX

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1, X2, X3

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)

            X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

            DO iNodeX2 = 1, nNodesX(2)

              X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                VolJacX(iNodeX,iX1,iX2,iX3) = d( [ X1, X2, X3 ] )

              END DO
            END DO
          END DO

          VolX(iX1,iX2,iX3) &
            = DOT_PRODUCT( WeightsGX(:), VolJacX(:,iX1,iX2,iX3) )

        END DO
      END DO
    END DO

    CALL ApplyBoundaryConditions_GeometryX

  END SUBROUTINE ComputeGeometryX


  SUBROUTINE ComputeGeometry( nX, nNodesX, swX, nE, nNodesE, swE )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, nNodesX, swX
    INTEGER,               INTENT(in) :: nE, nNodesE, swE

    INTEGER  :: iE, iX1, iX2, iX3
    INTEGER  :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: E, X1, X2, X3

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)

              X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

              DO iNodeX2 = 1, nNodesX(2)

                X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                DO iNodeX1 = 1, nNodesX(1)

                  X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    VolJac(iNode,iE,iX1,iX2,iX3) &
                      = d( [ X1, X2, X3 ] ) * E**2

                    VolJacE(iNode,iE,iX1,iX2,iX3) &
                      = d( [ X1, X2, X3 ] ) * E**3

                  END DO
                END DO
              END DO
            END DO

            Vol(iE,iX1,iX2,iX3) &
              = DOT_PRODUCT( WeightsG(:), VolJac(:,iE,iX1,iX2,iX3) )

          END DO
        END DO
      END DO
    END DO

    CALL ApplyBoundaryConditions_Geometry

  END SUBROUTINE ComputeGeometry


END MODULE GeometryComputationModule
