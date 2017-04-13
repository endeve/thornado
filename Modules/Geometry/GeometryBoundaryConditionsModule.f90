MODULE GeometryBoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, swX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    WeightsGX, VolX, VolJacX, &
    WeightsG,  Vol,  VolJac, d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_GeometryX
  PUBLIC :: ApplyBoundaryConditions_Geometry

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_GeometryX

    CALL ApplyBoundaryConditions_GeometryX_X1

  END SUBROUTINE ApplyBoundaryConditions_GeometryX


  SUBROUTINE ApplyBoundaryConditions_GeometryX_X1

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1, X2, X3

    IF( swX(1) < 1 ) RETURN

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        DO iNodeX3 = 1, nNodesX(3)

          X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

          DO iNodeX2 = 1, nNodesX(2)

            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              ! --- Inner Boundary ---

              X1 = NodeCoordinate( MeshX(1), 0, iNodeX1 )

              VolJacX(iNodeX,0,iX2,iX3) = d( [ ABS( X1 ), X2, X3 ] )

              ! --- Outer Boundary ---

              X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

              VolJacX(iNodeX,nX(1)+1,iX2,iX3) = d( [ X1, X2, X3 ] )

            END DO
          END DO
        END DO

        ! --- Inner Boundary ---

        VolX(0,      iX2,iX3) &
          = DOT_PRODUCT( WeightsGX(:), VolJacX(:,0,      iX2,iX3) )

        ! --- Outer Boundary ---

        VolX(nX(1)+1,iX2,iX3) &
          = DOT_PRODUCT( WeightsGX(:), VolJacX(:,nX(1)+1,iX2,iX3) )

      END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_GeometryX_X1


  SUBROUTINE ApplyBoundaryConditions_Geometry

    CALL ApplyBoundaryConditions_Geometry_X1

  END SUBROUTINE ApplyBoundaryConditions_Geometry


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1

    INTEGER  :: iE, iX1, iX2, iX3
    INTEGER  :: iNodeE, iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: E, X1_Inner, X1_Outer, X2, X3

    IF( swX(1) < 1 ) RETURN

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iE = 1, nE

          DO iNodeX3 = 1, nNodesX(3)

            X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

            DO iNodeX2 = 1, nNodesX(2)

              X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

              DO iNodeX1 = 1, nNodesX(1)

                X1_Inner = NodeCoordinate( MeshX(1), 0,       iNodeX1 )
                X1_Outer = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

                DO iNodeE = 1, nNodesE

                  E = NodeCoordinate( MeshE, iE, iNodeE )

                  iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                  ! --- Inner Boundary ---

                  VolJac(iNode,iE,0,      iX2,iX3) &
                    = d( [ ABS( X1_Inner ), X2, X3 ] ) * E**2

                  ! --- Outer Boundary ---

                  VolJac(iNode,iE,nX(1)+1,iX2,iX3) &
                    = d( [ X1_Outer,        X2, X3 ] ) * E**2

                END DO
              END DO
            END DO
          END DO

          ! --- Inner Boundary ---

          Vol(iE,0,      iX2,iX3) &
            = DOT_PRODUCT( WeightsG(:), VolJac(:,iE,0,      iX2,iX3) )

          ! --- Outer Boundary ---

          Vol(iE,nX(1)+1,iX2,iX3) &
            = DOT_PRODUCT( WeightsG(:), VolJac(:,iE,nX(1)+1,iX2,iX3) )

        END DO
      END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1


END MODULE GeometryBoundaryConditionsModule
