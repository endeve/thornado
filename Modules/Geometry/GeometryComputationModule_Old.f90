MODULE GeometryComputationModule

  USE KindModule, ONLY: &
    DP, One
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    uGF, nGF, iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE GeometryBoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_GeometryX, &
    ApplyBoundaryConditions_Geometry

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeGeometryX
  PUBLIC :: ComputeGeometry

CONTAINS


  SUBROUTINE ComputeGeometryX

    SELECT CASE ( TRIM( CoordinateSystem ) )

      CASE ( 'CARTESIAN' )

        CALL ComputeGeometryX_CARTESIAN

      CASE ( 'SPHERICAL' )

        CALL ComputeGeometryX_SPHERICAL

      CASE ( 'CYLINDRICAL' )

        CALL ComputeGeometryX_CYLINDRICAL

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

    CALL ApplyBoundaryConditions_GeometryX

  END SUBROUTINE ComputeGeometryX


  SUBROUTINE ComputeGeometryX_CARTESIAN

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                  = One
                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                  = One
                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                  = One

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CARTESIAN


  SUBROUTINE ComputeGeometryX_SPHERICAL

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1, X2

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)

            DO iNodeX2 = 1, nNodesX(2)

              X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                  = One
                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                  = X1**2
                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                  = ( X1 * SIN( X2 ) )**2

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeGeometryX_SPHERICAL


  SUBROUTINE ComputeGeometryX_CYLINDRICAL

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                  = One
                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                  = One
                uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) &
                  = X1**2

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeGeometryX_CYLINDRICAL


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

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

    CALL ApplyBoundaryConditions_Geometry

  END SUBROUTINE ComputeGeometry


END MODULE GeometryComputationModule
