MODULE GeometryBoundaryConditionsModule

  USE KindModule, ONLY: &
    DP, One
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, swX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Up
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    uGF, nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, iGF_Beta_1
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_GeometryX
  PUBLIC :: ApplyBoundaryConditions_Geometry
  PUBLIC :: ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting
  PUBLIC :: ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_GeometryX

    CALL ApplyBoundaryConditions_GeometryX_X1

  END SUBROUTINE ApplyBoundaryConditions_GeometryX


  SUBROUTINE ApplyBoundaryConditions_GeometryX_X1

    IF( swX(1) < 1 ) RETURN

    SELECT CASE ( TRIM( CoordinateSystem ) )

      CASE ( 'CARTESIAN' )

        CALL ApplyBoundaryConditions_GeometryX_X1_CARTESIAN

      CASE ( 'SPHERICAL' )

        CALL ApplyBoundaryConditions_GeometryX_X1_SPHERICAL

      CASE ( 'CYLINDRICAL' )

        CALL ApplyBoundaryConditions_GeometryX_X1_CYLINDRICAL

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_GeometryX_X1


  SUBROUTINE ApplyBoundaryConditions_GeometryX_X1_CARTESIAN

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              ! --- Inner Boundary ---

              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_11) &
                = One
              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_22) &
                = One
              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_33) &
                = One

              ! --- Outer Boundary ---

              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_11) &
                = One
              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_22) &
                = One
              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_33) &
                = One

            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_GeometryX_X1_CARTESIAN


  SUBROUTINE ApplyBoundaryConditions_GeometryX_X1_SPHERICAL

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1, X2

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)

            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              ! --- Inner Boundary ---

              X1 = NodeCoordinate( MeshX(1), 0, iNodeX1 )

              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_11) &
                = One
              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_22) &
                = X1**2
              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_33) &
                = ( X1 * SIN( X2 ) )**2

              ! --- Outer Boundary ---

              X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_11) &
                = One
              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_22) &
                = X1**2
              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_33) &
                = ( X1 * SIN( X2 ) )**2

            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_GeometryX_X1_SPHERICAL


  SUBROUTINE ApplyBoundaryConditions_GeometryX_X1_CYLINDRICAL

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              ! --- Inner Boundary ---

              X1 = NodeCoordinate( MeshX(1), 0, iNodeX1 )

              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_11) &
                = One
              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_22) &
                = One
              uGF(iNodeX,0,iX2,iX3,iGF_Gm_dd_33) &
                = X1**2

              ! --- Outer Boundary ---

              X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_11) &
                = One
              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_22) &
                = One
              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Gm_dd_33) &
                = X1**2

            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_GeometryX_X1_CYLINDRICAL


  SUBROUTINE ApplyBoundaryConditions_Geometry

    CALL ApplyBoundaryConditions_Geometry_X1

  END SUBROUTINE ApplyBoundaryConditions_Geometry


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1

    INTEGER  :: iE, iX2, iX3
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

                END DO
              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNX1, iNX2, iNX3, iNX
    INTEGER :: jNX1, jNX
    INTEGER :: iGF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1       , swX  (1)

      DO iNX3 = 1, nNodesX(3)
      DO iNX2 = 1, nNodesX(2)
      DO iNX1 = 1, nNodesX(1)

        jNX1 = ( nNodesX(1) - iNX1 ) + 1

        iNX = NodeNumberX( iNX1, iNX2, iNX3 )
        jNX = NodeNumberX( jNX1, iNX2, iNX3 )

        DO iGF = 1, nGF

          G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF) &
            = +G(jNX,iX_B0(1),iX2,iX3,iGF)

        END DO

        G(iNX,iX_B0(1)-iX1,iX2,iX3,iGF_Beta_1) &
          = -G(jNX,iX_B0(1),iX2,iX3,iGF_Beta_1)

      END DO
      END DO
      END DO

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX_B0(1)-iX1,iX2,iX3,:) )

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting


  SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    REAL(DP) :: G_F(nDOFX_X1,nGF)

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iGF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = 1       , swX  (1)

      DO iGF = 1, nGF

        G_F(:,iGF) = MATMUL( LX_X1_Up, G(:,iX_E0(1),iX2,iX3,iGF) )

        G(:,iX_E0(1)+iX1,iX2,iX3,iGF) = G_F(1,iGF)

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G(:,iX_E0(1)+iX1,iX2,iX3,:) )

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace


END MODULE GeometryBoundaryConditionsModule
