MODULE GravitySolutionModule_Newtonian_PointMass_Old

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    GravitationalConstant, &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Phi_N

  IMPLICIT NONE
  PRIVATE

  REAL(DP) :: PointMass

  PUBLIC :: InitializeGravitySolver_Newtonian_PointMass
  PUBLIC :: FinalizeGravitySolver_Newtonian_PointMass
  PUBLIC :: SolveGravity_Newtonian_PointMass

CONTAINS


  SUBROUTINE InitializeGravitySolver_Newtonian_PointMass &
               ( PointMass_Option )

    REAL(DP), INTENT(in), OPTIONAL :: PointMass_Option

    PointMass = 1.0_DP
    IF( PRESENT( PointMass_Option ) ) &
      PointMass = PointMass_Option

    WRITE(*,*)
    WRITE(*,'(A7,A12,ES10.4E2,A1,A)') &
      '', 'PointMass = ', PointMass / UnitsDisplay % MassUnit, &
      '', UnitsDisplay % MassLabel

  END SUBROUTINE InitializeGravitySolver_Newtonian_PointMass


  SUBROUTINE FinalizeGravitySolver_Newtonian_PointMass

    PointMass = 1.0_DP

  END SUBROUTINE FinalizeGravitySolver_Newtonian_PointMass


  SUBROUTINE SolveGravity_Newtonian_PointMass

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNode
    REAL(DP) :: X1

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNode = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uGF(iNode,iX1,iX2,iX3,iGF_Phi_N) &
                  = - GravitationalConstant * PointMass / X1

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

    CALL SetBoundaryConditions

  END SUBROUTINE SolveGravity_Newtonian_PointMass


  SUBROUTINE SetBoundaryConditions

    CALL SetBoundaryConditions_X1

  END SUBROUTINE SetBoundaryConditions


  SUBROUTINE SetBoundaryConditions_X1

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: X1

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

              ! --- Inner Boundary: Dirichlet ---

              X1 = NodeCoordinate( MeshX(1), 0, iNodeX1 )

              uGF(iNodeX,0,iX2,iX3,iGF_Phi_N) &
                = - PointMass / X1

              ! --- Outer Boundary: Dirichlet ---

              X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

              uGF(iNodeX,nX(1)+1,iX2,iX3,iGF_Phi_N) &
                = - PointMass / X1

            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE SetBoundaryConditions_X1


END MODULE GravitySolutionModule_Newtonian_PointMass_Old
