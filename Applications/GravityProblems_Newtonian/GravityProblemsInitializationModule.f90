MODULE GravityProblemsInitializationModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeHomogeneousSphere

CONTAINS


  SUBROUTINE InitializeHomogeneousSphere( SphereRadius_Option )

    REAL(DP), INTENT(in), OPTIONAL :: SphereRadius_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX
    REAL(DP) :: SphereRadius, X1

    SphereRadius = 1.0_DP
    IF( PRESENT( SphereRadius_Option ) ) &
      SphereRadius = SphereRadius_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A16,ES10.3E2)') &
      '', 'Sphere Radius = ', SphereRadius
    WRITE(*,*)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                IF( X1 <= SphereRadius )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0_DP

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 0.0_DP

                END IF

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeHomogeneousSphere


END MODULE GravityProblemsInitializationModule
