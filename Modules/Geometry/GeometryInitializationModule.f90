MODULE GeometryInitializationModule

  USE KindModule, ONLY: &
    DP
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    InitializeGeometryFields_CARTESIAN, &
    InitializeGeometryFields_SPHERICAL, &
    InitializeGeometryFields_CYLINDRICAL, &
    FinalizeGeometryFields
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX, &
    ComputeGeometry

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGeometry
  PUBLIC :: FinalizeGeometry

CONTAINS


  SUBROUTINE InitializeGeometry &
               ( nX, nNodesX, swX, nE, nNodesE, swE, CoordinateSystem_Option )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, nNodesX, swX
    INTEGER,               INTENT(in) :: nE, nNodesE, swE
    CHARACTER(LEN=*),      INTENT(in), OPTIONAL :: &
      CoordinateSystem_Option

    IF( PRESENT( CoordinateSystem_Option ) )THEN
      CoordinateSystem = TRIM( CoordinateSystem_Option )
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A19,A)') &
      '', 'Coordinate System: ', TRIM( CoordinateSystem )
    WRITE(*,'(A5,A19)') &
      '', '------------------ '

    SELECT CASE ( TRIM( CoordinateSystem ) )
      CASE ( 'CARTESIAN' )

        CALL InitializeGeometryFields_CARTESIAN

      CASE ( 'SPHERICAL' )

        CALL InitializeGeometryFields_SPHERICAL

      CASE ( 'CYLINDRICAL' )

        CALL InitializeGeometryFields_CYLINDRICAL

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

    CALL ComputeGeometryX &
           ( nX, nNodesX, swX )

    CALL ComputeGeometry &
           ( nX, nNodesX, swX, nE, nNodesE, swE )

  END SUBROUTINE InitializeGeometry


  SUBROUTINE FinalizeGeometry

    CALL FinalizeGeometryFields

  END SUBROUTINE FinalizeGeometry


END MODULE GeometryInitializationModule
