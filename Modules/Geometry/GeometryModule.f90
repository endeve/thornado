MODULE GeometryModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(16) :: &
    CoordinateSystem &
      = 'CARTESIAN'

  PUBLIC :: InitializeGeometry
  PUBLIC :: FinalizeGeometry

CONTAINS


  SUBROUTINE InitializeGeometry( CoordinateSystem_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: CoordinateSystem_Option

    IF( PRESENT( CoordinateSystem_Option ) )THEN
      CoordinateSystem = TRIM( CoordinateSystem_Option )
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A19,A)') &
      '', 'Coordinate System: ', TRIM( CoordinateSystem )
    WRITE(*,'(A5,A19)') &
      '', '------------------ '

  END SUBROUTINE InitializeGeometry


  SUBROUTINE FinalizeGeometry

  END SUBROUTINE FinalizeGeometry


END MODULE GeometryModule
