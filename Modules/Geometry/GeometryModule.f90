MODULE GeometryModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(16), PUBLIC :: &
    CoordinateSystem &
      = 'CARTESIAN'

  INTERFACE
    PURE REAL(DP) FUNCTION MetricFunction( X )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(3), INTENT(in) :: X
    END FUNCTION MetricFunction
  END INTERFACE

  PROCEDURE (MetricFunction), POINTER, PUBLIC :: a, b, c, d
  PROCEDURE (MetricFunction), POINTER, PUBLIC :: dlnadX1
  PROCEDURE (MetricFunction), POINTER, PUBLIC :: dlnbdX1
  PROCEDURE (MetricFunction), POINTER, PUBLIC :: dlncdX2


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

    SELECT CASE ( TRIM( CoordinateSystem ) )
      CASE ( 'CARTESIAN' )

        a => a_CARTESIAN
        b => b_CARTESIAN
        c => c_CARTESIAN
        d => SqrtDet_CARTESIAN

        dlnadX1 => dlnadX1_CARTESIAN
        dlnbdX1 => dlnbdX1_CARTESIAN
        dlncdX2 => dlncdX2_CARTESIAN

      CASE ( 'SPHERICAL' )

        a => a_SPHERICAL
        b => b_SPHERICAL
        c => c_SPHERICAL
        d => SqrtDet_SPHERICAL

        dlnadX1 => dlnadX1_SPHERICAL
        dlnbdX1 => dlnbdX1_SPHERICAL
        dlncdX2 => dlncdX2_SPHERICAL

      CASE ( 'CYLINDRICAL' )

        a => a_CYLINDRICAL
        b => b_CYLINDRICAL
        c => c_CYLINDRICAL
        d => SqrtDet_CYLINDRICAL

        dlnadX1 => dlnadX1_CYLINDRICAL
        dlnbdX1 => dlnbdX1_CYLINDRICAL
        dlncdX2 => dlncdX2_CYLINDRICAL

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A27,A)') &
          '', 'Invalid Coordinate System: ', TRIM( CoordinateSystem )
        STOP

    END SELECT

  END SUBROUTINE InitializeGeometry


  SUBROUTINE FinalizeGeometry

    NULLIFY( a, b, c, d )

  END SUBROUTINE FinalizeGeometry


  ! --- Cartesian Coordinates ---


  PURE REAL(DP) FUNCTION a_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    a_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION a_CARTESIAN


  PURE REAL(DP) FUNCTION b_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    b_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION b_CARTESIAN


  PURE REAL(DP) FUNCTION c_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    c_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION c_CARTESIAN


  PURE REAL(DP) FUNCTION SqrtDet_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    SqrtDet_CARTESIAN = 1.0_DP

    RETURN
  END FUNCTION SqrtDet_CARTESIAN


  PURE REAL(DP) FUNCTION dlnadX1_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnadX1_CARTESIAN = 0.0_DP

    RETURN
  END FUNCTION dlnadX1_CARTESIAN


  PURE REAL(DP) FUNCTION dlnbdX1_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnbdX1_CARTESIAN = 0.0_DP

    RETURN
  END FUNCTION dlnbdX1_CARTESIAN


  PURE REAL(DP) FUNCTION dlncdX2_CARTESIAN( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlncdX2_CARTESIAN = 0.0_DP

    RETURN
  END FUNCTION dlncdX2_CARTESIAN


  ! --- Spherical Coordinates ---


  PURE REAL(DP) FUNCTION a_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    a_SPHERICAL = X(1)

    RETURN
  END FUNCTION a_SPHERICAL


  PURE REAL(DP) FUNCTION b_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    b_SPHERICAL = X(1)

    RETURN
  END FUNCTION b_SPHERICAL


  PURE REAL(DP) FUNCTION c_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    c_SPHERICAL = SIN( X(2) )

    RETURN
  END FUNCTION c_SPHERICAL


  PURE REAL(DP) FUNCTION SqrtDet_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    SqrtDet_SPHERICAL = X(1)**2 * SIN( X(2) )

    RETURN
  END FUNCTION SqrtDet_SPHERICAL


  PURE REAL(DP) FUNCTION dlnadX1_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnadX1_SPHERICAL = 1.0_DP / X(1)

    RETURN
  END FUNCTION dlnadX1_SPHERICAL


  PURE REAL(DP) FUNCTION dlnbdX1_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnbdX1_SPHERICAL = 1.0_DP / X(1)

    RETURN
  END FUNCTION dlnbdX1_SPHERICAL


  PURE REAL(DP) FUNCTION dlncdX2_SPHERICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlncdX2_SPHERICAL = 1.0_DP / TAN( X(2) )

    RETURN
  END FUNCTION dlncdX2_SPHERICAL


  ! --- Cylindrical Coordinates ---


  PURE REAL(DP) FUNCTION a_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    a_CYLINDRICAL = 1.0_DP

    RETURN
  END FUNCTION a_CYLINDRICAL


  PURE REAL(DP) FUNCTION b_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    b_CYLINDRICAL = X(1)

    RETURN
  END FUNCTION b_CYLINDRICAL


  PURE REAL(DP) FUNCTION c_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    c_CYLINDRICAL = 1.0_DP

    RETURN
  END FUNCTION c_CYLINDRICAL


  PURE REAL(DP) FUNCTION SqrtDet_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    SqrtDet_CYLINDRICAL = X(1)

    RETURN
  END FUNCTION SqrtDet_CYLINDRICAL


  PURE REAL(DP) FUNCTION dlnadX1_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnadX1_CYLINDRICAL = 0.0_DP

    RETURN
  END FUNCTION dlnadX1_CYLINDRICAL


  PURE REAL(DP) FUNCTION dlnbdX1_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlnbdX1_CYLINDRICAL = 1.0_DP / X(1)

    RETURN
  END FUNCTION dlnbdX1_CYLINDRICAL


  PURE REAL(DP) FUNCTION dlncdX2_CYLINDRICAL( X )

    REAL(DP), DIMENSION(3), INTENT(in) :: X

    dlncdX2_CYLINDRICAL = 0.0_DP

    RETURN
  END FUNCTION dlncdX2_CYLINDRICAL


END MODULE GeometryModule
