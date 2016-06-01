MODULE ProgramHeaderModule

  USE KindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=32), PUBLIC :: ProgramName = 'thornado'

  CHARACTER(LEN=*), PUBLIC, PARAMETER :: &
    GeometryName = 'Cartesian'

  INTEGER,  DIMENSION(3), PUBLIC :: nX      = [ 8, 1, 1 ]
  INTEGER,  DIMENSION(3), PUBLIC :: nNodesX = [ 1, 1, 1 ]
  INTEGER,  DIMENSION(3), PUBLIC :: swX     = [ 1, 0, 0 ]
  INTEGER,  DIMENSION(3), PUBLIC :: bcX     = [ 1, 0, 0 ]
  REAL(DP), DIMENSION(3), PUBLIC :: xL      = [ 0.0d0, 0.0d0, 0.0d0 ]
  REAL(DP), DIMENSION(3), PUBLIC :: xR      = [ 1.0d0, 1.0d0, 1.0d0 ]
  REAL(DP), DIMENSION(3), PUBLIC :: zoomX   = [ 1.0d0, 1.0d0, 1.0d0 ]

  INTEGER,  PUBLIC :: nE      = 1
  INTEGER,  PUBLIC :: nNodesE = 1
  INTEGER,  PUBLIC :: swE     = 0
  INTEGER,  PUBLIC :: bcE     = 0
  REAL(DP), PUBLIC :: eL      = 0.0d0
  REAL(DP), PUBLIC :: eR      = 1.0d0
  REAL(DP), PUBLIC :: ZoomE   = 1.0d0

  INTEGER,  PUBLIC :: nNodes = 1

  INTEGER, PUBLIC :: nDOFX
  INTEGER, PUBLIC :: nDOFE
  INTEGER, PUBLIC :: nDOF
  INTEGER, PUBLIC :: nDimsX
  INTEGER, PUBLIC :: nDimsE
  INTEGER, PUBLIC :: nDims

  REAL(DP), PUBLIC, PARAMETER :: NewtonTol = 1.0d-12

  PUBLIC :: InitializeProgramHeader

CONTAINS


  SUBROUTINE InitializeProgramHeader &
               ( ProgramName_Option, nX_Option, swX_Option, bcX_Option, &
                 xL_Option, xR_Option, zoomX_Option, nE_Option, swE_Option, &
                 bcE_Option, eL_Option, eR_Option, zoomE_Option, nNodes_Option )

    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: ProgramName_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: nX_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: swX_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: bcX_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: xL_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: xR_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: zoomX_Option
    INTEGER,                INTENT(in), OPTIONAL :: nE_Option
    INTEGER,                INTENT(in), OPTIONAL :: swE_Option
    INTEGER,                INTENT(in), OPTIONAL :: bcE_Option
    REAL(DP),               INTENT(in), OPTIONAL :: eL_Option
    REAL(DP),               INTENT(in), OPTIONAL :: eR_Option
    REAL(DP),               INTENT(in), OPTIONAL :: zoomE_Option
    INTEGER,                INTENT(in), OPTIONAL :: nNodes_Option

    IF( PRESENT( ProgramName_Option ) )THEN
      ProgramName = TRIM( ProgramName_Option )
    END IF

    IF( PRESENT( nX_Option ) )THEN
      nX = nX_Option
    END IF

    IF( PRESENT( swX_Option ) )THEN
      swX = swX_Option
    END IF

    IF( PRESENT( bcX_Option ) )THEN
      bcX = bcX_Option
    END IF

    IF( PRESENT( xL_Option ) )THEN
      xL = xL_Option
    END IF

    IF( PRESENT( xR_Option ) )THEN
      xR = xR_Option
    END IF

    IF( PRESENT( zoomX_Option ) )THEN
      zoomX = zoomX_Option
    END IF

    IF( PRESENT( nE_Option ) )THEN
      nE = nE_Option
    END IF

    IF( PRESENT( swE_Option ) )THEN
      swE = swE_Option
    END IF

    IF( PRESENT( bcE_Option ) )THEN
      bcE = bcE_Option
    END IF

    IF( PRESENT( eL_Option ) )THEN
      eL = eL_Option
    END IF

    IF( PRESENT( eR_Option ) )THEN
      eR = eR_Option
    END IF

    IF( PRESENT( zoomE_Option ) )THEN
      zoomE = zoomE_Option
    END IF

    IF( PRESENT( nNodes_Option ) )THEN
      nNodes = nNodes_Option
    END IF

  END SUBROUTINE InitializeProgramHeader


END MODULE ProgramHeaderModule
