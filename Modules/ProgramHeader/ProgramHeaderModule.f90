MODULE ProgramHeaderModule

  USE KindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=32), PUBLIC :: ProgramName
  INTEGER,           PUBLIC :: nNodes

  ! --- Position Space ---

  INTEGER,  DIMENSION(3), PUBLIC :: nX
  INTEGER,  DIMENSION(3), PUBLIC :: nNodesX
  INTEGER,  DIMENSION(3), PUBLIC :: swX
  INTEGER,  DIMENSION(3), PUBLIC :: bcX
  INTEGER,  DIMENSION(3), PUBLIC :: iX_B0
  INTEGER,  DIMENSION(3), PUBLIC :: iX_B1
  INTEGER,  DIMENSION(3), PUBLIC :: iX_E0
  INTEGER,  DIMENSION(3), PUBLIC :: iX_E1
  REAL(DP), DIMENSION(3), PUBLIC :: xL
  REAL(DP), DIMENSION(3), PUBLIC :: xR
  REAL(DP), DIMENSION(3), PUBLIC :: zoomX

  INTEGER,                PUBLIC :: nDimsX
  INTEGER,                PUBLIC :: nDOFX

  ! --- Energy Space ---

  INTEGER,                PUBLIC :: nE
  INTEGER,                PUBLIC :: nNodesE
  INTEGER,                PUBLIC :: swE
  INTEGER,                PUBLIC :: bcE
  INTEGER,                PUBLIC :: iE_B0
  INTEGER,                PUBLIC :: iE_B1
  INTEGER,                PUBLIC :: iE_E0
  INTEGER,                PUBLIC :: iE_E1
  REAL(DP),               PUBLIC :: eL
  REAL(DP),               PUBLIC :: eR
  REAL(DP),               PUBLIC :: ZoomE

  INTEGER,                PUBLIC :: nDimsE
  INTEGER,                PUBLIC :: nDOFE

  ! --- Energy-Position Space ---

  INTEGER,  DIMENSION(4), PUBLIC :: nZ
  INTEGER,  DIMENSION(4), PUBLIC :: nNodesZ
  INTEGER,  DIMENSION(4), PUBLIC :: swZ
  INTEGER,  DIMENSION(4), PUBLIC :: bcZ
  INTEGER,  DIMENSION(4), PUBLIC :: iZ_B0
  INTEGER,  DIMENSION(4), PUBLIC :: iZ_B1
  INTEGER,  DIMENSION(4), PUBLIC :: iZ_E0
  INTEGER,  DIMENSION(4), PUBLIC :: iZ_E1
  REAL(DP), DIMENSION(4), PUBLIC :: zL
  REAL(DP), DIMENSION(4), PUBLIC :: zR
  REAL(DP), DIMENSION(4), PUBLIC :: zoomZ

  INTEGER,                PUBLIC :: nDims
  INTEGER,                PUBLIC :: nDOF

  REAL(DP), PUBLIC, PARAMETER :: NewtonTol = 1.0d-12

  PUBLIC :: InitializeProgramHeader

CONTAINS


  SUBROUTINE InitializeProgramHeader &
    ( ProgramName_Option, nNodes_Option, nX_Option, swX_Option,  &
      bcX_Option, xL_Option, xR_Option, zoomX_Option, nE_Option, &
      swE_Option, bcE_Option, eL_Option, eR_Option, zoomE_Option )

    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: ProgramName_Option
    INTEGER,                INTENT(in), OPTIONAL :: nNodes_Option
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

    INTEGER :: iDim

    ProgramName = 'thornado'
    IF( PRESENT( ProgramName_Option ) ) &
      ProgramName = TRIM( ProgramName_Option )

    nNodes = 1
    IF( PRESENT( nNodes_Option ) ) &
      nNodes = nNodes_Option

    ! --- Position Space ---

    nX = [ 1, 1, 1 ]
    IF( PRESENT( nX_Option ) ) &
      nX = nX_Option

    swX = [ 1, 0, 0 ]
    IF( PRESENT( swX_Option ) ) &
      swX = swX_Option

    bcX = [ 0, 0, 0 ]
    IF( PRESENT( bcX_Option ) )THEN
      bcX = bcX_Option
    END IF

    iX_B0 = 1;  iX_B1 = 1  - swX
    iX_E0 = nX; iX_E1 = nX + swX

    xL = [ 0.0d0, 0.0d0, 0.0d0 ]
    IF( PRESENT( xL_Option ) ) &
      xL = xL_Option

    xR = [ 1.0d0, 1.0d0, 1.0d0 ]
    IF( PRESENT( xR_Option ) ) &
      xR = xR_Option

    zoomX = [ 1.0d0, 1.0d0, 1.0d0 ]
    IF( PRESENT( zoomX_Option ) ) &
      zoomX = zoomX_Option

    nDimsX = 0
    DO iDim = 1, 3
      nNodesX(iDim) = 1
      IF( nX(iDim) > 1 )THEN
        nDimsX = nDimsX + 1
        nNodesX(iDim) = nNodes
      END IF
    END DO

    nDOFX = nNodes**nDimsX

    ! --- Energy Space ---

    nE = 1
    IF( PRESENT( nE_Option ) ) &
      nE = nE_Option

    swE = 0
    IF( PRESENT( swE_Option ) ) &
      swE = swE_Option

    bcE = 0
    IF( PRESENT( bcE_Option ) ) &
      bcE = bcE_Option

    iE_B0 = 1;  iE_B1 = 1  - swE
    iE_E0 = nE; iE_E1 = nE + swE

    eL = 0.0d0
    IF( PRESENT( eL_Option ) ) &
      eL = eL_Option

    eR = 1.0d0
    IF( PRESENT( eR_Option ) ) &
      eR = eR_Option

    zoomE = 1.0d0
    IF( PRESENT( zoomE_Option ) ) &
      zoomE = zoomE_Option

    nDimsE  = 0
    nNodesE = 1
    IF( nE > 1 )THEN
      nDimsE  = 1
      nNodesE = nNodes
    END IF

    nDOFE = nNodes**nDimsE

    ! --- Energy-Position Space ---

    nZ      = [ nE,      nX      ]
    nNodesZ = [ nNodesE, nNodesX ]
    swZ     = [ swE,     swX     ]
    bcZ     = [ bcE,     bcX     ]
    iZ_B0   = [ iE_B0,   iX_B0   ]
    iZ_B1   = [ iE_B1,   iX_B1   ]
    iZ_E0   = [ iE_E0,   iX_E0   ]
    iZ_E1   = [ iE_E1,   iX_E1   ]
    zL      = [ eL,      xL      ]
    zR      = [ eR,      xR      ]
    zoomZ   = [ zoomE,   zoomX   ]

    nDims = nDimsX + nDimsE
    nDOF  = nNodes**nDims

  END SUBROUTINE InitializeProgramHeader


END MODULE ProgramHeaderModule
