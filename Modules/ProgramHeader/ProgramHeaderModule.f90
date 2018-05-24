MODULE ProgramHeaderModule

  USE KindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=32), PUBLIC :: ProgramName
  INTEGER,           PUBLIC :: nNodes

  ! --- Position Space ---

  INTEGER,  PUBLIC :: nX(3)
  INTEGER,  PUBLIC :: swX(3)
  INTEGER,  PUBLIC :: bcX(3)
  INTEGER,  PUBLIC :: iX_B0(3)
  INTEGER,  PUBLIC :: iX_B1(3)
  INTEGER,  PUBLIC :: iX_E0(3)
  INTEGER,  PUBLIC :: iX_E1(3)
  REAL(DP), PUBLIC :: xL(3)
  REAL(DP), PUBLIC :: xR(3)
  REAL(DP), PUBLIC :: zoomX(3)

  INTEGER,  PUBLIC :: nNodesX(3)
  INTEGER,  PUBLIC :: nDimsX
  INTEGER,  PUBLIC :: nDOFX

  ! --- Energy Space ---

  INTEGER,  PUBLIC :: nE
  INTEGER,  PUBLIC :: swE
  INTEGER,  PUBLIC :: bcE
  INTEGER,  PUBLIC :: iE_B0
  INTEGER,  PUBLIC :: iE_B1
  INTEGER,  PUBLIC :: iE_E0
  INTEGER,  PUBLIC :: iE_E1
  REAL(DP), PUBLIC :: eL
  REAL(DP), PUBLIC :: eR
  REAL(DP), PUBLIC :: ZoomE

  INTEGER,  PUBLIC :: nNodesE
  INTEGER,  PUBLIC :: nDimsE
  INTEGER,  PUBLIC :: nDOFE

  ! --- Energy-Position Space ---

  INTEGER,  PUBLIC :: nZ(4)
  INTEGER,  PUBLIC :: swZ(4)
  INTEGER,  PUBLIC :: bcZ(4)
  INTEGER,  PUBLIC :: iZ_B0(4)
  INTEGER,  PUBLIC :: iZ_B1(4)
  INTEGER,  PUBLIC :: iZ_E0(4)
  INTEGER,  PUBLIC :: iZ_E1(4)
  REAL(DP), PUBLIC :: zL(4)
  REAL(DP), PUBLIC :: zR(4)
  REAL(DP), PUBLIC :: zoomZ(4)

  INTEGER,  PUBLIC :: nNodesZ(4)
  INTEGER,  PUBLIC :: nDims
  INTEGER,  PUBLIC :: nDOF

  REAL(DP), PUBLIC, PARAMETER :: NewtonTol = 1.0d-12

  PUBLIC :: InitializeProgramHeader

CONTAINS


  SUBROUTINE InitializeProgramHeader &
    ( ProgramName_Option, nNodes_Option, nX_Option, swX_Option,  &
      bcX_Option, xL_Option, xR_Option, zoomX_Option, nE_Option, &
      swE_Option, bcE_Option, eL_Option, eR_Option, zoomE_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ProgramName_Option
    INTEGER,          INTENT(in), OPTIONAL :: nNodes_Option
    INTEGER,          INTENT(in), OPTIONAL :: nX_Option(3)
    INTEGER,          INTENT(in), OPTIONAL :: swX_Option(3)
    INTEGER,          INTENT(in), OPTIONAL :: bcX_Option(3)
    REAL(DP),         INTENT(in), OPTIONAL :: xL_Option(3)
    REAL(DP),         INTENT(in), OPTIONAL :: xR_Option(3)
    REAL(DP),         INTENT(in), OPTIONAL :: zoomX_Option(3)
    INTEGER,          INTENT(in), OPTIONAL :: nE_Option
    INTEGER,          INTENT(in), OPTIONAL :: swE_Option
    INTEGER,          INTENT(in), OPTIONAL :: bcE_Option
    REAL(DP),         INTENT(in), OPTIONAL :: eL_Option
    REAL(DP),         INTENT(in), OPTIONAL :: eR_Option
    REAL(DP),         INTENT(in), OPTIONAL :: zoomE_Option

    INTEGER :: iDim

    IF( PRESENT( ProgramName_Option ) )THEN
      ProgramName = TRIM( ProgramName_Option )
    ELSE
      ProgramName = 'thornado'
    END IF

    IF( LEN_TRIM( ProgramName ) > 0 )THEN

      WRITE(*,*)
      WRITE(*,'(A2,A28,A)') &
        '', 'INFO: Initializing Program: ', TRIM( ProgramName )
      WRITE(*,*)

    END IF

    IF( PRESENT( nNodes_Option ) )THEN
      nNodes = nNodes_Option
    ELSE
      nNodes = 1
    END IF

    ! --- Position Space ---

    IF( PRESENT( nX_Option ) )THEN
      nX = nX_Option
    ELSE
      nX = [ 1, 1, 1 ]
    END IF

    IF( PRESENT( swX_Option ) )THEN
      swX = swX_Option
    ELSE
      swX = [ 1, 0, 0 ]
    END IF

    IF( PRESENT( bcX_Option ) )THEN
      bcX = bcX_Option
    ELSE
      bcX = [ 0, 0, 0 ]
    END IF

    iX_B0 = 1;  iX_B1 = 1  - swX
    iX_E0 = nX; iX_E1 = nX + swX

    IF( PRESENT( xL_Option ) )THEN
      xL = xL_Option
    ELSE
      xL = [ 0.0d0, 0.0d0, 0.0d0 ]
    END IF

    IF( PRESENT( xR_Option ) )THEN
      xR = xR_Option
    ELSE
      xR = [ 1.0d0, 1.0d0, 1.0d0 ]
    END IF

    IF( PRESENT( zoomX_Option ) )THEN
      zoomX = zoomX_Option
    ELSE
      zoomX = [ 1.0d0, 1.0d0, 1.0d0 ]
    END IF

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

    IF( PRESENT( nE_Option ) )THEN
      nE = nE_Option
    ELSE
      nE = 1
    END IF

    IF( PRESENT( swE_Option ) )THEN
      swE = swE_Option
    ELSE
      swE = 0
    END IF

    IF( PRESENT( bcE_Option ) )THEN
      bcE = bcE_Option
    ELSE
      bcE = 0
    END IF

    iE_B0 = 1;  iE_B1 = 1  - swE
    iE_E0 = nE; iE_E1 = nE + swE

    IF( PRESENT( eL_Option ) )THEN
      eL = eL_Option
    ELSE
      eL = 0.0d0
    END IF

    IF( PRESENT( eR_Option ) )THEN
      eR = eR_Option
    ELSE
      eR = 1.0d0
    END IF

    IF( PRESENT( zoomE_Option ) )THEN
      zoomE = zoomE_Option
    ELSE
      zoomE = 1.0d0
    END IF

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
