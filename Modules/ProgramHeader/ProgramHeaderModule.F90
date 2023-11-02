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
  INTEGER,  PUBLIC :: nDimsZ
  INTEGER,  PUBLIC :: nDOFZ

  INTEGER,  PUBLIC :: nDims
  INTEGER,  PUBLIC :: nDOF

  REAL(DP), PUBLIC, PARAMETER :: NewtonTol = 1.0d-12

  LOGICAL :: Verbose

  PUBLIC :: InitializeProgramHeader
  PUBLIC :: InitializeProgramHeaderX
  PUBLIC :: InitializeProgramHeaderE
  PUBLIC :: DescribeProgramHeader
  PUBLIC :: DescribeProgramHeaderX

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE &
  !$OMP TARGET( nE, nNodesE, nDimsE, nDOFE, swE, bcE, eL, eR, ZoomE, &
  !$OMP         nX, nNodesX, nDimsX, nDOFX, swX, bcX, xL, xR, zoomX, &
  !$OMP         nZ, nNodesZ, nDimsZ, nDOFZ, swZ, bcZ, zL, zR, zoomZ, &
  !$OMP         iE_B0, iE_E0, iE_B1, iE_E1, &
  !$OMP         iX_B0, iX_E0, iX_B1, iX_E1, &
  !$OMP         iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE &
  !$ACC CREATE( nE, nNodesE, nDimsE, nDOFE, swE, bcE, eL, eR, ZoomE, &
  !$ACC         nX, nNodesX, nDimsX, nDOFX, swX, bcX, xL, xR, zoomX, &
  !$ACC         nZ, nNodesZ, nDimsZ, nDOFZ, swZ, bcZ, zL, zR, zoomZ, &
  !$ACC         iE_B0, iE_E0, iE_B1, iE_E1, &
  !$ACC         iX_B0, iX_E0, iX_B1, iX_E1, &
  !$ACC         iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

CONTAINS


  SUBROUTINE InitializeProgramHeader &
    ( ProgramName_Option, nNodes_Option, nX_Option, swX_Option,  &
      bcX_Option, xL_Option, xR_Option, zoomX_Option, nE_Option, &
      swE_Option, bcE_Option, eL_Option, eR_Option, zoomE_Option, &
      Verbose_Option )

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
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    IF( PRESENT( ProgramName_Option ) )THEN
      ProgramName = TRIM( ProgramName_Option )
    ELSE
      ProgramName = 'thornado'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN

#if defined( THORNADO_GIT_HASH )

      WRITE(*,*)

      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git hash:   ', THORNADO_GIT_HASH
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git date:   ', THORNADO_GIT_DATE
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git branch: ', THORNADO_GIT_BRANCH
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git url:    ', THORNADO_GIT_URL

#endif

    END IF

    IF( ( LEN_TRIM( ProgramName ) > 0 ) .AND. Verbose )THEN

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

    CALL InitializeProgramHeaderX &
           ( nX_Option = nX_Option, swX_Option = swX_Option, &
             bcX_Option = bcX_Option, xL_Option = xL_Option, &
             xR_Option = xR_Option, zoomX_Option = zoomX_Option )

    ! --- Energy Space ---

    CALL InitializeProgramHeaderE &
           ( nE_Option = nE_Option, swE_Option = swE_Option, &
             bcE_Option = bcE_Option, eL_Option = eL_Option, &
             eR_Option = eR_Option, zoomE_Option = zoomE_Option )

    ! --- Energy-Position Space ---

    CALL InitializeProgramHeaderZ


  END SUBROUTINE InitializeProgramHeader


  SUBROUTINE InitializeProgramHeaderX &
    ( nX_Option, swX_Option, bcX_Option, xL_Option, xR_Option, zoomX_Option, &
      reinitializeZ_Option )

    INTEGER,  INTENT(in), OPTIONAL :: nX_Option(3)
    INTEGER,  INTENT(in), OPTIONAL :: swX_Option(3)
    INTEGER,  INTENT(in), OPTIONAL :: bcX_Option(3)
    REAL(DP), INTENT(in), OPTIONAL :: xL_Option(3)
    REAL(DP), INTENT(in), OPTIONAL :: xR_Option(3)
    REAL(DP), INTENT(in), OPTIONAL :: zoomX_Option(3)
    LOGICAL,  INTENT(in), OPTIONAL :: reinitializeZ_Option

    LOGICAL :: reinitializeZ
    INTEGER :: iDim

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

    IF( PRESENT( reinitializeZ_Option ) )THEN
      reinitializeZ = reinitializeZ_Option
    ELSE
      reinitializeZ = .FALSE.
    END IF

    IF( reinitializeZ )THEN
      CALL InitializeProgramHeaderZ
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE &
    !$OMP TO( nX, nNodesX, nDimsX, nDOFX, swX, bcX, xL, xR, zoomX, &
    !$OMP     iX_B0, iX_E0, iX_B1, iX_E1 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( nX, nNodesX, nDimsX, nDOFX, swX, bcX, xL, xR, zoomX, &
    !$ACC         iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

  END SUBROUTINE InitializeProgramHeaderX


  SUBROUTINE InitializeProgramHeaderE &
    ( nE_Option, swE_Option, bcE_Option, eL_Option, eR_Option, zoomE_Option, &
      reinitializeZ_Option )

    INTEGER,  INTENT(in), OPTIONAL :: nE_Option
    INTEGER,  INTENT(in), OPTIONAL :: swE_Option
    INTEGER,  INTENT(in), OPTIONAL :: bcE_Option
    REAL(DP), INTENT(in), OPTIONAL :: eL_Option
    REAL(DP), INTENT(in), OPTIONAL :: eR_Option
    REAL(DP), INTENT(in), OPTIONAL :: zoomE_Option
    LOGICAL,  INTENT(in), OPTIONAL :: reinitializeZ_Option

    LOGICAL :: reinitializeZ

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

    IF( PRESENT( reinitializeZ_Option ) )THEN
      reinitializeZ = reinitializeZ_Option
    ELSE
      reinitializeZ = .FALSE.
    END IF

    IF( reinitializeZ )THEN
      CALL InitializeProgramHeaderZ
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE &
    !$OMP TO( nE, nNodesE, nDimsE, nDOFE, swE, bcE, eL, eR, ZoomE, &
    !$OMP     iE_B0, iE_E0, iE_B1, iE_E1 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( nE, nNodesE, nDimsE, nDOFE, swE, bcE, eL, eR, ZoomE, &
    !$ACC         iE_B0, iE_E0, iE_B1, iE_E1 )
#endif

  END SUBROUTINE InitializeProgramHeaderE


  SUBROUTINE InitializeProgramHeaderZ

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

    nDimsZ = nDimsX + nDimsE
    nDOFZ  = nNodes**nDimsZ

    nDims = nDimsX + nDimsE
    nDOF  = nNodes**nDims

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE &
    !$OMP TO( nZ, nNodesZ, nDimsZ, nDOFZ, swZ, bcZ, zL, zR, zoomZ, &
    !$OMP     iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( nZ, nNodesZ, nDimsZ, nDOFZ, swZ, bcZ, zL, zR, zoomZ, &
    !$ACC         iZ_B0, iZ_E0, iZ_B1, iZ_E1 )
#endif

  END SUBROUTINE InitializeProgramHeaderZ


  SUBROUTINE DescribeProgramHeader

    CALL DescribeProgramHeaderX
    CALL DescribeProgramHeaderE
    CALL DescribeProgramHeaderZ

  END SUBROUTINE DescribeProgramHeader


  SUBROUTINE DescribeProgramHeaderX

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Program Header (X)'
    WRITE(*,*)
    WRITE(*,'(A6,A10,3I6.4)')     '',      'nX = ', nX
    WRITE(*,'(A6,A10,3I6.4)')     '',     'swX = ', swX
    WRITE(*,'(A6,A10,3I6.4)')     '',     'bcX = ', bcX
    WRITE(*,'(A6,A10,3I6.4)')     '',   'iX_B0 = ', iX_B0
    WRITE(*,'(A6,A10,3I6.4)')     '',   'iX_B1 = ', iX_B1
    WRITE(*,'(A6,A10,3I6.4)')     '',   'iX_E0 = ', iX_E0
    WRITE(*,'(A6,A10,3I6.4)')     '',   'iX_E1 = ', iX_E1
    WRITE(*,'(A6,A10,3ES11.2E3)') '',      'xL = ', xL
    WRITE(*,'(A6,A10,3ES11.2E3)') '',      'xR = ', xR
    WRITE(*,'(A6,A10,3ES11.2E3)') '',   'zoomX = ', zoomX
    WRITE(*,'(A6,A10,3I6.4)')     '', 'nNodesX = ', nNodesX
    WRITE(*,'(A6,A10,1I6.4)')     '',  'nDimsX = ', nDimsX
    WRITE(*,'(A6,A10,1I6.4)')     '',   'nDOFX = ', nDOFX
    WRITE(*,*)

  END SUBROUTINE DescribeProgramHeaderX


  SUBROUTINE DescribeProgramHeaderE

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Program Header (E)'
    WRITE(*,*)
    WRITE(*,'(A6,A10,1I6.4)')     '',      'nE = ', nE
    WRITE(*,'(A6,A10,1I6.4)')     '',     'swE = ', swE
    WRITE(*,'(A6,A10,1I6.4)')     '',     'bcE = ', bcE
    WRITE(*,'(A6,A10,1I6.4)')     '',   'iE_B0 = ', iE_B0
    WRITE(*,'(A6,A10,1I6.4)')     '',   'iE_B1 = ', iE_B1
    WRITE(*,'(A6,A10,1I6.4)')     '',   'iE_E0 = ', iE_E0
    WRITE(*,'(A6,A10,1I6.4)')     '',   'iE_E1 = ', iE_E1
    WRITE(*,'(A6,A10,1ES11.2E3)') '',      'eL = ', eL
    WRITE(*,'(A6,A10,1ES11.2E3)') '',      'eR = ', eR
    WRITE(*,'(A6,A10,1ES11.2E3)') '',   'zoomE = ', zoomE
    WRITE(*,'(A6,A10,1I6.4)')     '', 'nNodesE = ', nNodesE
    WRITE(*,'(A6,A10,1I6.4)')     '',  'nDimsE = ', nDimsE
    WRITE(*,'(A6,A10,1I6.4)')     '',   'nDOFE = ', nDOFE
    WRITE(*,*)

  END SUBROUTINE DescribeProgramHeaderE


  SUBROUTINE DescribeProgramHeaderZ

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Program Header (Z)'
    WRITE(*,*)
    WRITE(*,'(A6,A10,4I6.4)')     '',      'nZ = ', nZ
    WRITE(*,'(A6,A10,4I6.4)')     '',     'swZ = ', swZ
    WRITE(*,'(A6,A10,4I6.4)')     '',     'bcZ = ', bcZ
    WRITE(*,'(A6,A10,4I6.4)')     '',   'iZ_B0 = ', iZ_B0
    WRITE(*,'(A6,A10,4I6.4)')     '',   'iZ_B1 = ', iZ_B1
    WRITE(*,'(A6,A10,4I6.4)')     '',   'iZ_E0 = ', iZ_E0
    WRITE(*,'(A6,A10,4I6.4)')     '',   'iZ_E1 = ', iZ_E1
    WRITE(*,'(A6,A10,4ES11.2E3)') '',      'zL = ', zL
    WRITE(*,'(A6,A10,4ES11.2E3)') '',      'zR = ', zR
    WRITE(*,'(A6,A10,4ES11.2E3)') '',   'zoomZ = ', zoomZ
    WRITE(*,'(A6,A10,4I6.4)')     '', 'nNodesZ = ', nNodesZ
    WRITE(*,'(A6,A10,1I6.4)')     '',   'nDims = ', nDims
    WRITE(*,'(A6,A10,1I6.4)')     '',    'nDOF = ', nDOF
    WRITE(*,*)

  END SUBROUTINE DescribeProgramHeaderZ


END MODULE ProgramHeaderModule
