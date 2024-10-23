MODULE TimersModule_MHD

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64
  USE KindModule, Only: &
    DP, SqrtTiny

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PUBLIC :: TimeIt_MHD = .FALSE.

  REAL(DP), PUBLIC :: Timer_MHD_Program

  ! --- ApplicationDriver ---

  REAL(DP), PUBLIC :: Timer_MHD_Initialize
  REAL(DP), PUBLIC :: Timer_MHD_UpdateFluid
  REAL(DP), PUBLIC :: Timer_MHD_InputOutput
  REAL(DP), PUBLIC :: Timer_MHD_Finalize

  REAL(DP), PUBLIC :: Timer_MHD_ComputeTimeStep
  REAL(DP), PUBLIC :: Timer_MHD_CTS_ComputeTimeStep
  REAL(DP), PUBLIC :: Timer_MHD_CTS_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_CTS_CopyOut

  ! --- DG discretization ---

  REAL(DP), PUBLIC :: Timer_MHD_DG
  REAL(DP), PUBLIC :: Timer_MHD_Increment
  REAL(DP), PUBLIC :: Timer_MHD_Divergence
  REAL(DP), PUBLIC :: Timer_MHD_SurfaceTerm
  REAL(DP), PUBLIC :: Timer_MHD_VolumeTerm
  REAL(DP), PUBLIC :: Timer_MHD_Geometry
  REAL(DP), PUBLIC :: Timer_MHD_Gravity
  REAL(DP), PUBLIC :: Timer_MHD_DG_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_DG_CopyOut
  REAL(DP), PUBLIC :: Timer_MHD_DG_ComputePrimitive
  REAL(DP), PUBLIC :: Timer_MHD_DG_Permute
  REAL(DP), PUBLIC :: Timer_MHD_DG_Interpolate
  REAL(DP), PUBLIC :: Timer_MHD_DG_ErrorCheck

  ! --- Compute Primitive --

  REAL(DP), PUBLIC :: Timer_MHD_ComputePrimitive
  REAL(DP), PUBLIC :: Timer_MHD_CP_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_CP_CopyOut
  REAL(DP), PUBLIC :: Timer_MHD_CP_Permute
  REAL(DP), PUBLIC :: Timer_MHD_CP_GetBounds
  REAL(DP), PUBLIC :: Timer_MHD_CP_Bisection
  REAL(DP), PUBLIC :: Timer_MHD_CP_RecoverPrimitives

  ! --- Discontinuity Detection ---

  ! --- Troubled-Cell Indicator

  REAL(DP), PUBLIC :: Timer_MHD_DD_TCI
  REAL(DP), PUBLIC :: Timer_MHD_DD_TCI_DetectTroubledCells
  REAL(DP), PUBLIC :: Timer_MHD_DD_TCI_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_DD_TCI_CopyOut
  REAL(DP), PUBLIC :: Timer_MHD_DD_TCI_Permute
  REAL(DP), PUBLIC :: Timer_MHD_DD_TCI_Integrate

  ! --- Shock Detector ---

  REAL(DP), PUBLIC :: Timer_MHD_DD_SD
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_DetectShocks
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_CopyOut
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_Permute
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_Integrate
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_ComputePrimitive
  REAL(DP), PUBLIC :: Timer_MHD_DD_SD_ErrorCheck

  ! --- Limiter-specific ---

  ! --- Slope-Limiter ---

  REAL(DP), PUBLIC :: Timer_MHD_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_MHD_SL_LimitCells
  REAL(DP), PUBLIC :: Timer_MHD_SL_CharDecomp
  REAL(DP), PUBLIC :: Timer_MHD_SL_ConsCorr
  REAL(DP), PUBLIC :: Timer_MHD_SL_Mapping
  REAL(DP), PUBLIC :: Timer_MHD_SL_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_SL_CopyOut
  REAL(DP), PUBLIC :: Timer_MHD_SL_Permute
  REAL(DP), PUBLIC :: Timer_MHD_SL_Integrate

  ! --- Positivity-Limiter ---

  REAL(DP), PUBLIC :: Timer_MHD_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_MHD_PL_LimitCells
  REAL(DP), PUBLIC :: Timer_MHD_PL_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_PL_CopyOut
  REAL(DP), PUBLIC :: Timer_MHD_PL_Permute
  REAL(DP), PUBLIC :: Timer_MHD_PL_Integrate
  REAL(DP), PUBLIC :: Timer_MHD_PL_ErrorCheck

  ! --- Boundary Conditions ---

  REAL(DP), PUBLIC :: Timer_MHD_BoundaryConditions
  REAL(DP), PUBLIC :: Timer_MHD_BC_ApplyBC
  REAL(DP), PUBLIC :: Timer_MHD_BC_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_BC_CopyOut

  ! --- Compute From Conserved ---

  REAL(DP), PUBLIC :: Timer_MHD_ComputeFromConserved
  REAL(DP), PUBLIC :: Timer_MHD_CFC_ComputePrimitive
  REAL(DP), PUBLIC :: Timer_MHD_CFC_ErrorCheck
  REAL(DP), PUBLIC :: Timer_MHD_CFC_CopyIn
  REAL(DP), PUBLIC :: Timer_MHD_CFC_CopyOut

  ! --- Gravity Solver ---

  REAL(DP), PUBLIC :: Timer_GravitySolver
  REAL(DP), PUBLIC :: Timer_GS_ComputeSourceTerms

  PUBLIC :: InitializeTimers_MHD
  PUBLIC :: FinalizeTimers_MHD
  PUBLIC :: TimersStart_MHD
  PUBLIC :: TimersStop_MHD
  PUBLIC :: TimersWtime_MHD


CONTAINS


  SUBROUTINE InitializeTimers_MHD

    IF( .NOT. TimeIt_MHD ) RETURN

    Timer_MHD_Program = SqrtTiny

    CALL TimersStart_MHD( Timer_MHD_Program )

    Timer_MHD_Initialize  = SqrtTiny
    Timer_MHD_UpdateFluid = SqrtTiny
    Timer_MHD_InputOutput = SqrtTiny
    Timer_MHD_Finalize    = SqrtTiny

    Timer_MHD_ComputeTimeStep     = SqrtTiny
    Timer_MHD_CTS_ComputeTimeStep = SqrtTiny
    Timer_MHD_CTS_CopyIn          = SqrtTiny
    Timer_MHD_CTS_CopyOut         = SqrtTiny

    Timer_MHD_DG                  = SqrtTiny
    Timer_MHD_Increment           = SqrtTiny
    Timer_MHD_Divergence          = SqrtTiny
    Timer_MHD_SurfaceTerm         = SqrtTiny
    Timer_MHD_VolumeTerm          = SqrtTiny
    Timer_MHD_Geometry            = SqrtTiny
    Timer_MHD_Gravity             = SqrtTiny
    Timer_MHD_DG_CopyIn           = SqrtTiny
    Timer_MHD_DG_CopyOut          = SqrtTiny
    Timer_MHD_DG_ComputePrimitive = SqrtTiny
    Timer_MHD_DG_Permute          = SqrtTiny
    Timer_MHD_DG_Interpolate      = SqrtTiny
    Timer_MHD_DG_ErrorCheck       = SqrtTiny

    Timer_MHD_ComputePrimitive     = SqrtTiny
    TImer_MHD_CP_CopyIn            = SqrtTiny
    TImer_MHD_CP_CopyOut           = SqrtTiny
    TImer_MHD_CP_Permute           = SqrtTiny
    TImer_MHD_CP_GetBounds         = SqrtTiny
    TImer_MHD_CP_Bisection         = SqrtTiny
    TImer_MHD_CP_RecoverPrimitives = SqrtTiny

    Timer_MHD_DD_TCI                     = SqrtTiny
    Timer_MHD_DD_TCI_DetectTroubledCells = SqrtTiny
    Timer_MHD_DD_TCI_CopyIn              = SqrtTiny
    Timer_MHD_DD_TCI_CopyOut             = SqrtTiny
    Timer_MHD_DD_TCI_Permute             = SqrtTiny
    Timer_MHD_DD_TCI_Integrate           = SqrtTiny

    Timer_MHD_DD_SD                  = SqrtTiny
    Timer_MHD_DD_SD_DetectShocks     = SqrtTiny
    Timer_MHD_DD_SD_CopyIn           = SqrtTiny
    Timer_MHD_DD_SD_CopyOut          = SqrtTiny
    Timer_MHD_DD_SD_Permute          = SqrtTiny
    Timer_MHD_DD_SD_Integrate        = SqrtTiny
    Timer_MHD_DD_SD_ComputePrimitive = SqrtTiny
    Timer_MHD_DD_SD_ErrorCheck       = SqrtTiny

    Timer_MHD_SlopeLimiter  = SqrtTiny
    Timer_MHD_SL_LimitCells = SqrtTiny
    Timer_MHD_SL_CharDecomp = SqrtTiny
    Timer_MHD_SL_ConsCorr   = SqrtTiny
    Timer_MHD_SL_Mapping    = SqrtTiny
    Timer_MHD_SL_CopyIn     = SqrtTiny
    Timer_MHD_SL_CopyOut    = SqrtTiny
    Timer_MHD_SL_Permute    = SqrtTiny
    Timer_MHD_SL_Integrate  = SqrtTiny

    Timer_MHD_PositivityLimiter = SqrtTiny
    Timer_MHD_PL_LimitCells     = SqrtTiny
    Timer_MHD_PL_CopyIn         = SqrtTiny
    Timer_MHD_PL_CopyOut        = SqrtTiny
    Timer_MHD_PL_Permute        = SqrtTiny
    Timer_MHD_PL_Integrate      = SqrtTiny
    Timer_MHD_PL_ErrorCheck     = SqrtTiny

    Timer_MHD_BoundaryConditions = SqrtTiny
    Timer_MHD_BC_ApplyBC         = SqrtTiny
    Timer_MHD_BC_CopyIn          = SqrtTiny
    Timer_MHD_BC_CopyOut         = SqrtTiny

    Timer_MHD_ComputeFromConserved = SqrtTiny
    Timer_MHD_CFC_ComputePrimitive = SqrtTiny
    Timer_MHD_CFC_ErrorCheck       = SqrtTiny
    Timer_MHD_CFC_CopyIn           = SqrtTiny
    Timer_MHD_CFC_CopyOut          = SqrtTiny

    Timer_GravitySolver         = SqrtTiny
    Timer_GS_ComputeSourceTerms = SqrtTiny

  END SUBROUTINE InitializeTimers_MHD


  SUBROUTINE FinalizeTimers_MHD &
    ( Verbose_Option, SuppressApplicationDriver_Option, &
      WriteAtIntermediateTime_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL, INTENT(in), OPTIONAL :: SuppressApplicationDriver_Option
    LOGICAL, INTENT(in), OPTIONAL :: WriteAtIntermediateTime_Option

    INTEGER  :: iT, nT
    LOGICAL  :: Verbose, SuppressApplicationDriver, WriteAtIntermediateTime
    REAL(DP) :: TotalTime

    REAL(DP),      ALLOCATABLE :: Timers(:)
    CHARACTER(64), ALLOCATABLE :: Labels(:)

    CHARACTER(6)  :: Label = '(6x,A)'

    CHARACTER(64) :: TimeL1 = '(8x,A,ES10.3E3,A,ES10.3E3)'
    CHARACTER(64) :: TimeL2 = '(8x,A,ES10.3E3,A,ES10.3E3,A,ES10.3E3)'

    IF( .NOT. TimeIt_MHD ) RETURN

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    SuppressApplicationDriver = .FALSE.
    IF( PRESENT( SuppressApplicationDriver_Option ) ) &
      SuppressApplicationDriver = SuppressApplicationDriver_Option

    WriteAtIntermediateTime = .FALSE.
    IF( PRESENT( WriteAtIntermediateTime_Option ) ) &
      WriteAtIntermediateTime = WriteAtIntermediateTime_Option

    CALL TimersStop_MHD( Timer_MHD_Program )

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(2x,A)') 'Timers Summary (MHD)'
      WRITE(*,'(2x,A)') '----------------------'
      WRITE(*,*)

      WRITE(*,'(4x,A,ES13.6E3,A)') &
        'Total run-time = ', Timer_MHD_Program, ' s'

    END IF

    IF( .NOT. SuppressApplicationDriver )THEN

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Application Driver'
      WRITE(*,TRIM(Label)) '------------------'
      WRITE(*,*)

      nT = 5
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_Initialize, &
                 Timer_MHD_ComputeTimeStep, &
                 Timer_MHD_UpdateFluid, &
                 Timer_MHD_InputOutput, &
                 Timer_MHD_Finalize ]

      Labels = [ 'Initialize        : ', &
                 'Compute Time-Step : ', &
                 'Update Fluid      : ', &
                 'Input/Output      : ', &
                 'Finalize          : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,'(8x,A,ES13.6E3,A,F7.3,A)') &
        'Timers = ', TotalTime, ' s = ', &
        TotalTime / Timer_MHD_Program

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL1, Labels(iT), &
                 Timers(iT), Timer_MHD_Program )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'DG discretization'
      WRITE(*,TRIM(Label)) '-----------------'
      WRITE(*,*)

      nT = 6
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_Increment, &
                 Timer_MHD_Geometry, &
                 Timer_MHD_Gravity, &
                 Timer_MHD_Divergence, &
                 Timer_MHD_DG_CopyIn, &
                 Timer_MHD_DG_CopyOut ]

      Labels = [ 'Increment   : ', &
                 'Geometry    : ', &
                 'Gravity     : ', &
                 'Divergence  : ', &
                 'CopyIn      : ', &
                 'CopyOut     : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'dgDiscretization : ', &
        Timer_MHD_DG, ' s = ', &
        Timer_MHD_DG / Timer_MHD_Program, ' = ', &
        Timer_MHD_DG / TotalTime

      WRITE(*,*)

      DO iT = 1, nT - 2

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_DG )

      END DO

      WRITE(*,*)

      WRITE(*,TRIM(TimeL2)) &
        '  Surface Term       : ', &
        Timer_MHD_SurfaceTerm, ' s = ', &
        Timer_MHD_SurfaceTerm / Timer_MHD_Program, ' = ', &
        Timer_MHD_SurfaceTerm / Timer_MHD_Divergence

      WRITE(*,TRIM(TimeL2)) &
        '  Volume Term        : ', &
        Timer_MHD_VolumeTerm, ' s = ', &
        Timer_MHD_VolumeTerm / Timer_MHD_Program, ' = ', &
        Timer_MHD_VolumeTerm / Timer_MHD_Divergence

      WRITE(*,TRIM(TimeL2)) &
        '  Compute Primitive  : ', &
        Timer_MHD_DG_ComputePrimitive, ' s = ', &
        Timer_MHD_DG_ComputePrimitive / Timer_MHD_Program, ' = ', &
        Timer_MHD_DG_ComputePrimitive / Timer_MHD_Divergence

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      nT = 5
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_DG_CopyIn, &
                 Timer_MHD_DG_CopyOut, &
                 Timer_MHD_DG_Interpolate, &
                 Timer_MHD_DG_Permute, &
                 Timer_MHD_DG_ErrorCheck ]

      Labels = [ 'CopyIn      : ', &
                 'CopyOut     : ', &
                 'Interpolate : ', &
                 'Permute     : ', &
                 'Error Check : ' ]

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_DG )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Compute Primitive'
      WRITE(*,TRIM(Label)) '-----------------'
      WRITE(*,*)

      nT = 6
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_CP_CopyIn, &
                 Timer_MHD_CP_CopyOut, &
                 Timer_MHD_CP_Permute, &
                 Timer_MHD_CP_GetBounds, &
                 Timer_MHD_CP_Bisection, &
                 Timer_MHD_CP_RecoverPrimitives ]

      Labels = [ 'CopyIn            : ', &
                 'CopyOut           : ', &
                 'Permute           : ', &
                 'GetBounds         : ', &
                 'Bisection         : ', &
                 'RecoverPrimitives : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'ComputePrimitive : ', &
        Timer_MHD_ComputePrimitive, ' s = ', &
        Timer_MHD_ComputePrimitive / Timer_MHD_Program, ' = ', &
        Timer_MHD_ComputePrimitive / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_ComputePrimitive )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Troubled-Cell Indicator'
      WRITE(*,TRIM(Label)) '-----------------------'
      WRITE(*,*)

      nT = 5
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_DD_TCI_DetectTroubledCells, &
                 Timer_MHD_DD_TCI_CopyIn, &
                 Timer_MHD_DD_TCI_CopyOut, &
                 Timer_MHD_DD_TCI_Permute, &
                 Timer_MHD_DD_TCI_Integrate ]

      Labels = [ 'Detect Troubled Cells : ', &
                 'CopyIn                : ', &
                 'CopyOut               : ', &
                 'Permute               : ', &
                 'Integrate             : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'Troubled-Cell Indicator : ', &
        Timer_MHD_DD_TCI, ' s = ', &
        Timer_MHD_DD_TCI / Timer_MHD_Program, ' = ', &
        Timer_MHD_DD_TCI / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_DD_TCI )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Shock Detector'
      WRITE(*,TRIM(Label)) '--------------'
      WRITE(*,*)

      nT = 7
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_DD_SD_DetectShocks, &
                 Timer_MHD_DD_SD_CopyIn, &
                 Timer_MHD_DD_SD_CopyOut, &
                 Timer_MHD_DD_SD_Permute, &
                 Timer_MHD_DD_SD_Integrate, &
                 Timer_MHD_DD_SD_ComputePrimitive, &
                 Timer_MHD_DD_SD_ErrorCheck ]

      Labels = [ 'Detect Shocks     : ', &
                 'CopyIn            : ', &
                 'CopyOut           : ', &
                 'Permute           : ', &
                 'Integrate         : ', &
                 'Compute Primitive : ', &
                 'Error Check       : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'Shock Detector : ', &
        Timer_MHD_DD_SD, ' s = ', &
        Timer_MHD_DD_SD / Timer_MHD_Program, ' = ', &
        Timer_MHD_DD_SD / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_DD_SD )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Slope-Limiter'
      WRITE(*,TRIM(Label)) '-------------'
      WRITE(*,*)

      nT = 8
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_SL_LimitCells, &
                 Timer_MHD_SL_CharDecomp, &
                 Timer_MHD_SL_ConsCorr, &
                 Timer_MHD_SL_Mapping, &
                 Timer_MHD_SL_CopyIn, &
                 Timer_MHD_SL_CopyOut, &
                 Timer_MHD_SL_Permute, &
                 Timer_MHD_SL_Integrate ]

      Labels = [ 'Limit Cells                  : ', &
                 'Characteristic Decomposition : ', &
                 'Polynomial Mapping           : ', &
                 'Conservative Correction      : ', &
                 'CopyIn                       : ', &
                 'CopyOut                      : ', &
                 'Permute                      : ', &
                 'Integrate                    : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'Slope-Limiter : ', &
        Timer_MHD_SlopeLimiter, ' s = ', &
        Timer_MHD_SlopeLimiter / Timer_MHD_Program, ' = ', &
        Timer_MHD_SlopeLimiter / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_SlopeLimiter )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      nT = 6
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Positivity-Limiter'
      WRITE(*,TRIM(Label)) '------------------'
      WRITE(*,*)

      Timers = [ Timer_MHD_PL_LimitCells, &
                 Timer_MHD_PL_CopyIn, &
                 Timer_MHD_PL_CopyOut, &
                 Timer_MHD_PL_Permute, &
                 Timer_MHD_PL_Integrate, &
                 Timer_MHD_PL_ErrorCheck ]

      Labels = [ 'Limit Cells : ', &
                 'CopyIn      : ', &
                 'CopyOut     : ', &
                 'Permute     : ', &
                 'Integrate   : ', &
                 'Error Check : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'Positivity-Limiter : ', &
        Timer_MHD_PositivityLimiter, ' s = ', &
        Timer_MHD_PositivityLimiter / Timer_MHD_Program, ' = ', &
        Timer_MHD_PositivityLimiter / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_PositivityLimiter )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Compute From Conserved'
      WRITE(*,TRIM(Label)) '----------------------'
      WRITE(*,*)

      nT = 4
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      Timers = [ Timer_MHD_CFC_ComputePrimitive, &
                 Timer_MHD_CFC_ErrorCheck, &
                 Timer_MHD_CFC_CopyIn, &
                 Timer_MHD_CFC_CopyOut ]

      Labels = [ 'Compute Primitive : ', &
                 'ErrorCheck        : ', &
                 'CopyIn            : ', &
                 'CopyOut           : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'ComputeFromConserved : ', &
        Timer_MHD_ComputeFromConserved, ' s = ', &
        Timer_MHD_ComputeFromConserved / Timer_MHD_Program, ' = ', &
        Timer_MHD_ComputeFromConserved / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_ComputeFromConserved )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      nT = 3
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Compute Time Step'
      WRITE(*,TRIM(Label)) '-----------------'
      WRITE(*,*)

      Timers = [ Timer_MHD_CTS_ComputeTimeStep, &
                 Timer_MHD_CTS_CopyIn, &
                 Timer_MHD_CTS_CopyOut ]

      Labels = [ 'Compute TimeStep : ', &
                 'CopyIn           : ', &
                 'CopyOut          : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'ComputeTimeStep : ', &
        Timer_MHD_ComputeTimeStep, ' s = ', &
        Timer_MHD_ComputeTimeStep / Timer_MHD_Program, ' = ', &
        Timer_MHD_ComputeTimeStep / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_ComputeTimeStep )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      nT = 3
      ALLOCATE( Timers(nT) )
      ALLOCATE( Labels(nT) )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Boundary Conditions'
      WRITE(*,TRIM(Label)) '-------------------'
      WRITE(*,*)

      Timers = [ Timer_MHD_BC_ApplyBC, &
                 Timer_MHD_BC_CopyIn, &
                 Timer_MHD_BC_CopyOut ]

      Labels = [ 'Apply Boundary Conditions : ', &
                 'CopyIn                    : ', &
                 'CopyOut                   : ' ]

      TotalTime = SUM( Timers )

      WRITE(*,TRIM(TimeL2)) &
        'BoundaryConditions: ', &
        Timer_MHD_BoundaryConditions, ' s = ', &
        Timer_MHD_BoundaryConditions / Timer_MHD_Program, ' = ', &
        Timer_MHD_BoundaryConditions / TotalTime

      WRITE(*,*)

      DO iT = 1, nT

        CALL WriteTimer &
               ( TimeL2, Labels(iT), &
                 Timers(iT), Timer_MHD_Program, &
                 SubroutineRunTime_Option = Timer_MHD_BoundaryConditions )

      END DO

      WRITE(*,*)

      DEALLOCATE( Labels )
      DEALLOCATE( Timers )

      WRITE(*,*)
      WRITE(*,TRIM(Label)) 'Gravity Solver'
      WRITE(*,TRIM(Label)) '--------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'Solve Gravity : ', &
        Timer_GravitySolver, ' s = ', &
        Timer_GravitySolver / Timer_MHD_Program

      WRITE(*,TRIM(TimeL1)) &
        'ComputeSourceTerms : ', &
        Timer_GS_ComputeSourceTerms, ' s = ', &
        Timer_GS_ComputeSourceTerms / Timer_MHD_Program

    END IF

    IF( WriteAtIntermediateTime ) &
      CALL TimersStart_MHD( Timer_MHD_Program )

  END SUBROUTINE FinalizeTimers_MHD


  SUBROUTINE TimersStart_MHD( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_MHD ) RETURN

    Timer = Timer - TimersWtime_MHD()

    RETURN
  END SUBROUTINE TimersStart_MHD


  SUBROUTINE TimersStop_MHD( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_MHD ) RETURN

    Timer = Timer + TimersWtime_MHD()

    RETURN
  END SUBROUTINE TimersStop_MHD


  REAL(DP) FUNCTION TimersWtime_MHD()

    INTEGER(I8) :: clock_read
    INTEGER(I8) :: clock_rate
    INTEGER(I8) :: clock_max

    IF( .NOT. TimeIt_MHD )THEN

      TimersWtime_MHD = -HUGE( 1.0_DP )

      RETURN

    END IF

    CALL SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
    TimersWtime_MHD = REAL( clock_read, DP ) / REAL( clock_rate, DP )

    RETURN
  END FUNCTION TimersWtime_MHD


  SUBROUTINE WriteTimer &
    ( FMT, Label, Timer, ProgramRunTime, SubroutineRunTime_Option )

    CHARACTER(64), INTENT(in) :: FMT
    CHARACTER(64), INTENT(in) :: Label
    REAL(DP),      INTENT(in) :: Timer, ProgramRunTime
    REAL(DP),      INTENT(in), OPTIONAL :: SubroutineRunTime_Option

    IF( PRESENT( SubroutineRunTime_Option ) )THEN

      WRITE(*,TRIM( FMT )) &
        TRIM( Label ), &
        Timer, ' s = ', &
        Timer / ProgramRunTime, ' = ', &
        Timer / SubroutineRunTime_Option

    ELSE

      WRITE(*,TRIM( FMT )) &
        TRIM( Label ), &
        Timer, ' s = ', &
        Timer / ProgramRunTime

    END IF

  END SUBROUTINE WriteTimer


END MODULE TimersModule_MHD
