MODULE TimersModule_Euler

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64
  USE KindModule, Only: &
    DP, SqrtTiny

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PUBLIC :: TimeIt_Euler = .FALSE.

  REAL(DP), PUBLIC :: Timer_Euler_Program

  ! --- ApplicationDriver ---

  REAL(DP), PUBLIC :: Timer_Euler_Initialize
  REAL(DP), PUBLIC :: Timer_Euler_ComputeTimeStep
  REAL(DP), PUBLIC :: Timer_Euler_UpdateFluid
  REAL(DP), PUBLIC :: Timer_Euler_InputOutput
  REAL(DP), PUBLIC :: Timer_Euler_Finalize

  ! --- DG discretization ---

  REAL(DP), PUBLIC :: Timer_Euler_DG
  REAL(DP), PUBLIC :: Timer_Euler_Increment
  REAL(DP), PUBLIC :: Timer_Euler_Divergence
  REAL(DP), PUBLIC :: Timer_Euler_SurfaceTerm
  REAL(DP), PUBLIC :: Timer_Euler_VolumeTerm
  REAL(DP), PUBLIC :: Timer_Euler_ComputePrimitive
  REAL(DP), PUBLIC :: Timer_Euler_Geometry
  REAL(DP), PUBLIC :: Timer_Euler_Gravity
  REAL(DP), PUBLIC :: Timer_Euler_DG_CopyIn
  REAL(DP), PUBLIC :: Timer_Euler_DG_CopyOut
  REAL(DP), PUBLIC :: Timer_Euler_DG_Permute
  REAL(DP), PUBLIC :: Timer_Euler_DG_Interpolate
  REAL(DP), PUBLIC :: Timer_Euler_DG_ErrorCheck


  ! --- Data Manipulation ---

  REAL(DP), PUBLIC :: Timer_Euler_CopyIn
  REAL(DP), PUBLIC :: Timer_Euler_CopyOut
  REAL(DP), PUBLIC :: Timer_Euler_Permute
  REAL(DP), PUBLIC :: Timer_Euler_Interpolate

  ! --- Discontinuity Detection ---

  ! --- Troubled-Cell Indicator

  REAL(DP), PUBLIC :: Timer_Euler_DD_TCI
  REAL(DP), PUBLIC :: Timer_Euler_DD_TCI_DetectTroubledCells
  REAL(DP), PUBLIC :: Timer_Euler_DD_TCI_CopyIn
  REAL(DP), PUBLIC :: Timer_Euler_DD_TCI_CopyOut
  REAL(DP), PUBLIC :: Timer_Euler_DD_TCI_Permute
  REAL(DP), PUBLIC :: Timer_Euler_DD_TCI_Integrate

  ! --- Shock Detector ---

  REAL(DP), PUBLIC :: Timer_Euler_DD_SD
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_DetectShocks
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_CopyIn
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_CopyOut
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_Permute
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_Integrate
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_ComputePrimitive
  REAL(DP), PUBLIC :: Timer_Euler_DD_SD_ErrorCheck

  ! --- Limiter-specific ---

  ! --- Slope-Limiter ---

  REAL(DP), PUBLIC :: Timer_Euler_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_Euler_SL_LimitCells
  REAL(DP), PUBLIC :: Timer_Euler_SL_CharDecomp
  REAL(DP), PUBLIC :: Timer_Euler_SL_ConsCorr
  REAL(DP), PUBLIC :: Timer_Euler_SL_Mapping
  REAL(DP), PUBLIC :: Timer_Euler_SL_CopyIn
  REAL(DP), PUBLIC :: Timer_Euler_SL_CopyOut
  REAL(DP), PUBLIC :: Timer_Euler_SL_Permute
  REAL(DP), PUBLIC :: Timer_Euler_SL_Integrate

  ! --- Positivity-Limiter ---

  REAL(DP), PUBLIC :: Timer_Euler_PositivityLimiter

  ! --- Miscellaneous ---

  REAL(DP), PUBLIC :: Timer_Euler_BoundaryConditions
  REAL(DP), PUBLIC :: Timer_GravitySolver


  PUBLIC :: InitializeTimers_Euler
  PUBLIC :: FinalizeTimers_Euler
  PUBLIC :: TimersStart_Euler
  PUBLIC :: TimersStop_Euler
  PUBLIC :: TimersWtime_Euler

  REAL(DP), PARAMETER :: Hundred = 100.0_DP


CONTAINS


  SUBROUTINE InitializeTimers_Euler

    IF( .NOT. TimeIt_Euler ) RETURN

    Timer_Euler_Program = SqrtTiny

    CALL TimersStart_Euler( Timer_Euler_Program )

    Timer_Euler_Initialize      = SqrtTiny
    Timer_Euler_ComputeTimeStep = SqrtTiny
    Timer_Euler_UpdateFluid     = SqrtTiny
    Timer_Euler_InputOutput     = SqrtTiny
    Timer_Euler_Finalize        = SqrtTiny

    Timer_Euler_DG               = SqrtTiny
    Timer_Euler_Increment        = SqrtTiny
    Timer_Euler_Divergence       = SqrtTiny
    Timer_Euler_SurfaceTerm      = SqrtTiny
    Timer_Euler_VolumeTerm       = SqrtTiny
    Timer_Euler_ComputePrimitive = SqrtTiny
    Timer_Euler_Geometry         = SqrtTiny
    Timer_Euler_Gravity          = SqrtTiny
    Timer_Euler_DG_CopyIn        = SqrtTiny
    Timer_Euler_DG_CopyOut       = SqrtTiny
    Timer_Euler_DG_Permute       = SqrtTiny
    Timer_Euler_DG_Interpolate   = SqrtTiny
    Timer_Euler_DG_ErrorCheck    = SqrtTiny

    Timer_Euler_CopyIn      = SqrtTiny
    Timer_Euler_CopyOut     = SqrtTiny
    Timer_Euler_Permute     = SqrtTiny
    Timer_Euler_Interpolate = SqrtTiny

    Timer_Euler_DD_TCI                     = SqrtTiny
    Timer_Euler_DD_TCI_DetectTroubledCells = SqrtTiny
    Timer_Euler_DD_TCI_CopyIn              = SqrtTiny
    Timer_Euler_DD_TCI_CopyOut             = SqrtTiny
    Timer_Euler_DD_TCI_Permute             = SqrtTiny
    Timer_Euler_DD_TCI_Integrate           = SqrtTiny

    Timer_Euler_DD_SD                  = SqrtTiny
    Timer_Euler_DD_SD_DetectShocks     = SqrtTiny
    Timer_Euler_DD_SD_CopyIn           = SqrtTiny
    Timer_Euler_DD_SD_CopyOut          = SqrtTiny
    Timer_Euler_DD_SD_Permute          = SqrtTiny
    Timer_Euler_DD_SD_Integrate        = SqrtTiny
    Timer_Euler_DD_SD_ComputePrimitive = SqrtTiny
    Timer_Euler_DD_SD_ErrorCheck       = SqrtTiny

    Timer_Euler_SlopeLimiter  = SqrtTiny
    Timer_Euler_SL_LimitCells = SqrtTiny
    Timer_Euler_SL_CharDecomp = SqrtTiny
    Timer_Euler_SL_ConsCorr   = SqrtTiny
    Timer_Euler_SL_Mapping    = SqrtTiny
    Timer_Euler_SL_CopyIn     = SqrtTiny
    Timer_Euler_SL_CopyOut    = SqrtTiny
    Timer_Euler_SL_Integrate  = SqrtTiny
    Timer_Euler_SL_Permute    = SqrtTiny

    Timer_Euler_PositivityLimiter = SqrtTiny

    Timer_GravitySolver            = SqrtTiny
    Timer_Euler_BoundaryConditions = SqrtTiny

  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler &
    ( Verbose_Option, SuppressApplicationDriver_Option, &
      WriteAtIntermediateTime_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL, INTENT(in), OPTIONAL :: SuppressApplicationDriver_Option
    LOGICAL, INTENT(in), OPTIONAL :: WriteAtIntermediateTime_Option

    LOGICAL  :: Verbose, SuppressApplicationDriver, WriteAtIntermediateTime
    REAL(DP) :: TotalTime

    CHARACTER(6)  :: Label_Level1 = '(6x,A)'

    CHARACTER(64) :: TimeL3 = '(8x,A,ES10.3E3,A,ES10.3E3,A,ES10.3E3)'
    CHARACTER(64) :: TimeL1 = '(8x,A,ES10.3E3,A,F6.3,A)'
    CHARACTER(64) :: TimeL2 = '(10x,A,ES10.3E3,A,F6.3,A,F6.3,A)'

    IF( .NOT. TimeIt_Euler ) RETURN

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    SuppressApplicationDriver = .FALSE.
    IF( PRESENT( SuppressApplicationDriver_Option ) ) &
      SuppressApplicationDriver = SuppressApplicationDriver_Option

    WriteAtIntermediateTime = .FALSE.
    IF( PRESENT( WriteAtIntermediateTime_Option ) ) &
      WriteAtIntermediateTime = WriteAtIntermediateTime_Option

    CALL TimersStop_Euler( Timer_Euler_Program )

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(2x,A)') 'Timers Summary (Euler)'
      WRITE(*,'(2x,A)') '----------------------'
      WRITE(*,*)

      WRITE(*,'(4x,A,ES13.6E3,A)') &
        'Total run-time = ', Timer_Euler_Program, ' s'

    END IF

    IF( .NOT. SuppressApplicationDriver )THEN

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Application Driver'
      WRITE(*,TRIM(Label_Level1)) '------------------'
      WRITE(*,*)
      TotalTime = Timer_Euler_Initialize &
                    + Timer_Euler_ComputeTimeStep &
                    + Timer_Euler_UpdateFluid &
                    + Timer_Euler_InputOutput &
                    + Timer_Euler_Finalize

      WRITE(*,'(8x,A,ES13.6E3,A,F7.3,A)') &
        'Timers = ', TotalTime, ' s = ', &
        Hundred * TotalTime / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'Initialize:        ', &
        Timer_Euler_Initialize, ' s = ', &
        Hundred &
        * Timer_Euler_Initialize / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Compute Time-Step: ', &
        Timer_Euler_ComputeTimeStep, ' s = ', &
        Hundred &
        * Timer_Euler_ComputeTimeStep / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Update Fluid:      ', &
        Timer_Euler_UpdateFluid, ' s = ', &
        Hundred &
          * Timer_Euler_UpdateFluid / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Input/Output:      ', &
        Timer_Euler_InputOutput, ' s = ', &
        Hundred &
          * Timer_Euler_InputOutput / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Finalize:          ', &
        Timer_Euler_Finalize, ' s = ', &
        Hundred &
          * Timer_Euler_Finalize / Timer_Euler_Program, ' %'

    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'DG discretization'
      WRITE(*,TRIM(Label_Level1)) '-----------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'dgDiscretization: ', &
        Timer_Euler_DG, ' s = ', &
        Hundred &
          * Timer_Euler_DG / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        '  Increment:    ', &
        Timer_Euler_Increment, ' s = ', &
        Hundred &
          * Timer_Euler_Increment / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        '  Geometry:     ', &
        Timer_Euler_Geometry, ' s = ', &
        Hundred &
          * Timer_Euler_Geometry / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        '  Gravity:      ', &
        Timer_Euler_Gravity, ' s = ', &
        Hundred &
          * Timer_Euler_Gravity / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        '  Divergence:   ', &
        Timer_Euler_Divergence, ' s = ', &
        Hundred &
          * Timer_Euler_Divergence / Timer_Euler_Program, ' %'
      WRITE(*,*)

     WRITE(*,TRIM(TimeL2)) &
       '  Surface Term: ', &
       Timer_Euler_SurfaceTerm, ' s = ', &
       Hundred &
         * Timer_Euler_SurfaceTerm / Timer_Euler_Program, ' % = ', &
       Hundred &
         * Timer_Euler_SurfaceTerm / Timer_Euler_Divergence, ' %'

     WRITE(*,TRIM(TimeL2)) &
       '  Volume Term:  ', &
       Timer_Euler_VolumeTerm, ' s = ', &
       Hundred &
         * Timer_Euler_VolumeTerm / Timer_Euler_Program, ' % = ', &
       Hundred &
         * Timer_Euler_VolumeTerm / Timer_Euler_Divergence, ' %'

     WRITE(*,*)

     WRITE(*,TRIM(TimeL2)) &
       '    ComputePrimitive: ', &
       Timer_Euler_ComputePrimitive, ' s = ', &
       Hundred &
         * Timer_Euler_ComputePrimitive / Timer_Euler_Program, ' % = ', &
       Hundred &
         * Timer_Euler_ComputePrimitive &
         / ( Timer_Euler_SurfaceTerm + Timer_Euler_VolumeTerm ), ' %'

      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        '  CopyIn:       ', &
        Timer_Euler_DG_CopyIn, ' s = ', &
        Hundred &
          * Timer_Euler_DG_CopyIn / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        '  CopyOut:      ', &
        Timer_Euler_DG_CopyOut, ' s = ', &
        Hundred &
          * Timer_Euler_DG_CopyOut / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        '  Interpolate:  ', &
        Timer_Euler_DG_Interpolate, ' s = ', &
        Hundred &
          * Timer_Euler_DG_Interpolate / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        '  Permute:      ', &
        Timer_Euler_DG_Permute, ' s = ', &
        Hundred &
          * Timer_Euler_DG_Permute / Timer_Euler_Program, ' %'

      Timer_Euler_CopyIn &
        = Timer_Euler_CopyIn + Timer_Euler_DG_CopyIn

      Timer_Euler_CopyOut &
        = Timer_Euler_CopyOut + Timer_Euler_DG_CopyOut

      Timer_Euler_Interpolate &
        = Timer_Euler_Interpolate + Timer_Euler_DG_Interpolate

      Timer_Euler_Permute &
        = Timer_Euler_Permute + Timer_Euler_DG_Permute

      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        '  Error Check:  ', &
        Timer_Euler_DG_ErrorCheck, ' s = ', &
        Hundred &
          * Timer_Euler_DG_ErrorCheck / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Data Manipulation'
      WRITE(*,TRIM(Label_Level1)) '-----------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'CopyIn:      ', &
        Timer_Euler_CopyIn, ' s = ', &
        Hundred &
          * Timer_Euler_CopyIn / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'CopyOut:     ', &
        Timer_Euler_CopyOut, ' s = ', &
        Hundred &
          * Timer_Euler_CopyOut / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Permute:     ', &
        Timer_Euler_Permute, ' s = ', &
        Hundred &
          * Timer_Euler_Permute / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Interpolate: ', &
        Timer_Euler_Permute, ' s = ', &
        Hundred &
          * Timer_Euler_Permute / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Limiter-specific'
      WRITE(*,TRIM(Label_Level1)) '----------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'Positivity-Limiter: ', &
        Timer_Euler_PositivityLimiter, ' s = ', &
        Hundred &
          * Timer_Euler_PositivityLimiter / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Troubled-Cell Indicator'
      WRITE(*,TRIM(Label_Level1)) '-----------------------'
      WRITE(*,*)

      TotalTime = Timer_Euler_DD_TCI_DetectTroubledCells &
                    + Timer_Euler_DD_TCI_CopyIn &
                    + Timer_Euler_DD_TCI_CopyOut &
                    + Timer_Euler_DD_TCI_Permute &
                    + Timer_Euler_DD_TCI_Integrate

      WRITE(*,TRIM(TimeL3)) &
        'Troubled-Cell Indicator: ', &
        Timer_Euler_DD_TCI, ' s = ', &
        Timer_Euler_DD_TCI / Timer_Euler_Program, ' = ', &
        TotalTime / Timer_Euler_DD_TCI

      WRITE(*,*)

      WRITE(*,TRIM(TimeL3)) &
        '  Detect Troubled Cells: ', &
        Timer_Euler_DD_TCI_DetectTroubledCells, ' s = ', &
        Timer_Euler_DD_TCI_DetectTroubledCells / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_TCI_DetectTroubledCells / Timer_Euler_DD_TCI

      WRITE(*,TRIM(TimeL3)) &
        '  CopyIn:                ', &
        Timer_Euler_DD_TCI_CopyIn, ' s = ', &
        Timer_Euler_DD_TCI_CopyIn / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_TCI_CopyIn / Timer_Euler_DD_TCI

      WRITE(*,TRIM(TimeL3)) &
        '  CopyOut:               ', &
        Timer_Euler_DD_TCI_CopyOut, ' s = ', &
        Timer_Euler_DD_TCI_CopyOut / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_TCI_CopyOut / Timer_Euler_DD_TCI

      WRITE(*,TRIM(TimeL3)) &
        '  Permute:               ', &
        Timer_Euler_DD_TCI_Permute, ' s = ', &
        Timer_Euler_DD_TCI_Permute / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_TCI_Permute / Timer_Euler_DD_TCI

      WRITE(*,TRIM(TimeL3)) &
        '  Integrate:             ', &
        Timer_Euler_DD_TCI_Integrate, ' s = ', &
        Timer_Euler_DD_TCI_Integrate / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_TCI_Integrate / Timer_Euler_DD_TCI

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Shock Detector'
      WRITE(*,TRIM(Label_Level1)) '--------------'
      WRITE(*,*)

      TotalTime = Timer_Euler_DD_SD_DetectShocks &
                    + Timer_Euler_DD_SD_CopyIn &
                    + Timer_Euler_DD_SD_CopyOut &
                    + Timer_Euler_DD_SD_Permute &
                    + Timer_Euler_DD_SD_Integrate &
                    + Timer_Euler_DD_SD_ComputePrimitive &
                    + Timer_Euler_DD_SD_ErrorCheck

      WRITE(*,TRIM(TimeL3)) &
        'Shock Detector: ', &
        Timer_Euler_DD_SD, ' s = ', &
        Timer_Euler_DD_SD / Timer_Euler_Program, ' = ', &
        TotalTime / Timer_Euler_DD_SD

      WRITE(*,*)

      WRITE(*,TRIM(TimeL3)) &
        '  Detect Shocks:     ', &
        Timer_Euler_DD_SD_DetectShocks, ' s = ', &
        Timer_Euler_DD_SD_DetectShocks / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_DetectShocks / Timer_Euler_DD_SD

      WRITE(*,TRIM(TimeL3)) &
        '  CopyIn:            ', &
        Timer_Euler_DD_SD_CopyIn, ' s = ', &
        Timer_Euler_DD_SD_CopyIn / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_CopyIn / Timer_Euler_DD_SD

      WRITE(*,TRIM(TimeL3)) &
        '  CopyOut:           ', &
        Timer_Euler_DD_SD_CopyOut, ' s = ', &
        Timer_Euler_DD_SD_CopyOut / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_CopyOut / Timer_Euler_DD_SD

      WRITE(*,TRIM(TimeL3)) &
        '  Permute:           ', &
        Timer_Euler_DD_SD_Permute, ' s = ', &
        Timer_Euler_DD_SD_Permute / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_Permute / Timer_Euler_DD_SD

      WRITE(*,TRIM(TimeL3)) &
        '  Integrate:         ', &
        Timer_Euler_DD_SD_Integrate, ' s = ', &
        Timer_Euler_DD_SD_Integrate / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_Integrate / Timer_Euler_DD_SD

      WRITE(*,TRIM(TimeL3)) &
        '  Compute Primitive: ', &
        Timer_Euler_DD_SD_ComputePrimitive, ' s = ', &
        Timer_Euler_DD_SD_ComputePrimitive / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_ComputePrimitive / Timer_Euler_DD_SD

      WRITE(*,TRIM(TimeL3)) &
        '  Error Check:       ', &
        Timer_Euler_DD_SD_ErrorCheck, ' s = ', &
        Timer_Euler_DD_SD_ErrorCheck / Timer_Euler_Program, ' = ', &
        Timer_Euler_DD_SD_ErrorCheck / Timer_Euler_DD_SD

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Slope-Limiter'
      WRITE(*,TRIM(Label_Level1)) '-------------'
      WRITE(*,*)

      TotalTime = Timer_Euler_SL_LimitCells &
                    + Timer_Euler_SL_CharDecomp &
                    + Timer_Euler_SL_ConsCorr &
                    + Timer_Euler_SL_Mapping &
                    + Timer_Euler_SL_CopyIn &
                    + Timer_Euler_SL_CopyOut &
                    + Timer_Euler_SL_Integrate &
                    + Timer_Euler_SL_Permute

      WRITE(*,TRIM(TimeL3)) &
        'Slope-Limiter: ', &
        Timer_Euler_SlopeLimiter, ' s = ', &
        Timer_Euler_SlopeLimiter / Timer_Euler_Program, ' = ', &
        TotalTime / Timer_Euler_SlopeLimiter

      WRITE(*,*)

      WRITE(*,TRIM(TimeL3)) &
        '  Limit Cells:                  ', &
        Timer_Euler_SL_LimitCells, ' s = ', &
        Timer_Euler_SL_LimitCells / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_LimitCells / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  Characteristic Decomposition: ', &
        Timer_Euler_SL_CharDecomp, ' s = ', &
        Timer_Euler_SL_CharDecomp / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_CharDecomp / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  Polynomial Mapping:           ', &
        Timer_Euler_SL_Mapping, ' s = ', &
        Timer_Euler_SL_Mapping / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_Mapping / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  Conservative Correction:      ', &
        Timer_Euler_SL_ConsCorr, ' s = ', &
        Timer_Euler_SL_ConsCorr / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_ConsCorr / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  CopyIn:                       ', &
        Timer_Euler_SL_CopyIn, ' s = ', &
        Timer_Euler_SL_CopyIn / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_CopyIn / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  CopyOut:                      ', &
        Timer_Euler_SL_CopyOut, ' s = ', &
        Timer_Euler_SL_CopyOut / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_CopyOut / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  Permute:                      ', &
        Timer_Euler_SL_Permute, ' s = ', &
        Timer_Euler_SL_Permute / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_Permute / Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(TimeL3)) &
        '  Integrate:                    ', &
        Timer_Euler_SL_Integrate, ' s = ', &
        Timer_Euler_SL_Integrate / Timer_Euler_Program, ' = ', &
        Timer_Euler_SL_Integrate / Timer_Euler_SlopeLimiter

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Miscellaneous'
      WRITE(*,TRIM(Label_Level1)) '-------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'GravitySolver:       ', &
        Timer_GravitySolver, ' s = ', &
        Hundred &
          * Timer_GravitySolver / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Boundary Conditions: ', &
        Timer_Euler_BoundaryConditions, ' s = ', &
        Hundred &
          * Timer_Euler_BoundaryConditions / Timer_Euler_Program, ' %'

    END IF

    IF( WriteAtIntermediateTime ) &
      CALL TimersStart_Euler( Timer_Euler_Program )

  END SUBROUTINE FinalizeTimers_Euler


  SUBROUTINE TimersStart_Euler( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_Euler ) RETURN

    Timer = Timer - TimersWtime_Euler()

    RETURN
  END SUBROUTINE TimersStart_Euler


  SUBROUTINE TimersStop_Euler( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_Euler ) RETURN

    Timer = Timer + TimersWtime_Euler()

    RETURN
  END SUBROUTINE TimersStop_Euler


  REAL(DP) FUNCTION TimersWtime_Euler()

    INTEGER(I8) :: clock_read
    INTEGER(I8) :: clock_rate
    INTEGER(I8) :: clock_max

    IF( .NOT. TimeIt_Euler ) RETURN

    CALL SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
    TimersWtime_Euler = REAL( clock_read, DP ) / REAL( clock_rate, DP )

    RETURN
  END FUNCTION TimersWtime_Euler


END MODULE TimersModule_Euler
