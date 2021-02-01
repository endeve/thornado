MODULE TimersModule_Euler

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64
  USE KindModule, Only: &
    DP, Zero

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

  REAL(DP), PUBLIC :: Timer_Euler_Increment
  REAL(DP), PUBLIC :: Timer_Euler_Divergence
  REAL(DP), PUBLIC :: Timer_Euler_SurfaceTerm
  REAL(DP), PUBLIC :: Timer_Euler_VolumeTerm
  REAL(DP), PUBLIC :: Timer_Euler_Geometry
  REAL(DP), PUBLIC :: Timer_Euler_Gravity

  ! --- Data Manipulation ---

  REAL(DP), PUBLIC :: Timer_Euler_CopyIn
  REAL(DP), PUBLIC :: Timer_Euler_CopyOut
  REAL(DP), PUBLIC :: Timer_Euler_Permute
  REAL(DP), PUBLIC :: Timer_Euler_Interpolate

  ! --- Limiter-specific ---

  REAL(DP), PUBLIC :: Timer_Euler_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_Euler_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_Euler_TroubledCellIndicator
  REAL(DP), PUBLIC :: Timer_Euler_ShockDetector
  REAL(DP), PUBLIC :: Timer_Euler_CharacteristicDecomposition

  ! --- Miscellaneous ---

  REAL(DP), PUBLIC :: Timer_Euler_BoundaryConditions
  REAL(DP), PUBLIC :: Timer_GravitySolver


  PUBLIC :: InitializeTimers_Euler
  PUBLIC :: FinalizeTimers_Euler
  PUBLIC :: TimersStart_Euler
  PUBLIC :: TimersStop_Euler
  PUBLIC :: TimersWtime_Euler


CONTAINS


  SUBROUTINE InitializeTimers_Euler

    IF( .NOT. TimeIt_Euler ) RETURN

    Timer_Euler_Program         = Zero

    Timer_Euler_Initialize      = Zero
    Timer_Euler_ComputeTimeStep = Zero
    Timer_Euler_UpdateFluid     = Zero
    Timer_Euler_InputOutput     = Zero
    Timer_Euler_Finalize        = Zero

    Timer_Euler_Increment   = Zero
    Timer_Euler_Divergence  = Zero
    Timer_Euler_SurfaceTerm = Zero
    Timer_Euler_VolumeTerm  = Zero
    Timer_Euler_Geometry    = Zero
    Timer_Euler_Gravity     = Zero

    Timer_Euler_CopyIn      = Zero
    Timer_Euler_CopyOut     = Zero
    Timer_Euler_Permute     = Zero
    Timer_Euler_Interpolate = Zero

    Timer_Euler_PositivityLimiter           = Zero
    Timer_Euler_SlopeLimiter                = Zero
    Timer_Euler_TroubledCellIndicator       = Zero
    Timer_Euler_ShockDetector               = Zero
    Timer_Euler_CharacteristicDecomposition = Zero

    Timer_GravitySolver            = Zero
    Timer_Euler_BoundaryConditions = Zero

    CALL TimersStart_Euler( Timer_Euler_Program )

  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler &
    ( Verbose_Option, SuppressApplicationDriver_Option, &
      WriteAtIntermediateTime_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL, INTENT(in), OPTIONAL :: SuppressApplicationDriver_Option
    LOGICAL, INTENT(in), OPTIONAL :: WriteAtIntermediateTime_Option

    LOGICAL  :: Verbose, SuppressApplicationDriver, WriteAtIntermediateTime
    REAL(DP) :: TotalTime

    CHARACTER(6)  :: Label_Level1 = '(8x,A)'

    CHARACTER(64) :: TimeL1 = '(10x,A,ES13.6E3,A,F6.3,A)'
    CHARACTER(64) :: TimeL2 = '(12x,A,ES13.6E3,A,F6.3,A,F6.3,A)'

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
      WRITE(*,'(4x,A)') 'Timers Summary (Euler)'
      WRITE(*,'(4x,A)') '----------------------'
      WRITE(*,*)

      WRITE(*,'(6x,A,ES13.6E3,A)') &
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

      WRITE(*,'(10x,A,ES13.6E3,A,F7.3,A)') &
        'Timers = ', TotalTime, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'Initialize:        ', &
        Timer_Euler_Initialize, ' s = ', &
        100.0_DP &
        * Timer_Euler_Initialize / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Compute Time-Step: ', &
        Timer_Euler_ComputeTimeStep, ' s = ', &
        100.0_DP &
        * Timer_Euler_ComputeTimeStep / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Update Fluid:      ', &
        Timer_Euler_UpdateFluid, ' s = ', &
        100.0_DP &
          * Timer_Euler_UpdateFluid / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Input/Output:      ', &
        Timer_Euler_InputOutput, ' s = ', &
        100.0_DP &
          * Timer_Euler_InputOutput / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Finalize:          ', &
        Timer_Euler_Finalize, ' s = ', &
        100.0_DP &
          * Timer_Euler_Finalize / Timer_Euler_Program, ' %'

    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'DG discretization'
      WRITE(*,TRIM(Label_Level1)) '---------------- '
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'Divergence: ', &
        Timer_Euler_Divergence, ' s = ', &
        100.0_DP &
          * Timer_Euler_Divergence / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL2)) &
        'Surface Term: ', &
        Timer_Euler_SurfaceTerm, ' s = ', &
        100.0_DP &
          * Timer_Euler_SurfaceTerm / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL2)) &
        'Volume Term:  ', &
        Timer_Euler_VolumeTerm, ' s = ', &
        100.0_DP &
          * Timer_Euler_VolumeTerm / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(TimeL1)) &
        'Geometry: ', &
        Timer_Euler_Geometry, ' s = ', &
        100.0_DP &
          * Timer_Euler_Geometry / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Gravity:  ', &
        Timer_Euler_Gravity, ' s = ', &
        100.0_DP &
          * Timer_Euler_Gravity / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Data Manipulation'
      WRITE(*,TRIM(Label_Level1)) '-----------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'CopyIn:      ', &
        Timer_Euler_CopyIn, ' s = ', &
        100.0_DP &
          * Timer_Euler_CopyIn / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'CopyOut:     ', &
        Timer_Euler_CopyOut, ' s = ', &
        100.0_DP &
          * Timer_Euler_CopyOut / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Permute:     ', &
        Timer_Euler_Permute, ' s = ', &
        100.0_DP &
          * Timer_Euler_Permute / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Interpolate: ', &
        Timer_Euler_Permute, ' s = ', &
        100.0_DP &
          * Timer_Euler_Permute / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Limiter-specific'
      WRITE(*,TRIM(Label_Level1)) '----------------'

      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'Positivity-Limiter: ', &
        Timer_Euler_PositivityLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_PositivityLimiter / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Slope-Limiter:      ', &
        Timer_Euler_SlopeLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_SlopeLimiter / Timer_Euler_Program, ' %'

      WRITE(*,*)

      WRITE(*,TRIM(TimeL2)) &
        'Troubled-Cell Indicator:      ', &
        Timer_Euler_TroubledCellIndicator, ' s = ', &
        100.0_DP &
          * Timer_Euler_TroubledCellIndicator / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL2)) &
        'Shock Detector:               ', &
        Timer_Euler_ShockDetector, ' s = ', &
        100.0_DP &
          * Timer_Euler_ShockDetector / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL2)) &
        'Characteristic Decomposition: ', &
        Timer_Euler_CharacteristicDecomposition, ' s = ', &
        100.0_DP &
          * Timer_Euler_CharacteristicDecomposition / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Miscellaneous'
      WRITE(*,TRIM(Label_Level1)) '-------------'
      WRITE(*,*)

      WRITE(*,TRIM(TimeL1)) &
        'GravitySolver:       ', &
        Timer_GravitySolver, ' s = ', &
        100.0_DP &
          * Timer_GravitySolver / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeL1)) &
        'Boundary Conditions: ', &
        Timer_Euler_BoundaryConditions, ' s = ', &
        100.0_DP &
          * Timer_Euler_BoundaryConditions / Timer_Euler_Program, ' %'

    END IF

    IF( WriteAtIntermediateTime ) &
      CALL TimersStart_Euler( Timer_Euler_Program )

    RETURN

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
