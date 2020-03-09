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
  REAL(DP), PUBLIC :: Timer_Euler_dgDiscretization
  REAL(DP), PUBLIC :: Timer_Euler_Divergence
  REAL(DP), PUBLIC :: Timer_Euler_Geometry
  REAL(DP), PUBLIC :: Timer_Euler_Gravity
  REAL(DP), PUBLIC :: Timer_Euler_SurfaceTerm
  REAL(DP), PUBLIC :: Timer_Euler_NumericalFlux
  REAL(DP), PUBLIC :: Timer_Euler_VolumeTerm
  REAL(DP), PUBLIC :: Timer_Euler_ComputePrimitive

  ! --- CPU ---
  REAL(DP), PUBLIC :: Timer_Euler_MatrixVectorMultiply

  ! --- Limiter-specific ---
  REAL(DP), PUBLIC :: Timer_Euler_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_Euler_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_Euler_TroubledCellIndicator
  REAL(DP), PUBLIC :: Timer_Euler_ShockDetector
  REAL(DP), PUBLIC :: Timer_Euler_CharacteristicDecomposition

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

    Timer_Euler_dgDiscretization = Zero
    Timer_Euler_Divergence       = Zero
    Timer_Euler_Geometry         = Zero
    Timer_Euler_Gravity          = Zero
    Timer_Euler_SurfaceTerm      = Zero
    Timer_Euler_NumericalFlux    = Zero
    Timer_Euler_VolumeTerm       = Zero
    Timer_Euler_Increment        = Zero
    Timer_Euler_ComputePrimitive = Zero

    Timer_Euler_MatrixVectorMultiply = Zero

    Timer_Euler_PositivityLimiter           = Zero
    Timer_Euler_SlopeLimiter                = Zero
    Timer_Euler_CharacteristicDecomposition = Zero
    Timer_Euler_TroubledCellIndicator       = Zero
    Timer_Euler_ShockDetector               = Zero

    CALL TimersStart_Euler( Timer_Euler_Program )

    RETURN
  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler &
    ( Verbose_Option, SuppressApplicationDriver_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL, INTENT(in), OPTIONAL :: SuppressApplicationDriver_Option

    LOGICAL  :: Verbose, SuppressApplicationDriver
    REAL(DP) :: TotalTime

    CHARACTER(6)  :: Label_Level1  = '(8x,A)'
    CHARACTER(7)  :: Label_Level2  = '(10x,A)'

    CHARACTER(64) :: OverallTimeAD = '(10x,A,ES13.6E3,A,F7.3,A)'
    CHARACTER(64) :: TimeAD        = '(12x,A,ES13.6E3,A,F6.3,A)'

    CHARACTER(64) :: OverallTimeDG = '(10x,A,ES13.6E3,A,F6.3,A,F7.3,A)'
    CHARACTER(64) :: TimeDG        = '(12x,A,ES13.6E3,A,F6.3,A,F6.3,A)'

    CHARACTER(64) :: OverallTimeAux = '(12x,A,ES13.6E3,A,F6.3,A,F7.3,A)'
    CHARACTER(64) :: TimeAux        = '(12x,A,ES13.6E3,A,F6.3,A,F6.3,A)'

    CHARACTER(64) :: OverallTimeCPU = '(12x,A,ES13.6E3,A,F6.3,A,F7.3,A)'
    CHARACTER(64) :: TimeCPU        = '(12x,A,ES13.6E3,A,F6.3,A,F6.3,A)'

    CHARACTER(64) :: OverallTimeLim = '(10x,A,ES13.6E3,A,F6.3,A,F7.3,A)'
    CHARACTER(64) :: TimeLim        = '(12x,A,ES13.6E3,A,F6.3,A,F6.3,A)'
    CHARACTER(64) :: TimeLim2       = '(14x,A,ES13.6E3,A,F6.3,A,F6.3,A)'

    IF( .NOT. TimeIt_Euler ) RETURN

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    SuppressApplicationDriver = .FALSE.
    IF( PRESENT( SuppressApplicationDriver_Option ) ) &
      SuppressApplicationDriver = SuppressApplicationDriver_Option

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

      WRITE(*,TRIM(OverallTimeAD)) &
        'Timers = ', TotalTime, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(TimeAD)) &
        'Initialize:        ', &
        Timer_Euler_Initialize, ' s = ', &
        100.0_DP &
        * Timer_Euler_Initialize / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeAD)) &
        'Compute Time-Step: ', &
        Timer_Euler_ComputeTimeStep, ' s = ', &
        100.0_DP &
        * Timer_Euler_ComputeTimeStep / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeAD)) &
        'Update Fluid:      ', &
        Timer_Euler_UpdateFluid, ' s = ', &
        100.0_DP &
          * Timer_Euler_UpdateFluid / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeAD)) &
        'Input/Output:      ', &
        Timer_Euler_InputOutput, ' s = ', &
        100.0_DP &
          * Timer_Euler_InputOutput / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeAD)) &
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
      TotalTime = Timer_Euler_Divergence &
                    + Timer_Euler_Geometry &
                    + Timer_Euler_Gravity

      WRITE(*,'(49x,A)') '%DG'
      WRITE(*,TRIM(OverallTimeDG)) &
        'Timers-DG = ', Timer_Euler_dgDiscretization, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' % === ', &
        100.0_DP * TotalTime / Timer_Euler_dgDiscretization
      WRITE(*,*)

      WRITE(*,'(51x,A)') '%DG'
      WRITE(*,TRIM(TimeDG)) &
        'Divergence: ', &
        Timer_Euler_Divergence, ' s = ', &
        100.0_DP &
          * Timer_Euler_Divergence / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_Divergence / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,TRIM(TimeDG)) &
        'Geometry:   ', &
        Timer_Euler_Geometry, ' s = ', &
        100.0_DP &
          * Timer_Euler_Geometry / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_Geometry / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,TRIM(TimeDG)) &
        'Gravity:    ', &
        Timer_Euler_Gravity, ' s = ', &
        100.0_DP &
          * Timer_Euler_Gravity / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_Gravity / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level2)) 'DG discretization-auxiliary'
      WRITE(*,TRIM(Label_Level2)) '---------------------------'

      WRITE(*,'(58x,A)') '%DG'
      WRITE(*,TRIM(TimeAux)) &
        'Surface Term:      ', &
        Timer_Euler_SurfaceTerm, ' s = ', &
        100.0_DP &
          * Timer_Euler_SurfaceTerm / Timer_Euler_Program, ' % === ', &
        100.0_DP * Timer_Euler_SurfaceTerm / Timer_Euler_Divergence, ' %'

      WRITE(*,TRIM(TimeAux)) &
        'Numerical Flux:    ', &
        Timer_Euler_NumericalFlux, ' s = ', &
        100.0_DP &
          * Timer_Euler_NumericalFlux / Timer_Euler_Program, ' % === ', &
        100.0_DP * Timer_Euler_NumericalFlux / Timer_Euler_Divergence, ' %'

      WRITE(*,TRIM(TimeAux)) &
        'Volume Term:       ', &
        Timer_Euler_VolumeTerm, ' s = ', &
        100.0_DP &
          * Timer_Euler_VolumeTerm / Timer_Euler_Program, ' % === ', &
        100.0_DP * Timer_Euler_VolumeTerm / Timer_Euler_Divergence, ' %'

      WRITE(*,TRIM(TimeAux)) &
        'Compute Primitive: ', &
        Timer_Euler_ComputePrimitive, ' s = ', &
        100.0_DP &
          * Timer_Euler_ComputePrimitive / Timer_Euler_Program, ' % === ', &
        100.0_DP * Timer_Euler_ComputePrimitive &
                    / ( Timer_Euler_SurfaceTerm + Timer_Euler_VolumeTerm ), ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level2)) 'CPU-specific'
      WRITE(*,TRIM(Label_Level2)) '------------'
      TotalTime = Timer_Euler_MatrixVectorMultiply

      WRITE(*,'(63x,A)') '%DG'
      WRITE(*,TRIM(TimeCPU)) &
        'Matrix-Vector Multiply: ', &
        Timer_Euler_MatrixVectorMultiply, ' s = ', &
        100.0_DP &
          * Timer_Euler_MatrixVectorMultiply / Timer_Euler_Program, ' % === ', &
        100.0_DP * Timer_Euler_ComputePrimitive &
                    / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(Label_Level1)) 'Limiter-specific'
      WRITE(*,TRIM(Label_Level1)) '----------------'
      WRITE(*,*)
      TotalTime = Timer_Euler_PositivityLimiter &
                    + Timer_Euler_SlopeLimiter

      WRITE(*,TRIM(OverallTimeLim)) &
        'Timers-Limiters = ', TotalTime, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(TimeLim)) &
        'Positivity-Limiter: ', &
        Timer_Euler_PositivityLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_PositivityLimiter / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeLim)) &
        'Slope-Limiter:      ', &
        Timer_Euler_SlopeLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_SlopeLimiter / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(TimeLim2)) &
        'Troubled-Cell Indicator:      ', &
        Timer_Euler_TroubledCellIndicator, ' s = ', &
        100.0_DP &
          * Timer_Euler_TroubledCellIndicator / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(TimeLim2)) &
        'Shock Detector:               ', &
        Timer_Euler_ShockDetector, ' s = ', &
        100.0_DP &
          * Timer_Euler_ShockDetector / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(TimeLim2)) &
        'Characteristic Decomposition: ', &
        Timer_Euler_CharacteristicDecomposition, ' s = ', &
        100.0_DP &
          * Timer_Euler_CharacteristicDecomposition / Timer_Euler_Program, ' %'

    END IF

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
