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

  ! --- DG discretization-specific ---
  REAL(DP), PUBLIC :: Timer_Euler_dgDiscretization
  REAL(DP), PUBLIC :: Timer_Euler_Divergence
  REAL(DP), PUBLIC :: Timer_Euler_Geometry
  REAL(DP), PUBLIC :: Timer_Euler_Gravity
  REAL(DP), PUBLIC :: Timer_Euler_MV     ! --- Matrix-Vector multiplies ---
  REAL(DP), PUBLIC :: Timer_Euler_RS     ! --- Riemann solvers ---
  REAL(DP), PUBLIC :: Timer_Euler_ComputePrimitive

  ! --- Limiter-specific ---
  REAL(DP), PUBLIC :: Timer_Euler_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_Euler_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_Euler_CharacteristicDecomposition
  REAL(DP), PUBLIC :: Timer_Euler_TroubledCellIndicator

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
    Timer_Euler_MV               = Zero
    Timer_Euler_RS               = Zero
    Timer_Euler_ComputePrimitive = Zero

    Timer_Euler_PositivityLimiter           = Zero
    Timer_Euler_SlopeLimiter                = Zero
    Timer_Euler_CharacteristicDecomposition = Zero
    Timer_Euler_TroubledCellIndicator       = Zero

    CALL TimersStart_Euler( Timer_Euler_Program )

    RETURN
  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler &
    ( Verbose_Option, SuppressApplicationDriver_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL, INTENT(in), OPTIONAL :: SuppressApplicationDriver_Option

    LOGICAL  :: Verbose, SuppressApplicationDriver
    REAL(DP) :: TotalTime

    CHARACTER(64) :: OutputFMT1 = '(14x,A19,ES13.6E3,A,F6.3,A)'
    CHARACTER(64) :: OutputFMT2 = '(14x,A24,ES13.6E3,A,F6.3,A,F6.3,A)'
    CHARACTER(64) :: OutputFMT3 = '(14x,A20,ES13.6E3,A,F6.3,A)'
    CHARACTER(64) :: OutputFMT4 = '(16x,A30,ES13.6E3,A,F6.3,A)'

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
      WRITE(*,'(6x,A)') 'Timers Summary (Euler)'
      WRITE(*,'(6x,A)') '----------------------'
      WRITE(*,*)

      WRITE(*,'(8x,A,ES13.6E3,A)') &
        'Total run-time = ', Timer_Euler_Program, ' s'

    END IF

    IF( .NOT. SuppressApplicationDriver )THEN

      WRITE(*,*)
      WRITE(*,'(10x,A)') 'Application Driver'
      WRITE(*,'(10x,A)') '------------------'
      WRITE(*,*)
      TotalTime = Timer_Euler_Initialize &
                    + Timer_Euler_ComputeTimeStep &
                    + Timer_Euler_UpdateFluid &
                    + Timer_Euler_InputOutput &
                    + Timer_Euler_Finalize

      WRITE(*,'(12x,A,ES13.6E3,A,F7.3,A)') &
        'Timers = ', TotalTime, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(OutputFMT1)) &
        'Initialize: ', &
        Timer_Euler_Initialize, ' s = ', &
        100.0_DP &
        * Timer_Euler_Initialize / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT1)) &
        'Compute Time-Step: ', &
        Timer_Euler_ComputeTimeStep, ' s = ', &
        100.0_DP &
        * Timer_Euler_ComputeTimeStep / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT1)) &
        'Update Fluid: ', &
        Timer_Euler_UpdateFluid, ' s = ', &
        100.0_DP &
          * Timer_Euler_UpdateFluid / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT1)) &
        'Input/Output: ', &
        Timer_Euler_InputOutput, ' s = ', &
        100.0_DP &
          * Timer_Euler_InputOutput / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT1)) &
        'Finalize: ', &
        Timer_Euler_Finalize, ' s = ', &
        100.0_DP &
          * Timer_Euler_Finalize / Timer_Euler_Program, ' %'

    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(10x,A)') 'DG discretization-specific'
      WRITE(*,'(10x,A)') '--------------------------'
      WRITE(*,*)
      TotalTime = Timer_Euler_Divergence &
                    + Timer_Euler_Geometry &
                    + Timer_Euler_Gravity

      WRITE(*,'(51x,A)') '%DG'
      WRITE(*,'(12x,A,ES13.6E3,A,F6.3,A,F7.3,A)') &
        'Timers-DG = ', Timer_Euler_dgDiscretization, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' % === ', &
        100.0_DP * TotalTime / Timer_Euler_dgDiscretization
      WRITE(*,*)

      WRITE(*,'(65x,A)') '%DG'
      WRITE(*,TRIM(OutputFMT2)) &
        'Divergence: ', &
        Timer_Euler_Divergence, ' s = ', &
        100.0_DP &
          * Timer_Euler_Divergence / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_Divergence / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,TRIM(OutputFMT2)) &
        'Geometry: ', &
        Timer_Euler_Geometry, ' s = ', &
        100.0_DP &
          * Timer_Euler_Geometry / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_Geometry / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,TRIM(OutputFMT2)) &
        'Gravity: ', &
        Timer_Euler_Gravity, ' s = ', &
        100.0_DP &
          * Timer_Euler_Gravity / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_Gravity / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,*)
      WRITE(*,TRIM(OutputFMT2)) &
        'Compute Primitive: ', &
        Timer_Euler_ComputePrimitive, ' s = ', &
        100.0_DP &
          * Timer_Euler_ComputePrimitive / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_ComputePrimitive / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,TRIM(OutputFMT2)) &
        'Matrix-Vector Multiply: ', &
        Timer_Euler_MV, ' s = ', &
        100.0_DP &
          * Timer_Euler_MV / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_MV / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,TRIM(OutputFMT2)) &
        'Riemann Solver: ', &
        Timer_Euler_RS, ' s = ', &
        100.0_DP &
          * Timer_Euler_RS / Timer_Euler_Program, ' % === ', &
        100.0_DP &
          * Timer_Euler_RS / Timer_Euler_dgDiscretization, ' %'

      WRITE(*,*)
      WRITE(*,'(10x,A)') 'Limiter-specific'
      WRITE(*,'(10x,A)') '----------------'
      WRITE(*,*)
      TotalTime = Timer_Euler_PositivityLimiter &
                  +  Timer_Euler_SlopeLimiter

      WRITE(*,'(12x,A,ES13.6E3,A,F6.3,A)') &
        'Timers-Limiters = ', TotalTime, ' s = ', &
        100.0_DP * TotalTime / Timer_Euler_Program, ' %'
      WRITE(*,*)

      WRITE(*,TRIM(OutputFMT3)) &
        'Positivity-Limiter: ', &
        Timer_Euler_PositivityLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_PositivityLimiter / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT3)) &
        'Slope-Limiter: ', &
        Timer_Euler_SlopeLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_SlopeLimiter / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT4)) &
        'Troubled-Cell Indicator: ', &
        Timer_Euler_TroubledCellIndicator, ' s = ', &
        100.0_DP &
          * Timer_Euler_TroubledCellIndicator / Timer_Euler_Program, ' %'

      WRITE(*,TRIM(OutputFMT4)) &
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
