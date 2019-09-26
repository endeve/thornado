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
  REAL(DP), PUBLIC :: Timer_Euler_Inc
  REAL(DP), PUBLIC :: Timer_Euler_MV     ! --- Matrix-Vector multiplies ---
  REAL(DP), PUBLIC :: Timer_Euler_RS     ! --- Riemann solvers ---
  REAL(DP), PUBLIC :: Timer_Euler_Div_X1
  REAL(DP), PUBLIC :: Timer_Euler_Div_X2
  REAL(DP), PUBLIC :: Timer_Euler_Div_X3
  REAL(DP), PUBLIC :: Timer_Euler_Geom
  REAL(DP), PUBLIC :: Timer_Euler_Grav
  REAL(DP), PUBLIC :: Timer_Euler_CompPrim

  ! --- Limiter-specific ---
  REAL(DP), PUBLIC :: Timer_Euler_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_Euler_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_Euler_TroubledCellIndicator

  CHARACTER(24) :: OutputFMT = '(7x,A,ES13.6E3,A,F6.3,A)'

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

    Timer_Euler_Inc      = Zero
    Timer_Euler_MV       = Zero
    Timer_Euler_RS       = Zero
    Timer_Euler_Div_X1   = Zero
    Timer_Euler_Div_X2   = Zero
    Timer_Euler_Div_X3   = Zero
    Timer_Euler_Geom     = Zero
    Timer_Euler_Grav     = Zero
    Timer_Euler_CompPrim = Zero

    Timer_Euler_PositivityLimiter     = Zero
    Timer_Euler_SlopeLimiter          = Zero
    Timer_Euler_TroubledCellIndicator = Zero

    CALL TimersStart_Euler( Timer_Euler_Program )

    RETURN
  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler &
    ( Verbose_Option, SuppressApplicationDriver_Option )

    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL, INTENT(in), OPTIONAL :: SuppressApplicationDriver_Option

    LOGICAL  :: Verbose, SuppressApplicationDriver
    REAL(DP) :: TotalTime

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
      WRITE(*,'(5x,A)') 'Timers (Euler) Summary'
      WRITE(*,'(5x,A)') '---------------------'
      WRITE(*,*)

      WRITE(*,'(7x,A,ES13.6E3,A)') &
        'TotalRunTime = ', Timer_Euler_Program, ' s'

    END IF

    IF( .NOT. SuppressApplicationDriver )THEN


      WRITE(*,*)
      WRITE(*,'(5x,A)') '  ApplicationDriver'
      WRITE(*,'(5x,A)') '  -----------------'
      WRITE(*,*)
      TotalTime = Zero

      WRITE(*,OutputFMT) &
        '  Initialize:        ', &
        Timer_Euler_Initialize, ' s = ', &
        100.0_DP &
        * Timer_Euler_Initialize / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Initialize

      WRITE(*,OutputFMT) &
        '  Compute Time-Step: ', &
        Timer_Euler_ComputeTimeStep, ' s = ', &
        100.0_DP &
        * Timer_Euler_ComputeTimeStep / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_ComputeTimeStep

      WRITE(*,OutputFMT) &
        '  Update Fluid:      ', &
        Timer_Euler_UpdateFluid, ' s = ', &
        100.0_DP &
          * Timer_Euler_UpdateFluid / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_UpdateFluid

      WRITE(*,OutputFMT) &
        '  Input/Output:      ', &
        Timer_Euler_InputOutput, ' s = ', &
        100.0_DP &
          * Timer_Euler_InputOutput / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_InputOutput

      WRITE(*,OutputFMT) &
        '  Finalize:          ', &
        Timer_Euler_Finalize, ' s = ', &
        100.0_DP &
          * Timer_Euler_Finalize / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Finalize

      WRITE(*,*)
      WRITE(*,'(7x,A,ES13.6E3,A)') &
        '  Timers = ', TotalTime, ' s'
      WRITE(*,*)
      WRITE(*,'(7x,A,F7.3,A)') &
        '  Timers / TotalRunTime = ', &
        100.0_DP &
        * TotalTime / Timer_Euler_Program, ' %'

    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(7x,A)') '  DG discretization-specific'
      WRITE(*,'(7x,A)') '  --------------------------'
      WRITE(*,*)
      TotalTime = Zero

      WRITE(*,OutputFMT) &
        '    Increment:  ', &
        Timer_Euler_Inc, ' s = ', &
        100.0_DP &
          * Timer_Euler_Inc / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Inc

      WRITE(*,OutputFMT) &
        '      Div X1:   ', &
        Timer_Euler_Div_X1, ' s = ', &
        100.0_DP &
          * Timer_Euler_Div_X1 / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Div_X1

      WRITE(*,OutputFMT) &
        '      Div X2:   ', &
        Timer_Euler_Div_X2, ' s = ', &
        100.0_DP &
          * Timer_Euler_Div_X2 / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Div_X2

      WRITE(*,OutputFMT) &
        '      Div X3:   ', &
        Timer_Euler_Div_X3, ' s = ', &
        100.0_DP &
          * Timer_Euler_Div_X3 / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Div_X3

      WRITE(*,OutputFMT) &
        '      Geometry: ', &
        Timer_Euler_Geom, ' s = ', &
        100.0_DP &
          * Timer_Euler_Geom / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Geom

      WRITE(*,OutputFMT) &
        '      Gravity:  ', &
        Timer_Euler_Grav, ' s = ', &
        100.0_DP &
          * Timer_Euler_Grav / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_Grav

      WRITE(*,*)
      WRITE(*,'(7x,A)') '  Auxiliary'
      WRITE(*,'(7x,A)') '  ---------'
      WRITE(*,*)
      WRITE(*,OutputFMT) &
        '      Compute Primitive:      ', &
        Timer_Euler_CompPrim, ' s = ', &
        100.0_DP &
          * Timer_Euler_CompPrim / Timer_Euler_Program, ' %'

      WRITE(*,OutputFMT) &
        '      Matrix-Vector Multiply: ', &
        Timer_Euler_MV, ' s = ', &
        100.0_DP &
          * Timer_Euler_MV / Timer_Euler_Program, ' %'

      WRITE(*,OutputFMT) &
        '      Riemann Solver:         ', &
        Timer_Euler_RS, ' s = ', &
        100.0_DP &
          * Timer_Euler_RS / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,'(7x,A)') '  Limiter-specific'
      WRITE(*,'(7x,A)') '  ----------------'
      WRITE(*,*)
      TotalTime = Zero

      WRITE(*,OutputFMT) &
        '    Positivity-Limiter:        ', &
        Timer_Euler_PositivityLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_PositivityLimiter / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_PositivityLimiter

      WRITE(*,OutputFMT) &
        '    Slope-Limiter:             ', &
        Timer_Euler_SlopeLimiter, ' s = ', &
        100.0_DP &
          * Timer_Euler_SlopeLimiter / Timer_Euler_Program, ' %'
      TotalTime = TotalTime + Timer_Euler_SlopeLimiter

      WRITE(*,OutputFMT) &
        '      Troubled-Cell Indicator: ', &
        Timer_Euler_TroubledCellIndicator, ' s = ', &
        100.0_DP &
          * Timer_Euler_TroubledCellIndicator / Timer_Euler_Program, ' %'

      WRITE(*,*)
      WRITE(*,'(7x,A,ES13.6E3,A)') &
        '    Timers-Limiters = ', TotalTime, ' s'
      WRITE(*,*)
      WRITE(*,'(7x,A,F6.3,A)') &
        '    Timers-Limiters / TotalRunTime = ', &
        100.0_DP &
          * TotalTime / Timer_Euler_Program, ' %'

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
