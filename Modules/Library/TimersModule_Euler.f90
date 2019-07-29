MODULE TimersModule_Euler

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64
  USE KindModule, Only: &
    DP, Zero

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PUBLIC :: TimeIt_Euler

  ! --- ApplicationDriver (Relativistic) ---
  REAL(DP), PUBLIC :: Timer_Euler_Initialize
  REAL(DP), PUBLIC :: Timer_Euler_Finalize
  REAL(DP), PUBLIC :: Timer_Euler_InputOutput
  REAL(DP), PUBLIC :: Timer_Euler_UpdateFluid
  REAL(DP), PUBLIC :: Timer_Euler_ComputeTimeStep
  REAL(DP), PUBLIC :: Timer_Euler_Program

  ! --- Euler_dgDiscretizationModule ---
  REAL(DP), PUBLIC :: Timer_Euler_Inc
  REAL(DP), PUBLIC :: Timer_Euler_Div_X1
  REAL(DP), PUBLIC :: Timer_Euler_Div_X2
  REAL(DP), PUBLIC :: Timer_Euler_Div_X3
  REAL(DP), PUBLIC :: Timer_Euler_Geom
  REAL(DP), PUBLIC :: Timer_Euler_Grav
  REAL(DP), PUBLIC :: Timer_Euler_CompPrim
  REAL(DP), PUBLIC :: Timer_Euler_MV       ! --- Matrix-Vector multiplies ---
  REAL(DP), PUBLIC :: Timer_Euler_RS       ! --- Riemann solvers ---

  ! --- Limiters ---
  REAL(DP), PUBLIC :: Timer_Euler_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_Euler_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_Euler_TroubledCellIndicator

  PUBLIC :: InitializeTimers_Euler
  PUBLIC :: FinalizeTimers_Euler
  PUBLIC :: TimersStart_Euler
  PUBLIC :: TimersStop_Euler
  PUBLIC :: TimersWtime_Euler

  CHARACTER(24) :: OutputFMT = '(7x,A,ES13.6E3,A,F6.3,A)'

CONTAINS


  SUBROUTINE InitializeTimers_Euler( TimeIt_Euler_Option )

    LOGICAL, INTENT(in), OPTIONAL :: TimeIt_Euler_Option

    IF( PRESENT( TimeIt_Euler_Option ) )THEN
      TimeIt_Euler = TimeIt_Euler_Option
    END IF

    IF( .NOT. TimeIt_Euler ) RETURN

    Timer_Euler_Initialize      = Zero
    Timer_Euler_Finalize        = Zero
    Timer_Euler_InputOutput     = Zero
    Timer_Euler_UpdateFluid     = Zero
    Timer_Euler_ComputeTimeStep = Zero
    Timer_Euler_Program         = Zero

    Timer_Euler_Inc    = Zero
    Timer_Euler_Div_X1 = Zero
    Timer_Euler_Div_X2 = Zero
    Timer_Euler_Div_X3 = Zero
    Timer_Euler_Geom   = Zero
    Timer_Euler_Grav   = Zero
    Timer_Euler_MV     = Zero
    Timer_Euler_RS     = Zero

    Timer_Euler_PositivityLimiter     = Zero
    Timer_Euler_SlopeLimiter          = Zero
    Timer_Euler_TroubledCellIndicator = Zero

    RETURN
  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler

    REAL(DP) :: TotalTime

    IF( .NOT. TimeIt_Euler ) RETURN

    WRITE(*,'(5x,A)') 'Timers Summary'
    WRITE(*,'(5x,A)') '--------------'
    WRITE(*,*)

    WRITE(*,'(5x,A)') 'ApplicationDriver'
    WRITE(*,'(5x,A)') '-----------------'
    WRITE(*,OutputFMT) &
      'Euler_Initialize:      ', Timer_Euler_Initialize, ' s = ', &
      100.0d0 * Timer_Euler_Initialize / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      'Euler_ComputeTimeStep: ', Timer_Euler_ComputeTimeStep, ' s = ', &
      100.0d0 * Timer_Euler_ComputeTimeStep / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      'Euler_UpdateFluid:     ', Timer_Euler_UpdateFluid, ' s = ', &
      100.0d0 * Timer_Euler_UpdateFluid / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      'Euler_Finalize:        ', Timer_Euler_Finalize, ' s = ', &
      100.0d0 * Timer_Euler_Finalize / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      'InputOutput:           ', Timer_Euler_InputOutput, ' s = ', &
      100.0d0 * Timer_Euler_InputOutput / Timer_Euler_Program, ' %'
    WRITE(*,*)

    TotalTime = Timer_Euler_Initialize + Timer_Euler_Finalize &
                  + Timer_Euler_UpdateFluid + Timer_Euler_ComputeTimeStep &
                  + Timer_Euler_Finalize + Timer_Euler_InputOutput
    WRITE(*,'(5x,A,ES14.6E3,A,F7.3,A)') &
      'Total:                  ', TotalTime, ' s = ', &
      100.0d0 * TotalTime / Timer_Euler_Program, ' %'
    WRITE(*,*)

    WRITE(*,'(5x,A)') 'Euler_dgDiscretizationModule'
    WRITE(*,'(5x,A)') '----------------------------'
    WRITE(*,OutputFMT) &
      'Increment:          ', Timer_Euler_Inc, ' s = ', &
      100.0d0 * Timer_Euler_Inc / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '  Euler_Div_X1:     ', Timer_Euler_Div_X1, ' s = ', &
      100.0d0 * Timer_Euler_Div_X1 / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '  Euler_Div_X2:     ', Timer_Euler_Div_X2, ' s = ', &
      100.0d0 * Timer_Euler_Div_X2 / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '  Euler_Div_X3:     ', Timer_Euler_Div_X3, ' s = ', &
      100.0d0 * Timer_Euler_Div_X3 / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '  Euler_Geom:       ', Timer_Euler_Geom, ' s = ', &
      100.0d0 * Timer_Euler_Geom / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '  Euler_Grav:       ', Timer_Euler_Grav, ' s = ', &
      100.0d0 * Timer_Euler_Grav / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '    Euler_CompPrim: ', Timer_Euler_CompPrim, ' s = ', &
      100.0d0 * Timer_Euler_CompPrim / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '    Euler_MV:       ', Timer_Euler_MV, ' s = ', &
      100.0d0 * Timer_Euler_MV / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '    Euler_RS:       ', Timer_Euler_RS, ' s = ', &
      100.0d0 * Timer_Euler_RS / Timer_Euler_Program, ' %'
    WRITE(*,*)

    WRITE(*,'(5x,A)') 'Limiters'
    WRITE(*,'(5x,A)') '--------'
    WRITE(*,OutputFMT) &
      'PositivityLimiter:       ', Timer_Euler_PositivityLimiter, ' s = ', &
      100.0d0 * Timer_Euler_PositivityLimiter / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      'SlopeLimiter:            ', Timer_Euler_SlopeLimiter, ' s = ', &
      100.0d0 * Timer_Euler_SlopeLimiter / Timer_Euler_Program, ' %'
    WRITE(*,OutputFMT) &
      '  TroubledCellIndicator: ', Timer_Euler_TroubledCellIndicator, ' s = ', &
      100.0d0 * Timer_Euler_TroubledCellIndicator / Timer_Euler_Program, ' %'
    WRITE(*,*)

    WRITE(*,'(7x,A,ES13.6E3,A)') &
      'Program: ', Timer_Euler_Program, ' s'
    WRITE(*,*)

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
