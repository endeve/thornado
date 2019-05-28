MODULE TimersModule_Euler

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64
  USE KindModule, Only: &
    DP, Zero

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: Timer_Euler_InputOutput
  REAL(DP), PUBLIC :: Timer_Euler_Inc   
  REAL(DP), PUBLIC :: Timer_Euler_Div_X1
  REAL(DP), PUBLIC :: Timer_Euler_Div_X2
  REAL(DP), PUBLIC :: Timer_Euler_Div_X3
  REAL(DP), PUBLIC :: Timer_Euler_Geom
  REAL(DP), PUBLIC :: Timer_Euler_Grav
  REAL(DP), PUBLIC :: Timer_Euler_CompPrim
  REAL(DP), PUBLIC :: Timer_Euler_MV

  PUBLIC :: InitializeTimers_Euler
  PUBLIC :: FinalizeTimers_Euler
  PUBLIC :: TimersStart_Euler
  PUBLIC :: TimersStop_Euler
  PUBLIC :: TimersWtime_Euler

CONTAINS


  SUBROUTINE InitializeTimers_Euler

    Timer_Euler_Inc         = Zero
    Timer_Euler_InputOutput = Zero
    Timer_Euler_Inc         = Zero
    Timer_Euler_Div_X1      = Zero
    Timer_Euler_Div_X2      = Zero
    Timer_Euler_Div_X3      = Zero
    Timer_Euler_Geom        = Zero
    Timer_Euler_Grav        = Zero
    Timer_Euler_MV          = Zero

    RETURN
  END SUBROUTINE InitializeTimers_Euler


  SUBROUTINE FinalizeTimers_Euler

    WRITE(*,'(5x,A)') 'Timers Summary'
    WRITE(*,'(5x,A)') '--------------'
    WRITE(*,*)
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      'Increment:              ', Timer_Euler_Inc,         ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      'InputOutput:            ', Timer_Euler_InputOutput, ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '    Euler_Div:          ', Timer_Euler_Inc,         ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '      Euler_Div_X1:     ', Timer_Euler_Div_X1,      ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '      Euler_Div_X2:     ', Timer_Euler_Div_X2,      ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '      Euler_Div_X3:     ', Timer_Euler_Div_X3,      ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '      Euler_Geom:       ', Timer_Euler_Geom,       ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '      Euler_Grav:       ', Timer_Euler_Grav,       ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '        Euler_CompPrim: ', Timer_Euler_CompPrim,   ' s'
    WRITE(*,'(7x,A,ES13.6E3,A)') &
      '        Euler_MV:       ', Timer_Euler_MV,         ' s'
    WRITE(*,*)

    RETURN
  END SUBROUTINE FinalizeTimers_Euler


  SUBROUTINE TimersStart_Euler( Timer )

    REAL(DP), INTENT(inout) :: Timer

    Timer = Timer - TimersWtime_Euler()

    RETURN
  END SUBROUTINE TimersStart_Euler


  SUBROUTINE TimersStop_Euler( Timer )

    REAL(DP), INTENT(inout) :: Timer

    Timer = Timer + TimersWtime_Euler()

    RETURN
  END SUBROUTINE TimersStop_Euler


  REAL(DP) FUNCTION TimersWtime_Euler()

    INTEGER(I8) :: clock_read
    INTEGER(I8) :: clock_rate
    INTEGER(I8) :: clock_max

    CALL SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
    TimersWtime_Euler = REAL( clock_read, DP ) / REAL( clock_rate, DP )

    RETURN
  END FUNCTION TimersWtime_Euler


END MODULE TimersModule_Euler
