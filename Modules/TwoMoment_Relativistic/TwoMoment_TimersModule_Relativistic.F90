MODULE TwoMoment_TimersModule_Relativistic

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8 => INT64
  USE KindModule, Only: &
    DP, Zero, SqrtTiny

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: Timer_Total
  REAL(DP), PUBLIC :: Timer_ComputePrimative
  REAL(DP), PUBLIC :: Timer_NestedSolve

  PUBLIC :: InitializeTimers
  PUBLIC :: FinalizeTimers
  PUBLIC :: TimersStart
  PUBLIC :: TimersStop
  PUBLIC :: TimersWtime

CONTAINS


  SUBROUTINE InitializeTimers

    Timer_Total                          = Zero

    Timer_ComputePrimative               = Zero

    Timer_NestedSolve                    = Zero

    CALL TimersStart( Timer_Total )

  END SUBROUTINE InitializeTimers


  SUBROUTINE FinalizeTimers

    CALL TimersStop( Timer_Total )

    WRITE(*,*)
    WRITE(*,'(5X,A)') 'Timers Summary'
    WRITE(*,'(5X,A)') '--------------'
    WRITE(*,*)
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      'Timer_Total                              :', Timer_Total                         , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_ComputePrimative                 :', Timer_ComputePrimative              , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_NestedSolve                      :', Timer_NestedSolve                   , ' s'
    WRITE(*,*)

  END SUBROUTINE FinalizeTimers


  SUBROUTINE TimersStart( Timer )

    REAL(DP), INTENT(inout) :: Timer

    Timer = Timer - TimersWtime()

    RETURN
  END SUBROUTINE TimersStart


  SUBROUTINE TimersStop( Timer )

    REAL(DP), INTENT(inout) :: Timer

    Timer = Timer + TimersWtime()

    RETURN
  END SUBROUTINE TimersStop


  REAL(DP) FUNCTION TimersWtime()

    INTEGER(I8) :: clock_read
    INTEGER(I8) :: clock_rate
    INTEGER(I8) :: clock_max

    CALL SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
    TimersWtime = REAL( clock_read, DP ) / REAL( clock_rate, DP )

    RETURN
  END FUNCTION TimersWtime


END MODULE TwoMoment_TimersModule_Relativistic
