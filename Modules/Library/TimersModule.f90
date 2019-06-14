MODULE TimersModule

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64
  USE KindModule, Only: &
    DP, Zero

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: Timer_Initialize
  REAL(DP), PUBLIC :: Timer_InputOutput
  REAL(DP), PUBLIC :: Timer_Evolve
  REAL(DP), PUBLIC :: Timer_AddFieldsF
  REAL(DP), PUBLIC :: Timer_AddFieldsR
  REAL(DP), PUBLIC :: Timer_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_PL_In
  REAL(DP), PUBLIC :: Timer_PL_P
  REAL(DP), PUBLIC :: Timer_PL_K
  REAL(DP), PUBLIC :: Timer_PL_Theta_1
  REAL(DP), PUBLIC :: Timer_PL_Theta_2
  REAL(DP), PUBLIC :: Timer_PL_Out
  REAL(DP), PUBLIC :: Timer_Explicit
  REAL(DP), PUBLIC :: Timer_Ex_In
  REAL(DP), PUBLIC :: Timer_Ex_Div
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_In
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_G
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_U
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_S
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_V
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_dU
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_Out
  REAL(DP), PUBLIC :: Timer_Ex_Div_X1_MM
  REAL(DP), PUBLIC :: Timer_Ex_Div_X2
  REAL(DP), PUBLIC :: Timer_Ex_Div_X3
  REAL(DP), PUBLIC :: Timer_Ex_Out
  REAL(DP), PUBLIC :: Timer_Implicit
  REAL(DP), PUBLIC :: Timer_Im_In
  REAL(DP), PUBLIC :: Timer_Im_ComputeTS_Aux
  REAL(DP), PUBLIC :: Timer_Im_ComputeOpacity
  REAL(DP), PUBLIC :: Timer_Im_MapForward
  REAL(DP), PUBLIC :: Timer_Im_Solve
  REAL(DP), PUBLIC :: Timer_Im_CoupledAA
  REAL(DP), PUBLIC :: Timer_Im_NestedAA
  REAL(DP), PUBLIC :: Timer_Im_NestedNewton
  REAL(DP), PUBLIC :: Timer_Im_Newton
  REAL(DP), PUBLIC :: Timer_Im_Out
  REAL(DP), PUBLIC :: Timer_Im_ComputeTS_Prim
  REAL(DP), PUBLIC :: Timer_Im_Increment
  REAL(DP), PUBLIC :: Timer_Im_MapBackward

  PUBLIC :: InitializeTimers
  PUBLIC :: FinalizeTimers
  PUBLIC :: TimersStart
  PUBLIC :: TimersStop
  PUBLIC :: TimersWtime

CONTAINS


  SUBROUTINE InitializeTimers

    Timer_Initialize        = Zero
    Timer_InputOutput       = Zero
    Timer_Evolve            = Zero
    Timer_AddFieldsF        = Zero
    Timer_AddFieldsR        = Zero
    Timer_PositivityLimiter = Zero
    Timer_PL_In             = Zero
    Timer_PL_P              = Zero
    Timer_PL_K              = Zero
    Timer_PL_Theta_1        = Zero
    Timer_PL_Theta_2        = Zero
    Timer_PL_Out            = Zero
    Timer_Explicit          = Zero
    Timer_Ex_In             = Zero
    Timer_Ex_Div            = Zero
    Timer_Ex_Div_X1         = Zero
    Timer_Ex_Div_X1_In      = Zero
    Timer_Ex_Div_X1_G       = Zero
    Timer_Ex_Div_X1_U       = Zero
    Timer_Ex_Div_X1_S       = Zero
    Timer_Ex_Div_X1_V       = Zero
    Timer_Ex_Div_X1_dU      = Zero
    Timer_Ex_Div_X1_Out     = Zero
    Timer_Ex_Div_X1_MM      = Zero
    Timer_Ex_Div_X2         = Zero
    Timer_Ex_Div_X3         = Zero
    Timer_Ex_Out            = Zero
    Timer_Implicit          = Zero
    Timer_Im_In             = Zero
    Timer_Im_ComputeTS_Aux  = Zero
    Timer_Im_ComputeOpacity = Zero
    Timer_Im_MapForward     = Zero
    Timer_Im_Solve          = Zero
    Timer_Im_CoupledAA      = Zero
    Timer_Im_NestedAA       = Zero
    Timer_Im_NestedNewton   = Zero
    Timer_Im_Newton         = Zero
    Timer_Im_Out            = Zero
    Timer_Im_ComputeTS_Prim = Zero
    Timer_Im_Increment      = Zero
    Timer_Im_MapBackward    = Zero

    RETURN
  END SUBROUTINE InitializeTimers


  SUBROUTINE FinalizeTimers

    WRITE(*,'(5X,A)') 'Timers Summary'
    WRITE(*,'(5X,A)') '--------------'
    WRITE(*,*)
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') 'Initialize              :', Timer_Initialize       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') 'InputOutput             :', Timer_InputOutput      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') 'Evolve                  :', Timer_Evolve           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  AddFieldsF            :', Timer_AddFieldsF       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  AddFieldsR            :', Timer_AddFieldsR       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  PositivityLimiter     :', Timer_PositivityLimiter, ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_In               :', Timer_PL_In            , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_P                :', Timer_PL_P             , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_K                :', Timer_PL_K             , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Theta_1          :', Timer_PL_Theta_1       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Theta_2          :', Timer_PL_Theta_2       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Out              :', Timer_PL_Out           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  Explicit              :', Timer_Explicit         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_In               :', Timer_Ex_In            , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Div              :', Timer_Ex_Div           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Ex_Div_X1         :', Timer_Ex_Div_X1        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_In    :', Timer_Ex_Div_X1_In     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_G     :', Timer_Ex_Div_X1_G      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_U     :', Timer_Ex_Div_X1_U      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_S     :', Timer_Ex_Div_X1_S      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_V     :', Timer_Ex_Div_X1_V      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_dU    :', Timer_Ex_Div_X1_dU     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_Out   :', Timer_Ex_Div_X1_Out    , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '        Ex_Div_X1_MM    :', Timer_Ex_Div_X1_MM     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Ex_Div_X2         :', Timer_Ex_Div_X2        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Ex_Div_X3         :', Timer_Ex_Div_X3        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Out              :', Timer_Ex_Out           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  Implicit              :', Timer_Implicit         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_In               :', Timer_Im_In            , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_ComputeTS_Aux  :', Timer_Im_ComputeTS_Aux , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_ComputeOpacity :', Timer_Im_ComputeOpacity, ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_MapForward     :', Timer_Im_MapForward    , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_Solve            :', Timer_Im_Solve         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    CoupledAA           :', Timer_Im_CoupledAA     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    NestedAA            :', Timer_Im_NestedAA      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    NestedNewton        :', Timer_Im_NestedNewton  , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Newton              :', Timer_Im_Newton        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_Out              :', Timer_Im_Out           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_ComputeTS_Prim :', Timer_Im_ComputeTS_Prim, ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_Increment      :', Timer_Im_Increment     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_MapBackward    :', Timer_Im_MapBackward   , ' s'
    WRITE(*,*)

    RETURN
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


END MODULE TimersModule
