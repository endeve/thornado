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
  REAL(DP), PUBLIC :: Timer_PL_Points
  REAL(DP), PUBLIC :: Timer_PL_CellAverage
  REAL(DP), PUBLIC :: Timer_PL_Theta_1
  REAL(DP), PUBLIC :: Timer_PL_Theta_2
  REAL(DP), PUBLIC :: Timer_PL_Out
  REAL(DP), PUBLIC :: Timer_Explicit
  REAL(DP), PUBLIC :: Timer_Ex_In
  REAL(DP), PUBLIC :: Timer_Ex_Div
  REAL(DP), PUBLIC :: Timer_Ex_Geometry
  REAL(DP), PUBLIC :: Timer_Ex_Permute
  REAL(DP), PUBLIC :: Timer_Ex_Interpolate
  REAL(DP), PUBLIC :: Timer_Ex_Flux
  REAL(DP), PUBLIC :: Timer_Ex_Increment
  REAL(DP), PUBLIC :: Timer_Ex_Out
  REAL(DP), PUBLIC :: Timer_Implicit
  REAL(DP), PUBLIC :: Timer_Im_In
  REAL(DP), PUBLIC :: Timer_Im_MapForward
  REAL(DP), PUBLIC :: Timer_Im_EosIn
  REAL(DP), PUBLIC :: Timer_Im_Solve
  REAL(DP), PUBLIC :: Timer_Im_ComputeOpacity
  REAL(DP), PUBLIC :: Timer_Im_ComputeRate
  REAL(DP), PUBLIC :: Timer_Im_ComputeLS
  REAL(DP), PUBLIC :: Timer_Im_UpdateFP
  REAL(DP), PUBLIC :: Timer_Im_CoupledAA
  REAL(DP), PUBLIC :: Timer_Im_NestedAA
  REAL(DP), PUBLIC :: Timer_Im_NestedNewton
  REAL(DP), PUBLIC :: Timer_Im_Newton
  REAL(DP), PUBLIC :: Timer_Im_NestInner
  REAL(DP), PUBLIC :: Timer_Im_EmAb_FP
  REAL(DP), PUBLIC :: Timer_Im_Presolve
  REAL(DP), PUBLIC :: Timer_Im_Increment
  REAL(DP), PUBLIC :: Timer_Im_EosOut
  REAL(DP), PUBLIC :: Timer_Im_MapBackward
  REAL(DP), PUBLIC :: Timer_Im_Out

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
    Timer_PL_Points         = Zero
    Timer_PL_CellAverage    = Zero
    Timer_PL_Theta_1        = Zero
    Timer_PL_Theta_2        = Zero
    Timer_PL_Out            = Zero
    Timer_Explicit          = Zero
    Timer_Ex_In             = Zero
    Timer_Ex_Div            = Zero
    Timer_Ex_Geometry       = Zero
    Timer_Ex_Permute        = Zero
    Timer_Ex_Interpolate    = Zero
    Timer_Ex_Flux           = Zero
    Timer_Ex_Increment      = Zero
    Timer_Ex_Out            = Zero
    Timer_Implicit          = Zero
    Timer_Im_In             = Zero
    Timer_Im_MapForward     = Zero
    Timer_Im_EosIn          = Zero
    Timer_Im_Solve          = Zero
    Timer_Im_ComputeOpacity = Zero
    Timer_Im_ComputeRate    = Zero
    Timer_Im_ComputeLS      = Zero
    Timer_Im_UpdateFP       = Zero
    Timer_Im_CoupledAA      = Zero
    Timer_Im_NestedAA       = Zero
    Timer_Im_NestedNewton   = Zero
    Timer_Im_Newton         = Zero
    Timer_Im_NestInner      = Zero
    Timer_Im_EmAb_FP        = Zero
    Timer_Im_Presolve       = Zero
    Timer_Im_Increment      = Zero
    Timer_Im_EosOut         = Zero
    Timer_Im_MapBackward    = Zero
    Timer_Im_Out            = Zero

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
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Points           :', Timer_PL_Points        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_CellAverage      :', Timer_PL_CellAverage   , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Theta_1          :', Timer_PL_Theta_1       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Theta_2          :', Timer_PL_Theta_2       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    PL_Out              :', Timer_PL_Out           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  Explicit              :', Timer_Explicit         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_In               :', Timer_Ex_In            , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Div              :', Timer_Ex_Div           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Geometry         :', Timer_Ex_Geometry      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Permute          :', Timer_Ex_Permute       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Interpolate      :', Timer_Ex_Interpolate   , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Flux             :', Timer_Ex_Flux          , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Increment        :', Timer_Ex_Increment     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Ex_Out              :', Timer_Ex_Out           , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '  Implicit              :', Timer_Implicit         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_In               :', Timer_Im_In            , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_MapForward       :', Timer_Im_MapForward    , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_EosIn            :', Timer_Im_EosIn         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_Solve            :', Timer_Im_Solve         , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_ComputeOpacity :', Timer_Im_ComputeOpacity, ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_ComputeRate    :', Timer_Im_ComputeRate   , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_ComputeLS      :', Timer_Im_ComputeLS     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '      Im_UpdateFP       :', Timer_Im_UpdateFP      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    CoupledAA           :', Timer_Im_CoupledAA     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    NestedAA            :', Timer_Im_NestedAA      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    NestedNewton        :', Timer_Im_NestedNewton  , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Newton              :', Timer_Im_Newton        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Nested Inner loop   :', Timer_Im_NestInner     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    EmAb_FP precond     :', Timer_Im_EmAb_FP       , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_Pair_Presolve    :', Timer_Im_Presolve      , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_Increment        :', Timer_Im_Increment     , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_EosOut           :', Timer_Im_EosOut        , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_MapBackward      :', Timer_Im_MapBackward   , ' s'
    WRITE(*,'(7X,A,5x,ES12.6E2,A)') '    Im_Out              :', Timer_Im_Out           , ' s'
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
