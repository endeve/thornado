MODULE TwoMoment_TimersModule_OrderV

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8 => INT64
  USE KindModule, Only: &
    DP, Zero, SqrtTiny

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: Timer_Total

  REAL(DP), PUBLIC :: Timer_Streaming
  REAL(DP), PUBLIC :: Timer_Streaming_Permute
  REAL(DP), PUBLIC :: Timer_Streaming_LinearAlgebra
  REAL(DP), PUBLIC :: Timer_Streaming_BCs
  REAL(DP), PUBLIC :: Timer_Streaming_Zero
  REAL(DP), PUBLIC :: Timer_Streaming_Divergence_X1
  REAL(DP), PUBLIC :: Timer_Streaming_Divergence_X2
  REAL(DP), PUBLIC :: Timer_Streaming_Divergence_X3
  REAL(DP), PUBLIC :: Timer_Streaming_ObserverCorrections
  REAL(DP), PUBLIC :: Timer_Streaming_Derivatives_X1
  REAL(DP), PUBLIC :: Timer_Streaming_Derivatives_X2
  REAL(DP), PUBLIC :: Timer_Streaming_Derivatives_X3
  REAL(DP), PUBLIC :: Timer_Streaming_InverseMassMatrix
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_InOut
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_RHS
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_LS
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_Update
  REAL(DP), PUBLIC :: Timer_Streaming_Sources
  REAL(DP), PUBLIC :: Timer_Collisions
  REAL(DP), PUBLIC :: Timer_Collisions_Permute
  REAL(DP), PUBLIC :: Timer_Collisions_PrimitiveFluid
  REAL(DP), PUBLIC :: Timer_Collisions_Solve
  REAL(DP), PUBLIC :: Timer_TCI
  REAL(DP), PUBLIC :: Timer_SlopeLimiter
  REAL(DP), PUBLIC :: Timer_PositivityLimiter
  REAL(DP), PUBLIC :: Timer_TimeStepper

  PUBLIC :: InitializeTimers
  PUBLIC :: FinalizeTimers
  PUBLIC :: TimersStart
  PUBLIC :: TimersStop
  PUBLIC :: TimersWtime

CONTAINS


  SUBROUTINE InitializeTimers

    Timer_Total                         = Zero

    Timer_Streaming                     = Zero
    Timer_Streaming_Permute             = Zero
    Timer_Streaming_LinearAlgebra       = Zero
    Timer_Streaming_BCs                 = Zero
    Timer_Streaming_Zero                = Zero
    Timer_Streaming_Divergence_X1       = Zero
    Timer_Streaming_Divergence_X2       = Zero
    Timer_Streaming_Divergence_X3       = Zero
    Timer_Streaming_ObserverCorrections = Zero
    Timer_Streaming_Derivatives_X1      = Zero
    Timer_Streaming_Derivatives_X2      = Zero
    Timer_Streaming_Derivatives_X3      = Zero
    Timer_Streaming_InverseMassMatrix   = Zero
    Timer_Streaming_NumericalFlux       = Zero
    Timer_Streaming_NumericalFlux_InOut = Zero
    Timer_Streaming_NumericalFlux_RHS   = Zero
    Timer_Streaming_NumericalFlux_LS    = Zero
    Timer_Streaming_NumericalFlux_Update= Zero
    Timer_Streaming_Sources             = Zero

    Timer_Collisions                    = Zero
    Timer_Collisions_Permute            = Zero
    Timer_Collisions_PrimitiveFluid     = Zero
    Timer_Collisions_Solve              = Zero

    Timer_TCI                           = Zero

    Timer_SlopeLimiter                  = Zero

    Timer_PositivityLimiter             = Zero

    Timer_TimeStepper                   = Zero

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
      '  Timer_Streaming                        :', Timer_Streaming                     , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_BCs                  :', Timer_Streaming_BCs                 , ' s ', &
      Timer_Streaming_BCs / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_Zero                 :', Timer_Streaming_Zero                , ' s ', &
      Timer_Streaming_Zero / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_Divergence_X1        :', Timer_Streaming_Divergence_X1       , ' s ', &
      Timer_Streaming_Divergence_X1 / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_Divergence_X2        :', Timer_Streaming_Divergence_X2       , ' s ', &
      Timer_Streaming_Divergence_X2 / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_Divergence_X3        :', Timer_Streaming_Divergence_X3       , ' s ', &
      Timer_Streaming_Divergence_X3 / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_ObserverCorrections  :', Timer_Streaming_ObserverCorrections , ' s ', &
      Timer_Streaming_ObserverCorrections / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_InverseMassMatrix    :', Timer_Streaming_InverseMassMatrix   , ' s ', &
      Timer_Streaming_InverseMassMatrix / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Derivatives_X1       :', Timer_Streaming_Derivatives_X1      , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Derivatives_X2       :', Timer_Streaming_Derivatives_X2      , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Derivatives_X3       :', Timer_Streaming_Derivatives_X3      , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_NumericalFlux        :', Timer_Streaming_NumericalFlux       , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_NumericalFlux_InOut  :', Timer_Streaming_NumericalFlux_InOut , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_NumericalFlux_RHS    :', Timer_Streaming_NumericalFlux_RHS   , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_NumericalFlux_LS     :', Timer_Streaming_NumericalFlux_LS    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_NumericalFlux_Update :', Timer_Streaming_NumericalFlux_Update, ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Sources              :', Timer_Streaming_Sources             , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Permute              :', Timer_Streaming_Permute             , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_LinearAlgebra        :', Timer_Streaming_LinearAlgebra       , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_Collisions                       :', Timer_Collisions                    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_Collisions_Permute               :', Timer_Collisions_Permute            , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_Collisions_PrimitiveFluid        :', Timer_Collisions_PrimitiveFluid     , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_Collisions_Solve                 :', Timer_Collisions_Solve              , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_TCI                              :', Timer_TCI                           , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_SlopeLimiter                     :', Timer_SlopeLimiter                  , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_PositivityLimiter                :', Timer_PositivityLimiter             , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_TimeStepper                      :', Timer_TimeStepper                   , ' s'

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


END MODULE TwoMoment_TimersModule_OrderV
