MODULE TwoMoment_TimersModule_OrderV

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8 => INT64
  USE KindModule, Only: &
    DP, Zero, SqrtTiny

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PUBLIC :: Timer_Total
  REAL(DP), PUBLIC :: Timer_IMEX
  REAL(DP), PUBLIC :: Timer_Streaming
  REAL(DP), PUBLIC :: Timer_Streaming_LinearAlgebra
  REAL(DP), PUBLIC :: Timer_Streaming_Divergence
  REAL(DP), PUBLIC :: Timer_Streaming_ObserverCorrections
  REAL(DP), PUBLIC :: Timer_Streaming_Derivatives
  REAL(DP), PUBLIC :: Timer_Streaming_Eigenvalues
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_InOut
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_RHS
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_LS
  REAL(DP), PUBLIC :: Timer_Streaming_NumericalFlux_Update
  REAL(DP), PUBLIC :: Timer_Streaming_PrimitiveTwoMoment
  REAL(DP), PUBLIC :: Timer_Streaming_Sources
  REAL(DP), PUBLIC :: Timer_Collisions
  REAL(DP), PUBLIC :: Timer_Collisions_PrimitiveFluid
  REAL(DP), PUBLIC :: Timer_Collisions_PrimitiveTwoMoment
  REAL(DP), PUBLIC :: Timer_Collisions_Solve
  REAL(DP), PUBLIC :: Timer_Collisions_OuterLoop
  REAL(DP), PUBLIC :: Timer_Collisions_InnerLoop
  REAL(DP), PUBLIC :: Timer_Collisions_ComputeOpacity
  REAL(DP), PUBLIC :: Timer_Collisions_ComputeRates
  REAL(DP), PUBLIC :: Timer_Collisions_InitializeRHS
  REAL(DP), PUBLIC :: Timer_Collisions_NeutrinoRHS
  REAL(DP), PUBLIC :: Timer_Collisions_MatterRHS
  REAL(DP), PUBLIC :: Timer_Collisions_SolveLS
  REAL(DP), PUBLIC :: Timer_Collisions_UpdateFP
  REAL(DP), PUBLIC :: Timer_Collisions_CheckOuter
  REAL(DP), PUBLIC :: Timer_Collisions_CheckInner
  REAL(DP), PUBLIC :: Timer_TCI
  REAL(DP), PUBLIC :: Timer_TCI_Permute
  REAL(DP), PUBLIC :: Timer_TCI_LinearAlgebra
  REAL(DP), PUBLIC :: Timer_TCI_Compute
  REAL(DP), PUBLIC :: Timer_SL
  REAL(DP), PUBLIC :: Timer_SL_Permute
  REAL(DP), PUBLIC :: Timer_SL_LinearAlgebra
  REAL(DP), PUBLIC :: Timer_SL_MinMod
  REAL(DP), PUBLIC :: Timer_SL_ReplaceSlopes
  REAL(DP), PUBLIC :: Timer_SL_Correction
  REAL(DP), PUBLIC :: Timer_PL
  REAL(DP), PUBLIC :: Timer_PL_Permute
  REAL(DP), PUBLIC :: Timer_PL_PointValues
  REAL(DP), PUBLIC :: Timer_PL_CellAverage
  REAL(DP), PUBLIC :: Timer_PL_Theta_1
  REAL(DP), PUBLIC :: Timer_PL_Theta_2
  REAL(DP), PUBLIC :: Timer_PL_EnergyLimiter
  REAL(DP), PUBLIC :: Timer_TimeStepper

  PUBLIC :: InitializeTimers
  PUBLIC :: FinalizeTimers
  PUBLIC :: TimersStart
  PUBLIC :: TimersStop
  PUBLIC :: TimersWtime

CONTAINS


  SUBROUTINE InitializeTimers

    Timer_Total                          = Zero
    Timer_IMEX                           = Zero

    Timer_Streaming                      = Zero
    Timer_Streaming_LinearAlgebra        = Zero
    Timer_Streaming_Divergence           = Zero
    Timer_Streaming_ObserverCorrections  = Zero
    Timer_Streaming_Derivatives          = Zero
    Timer_Streaming_Eigenvalues          = Zero
    Timer_Streaming_NumericalFlux        = Zero
    Timer_Streaming_NumericalFlux_InOut  = Zero
    Timer_Streaming_NumericalFlux_RHS    = Zero
    Timer_Streaming_NumericalFlux_LS     = Zero
    Timer_Streaming_NumericalFlux_Update = Zero
    Timer_Streaming_PrimitiveTwoMoment   = Zero
    Timer_Streaming_Sources              = Zero

    Timer_Collisions                     = Zero
    Timer_Collisions_PrimitiveFluid      = Zero
    Timer_Collisions_PrimitiveTwoMoment  = Zero
    Timer_Collisions_Solve               = Zero
    Timer_Collisions_OuterLoop           = Zero
    Timer_Collisions_InnerLoop           = Zero
    Timer_Collisions_ComputeOpacity      = Zero
    Timer_Collisions_ComputeRates        = Zero
    Timer_Collisions_InitializeRHS       = Zero
    Timer_Collisions_NeutrinoRHS         = Zero
    Timer_Collisions_MatterRHS           = Zero
    Timer_Collisions_SolveLS             = Zero
    Timer_Collisions_UpdateFP            = Zero
    Timer_Collisions_CheckOuter          = Zero
    Timer_Collisions_CheckInner          = Zero

    Timer_TCI                            = Zero
    Timer_TCI_Permute                    = Zero
    Timer_TCI_LinearAlgebra              = Zero
    Timer_TCI_Compute                    = Zero

    Timer_SL                             = Zero
    Timer_SL_Permute                     = Zero
    Timer_SL_LinearAlgebra               = Zero
    Timer_SL_MinMod                      = Zero
    Timer_SL_ReplaceSlopes               = Zero
    Timer_SL_Correction                  = Zero

    Timer_PL                             = Zero
    Timer_PL_Permute                     = Zero
    Timer_PL_PointValues                 = Zero
    Timer_PL_CellAverage                 = Zero
    Timer_PL_Theta_1                     = Zero
    Timer_PL_Theta_2                     = Zero
    Timer_PL_EnergyLimiter               = Zero

    Timer_TimeStepper                    = Zero

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
      '  Timer_IMEX                             :', Timer_IMEX                          , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_Streaming                        :', Timer_Streaming                     , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_Divergence           :', Timer_Streaming_Divergence          , ' s ', &
      Timer_Streaming_Divergence / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A,ES12.6E2)') &
      '    Timer_Streaming_ObserverCorrections  :', Timer_Streaming_ObserverCorrections , ' s ', &
      Timer_Streaming_ObserverCorrections / MAX( Timer_Streaming, SqrtTiny )
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Derivatives          :', Timer_Streaming_Derivatives         , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Eigenvalues          :', Timer_Streaming_Eigenvalues         , ' s'
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
      '    Timer_Streaming_PrimitiveTwoMoment   :', Timer_Streaming_PrimitiveTwoMoment  , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_Sources              :', Timer_Streaming_Sources             , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Streaming_LinearAlgebra        :', Timer_Streaming_LinearAlgebra       , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_Collisions                       :', Timer_Collisions                    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_PrimitiveFluid      :', Timer_Collisions_PrimitiveFluid     , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_PrimitiveTwoMoment  :', Timer_Collisions_PrimitiveTwoMoment , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_Solve               :', Timer_Collisions_Solve              , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_OuterLoop           :', Timer_Collisions_OuterLoop          , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_InnerLoop           :', Timer_Collisions_InnerLoop          , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_ComputeOpacity      :', Timer_Collisions_ComputeOpacity     , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_ComputeRates        :', Timer_Collisions_ComputeRates       , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_InitializeRHS       :', Timer_Collisions_InitializeRHS      , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_NeutrinoRHS         :', Timer_Collisions_NeutrinoRHS        , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_MatterRHS           :', Timer_Collisions_MatterRHS          , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_SolveLS             :', Timer_Collisions_SolveLS            , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_UpdateFP            :', Timer_Collisions_UpdateFP           , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_CheckOuter          :', Timer_Collisions_CheckOuter         , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_Collisions_CheckInner          :', Timer_Collisions_CheckInner         , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_TCI                              :', Timer_TCI                           , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_TCI_Permute                    :', Timer_TCI_Permute                   , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_TCI_LinearAlgebra              :', Timer_TCI_LinearAlgebra             , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_TCI_Compute                    :', Timer_TCI_Compute                   , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_SL                               :', Timer_SL                            , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_SL_Permute                     :', Timer_SL_Permute                    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_SL_LinearAlgebra               :', Timer_SL_LinearAlgebra              , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_SL_MinMod                      :', Timer_SL_MinMod                     , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_SL_ReplaceSlopes               :', Timer_SL_ReplaceSlopes              , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_SL_Correction                  :', Timer_SL_Correction                 , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '  Timer_PL                               :', Timer_PL                            , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_PL_Permute                     :', Timer_PL_Permute                    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_PL_CellAverage                 :', Timer_PL_CellAverage                , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_PL_PointValues                 :', Timer_PL_PointValues                , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_PL_Theta_1                     :', Timer_PL_Theta_1                    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_PL_Theta_2                     :', Timer_PL_Theta_2                    , ' s'
    WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
      '    Timer_PL_EnergyLimiter               :', Timer_PL_EnergyLimiter              , ' s'
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
