MODULE MF_TimersModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor,  &
    amrex_parallel_reduce_min, &
    amrex_parallel_reduce_max, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_nprocs

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_Euler_TimersModule, ONLY: &
    TimeIt_AMReX_Euler, &
    InitializeTimers_AMReX_Euler, &
    FinalizeTimers_AMReX_Euler

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeTimers_AMReX
  PUBLIC :: FinalizeTimers_AMReX
  PUBLIC :: TimersStart_AMReX
  PUBLIC :: TimersStop_AMReX
  PUBLIC :: TimersWtime_AMReX

  LOGICAL,  PUBLIC :: TimeIt_AMReX

  REAL(DP), PUBLIC :: Timer_AMReX_Program;               INTEGER :: iT_P   = 1

  REAL(DP), PUBLIC :: Timer_AMReX_Initialize;            INTEGER :: iT_I   = 2
  REAL(DP), PUBLIC :: Timer_AMReX_ComputeTimeStep_Euler; INTEGER :: iT_CTE = 3
  REAL(DP), PUBLIC :: Timer_AMReX_UpdateFluid;           INTEGER :: iT_UF  = 4
  REAL(DP), PUBLIC :: Timer_AMReX_InputOutput;           INTEGER :: iT_IO  = 5
  REAL(DP), PUBLIC :: Timer_AMReX_Finalize;              INTEGER :: iT_F   = 6

  REAL(DP), PUBLIC :: Timer_AMReX_CopyMultiFab;          INTEGER :: iT_CMF = 7
  REAL(DP), PUBLIC :: Timer_AMReX_Allocate_X;            INTEGER :: iT_ALX = 8
  REAL(DP), PUBLIC :: Timer_AMReX_Allocate_Z;            INTEGER :: iT_ALZ = 9
  REAL(DP), PUBLIC :: Timer_AMReX_PermuteData_X;         INTEGER :: iT_PDX = 10
  REAL(DP), PUBLIC :: Timer_AMReX_PermuteData_Z;         INTEGER :: iT_PDZ = 11
  REAL(DP), PUBLIC :: Timer_AMReX_FillPatch;             INTEGER :: iT_FP  = 12
  REAL(DP), PUBLIC :: Timer_AMReX_AverageDown;           INTEGER :: iT_AD  = 13

  REAL(DP), PUBLIC :: Timer_AMReX_GravitySolve;          INTEGER :: iT_GS  = 14

  INTEGER, PARAMETER :: nTimers = 14

  LOGICAL :: WriteMMS

  INTEGER, PARAMETER :: nChar = 48

CONTAINS


  SUBROUTINE InitializeTimers_AMReX( WriteMMS_Option )

    LOGICAL, INTENT(in), OPTIONAL :: WriteMMS_Option

    WriteMMS = .FALSE.
    IF( PRESENT( WriteMMS_Option ) ) &
      WriteMMS = WriteMMS_Option

    IF( .NOT. TimeIt_AMReX ) RETURN

    Timer_AMReX_Program = Zero

    CALL TimersStart_AMReX( Timer_AMReX_Program )

    Timer_AMReX_Initialize            = Zero
    Timer_AMReX_ComputeTimeStep_Euler = Zero
    Timer_AMReX_UpdateFluid           = Zero
    Timer_AMReX_InputOutput           = Zero
    Timer_AMReX_Finalize              = Zero

    Timer_AMReX_CopyMultiFab  = Zero
    Timer_AMReX_Allocate_X    = Zero
    Timer_AMReX_Allocate_Z    = Zero
    Timer_AMReX_PermuteData_X = Zero
    Timer_AMReX_PermuteData_Z = Zero
    Timer_AMReX_FillPatch     = Zero
    Timer_AMReX_AverageDown   = Zero

    Timer_AMReX_GravitySolve  = Zero

    IF( TimeIt_AMReX_Euler ) &
      CALL InitializeTimers_AMReX_Euler &
             ( WriteMMS_Option = WriteMMS )

  END SUBROUTINE InitializeTimers_AMReX


  SUBROUTINE FinalizeTimers_AMReX( RestartProgramTimer_Option )

    LOGICAL, INTENT(in), OPTIONAL :: RestartProgramTimer_Option

    REAL(DP) :: Timer(nTimers), TotalTime
    REAL(DP) :: TimerSum(nTimers), TimerMin(nTimers), &
                TimerMax(nTimers), TimerAve(nTimers)
    CHARACTER(nChar) :: Label(nTimers)

    INTEGER :: iT

    CHARACTER(128) :: &
      OutFMT = '(8x,A32,A,ES11.3E3,A,ES11.3E3)'
    CHARACTER(128) :: &
      OutMMS = '(10x,A,ES13.6E3,A,A,ES13.6E3,A,A,ES13.6E3,A)'

    LOGICAL :: RestartProgramTimer

    RestartProgramTimer = .FALSE.
    IF( PRESENT( RestartProgramTimer_Option ) ) &
      RestartProgramTimer = RestartProgramTimer_Option

    IF( .NOT. TimeIt_AMReX ) RETURN

    IF( RestartProgramTimer ) &
      CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

    CALL TimersStop_AMReX( Timer_AMReX_Program )

    Timer(iT_P  ) = Timer_AMReX_Program

    Timer(iT_I  ) = Timer_AMReX_Initialize
    Timer(iT_CTE) = Timer_AMReX_ComputeTimeStep_Euler
    Timer(iT_UF ) = Timer_AMReX_UpdateFluid
    Timer(iT_IO ) = Timer_AMReX_InputOutput
    Timer(iT_F  ) = Timer_AMReX_Finalize

    Timer(iT_CMF) = Timer_AMReX_CopyMultiFab
    Timer(iT_ALX) = Timer_AMReX_Allocate_X
    Timer(iT_ALZ) = Timer_AMReX_Allocate_Z
    Timer(iT_PDX) = Timer_AMReX_PermuteData_X
    Timer(iT_PDZ) = Timer_AMReX_PermuteData_Z
    Timer(iT_FP ) = Timer_AMReX_FillPatch
    Timer(iT_AD ) = Timer_AMReX_AverageDown

    Timer(iT_GS ) = Timer_AMReX_GravitySolve

    DO iT = 1, nTimers

      CALL SumMinMaxAve( Timer(iT), TimerSum(iT), &
                         TimerMin(iT), TimerMax(iT), TimerAve(iT) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') 'Timers Summary (AMReX)'
      WRITE(*,'(4x,A)') '----------------------'
      WRITE(*,*)

      WRITE(*,'(6x,A23,ES13.6E3,A2)') &
        'Total run-time (ave) = ', TimerAve(iT_P), ' s'
      WRITE(*,*)

      TotalTime = Zero
      DO iT = iT_I, iT_F
        TotalTime = TotalTime + TimerAve(iT)
      END DO

      WRITE(*,'(6x,A26,ES11.3E3)') &
        'Timers / ProgramRunTime = ', &
        TotalTime / TimerAve(iT_P)

      Label(iT_I  ) = 'Initialize'
      Label(iT_CTE) = 'Compute Timestep (Euler)'
      Label(iT_UF ) = 'Update Fields'
      Label(iT_IO ) = 'Input/Output'
      Label(iT_F  ) = 'Finalize'

      Label(iT_CMF) = 'Copy MultiFab'
      Label(iT_ALX) = 'Allocate/Deallocate (X)'
      Label(iT_ALZ) = 'Allocate/Deallocate (Z)'
      Label(iT_PDX) = 'Permute Data (X)'
      Label(iT_PDZ) = 'Permute Data (Z)'
      Label(iT_FP ) = 'Fill Patch'
      Label(iT_AD ) = 'Average Down'

      Label(iT_GS ) = 'Gravity Solve'

      WRITE(*,*)
      WRITE(*,'(6x,A)') 'AMReX-specific'
      WRITE(*,'(6x,A)') '--------------'
      WRITE(*,*)

      DO iT = iT_I, iT_F

        CALL WriteTimer &
               ( OutFMT, OutMMS, TimerAve(iT_P), &
                 Label(iT), &
                 TimerAve(iT), TimerMin(iT), TimerMax(iT), TimerSum(iT) )

      END DO
      WRITE(*,*)

      DO iT = iT_F+1, nTimers

        CALL WriteTimer &
               ( OutFMT, OutMMS, TimerAve(iT_P), &
                 Label(iT), &
                 TimerAve(iT), TimerMin(iT), TimerMax(iT), TimerSum(iT) )

      END DO
      WRITE(*,*)

    END IF

    IF( TimeIt_AMReX_Euler ) &
      CALL FinalizeTimers_AMReX_Euler( TimerAve(iT_P) )

    IF( RestartProgramTimer )THEN

      CALL TimersStart_AMReX( Timer_AMReX_Program )
      CALL TimersStart_AMReX( Timer_AMReX_InputOutput )

    END IF

  END SUBROUTINE FinalizeTimers_AMReX


  SUBROUTINE TimersStart_AMReX( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX ) RETURN

    Timer = Timer - TimersWtime_AMReX()

  END SUBROUTINE TimersStart_AMReX


  SUBROUTINE TimersStop_AMReX( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX ) RETURN

    Timer = Timer + TimersWtime_AMReX()

  END SUBROUTINE TimersStop_AMReX


  REAL(DP) FUNCTION TimersWtime_AMReX() RESULT( Time )

    IF( .NOT. TimeIt_AMReX ) RETURN

    Time = MPI_WTIME()

    RETURN
  END FUNCTION TimersWtime_AMReX


  SUBROUTINE SumMinMaxAve &
    ( Timer, TimerSum, TimerMin, TimerMax, TimerAve )

    REAL(DP), INTENT(in)    :: Timer
    REAL(DP), INTENT(inout) :: TimerSum, TimerMin, TimerMax, TimerAve

    TimerSum = Timer
    CALL amrex_parallel_reduce_sum( TimerSum )

    TimerMin = Timer
    CALL amrex_parallel_reduce_min( TimerMin )

    TimerMax = Timer
    CALL amrex_parallel_reduce_max( TimerMax )

    TimerAve = Timer
    CALL amrex_parallel_reduce_sum( TimerAve )
    TimerAve = TimerAve / amrex_parallel_nprocs()

  END SUBROUTINE SumMinMaxAve


  SUBROUTINE WriteTimer &
    ( OutFMT, OutMMS, TimerAve_P, &
      Label, &
      TimerAve, TimerMin, TimerMax, TimerSum )

    CHARACTER(nChar), INTENT(in) :: OutFMT, OutMMS, Label
    REAL(DP), INTENT(in) :: TimerAve_P
    REAL(DP), INTENT(in) :: TimerAve, TimerMin, TimerMax, TimerSum

    WRITE(*,TRIM(OutFMT)) &
      Label, ' : ', TimerAve, ' s = ', TimerAve / TimerAve_P

    IF( WriteMMS )THEN

      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin, ' s, ', &
        'Max = ', TimerMax, ' s, ', &
        'Sum = ', TimerSum, ' s'

    END IF

  END SUBROUTINE WriteTimer


END MODULE MF_TimersModule
