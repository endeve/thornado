MODULE TimersModule_AMReX_Euler

  USE amrex_fort_module, ONLY: &
    AR => amrex_real
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_min, &
    amrex_parallel_reduce_max, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_nprocs
  USE MyAmrModule, ONLY: &
    nLevels

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  LOGICAL,  PUBLIC :: TimeIt_AMReX_Euler = .FALSE.

  REAL(AR), PUBLIC :: Timer_AMReX_Euler_Program;          INTEGER :: iT_P   = 1

  ! --- fmain ---
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_Initialize;       INTEGER :: iT_I   = 2
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_MPI_Barrier;      INTEGER :: iT_B   = 3
  REAL(AR), PUBLIC :: Timer_AMReX_ComputeTimeStep_Euler;  INTEGER :: iT_CTS = 4
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_UpdateFluid;      INTEGER :: iT_UF  = 5
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_InputOutput;      INTEGER :: iT_IO  = 6
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_Finalize;         INTEGER :: iT_F   = 7

  ! --- AMReX-specific ---
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_DataTransfer;     INTEGER :: iT_DT  = 8
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_InteriorBC;       INTEGER :: iT_IBC = 9
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_CopyMultiFab;     INTEGER :: iT_CMF = 10
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_ConstructEdgeMap; INTEGER :: iT_CEM = 11
  REAL(AR), PUBLIC :: Timer_AMReX_Euler_GetBC;            INTEGER :: iT_GBC = 12

  INTEGER :: nT = 12

  PUBLIC :: InitializeTimers_AMReX_Euler
  PUBLIC :: FinalizeTimers_AMReX_Euler
  PUBLIC :: TimersStart_AMReX_Euler
  PUBLIC :: TimersStop_AMReX_Euler
  PUBLIC :: TimersWtime_AMReX

  INTEGER :: nProcs

CONTAINS


  SUBROUTINE InitializeTimers_AMReX_Euler

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    nProcs = amrex_parallel_nprocs()

    Timer_AMReX_Euler_Program          = 0.0_AR

    Timer_AMReX_Euler_Initialize       = 0.0_AR
    Timer_AMReX_Euler_MPI_Barrier      = 0.0_AR
    Timer_AMReX_ComputeTimeStep_Euler  = 0.0_AR
    Timer_AMReX_Euler_UpdateFluid      = 0.0_AR
    Timer_AMReX_Euler_InputOutput      = 0.0_AR
    Timer_AMReX_Euler_Finalize         = 0.0_AR

    Timer_AMReX_Euler_DataTransfer     = 0.0_AR
    Timer_AMReX_Euler_InteriorBC       = 0.0_AR
    Timer_AMReX_Euler_CopyMultiFab     = 0.0_AR
    Timer_AMReX_Euler_ConstructEdgeMap = 0.0_AR
    Timer_AMReX_Euler_GetBC            = 0.0_AR

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Program )

    RETURN
  END SUBROUTINE InitializeTimers_AMReX_Euler


  SUBROUTINE FinalizeTimers_AMReX_Euler

    REAL(AR) :: Timer(nT), TotalTime

    REAL(AR) :: TimerSum(nT), TimerMin(nT), TimerMax(nT), TimerAve(nT)

    CHARACTER(32) :: OutFMT = '(8x,A,ES13.6E3,A5,F7.3,A2)'
    CHARACTER(64) :: OutMMA = '(10x,A6,ES13.6E3,A4,A6,ES13.6E3,A4,A6,ES13.6E3,A2)'
    INTEGER       :: iT

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Program )

    Timer(iT_P  ) = Timer_AMReX_Euler_Program

    Timer(iT_I  ) = Timer_AMReX_Euler_Initialize
    Timer(iT_B  ) = Timer_AMReX_Euler_MPI_Barrier
    Timer(iT_CTS) = Timer_AMReX_ComputeTimeStep_Euler
    Timer(iT_UF ) = Timer_AMReX_Euler_UpdateFluid
    Timer(iT_IO ) = Timer_AMReX_Euler_InputOutput
    Timer(iT_F  ) = Timer_AMReX_Euler_Finalize

    Timer(iT_DT ) = Timer_AMReX_Euler_DataTransfer
    Timer(iT_IBC) = Timer_AMReX_Euler_InteriorBC
    Timer(iT_CMF) = Timer_AMReX_Euler_CopyMultiFab
    Timer(iT_CEM) = Timer_AMReX_Euler_ConstructEdgeMap
    Timer(iT_GBC) = Timer_AMReX_Euler_GetBC

    DO iT = 1, nT

      CALL SumMinMaxAve( Timer(iT), TimerSum(iT), &
                         TimerMin(iT), TimerMax(iT), TimerAve(iT) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(6x,A28)') 'Timers (AMReX-Euler) Summary'
      WRITE(*,'(6x,A28)') '----------------------------'
      WRITE(*,*)

      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (sum) = ', TimerSum(iT_P), ' s'
      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (min) = ', TimerMin(iT_P), ' s'
      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (max) = ', TimerMax(iT_P), ' s'
      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (ave) = ', TimerAve(iT_P), ' s'
      WRITE(*,*)

      TotalTime = 0.0_AR
      DO iT = iT_I, iT_F
        TotalTime = TotalTime + TimerSum(iT)
      END DO

      WRITE(*,'(8x,A26,F7.3,A2)') &
        'Timers / ProgramRunTime = ', &
        100.0_AR * TotalTime / TimerSum(iT_P), ' %'

      WRITE(*,*)
      WRITE(*,'(6x,A5)') 'fmain'
      WRITE(*,'(6x,A5)') '-----'
      WRITE(*,*)

      iT = iT_I
      WRITE(*,TRIM(OutFMT)) &
        'Initialize:        ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_CTS
      WRITE(*,TRIM(OutFMT)) &
        'Compute Time-Step: ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_UF
      WRITE(*,TRIM(OutFMT)) &
        'Update Fluid:      ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_IO
      WRITE(*,TRIM(OutFMT)) &
        'Input/Output:      ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_F
      WRITE(*,TRIM(OutFMT)) &
        'Finalize:          ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      WRITE(*,*)
      WRITE(*,'(6x,A)') 'AMReX-specific'
      WRITE(*,'(6x,A)') '--------------'
      WRITE(*,*)

      iT = iT_DT
      WRITE(*,TRIM(OutFMT)) &
        'Data Transfer:                ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_IBC
      WRITE(*,TRIM(OutFMT)) &
        'Interior Boundary Conditions: ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_CMF
      WRITE(*,TRIM(OutFMT)) &
        'Copy MultiFab:                ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_CEM
      WRITE(*,TRIM(OutFMT)) &
        'Construct Edge-Map:           ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

      iT = iT_GBC
      WRITE(*,TRIM(OutFMT)) &
        'Get Boundary-Conditions:      ', TimerSum(iT), ' s = ', &
        100.0_AR * TimerSum(iT) / TimerSum(iT_P), ' %'
      WRITE(*,TRIM(OutMMA)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Ave = ', TimerAve(iT), ' s'

      WRITE(*,*)

    END IF

    RETURN
  END SUBROUTINE FinalizeTimers_AMReX_Euler


  SUBROUTINE TimersStart_AMReX_Euler( Timer )

    REAL(AR), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    Timer = Timer - TimersWtime_AMReX()

    RETURN
  END SUBROUTINE TimersStart_AMReX_Euler


  SUBROUTINE TimersStop_AMReX_Euler( Timer )

    REAL(AR), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    Timer = Timer + TimersWtime_AMReX()

    RETURN
  END SUBROUTINE TimersStop_AMReX_Euler


  REAL(AR) FUNCTION TimersWtime_AMReX()

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    TimersWtime_AMReX = MPI_WTIME()

    RETURN
  END FUNCTION TimersWtime_AMReX


  SUBROUTINE SumMinMaxAve &
    ( Timer, TimerSum, TimerMin, TimerMax, TimerAve )

    REAL(AR), INTENT(in)    :: Timer
    REAL(AR), INTENT(inout) :: TimerSum, TimerMin, TimerMax, TimerAve

    TimerSum = Timer
    CALL amrex_parallel_reduce_sum( TimerSum )

    TimerMin = Timer
    CALL amrex_parallel_reduce_min( TimerMin )

    TimerMax = Timer
    CALL amrex_parallel_reduce_max( TimerMax )

    TimerAve = Timer
    CALL amrex_parallel_reduce_sum( TimerAve )
    TimerAve = TimerAve / nProcs

  END SUBROUTINE SumMinMaxAve


END MODULE TimersModule_AMReX_Euler
