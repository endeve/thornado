MODULE MF_Euler_TimersModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor,  &
    amrex_parallel_reduce_min,   &
    amrex_parallel_reduce_max,   &
    amrex_parallel_reduce_sum,   &
    amrex_parallel_nprocs

  ! --- Thornado Modules ---

  USE TimersModule_Euler,    ONLY: &
    InitializeTimers_Euler,             &
    FinalizeTimers_Euler,               &
    Timer_Euler_DG,                     &
    Timer_Euler_Increment,              &
    Timer_Euler_Divergence,             &
    Timer_Euler_SurfaceTerm,            &
    Timer_Euler_VolumeTerm,             &
    Timer_Euler_Geometry,               &
    Timer_Euler_DD_TCI,                 &
    Timer_Euler_DD_SD,                  &
    Timer_Euler_SlopeLimiter,           &
    Timer_Euler_PositivityLimiter,      &
    Timer_Euler_BoundaryConditions,     &
    Timer_Euler_SL_CharDecomp

  ! --- Local Modules ---

  USE MF_KindModule,         ONLY: &
    DP, &
    Zero

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeTimers_AMReX_Euler
  PUBLIC :: FinalizeTimers_AMReX_Euler
  PUBLIC :: TimersStart_AMReX_Euler
  PUBLIC :: TimersStop_AMReX_Euler
  PUBLIC :: TimersWtime_AMReX

  LOGICAL,  PUBLIC :: TimeIt_AMReX_Euler = .FALSE.

  REAL(DP), PUBLIC :: Timer_AMReX_Euler_Program;          INTEGER :: iT_P   = 1

  ! --- ApplicationDriver ---

  REAL(DP), PUBLIC :: Timer_AMReX_Euler_Initialize;       INTEGER :: iT_I   = 2
  REAL(DP), PUBLIC :: Timer_AMReX_ComputeTimeStep_Euler;  INTEGER :: iT_CTS = 3
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_UpdateFluid;      INTEGER :: iT_UF  = 4
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_InputOutput;      INTEGER :: iT_IO  = 5
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_Finalize;         INTEGER :: iT_F   = 6

  ! --- AMReX-specific ---

  REAL(DP), PUBLIC :: Timer_AMReX_Euler_DataTransfer;     INTEGER :: iT_DT  = 7
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_Allocate;         INTEGER :: iT_AL  = 8
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_InteriorBC;       INTEGER :: iT_IBC = 9
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_CopyMultiFab;     INTEGER :: iT_CMF = 10
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_ConstructEdgeMap; INTEGER :: iT_CEM = 11
  REAL(DP), PUBLIC :: Timer_AMReX_Euler_GetBC;            INTEGER :: iT_GBC = 12

  ! --- Thornado Modules ---

  INTEGER :: iT_DG  = 13
  INTEGER :: iT_INC = 14
  INTEGER :: iT_DIV = 15
  INTEGER :: iT_SUR = 16
  INTEGER :: iT_VOL = 17
  INTEGER :: iT_GEO = 18
  INTEGER :: iT_TCI = 19
  INTEGER :: iT_SD  = 20
  INTEGER :: iT_SL  = 21
  INTEGER :: iT_PL  = 22
  INTEGER :: iT_BC  = 23
  INTEGER :: iT_CD  = 24

  INTEGER :: nTimers = 24

  REAL(DP), PARAMETER :: Hundred = 100.0_DP

CONTAINS


  SUBROUTINE InitializeTimers_AMReX_Euler

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    Timer_AMReX_Euler_Program          = Zero

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Program )

    Timer_AMReX_Euler_Initialize       = Zero
    Timer_AMReX_ComputeTimeStep_Euler  = Zero
    Timer_AMReX_Euler_UpdateFluid      = Zero
    Timer_AMReX_Euler_InputOutput      = Zero
    Timer_AMReX_Euler_Finalize         = Zero

    Timer_AMReX_Euler_DataTransfer     = Zero
    Timer_AMReX_Euler_Allocate         = Zero
    Timer_AMReX_Euler_InteriorBC       = Zero
    Timer_AMReX_Euler_CopyMultiFab     = Zero
    Timer_AMReX_Euler_ConstructEdgeMap = Zero
    Timer_AMReX_Euler_GetBC            = Zero

    CALL InitializeTimers_Euler

  END SUBROUTINE InitializeTimers_AMReX_Euler


  SUBROUTINE FinalizeTimers_AMReX_Euler( WriteAtIntermediateTime_Option )

    LOGICAL, INTENT(in), OPTIONAL :: WriteAtIntermediateTime_Option

    REAL(DP) :: Timer(nTimers), TotalTime
    REAL(DP) :: TimerSum(nTimers), TimerMin(nTimers), &
                TimerMax(nTimers), TimerAve(nTimers)

    CHARACTER(32) :: &
      OutFMT = '(8x,A,ES13.6E3,A,F7.3,A)'
    CHARACTER(64) :: &
      OutMMS = '(10x,A,ES13.6E3,A,A,ES13.6E3,A,A,ES13.6E3,A)'
    INTEGER       :: iT

    LOGICAL :: WriteAtIntermediateTime

    WriteAtIntermediateTime = .FALSE.
    IF( PRESENT( WriteAtIntermediateTime_Option ) ) &
      WriteAtIntermediateTime = WriteAtIntermediateTime_Option

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    CALL FinalizeTimers_Euler &
           ( Verbose_Option = .FALSE., &
             SuppressApplicationDriver_Option = .TRUE. )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Program )

    Timer(iT_P  ) = Timer_AMReX_Euler_Program

    Timer(iT_I  ) = Timer_AMReX_Euler_Initialize
    Timer(iT_CTS) = Timer_AMReX_ComputeTimeStep_Euler
    Timer(iT_UF ) = Timer_AMReX_Euler_UpdateFluid
    Timer(iT_IO ) = Timer_AMReX_Euler_InputOutput
    Timer(iT_F  ) = Timer_AMReX_Euler_Finalize

    Timer(iT_DT ) = Timer_AMReX_Euler_DataTransfer
    Timer(iT_AL ) = Timer_AMReX_Euler_Allocate
    Timer(iT_IBC) = Timer_AMReX_Euler_InteriorBC
    Timer(iT_CMF) = Timer_AMReX_Euler_CopyMultiFab
    Timer(iT_CEM) = Timer_AMReX_Euler_ConstructEdgeMap
    Timer(iT_GBC) = Timer_AMReX_Euler_GetBC

    Timer(iT_DG ) = Timer_Euler_DG
    Timer(iT_INC) = Timer_Euler_Increment
    Timer(iT_DIV) = Timer_Euler_Divergence
    Timer(iT_SUR) = Timer_Euler_SurfaceTerm
    Timer(iT_VOL) = Timer_Euler_VolumeTerm
    Timer(iT_GEO) = Timer_Euler_Geometry
    Timer(iT_TCI) = Timer_Euler_DD_TCI
    Timer(iT_SD ) = Timer_Euler_DD_SD
    Timer(iT_SL ) = Timer_Euler_SlopeLimiter
    Timer(iT_PL ) = Timer_Euler_PositivityLimiter
    Timer(iT_BC ) = Timer_Euler_BoundaryConditions
    Timer(iT_CD ) = Timer_Euler_SL_CharDecomp

    DO iT = 1, nTimers

      CALL SumMinMaxAve( Timer(iT), TimerSum(iT), &
                         TimerMin(iT), TimerMax(iT), TimerAve(iT) )

    END DO

    IF( amrex_parallel_ioprocessor() .AND. WriteAtIntermediateTime )THEN

      WRITE(*,*)
      WRITE(*,'(6x,A28)') 'Timers (AMReX-Euler) Summary'
      WRITE(*,'(6x,A28)') '----------------------------'
      WRITE(*,*)

      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (ave) = ', TimerAve(iT_P), ' s'
      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (min) = ', TimerMin(iT_P), ' s'
      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (max) = ', TimerMax(iT_P), ' s'
      WRITE(*,'(8x,A23,ES13.6E3,A2)') &
        'Total run-time (sum) = ', TimerSum(iT_P), ' s'
      WRITE(*,*)

      TotalTime = Zero
      DO iT = iT_I, iT_F
        TotalTime = TotalTime + TimerAve(iT)
      END DO

      WRITE(*,'(8x,A26,F7.3,A2)') &
        'Timers / ProgramRunTime = ', &
        Hundred * TotalTime / TimerAve(iT_P), ' %'

      WRITE(*,*)
      WRITE(*,'(6x,A5)') 'fmain'
      WRITE(*,'(6x,A5)') '-----'
      WRITE(*,*)

      iT = iT_I
      WRITE(*,TRIM(OutFMT)) &
        'Initialize:        ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_CTS
      WRITE(*,TRIM(OutFMT)) &
        'Compute Time-Step: ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_UF
      WRITE(*,TRIM(OutFMT)) &
        'Update Fluid:      ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_IO
      WRITE(*,TRIM(OutFMT)) &
        'Input/Output:      ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_F
      WRITE(*,TRIM(OutFMT)) &
        'Finalize:          ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Native thornado'
      WRITE(*,'(6x,A)') '---------------'
      WRITE(*,*)

      iT = iT_DG
      WRITE(*,TRIM(OutFMT)) &
        'DG Discretization:      ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_INC
      WRITE(*,TRIM(OutFMT)) &
        'Increment:              ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_DIV
      WRITE(*,TRIM(OutFMT)) &
        'Divergence:             ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_SUR
      WRITE(*,TRIM(OutFMT)) &
        '  SurfaceTerm:          ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
          'Min = ', TimerMin(iT), ' s, ', &
          'Max = ', TimerMax(iT), ' s, ', &
          'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_VOL
      WRITE(*,TRIM(OutFMT)) &
        '  VolumeTerm:           ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
          'Min = ', TimerMin(iT), ' s, ', &
          'Max = ', TimerMax(iT), ' s, ', &
          'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_GEO
      WRITE(*,TRIM(OutFMT)) &
        'Geometry:               ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_TCI
      WRITE(*,TRIM(OutFMT)) &
        'Troubled-Cell Indicator:', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_SD
      WRITE(*,TRIM(OutFMT)) &
        'Shock Detector:         ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_SL
      WRITE(*,TRIM(OutFMT)) &
        'SlopeLimiter:           ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_PL
      WRITE(*,TRIM(OutFMT)) &
        'PositivityLimiter:      ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_BC
      WRITE(*,TRIM(OutFMT)) &
        'BoundaryConditions:     ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_CD
      WRITE(*,TRIM(OutFMT)) &
        'Characteristic Decomp:  ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      WRITE(*,*)
      WRITE(*,'(6x,A)') 'AMReX-specific'
      WRITE(*,'(6x,A)') '--------------'
      WRITE(*,*)

      iT = iT_DT
      WRITE(*,TRIM(OutFMT)) &
        'Data Transfer:                ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_AL
      WRITE(*,TRIM(OutFMT)) &
        'Allocate/Deallocate:          ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_IBC
      WRITE(*,TRIM(OutFMT)) &
        'Interior Boundary Conditions: ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_CMF
      WRITE(*,TRIM(OutFMT)) &
        'Copy MultiFab:                ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_CEM
      WRITE(*,TRIM(OutFMT)) &
        'Construct Edge-Map:           ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

      iT = iT_GBC
      WRITE(*,TRIM(OutFMT)) &
        'Get Boundary-Conditions:      ', TimerAve(iT), ' s = ', &
        Hundred * TimerAve(iT) / TimerAve(iT_P), ' %'
      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin(iT), ' s, ', &
        'Max = ', TimerMax(iT), ' s, ', &
        'Sum = ', TimerSum(iT), ' s'
      WRITE(*,*)

    END IF

    IF( WriteAtIntermediateTime ) &
      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Program )

  END SUBROUTINE FinalizeTimers_AMReX_Euler


  SUBROUTINE TimersStart_AMReX_Euler( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    Timer = Timer - TimersWtime_AMReX()

  END SUBROUTINE TimersStart_AMReX_Euler


  SUBROUTINE TimersStop_AMReX_Euler( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    Timer = Timer + TimersWtime_AMReX()

  END SUBROUTINE TimersStop_AMReX_Euler


  REAL(DP) FUNCTION TimersWtime_AMReX() RESULT( Time )

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

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


END MODULE MF_Euler_TimersModule
