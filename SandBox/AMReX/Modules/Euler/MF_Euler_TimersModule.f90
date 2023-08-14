MODULE MF_Euler_TimersModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor,  &
    amrex_parallel_reduce_min, &
    amrex_parallel_reduce_max, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_nprocs

  ! --- Thornado Modules ---

  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    InitializeTimers_Euler, &
    FinalizeTimers_Euler, &
    Timer_Euler_DG, &
    Timer_Euler_Increment, &
    Timer_Euler_Divergence, &
    Timer_Euler_SurfaceTerm, &
    Timer_Euler_VolumeTerm, &
    Timer_Euler_Geometry, &
    Timer_Euler_DD_TCI, &
    Timer_Euler_DD_SD, &
    Timer_Euler_SlopeLimiter, &
    Timer_Euler_PositivityLimiter, &
    Timer_Euler_BoundaryConditions, &
    Timer_Euler_SL_CharDecomp

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
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

  INTEGER :: iT_DG  = 1
  INTEGER :: iT_INC = 2
  INTEGER :: iT_DIV = 3
  INTEGER :: iT_SUR = 4
  INTEGER :: iT_VOL = 5
  INTEGER :: iT_GEO = 6
  INTEGER :: iT_TCI = 7
  INTEGER :: iT_SD  = 8
  INTEGER :: iT_SL  = 9
  INTEGER :: iT_CD  = 10
  INTEGER :: iT_PL  = 11
  INTEGER :: iT_BC  = 12

  INTEGER :: nTimers = 12

  LOGICAL :: WriteMMS

  INTEGER, PARAMETER :: nChar = 48

CONTAINS


  SUBROUTINE InitializeTimers_AMReX_Euler( WriteMMS_Option )

    LOGICAL, INTENT(in), OPTIONAL :: WriteMMS_Option

    WriteMMS = .FALSE.
    IF( PRESENT( WriteMMS_Option ) ) &
      WriteMMS = WriteMMS_Option

    TimeIt_Euler = .FALSE.

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    TimeIt_Euler = .TRUE.

    CALL InitializeTimers_Euler

  END SUBROUTINE InitializeTimers_AMReX_Euler


  SUBROUTINE FinalizeTimers_AMReX_Euler( Timer_ProgramAve )

    REAL(DP), INTENT(in) :: Timer_ProgramAve

    REAL(DP) :: Timer(nTimers)
    REAL(DP) :: TimerSum(nTimers), TimerMin(nTimers), &
                TimerMax(nTimers), TimerAve(nTimers)
    CHARACTER(nChar) :: Label(nTimers)

    INTEGER :: iT

    CHARACTER(128) :: &
      OutFMT = '(8x,A32,A,ES11.3E3,A,ES11.3E3)'
    CHARACTER(128) :: &
      OutMMS = '(10x,A,ES13.6E3,A,A,ES13.6E3,A,A,ES13.6E3,A)'

    IF( .NOT. TimeIt_AMReX_Euler ) RETURN

    CALL FinalizeTimers_Euler &
           ( Verbose_Option = .FALSE., &
             SuppressApplicationDriver_Option = .TRUE. )

    Timer(iT_DG ) = Timer_Euler_DG
    Timer(iT_INC) = Timer_Euler_Increment
    Timer(iT_DIV) = Timer_Euler_Divergence
    Timer(iT_SUR) = Timer_Euler_SurfaceTerm
    Timer(iT_VOL) = Timer_Euler_VolumeTerm
    Timer(iT_GEO) = Timer_Euler_Geometry
    Timer(iT_TCI) = Timer_Euler_DD_TCI
    Timer(iT_SD ) = Timer_Euler_DD_SD
    Timer(iT_SL ) = Timer_Euler_SlopeLimiter
    Timer(iT_CD ) = Timer_Euler_SL_CharDecomp
    Timer(iT_PL ) = Timer_Euler_PositivityLimiter
    Timer(iT_BC ) = Timer_Euler_BoundaryConditions

    DO iT = 1, nTimers

      CALL SumMinMaxAve( Timer(iT), TimerSum(iT), &
                         TimerMin(iT), TimerMax(iT), TimerAve(iT) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(4x,A)') 'Timers Summary (Euler)'
      WRITE(*,'(4x,A)') '----------------------'
      WRITE(*,*)

      Label(iT_DG ) = 'DG Discretization'
      Label(iT_INC) = '  Increment'
      Label(iT_DIV) = '  Divergence'
      Label(iT_SUR) = '  Surface Term'
      Label(iT_VOL) = '  Volume Term'
      Label(iT_GEO) = '  Geometry'
      Label(iT_TCI) = 'Troubled-Cell Indicator'
      Label(iT_SD ) = 'Shock Detector'
      Label(iT_SL ) = 'Slope Limiter'
      Label(iT_CD ) = '  Characteristic Decomposition'
      Label(iT_PL ) = 'Positivity Limiter'
      Label(iT_BC ) = 'Boundary Conditions'

      WRITE(*,'(6x,A)') 'Native thornado'
      WRITE(*,'(6x,A)') '---------------'
      WRITE(*,*)

      DO iT = iT_DG, nTimers

        CALL WriteTimer &
               ( OutFMT, OutMMS, Timer_ProgramAve, &
                 Label(iT), &
                 TimerAve(iT), TimerMin(iT), TimerMax(iT), TimerSum(iT) )

      END DO
      WRITE(*,*)

    END IF

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


END MODULE MF_Euler_TimersModule
