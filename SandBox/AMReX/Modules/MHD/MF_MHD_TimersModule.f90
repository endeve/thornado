MODULE MF_MHD_TimersModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor,  &
    amrex_parallel_reduce_min, &
    amrex_parallel_reduce_max, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_nprocs

  ! --- Thornado Modules ---

  USE TimersModule_MHD, ONLY: &
    TimeIt_MHD, &
    InitializeTimers_MHD, &
    FinalizeTimers_MHD, &
    Timer_MHD_DG, &
    Timer_MHD_Increment, &
    Timer_MHD_Divergence, &
    Timer_MHD_SurfaceTerm, &
    Timer_MHD_VolumeTerm, &
    Timer_MHD_Geometry, &
    Timer_MHD_DD_TCI, &
    Timer_MHD_DD_SD, &
    Timer_MHD_SlopeLimiter, &
    Timer_MHD_SL_CharDecomp, &
    Timer_MHD_PositivityLimiter, &
    Timer_MHD_BoundaryConditions

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeTimers_AMReX_MHD
  PUBLIC :: FinalizeTimers_AMReX_MHD
  PUBLIC :: TimersStart_AMReX_MHD
  PUBLIC :: TimersStop_AMReX_MHD
  PUBLIC :: TimersWtime_AMReX

  LOGICAL, PUBLIC :: TimeIt_AMReX_MHD = .FALSE.

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


  SUBROUTINE InitializeTimers_AMReX_MHD( WriteMMS_Option )

    LOGICAL, INTENT(in), OPTIONAL :: WriteMMS_Option

    WriteMMS = .FALSE.
    IF( PRESENT( WriteMMS_Option ) ) &
      WriteMMS = WriteMMS_Option

    TimeIt_MHD = .FALSE.

    IF( .NOT. TimeIt_AMReX_MHD ) RETURN

    TimeIt_MHD = .TRUE.

    CALL InitializeTimers_MHD

  END SUBROUTINE InitializeTimers_AMReX_MHD


  SUBROUTINE FinalizeTimers_AMReX_MHD( Timer_ProgramAve, Verbose_Option )

    REAL(DP), INTENT(in) :: Timer_ProgramAve
    LOGICAL , INTENT(in), OPTIONAL :: Verbose_Option

    REAL(DP) :: Timer(nTimers)
    REAL(DP) :: TimerSum(nTimers), TimerMin(nTimers), &
                TimerMax(nTimers), TimerAve(nTimers)
    CHARACTER(nChar) :: Label(nTimers)
    LOGICAL          :: Verbose

    INTEGER :: iT

    CHARACTER(128) :: &
      OutFMT = '(8x,A32,A,ES11.3E3,A,ES11.3E3)'
    CHARACTER(128) :: &
      OutMMS = '(10x,A,ES13.6E3,A,A,ES13.6E3,A,A,ES13.6E3,A)'

    IF( .NOT. TimeIt_AMReX_MHD ) RETURN

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL FinalizeTimers_MHD &
           ( Verbose_Option = .FALSE., &
             SuppressApplicationDriver_Option = .TRUE. )

    Timer(iT_DG ) = Timer_MHD_DG
    Timer(iT_INC) = Timer_MHD_Increment
    Timer(iT_DIV) = Timer_MHD_Divergence
    Timer(iT_SUR) = Timer_MHD_SurfaceTerm
    Timer(iT_VOL) = Timer_MHD_VolumeTerm
    Timer(iT_GEO) = Timer_MHD_Geometry
    Timer(iT_TCI) = Timer_MHD_DD_TCI
    Timer(iT_SD ) = Timer_MHD_DD_SD
    Timer(iT_SL ) = Timer_MHD_SlopeLimiter
    Timer(iT_CD ) = Timer_MHD_SL_CharDecomp
    Timer(iT_PL ) = Timer_MHD_PositivityLimiter
    Timer(iT_BC ) = Timer_MHD_BoundaryConditions

    DO iT = 1, nTimers

      CALL SumMinMaxAve( Timer(iT), TimerSum(iT), &
                         TimerMin(iT), TimerMax(iT), TimerAve(iT) )

    END DO

    IF( Verbose )THEN

      WRITE(*,'(4x,A)') 'Timers Summary (MHD)'
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

  END SUBROUTINE FinalizeTimers_AMReX_MHD


  SUBROUTINE TimersStart_AMReX_MHD( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX_MHD ) RETURN

    Timer = Timer - TimersWtime_AMReX()

  END SUBROUTINE TimersStart_AMReX_MHD


  SUBROUTINE TimersStop_AMReX_MHD( Timer )

    REAL(DP), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX_MHD ) RETURN

    Timer = Timer + TimersWtime_AMReX()

  END SUBROUTINE TimersStop_AMReX_MHD


  REAL(DP) FUNCTION TimersWtime_AMReX() RESULT( Time )

    IF( .NOT. TimeIt_AMReX_MHD ) RETURN

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
    REAL(DP)        , INTENT(in) :: TimerAve_P
    REAL(DP)        , INTENT(in) :: TimerAve, TimerMin, TimerMax, TimerSum

    WRITE(*,TRIM(OutFMT)) &
      Label, ' : ', TimerAve, ' s = ', TimerAve / TimerAve_P

    IF( WriteMMS )THEN

      WRITE(*,TRIM(OutMMS)) &
        'Min = ', TimerMin, ' s, ', &
        'Max = ', TimerMax, ' s, ', &
        'Sum = ', TimerSum, ' s'

    END IF

  END SUBROUTINE WriteTimer


END MODULE MF_MHD_TimersModule
