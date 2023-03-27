PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator, &
    amrex_parallel_myproc
  USE amrex_amrcore_module, ONLY: &
    amrex_regrid, &
    amrex_get_numlevels

  ! --- thornado Modules ---

  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE MemoryProfilingModule, ONLY: &
    WriteMemoryUsage

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeTimeStep_Euler_MF, &
    ComputeFromConserved_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    ComputeTally_Euler_MF
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    UpdateFluid_SSPRK_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    t_end, &
    t_new, &
    t_old, &
    dt, &
    CFL, &
    iCycleD, &
    iCycleW, &
    iCycleChk, &
    t_wrt, &
    t_chk, &
    dt_wrt, &
    dt_chk, &
    UseAMR, &
    iReGrid, &
    DEBUG
  USE MF_Euler_TimersModule, ONLY: &
    TimeIt_AMReX_Euler
  USE MF_TimersModule, ONLY: &
    TimeIt_AMReX, &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_InputOutput, &
    FinalizeTimers_AMReX

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: iLevel, iErr
  LOGICAL  :: wrt, chk
  REAL(DP) :: Timer_Evolution

  CHARACTER(128) :: MemFileName
  INTEGER :: UnitNo

  TimeIt_AMReX       = .TRUE.
  TimeIt_AMReX_Euler = .TRUE.

  wrt = .FALSE.
  chk = .FALSE.

  CALL InitializeProgram

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

!  UnitNo = 100 + amrex_parallel_myproc()
!  WRITE( MemFileName, '(A,I2.2,A)' ) &
!    'MemoryUsage_iProc', amrex_parallel_myproc(), '.txt'
!  OPEN( UNIT = UnitNo, FILE = TRIM( MemFileName ) )
!  WRITE( UnitNo, '(A)' ) &
!    '# StepNo, Current Time [ms], Current Memory Usage [kB], Maximum Memory Usage [kB]'
!  CLOSE( UnitNo )
!  CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

  ! --- Begin evolution ---

  DO WHILE( MAXVAL( t_new ) .LT. t_end )

    StepNo = StepNo + 1

    t_old = t_new

!    UnitNo = 100 + amrex_parallel_myproc()
!    WRITE( MemFileName, '(A,I2.2,A)' ) &
!      'MemoryUsage_iProc', amrex_parallel_myproc(), '.txt'
!    OPEN( UNIT = UnitNo, FILE = TRIM( MemFileName ), POSITION = 'APPEND' )
!    CALL WriteMemoryUsage &
!           ( UnitNo, 'Before call', StepNo(0), &
!             t_new(0) / UnitsDisplay % TimeUnit )
!    ! CALL YourFavoriteSubroutine
!    CALL WriteMemoryUsage &
!           ( UnitNo, 'Before call', StepNo(0), &
!             t_new(0) / UnitsDisplay % TimeUnit )
!    CLOSE( UnitNo )

    CALL ReGrid

    CALL ComputeTimeStep_Euler_MF( MF_uGF, MF_uCF, CFL, dt )

    dt = MINVAL( dt )

    IF( MAXVAL( t_old + dt ) .LT. t_end )THEN

      t_new = t_old + dt

    ELSE

      dt = t_end - t_old

      t_new = t_end

    END IF

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,'(A)') 'CALL UpdateFluid_SSPRK_MF'

    END IF

    CALL UpdateFluid_SSPRK_MF

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(A)') 'CALL ComputeFromConserved_Euler_MF'
        WRITE(*,*)

      END IF

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      IF( ( MOD( StepNo(0), iCycleD ) .EQ. 0 ) .OR. DEBUG )THEN

        WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
          'StepNo: ', StepNo(0), ' t = ', t_new(0) / UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel ), &
          ' dt = ', dt(0) /  UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel )

     END IF

    END IF

    CALL WritePlotFile

    CALL WriteCheckpointFile

  END DO

  ! --- END of evolution ---

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(2x,A,ES13.6E3,A)') &
      'Total evolution time: ', MPI_WTIME() - Timer_Evolution, ' s'

  END IF

  StepNo = StepNo + 1

  CALL FinalizeProgram

CONTAINS


  SUBROUTINE WritePlotFile

    CALL TimersStart_AMReX( Timer_AMReX_InputOutput )

    IF( iCycleW .GT. 0 )THEN

      IF( MOD( StepNo(0), iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.

    ELSE

      IF( ALL( t_new + dt .GT. t_wrt ) )THEN

        t_wrt = t_wrt + dt_wrt
        wrt   = .TRUE.

      END IF

    END IF

    IF( wrt )THEN

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() ) &
          WRITE(*,'(A)') 'CALL WritePlotFile'

      END IF

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_PlotFile &
             ( t_new(0), StepNo, MF_uGF, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF, &
               MF_uDF_Option = MF_uDF )

      CALL ComputeTally_Euler_MF &
             ( t_new, MF_uGF, MF_uCF, Verbose_Option = .TRUE. )

      wrt = .FALSE.

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

  END SUBROUTINE WritePlotFile


  SUBROUTINE WriteCheckpointFile

    CALL TimersStart_AMReX( Timer_AMReX_InputOutput )

    IF( iCycleChk .GT. 0 )THEN

      IF( MOD( StepNo(0), iCycleChk ) .EQ. 0 ) &
        chk = .TRUE.

    ELSE

      IF( ALL( t_new + dt .GT. t_chk ) )THEN

        t_chk = t_chk + dt_chk
        chk   = .TRUE.

      END IF

    END IF

    IF( chk )THEN

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() ) &
          WRITE(*,'(A)') 'CALL WriteCheckpointFile'

      END IF

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t_new, &
               MF_uGF % BA % P, &
               iWriteFields_uGF = 1, &
               iWriteFields_uCF = 1, &
               iWriteFields_uCR = 0, &
               pMF_uGF_Option = MF_uGF % P, &
               pMF_uCF_Option = MF_uCF % P )

      CALL FinalizeTimers_AMReX &
             ( RestartProgramTimer_Option = .TRUE. )

      chk = .FALSE.

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

  END SUBROUTINE WriteCheckpointFile


  SUBROUTINE ReGrid

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

    END IF

    IF( UseAMR )THEN

      IF( MOD( StepNo(0), iReGrid ) .EQ. 0 )THEN

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(6x,A,I2.2)') 'nLevels (before regrid): ', nLevels
            WRITE(*,'(6x,A)') 'Regridding'

          END IF

        END IF

        DO iLevel = 0, nLevels

          IF( iLevel .LT. nLevels-1 ) &
            CALL amrex_regrid( iLevel, t_new(iLevel) )

        END DO

        nLevels = amrex_get_numlevels()

        ! --- nLevels <= nMaxLevels; entire arrays t_old(0:nMaxLevels-1) and
        !     t_new(0:nMaxLevels-1) must have valid data ---
        t_old = t_old(0)
        t_new = t_new(0)

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,'(6x,A,I2.2)') 'nLevels (after regrid): ', nLevels
            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ApplyBoundaryConditions_Geometry_MF'

          END IF

        END IF

        CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ApplyPositivityLimiter_Euler_MF'
            WRITE(*,*)

          END IF

        END IF

        ! --- Regridding may cause some cells to be un-physical ---
        CALL ApplyPositivityLimiter_Euler_MF( MF_uGF, MF_uCF, MF_uDF )

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ComputeFromConserved_Euler_MF'
            WRITE(*,*)

          END IF

          CALL ComputeFromConserved_Euler_MF &
                 ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

        END IF

      END IF ! MOD( StepNo(0), 10 ) .EQ. 0

    END IF ! UseAMR

  END SUBROUTINE ReGrid


END PROGRAM main
