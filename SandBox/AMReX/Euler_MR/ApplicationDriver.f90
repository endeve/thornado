PROGRAM ApplicationDriver

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    FinalizeTimers_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Two
  USE MF_FieldsModule, ONLY: &
    MF_uGF, &
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
    dt_chk
  USE MF_Euler_TimersModule, ONLY: &
    TimeIt_AMReX_Euler, &
    FinalizeTimers_AMReX_Euler, &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_InputOutput

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL  :: wrt, chk
  REAL(DP) :: Timer_Evolution

  TimeIt_AMReX_Euler = .TRUE.

  TimeIt_Euler = .TRUE.

  wrt = .FALSE.
  chk = .FALSE.

  CALL InitializeProgram

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

  ! --- Begin evolution ---

  DO WHILE( MAXVAL( t_new ) .LT. t_end )

    StepNo = StepNo + 1

    t_old = t_new

    CALL ComputeTimeStep_Euler_MF( MF_uGF, MF_uCF, CFL, dt )

    dt = MINVAL( dt )

    IF( MAXVAL( t_old + dt ) .LT. t_end )THEN

      t_new = t_old + dt

    ELSE

      dt = t_end - t_old

      t_new = t_end

    END IF

    CALL UpdateFluid_SSPRK_MF( t_new, dt, MF_uGF, MF_uCF, MF_uDF )

    IF( amrex_parallel_ioprocessor() )THEN

      IF( MOD( StepNo(0), iCycleD ) .EQ. 0 )THEN

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

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

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

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

  END SUBROUTINE WritePlotFile


  SUBROUTINE WriteCheckpointFile

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

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

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t_new, &
               MF_uGF % BA % P, &
               MF_uGF % P, &
               MF_uCF % P )

      CALL FinalizeTimers_Euler &
             !( Verbose_Option = amrex_parallel_ioprocessor(), &
             ( Verbose_Option = .FALSE., &
               SuppressApplicationDriver_Option = .TRUE., &
               WriteAtIntermediateTime_Option = .FALSE. )

      CALL FinalizeTimers_AMReX_Euler &
             ( WriteAtIntermediateTime_Option = .FALSE. )

      chk = .FALSE.

    END IF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

  END SUBROUTINE WriteCheckpointFile

END PROGRAM ApplicationDriver
