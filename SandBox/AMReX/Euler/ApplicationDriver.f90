PROGRAM ApplicationDriver

  ! --- AMReX Modules ---

  USE amrex_amrcore_module, ONLY: &
    amrex_max_level
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
    DP
  USE MF_FieldsModule, ONLY: &
    MF_uGF_new, &
    MF_uCF_new, &
    MF_uPF_new, &
    MF_uAF_new, &
    MF_uDF_new
  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeTimeStep_Euler_MF, &
    ComputeFromConserved_Euler_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile
  USE MF_Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_MF
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    UpdateFluid_SSPRK_MF
  USE InputParsingModule, ONLY: &
    StepNo, &
    t_end, &
    t_new, &
    t_old, &
    dt, &
    CFL, &
    iCycleD, &
    iCycleW, &
    iCycleChk, &
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

  INTEGER  :: iLevel
  LOGICAL  :: wrt, chk
  REAL(DP) :: t_wrt, t_chk
  REAL(DP) :: Timer_Evolution

  TimeIt_AMReX_Euler = .TRUE.

  TimeIt_Euler = .TRUE.

  wrt = .FALSE.
  chk = .FALSE.

  t_wrt = dt_wrt
  t_chk = dt_chk

  CALL InitializeProgram

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

  DO WHILE( ALL( t_new .LT. t_end ) )

    StepNo = StepNo + 1

    CALL ComputeTimeStep_Euler_MF( MF_uGF_new, MF_uCF_new, CFL, dt )

    IF( ANY( t_new + dt .LE. t_end ) )THEN

      t_new = t_new + dt

    ELSE

      DO iLevel = 0, amrex_max_level

        dt   (iLevel) = t_end - t_new(iLevel)
        t_new(iLevel) = t_end

      END DO

    END IF

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    IF( amrex_parallel_ioprocessor() )THEN

      IF( MOD( StepNo(0), iCycleD ) .EQ. 0 )THEN

        WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
          'StepNo: ', StepNo(0), ' t = ', t_new(0) / UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel ), &
          ' dt = ', dt(0) /  UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel )

     END IF

    END IF

    IF( iCycleChk .GT. 0 )THEN

      IF( MOD( StepNo(0), iCycleChk ) .EQ. 0 ) &
        chk = .TRUE.

    ELSE

      IF( ALL( t_new + dt .GT. t_chk ) )THEN

        t_chk = t_chk + dt_chk
        chk   = .TRUE.

      END IF

    END IF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL UpdateFluid_SSPRK_MF &
          ( t_new, dt, MF_uGF_new, MF_uCF_new, MF_uDF_new, &
            ComputeIncrement_Euler_MF )

    IF( chk )THEN

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF_new, MF_uCF_new, MF_uPF_new, MF_uAF_new )

!      CALL WriteFieldsAMReX_Checkpoint &
!             ( StepNo, nLevels, dt, t, &
!               MF_uGF % BA % P, &
!               MF_uGF % P, &
!               MF_uCF % P )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

      CALL FinalizeTimers_Euler &
             ( Verbose_Option = amrex_parallel_ioprocessor(), &
               SuppressApplicationDriver_Option = .TRUE., &
               WriteAtIntermediateTime_Option = .TRUE. )

      CALL FinalizeTimers_AMReX_Euler &
             ( WriteAtIntermediateTime_Option = .TRUE. )

      chk = .FALSE.

    END IF

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
             ( MF_uGF_new, MF_uCF_new, MF_uPF_new, MF_uAF_new )

      CALL WriteFieldsAMReX_PlotFile &
             ( t_new(0), StepNo, MF_uGF_new, &
               MF_uGF_Option = MF_uGF_new, &
               MF_uCF_Option = MF_uCF_new, &
               MF_uPF_Option = MF_uPF_new, &
               MF_uAF_Option = MF_uAF_new, &
               MF_uDF_Option = MF_uDF_new )

      wrt = .FALSE.

    END IF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

  END DO

  ! --- END of evolution ---

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(2x,A,ES13.6E3,A)') &
      'Total evolution time: ', MPI_WTIME() - Timer_Evolution, ' s'

  END IF

  StepNo = StepNo + 1

  CALL FinalizeProgram

END PROGRAM ApplicationDriver
