PROGRAM main

  ! --- AMReX Modules ---
  USE amrex_fort_module, ONLY: &
    amrex_real
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- thornado Modules ---
  USE MeshModule,                       ONLY: &
    MeshX, DestroyMesh
  USE InputOutputModuleAMReX,           ONLY: &
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile!, &
!!$    MakeMF_Diff

  ! --- Local Modules ---
  USE MF_Euler_UtilitiesModule,         ONLY: &
    MF_ComputeFromConserved, &
    MF_ComputeTimeStep
  USE MF_Euler_SlopeLimiterModule,      ONLY: &
    MF_Euler_ApplySlopeLimiter
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    MF_Euler_ApplyPositivityLimiter
  USE MF_Euler_dgDiscretizationModule,  ONLY: &
    MF_Euler_ComputeIncrement
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_UpdateFluid_SSPRK
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE MF_UtilitiesModule,               ONLY: &
    ShowVariableFromMultifab
  USE MyAmrDataModule,                  ONLY: &
    MF_uGF, MF_uCF, MF_uPF, MF_uAF
  USE InitializationModule
  USE MyAmrModule

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER          :: iErr
  REAL(amrex_real) :: Timer_Evolution

!!$  CALL MakeMF_Diff( 0, 2929 )

  ! --- Argument is integer corresponding to
  !     checkpoint file from which to restart ---
  CALL InitializeProblem()

  IF( amrex_parallel_ioprocessor() ) &
    Timer_Evolution = MPI_WTIME()

  CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

  DO WHILE( ALL( t .LT. t_end ) )

    StepNo = StepNo + 1

    CALL MF_ComputeTimeStep( MF_uGF, MF_uCF, CFL, dt )

    IF( ALL( t + dt .LE. t_end ) )THEN
      t = t + dt
    ELSE
      dt = t_end - [t]
      t  = [t_end]
    END IF

    IF( amrex_parallel_ioprocessor() )THEN
      IF( MOD( StepNo(0), iCycleD ) .EQ. 0 ) &
        WRITE(*,'(A5,A,I6.6,A,ES13.6E3,A,ES13.6E3)') &
          '', 'StepNo: ', StepNo(0), ', t = ', t, ', dt = ', dt(0)
    END IF

    CALL MF_UpdateFluid_SSPRK &
           ( t, dt, MF_uGF, MF_uCF, &
             GEOM, MF_Euler_ComputeIncrement )

    IF( iCycleW .GT. 0 )THEN
      IF( MOD( StepNo(0), iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.
    ELSE
      IF( ALL( t + dt .GT. t_wrt ) )THEN
        t_wrt = t_wrt + dt_wrt
        wrt   = .TRUE.
      END IF
    END IF

    IF( wrt )THEN

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF )

      wrt = .FALSE.

    END IF

    IF( iCycleChk .GT. 0 )THEN
      IF( MOD( StepNo(0), iCycleChk ) .EQ. 0 ) &
        chk = .TRUE.
    ELSE
      IF( ALL( t + dt .GT. t_chk ) )THEN
        t_chk = t_chk + dt_chk
        chk   = .TRUE.
      END IF
    END IF

    IF( chk )THEN

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t, &
               MF_uGF % BA % P, &
               MF_uGF % P, &
               MF_uCF % P, &
               MF_uPF % P, &
               MF_uAF % P )

      chk = .FALSE.
     
    END IF

  END DO

  ! --- END of evolution ---

  IF( amrex_parallel_ioprocessor() )THEN
    WRITE(*,*)
    WRITE(*,'(A,ES13.6E3,A)') &
      'Total evolution time: ', MPI_WTIME() - Timer_Evolution, ' s'
  END IF

  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

  StepNo = StepNo + 1
  CALL WriteFieldsAMReX_PlotFile &
         ( t(0), StepNo, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t, &
           MF_uGF % BA % P, &
           MF_uGF % P, &
           MF_uCF % P, &
           MF_uPF % P, &
           MF_uAF % P )

  ! --- Finalize everything ---

  CALL FinalizeProgram( GEOM, MeshX )


END PROGRAM main

