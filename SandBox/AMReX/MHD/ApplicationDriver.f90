PROGRAM ApplicationDriver

  ! --- AMReX Modules ---

  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- thornado Modules ---

  USE InputOutputModuleAMReX_MHD,           ONLY: &
    WriteFieldsAMReX_Checkpoint, &
    WriteFieldsAMReX_PlotFile
  USE UnitsModule,                      ONLY: &
    UnitsDisplay

  ! --- Local Modules ---

  USE MF_KindModule,                    ONLY: &
    DP
  USE MF_UtilitiesModule,               ONLY: &
    WriteNodalDataToFile
  USE MF_MHD_UtilitiesModule,           ONLY: &
    MF_ComputeFromConserved, &
    MF_ComputeTimeStep
  USE MF_MHD_dgDiscretizationModule,    ONLY: &
    MF_ComputeIncrement_MHD
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_UpdateMagnetofluid_SSPRK
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE MF_FieldsModule_MHD,           ONLY: &
    MF_uGF, &
    MF_uCM, &
    MF_uPM, &
    MF_uAM, &
    MF_uDM
  USE InitializationModule,             ONLY: &
    InitializeProgram, &
    chk,               &
    wrt
  USE InputParsingModule,               ONLY: &
    nLevels,   &
    StepNo,    &
    t,         &
    dt,        &
    t_end,     &
    CFL,       &
    t_wrt,     &
    dt_wrt,    &
    t_chk,     &
    dt_chk,    &
    iCycleD,   &
    iCycleW,   &
    iCycleChk, &
    GEOM,      &
    WriteNodalData

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: iErr
  REAL(DP) :: Timer_Evolution

  CHARACTER(LEN=6) :: NodalFileBaseName

  CALL InitializeProgram

  IF( WriteNodalData )THEN

    WRITE( NodalFileBaseName, '(I6.6)' ) StepNo

    CALL WriteNodalDataToFile( GEOM, MF_uGF, MF_uCM, MF_uDM, NodalFileBaseName )

  END IF

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

  CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

  DO WHILE( ALL( t .LT. t_end ) )

    StepNo = StepNo + 1

    CALL MF_ComputeTimeStep( MF_uGF, MF_uCM, CFL, dt )

    IF( ALL( t + dt .LE. t_end ) )THEN

      t = t + dt

    ELSE

      dt = t_end - [t]
      t  = [t_end]

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      IF( MOD( StepNo(0), iCycleD ) .EQ. 0 )THEN

        WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
          'StepNo: ', StepNo(0), ' t = ', t / UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel ), &
          ' dt = ', dt(0) /  UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel )

     END IF

    END IF

    CALL MF_UpdateMagnetofluid_SSPRK &
           ( t, dt, MF_uGF, MF_uCM, MF_uDM, GEOM, MF_ComputeIncrement_MHD )

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

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t, &
               MF_uGF % BA % P, &
               MF_uGF % P, &
               MF_uCM % P )

      chk = .FALSE.

    END IF

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

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCM_Option = MF_uCM, &
               MF_uPM_Option = MF_uPM, &
               MF_uAM_Option = MF_uAM, &
               MF_uDM_Option = MF_uDM )

      IF( WriteNodalData )THEN

        WRITE( NodalFileBaseName, '(I6.6)' ) StepNo

        CALL WriteNodalDataToFile( GEOM, MF_uGF, MF_uCM, MF_uDM, NodalFileBaseName )

      END IF

      wrt = .FALSE.

    END IF

  END DO

  ! --- END of evolution ---

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(2x,A,ES13.6E3,A)') &
      'Total evolution time: ', MPI_WTIME() - Timer_Evolution, ' s'

  END IF

  StepNo = StepNo + 1

  CALL MF_ComputeFromConserved( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t, &
           MF_uGF % BA % P, &
           MF_uGF % P, &
           MF_uCM % P )

  CALL WriteFieldsAMReX_PlotFile &
         ( t(0), StepNo,           &
           MF_uGF_Option = MF_uGF, &
           MF_uCM_Option = MF_uCM, &
           MF_uPM_Option = MF_uPM, &
           MF_uAM_Option = MF_uAM, &
           MF_uDM_Option = MF_uDM )

  IF( WriteNodalData )THEN

    WRITE( NodalFileBaseName, '(I6.6)' ) StepNo

    CALL WriteNodalDataToFile( GEOM, MF_uGF, MF_uCM, MF_uDM, NodalFileBaseName )

  END IF

  ! --- Finalize everything ---

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(2x,A)') 'git info'
    WRITE(*,'(2x,A)') '--------'
    WRITE(*,*)
    WRITE(*,'(2x,A)') 'git branch:'
    CALL EXECUTE_COMMAND_LINE( 'git branch' )
    WRITE(*,*)
    WRITE(*,'(2x,A)') 'git describe --tags:'
    CALL EXECUTE_COMMAND_LINE( 'git describe --tags' )
    WRITE(*,*)
    WRITE(*,'(2x,A)') 'git rev-parse HEAD:'
    CALL EXECUTE_COMMAND_LINE( 'git rev-parse HEAD' )
    WRITE(*,*)
    WRITE(*,'(2x,A)') 'date:'
    CALL EXECUTE_COMMAND_LINE( 'date' )
    WRITE(*,*)

  END IF

  CALL FinalizeProgram( GEOM )

END PROGRAM ApplicationDriver
