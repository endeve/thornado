PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- thornado Modules ---

  USE UnitsModule, ONLY: &
    UnitsDisplay

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_MHD, ONLY: &
    MF_uCM, &
    MF_uPM, &
    MF_uAM, &
    MF_uDM
  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram
  USE MF_MHD_UtilitiesModule, ONLY: &
    ComputeTimeStep_MHD_MF, &
    ComputeFromConserved_MHD_MF
  USE InputOutputModuleAMReX_MHD, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_MHD_TallyModule, ONLY: &
    ComputeTally_MHD_MF, &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    MHDMomentumX1_Initial, &
    MHDMomentumX1_OffGrid, &
    MHDMomentumX2_Initial, &
    MHDMomentumX2_OffGrid, &
    MHDMomentumX3_Initial, &
    MHDMomentumX3_OffGrid, &
    MHDEnergy_Initial, &
    MHDEnergy_OffGrid, &
    ElectronNumber_Initial, &
    ElectronNumber_OffGrid, &
    MHDMagFieldX1_Initial, &
    MHDMagFieldX1_OffGrid, &
    MHDMagFieldX2_Initial, &
    MHDMagFieldX2_OffGrid, &
    MHDMagFieldX3_Initial, &
    MHDMagFieldX3_OffGrid, &
    MHDCleaningField_Initial, &
    MHDCleaningField_OffGrid, &
    ADMMass_Initial, &
    ADMMass_OffGrid
  USE MF_TimeSteppingModule_SSPRK_MHD, ONLY: &
    UpdateFluid_SSPRK_MF, &
    CFL
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    t_end, &
    t_new, &
    t_old, &
    dt, &
    iCycleD, &
    iCycleW, &
    iCycleChk, &
    t_wrt, &
    t_chk, &
    dt_wrt, &
    dt_chk, &
    DEBUG
  USE MF_TimersModule_MHD, ONLY: &
    TimeIt_AMReX, &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_InputOutput, &
    FinalizeTimers_AMReX
  USE MF_MHD_TimersModule, ONLY: &
    TimeIt_AMReX_MHD
  USE ReGridModule_MHD, ONLY: &
    ReGrid

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: iErr
  LOGICAL  :: wrt, chk
  REAL(DP) :: Timer_Evolution

  TimeIt_AMReX       = .TRUE.
  TimeIt_AMReX_MHD = .TRUE.

  wrt = .FALSE.
  chk = .FALSE.

  CALL InitializeProgram

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

  ! --- Begin evolution ---

  DO WHILE( MAXVAL( t_new ) .LT. t_end )

    StepNo = StepNo + 1

    t_old = t_new

    CALL ReGrid

    CALL ComputeTimeStep_MHD_MF( MF_uGF, MF_uCM, CFL, dt )

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
        WRITE(*,'(A)') 'CALL ComputeFromConserved_MHD_MF'
        WRITE(*,*)

      END IF

      CALL ComputeFromConserved_MHD_MF &
             ( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

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

      CALL ComputeFromConserved_MHD_MF &
             ( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

      CALL WriteFieldsAMReX_PlotFile &
             ( t_new(0), StepNo, MF_uGF, &
               MF_uGF_Option = MF_uGF, &
               MF_uCM_Option = MF_uCM, &
               MF_uPM_Option = MF_uPM, &
               MF_uAM_Option = MF_uAM, &
               MF_uDM_Option = MF_uDM )

      CALL ComputeTally_MHD_MF &
             ( t_new, MF_uGF, MF_uCM, MF_uAM, Verbose_Option = .TRUE. )

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

      CALL ComputeFromConserved_MHD_MF &
             ( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t_new, &
               [ BaryonicMass_Initial   , BaryonicMass_OffGrid    ], &
               [ MHDMomentumX1_Initial, MHDMomentumX1_OffGrid ], &
               [ MHDMomentumX2_Initial, MHDMomentumX2_OffGrid ], &
               [ MHDMomentumX3_Initial, MHDMomentumX3_OffGrid ], &
               [ MHDEnergy_Initial    , MHDEnergy_OffGrid     ], &
               [ ElectronNumber_Initial , ElectronNumber_OffGrid  ], &
               [ MHDMagFieldX1_Initial, MHDMagFieldX1_OffGrid ], &
               [ MHDMagFieldX2_Initial, MHDMagFieldX2_OffGrid ], &
               [ MHDMagFieldX3_Initial, MHDMagFieldX3_OffGrid ], &
               [ MHDCleaningField_Initial, MHDCleaningField_OffGrid ], &
               [ ADMMass_Initial        , ADMMass_OffGrid         ], &
               MF_uGF % BA % P, &
               iWriteFields_uGF = 1, &
               iWriteFields_uCM = 1, &
               iWriteFields_uCR = 0, &
               pMF_uGF_Option = MF_uGF % P, &
               pMF_uCM_Option = MF_uCM % P )

      CALL FinalizeTimers_AMReX &
             ( RestartProgramTimer_Option = .TRUE., &
               Verbose_Option = amrex_parallel_ioprocessor() )

      chk = .FALSE.

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

  END SUBROUTINE WriteCheckpointFile


END PROGRAM main
