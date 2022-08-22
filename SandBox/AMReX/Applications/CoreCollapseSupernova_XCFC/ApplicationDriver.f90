PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_amrcore_module, ONLY: &
   amrex_regrid, &
   amrex_get_numlevels, &
   amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodes
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE TimersModule_Euler, ONLY: &
    TimeIt_Euler, &
    FinalizeTimers_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Two
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR, &
    MF_uPR
  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeTimeStep_Euler_MF, &
    ComputeFromConserved_Euler_MF
  USE MF_TwoMoment_UtilitiesModule, ONLY: &
    ComputeTimeStep_TwoMoment_Fancy_MF, &
    ComputeFromConserved_TwoMoment_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_TimeSteppingModule, ONLY: &
    Update_IMEX_RK_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    ComputeTally_Euler_MF
  USE MF_TwoMoment_TallyModule, ONLY: &
    ComputeTally_TwoMoment_MF
  USE AverageDownModule, ONLY: &
    AverageDownTo
  USE InputParsingModule, ONLY: &
    nLevels, &
    StepNo, &
    t_end, &
    t_new, &
    t_old, &
    dt, &
    dt_TM, &
    CFL, &
    iCycleD, &
    iCycleW, &
    iCycleChk, &
    t_wrt, &
    t_chk, &
    dt_wrt, &
    dt_chk, &
    UseAMR, &
    DEBUG, &
    nX, &
    xL, &
    xR
  USE MF_Euler_TimersModule, ONLY: &
    TimeIt_AMReX_Euler, &
    FinalizeTimers_AMReX_Euler, &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_InputOutput

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: iLevel, iErr
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

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

    END IF

    IF( UseAMR )THEN

      IF( MOD( StepNo(0), 10 ) .EQ. 0 )THEN

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

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,'(A)') 'CALL ComputeTimeStep_Euler_MF'

    END IF

    CALL ComputeTimeStep_Euler_MF( MF_uGF, MF_uCF, CFL, dt )
    CALL ComputeTimeStep_TwoMoment_Fancy_MF &
           ( MF_uGF, nX, nNodes, xR, xL, CFL, dt_TM )

    dt = MIN( MINVAL( dt ), MINVAL( dt_TM ) )

    IF( MAXVAL( t_old + dt ) .LT. t_end )THEN

      t_new = t_old + dt

    ELSE

      dt = t_end - t_old

      t_new = t_end

    END IF

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,'(A)') 'CALL UpdateFluid_IMEX_RK_MF'

    END IF

    CALL Update_IMEX_RK_MF

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

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() ) &
          WRITE(*,'(A)') 'CALL WritePlotFile'

      END IF

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL ComputeFromConserved_TwoMoment_MF &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

      CALL WriteFieldsAMReX_PlotFile &
             ( t_new(0), StepNo, MF_uGF, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF, &
               MF_uDF_Option = MF_uDF, &
               MF_uCR_Option = MF_uCR, &
               MF_uPR_Option = MF_uPR )

      CALL ComputeTally_Euler_MF &
             ( t_new, MF_uGF, MF_uCF, &
               Verbose_Option = .FALSE. )
               !Verbose_Option = amrex_parallel_ioprocessor() )

      CALL ComputeTally_TwoMoment_MF &
             ( amrex_geom, MF_uGF, MF_uCF, MF_uCR, t_new(0), &
               Verbose_Option = .FALSE. )
               !Verbose_Option = amrex_parallel_ioprocessor() )

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

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() ) &
          WRITE(*,'(A)') 'CALL WriteCheckpointFile'

      END IF

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL ComputeFromConserved_TwoMoment_MF &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t_new, &
               MF_uGF % BA % P, &
               iWriteFields_uGF = 1, &
               iWriteFields_uCF = 1, &
               iWriteFields_uCR = 1, &
               pMF_uGF_Option = MF_uGF % P, &
               pMF_uCF_Option = MF_uCF % P, &
               pMF_uCR_Option = MF_uCR % P )

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


END PROGRAM main