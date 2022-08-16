PROGRAM main

  ! --- AMReX Modules ---
  USE amrex_fort_module,                ONLY: &
    amrex_real
  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- Local Modules ---
  USE MF_TwoMoment_UtilitiesModule,     ONLY: &
    MF_ComputeTimeStep,                &
    MF_ComputeTimeStep_Fancy,          &
    MF_ComputeFromConserved,           &
    MF_ComputeFromConserved_Euler
  USE MF_UtilitiesModule,     ONLY: &
    ShowVariableFromMultiFab
  USE MF_FieldsModule_Geometry,                  ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler,                  ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE MF_FieldsModule_TwoMoment,                  ONLY: &
    MF_uPR, &
    MF_uCR
  USE MF_TwoMoment_TallyModule,         ONLY: &
    MF_ComputeTally_TwoMoment
  USE MF_Euler_TallyModule,         ONLY: &
    ComputeTally_Euler_MF
  USE InitializationModule,             ONLY: &
    InitializeProgram
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE InputParsingModule,                      ONLY: &
    nLevels,   &
    StepNo,    &
    t_new,     &
    dt,        &
    nX,        &
    xR,        &
    xL,        &
    nNodes,    &
    t_end,     &
    CFL,       &
    t_wrt,     &
    dt_wrt,    &
    dt_rel
  USE ProgramHeaderModule,  ONLY: &
    nDOFZ
  USE MF_TwoMoment_TimeSteppingModule_Relativistic,      ONLY: &
    MF_Update_IMEX_RK

  ! --- thornado Modules ---
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint, &
    ReadCheckpointFile
!!$  USE InputOutputEuler,           ONLY: &
!!$     WriteFieldsAMReX_PlotFile_Euler
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE UnitsModule,            ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay, &
    Millisecond, &
    Kilometer

  IMPLICIT NONE

  INTEGER :: num
  LOGICAL  :: wrt, chk
  REAL(amrex_real) :: n

  n = 1.0_amrex_real
  CALL InitializeProgram

num = 1
  DO WHILE( ALL( t_new .LT. t_end ) )

    StepNo = StepNo + 1
    IF ( dt_rel .NE. 0.0_amrex_real ) THEN

      dt = dt_rel

    ELSE

      CALL MF_ComputeTimeStep_Fancy( MF_uGF, nX, nNodes, xR, xL, CFL, dt )

    END IF

    IF( ALL( t_new + dt .LE. t_end ) )THEN
      t_new = t_new + dt
    ELSE
      dt = t_end - [t_new]
      t_new  = [t_end]
    END IF
    IF( amrex_parallel_ioprocessor() )THEN
      !WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
       print*,  'StepNo: ', StepNo(0), ' t = ', t_new/ UnitsDisplay % TimeUnit , &
       TRIM( UnitsDisplay % TimeLabel ), ' dt = ', dt(0) / UnitsDisplay % TimeUnit, &
       TRIM( UnitsDisplay % TimeLabel )
    END IF
    CALL MF_Update_IMEX_RK &
           ( t_new, dt, uGE, MF_uGF, MF_uCF, MF_uCR, amrex_geom, &
            Verbose_Option = amrex_parallel_ioprocessor()  )

      CALL ComputeTally_Euler_MF &
             ( t_new, MF_uGF, MF_uCF, Verbose_Option = .FALSE. )

      CALL MF_ComputeTally_TwoMoment( amrex_geom, MF_uGF, MF_uCF, MF_uCR, &
                                    t_new(0), Verbose_Option = .FALSE. )
    IF( ALL( t_new + dt .GT. t_wrt ) )THEN
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( wrt )THEN

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

      CALL MF_ComputeFromConserved_Euler( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_PlotFile &
               ( t_new(0), StepNo, MF_uGF, &
                 MF_uCR_Option = MF_uCR, &
                 MF_uPR_Option = MF_uPR, &
                 PlotFileNumber_Option = num )




!!$      CALL WriteFieldsAMReX_PlotFile_Euler &
!!$             ( t_new(0), StepNo, &
!!$               MF_uGF_Option = MF_uGF, &
!!$               MF_uCF_Option = MF_uCF, &
!!$               MF_uPF_Option = MF_uPF, &
!!$               MF_uAF_Option = MF_uAF, &
!!$               num_Option = num )

!      CALL ComputeTally_Euler_MF &
!             ( t_new, MF_uGF, MF_uCF, Verbose_Option = .FALSE. )

!      CALL MF_ComputeTally_TwoMoment( amrex_geom, MF_uGF, MF_uCF, MF_uCR, &
!                                    t_new(0), Verbose_Option = .FALSE. )

      num = num + 1
      wrt = .FALSE.
    END IF


  END DO

  IF (nDOFZ .GT. 1) THEN

!!$  CALL ShowVariableFromMultiFab_Vector &
!!$    ( MF_uCR, 1, WriteToFile_Option = .TRUE., &
!!$       FileName_Option = 'thornado_uCR_' )

  END IF

  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t_new, &
           MF_uGF % BA % P, &
           iWriteFields_uGF = 1, &
           iWriteFields_uCF = 0, &
           iWriteFields_uCR = 1, &
           pMF_uGF_Option = MF_uGF % P, &
           pMF_uCR_Option = MF_uCR % P )

  CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uCR_Option = MF_uCR, &
             MF_uPR_Option = MF_uPR, &
             PlotFileNumber_Option = num )

  CALL MF_ComputeFromConserved_Euler( MF_uGF, MF_uCF, MF_uPF, MF_uAF )


!!$  CALL WriteFieldsAMReX_PlotFile_Euler &
!!$             ( t_new(0), StepNo, MF_uGF, &
!!$               MF_uGF_Option = MF_uGF, &
!!$               MF_uCF_Option = MF_uCF, &
!!$               MF_uPF_Option = MF_uPF, &
!!$               MF_uAF_Option = MF_uAF, &
!!$               num_Option = num )

  CALL ComputeTally_Euler_MF &
         ( t_new, MF_uGF, MF_uCF, Verbose_Option = .FALSE. )

  CALL MF_ComputeTally_TwoMoment( amrex_geom, MF_uGF, MF_uCF, MF_uCR, &
                                    t_new(0), Verbose_Option = .FALSE. )
  CALL FinalizeProgram


END PROGRAM main
