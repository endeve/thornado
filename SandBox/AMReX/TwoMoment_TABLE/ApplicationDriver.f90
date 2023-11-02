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
    ComputeTimeStep_TwoMoment_MF,                &
    ComputeTimeStep_TwoMoment_Fancy_MF,          &
    ComputeFromConserved_TwoMoment_MF,           &
    ComputeGray_TwoMoment_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
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
    MF_uCR, &
    MF_uGR
  USE MF_TwoMoment_TallyModule,         ONLY: &
    ComputeTally_TwoMoment_MF
  USE MF_Euler_TallyModule,         ONLY: &
    ComputeTally_Euler_MF, &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    Energy_Initial, &
    Energy_OffGrid, &
    ElectronNumber_Initial, &
    ElectronNumber_OffGrid, &
    ADMMass_Initial, &
    ADMMass_OffGrid
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
    Update_IMEX_RK_MF

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

      CALL ComputeTimeStep_TwoMoment_Fancy_MF &
             ( MF_uGF, nX, nNodes, xR, xL, CFL, dt )

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
    CALL Update_IMEX_RK_MF &
           ( t_new, dt, uGE, MF_uGF, MF_uCF, MF_uCR, amrex_geom, &
            Verbose_Option = amrex_parallel_ioprocessor()  )

      CALL ComputeTally_Euler_MF &
             ( t_new, MF_uGF, MF_uCF, Verbose_Option = .FALSE. )

      CALL ComputeTally_TwoMoment_MF( amrex_geom, MF_uGF, MF_uCF, MF_uCR, &
                                    t_new(0), Verbose_Option = .FALSE. )
    IF( ALL( t_new + dt .GT. t_wrt ) )THEN
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( wrt )THEN


     CALL ComputeFromConserved_TwoMoment_MF &
            ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

     CALL ComputeFromConserved_Euler_MF &
            ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

     CALL ComputeGray_TwoMoment_MF &
            ( MF_uGF, MF_uPF, MF_uCR, MF_uPR, MF_uGR )

     CALL WriteFieldsAMReX_Checkpoint &
            ( StepNo, nLevels, dt, t_new, &
              [ BaryonicMass_Initial  , BaryonicMass_OffGrid   ], &
              [ Energy_Initial        , Energy_OffGrid         ], &
              [ ElectronNumber_Initial, ElectronNumber_OffGrid ], &
              [ ADMMass_Initial       , ADMMass_OffGrid        ], &
              MF_uGF % BA % P, &
              iWriteFields_uGF = 1, &
              iWriteFields_uCF = 1, &
              iWriteFields_uCR = 1, &
              pMF_uGF_Option = MF_uGF % P, &
              pMF_uCF_Option = MF_uCF % P, &
              pMF_uCR_Option = MF_uCR % P )

     CALL WriteFieldsAMReX_PlotFile &
            ( t_new(0), StepNo, MF_uGF, &
              MF_uGF_Option = MF_uGF, &
              MF_uCF_Option = MF_uCF, &
              MF_uPF_Option = MF_uPF, &
              MF_uAF_Option = MF_uAF, &
              MF_uDF_Option = MF_uDF, &
              MF_uPR_Option = MF_uPR, &
              MF_uCR_Option = MF_uCR, &
              MF_uGR_Option = MF_uGR )

      num = num + 1
      wrt = .FALSE.
    END IF


  END DO

  IF (nDOFZ .GT. 1) THEN

!!$  CALL ShowVariableFromMultiFab_Vector &
!!$    ( MF_uCR, 1, WriteToFile_Option = .TRUE., &
!!$       FileName_Option = 'thornado_uCR_' )

  END IF


  CALL ComputeFromConserved_TwoMoment_MF( MF_uGF, MF_uCF, MF_uCR, MF_uPR )


  CALL ComputeFromConserved_Euler_MF( MF_uGF, MF_uCF, MF_uPF, MF_uAF )



  CALL ComputeGray_TwoMoment_MF &
         ( MF_uGF, MF_uPF, MF_uCR, MF_uPR, MF_uGR )

  CALL WriteFieldsAMReX_PlotFile &
         ( t_new(0), StepNo, MF_uGF, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF, &
           MF_uDF_Option = MF_uDF, &
           MF_uPR_Option = MF_uPR, &
           MF_uCR_Option = MF_uCR, &
           MF_uGR_Option = MF_uGR )


  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t_new, &
           BaryonicMass_Initial, &
           Energy_Initial, &
           ElectronNumber_Initial, &
           ADMMass_Initial, &
           MF_uGF % BA % P, &
           iWriteFields_uGF = 1, &
           iWriteFields_uCF = 1, &
           iWriteFields_uCR = 1, &
           pMF_uGF_Option = MF_uGF % P, &
           pMF_uCF_Option = MF_uCF % P, &
           pMF_uCR_Option = MF_uCR % P )



  CALL ComputeTally_Euler_MF &
         ( t_new, MF_uGF, MF_uCF, Verbose_Option = .FALSE. )

  CALL ComputeTally_TwoMoment_MF( amrex_geom, MF_uGF, MF_uCF, MF_uCR, &
                                    t_new(0), Verbose_Option = .FALSE. )
  CALL FinalizeProgram


END PROGRAM main
