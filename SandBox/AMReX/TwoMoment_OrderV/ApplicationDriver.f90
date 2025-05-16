PROGRAM main

  ! --- AMReX Modules ---
  USE amrex_fort_module,                ONLY: &
    amrex_real
  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---
  USE UnitsModule,            ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay
  USE GeometryFieldsModuleE, ONLY: &
    uGE

  ! --- Local Modules ---
  USE MF_TwoMoment_UtilitiesModule,     ONLY: &
    ComputeFromConserved_TwoMoment_MF, &
    ComputeTimeStep_TwoMoment_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_FieldsModule_Geometry,                  ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler,                  ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE MF_FieldsModule_TwoMoment,                  ONLY: &
    MF_uCR, &
    MF_uPR, &
    MF_uAR, &
    MF_uGR
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
    t_wrt,     &
    dt_wrt,    &
    dt_rel
  USE MF_Euler_TallyModule, ONLY: &
    ComputeTally_Euler_MF, &
    FinalizeTally_Euler_MF, &
    BaryonicMass_Initial, &
    BaryonicMass_OffGrid, &
    EulerMomentumX1_Initial, &
    EulerMomentumX1_OffGrid, &
    EulerMomentumX2_Initial, &
    EulerMomentumX2_OffGrid, &
    EulerMomentumX3_Initial, &
    EulerMomentumX3_OffGrid, &
    EulerEnergy_Initial, &
    EulerEnergy_OffGrid, &
    ElectronNumber_Initial, &
    ElectronNumber_OffGrid, &
    ADMMass_Initial, &
    ADMMass_OffGrid
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Update_IMEX_RK_MF, &
    CFL
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab

  IMPLICIT NONE

  LOGICAL  :: wrt, chk
  REAL(amrex_real) :: n

  n = 1.0_amrex_real
  CALL InitializeProgram

  CALL ShowVariableFromMultifab(MF_uPR, 1, writetofile_option=.TRUE., FileNameBase_Option ='Primitive_Variables')

  CALL ComputeFromConserved_TwoMoment_MF(  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )

  CALL ShowVariableFromMultifab(MF_uCR, 1, writetofile_option=.TRUE., FileNameBase_Option ='Conserved_Variables')

  DO WHILE( ALL( t_new .LT. t_end ) )

    StepNo = StepNo + 1
    IF ( dt_rel .NE. 0.0_amrex_real ) THEN

      dt = dt_rel

    ELSE

      CALL ComputeTimeStep_TwoMoment_MF( MF_uGF, CFL, dt )
  PRINT *, 'MULTIFAB DATA START FROM Application driver'
                    CALL Showvariablefrommultifab(MF_uGF,1)
PRINT *, 'MULTIFAB DATA END'
      dt = MINVAL( dt )
    END IF


    IF( ALL( t_new + dt .LE. t_end ) )THEN
      t_new = t_new + dt
    ELSE
      dt = t_end - t_new
      t_new  = t_end
    END IF
    IF( amrex_parallel_ioprocessor() )THEN
      !WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
       print*,  'StepNo: ', StepNo(0), ' t = ', t_new/ UnitsDisplay % TimeUnit , &
       TRIM( UnitsDisplay % TimeLabel ), ' dt = ', dt(0) / UnitsDisplay % TimeUnit, &
       TRIM( UnitsDisplay % TimeLabel )
    END IF


    CALL Update_IMEX_RK_MF

    IF( ALL( t_new + dt .GT. t_wrt ) )THEN
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF
      !wrt  = .TRUE.
    IF( wrt )THEN
        CALL ComputeFromConserved_TwoMoment_MF(  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )

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
        
        CALL ShowVariableFromMultifab(MF_uPR, 1, writetofile_option=.TRUE., FileNameBase_Option ='Primitive_Variables')
        CALL ShowVariableFromMultifab(MF_uCR, 1, writetofile_option=.TRUE., FileNameBase_Option ='Conserved_Variables')

    CALL WriteFieldsAMReX_Checkpoint &
           ( StepNo, nLevels, dt, t_new, &
             [ BaryonicMass_Initial   , BaryonicMass_OffGrid    ], &
             [ EulerMomentumX1_Initial, EulerMomentumX1_OffGrid ], &
             [ EulerMomentumX2_Initial, EulerMomentumX2_OffGrid ], &
             [ EulerMomentumX3_Initial, EulerMomentumX3_OffGrid ], &
             [ EulerEnergy_Initial    , EulerEnergy_OffGrid     ], &
             [ ElectronNumber_Initial , ElectronNumber_OffGrid  ], &
             [ ADMMass_Initial        , ADMMass_OffGrid         ], &
             MF_uGF % BA % P, &
             iWriteFields_uGF = 1, &
             iWriteFields_uCF = 1, &
             iWriteFields_uCR = 0, &
             pMF_uGF_Option = MF_uGF % P, &
             pMF_uCF_Option = MF_uCF % P)
             !pMF_uCR_Option = MF_uCR % P )
    END IF

  END DO
  


  CALL FinalizeProgram


END PROGRAM main
