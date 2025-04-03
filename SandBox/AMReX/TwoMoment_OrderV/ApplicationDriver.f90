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
    MF_uAF, &
    MF_uCF
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
  USE MF_TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Update_IMEX_RK_MF, &
    CFL
    !Initialize_IMEX_RK_MF, &
    !Finalize_IMEX_RK_MF

  USE MF_UtilitiesModule, ONLY: &
  ShowVariableFromMultiFab


  IMPLICIT NONE

  LOGICAL  :: wrt, chk
  REAL(amrex_real) :: n

  n = 1.0_amrex_real
  CALL InitializeProgram

  CALL ComputeFromConserved_TwoMoment_MF(  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )

  DO WHILE( ALL( t_new .LT. t_end ) )

    StepNo = StepNo + 1
    IF ( dt_rel .NE. 0.0_amrex_real ) THEN

      dt = dt_rel

    ELSE

    PRINT *, CFL

    PRINT *, 'CALLING TIME STEPPER'

      CALL ComputeTimeStep_TwoMoment_MF( MF_uGF, CFL, dt )
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

  END DO
  


  CALL FinalizeProgram


END PROGRAM main
