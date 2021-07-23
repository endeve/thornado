PROGRAM main

  ! --- AMReX Modules ---
  USE amrex_fort_module,                ONLY: &
    amrex_real
  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- Local Modules ---
  USE MF_TwoMoment_UtilitiesModule,     ONLY: & 
    MF_ComputeTimeStep,                &
    MF_ComputeTimeStep_Fancy,                &
    MF_ComputeFromConserved
  USE MF_UtilitiesModule,     ONLY: & 
    WriteNodalDataToFile
  USE MyAmrDataModule,                  ONLY: &
    MF_uCR, &
    MF_uPR, &
    MF_uCF, &
    MF_uGF
  USE InitializationModule,             ONLY: &
    wrt,              &
    InitializeProgram
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE MyAmrModule,                      ONLY: &
    nLevels,   &
    StepNo,    &
    t,         &
    dt,        &
    nX,        &
    xR,        &
    xL,        &
    nNodes,    &
    t_end,     &
    CFL,       &
    t_wrt,     &
    dt_wrt,    &
    t_chk,     &
    dt_chk,    &
    iCycleD,   &
    iCycleW,   &
    iCycleChk, &
    BA,        &
    GEOM
  USE ProgramHeaderModule,  ONLY: &
    nDOFZ
  USE MF_TwoMoment_TimeSteppingModule_Relativistic,      ONLY: &
    MF_Update_IMEX_RK

  ! --- thornado Modules ---
  USE InputOutput,           ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint, &
    ReadCheckpointFile
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE UnitsModule,            ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay, &
    Millisecond, &
    Kilometer

  IMPLICIT NONE

  REAL(amrex_real) :: n, m, dt_Fancy
  
  n = 1.0_amrex_real
  CALL InitializeProgram

!  CALL WriteFieldsAMReX_Checkpoint & 
!      ( StepNo, nLevels, dt, t, t_wrt, BA % P, &
!        MF_uCR % P,  &
!        MF_uPR % P  )
  DO WHILE( ALL( t .LT. t_end ) )
    
    StepNo = StepNo + 1
 
    CALL MF_ComputeTimeStep_Fancy( MF_uGF, nX, nNodes, xR, xL, CFL, dt )
 !   dt_Fancy = dt(0)
 !   CALL MF_ComputeTimeStep( nX, xR, xL, nNodes, CFL, dt )

 !   print*, dt(0) / UnitsDisplay % TimeUnit, dt_Fancy / UnitsDisplay % TimeUnit
!STOP
    IF( ALL( t + dt .LE. t_end ) )THEN
      t = t + dt
    ELSE
      dt = t_end - [t]
      t  = [t_end]
    END IF
    IF( amrex_parallel_ioprocessor() )THEN
      !WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
       print*,  'StepNo: ', StepNo(0), ' t = ', t / UnitsDisplay % TimeUnit , &
       TRIM( UnitsDisplay % TimeLabel ), ' dt = ', dt(0) / UnitsDisplay % TimeUnit, &
       TRIM( UnitsDisplay % TimeLabel )
    END IF
    !this is where the issue is
    CALL MF_Update_IMEX_RK &
           ( t, dt, uGE, MF_uGF, MF_uCF, MF_uCR, GEOM, &
            Verbose_Option = amrex_parallel_ioprocessor()  )

    IF( ALL( t + dt .GT. t_wrt ) )THEN
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( wrt )THEN

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

      CALL WriteFieldsAMReX_PlotFile &
               ( t(0), StepNo, &
                 MF_uCR_Option = MF_uCR, &
                 MF_uPR_Option = MF_uPR )
      wrt = .FALSE.
    END IF


!    IF (t(0) .GE. n) THEN
!
!  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )
!
!   CALL WriteFieldsAMReX_PlotFile &
!         ( t(0), StepNo, &
!             MF_uCR_Option = MF_uCR, &
!             MF_uPR_Option = MF_uPR )
!
!             n = n+1.0_amrex_real
!        
!    END IF
!
  END DO
 
  IF (nDOFZ .GT. 1) THEN

    CALL WriteNodalDataToFile( GEOM, MF_uGF, MF_uCF, MF_uCR, 'thornado_')

  END IF

  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

  CALL WriteFieldsAMReX_Checkpoint & 
      ( StepNo, nLevels, dt, t, t_wrt, BA % P, &
        MF_uCR % P,  &
        MF_uPR % P  )

  CALL WriteFieldsAMReX_PlotFile &
           ( t(0), StepNo, &
             MF_uCR_Option = MF_uCR, &
             MF_uPR_Option = MF_uPR )


  CALL FinalizeProgram( GEOM )
  

END PROGRAM main
