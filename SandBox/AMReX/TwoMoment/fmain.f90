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
    MF_ComputeTimeStep_Fancy,          &
    MF_ComputeFromConserved,           &      
    MF_ComputeFromConserved_Euler
  USE MF_UtilitiesModule,     ONLY: & 
    WriteNodalDataToFile, &
    WriteEulerToFile
  USE MyAmrDataModule,                  ONLY: &
    MF_uCR, &
    MF_uPR, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uGF
  USE MF_TwoMoment_TallyModule,         ONLY: &
    MF_ComputeTally_TwoMoment
  USE MF_Euler_TallyModule,         ONLY: &
    MF_ComputeTally_Euler
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
    dt_rel,    &
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
  USE InputOutputEuler,           ONLY: &
     WriteFieldsAMReX_PlotFile_Euler
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
  REAL(amrex_real) :: n
  
  n = 1.0_amrex_real
  CALL InitializeProgram
  
num = 1
  DO WHILE( ALL( t .LT. t_end ) )
    
    StepNo = StepNo + 1
    IF ( dt_rel .NE. 0.0_amrex_real ) THEN

      dt = dt_rel

    ELSE

      CALL MF_ComputeTimeStep_Fancy( MF_uGF, nX, nNodes, xR, xL, CFL, dt )

    END IF

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
    CALL MF_Update_IMEX_RK &
           ( t, dt, uGE, MF_uGF, MF_uCF, MF_uCR, GEOM, &
            Verbose_Option = amrex_parallel_ioprocessor()  )

      CALL MF_ComputeTally_Euler( GEOM, MF_uGF, MF_uCF, t(0), &
                                Verbose_Option = .FALSE. )

      CALL MF_ComputeTally_TwoMoment( GEOM, MF_uGF, MF_uCF, MF_uCR, &
                                    t(0), Verbose_Option = .FALSE. )
    IF( ALL( t + dt .GT. t_wrt ) )THEN
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( wrt )THEN
      
      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

      CALL WriteFieldsAMReX_PlotFile &
               ( t(0), StepNo, &
                 MF_uCR_Option = MF_uCR, &
                 MF_uPR_Option = MF_uPR, &
                 num_Option = num )


      CALL MF_ComputeFromConserved_Euler( MF_uGF, MF_uCF, MF_uPF, MF_uAF )


      CALL WriteFieldsAMReX_PlotFile_Euler &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF, &
               num_Option = num )

!      CALL MF_ComputeTally_Euler( GEOM, MF_uGF, MF_uCF, t(0), &
!                                Verbose_Option = .FALSE. )

!      CALL MF_ComputeTally_TwoMoment( GEOM, MF_uGF, MF_uCF, MF_uCR, &
!                                    t(0), Verbose_Option = .FALSE. )

      num = num + 1
      wrt = .FALSE.
    END IF


  END DO
 
  IF (nDOFZ .GT. 1) THEN

!    CALL WriteNodalDataToFile( GEOM, MF_uGF, MF_uCF, MF_uCR, 'thornado_')

  END IF

  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

  CALL WriteFieldsAMReX_Checkpoint & 
      ( StepNo, nLevels, dt, t, t_wrt, BA % P, &
        MF_uCR % P,  &
        MF_uPR % P  )

  CALL WriteFieldsAMReX_PlotFile &
           ( t(0), StepNo, &
             MF_uCR_Option = MF_uCR, &
             MF_uPR_Option = MF_uPR, &
             num_Option = num )

  CALL MF_ComputeFromConserved_Euler( MF_uGF, MF_uCF, MF_uPF, MF_uAF )


  CALL WriteFieldsAMReX_PlotFile_Euler &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF, &
               num_Option = num )

  CALL MF_ComputeTally_Euler( GEOM, MF_uGF, MF_uCF, t(0), &
                                Verbose_Option = .FALSE. )

  CALL MF_ComputeTally_TwoMoment( GEOM, MF_uGF, MF_uCF, MF_uCR, &
                                    t(0), Verbose_Option = .FALSE. )
  CALL FinalizeProgram( GEOM )
  

END PROGRAM main
