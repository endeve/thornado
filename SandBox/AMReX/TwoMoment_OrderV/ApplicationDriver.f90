PROGRAM main

  ! --- AMReX Modules ---
  USE amrex_fort_module,                ONLY: &
    amrex_real
  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---
  USE UnitsModule,            ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE KindModule, ONLY: &
    DP, Zero, One, Two

  ! --- Local Modules ---
  USE MF_TwoMoment_UtilitiesModule_OrderV,     ONLY: &
    ComputeFromConserved_TwoMoment_MF, &
    ComputeTimeStep_TwoMoment_MF, &
    ComputeTimeStep_TwoMoment_Realizability_MF
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
    t_chk,     &
    dt_wrt,    &
    dt_chk,    &
    iCycleChk, &
    iCycleW,   &
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
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
  ApplyBoundaryConditions_Euler_MF
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab
  USE MF_TimersModule, ONLY: &
    TimeIt_AMReX, &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_InputOutput, &
    FinalizeTimers_AMReX
  USE MF_TwoMoment_TallyModule
    
  IMPLICIT NONE

  LOGICAL  :: wrt, chk, UseRealizabilityTimeStep=.FALSE.
  TYPE(amrex_parmparse) :: PP
  CHARACTER(:), ALLOCATABLE :: ProgramName
  REAL(DP) :: N_before

  CALL InitializeProgram

  !CALL ShowVariableFromMultiFab( MF_uCR, 1, &
  !                                 WriteToFile_Option = .TRUE., &
  !                                 FileNameBase_Option = 'MF_uCR' )

  CALL ComputeFromConserved_TwoMoment_MF(  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )

  !CALL ShowVariableFromMultifab(MF_uCR, 1, writetofile_option=.TRUE., FileNameBase_Option ='Conserved_Variables')

  CALL amrex_parmparse_build( PP, 'thornado' )
    CALL PP % query ( 'UseRealizabilityTimeStep', &
                       UseRealizabilityTimeStep )
    CALL PP % query ('dt_wrt', dt_wrt)
    CALL PP % query ('dt_chk', dt_chk)
    CALL PP % query ('iCycleChk', iCycleChk)
    CALL PP % query ('iCycleW', iCycleW)
  CALL amrex_parmparse_destroy( PP )

  wrt = .FALSE.

  chk = .FALSE.
  
  DO WHILE( ALL( t_new .LT. t_end ) )

    StepNo = StepNo + 1

    CALL ComputeTally_TwoMoment_MF &
           ( t_new, MF_uGF, MF_uCF, MF_uCR, &
             WriteTally_Option = .FALSE., Verbose_Option = .FALSE. )
    N_before = NeutrinoLeptonNumber_Interior

    CALL ReGrid

    CALL ComputeTally_TwoMoment_MF &
           ( t_new, MF_uGF, MF_uCF, MF_uCR, Verbose_Option = .FALSE. )

    IF( amrex_parallel_ioprocessor() .AND. N_before .NE. 0.0_DP )THEN
      IF( ABS( NeutrinoLeptonNumber_Interior - N_before ) &
            .GT. 1.0e-13_DP * ABS( N_before ) )THEN
        WRITE(*,'(A,ES13.6,A,ES13.6,A,I2)') '  REGRID dN/N = ', &
          ( NeutrinoLeptonNumber_Interior - N_before ) / N_before, &
          '  at t = ', t_new(0), '  nLevels = ', nLevels
      END IF
    END IF

    IF ( dt_rel .NE. 0.0_amrex_real ) THEN

      dt = dt_rel

    ELSE

      IF ( UseRealizabilityTimeStep ) THEN

        CALL ApplyBoundaryConditions_Euler_MF (MF_uCF)

        CALL ComputeTimeStep_TwoMoment_Realizability_MF (MF_uGF, MF_uCF, One, dt)

      ELSE

        CALL ComputeTimeStep_TwoMoment_MF( MF_uGF, CFL, dt )

      END IF

      dt = MINVAL( dt(0:nLevels-1) ) !#MINVAL( dt )
    END IF


    IF( ALL( t_new + dt .LE. t_end ) )THEN
      t_new = t_new + dt
    ELSE
      dt = t_end - t_new
      t_new  = t_end
    END IF
    IF( amrex_parallel_ioprocessor() )THEN
        WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
          'StepNo: ', StepNo(0), ' t = ', t_new(0) / UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel ), &
          ' dt = ', dt(0) /  UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel )
    END IF


    CALL Update_IMEX_RK_MF

    !CALL WritePlotFile

    !CALL WriteCheckpointFile

    CALL TimersStart_AMReX( Timer_AMReX_InputOutput )

    IF( iCycleW .GT. 0 )THEN
        IF( MOD( StepNo(0), iCycleW ) .EQ. 0 ) wrt = .TRUE.
      ELSE

        IF( ALL( t_new + dt .GT. t_wrt ) )THEN
          t_wrt = t_wrt + dt_wrt
          wrt   = .TRUE.

      END IF

    END IF

      !wrt  = .TRUE.
    IF( wrt )THEN
        CALL ComputeFromConserved_TwoMoment_MF(  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )

        CALL ComputeFromConserved_Euler_MF(  MF_uGF, MF_uCF, MF_uPF, MF_uAF )

        CALL ShowVariableFromMultiFab( MF_uPR, 1, &
                                  WriteToFile_Option = .TRUE., &
                                  FileNameBase_Option = 'MF_uPR' )
        !CALL ShowVariableFromMultiFab( MF_uPR, 2, &
        !                          WriteToFile_Option = .TRUE., &
        !                          FileNameBase_Option = 'MF_uI1' )
        !CALL ShowVariableFromMultiFab( MF_uPR, 3, &
        !                          WriteToFile_Option = .TRUE., &
        !                          FileNameBase_Option = 'MF_uI2' )

        CALL ShowVariableFromMultiFab( MF_uCR, 1, &
                                  WriteToFile_Option = .TRUE., &
                                  FileNameBase_Option = 'MF_uCR' )

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
             
        CALL ComputeTally_Euler_MF &
           ( t_new, MF_uGF, MF_uCF, Verbose_Option = .FALSE. )
        
        wrt = .FALSE.
        
        !CALL WriteFieldsAMReX_Checkpoint
    END IF

    CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

  CALL TimersStart_AMReX( Timer_AMReX_InputOutput )

  IF( iCycleChk .GT. 0 )THEN
    IF( MOD( StepNo(0), iCycleChk ) .EQ. 0 ) chk = .TRUE.
  ELSE
    IF( ALL( t_new + dt .GT. t_chk ) )THEN
      t_chk = t_chk + dt_chk
      chk   = .TRUE.
    END IF
  END IF

  IF( chk )THEN
    CALL ComputeFromConserved_TwoMoment_MF( MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )
    CALL ComputeFromConserved_Euler_MF( MF_uGF, MF_uCF, MF_uPF, MF_uAF )
    CALL WriteFieldsAMReX_Checkpoint
    chk = .FALSE.
  END IF

  CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

  END DO
  


  CALL FinalizeProgram



END PROGRAM main

SUBROUTINE ReGrid

  USE InputParsingModule, ONLY: &
    InitializeParameters, &
    nLevels, &
    nMaxLevels, &
    StepNo, &
    iRestart, &
    iReGrid, &
    t_old, &
    t_new, &
    UseTiling, &
    TagCriteria, &
    RefinementScheme, &
    DescribeProgramHeader_AMReX, &
    DEBUG, &
    UseAMR
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_init_module, ONLY: &
    amrex_init
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_regrid, &
    amrex_get_numlevels, &
    amrex_amrcore_init, &
    amrex_init_virtual_functions, &
    amrex_init_from_scratch
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
      ApplyPositivityLimiter_TwoMoment_MF
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF

    !TYPE(amrex_multifab) :: MF_uGS(0:nMaxLevels-1)
    !TYPE(amrex_multifab) :: MF_uMF(0:nMaxLevels-1)
    INTEGER :: iLevel, iErr, nLevelsOld
    

    IF( .NOT. UseAMR ) RETURN
    
    IF( DEBUG ) THEN
      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) THEN

        WRITE(*,*)
        WRITE(*,'(6x,A,I2.2)') 'nLevels (before regrid): ', nLevels
        WRITE(*,'(6x,A)') 'Regridding'

      END IF

    END IF

    nLevels = amrex_get_numlevels()
    !nLevelsOld = nLevels

    DO iLevel = 0, nMaxLevels-1

      CALL amrex_regrid( iLevel, t_new(iLevel) )

      nLevels = amrex_get_numlevels()

      IF( iLevel .GE. nLevels-1 ) EXIT

    END DO

    !DO iLevel = 0, nMaxLevels-1
      !PRINT *, 't_new(iLevel):',t_new(iLevel)
      !PRINT *, 'iLevel:', iLevel
      !CALL amrex_regrid( iLevel, t_new(iLevel) )
    !END DO


    !PRINT *, 'CALLED AMR regrid with levels:', nLevels
    
      CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )
      CALL AverageDown( MF_uGF, MF_uCF )
      CALL AverageDown( MF_uGF, MF_uCR)

      CALL ApplyPositivityLimiter_TwoMoment_MF ( MF_uGF, MF_uCF, MF_uCR )

      CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

      !PRINT *, 'CALLED AVG DWN from Regrid'
    
    t_old = t_old(0)
    t_new = t_new(0)
    
    IF( DEBUG ) THEN
      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )
      IF( amrex_parallel_ioprocessor() ) THEN
        WRITE(*,'(6x,A,I2.2)') 'nLevels (after regrid): ', nLevels
        WRITE(*,*)
      END IF
    END IF
    
END SUBROUTINE ReGrid
