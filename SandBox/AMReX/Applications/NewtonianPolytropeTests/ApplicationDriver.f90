PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator, &
    amrex_parallel_reduce_max

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Gram, &
    Centimeter
  USE FluidFieldsModule, ONLY: &
    iPF_D, &
    iCF_D

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE GeometryFieldsModule, ONLY: &
    iGF_Phi_N
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE InitializationModule, ONLY: &
    InitializeProgram
  USE FinalizationModule, ONLY: &
    FinalizeProgram
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeTimeStep_Euler_MF, &
    ComputeFromConserved_Euler_MF
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint
  USE MF_Euler_TallyModule, ONLY: &
    ComputeTally_Euler_MF, &
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
    ADMMass_OffGrid, &
    ADMMass_Interior
  USE MF_TimeSteppingModule_SSPRK_Newtonian, ONLY: &
    UpdateFluid_SSPRK_Newtonian_MF, &
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
    DEBUG, &
    UseTiling
  USE MF_Euler_TimersModule, ONLY: &
    TimeIt_AMReX_Euler
  USE MF_TimersModule, ONLY: &
    TimeIt_AMReX, &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_InputOutput, &
    FinalizeTimers_AMReX
  USE ReGridModule, ONLY: &
    ReGrid
  USE MF_UtilitiesModule

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: iErr, iLevel
  LOGICAL  :: wrt, chk
  REAL(DP) :: Timer_Evolution, CentralDensity
  REAL(DP), PARAMETER :: MaxCentralDensity = 1.0e15_DP * Gram / Centimeter**3

  TimeIt_AMReX       = .TRUE.
  TimeIt_AMReX_Euler = .TRUE.

  wrt = .FALSE.
  chk = .FALSE.

  CALL InitializeProgram
  CALL ShowVariableFromMultifab(MF_uCF, iCF_D, writetofile_option=.TRUE., FileNameBase_Option ='Density')
  CALL ShowVariableFromMultifab(MF_uGF, iGF_Phi_N, swXX_Option=[0,0,0], writetofile_option=.TRUE., FileNameBase_Option ='Potential')

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

  ! --- Begin evolution ---

  DO WHILE( MAXVAL( t_new ) .LT. t_end )

    StepNo = StepNo + 1

    t_old = t_new

    CALL ReGrid

    CALL ComputeTimeStep_Euler_MF( MF_uGF, MF_uCF, CFL, dt )

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
        WRITE(*,'(A)') 'CALL UpdateFluid_SSPRK_Newtonian_MF'

    END IF

    CALL UpdateFluid_SSPRK_Newtonian_MF

    CALL ShowVariableFromMultifab(MF_uCF, iCF_D, writetofile_option=.TRUE., FileNameBase_Option ='Density')
    CALL ShowVariableFromMultifab(MF_uGF, iGF_Phi_N, swXX_Option=[0,0,0], writetofile_option=.TRUE., FileNameBase_Option ='Potential')

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

    !CALL GetCentralDensity( CentralDensity )

    !IF( CentralDensity .GE. MaxCentralDensity ) EXIT


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

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_PlotFile &
             ( t_new(0), StepNo, MF_uGF, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF, &
               MF_uDF_Option = MF_uDF )
      !CALL ShowVariableFromMultiFab (MF_uGF, iGF_Phi_N, writetofile_option=.TRUE., FileNameBase_Option= 'Potential')
      CALL ComputeTally_Euler_MF &
             ( t_new, MF_uGF, MF_uCF, Verbose_Option = .TRUE. )

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

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t_new, &
               [ BaryonicMass_Initial   , BaryonicMass_OffGrid    ], &
               [ EulerMomentumX1_Initial, EulerMomentumX1_OffGrid ], &
               [ EulerMomentumX2_Initial, EulerMomentumX2_OffGrid ], &
               [ EulerMomentumX3_Initial, EulerMomentumX3_OffGrid ], &
               [ EulerEnergy_Initial    , EulerEnergy_OffGrid     ], &
               [ ElectronNumber_Initial , ElectronNumber_OffGrid  ], &
               [ ADMMass_Initial        , ADMMass_OffGrid, &
                 ADMMass_Interior ], &
               MF_uGF % BA % P, &
               iWriteFields_uGF = 1, &
               iWriteFields_uCF = 1, &
               iWriteFields_uCR = 0, &
               pMF_uGF_Option = MF_uGF % P, &
               pMF_uCF_Option = MF_uCF % P )

      CALL FinalizeTimers_AMReX &
             ( RestartProgramTimer_Option = .TRUE. )

      chk = .FALSE.

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_InputOutput )

  END SUBROUTINE WriteCheckpointFile


  SUBROUTINE GetCentralDensity( CentralDensity )

    REAL(DP), INTENT(out) :: CentralDensity

    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_box)       :: BX
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX2, iX3

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CentralDensity = -HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uPF % BA, MF_uPF % DM )

      CALL amrex_mfiter_build( MFI, MF_uPF(iLevel), tiling = UseTiling  )

      DO WHILE( MFI % next() )

        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % TileBox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          IF( IsNotLeafElement( uFM(iX_B0(1),iX_B0(2),iX_B0(3),1) ) ) CYCLE

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)

            CentralDensity &
              = MAX( CentralDensity, uPF(iX_B0(1),iX2,iX3,nDOFX*(iPF_D-1)+1) )

          END DO
          END DO

        END IF

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyFineMask( iMF_FineMask )

    END DO

    CALL amrex_parallel_reduce_max( CentralDensity )



  END SUBROUTINE GetCentralDensity


END PROGRAM main
