PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_base_module
  USE amrex_fort_module

  ! --- thornado Modules ---

  USE KindModule,                       ONLY: &
    DP
  USE ProgramHeaderModule,              ONLY: &
    InitializeProgramHeader,                  &
    DescribeProgramHeaderX,                   &
    nDOFX, nNodesX
  USE PolynomialBasisModuleX_Lagrange,  ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre,  ONLY: &
    InitializePolynomialBasisX_Legendre
  USE ReferenceElementModuleX,          ONLY: &
    InitializeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE MeshModule,                       ONLY: &
    MeshX, CreateMesh, DestroyMesh
  USE EquationOfStateModule,            ONLY: &
    InitializeEquationOfState
  USE GeometryFieldsModule,             ONLY: &
    nGF, CoordinateSystem
  USE FluidFieldsModule,                ONLY: &
    nCF, nPF, nAF
  USE InputOutputModuleAMReX,           ONLY: &
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile

  ! --- Local Modules ---

  USE MF_GeometryModule,                ONLY: &
    MF_ComputeGeometryX
  USE MF_InitializationModule,          ONLY: &
    MF_InitializeFields
  USE MF_Euler_UtilitiesModule,         ONLY: &
    MF_ComputeFromConserved, &
    MF_ComputeTimeStep
  USE MF_Euler_SlopeLimiterModule,      ONLY: &
    MF_Euler_ApplySlopeLimiter
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    MF_Euler_ApplyPositivityLimiter
  USE MF_Euler_dgDiscretizationModule,  ONLY: &
    MF_Euler_ComputeIncrement
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_InitializeFluid_SSPRK, &
    MF_UpdateFluid_SSPRK
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE MF_UtilitiesModule,               ONLY: &
    MakeMF_Diff

  ! --- Checkpoint ---
  USE InputOutputModuleAMReX
  USE amrex_amr_module, ONLY: &
    amrex_amrcore_init
  USE MyAmrDataModule
  USE MyAmrModule

  ! --- For slope limiter ---
  USE Euler_SlopeLimiterModule,       ONLY: &
    Euler_InitializeSlopeLimiter
  USE FluidFieldsModule,              ONLY: &
    Shock
  USE PolynomialBasisMappingModule,   ONLY: &
    InitializePolynomialBasisMapping
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre, ONLY: &
    InitializePolynomialBasis_Legendre

  ! --- For positivity limiter ---
  USE Euler_PositivityLimiterModule, ONLY: &
    Euler_InitializePositivityLimiter

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER                            :: iLevel, iDim
  TYPE(amrex_box)                    :: BX
  TYPE(amrex_boxarray),  ALLOCATABLE :: BA(:)
  TYPE(amrex_distromap), ALLOCATABLE :: DM(:)
  TYPE(amrex_geometry),  ALLOCATABLE :: GEOM(:)

  REAL(amrex_real) :: Timer_Evolution

  CALL MakeMF_Diff( 2, 4 )

  ! --- Initialize AMReX ---
  CALL amrex_init()

  CALL amrex_amrcore_init()

  ! --- Parse parameter file ---
  CALL MyAmrInit

  BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

  ALLOCATE( BA(0:nLevels) )
  DO iLevel = 0, nLevels
    CALL amrex_boxarray_build( BA(iLevel), BX )
  END DO

  DO iLevel = 0, nLevels
    CALL BA(iLevel) % maxSize( MaxGridSize )
  END DO

  ALLOCATE( GEOM(0:nLevels) )
  ALLOCATE( DM  (0:nLevels) )

  DO iLevel = 0, nLevels
    CALL amrex_geometry_build ( GEOM(iLevel), BX )
    CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )
  END DO

  ! -- (Almost) end of initializing AMReX ---

  ! --- Initialize thornado ---

  CALL InitializeProgramHeader &
         ( ProgramName_Option = TRIM( ProgramName ), &
           nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, &
           xL_Option = xL, xR_Option = xR, bcX_Option = bcX, &
           Verbose_Option = amrex_parallel_ioprocessor() )

  DO iLevel = 0, nLevels
    CALL amrex_multifab_build &
           ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX(1) )
    CALL amrex_multifab_build &
           ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX(1) )
    CALL amrex_multifab_build &
           ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, swX(1) )
    CALL amrex_multifab_build &
           ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX(1) )
  END DO

  CoordinateSystem = TRIM( CoordSys )

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(A4,A6,A)')         '', 'Name: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A24,ES10.3E2)') '', 't_end   =', t_end
    WRITE(*,'(A4,A24,ES10.3E2)') '', 'dt_wrt  =', dt_wrt
    WRITE(*,'(A4,A24,I7.6)')     '', 'nNodes  =', nNodes
    WRITE(*,'(A4,A24,I7.6)')     '', 'nStages =', nStages
    WRITE(*,'(A4,A24,I3.2)')     '', 'nDimsX  =', amrex_spacedim
    WRITE(*,'(A4,A24,ES10.3E2)') '', 'Gamma   =', Gamma_IDEAL
    WRITE(*,'(A5,A24,A)')        '', 'CoordinateSystem = ', TRIM( CoordSys )
    WRITE(*,'(A4,A24,3I7.6)')    '', 'nX          =', nX
    WRITE(*,'(A4,A24,3I7.6)')    '', 'swX         =', swX
    WRITE(*,'(A4,A24,3I7.6)')    '', 'bcX         =', bcX
    WRITE(*,'(A4,A24,3I7.6)')    '', 'MaxGridSize =', MaxGridSize

    CALL DescribeProgramHeaderX

  END IF

  DO iDim = 1, 3
    CALL CreateMesh &
           ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
             amrex_problo(iDim), amrex_probhi(iDim) )
  END DO

  CALL InitializePolynomialBasisX_Lagrange
  CALL InitializePolynomialBasisX_Legendre

  CALL InitializePolynomialBasis_Lagrange
  CALL InitializePolynomialBasis_Legendre

  CALL InitializePolynomialBasisMapping &
         ( [0.0d0], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

  CALL InitializeReferenceElementX
  CALL InitializeReferenceElementX_Lagrange

  DO iLevel = 0, nLevels
    CALL MF_ComputeGeometryX( MF_uGF(iLevel) )
  END DO

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL )

  CALL Euler_InitializeSlopeLimiter &
         ( BetaTVD_Option &
             = BetaTVD, &
           BetaTVB_Option &
             = BetaTVB, &
           SlopeTolerance_Option &
             = SlopeTolerance, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter, &
           Verbose_Option &
             = amrex_parallel_ioprocessor() )

  CALL Euler_InitializePositivityLimiter &
         ( Min_1_Option = 1.0d-12, &
           Min_2_Option = 1.0d-12, &
           UsePositivityLimiter_Option = .TRUE., &
           Verbose_Option = amrex_parallel_ioprocessor() )

  DO iLevel = 0, nLevels
    CALL MF_InitializeFields &
           ( TRIM( ProgramName ), MF_uGF(iLevel), MF_uCF(iLevel) )
  END DO

  ALLOCATE( Shock(1:nX(1),1:nX(2),1:nX(3)) )
  CALL MF_Euler_ApplySlopeLimiter &
         ( nLevels, MF_uGF, MF_uCF, GEOM )
  CALL MF_Euler_ApplyPositivityLimiter &
         ( nLevels, MF_uGF, MF_uCF )

  DO iLevel = 0, nLevels
    CALL MF_ComputeFromConserved &
           ( MF_uGF(iLevel), MF_uCF(iLevel), &
             MF_uPF(iLevel), MF_uAF(iLevel) )
  END DO

  CALL MF_InitializeFluid_SSPRK &
         ( nLevels, nStages, BA, DM, &
           Verbose_Option = amrex_parallel_ioprocessor() )

  DO iLevel = 0, nLevels
    CALL amrex_distromap_destroy( DM(iLevel) )
    CALL amrex_boxarray_destroy ( BA(iLevel) )
  END DO

  IF( amrex_parallel_ioprocessor() )THEN
    WRITE(*,*)
    WRITE(*,'(A,ES13.6E3)') 't = ', 0.0_amrex_real
  END IF

  CALL WriteFieldsAMReX_PlotFile &
         ( 0.0e0_amrex_real, nLevels, GEOM, StepNo, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  ! --- Beginning of evolution ---

  t  = 0.0_amrex_real

  IF( amrex_parallel_ioprocessor() ) &
    Timer_Evolution = MPI_WTIME()

  DO WHILE( ALL( t .LT. t_end ) )

    StepNo = StepNo + 1

    CALL MF_ComputeTimeStep( nLevels, MF_uGF, MF_uCF, CFL, dt )
    t = t + dt

    IF( amrex_parallel_ioprocessor() )THEN
      IF( MOD( StepNo(0), iCycleD ) .EQ. 0 ) &
        WRITE(*,'(A5,A,I6.6,A,ES13.6E3,A,ES13.6E3)') &
          '', 'StepNo: ', StepNo(0), ', t = ', t, ', dt = ', dt(0)
    END IF

    CALL MF_UpdateFluid_SSPRK &
           ( nLevels, t, dt, MF_uGF, MF_uCF, &
             GEOM, MF_Euler_ComputeIncrement )

    IF( MOD( StepNo(0), iCycleW ) .EQ. 0 )THEN

      DO iLevel = 0, nLevels
        CALL MF_ComputeFromConserved &
               ( MF_uGF(iLevel), MF_uCF(iLevel), &
                 MF_uPF(iLevel), MF_uAF(iLevel) )
      END DO

      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), nLevels, GEOM, StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF )
    END IF

    IF( MOD( StepNo(0), iCycleChk ) .EQ. 0 )THEN

      DO iLevel = 0, nLevels
        CALL MF_ComputeFromConserved &
               ( MF_uGF(iLevel), MF_uCF(iLevel), &
                 MF_uPF(iLevel), MF_uAF(iLevel) )
      END DO

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t, &
               MF_uGF % BA % P, &
               MF_uGF % P, &
               MF_uCF % P, &
               MF_uPF % P, &
               MF_uAF % P )
     
    END IF

  END DO

  ! --- END of evolution ---

  IF( amrex_parallel_ioprocessor() )THEN
    WRITE(*,*)
    WRITE(*,'(A,ES13.6E3,A)') &
      'Total evolution time: ', MPI_WTIME() - Timer_Evolution, ' s'
  END IF

  DO iLevel = 0, nLevels
    CALL MF_ComputeFromConserved &
           ( MF_uGF(iLevel), MF_uCF(iLevel), &
             MF_uPF(iLevel), MF_uAF(iLevel) )
  END DO

  StepNo = StepNo + 1
  CALL WriteFieldsAMReX_PlotFile &
         ( t(0), nLevels, GEOM, StepNo, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t, &
           MF_uGF % BA % P, &
           MF_uGF % P, &
           MF_uCF % P, &
           MF_uPF % P, &
           MF_uAF % P )

  ! --- Finalize everything ---

  CALL FinalizeProgram( nLevels, GEOM, MeshX )

  DEALLOCATE( Shock )
  DEALLOCATE( GEOM )
  DEALLOCATE( BA )
  DEALLOCATE( DM )


END PROGRAM main

