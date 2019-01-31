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

  USE MF_GeometryModule,        ONLY: &
    MF_ComputeGeometryX
  USE MF_InitializationModule,  ONLY: &
    MF_InitializeFields
  USE MF_Euler_UtilitiesModule, ONLY: &
    MF_ComputeFromConserved
  USE MF_SlopeLimiterModule_Euler, ONLY: &
    MF_ApplySlopeLimiter_Euler
  USE FinalizationModule, ONLY: &
    FinalizeProgram

  ! --- Checkpoint ---
  USE InputOutputModuleAMReX
  USE amrex_amr_module ! To call amrex_amrcore_init
  USE MyAmrDataModule
  USE MyAmrModule

  ! --- For slope limiter ---
  USE Euler_SlopeLimiterModule,       ONLY: &
    InitializeSlopeLimiter_Euler
  USE FluidFieldsModule,              ONLY: &
    Shock
  USE PolynomialBasisMappingModule,   ONLY: &
    InitializePolynomialBasisMapping
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre, ONLY: &
    InitializePolynomialBasis_Legendre

  IMPLICIT NONE

  INTEGER :: iLevel, iDim
  TYPE(amrex_box)                    :: BX
  TYPE(amrex_boxarray),  ALLOCATABLE :: BA(:)
  TYPE(amrex_distromap), ALLOCATABLE :: DM(:)
  TYPE(amrex_geometry),  ALLOCATABLE :: GEOM(:)

  ! --- Initialize AMReX ---

  CALL amrex_init()

  CALL amrex_amrcore_init() ! Gets refinement ratio

  ! --- Parse parameter file ---
  CALL MyAmrInit

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

  END IF

  CALL InitializeProgramHeader &
         ( ProgramName_Option = TRIM( ProgramName ), &
           nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, &
           xL_Option = xL, xR_Option = xR, bcX_Option = bcX )

  IF( amrex_parallel_ioprocessor() )THEN

    CALL DescribeProgramHeaderX

  END IF

  CALL InitializePolynomialBasisX_Lagrange
  CALL InitializePolynomialBasisX_Legendre

  CALL InitializePolynomialBasis_Lagrange
  CALL InitializePolynomialBasis_Legendre

  CALL InitializeReferenceElementX
  CALL InitializeReferenceElementX_Lagrange

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL )

  CALL InitializeSlopeLimiter_Euler &
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
             = LimiterThresholdParameter )

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

  DO iLevel = 0, nLevels
    CALL amrex_distromap_destroy( DM(iLevel) )
    CALL amrex_boxarray_destroy ( BA(iLevel) )
  END DO

  DO iDim = 1, 3
    CALL CreateMesh &
           ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
             amrex_problo(iDim), amrex_probhi(iDim) )
  END DO

  CALL InitializePolynomialBasisMapping &
    ( [0.0d0], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

  DO iLevel = 0, nLevels
    CALL MF_ComputeGeometryX( MF_uGF(iLevel) )
    CALL MF_InitializeFields &
           ( TRIM( ProgramName ), MF_uGF(iLevel), MF_uCF(iLevel) )
    CALL MF_ComputeFromConserved &
           ( MF_uGF(iLevel), MF_uCF(iLevel), MF_uPF(iLevel), MF_uAF(iLevel) )
  END DO

  CALL WriteFieldsAMReX_PlotFile &
         ( 0.0_DP, nLevels, GEOM, StepNo, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t_new, &
           MF_uGF % BA % P, &
           MF_uGF % P, &
           MF_uCF % P, &
           MF_uPF % P, &
           MF_uAF % P )

  ! --- START of evolution ---

  ALLOCATE( Shock(1:nX(1),1:nX(2),1:nX(3)) )
  CALL MF_ApplySlopeLimiter_Euler( nLevels, MF_uGF, MF_uCF )
  DEALLOCATE( Shock )
  StepNo(0) = StepNo(0) + 1

  CALL WriteFieldsAMReX_PlotFile &
         ( 0.1_DP, nLevels, GEOM, StepNo, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  ! --- END of evolution ---

  CALL MyAmrFinalize
  CALL ReadCheckpointFile

  ! --- Finalize everything ---

  CALL FinalizeProgram( nLevels, GEOM, MeshX )

  DEALLOCATE( GEOM )
  DEALLOCATE( BA )
  DEALLOCATE( DM )


END PROGRAM main

