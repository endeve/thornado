MODULE ProblemInitializationModule

  ! --- AMReX Modules ---
  USE amrex_fort_module, ONLY: &
    amrex_real, &
    amrex_spacedim
  USE amrex_amr_module, ONLY: &
    amrex_init, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray, &
    amrex_boxarray_build, &
    amrex_boxarray_destroy, &
    amrex_print
  USE amrex_distromap_module, ONLY: &
    amrex_distromap, &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse,       &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,              ONLY: &
    DescribeProgramHeaderX, &
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
    nCF, nPF, nAF, &
    CreateFluidFields
  USE Euler_SlopeLimiterModule,         ONLY: &
    Euler_InitializeSlopeLimiter
  USE PolynomialBasisMappingModule,     ONLY: &
    InitializePolynomialBasisMapping
  USE PolynomialBasisModule_Lagrange,   ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre,   ONLY: &
    InitializePolynomialBasis_Legendre
  USE Euler_PositivityLimiterModule,    ONLY: &
    Euler_InitializePositivityLimiter
  USE InputOutputModuleAMReX

  ! --- Local modules ---
  USE MF_Euler_UtilitiesModule,         ONLY: &
    MF_ComputeFromConserved
  USE MF_GeometryModule,                ONLY: &
    MF_ComputeGeometryX, &
    MF_ComputeGravitationalPotential
  USE MF_Euler_SlopeLimiterModule,      ONLY: &
    MF_Euler_ApplySlopeLimiter
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    MF_Euler_ApplyPositivityLimiter
  USE MF_InitializationModule,          ONLY: &
    MF_InitializeFields
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_InitializeFluid_SSPRK
  USE MyAmrDataModule
  USE MyAmrModule

  IMPLICIT NONE

  PUBLIC :: InitializeProblem

  INTEGER,                            PUBLIC :: iLevel, iDim
  TYPE(amrex_parmparse),              PUBLIC :: PP
  TYPE(amrex_box),                    PUBLIC :: BX
  REAL(amrex_real),                   PUBLIC :: Mass
  REAL(amrex_real),                   PUBLIC :: t_wrt, t_chk
  LOGICAL,                            PUBLIC :: wrt, chk


CONTAINS


  SUBROUTINE InitializeProblem( CheckpointRestart_Option )

    INTEGER, INTENT(in), OPTIONAL :: CheckpointRestart_Option

    ! --- Initialize AMReX ---
    CALL amrex_init()

    CALL amrex_amrcore_init()

    ! --- Parse parameter file ---
    CALL MyAmrInit

    IF( .NOT. PRESENT( CheckpointRestart_Option ) )THEN

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
        CALL MF_uGF(iLevel) % SetVal( 0.0_amrex_real )
        CALL amrex_multifab_build &
               ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX(1) )
        CALL MF_uCF(iLevel) % SetVal( 0.0_amrex_real )
        CALL amrex_multifab_build &
               ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, swX(1) )
        CALL MF_uPF(iLevel) % SetVal( 0.0_amrex_real )
        CALL amrex_multifab_build &
               ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX(1) )
        CALL MF_uAF(iLevel) % SetVal( 0.0_amrex_real )
      END DO

    ELSE

      CALL ReadCheckpointFile( CheckpointRestart_Option )

    END IF

    ! -- End of initializing AMReX ---

    ! --- Initialize thornado ---

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
      WRITE(*,'(A5,A24,A)')        '', 'CoordinateSystem = ', CoordinateSystem
      WRITE(*,'(A4,A24,3I7.6)')    '', 'nX          =', nX
      WRITE(*,'(A4,A24,3I7.6)')    '', 'swX         =', swX
      WRITE(*,'(A4,A24,3I7.6)')    '', 'bcX         =', bcX
      WRITE(*,'(A4,A24,3I7.6)')    '', 'MaxGridSize =', MaxGridSize

      CALL DescribeProgramHeaderX

    END IF

    DO iDim = 1, 3
      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )
    END DO

    CALL InitializePolynomialBasisX_Lagrange
    CALL InitializePolynomialBasisX_Legendre

    CALL InitializePolynomialBasis_Lagrange
    CALL InitializePolynomialBasis_Legendre

    CALL InitializePolynomialBasisMapping &
           ( [0.0d0], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL MF_ComputeGeometryX( MF_uGF )

    Mass = 0.0_amrex_real
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'Mass', Mass )
    CALL amrex_parmparse_destroy( PP )

    CALL MF_ComputeGravitationalPotential( MF_uGF, Mass )

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
           ( Min_1_Option = Min_1, &
             Min_2_Option = Min_2, &
             UsePositivityLimiter_Option = UsePositivityLimiter, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL MF_InitializeFluid_SSPRK &
           ( nStages, BA, DM, &
             Verbose_Option = amrex_parallel_ioprocessor() )
    IF( amrex_parallel_ioprocessor() ) WRITE(*,'(A6,A,ES11.3E3)') &
      '', 'CFL: ', &
      CFL * ( amrex_spacedim * ( 2.0_amrex_real * nNodes - 1.0_amrex_real ) )

    CALL MF_InitializeFields( TRIM( ProgramName ), MF_uGF, MF_uCF )
    CALL CreateFluidFields( nX, swX, amrex_parallel_ioprocessor() )

    CALL MF_Euler_ApplySlopeLimiter     ( MF_uGF, MF_uCF, GEOM )
    CALL MF_Euler_ApplyPositivityLimiter( MF_uGF, MF_uCF )

    CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    DO iLevel = 0, nLevels
      CALL amrex_distromap_destroy( DM(iLevel) )
      CALL amrex_boxarray_destroy ( BA(iLevel) )
    END DO

    ! --- Beginning of evolution ---

    t     = 0.0_amrex_real
    t_wrt = dt_wrt
    t_chk = dt_chk
    wrt   = .FALSE.
    chk   = .FALSE.

    IF( .NOT. PRESENT( CheckpointRestart_Option ) )THEN
      CALL WriteFieldsAMReX_PlotFile &
             ( 0.0e0_amrex_real, StepNo, &
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
    END IF

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A,ES13.6E3)') 't = ', t
    END IF

  END SUBROUTINE InitializeProblem


END MODULE ProblemInitializationModule
