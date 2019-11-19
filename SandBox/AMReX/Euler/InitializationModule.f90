MODULE InitializationModule

  ! --- AMReX Modules ---
  USE amrex_fort_module, ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_amr_module, ONLY: &
    amrex_init, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray,         &
    amrex_boxarray_build,   &
    amrex_boxarray_destroy, &
    amrex_print
  USE amrex_distromap_module, ONLY: &
    amrex_distromap,       &
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
    nDOFX,                  &
    nNodesX
  USE PolynomialBasisModuleX_Lagrange,  ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre,  ONLY: &
    InitializePolynomialBasisX_Legendre
  USE ReferenceElementModuleX,          ONLY: &
    InitializeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE MeshModule,                       ONLY: &
    MeshX,      &
    CreateMesh, &
    DestroyMesh
  USE EquationOfStateModule,            ONLY: &
    InitializeEquationOfState
  USE EquationOfStateModule_TABLE,      ONLY: &
    MinD, &
    MaxD, &
    MinT, &
    MaxT, &
    MinY, &
    MaxY
  USE GeometryFieldsModule,             ONLY: &
    nGF, &
    CoordinateSystem
  USE FluidFieldsModule,                ONLY: &
    nCF, &
    nPF, &
    nAF, &
    CreateFluidFields
  USE Euler_SlopeLimiterModule,         ONLY: &
    InitializeSlopeLimiter_Euler
  USE PolynomialBasisMappingModule,     ONLY: &
    InitializePolynomialBasisMapping
  USE PolynomialBasisModule_Lagrange,   ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre,   ONLY: &
    InitializePolynomialBasis_Legendre
  USE Euler_PositivityLimiterModule,    ONLY: &
    InitializePositivityLimiter_Euler
  USE InputOutputModuleAMReX,           ONLY: &
    ReadCheckpointFile,          &
    WriteFieldsAMReX_Checkpoint, &
    WriteFieldsAMReX_PlotFile
  USE UnitsModule,                      ONLY: &
    SolarMass, &
    UnitsDisplay

  ! --- Local modules ---
  USE MF_Euler_UtilitiesModule,         ONLY: &
    MF_ComputeFromConserved
  USE MF_GeometryModule,                ONLY: &
    MF_ComputeGeometryX, &
    MF_ComputeGravitationalPotential
  USE MF_Euler_SlopeLimiterModule,      ONLY: &
    MF_ApplySlopeLimiter_Euler
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    MF_ApplyPositivityLimiter_Euler
  USE MF_InitializationModule,          ONLY: &
    MF_InitializeFields
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_InitializeFluid_SSPRK
  USE MyAmrDataModule,                  ONLY: &
    MF_uGF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF
  USE MyAmrModule,                      ONLY: &
    t_end,                     &
    t,                         &
    dt,                        &
    t_wrt,                     &
    dt_wrt,                    &
    t_chk,                     &
    dt_chk,                    &
    CFL,                       &
    nNodes,                    &
    nStages,                   &
    nX,                        &
    swX,                       &
    bcX,                       &
    xL,                        &
    xR,                        &
    ProgramName,               &
    CoordSys,                  &
    UseSlopeLimiter,           &
    UseCharacteristicLimiting, &
    UseTroubledCellIndicator,  &
    SlopeTolerance,            &
    BetaTVD,                   &
    BetaTVB,                   &
    LimiterThresholdParameter, &
    UseConservativeCorrection, &
    UsePositivityLimiter,      &
    Min_1,                     &
    Min_2,                     &
    Gamma_IDEAL,               &
    EquationOfState,           &
    EosTableName,              &
    StepNo,                    &
    nLevels,                   &
    iRestart,                  &
    MaxGridSizeX,              &
    BA,                        &
    DM,                        &
    GEOM,                      &
    MyAmrInit
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler,      &
    TimersStop_AMReX_Euler,       &
    Timer_AMReX_Euler_Initialize, &
    Timer_AMReX_Euler_InputOutput

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

  LOGICAL, PUBLIC :: wrt, chk

  REAL(AR), PARAMETER :: Zero = 0.0_AR
  REAL(AR), PARAMETER :: One  = 1.0_AR
  REAL(AR), PARAMETER :: Two  = 2.0_AR

CONTAINS


  SUBROUTINE InitializeProgram

    INTEGER               :: iLevel, iDim
    TYPE(amrex_parmparse) :: PP
    TYPE(amrex_box)       :: BX
    REAL(AR)              :: Mass

    ! --- Initialize AMReX ---
    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    ! --- Parse parameter file ---
    CALL MyAmrInit

    IF( iRestart .LT. 0 )THEN

      BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

      ALLOCATE( BA(0:nLevels-1) )
      DO iLevel = 0, nLevels-1
        CALL amrex_boxarray_build( BA(iLevel), BX )
      END DO

      DO iLevel = 0, nLevels-1
        CALL BA(iLevel) % maxSize( MaxGridSizeX )
      END DO

      ALLOCATE( GEOM(0:nLevels-1) )
      ALLOCATE( DM  (0:nLevels-1) )

      DO iLevel = 0, nLevels-1
        CALL amrex_geometry_build ( GEOM(iLevel), BX )
        CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )
      END DO

      DO iLevel = 0, nLevels-1
        CALL amrex_multifab_build &
               ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX(1) )
        CALL MF_uGF(iLevel) % SetVal( Zero )
        CALL amrex_multifab_build &
               ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX(1) )
        CALL MF_uCF(iLevel) % SetVal( Zero )
        CALL amrex_multifab_build &
               ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, swX(1) )
        CALL MF_uPF(iLevel) % SetVal( Zero )
        CALL amrex_multifab_build &
               ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX(1) )
        CALL MF_uAF(iLevel) % SetVal( Zero )
      END DO

      t     = Zero
      dt    = Zero
      t_wrt = dt_wrt
      t_chk = dt_chk

    ELSE

      CALL ReadCheckpointFile( iRestart )
      t_chk = t(0) + dt_chk

    END IF

    wrt = .FALSE.
    chk = .FALSE.

    ! -- End of initializing AMReX ---

    ! --- Initialize thornado ---

    CoordinateSystem = TRIM( CoordSys )

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A6,A)')             'Name: ', TRIM( ProgramName )
      WRITE(*,*)
      WRITE(*,'(4x,A24,ES10.3E2,A,A)') 't_end   =', &
        t_end  / UnitsDisplay % TimeUnit, ' ', TRIM( UnitsDisplay  % TimeLabel )
      WRITE(*,'(4x,A24,ES10.3E2,A,A)') 'dt_wrt  =', &
        dt_wrt / UnitsDisplay % TimeUnit, ' ', TRIM( UnitsDisplay  % TimeLabel )
      WRITE(*,'(4x,A24,ES10.3E2,A,A)') 'dt_chk  =', &
        dt_chk / UnitsDisplay % TimeUnit, ' ', TRIM( UnitsDisplay  % TimeLabel )
      WRITE(*,'(4x,A24,I3.2)')         'nNodes  =', nNodes
      WRITE(*,'(4x,A24,I3.2)')         'nStages =', nStages
      WRITE(*,'(4x,A24,I3.2)')         'nDimsX  =', amrex_spacedim
      WRITE(*,'(4x,A24,ES10.3E2)')     'Gamma   =', Gamma_IDEAL
      WRITE(*,'(5x,A24,A)')            'CoordinateSystem = ', CoordinateSystem
      WRITE(*,'(4x,A24,3I7.6)')        'nX           =', nX
      WRITE(*,'(4x,A24,3I7.6)')        'swX          =', swX
      WRITE(*,'(4x,A24,3I7.6)')        'bcX          =', bcX
      WRITE(*,'(4x,A24,3I7.6)')        'MaxGridSizeX =', MaxGridSizeX

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
           ( [Zero], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    Mass = Zero
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'Mass', Mass )
    CALL amrex_parmparse_destroy( PP )

    IF( ProgramName .EQ. 'StandingAccretionShock_Relativistic' )THEN
      Mass = Mass * SolarMass
      CALL MF_ComputeGeometryX( MF_uGF, Mass )
    ELSE
      CALL MF_ComputeGeometryX( MF_uGF, 0.0_AR )
    END IF

    IF( ProgramName .EQ. 'StandingAccretionShock' ) &
      CALL MF_ComputeGravitationalPotential( MF_uGF, Mass )

    IF( EquationOfState .EQ. 'TABLE' )THEN

      CALL InitializeEquationOfState &
             ( EquationOfState_Option = EquationOfState, &
               EquationOfStateTableName_Option &
                 = EosTableName )

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = ( One + EPSILON(One) ) * MinD, &
               Min_2_Option = ( One + EPSILON(One) ) * MinT, &
               Min_3_Option = ( One + EPSILON(One) ) * MinY, &
               Max_1_Option = ( One - EPSILON(One) ) * MaxD, &
               Max_2_Option = ( One - EPSILON(One) ) * MaxT, &
               Max_3_Option = ( One - EPSILON(One) ) * MaxY )

    ELSE

      CALL InitializeEquationOfState &
             ( EquationOfState_Option = EquationOfState, &
               Gamma_IDEAL_Option = Gamma_IDEAL )

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = Min_1, &
               Min_2_Option = Min_2 )

    END IF

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
               = LimiterThresholdParameter, &
             UseConservativeCorrection_Option &
               = UseConservativeCorrection, &
             Verbose_Option &
               = amrex_parallel_ioprocessor() )

    CALL MF_InitializeFluid_SSPRK &
           ( nStages, BA, DM, &
             Verbose_Option = amrex_parallel_ioprocessor() )
    IF( amrex_parallel_ioprocessor() ) WRITE(*,'(A6,A,ES11.3E3)') &
      '', 'CFL: ', &
      CFL * ( amrex_spacedim * ( Two * nNodes - One ) )

    ! --- Allocates 'Shock' and sets units for fluid fields ---
    CALL CreateFluidFields( nX, swX, amrex_parallel_ioprocessor() )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    IF( iRestart .LT. 0 )THEN

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Initialize )

      CALL MF_InitializeFields( TRIM( ProgramName ), MF_uGF, MF_uCF )

      CALL MF_ApplySlopeLimiter_Euler     ( MF_uGF, MF_uCF, GEOM )
      CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_uCF )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )
      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    END IF

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    DO iLevel = 0, nLevels-1
      CALL amrex_distromap_destroy( DM(iLevel) )
      CALL amrex_boxarray_destroy ( BA(iLevel) )
    END DO

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(A)') '  Evolving fields...'
      WRITE(*,*)
    END IF

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

  END SUBROUTINE InitializeProgram


END MODULE InitializationModule
