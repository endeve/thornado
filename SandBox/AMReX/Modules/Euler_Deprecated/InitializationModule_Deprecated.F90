MODULE InitializationModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,      ONLY: &
    amrex_spacedim
  USE amrex_init_module,      ONLY: &
    amrex_init
  USE amrex_amrcore_module,   ONLY: &
    amrex_amrcore_init
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray,         &
    amrex_boxarray_build,   &
    amrex_boxarray_destroy
  USE amrex_distromap_module, ONLY: &
    amrex_distromap,       &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry, &
    amrex_geometry_build
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, &
    amrex_multifab_build
  USE amrex_parallel_module,  ONLY: &
    amrex_parallel_ioprocessor
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
    CreateMesh
  USE EquationOfStateModule,            ONLY: &
    InitializeEquationOfState
  USE EquationOfStateModule_TABLE,      ONLY: &
    Min_D, &
    Max_D, &
    Min_T, &
    Max_T, &
    Min_Y, &
    Max_Y
  USE GeometryFieldsModule,             ONLY: &
    nGF,                     &
    DescribeGeometryFields,  &
    SetUnitsGeometryFields,  &
    CoordinateSystem
  USE FluidFieldsModule,                ONLY: &
    nCF,                            &
    nPF,                            &
    nAF,                            &
    nDF,                            &
    DescribeFluidFields_Primitive,  &
    DescribeFluidFields_Conserved,  &
    DescribeFluidFields_Auxiliary,  &
    DescribeFluidFields_Diagnostic, &
    SetUnitsFluidFields
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
    ReadCheckpointFile, &
    WriteFieldsAMReX_PlotFile
  USE UnitsModule,                      ONLY: &
    SolarMass, &
    UnitsDisplay

  ! --- Local modules ---

  USE MF_KindModule,                    ONLY: &
    DP, &
    Zero, &
    One, &
    Two
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
  USE MF_FieldsModule,                  ONLY: &
    MF_uGF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE InputParsingModule,               ONLY: &
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
    UsePhysicalUnits,          &
    UseSlopeLimiter,           &
    SlopeLimiterMethod,        &
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
    InitializeParameters
  USE MF_Euler_TallyModule,     ONLY: &
    MF_InitializeTally_Euler, &
    MF_ComputeTally_Euler
  USE TimersModule_AMReX_Euler, ONLY: &
    InitializeTimers_AMReX_Euler,  &
    TimersStart_AMReX_Euler,      &
    TimersStop_AMReX_Euler,       &
    Timer_AMReX_Euler_Initialize, &
    Timer_AMReX_Euler_InputOutput

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

  LOGICAL, PUBLIC :: wrt, chk


CONTAINS


  SUBROUTINE InitializeProgram

    INTEGER               :: iLevel, iDim
    TYPE(amrex_parmparse) :: PP
    TYPE(amrex_box)       :: BX

    LOGICAL               :: SolveGravity
    REAL(DP)              :: Mass

    ! --- Initialize AMReX ---

    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL InitializeTimers_AMReX_Euler

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    ! --- Parse parameter file ---

    CALL InitializeParameters

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
               ( MF_uGF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nGF, swX )
        CALL MF_uGF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )
        CALL MF_uCF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPF, swX )
        CALL MF_uPF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uAF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAF, swX )
        CALL MF_uAF(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uDF(iLevel), BA(iLevel), DM(iLevel), nDOFX * nDF, swX )
        CALL MF_uDF(iLevel) % SetVal( Zero )

      END DO

      t     = Zero
      dt    = Zero
      t_wrt = dt_wrt
      t_chk = dt_chk

    ELSE

      CALL ReadCheckpointFile( iRestart )
      t_chk = t(0) + dt_chk
      t_wrt = t(0) + dt_wrt

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

    SolveGravity = .FALSE.
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'SolveGravity', SolveGravity )
    CALL amrex_parmparse_destroy( PP )

    Mass = Zero
    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % query( 'Mass', Mass )
    CALL amrex_parmparse_destroy( PP )

    IF( UsePhysicalUnits ) &
      Mass = Mass * SolarMass

#if defined HYDRO_RELATIVISTIC

      CALL MF_ComputeGeometryX( MF_uGF, Mass )

#else

      CALL MF_ComputeGeometryX( MF_uGF, Zero )

      IF( SolveGravity ) THEN

        CALL MF_ComputeGravitationalPotential( MF_uGF, Mass )

      END IF

#endif

    CALL SetUnitsGeometryFields

    CALL DescribeFluidFields_Conserved( amrex_parallel_ioprocessor() )

    CALL DescribeFluidFields_Primitive( amrex_parallel_ioprocessor() )

    CALL DescribeFluidFields_Auxiliary( amrex_parallel_ioprocessor() )

    CALL DescribeFluidFields_Diagnostic( amrex_parallel_ioprocessor() )

    CALL SetUnitsFluidFields( TRIM( CoordinateSystem ), &
                              Verbose_Option = amrex_parallel_ioprocessor() )

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
             SlopeLimiterMethod_Option &
               = SlopeLimiterMethod, &
             LimiterThresholdParameter_Option &
               = LimiterThresholdParameter, &
             UseConservativeCorrection_Option &
               = UseConservativeCorrection, &
             Verbose_Option &
               = amrex_parallel_ioprocessor() )

    IF( EquationOfState .EQ. 'TABLE' )THEN

      CALL InitializeEquationOfState &
             ( EquationOfState_Option = EquationOfState, &
               EquationOfStateTableName_Option &
                 = EosTableName )

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = ( One + EPSILON(One) ) * Min_D, &
               Min_2_Option = ( One + EPSILON(One) ) * Min_T, &
               Min_3_Option = ( One + EPSILON(One) ) * Min_Y, &
               Max_1_Option = ( One - EPSILON(One) ) * Max_D, &
               Max_2_Option = ( One - EPSILON(One) ) * Max_T, &
               Max_3_Option = ( One - EPSILON(One) ) * Max_Y )

    ELSE

      CALL InitializeEquationOfState &
             ( EquationOfState_Option = EquationOfState, &
               Gamma_IDEAL_Option = Gamma_IDEAL, &
               Verbose_Option = amrex_parallel_ioprocessor() )

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = Min_1, &
               Min_2_Option = Min_2 )

    END IF

    CALL MF_InitializeFluid_SSPRK &
           ( nStages, BA, DM, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    IF( amrex_parallel_ioprocessor() ) WRITE(*,'(A6,A,ES11.3E3)') &
      '', 'CFL: ', &
      CFL * ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

    CALL MF_InitializeTally_Euler

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    IF( iRestart .LT. 0 )THEN

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Initialize )

      CALL MF_InitializeFields( TRIM( ProgramName ), MF_uGF, MF_uCF, GEOM )

      CALL MF_ApplySlopeLimiter_Euler( MF_uGF, MF_uCF, MF_uDF, GEOM )

      CALL MF_ApplyPositivityLimiter_Euler( MF_uGF, MF_uCF, MF_uDF )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF, &
               MF_uDF_Option = MF_uDF )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    END IF

    CALL MF_ComputeTally_Euler &
           ( GEOM, MF_uGF, MF_uCF, t(0), &
             SetInitialValues_Option = .TRUE., Verbose_Option = .FALSE. )

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
