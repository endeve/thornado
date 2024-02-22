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
  USE GeometryFieldsModule,             ONLY: &
    nGF,                     &
    DescribeGeometryFields,  &
    SetUnitsGeometryFields,  &
    CoordinateSystem
  USE MagnetofluidfieldsModule,                ONLY: &
    nCM,                            &
    nPM,                            &
    nAM,                            &
    nDM,                            &
    DescribeFields_Primitive,  &
    DescribeFields_Conserved,  &
    DescribeFields_Auxiliary,  &
    DescribeFields_Diagnostic, &
    SetUnitsFields
  USE MHD_SlopeLimiterModule,           ONLY: &
    InitializeSlopeLimiter_MHD
  USE PolynomialBasisMappingModule,     ONLY: &
    InitializePolynomialBasisMapping
  USE PolynomialBasisModule_Lagrange,   ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre,   ONLY: &
    InitializePolynomialBasis_Legendre
  USE InputOutputModuleAMReX_MHD,           ONLY: &
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
  USE MF_MHD_UtilitiesModule,         ONLY: &
    MF_ComputeFromConserved
  USE MF_GeometryModule,                ONLY: &
    MF_ComputeGeometryX
  USE MF_MHD_SlopeLimiterModule, ONLY: &
    MF_ApplySlopeLimiter_MHD
  USE MF_InitializationModule,          ONLY: &
    MF_InitializeFields
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_InitializeMagnetofluid_SSPRK
  USE MF_FieldsModule_MHD,                  ONLY: &
    MF_uGF, &
    MF_uCM, &
    MF_uPM, &
    MF_uAM, &
    MF_uDM
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
    SlopeTolerance,            &
    BetaTVD,                   &
    BetaTVB,                   &
    UseConservativeCorrection, &
    Gamma_IDEAL,               &
    EquationOfState,           &
    StepNo,                    &
    nLevels,                   &
    iRestart,                  &
    MaxGridSizeX,              &
    BA,                        &
    DM,                        &
    GEOM,                      &
    InitializeParameters

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

  LOGICAL, PUBLIC :: wrt, chk


CONTAINS


  SUBROUTINE InitializeProgram

    INTEGER               :: iLevel, iDim
    TYPE(amrex_parmparse) :: PP
    TYPE(amrex_box)       :: BX

    ! --- Initialize AMReX ---

    CALL amrex_init()

    CALL amrex_amrcore_init()

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
               ( MF_uCM(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCM, swX )
        CALL MF_uCM(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uPM(iLevel), BA(iLevel), DM(iLevel), nDOFX * nPM, swX )
        CALL MF_uPM(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uAM(iLevel), BA(iLevel), DM(iLevel), nDOFX * nAM, swX )
        CALL MF_uAM(iLevel) % SetVal( Zero )

        CALL amrex_multifab_build &
               ( MF_uDM(iLevel), BA(iLevel), DM(iLevel), nDOFX * nDM, swX )
        CALL MF_uDM(iLevel) % SetVal( Zero )

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

#if defined HYDRO_RELATIVISTIC

      CALL MF_ComputeGeometryX( MF_uGF, Zero )

#else

      CALL MF_ComputeGeometryX( MF_uGF, Zero )

      IF( SolveGravity ) THEN

        CALL MF_ComputeGravitationalPotential( MF_uGF, Mass )

      END IF

#endif

    CALL SetUnitsGeometryFields

    CALL DescribeFields_Conserved( amrex_parallel_ioprocessor() )

    CALL DescribeFields_Primitive( amrex_parallel_ioprocessor() )

    CALL DescribeFields_Auxiliary( amrex_parallel_ioprocessor() )

    CALL DescribeFields_Diagnostic( amrex_parallel_ioprocessor() )

    CALL SetUnitsFields( TRIM( CoordinateSystem ), &
                              Verbose_Option = amrex_parallel_ioprocessor() )

    CALL InitializeSlopeLimiter_MHD &
           ( BetaTVD_Option &
               = BetaTVD, &
             BetaTVB_Option &
               = BetaTVB, &
             SlopeTolerance_Option &
               = SlopeTolerance, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter, &
             SlopeLimiterMethod_Option &
               = SlopeLimiterMethod, &
             UseConservativeCorrection_Option &
               = UseConservativeCorrection, &
             Verbose_Option &
               = amrex_parallel_ioprocessor() )

    CALL InitializeEquationOfState &
           ( EquationOfState_Option = EquationOfState, &
             Gamma_IDEAL_Option = Gamma_IDEAL, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL MF_InitializeMagnetofluid_SSPRK &
           ( nStages, BA, DM, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    IF( amrex_parallel_ioprocessor() ) WRITE(*,'(A6,A,ES11.3E3)') &
      '', 'CFL: ', &
      CFL * ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

    IF( iRestart .LT. 0 )THEN

      CALL MF_ApplySlopeLimiter_MHD( MF_uGF, MF_uCM, MF_uDM, GEOM )

      CALL MF_InitializeFields( TRIM( ProgramName ), MF_uGF, MF_uCM, GEOM )

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCM_Option = MF_uCM, &
               MF_uPM_Option = MF_uPM, &
               MF_uAM_Option = MF_uAM, &
               MF_uDM_Option = MF_uDM )

    END IF

    DO iLevel = 0, nLevels-1

      CALL amrex_distromap_destroy( DM(iLevel) )
      CALL amrex_boxarray_destroy ( BA(iLevel) )

    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A)') '  Evolving fields...'
      WRITE(*,*)

    END IF

  END SUBROUTINE InitializeProgram


END MODULE InitializationModule
