MODULE InitializationModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_init
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init, &
    amrex_init_virtual_functions, &
    amrex_init_from_scratch, &
    amrex_ref_ratio, &
    amrex_get_numlevels
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE thornado_amrex_fluxregister_module, ONLY: &
    amrex_fluxregister_build, &
    amrex_fluxregister_destroy
  USE amrex_tagbox_module, ONLY: &
    amrex_tagboxarray

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    DescribeProgramHeaderX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre, ONLY: &
    InitializePolynomialBasis_Legendre
  USE PolynomialBasisModuleX_Lagrange, ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre, ONLY: &
    InitializePolynomialBasisX_Legendre
  USE PolynomialBasisMappingModule, ONLY: &
    InitializePolynomialBasisMapping
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE UnitsModule, ONLY: &
    DescribeUnitsDisplay
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    nPM, &
    nAM, &
    nDM

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_FieldsModule_Geometry, ONLY: &
    CreateFields_Geometry_MF, &
    MF_uGF
  USE MF_GeometryModule, ONLY: &
    ComputeGeometryX_MF
  USE MF_FieldsModule_MHD, ONLY: &
    CreateFields_MHD_MF, &
    MF_uCM, &
    MF_uPM, &
    MF_uAM, &
    MF_uDM, &
    FluxRegister_MHD
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD_MF
  USE MF_EquationOfStateModule_MHD, ONLY: &
    InitializeEquationOfState_MF
  USE MF_MHD_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_MHD_MF, &
    ApplySlopeLimiter_MHD_MF
  USE MF_MHD_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_MHD_MF, &
    ApplyPositivityLimiter_MHD_MF
  USE MF_TimeSteppingModule_SSPRK_MHD, ONLY: &
    InitializeFluid_SSPRK_MF
  USE MF_InitializationModule, ONLY: &
    InitializeFields_MF
  USE MF_MHD_UtilitiesModule, ONLY: &
    ComputeFromConserved_MHD_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_MHD_TallyModule, ONLY: &
    InitializeTally_MHD_MF, &
    ComputeTally_MHD_MF
  USE FillPatchModule_MHD, ONLY: &
    FillPatch, &
    FillCoarsePatch
  USE MF_XCFC_UtilitiesModule, ONLY: &
    MultiplyWithPsi6_MF
  USE TaggingModule, ONLY: &
    TagElements_Advection1D, &
    TagElements_RiemannProblem1D, &
    TagElements_RiemannProblem2D, &
    TagElements_Advection2D, &
    TagElements_KelvinHelmholtz2D, &
    TagElements_Advection3D, &
    TagElements_uCM
  USE InputParsingModule, ONLY: &
    InitializeParameters, &
    nLevels, &
    nMaxLevels, &
    swX, &
    StepNo, &
    iRestart, &
    dt, &
    t_old, &
    t_new, &
    t_wrt, &
    t_chk, &
    dt_wrt, &
    dt_chk, &
    UseTiling, &
    UseFluxCorrection_MHD, &
    TagCriteria, &
    DescribeProgramHeader_AMReX
  USE InputOutputModuleAMReX_MHD, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile
  USE AverageDownModule_MHD, ONLY: &
    AverageDown
  USE MHD_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_MHD
  USE MF_TimersModule_MHD, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    InitializeTimers_AMReX, &
    Timer_AMReX_Initialize

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

CONTAINS


  SUBROUTINE InitializeProgram

    LOGICAL :: SetInitialValues

    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL InitializeTimers_AMReX

    CALL TimersStart_AMReX( Timer_AMReX_Initialize )

    CALL InitializeParameters

    IF( amrex_parallel_ioprocessor() )THEN

      CALL DescribeUnitsDisplay
      CALL DescribeProgramHeaderX

    END IF

    CALL CreateFields_Geometry_MF
    CALL CreateFields_MHD_MF

    CALL InitializePolynomialBasisX_Lagrange
    CALL InitializePolynomialBasisX_Legendre

    CALL InitializePolynomialBasis_Lagrange
    CALL InitializePolynomialBasis_Legendre

    CALL CreateMesh_MF( 0, MeshX )

    CALL InitializePolynomialBasisMapping &
           ( [Zero], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    CALL DestroyMesh_MF( MeshX )

    ! --- Ordering of calls is important here ---
    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL InitializeMeshRefinement_MHD

    CALL InitializeEquationOfState_MF

    CALL InitializePositivityLimiter_MHD_MF

    CALL InitializeSlopeLimiter_MHD_MF

    CALL amrex_init_virtual_functions &
           ( MakeNewLevelFromScratch, &
             MakeNewLevelFromCoarse, &
             RemakeLevel, &
             ClearLevel, &
             ErrorEstimate )

    ALLOCATE( StepNo(0:nMaxLevels-1) )
    ALLOCATE( dt    (0:nMaxLevels-1) )
    ALLOCATE( t_old (0:nMaxLevels-1) )
    ALLOCATE( t_new (0:nMaxLevels-1) )

    StepNo = 0
    dt     = 0.0_DP
    t_new  = 0.0_DP

    IF( iRestart .LT. 0 )THEN

      CALL amrex_init_from_scratch( Zero )
      nLevels = amrex_get_numlevels()

      SetInitialValues = .TRUE.

      CALL InitializeTally_MHD_MF

      CALL ApplySlopeLimiter_MHD_MF &
             ( t_new, MF_uGF, MF_uCM, MF_uDM )

      CALL ApplyPositivityLimiter_MHD_MF &
             ( t_new, MF_uGF, MF_uCM, MF_uDM )

    ELSE

      CALL ReadCheckpointFile( ReadFields_uCM_Option = .TRUE. )

      SetInitialValues = .FALSE.

      CALL InitializeTally_MHD_MF &
             ( InitializeFromCheckpoint_Option = .TRUE. )

    END IF

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )
    CALL AverageDown( MF_uGF, MF_uCM )
    CALL ApplyPositivityLimiter_MHD_MF &
           ( t_new, MF_uGF, MF_uCM, MF_uDM )

    t_old = t_new
    t_chk = t_new(0) + dt_chk
    t_wrt = t_new(0) + dt_wrt

    CALL InitializeFluid_SSPRK_MF &
           ( Verbose_Option = amrex_parallel_ioprocessor() )

    CALL DescribeProgramHeader_AMReX

    CALL ComputeFromConserved_MHD_MF &
           ( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCM_Option = MF_uCM, &
             MF_uPM_Option = MF_uPM, &
             MF_uAM_Option = MF_uAM, &
             MF_uDM_Option = MF_uDM )

    CALL ComputeTally_MHD_MF &
           ( t_new, MF_uGF, MF_uCM, MF_uAM, &
             SetInitialValues_Option = SetInitialValues, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL TimersStop_AMReX( Timer_AMReX_Initialize )

  END SUBROUTINE InitializeProgram


  SUBROUTINE MakeNewLevelFromScratch( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM

    BA = pBA
    DM = pDM

    t_new(iLevel) = Time
    t_old(iLevel) = Time - 1.0e200_DP

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL MF_uGF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uCM(iLevel), BA, DM, nDOFX * nCM, swX )
    CALL MF_uCM(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uPM(iLevel), BA, DM, nDOFX * nPM, swX )
    CALL MF_uPM(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uAM(iLevel), BA, DM, nDOFX * nAM, swX )
    CALL MF_uAM(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uDM(iLevel), BA, DM, nDOFX * nDM, swX )
    CALL MF_uDM(iLevel) % SetVal( Zero )

    ! Assume nDOFX_X2 = nDOFX_X3 = nDOFX_X1
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_MHD ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_MHD(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nDOFX_X1*nCM )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )

    CALL InitializeFields_MF( iLevel, MF_uGF(iLevel), MF_uCM(iLevel), MF_uDM(iLevel) )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE MakeNewLevelFromScratch


  SUBROUTINE MakeNewLevelFromCoarse( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM

    BA = pBA
    DM = pDM

    CALL ClearLevel( iLevel )

    t_new( iLevel ) = Time
    t_old( iLevel ) = Time - 1.0e200_DP

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCM(iLevel), BA, DM, nDOFX * nCM, swX )
    CALL amrex_multifab_build( MF_uDM(iLevel), BA, DM, nDOFX * nDM, swX )
    CALL amrex_multifab_build( MF_uPM(iLevel), BA, DM, nDOFX * nPM, swX )
    CALL amrex_multifab_build( MF_uAM(iLevel), BA, DM, nDOFX * nAM, swX )

    CALL MF_uGF(iLevel) % SetVal( Zero )
    CALL MF_uCM(iLevel) % SetVal( Zero )
    CALL MF_uDM(iLevel) % SetVal( Zero )
    CALL MF_uPM(iLevel) % SetVal( Zero )
    CALL MF_uAM(iLevel) % SetVal( Zero )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_MHD ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_MHD(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCM )

    CALL FillCoarsePatch( iLevel, MF_uGF, &
                          ApplyBoundaryConditions_Geometry_Option = .TRUE., &
                          UpdateSpatialMetric_Option = .TRUE. )

    CALL FillCoarsePatch( iLevel, MF_uDM )

    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCM )

    CALL ApplyBoundaryConditions_MHD_MF( t_new(iLevel), iLevel, MF_uCM(iLevel), MF_uDM(iLevel) )

    CALL MultiplyWithPsi6_MF( MF_uGF(iLevel), MF_uCM(iLevel), -1 )

    CALL ApplyPositivityLimiter_MHD_MF &
           ( Time, iLevel, MF_uGF(iLevel), MF_uCM(iLevel), MF_uDM(iLevel) )

    CALL MultiplyWithPsi6_MF( MF_uGF(iLevel), MF_uCM(iLevel), +1 )

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uDM(iLevel) )
    CALL amrex_multifab_destroy( MF_uAM(iLevel) )
    CALL amrex_multifab_destroy( MF_uPM(iLevel) )
    CALL amrex_multifab_destroy( MF_uCM(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_MHD ) &
      CALL amrex_fluxregister_destroy( FluxRegister_MHD(iLevel) )

  END SUBROUTINE ClearLevel


  SUBROUTINE RemakeLevel( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER    , INTENT(in), VALUE :: iLevel
    REAL(DP)   , INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM
    TYPE(amrex_multifab)  :: MF_uGF_tmp, MF_uCM_tmp, MF_uDM_tmp

    BA = pBA
    DM = pDM

    CALL amrex_multifab_build( MF_uGF_tmp, BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCM_tmp, BA, DM, nDOFX * nCM, swX )
    CALL amrex_multifab_build( MF_uDM_tmp, BA, DM, nDOFX * nDM, swX )

    CALL MF_uGF_tmp % SetVal( Zero )
    CALL MF_uCM_tmp % SetVal( Zero )
    CALL MF_uDM_tmp % SetVal( Zero )

    CALL FillPatch &
           ( iLevel, MF_uGF, MF_uGF_tmp, &
             ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL FillPatch( iLevel, MF_uDM, MF_uDM_tmp )

    CALL FillPatch &
           ( iLevel, MF_uGF, MF_uGF_tmp, MF_uCM, MF_uCM_tmp )

    CALL ApplyBoundaryConditions_MHD_MF( Time, iLevel, MF_uCM_tmp, MF_uDM_tmp )

    CALL MultiplyWithPsi6_MF( MF_uGF_tmp, MF_uCM_tmp, -1, swX_Option = swX )

    CALL ApplyPositivityLimiter_MHD_MF &
           ( Time, iLevel, MF_uGF_tmp, MF_uCM_tmp, MF_uDM_tmp, swX_Option = swX )

    CALL MultiplyWithPsi6_MF( MF_uGF_tmp, MF_uCM_tmp, +1, swX_Option = swX )

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCM(iLevel), BA, DM, nDOFX * nCM, swX )
    CALL amrex_multifab_build( MF_uDM(iLevel), BA, DM, nDOFX * nDM, swX )
    CALL amrex_multifab_build( MF_uPM(iLevel), BA, DM, nDOFX * nPM, swX )
    CALL amrex_multifab_build( MF_uAM(iLevel), BA, DM, nDOFX * nAM, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_MHD ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_MHD(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCM )

    CALL MF_uGF(iLevel) % COPY( MF_uGF_tmp, 1, 1, nDOFX * nGF, swX )
    CALL MF_uCM(iLevel) % COPY( MF_uCM_tmp, 1, 1, nDOFX * nCM, swX )
    CALL MF_uDM(iLevel) % COPY( MF_uDM_tmp, 1, 1, nDOFX * nDM, swX )

    CALL amrex_multifab_destroy( MF_uDM_tmp )
    CALL amrex_multifab_destroy( MF_uCM_tmp )
    CALL amrex_multifab_destroy( MF_uGF_tmp )

  END SUBROUTINE RemakeLevel


  SUBROUTINE ErrorEstimate( iLevel, cp, Time, SetTag, ClearTag ) BIND(c)

    INTEGER,                INTENT(in), VALUE :: iLevel
    TYPE(c_ptr),            INTENT(in), VALUE :: cp
    REAL(DP),               INTENT(in), VALUE :: Time
    CHARACTER(KIND=c_char), INTENT(in), VALUE :: SetTag, ClearTag

    TYPE(amrex_parmparse)   :: PP
    TYPE(amrex_tagboxarray) :: Tag
    TYPE(amrex_mfiter)      :: MFI
    TYPE(amrex_box)         :: BX
    REAL(DP),               CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    CHARACTER(KIND=c_char), CONTIGUOUS, POINTER :: TagArr(:,:,:,:)

    IF( .NOT. ALLOCATED( TagCriteria ) )THEN

       CALL amrex_parmparse_build( PP, "amr" )

         CALL PP % getarr( "TagCriteria", TagCriteria )

       CALL amrex_parmparse_destroy( PP )

    END IF

    Tag = cp

    CALL CreateMesh_MF( iLevel, MeshX )

    !$OMP PARALLEL PRIVATE( MFI, BX, uCM, TagArr )
    CALL amrex_mfiter_build( MFI, MF_uCM( iLevel ), Tiling = UseTiling )

    DO WHILE( MFI % next() )

      BX = MFI % TileBox()

      uCM    => MF_uCM( iLevel ) % DataPtr( MFI )
      TagArr => Tag              % DataPtr( MFI )

      ! TagCriteria(iLevel+1) because iLevel starts at 0 but
      ! TagCriteria starts with 1

      SELECT CASE( TRIM( ProgramName ) )

        CASE( 'Advection1D' )

          CALL TagElements_Advection1D &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'RiemannProblem1D' )

          CALL TagElements_RiemannProblem1D &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'RiemannProblem2D' )

          CALL TagElements_RiemannProblem2D &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'Advection2D' )

          CALL TagElements_Advection2D &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'KelvinHelmholtz2D' )

          CALL TagElements_KelvinHelmholtz2D &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'Advection3D' )

          CALL TagElements_Advection3D &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE DEFAULT

          CALL TagElements_uCM &
                 ( iLevel, BX % lo, BX % hi, LBOUND( uCM ), UBOUND( uCM ), &
                   uCM, TagCriteria(iLevel+1), SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      END SELECT

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate


END MODULE InitializationModule
