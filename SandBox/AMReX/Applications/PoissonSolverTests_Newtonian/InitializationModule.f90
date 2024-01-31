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
  USE amrex_fluxregister_module, ONLY: &
    amrex_fluxregister_build, &
    amrex_fluxregister_destroy
  USE amrex_tagbox_module, ONLY: &
    amrex_tagboxarray

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX, &
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
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nPF, &
    nAF, &
    nDF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_EquationOfStateModule, ONLY: &
    InitializeEquationOfState_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    CreateFields_Geometry_MF, &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    CreateFields_Euler_MF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF, &
    FluxRegister_Euler
  USE MF_GeometryModule, ONLY: &
    ComputeGeometryX_MF
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_Euler_MF, &
    ApplySlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_Euler_MF, &
    ApplyPositivityLimiter_Euler_MF
  USE MF_InitializationModule, ONLY: &
    InitializeFields_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_TallyModule, ONLY: &
    InitializeTally_Euler_MF, &
    ComputeTally_Euler_MF
  USE FillPatchModule, ONLY: &
    FillPatch, &
    FillCoarsePatch
  USE TaggingModule, ONLY: &
    TagElements, &
    TagElements_Density
  USE InputParsingModule, ONLY: &
    InitializeParameters, &
    nLevels, &
    nMaxLevels, &
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
    UseFluxCorrection_Euler, &
    TagCriteria, &
    RefinementScheme, &
    DescribeProgramHeader_AMReX
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile
  USE AverageDownModule, ONLY: &
    AverageDown
  USE Euler_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_Euler
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    InitializeTimers_AMReX, &
    Timer_AMReX_Initialize
  USE MF_GravitySolutionModule, ONLY: &
    InitializeGravitySolver_MF
  USE MF_GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    ComputeGravitationalPotential_Newtonian_MF_Poseidon

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
    CALL CreateFields_Euler_MF

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

    CALL InitializeMeshRefinement_Euler

    CALL InitializeEquationOfState_MF

    CALL InitializePositivityLimiter_Euler_MF

    CALL InitializeSlopeLimiter_Euler_MF

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

    CALL amrex_init_from_scratch( t_new(0) )
    nLevels = amrex_get_numlevels()

    SetInitialValues   = .TRUE.

    CALL InitializeTally_Euler_MF

    CALL ApplySlopeLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL CreateMesh_MF( 0, MeshX )

    CALL InitializeGravitySolver_MF &
           ( Verbose_Option = amrex_parallel_ioprocessor() )

    CALL DestroyMesh_MF( MeshX )

    ! Call Newtonian gravity solve here

    Call ComputeGravitationalPotential_Newtonian_MF_Poseidon( MF_uCF, MF_uGF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )
    CALL AverageDown( MF_uGF, MF_uCF )
    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    t_old = t_new
    t_chk = t_new(0) + dt_chk
    t_wrt = t_new(0) + dt_wrt

    CALL DescribeProgramHeader_AMReX

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCF_Option = MF_uCF, &
             MF_uPF_Option = MF_uPF, &
             MF_uAF_Option = MF_uAF, &
             MF_uDF_Option = MF_uDF )

    CALL ComputeTally_Euler_MF &
           ( t_new, MF_uGF, MF_uCF, &
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

    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL MF_uCF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL MF_uPF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL MF_uAF(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL MF_uDF(iLevel) % SetVal( Zero )

    ! Assume nDOFX_X2 = nDOFX_X3 = nDOFX_X1
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nDOFX_X1*nCF )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )

    CALL InitializeFields_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), &
                     MF_uPF(iLevel), MF_uAF(iLevel) )

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
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )

    CALL MF_uGF(iLevel) % SetVal( Zero )
    CALL MF_uCF(iLevel) % SetVal( Zero )
    CALL MF_uDF(iLevel) % SetVal( Zero )
    CALL MF_uPF(iLevel) % SetVal( Zero )
    CALL MF_uAF(iLevel) % SetVal( Zero )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCF )

    CALL FillCoarsePatch( iLevel, MF_uGF, &
                          ApplyBoundaryConditions_Geometry_Option = .TRUE., &
                          UpdateSpatialMetric_Option = .TRUE. )

    CALL FillCoarsePatch( iLevel, MF_uDF )

    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCF, &
                          ApplyBoundaryConditions_Euler_Option = .TRUE. )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uDF(iLevel) )

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uDF(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_destroy( FluxRegister_Euler(iLevel) )

  END SUBROUTINE ClearLevel


  SUBROUTINE RemakeLevel( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER    , INTENT(in), VALUE :: iLevel
    REAL(DP)   , INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM
    TYPE(amrex_multifab)  :: MF_uGF_tmp, MF_uCF_tmp, MF_uDF_tmp

    BA = pBA
    DM = pDM

    CALL amrex_multifab_build( MF_uGF_tmp, BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF_tmp, BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uDF_tmp, BA, DM, nDOFX * nDF, swX )

    CALL MF_uGF_tmp % SetVal( Zero )
    CALL MF_uCF_tmp % SetVal( Zero )
    CALL MF_uDF_tmp % SetVal( Zero )

    CALL FillPatch &
           ( iLevel, MF_uGF, MF_uGF_tmp, &
             ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL FillPatch( iLevel, MF_uDF, MF_uDF_tmp )

    CALL FillPatch &
           ( iLevel, MF_uGF, MF_uGF_tmp, MF_uCF, MF_uCF_tmp, &
             ApplyBoundaryConditions_Euler_Option = .TRUE. )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( iLevel, MF_uGF_tmp, MF_uCF_tmp, MF_uDF_tmp, swX_Option = swX )

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCF )

    CALL MF_uGF(iLevel) % COPY( MF_uGF_tmp, 1, 1, nDOFX * nGF, swX )
    CALL MF_uCF(iLevel) % COPY( MF_uCF_tmp, 1, 1, nDOFX * nCF, swX )
    CALL MF_uDF(iLevel) % COPY( MF_uDF_tmp, 1, 1, nDOFX * nDF, swX )

    CALL amrex_multifab_destroy( MF_uDF_tmp )
    CALL amrex_multifab_destroy( MF_uCF_tmp )
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
    REAL(DP),               CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    CHARACTER(KIND=c_char), CONTIGUOUS, POINTER :: TagArr(:,:,:,:)

    IF( .NOT. ALLOCATED( TagCriteria ) )THEN

       CALL amrex_parmparse_build( PP, "amr" )

         CALL PP % getarr( "TagCriteria", TagCriteria )

       CALL amrex_parmparse_destroy( PP )

    END IF

    Tag = cp

    CALL CreateMesh_MF( iLevel, MeshX )

    !$OMP PARALLEL PRIVATE( MFI, BX, uCF, TagArr )
    CALL amrex_mfiter_build( MFI, MF_uCF( iLevel ), Tiling = UseTiling )

    DO WHILE( MFI % next() )

      BX = MFI % TileBox()

      uCF    => MF_uCF( iLevel ) % DataPtr( MFI )
      TagArr => Tag              % DataPtr( MFI )

      ! TagCriteria(iLevel+1) because iLevel starts at 0 but
      ! TagCriteria starts with 1

      IF( TRIM( RefinementScheme ) .EQ. 'Density' )THEN

        CALL TagElements_Density &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCF ), UBOUND( uCF ), &
                 uCF, TagCriteria(iLevel+1), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) .EQ. 'Mesh' )THEN

        CALL TagElements &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCF ), UBOUND( uCF ), &
                 uCF, TagCriteria(iLevel+1), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE

        CALL TagElements &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCF ), UBOUND( uCF ), &
                 uCF, TagCriteria(iLevel+1), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      END IF

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate


END MODULE InitializationModule
