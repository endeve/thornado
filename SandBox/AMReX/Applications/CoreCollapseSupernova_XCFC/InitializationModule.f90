MODULE InitializationModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_init_module, ONLY: &
    amrex_init
  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init, &
    amrex_init_virtual_functions, &
    amrex_init_from_scratch, &
    amrex_ref_ratio, &
    amrex_get_numlevels, &
    amrex_geom
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
    nDOFZ, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1, &
    nNodesX, &
    nNodesE, &
    DescribeProgramHeaderX
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    InitializeNeutrinoMatterSolverParameters
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
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange
  USE ReferenceElementModuleZ, ONLY: &
    InitializeReferenceElementZ, &
    nDOFZ_Z2
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange
  USE UnitsModule, ONLY: &
    DescribeUnitsDisplay
  USE MeshModule, ONLY: &
    CreateMesh, &
    MeshX, &
    MeshE
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    CoordinateSystem, &
    DescribeGeometryFields, &
    SetUnitsGeometryFields
  USE GeometryFieldsModuleE, ONLY: &
    CreateGeometryFieldsE, &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nPF, &
    nAF, &
    nDF, &
    DescribeFluidFields_Primitive, &
    DescribeFluidFields_Conserved, &
    DescribeFluidFields_Auxiliary, &
    DescribeFluidFields_Diagnostic, &
    SetUnitsFluidFields
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nPR, &
    nGR, &
    DescribeRadiationFields_Primitive, &
    DescribeRadiationFields_Conserved, &
    SetUnitsRadiationFields
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE OpacityModule_Table, ONLY: &
    InitializeOpacities_TABLE
  USE TwoMoment_TimersModule, ONLY: &
    InitializeTimers

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
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
  USE MF_FieldsModule_TwoMoment, ONLY: &
    CreateFields_TwoMoment_MF, &
    MF_uCR, &
    MF_uPR, &
    MF_uGR, &
    FluxRegister_TwoMoment
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_Euler_MF, &
    ApplySlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_Euler_MF, &
    ApplyPositivityLimiter_Euler_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_TwoMoment_MF, &
    ApplySlopeLimiter_TwoMoment_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment_MF, &
    ApplyPositivityLimiter_TwoMoment_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_TwoMoment_UtilitiesModule, ONLY: &
    ComputeGray_TwoMoment_MF, &
    ComputeFromConserved_TwoMoment_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_TallyModule, ONLY: &
    InitializeTally_Euler_MF, &
    ComputeTally_Euler_MF
  USE MF_TwoMoment_TallyModule, ONLY: &
    InitializeTally_TwoMoment_MF, &
    ComputeTally_TwoMoment_MF
  USE MF_TimeSteppingModule_IMEX, ONLY: &
    Initialize_IMEX_RK_MF
  USE FillPatchModule, ONLY: &
    FillPatch, &
    FillCoarsePatch
  USE InputParsingModule, ONLY: &
    InitializeParameters, &
    nLevels, &
    nMaxLevels, &
    swX, &
    swE, &
    zoomE, &
    StepNo, &
    iRestart, &
    dt, &
    dt_TM, &
    t_old, &
    t_new, &
    t_wrt, &
    t_chk, &
    dt_wrt, &
    dt_chk, &
    UseTiling, &
    UseFluxCorrection_Euler, &
    UseFluxCorrection_TwoMoment, &
    MaxGridSizeX, &
    BlockingFactor, &
    xL, &
    xR, &
    eL, &
    eR, &
    nE, &
    OpacityTableName_AbEm, &
    OpacityTableName_Iso, &
    OpacityTableName_NES, &
    OpacityTableName_Pair, &
    OpacityTableName_Brem, &
    M_outer, &
    MaxIter_outer, &
    Rtol_outer, &
    M_inner, &
    MaxIter_inner, &
    Rtol_inner, &
    Include_NES, &
    Include_Pair, &
    Include_Brem, &
    Include_LinCorr, &
    wMatterRHS, &
    Scheme, &
    nSpecies, &
    EquationOfState, &
    Gamma_IDEAL, &
    EosTableName, &
    ProgramName, &
    TagCriteria, &
    nRefinementBuffer, &
    UseAMR, &
    DescribeProgramHeader_AMReX
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile
  USE AverageDownModule, ONLY: &
    AverageDown
  USE Euler_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_Euler
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    InitializeGravitySolver_XCFC_MF, &
    InitializeMetric_TwoMoment_MF
  USE MF_TimersModule, ONLY: &
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
    CALL CreateFields_Euler_MF
    CALL CreateFields_TwoMoment_MF

    CALL InitializePolynomialBasisX_Lagrange
    CALL InitializePolynomialBasisX_Legendre

    CALL InitializePolynomialBasis_Lagrange
    CALL InitializePolynomialBasis_Legendre

    CALL CreateMesh_MF( 0, MeshX )

    CALL CreateMesh &
           ( MeshE, nE, nNodesE, swE, eL, eR, zoomOption = zoomE )

    CALL InitializePolynomialBasisMapping &
           ( MeshE % Nodes, &
             MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    CALL DestroyMesh_MF( MeshX )

    ! --- Ordering of calls is important here ---
    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL InitializeReferenceElementE
    CALL InitializeReferenceElementE_Lagrange

    CALL InitializeReferenceElementZ

    CALL InitializeReferenceElement
    CALL InitializeReferenceElement_Lagrange

    CALL InitializeMeshRefinement_Euler

    CALL SetUnitsGeometryFields

    CALL DescribeFluidFields_Conserved ( amrex_parallel_ioprocessor() )

    CALL DescribeFluidFields_Primitive ( amrex_parallel_ioprocessor() )

    CALL DescribeFluidFields_Auxiliary ( amrex_parallel_ioprocessor() )

    CALL DescribeFluidFields_Diagnostic( amrex_parallel_ioprocessor() )

    CALL SetUnitsFluidFields &
           ( TRIM( CoordinateSystem ), &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL DescribeRadiationFields_Conserved( amrex_parallel_ioprocessor() )

    CALL DescribeRadiationFields_Primitive( amrex_parallel_ioprocessor() )

    CALL SetUnitsRadiationFields

    CALL CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = amrex_parallel_ioprocessor() )

    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    CALL InitializeEquationOfState &
           ( EquationOfState_Option = EquationOfState, &
             EquationOfStateTableName_Option = EosTableName, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL InitializeOpacities_TABLE &
           ( OpacityTableName_EmAb_Option = OpacityTableName_AbEm, &
             OpacityTableName_Iso_Option  = OpacityTableName_Iso,  &
             OpacityTableName_NES_Option  = OpacityTableName_NES,  &
             OpacityTableName_Pair_Option = OpacityTableName_Pair, &
             OpacityTableName_Brem_Option = OpacityTableName_Brem, &
             EquationOfStateTableName_Option = EosTableName, &
             Verbose_Option = amrex_parallel_ioprocessor())

    CALL InitializeNeutrinoMatterSolverParameters &
           ( M_outer_Option &
               = M_outer, &
             M_inner_Option &
               = M_inner, &
             MaxIter_outer_Option &
               = MaxIter_outer, &
             MaxIter_inner_Option &
               = MaxIter_inner, &
             Rtol_inner_Option &
               = Rtol_inner, &
             Rtol_outer_Option &
               = Rtol_outer, &
             Include_NES_Option &
               = Include_NES, &
             Include_Pair_Option &
               = Include_Pair, &
             Include_Brem_Option &
               = Include_Brem, &
             Include_LinCorr_Option &
               = Include_LinCorr, &
             wMatrRHS_Option &
               = wMatterRHS, &
             Verbose_Option &
               = amrex_parallel_ioprocessor() )

    CALL InitializeClosure_TwoMoment &
           ( Verbose_Option = amrex_parallel_ioprocessor() )

    CALL InitializePositivityLimiter_Euler_MF

    CALL InitializeSlopeLimiter_Euler_MF

    CALL InitializePositivityLimiter_TwoMoment_MF

    CALL InitializeSlopeLimiter_TwoMoment_MF

    CALL amrex_init_virtual_functions &
           ( MakeNewLevelFromScratch, &
             MakeNewLevelFromCoarse, &
             RemakeLevel, &
             ClearLevel, &
             ErrorEstimate )

    ALLOCATE( StepNo(0:nMaxLevels-1) )
    ALLOCATE( dt    (0:nMaxLevels-1) )
    ALLOCATE( dt_TM (0:nMaxLevels-1) )
    ALLOCATE( t_old (0:nMaxLevels-1) )
    ALLOCATE( t_new (0:nMaxLevels-1) )

    StepNo = 0
    dt     = 0.0_DP
    dt_TM  = 0.0_DP
    t_new  = 0.0_DP

    CALL InitializeTally_TwoMoment_MF

    CALL Initialize_IMEX_RK_MF &
           ( Scheme, &
             EvolveEuler_Option     = .TRUE., &
             EvolveTwoMoment_Option = .TRUE., &
             Verbose_Option         = amrex_parallel_ioprocessor() )

    IF( iRestart .LT. 0 )THEN

      CALL amrex_init_from_scratch( 0.0_DP )
      nLevels = amrex_get_numlevels()

      SetInitialValues = .TRUE.

      CALL InitializeTally_Euler_MF

      CALL ApplySlopeLimiter_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uDF )

      CALL ApplyPositivityLimiter_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uDF )

      CALL ApplySlopeLimiter_TwoMoment_MF &
             ( amrex_geom, MF_uGF, MF_uCF, MF_uCR )

      CALL ApplyPositivityLimiter_TwoMoment_MF &
             ( amrex_geom, MF_uGF, MF_uCF, MF_uCR )

      CALL CreateMesh_MF( 0, MeshX )

      CALL InitializeGravitySolver_XCFC_MF

      CALL DestroyMesh_MF( MeshX )

      CALL ComputeFromConserved_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL ComputeFromConserved_TwoMoment_MF &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

      CALL InitializeMetric_TwoMoment_MF &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

      CALL ApplySlopeLimiter_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uDF )

      CALL ApplyPositivityLimiter_Euler_MF &
             ( MF_uGF, MF_uCF, MF_uDF )

      CALL ApplySlopeLimiter_TwoMoment_MF &
             ( amrex_geom, MF_uGF, MF_uCF, MF_uCR )

      CALL ApplyPositivityLimiter_TwoMoment_MF &
             ( amrex_geom, MF_uGF, MF_uCF, MF_uCR )

    ELSE

      CALL ReadCheckpointFile &
             ( ReadFields_uCF_Option = .TRUE., &
               ReadFields_uCR_Option = .TRUE. )

      SetInitialValues = .FALSE.

      CALL InitializeTally_Euler_MF &
             ( InitializeFromCheckpoint_Option = .TRUE. )

      CALL CreateMesh_MF( 0, MeshX )

      CALL InitializeGravitySolver_XCFC_MF( MF_uGF, MF_uCF )

      CALL DestroyMesh_MF( MeshX )

    END IF

    CALL AverageDown( MF_uGF )
    CALL AverageDown &
           ( MF_uGF, MF_uCF, MF_uDF, ApplyPositivityLimiter_Option = .TRUE. )
    !!$CALL AverageDown( MF_uGF, MF_uCR )

    t_old = t_new
    t_chk = t_new(0) + dt_chk
    t_wrt = t_new(0) + dt_wrt

    CALL DescribeProgramHeader_AMReX

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CALL ComputeFromConserved_TwoMoment_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

    CALL ComputeGray_TwoMoment_MF &
           ( MF_uGF, MF_uPF, MF_uCR, MF_uPR, MF_uGR )

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
           ( t_new, MF_uGF, MF_uCF, &
             SetInitialValues_Option = SetInitialValues, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL ComputeTally_TwoMoment_MF &
           ( amrex_geom, MF_uGF, MF_uCF, MF_uCR, t_new(0), &
             SetInitialValues_Option = .TRUE., &
             Verbose_Option = amrex_parallel_ioprocessor() )

    CALL TimersStop_AMReX( Timer_AMReX_Initialize )

  END SUBROUTINE InitializeProgram


  SUBROUTINE MakeNewLevelFromScratch( iLevel, Time, pBA, pDM ) BIND(c)

    USE MF_GeometryModule, ONLY: &
      ComputeGeometryX_MF

    USE MF_InitializationModule, ONLY: &
      InitializeFields_MF

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

    CALL amrex_multifab_build &
           ( MF_uCR(iLevel), BA, DM, nDOFZ * nCR * nE * nSpecies, swX )
    CALL MF_uCR(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, nDOFZ * nPR * nE * nSpecies, swX )
    CALL MF_uPR(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uGR(iLevel), BA, DM, nDOFX * nGR * nSpecies, swX )
    CALL MF_uGR(iLevel) % SetVal( Zero )

    ! Assume nDOFZ_Z3 = nDOFZ_Z4 = nDOFZ_Z2
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nDOFZ_Z2*nCR*nE*nSpecies )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )

    CALL InitializeFields_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uCR(iLevel) )

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
    CALL amrex_multifab_build &
           ( MF_uCR(iLevel), BA, DM, nDOFZ * nCR * nE * nSpecies, swX )
    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, nDOFZ * nPR * nE * nSpecies, swX )
    CALL amrex_multifab_build( MF_uGR(iLevel), BA, DM, nDOFX * nGR * nSpecies, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCF )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFZ_Z2 * nCR * nE * nSpecies )

    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uGF )
    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCF )
    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uDF )
    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCR )

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uCR(iLevel) )
    CALL amrex_multifab_destroy( MF_uPR(iLevel) )
    CALL amrex_multifab_destroy( MF_uGR(iLevel) )
    CALL amrex_multifab_destroy( MF_uDF(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_destroy( FluxRegister_Euler(iLevel) )

  END SUBROUTINE ClearLevel


  SUBROUTINE RemakeLevel( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM
    TYPE(amrex_multifab)  :: MF_uGF_tmp, MF_uCF_tmp, MF_uPF_tmp, &
                             MF_uAF_tmp, MF_uDF_tmp, MF_uCR_tmp, MF_uPR_tmp

    BA = pBA
    DM = pDM

    CALL amrex_multifab_build( MF_uGF_tmp, BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF_tmp, BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF_tmp, BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF_tmp, BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF_tmp, BA, DM, nDOFX * nDF, swX )
    CALL amrex_multifab_build &
           ( MF_uCR_tmp, BA, DM, nDOFZ * nCR * nE * nSpecies, swX )
    CALL amrex_multifab_build &
           ( MF_uPR_tmp, BA, DM, nDOFZ * nPR * nE * nSpecies, swX )

    CALL FillPatch( iLevel, MF_uGF, MF_uGF, MF_uGF_tmp )
    CALL FillPatch( iLevel, MF_uGF, MF_uCF, MF_uCF_tmp )
    CALL FillPatch( iLevel, MF_uGF, MF_uDF, MF_uDF_tmp )
    CALL FillPatch( iLevel, MF_uGF, MF_uCR, MF_uCR_tmp )

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL amrex_multifab_build &
           ( MF_uCR(iLevel), BA, DM, nDOFZ * nCF * nE * nSpecies, swX )
    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, nDOFZ * nPF * nE * nSpecies, swX )
    CALL amrex_multifab_build( MF_uGR(iLevel), BA, DM, nDOFX * nGR * nSpecies, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_Euler ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_Euler(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFX_X1 * nCF )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), &
               iLevel, nDOFZ_Z2 * nCR * nE * nSpecies )

    CALL MF_uGF(iLevel) % COPY( MF_uGF_tmp, 1, 1, nDOFX * nGF, swX )
    CALL MF_uCF(iLevel) % COPY( MF_uCF_tmp, 1, 1, nDOFX * nCF, swX )
    CALL MF_uDF(iLevel) % COPY( MF_uDF_tmp, 1, 1, nDOFX * nDF, swX )
    CALL MF_uCR(iLevel) &
           % COPY( MF_uCR_tmp, 1, 1, nDOFZ * nCR * nE * nSpecies, swX )

    CALL amrex_multifab_destroy( MF_uDF_tmp )
    CALL amrex_multifab_destroy( MF_uAF_tmp )
    CALL amrex_multifab_destroy( MF_uPF_tmp )
    CALL amrex_multifab_destroy( MF_uCF_tmp )
    CALL amrex_multifab_destroy( MF_uGF_tmp )
    CALL amrex_multifab_destroy( MF_uPR_tmp )
    CALL amrex_multifab_destroy( MF_uCR_tmp )

  END SUBROUTINE RemakeLevel


  SUBROUTINE ErrorEstimate( iLevel, cp, Time, SetTag, ClearTag ) BIND(c)

    USE TaggingModule, ONLY: &
      TagElements

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

      CALL TagElements &
             ( iLevel, BX % lo, BX % hi, LBOUND( uCF ), UBOUND( uCF ), &
               uCF, TagCriteria(iLevel+1), SetTag, ClearTag, &
               LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate

END MODULE InitializationModule
