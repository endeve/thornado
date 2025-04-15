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
  USE thornado_amrex_fluxregister_module, ONLY: &
    amrex_fluxregister_build, &
    amrex_fluxregister_destroy
  USE amrex_tagbox_module, ONLY: &
    amrex_tagboxarray

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    swX, &
    swE, &
    zoomE, &
    nDOFX, &
    nDOFZ, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1, &
    iZ_B1, &
    iZ_E1, &
    iZ_B0, &
    iZ_E0, &
    nNodesE, &
    eL, &
    eR, &
    nE, &
    nX, &
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
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    nDOF_X1
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    NodesX1, NodesX2, NodesX3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    NodesE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange
  USE UnitsModule, ONLY: &
    DescribeUnitsDisplay, &
    Centimeter, &
    UnitsDisplay
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    CreateMesh
  USE EquationOfStateModule, ONLY: &
    EquationOfState
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    CreateGeometryFieldsE, &
    uGE
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nPF, &
    nAF, &
    nDF
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nPR, &
    nAR, &
    nGR, &
    nSpecies
  USE TwoMoment_OpacityModule, ONLY: &
    CreateOpacities, &
    SetOpacities
  USE OpacityModule_Table, ONLY:   &
    InitializeOpacities_TABLE
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_TimersModule, ONLY: &
    InitializeTimers
  USE Euler_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_Euler
  USE TwoMoment_UtilitiesModule
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_EquationOfStateModule, ONLY: &
    InitializeEquationOfState_MF, &
    EosTableName
  USE MF_FieldsModule_Geometry, ONLY: &
    CreateFields_Geometry_MF, &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    CreateFields_Euler_MF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    CreateFields_TwoMoment_MF, &
    MF_uCR, &
    MF_Permute, &
    MF_uPR, &
    MF_uAR, &
    MF_uGR, &
    FluxRegister_TwoMoment
  USE MF_EquationOfStateModule, ONLY: &
    InitializeEquationOfState_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_TwoMoment_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_UtilitiesModule, ONLY: &
    ComputeFromConserved_TwoMoment_MF

  USE MF_TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Initialize_IMEX_RK_MF

  USE FillPatchModule, ONLY: &
    FillPatch, &
    FillCoarsePatch
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
    UseFluxCorrection_TwoMoment, &
    TagCriteria, &
    OpacityTableName_AbEm, &
    OpacityTableName_Iso, &
    OpacityTableName_NES, &
    OpacityTableName_Pair, &
    IOS_CPP,               &
    DescribeProgramHeader_AMReX
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile
  USE AverageDownModule, ONLY: &
    AverageDown

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

CONTAINS


  SUBROUTINE InitializeProgram


    LOGICAL :: SetInitialValues

    TYPE(amrex_parmparse) :: PP

    REAL(DP) :: R0, kT, Mu0, E0, D_0, Chi, Sigma

    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL InitializeTimers

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

    CALL CreateMesh &
           ( MeshE, nE, nNodesE, swE, eL, eR, zoomOption = zoomE )

    ! --- Ordering of calls is important here ---
    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL InitializeReferenceElementE
    CALL InitializeReferenceElementE_Lagrange

    CALL InitializeReferenceElement
    CALL InitializeReferenceElement_Lagrange

    CALL InitializePolynomialBasisMapping &
           ( NodesE, NodesX1, NodesX2, NodesX3 )

    CALL InitializeMeshRefinement_Euler

    CALL CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = amrex_parallel_ioprocessor() )

    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    CALL InitializeEquationOfState_MF

    !CALL InitializePositivityLimiter_TwoMoment_MF

    !CALL InitializeSlopeLimiter_TwoMoment_MF


    IF( TRIM( EquationOfState ) .EQ. 'TABLE' )THEN

      CALL InitializeOpacities_TABLE &
             ( OpacityTableName_EmAb_Option = OpacityTableName_AbEm, &
               OpacityTableName_Iso_Option  = OpacityTableName_Iso,  &
               OpacityTableName_NES_Option  = OpacityTableName_NES,  &
               OpacityTableName_Pair_Option = OpacityTableName_Pair, &
               EquationOfStateTableName_Option = EosTableName, &
               Verbose_Option = amrex_parallel_ioprocessor())

    ELSE

      CALL CreateMesh_MF( 0, MeshX )

      CALL CreateOpacities &
             ( nX, swX, nE, swE, &
               Verbose_Option = amrex_parallel_ioprocessor() )

      R0    = Zero
      E0    = Zero
      Mu0   = Zero
      kT    = Zero
      D_0   = Zero
      Chi   = Zero
      Sigma = Zero
      CALL amrex_parmparse_build( PP, 'ST' )
        CALL PP % query( 'R0'   , R0    )
        CALL PP % query( 'Mu0'  , Mu0   )
        CALL PP % query( 'E0'   , E0    )
        CALL PP % query( 'kT'   , kT    )
        CALL PP % query( 'D_0'  , D_0   )
        CALL PP % query( 'Chi'  , Chi   )
        CALL PP % query( 'Sigma', Sigma )
      CALL amrex_parmparse_destroy( PP )
      Chi  = Chi  * ( One / Centimeter )
      E0   = E0   * UnitsDisplay % EnergyUnit
      mu0  = mu0  * UnitsDisplay % EnergyUnit
      kT   = kT   * UnitsDisplay % EnergyUnit
      R0   = R0   * UnitsDisplay % LengthX1Unit

      CALL SetOpacities &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D_0, Chi, Sigma, Verbose_Option = amrex_parallel_ioprocessor() )

      CALL DestroyMesh_MF( MeshX )

    END IF





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

    CALL Initialize_IMEX_RK_MF &
           ( Verbose_Option = amrex_parallel_ioprocessor() )

    StepNo = 0
    dt     = 0.0_DP
    t_new  = 0.0_DP
    IF( iRestart .LT. 0 )THEN

      CALL amrex_init_from_scratch( 0.0_DP )

      nLevels = amrex_get_numlevels()

      SetInitialValues = .TRUE.

    ELSE

      CALL ReadCheckpointFile &
             ( ReadFields_uCF_Option = .TRUE., &
               ReadFields_uCR_Option = .TRUE. )

      SetInitialValues = .FALSE.

    END IF

    t_old = t_new
    t_chk = t_new(0) + dt_chk
    t_wrt = t_new(0) + dt_wrt

    CALL DescribeProgramHeader_AMReX

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CALL ComputeFromConserved_TwoMoment_MF &
           (  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )

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
    CALL ShowVariableFromMultifab(MF_uGR, 1, writetofile_option=.TRUE., FileNameBase_Option ='Gray_Variables')
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

    INTEGER :: iLo_MF(4)
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

    CALL amrex_multifab_build &
           ( MF_uCR(iLevel), BA, DM, &
             nDOFZ * nCR * ( iE_E0 - iE_B0 + 1 ) * nSpecies, swX )
    CALL MF_uCR(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, &
             nDOFZ * nPR * ( iE_E0 - iE_B0 + 1 ) * nSpecies, swX )
    CALL MF_uPR(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build &
           ( MF_uAR(iLevel), BA, DM, &
             nDOFZ * nAR * ( iE_E0 - iE_B0 + 1 ) * nSpecies, swX )
    CALL MF_uAR(iLevel) % SetVal ( Zero )

    CALL amrex_multifab_build &
           ( MF_Permute(iLevel), BA, DM, &
             nDOFZ * nCR * ( iE_E0 - iE_B0 + 1 ) * nSpecies, swX )
    CALL MF_Permute(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uGR(iLevel), BA, DM, nDOFX * nGR * nSpecies, swX )
    CALL MF_uGR(iLevel) % SetVal( Zero )

    ! Assume nDOF_X1 = nDOF_X2 = nDOFX3
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nDOF_X1*nCR*nE*nSpecies )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )

    CALL InitializeFields_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCR(iLevel), MF_uCF(iLevel) )

    CALL FillPatch( iLevel, MF_uGF )
    CALL FillPatch( iLevel, MF_uGF, MF_uCF )
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
    CALL amrex_multifab_build &
           ( MF_Permute(iLevel), BA, DM, nDOFZ * nCR * nE * nSpecies, swX )
    CALL amrex_multifab_build( MF_uGR(iLevel), BA, DM, nDOFX * nGR * nSpecies, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), &
               iLevel, nDOF_X1 * nCR * nE * nSpecies )

    CALL FillCoarsePatch( iLevel, MF_uGF, &
                          ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL FillCoarsePatch( iLevel, MF_uDF )

    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCF, &
                          ApplyBoundaryConditions_Euler_Option = .TRUE. )

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uPR(iLevel) )
    CALL amrex_multifab_destroy( MF_uCR(iLevel) )
    CALL amrex_multifab_destroy( MF_uGR(iLevel) )
    CALL amrex_multifab_destroy( MF_uDF(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_destroy( FluxRegister_TwoMoment(iLevel) )

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

    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, &
                    ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL FillPatch( iLevel, MF_uDF, MF_uDF_tmp )

    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, MF_uCF, MF_uCF_tmp, &
                    ApplyBoundaryConditions_Euler_Option = .TRUE. )


    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uDF(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL amrex_multifab_build &
           ( MF_uCR(iLevel), BA, DM, nDOFZ * nCR * nE * nSpecies, swX )
    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, nDOFZ * nPR * nE * nSpecies, swX )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), &
               iLevel, nDOF_X1 * nCR * nE * nSpecies )

    CALL MF_uGF(iLevel) % COPY( MF_uGF_tmp, 1, 1, nDOFX * nGF, swX )
    CALL MF_uCF(iLevel) % COPY( MF_uCF_tmp, 1, 1, nDOFX * nCF, swX )
    CALL MF_uDF(iLevel) % COPY( MF_uDF_tmp, 1, 1, nDOFX * nDF, swX )


    CALL amrex_multifab_destroy( MF_uPR_tmp )
    CALL amrex_multifab_destroy( MF_uCR_tmp )
    CALL amrex_multifab_destroy( MF_uDF_tmp )
    CALL amrex_multifab_destroy( MF_uAF_tmp )
    CALL amrex_multifab_destroy( MF_uPF_tmp )
    CALL amrex_multifab_destroy( MF_uCF_tmp )
    CALL amrex_multifab_destroy( MF_uGF_tmp )

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
    CHARACTER(KIND=c_char), CONTIGUOUS, POINTER :: TagArr(:,:,:,:)
    IF( .NOT. ALLOCATED( TagCriteria ) )THEN

       CALL amrex_parmparse_build( PP, "amr" )

         CALL PP % getarr( "TagCriteria", TagCriteria )

       CALL amrex_parmparse_destroy( PP )

    END IF

    Tag = cp

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uCF( iLevel ), Tiling = UseTiling )

    DO WHILE( MFI % next() )

      BX = MFI % TileBox()

      TagArr => Tag              % DataPtr( MFI )

      ! TagCriteria(iLevel+1) because iLevel starts at 0 but
      ! TagCriteria starts with 1

      CALL TagElements &
             ( iLevel, BX % lo, BX % hi, &
               TagCriteria(iLevel+1), SetTag, ClearTag, &
               LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate

END MODULE InitializationModule
