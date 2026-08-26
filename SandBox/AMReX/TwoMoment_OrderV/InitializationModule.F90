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
  USE MF_Euler_TallyModule, ONLY: &
    InitializeTally_Euler_MF, &
    ComputeTally_Euler_MF
  USE MF_TwoMoment_TallyModule

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    swX, &
    swE, &
    zoomE, &
    nDOFX, &
    nDOFZ, &
    nDOFE, &
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
    NodesX1, NodesX2, NodesX3, nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    NodesE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange
  USE ReferenceElementModuleZ, ONLY: &
    InitializeReferenceElementZ
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
    InitializeMeshRefinement_Euler !, &
    !InitializeMeshRefinement_Euler_Aniso
  USE MF_TwoMoment_OpacityModule, ONLY: &
    InitializeOpacities_MF, &
    FinalizeOpacities_MF, &
    BuildOpacities_MF_Level, &
    ClearOpacities_MF_Level
    !SetOpacities_TwoMoment_MF
  USE amrex_amrcore_module, ONLY: &
    GEOM => amrex_geom
  USE TwoMoment_UtilitiesModule
  USE MF_UtilitiesModule, ONLY: &
    ShowVariableFromMultiFab
  !USE AnisotropicRefinementModule, ONLY: &
  !  UseAnisotropicRefinement, &
  !  RefRatioVect, &
  !  ParseAnisotropicRefinement

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
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment_MF
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
  USE MF_TwoMoment_UtilitiesModule_OrderV, ONLY: &
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
    TagCriteria, &
    RefinementScheme, &
    UseFluxCorrection_Euler, &
    UseFluxCorrection_TwoMoment, &
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
  USE MF_GeometryModule, ONLY: &
    ComputeGeometryX_MF
  USE MF_InitializationModule, ONLY: &
      InitializeFields_MF

  !IMPLICIT NONE
  !PRIVATE

IMPLICIT NONE
  PRIVATE
 
  PUBLIC :: InitializeProgram
 
CONTAINS
 
 
  SUBROUTINE InitializeProgram
 
    LOGICAL :: SetInitialValues
 
    TYPE(amrex_parmparse) :: PP

    INTEGER :: iLevel
 
    REAL(DP) :: R0, kT, Mu0, E0, D_0, Chi, Sigma

    CALL amrex_init()
 
    CALL amrex_amrcore_init()
 
    CALL InitializeTimers
 
    CALL InitializeParameters

    !CALL ParseAnisotropicRefinement &
    !       ( UseFluxCorrection_Euler, UseFluxCorrection_TwoMoment )
 
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
 
    CALL InitializeReferenceElementZ
 
    CALL InitializeReferenceElement
    CALL InitializeReferenceElement_Lagrange
 
    CALL InitializePolynomialBasisMapping &
           ( NodesE, NodesX1, NodesX2, NodesX3 )
 
    CALL InitializeMeshRefinement_Euler

    !IF( UseAnisotropicRefinement ) &
    !  CALL InitializeMeshRefinement_Euler_Aniso( RefRatioVect )
 
    CALL CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = amrex_parallel_ioprocessor() )
 
    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )
 
    CALL InitializeEquationOfState_MF
 
    CALL InitializePositivityLimiter_TwoMoment_MF
 
    CALL InitializeSlopeLimiter_TwoMoment_MF
 
    IF( TRIM( EquationOfState ) .EQ. 'TABLE' )THEN
 
      CALL InitializeOpacities_TABLE &
             ( OpacityTableName_EmAb_Option = OpacityTableName_AbEm, &
               OpacityTableName_Iso_Option  = OpacityTableName_Iso,  &
               OpacityTableName_NES_Option  = OpacityTableName_NES,  &
               OpacityTableName_Pair_Option = OpacityTableName_Pair, &
               EquationOfStateTableName_Option = EosTableName, &
               Verbose_Option = amrex_parallel_ioprocessor() )
 
    ELSE
 
      R0    = Zero
      E0    = Zero
      Mu0   = Zero
      kT    = Zero
      D_0   = Zero
      Chi   = Zero
      Sigma = Zero
      CALL amrex_parmparse_build( PP, 'thornado' )
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
      Mu0  = Mu0  * UnitsDisplay % EnergyUnit
      kT   = kT   * UnitsDisplay % EnergyUnit
      R0   = R0   * UnitsDisplay % LengthX1Unit
 
      CALL InitializeOpacities_MF( D_0, Chi, Sigma )
 
    END IF
    ! ================================================================
 
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
 
      CALL amrex_init_from_scratch( 0.0_DP )
 
      nLevels = amrex_get_numlevels()
 
      SetInitialValues = .TRUE.
 
      CALL InitializeTally_Euler_MF
      CALL InitializeTally_TwoMoment_MF
 
    ELSE
 
      CALL ReadCheckpointFile
      SetInitialValues = .FALSE.
 
      CALL InitializeTally_Euler_MF &
         ( InitializeFromCheckpoint_Option = .TRUE. )

      CALL InitializeTally_TwoMoment_MF &
         ( InitializeFromCheckpoint_Option = .TRUE. )
      
      DO iLevel = 0, nLevels - 1
    
         CALL BuildOpacities_MF_Level( iLevel, MF_uGF( iLevel ) % BA, MF_uGF( iLevel ) % DM )
      
      END DO
 
    END IF
 
    t_old = t_new
    t_chk = t_new(0) + dt_chk
    t_wrt = t_new(0) + dt_wrt
 
    CALL Initialize_IMEX_RK_MF &
           ( MF_uGF % BA, MF_uGF % DM, &
             Verbose_Option = amrex_parallel_ioprocessor() )
 
    CALL DescribeProgramHeader_AMReX

    DO iLevel = 0, nLevels-1

      CALL FillPatch( iLevel, MF_uGF, &
                      ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    END DO
 
    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )
 
    CALL ComputeFromConserved_TwoMoment_MF &
           (  MF_uGF, MF_uCF, MF_uCR, MF_uPR, MF_uAR, MF_uGR )
 
    CALL ComputeTally_Euler_MF &
       ( t_new, MF_uGF, MF_uCF, &
         SetInitialValues_Option = SetInitialValues, &
         Verbose_Option = amrex_parallel_ioprocessor() )

    CALL ComputeTally_TwoMoment_MF &
       ( t_new, MF_uGF, MF_uCF, MF_uCR, &
         SetInitialValues_Option = SetInitialValues, &
         Verbose_Option = amrex_parallel_ioprocessor() )
 
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
 
    CALL ShowVariableFromMultiFab( MF_uPR, 1, &
                             WriteToFile_Option  = .TRUE., &
                             FileNameBase_Option = 'MF_uPR' )
 
    CALL ShowVariableFromMultiFab( MF_uCR, 1, &
                             WriteToFile_Option  = .TRUE., &
                             FileNameBase_Option = 'MF_uCR' )
 
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

    CALL amrex_multifab_build &
           ( MF_uCR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uCR(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uPR(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build &
           ( MF_uAR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uAR(iLevel) % SetVal ( Zero )

    CALL amrex_multifab_build &
           ( MF_Permute(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_Permute(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uGR(iLevel), BA, DM, &
             nDOFX *  nGR * nSpecies, swX )
    CALL MF_uGR(iLevel) % SetVal( Zero )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment )THEN
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, &
               nDOFX_X1 * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies )
    END IF

    IF( TRIM( EquationOfState ) .NE. 'TABLE' )THEN
      CALL BuildOpacities_MF_Level( iLevel, BA, DM )
    END IF

    CALL CreateMesh_MF( iLevel, MeshX )
    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )
    CALL FillPatch( iLevel, MF_uGF, &
                    ApplyBoundaryConditions_Geometry_Option = .TRUE. )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL InitializeFields_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCR(iLevel), MF_uCF(iLevel) )

    CALL FillPatch( iLevel, MF_uGF, MF_uCF, &
                    ApplyBoundaryConditions_Euler_Option = .TRUE. )
    CALL FillPatch( iLevel, MF_uGF, MF_uCR, &
                    ApplyBoundaryConditions_TwoMoment_Option = .TRUE. )

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
           ( MF_uCR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uGF(iLevel) % SetVal( Zero )
    CALL MF_uCF(iLevel) % SetVal( Zero )
    CALL MF_uPF(iLevel) % SetVal( Zero )
    CALL MF_uAF(iLevel) % SetVal( Zero )
    CALL MF_uDF(iLevel) % SetVal( Zero )
    CALL MF_uCR(iLevel) % SetVal( Zero )
    CALL MF_uPR(iLevel) % SetVal( Zero )
    CALL amrex_multifab_build &
           ( MF_uGR(iLevel), BA, DM, &
             nDOFX * nGR * nSpecies, swX )
    CALL MF_uGR(iLevel) % SetVal( Zero )
    CALL amrex_multifab_build &
           ( MF_uAR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uAR(iLevel) % SetVal( Zero )
    CALL amrex_multifab_build &
           ( MF_Permute(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_Permute(iLevel) % SetVal( Zero )
 
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, &
               nDOFX_X1 * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies )
 
    IF( TRIM( EquationOfState ) .NE. 'TABLE' )THEN
      CALL BuildOpacities_MF_Level( iLevel, BA, DM )
    END IF
 
    CALL FillCoarsePatch( iLevel, MF_uGF, &
                          ApplyBoundaryConditions_Geometry_Option = .TRUE. )
 
    CALL CreateMesh_MF( iLevel, MeshX )
 
    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )
 
    CALL FillCoarsePatch( iLevel, MF_uDF )
 
    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCF, &
                          ApplyBoundaryConditions_Euler_Option = .TRUE. )
    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uCR, &
                    ApplyBoundaryConditions_TwoMoment_Option = .TRUE. ) 
 
    CALL ApplyPositivityLimiter_TwoMoment_MF ( MF_uGF, MF_uCF, MF_uCR )
 
    CALL FillCoarsePatch( iLevel, MF_uGF, MF_uPR )

    CALL DestroyMesh_MF( MeshX )
 
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
    CALL amrex_multifab_destroy( MF_uAR    (iLevel) )
    CALL amrex_multifab_destroy( MF_Permute(iLevel) )

    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_destroy( FluxRegister_TwoMoment(iLevel) )
    IF( TRIM( EquationOfState ) .NE. 'TABLE' )THEN
      CALL ClearOpacities_MF_Level( iLevel )
    END IF
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
    CALL MF_uGF_tmp % SetVal( Zero )
 
    CALL amrex_multifab_build( MF_uCF_tmp, BA, DM, nDOFX * nCF, swX )
    CALL MF_uCF_tmp % SetVal( Zero )
 
    CALL amrex_multifab_build( MF_uPF_tmp, BA, DM, nDOFX * nPF, swX )
    CALL MF_uPF_tmp % SetVal( Zero )
 
    CALL amrex_multifab_build( MF_uAF_tmp, BA, DM, nDOFX * nAF, swX )
    CALL MF_uAF_tmp % SetVal( Zero )
 
    CALL amrex_multifab_build( MF_uDF_tmp, BA, DM, nDOFX * nDF, swX )
    CALL MF_uDF_tmp % SetVal( Zero )
 
    CALL amrex_multifab_build &
           ( MF_uCR_tmp, BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uCR_tmp % SetVal( Zero )
 
    CALL amrex_multifab_build &
           ( MF_uPR_tmp, BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nPR * nSpecies, swX )
    CALL MF_uPR_tmp % SetVal( Zero )
 
    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, &
                    ApplyBoundaryConditions_Geometry_Option = .TRUE. )
 
    CALL FillPatch( iLevel, MF_uDF, MF_uDF_tmp )
 
    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, MF_uCF, MF_uCF_tmp, &
                    ApplyBoundaryConditions_Euler_Option = .TRUE. )
 
    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, MF_uCR, MF_uCR_tmp, &
                    ApplyBoundaryConditions_TwoMoment_Option = .TRUE. )
 
    CALL FillPatch( iLevel, MF_uGF, MF_uGF_tmp, MF_uPR, MF_uPR_tmp )
 
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
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uCR(iLevel) % SetVal( Zero )
 
    CALL amrex_multifab_build &
           ( MF_uPR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nPR * nSpecies, swX )
    CALL MF_uPR(iLevel) % SetVal( Zero )
 
    CALL amrex_multifab_build &
           ( MF_uAR(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uAR(iLevel) % SetVal( Zero )
 
    CALL amrex_multifab_build &
           ( MF_Permute(iLevel), BA, DM, &
             nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_Permute(iLevel) % SetVal( Zero )
 
    CALL amrex_multifab_build &
           ( MF_uGR(iLevel), BA, DM, &
             nDOFX * nGR * nSpecies, swX )
    CALL MF_uGR(iLevel) % SetVal( Zero )
 
    IF( iLevel .GT. 0 .AND. UseFluxCorrection_TwoMoment ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister_TwoMoment(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, &
               nDOFX_X1 * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies )
 
    IF( TRIM( EquationOfState ) .NE. 'TABLE' )THEN
      CALL BuildOpacities_MF_Level( iLevel, BA, DM )
    END IF
 
    CALL CreateMesh_MF( iLevel, MeshX )
 
    CALL MF_uGF(iLevel) % COPY( MF_uGF_tmp, 1, 1, nDOFX * nGF, swX )
    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )
 
    CALL MF_uCF(iLevel) % COPY( MF_uCF_tmp, 1, 1, nDOFX * nCF, swX )
    CALL MF_uDF(iLevel) % COPY( MF_uDF_tmp, 1, 1, nDOFX * nDF, swX )
    CALL MF_uCR(iLevel) % COPY( MF_uCR_tmp, 1, 1, &
        nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )
    CALL MF_uPR(iLevel) % COPY( MF_uPR_tmp, 1, 1, &
        nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nPR * nSpecies, swX )

 
    !CALL ApplyPositivityLimiter_TwoMoment_MF( MF_uGF, MF_uCF, MF_uCR )
    CALL ApplyPositivityLimiter_TwoMoment_MF &
           ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uCR(iLevel) )
 
    CALL amrex_multifab_destroy( MF_uPR_tmp )
    CALL amrex_multifab_destroy( MF_uCR_tmp )
    CALL amrex_multifab_destroy( MF_uDF_tmp )
    CALL amrex_multifab_destroy( MF_uAF_tmp )
    CALL amrex_multifab_destroy( MF_uPF_tmp )
    CALL amrex_multifab_destroy( MF_uCF_tmp )
    CALL amrex_multifab_destroy( MF_uGF_tmp )

    CALL DestroyMesh_MF( MeshX )
 
  END SUBROUTINE RemakeLevel

  SUBROUTINE ErrorEstimate( iLevel, cp, Time, SetTag, ClearTag ) BIND(c)

    USE TaggingModule, ONLY: &
      TagElements, TagElements_Density, TagElements_ShadowCasting, &
      TagElements_TransparentVortex_Spherical, TagElements_Exterior, TagElements_TVSD, TagElements_TransparentVortex_patch, &
      TagElements_TransparentVortex_MovingWedge, TagElements_TransparentVortex_Wedge
    USE amrex_parallel_module, ONLY: &
      amrex_parallel_ioprocessor

    INTEGER,                INTENT(in), VALUE :: iLevel
    TYPE(c_ptr),            INTENT(in), VALUE :: cp
    REAL(DP),               INTENT(in), VALUE :: Time
    CHARACTER(KIND=c_char), INTENT(in), VALUE :: SetTag, ClearTag

    TYPE(amrex_parmparse)   :: PP
    TYPE(amrex_tagboxarray) :: Tag
    TYPE(amrex_mfiter)      :: MFI
    TYPE(amrex_box)         :: BX
    REAL(DP),               CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    CHARACTER(KIND=c_char), CONTIGUOUS, POINTER :: TagArr(:,:,:,:)
    LOGICAL, SAVE :: PrintedScheme = .FALSE.
    INTEGER       :: iTC

    IF( .NOT. ALLOCATED( TagCriteria ) )THEN
      CALL amrex_parmparse_build( PP, "amr" )
        CALL PP % queryarr( "TagCriteria", TagCriteria )
      CALL amrex_parmparse_destroy( PP )
      IF( .NOT. ALLOCATED( TagCriteria ) )THEN
        ALLOCATE( TagCriteria(1) )
        TagCriteria = Zero
      END IF
    END IF

    iTC = MIN( iLevel + 1, SIZE( TagCriteria ) )

    CALL amrex_parmparse_build( PP, "amr" )
      CALL PP % get( "RefinementScheme", RefinementScheme )
    CALL amrex_parmparse_destroy( PP )

    IF( .NOT. PrintedScheme .AND. amrex_parallel_ioprocessor() )THEN
      WRITE(*,'(A)') 'ErrorEstimate: RefinementScheme = |' &
        // TRIM( RefinementScheme ) // '|'
      PrintedScheme = .TRUE.
    END IF

    Tag = cp

    CALL CreateMesh_MF( iLevel, MeshX )

    !$OMP PARALLEL PRIVATE( MFI, BX, uCR, TagArr )
    CALL amrex_mfiter_build( MFI, MF_uCR( iLevel ), Tiling = UseTiling )

    DO WHILE( MFI % next() )

      BX     =  MFI % TileBox()
      uCR    => MF_uCR( iLevel ) % DataPtr( MFI )
      TagArr => Tag              % DataPtr( MFI )

      IF( TRIM( RefinementScheme ) == "Density" )THEN

        CALL TagElements_Density &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) == "ShadowCasting" )THEN

        CALL TagElements_ShadowCasting &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) == "TransparentVortex_Spherical" )THEN

        CALL TagElements_TransparentVortex_Spherical &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )


      ELSE IF( TRIM( RefinementScheme ) == "TransparentVortex_patch" )THEN

        CALL TagElements_TransparentVortex_patch &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) == "TransparentVortex_MovingWedge" )THEN

        CALL TagElements_TransparentVortex_MovingWedge &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) == "TransparentVortex_Wedge" )THEN

        CALL TagElements_TransparentVortex_Wedge &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) == "TVSD" )THEN

        CALL TagElements_TVSD &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE IF( TRIM( RefinementScheme ) == "Exterior" )THEN

        CALL TagElements_Exterior &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      ELSE

        CALL TagElements &
               ( iLevel, BX % lo, BX % hi, LBOUND( uCR ), UBOUND( uCR ), &
                 uCR, TagCriteria(iTC), SetTag, ClearTag, &
                 LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      END IF

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate
 
 
END MODULE InitializationModule
