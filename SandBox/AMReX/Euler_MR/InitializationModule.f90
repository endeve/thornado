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
    amrex_ref_ratio
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap, &
    amrex_distromap_build
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
  USE amrex_bc_types_module, ONLY: &
    amrex_bc_foextrap, &
    amrex_bc_bogus

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
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
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    CoordinateSystem
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nPF, &
    nAF, &
    nDF
  USE Euler_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_Euler
  USE Euler_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_Euler
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, &
    Max_D, &
    Min_T, &
    Max_T, &
    Min_Y, &
    Max_Y

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_FieldsModule, ONLY: &
    CreateFields_MF, &
    MF_uGF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF, &
    FluxRegister
  USE MF_Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler_MF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE InputParsingModule, ONLY: &
    InitializeParameters, &
    nLevels, &
    swX, &
    StepNo, &
    dt, &
    t_old, &
    t_new, &
    xL, &
    UseTiling, &
    do_reflux, &
    MaxGridSizeX, &
    xL, &
    xR, &
    CoordSys, &
    EquationOfState, &
    Gamma_IDEAL, &
    EosTableName, &
    UseSlopeLimiter, &
    BetaTVD, &
    BetaTVB, &
    SlopeTolerance, &
    UseCharacteristicLimiting, &
    UseTroubledCellIndicator, &
    SlopeLimiterMethod, &
    LimiterThresholdParameter, &
    UseConservativeCorrection, &
    UsePositivityLimiter, &
    Min_1, &
    Min_2, &
    lo_bc, &
    hi_bc, &
    lo_bc_uCF, &
    hi_bc_uCF, &
    ProgramName
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF
  USE AverageDownModule, ONLY: &
    AverageDownTo
  USE Euler_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_Euler
  USE MF_Euler_TimersModule, ONLY: &
    InitializeTimers_AMReX_Euler, &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Initialize, &
    Timer_AMReX_Euler_InputOutput

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

CONTAINS


  SUBROUTINE InitializeProgram

    CALL amrex_init()

    CALL amrex_amrcore_init()

    CALL InitializeTimers_AMReX_Euler

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    CALL InitializeParameters

    IF     ( CoordSys .EQ. 0 )THEN

      CoordinateSystem = 'CARTESIAN'

    ELSE IF( CoordSys .EQ. 1 )THEN

      CoordinateSystem = 'CYLINDRICAL'

    ELSE IF( CoordSys .EQ. 2 )THEN

      CoordinateSystem = 'SPHERICAL'

    ELSE

      CALL DescribeError_Euler_MF &
             ( 02, Message_Option = 'Invalid CoordSys:', &
                   Int_Option = [ CoordSys ] )

    END IF

    IF( amrex_parallel_ioprocessor() ) CALL DescribeProgramHeaderX

    CALL CreateFields_MF

    ALLOCATE( lo_bc    (1:amrex_spacedim,1) )
    ALLOCATE( hi_bc    (1:amrex_spacedim,1) )
    ALLOCATE( lo_bc_uCF(1:amrex_spacedim,1:nDOFX*nCF) )
    ALLOCATE( hi_bc_uCF(1:amrex_spacedim,1:nDOFX*nCF) )

    lo_bc     = amrex_bc_bogus
    hi_bc     = amrex_bc_bogus
    lo_bc_uCF = amrex_bc_foextrap
    hi_bc_uCF = amrex_bc_foextrap

    CALL InitializePolynomialBasisX_Lagrange
    CALL InitializePolynomialBasisX_Legendre

    CALL InitializePolynomialBasis_Lagrange
    CALL InitializePolynomialBasis_Legendre

    CALL CreateMesh_MF( 0, MeshX )

    CALL InitializePolynomialBasisMapping &
           ( [Zero], MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    CALL DestroyMesh_MF( MeshX )

    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    CALL InitializeMeshRefinement_Euler

    IF( EquationOfState .EQ. 'TABLE' )THEN

        CALL InitializeEquationOfState &
               ( EquationOfState_Option = EquationOfState, &
                 EquationOfStateTableName_Option = EosTableName, &
                 Verbose_Option = amrex_parallel_ioprocessor() )

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

    CALL amrex_init_virtual_functions &
           ( MakeNewLevelFromScratch, &
             MakeNewLevelFromCoarse, &
             RemakeLevel, &
             ClearLevel, &
             ErrorEstimate )

    ALLOCATE( StepNo(0:nLevels-1) )
    ALLOCATE( dt    (0:nLevels-1) )
    ALLOCATE( t_old (0:nLevels-1) )
    ALLOCATE( t_new (0:nLevels-1) )

    StepNo = 0
    dt     = 0.0_DP
    t_old  = 0.0_DP
    t_new  = 0.0_DP

    CALL amrex_init_from_scratch( 0.0_DP )

    CALL InitializeFluid_SSPRK_MF &
           ( Verbose_Option = amrex_parallel_ioprocessor() )

    CALL ApplySlopeLimiter_Euler_MF &
           ( t_new, MF_uGF, MF_uCF, MF_uDF )

    CALL ApplyPositivityLimiter_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uDF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF, &
             MF_uGF_Option = MF_uGF, &
             MF_uCF_Option = MF_uCF, &
             MF_uPF_Option = MF_uPF, &
             MF_uAF_Option = MF_uAF, &
             MF_uDF_Option = MF_uDF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

  END SUBROUTINE InitializeProgram


  SUBROUTINE MakeNewLevelFromScratch( iLevel, Time, pBA, pDM ) BIND(c)

    USE MF_GeometryModule, ONLY: &
      ComputeGeometryX_MF

    USE MF_InitializationModule, ONLY: &
      InitializeFields_Euler_MF

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
    IF( iLevel .GT. 0 .AND. do_reflux ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nDOFX_X1*nCF )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF(iLevel) )

    CALL InitializeFields_Euler_MF( iLevel, MF_uGF(iLevel), MF_uCF(iLevel) )

    IF( iLevel .GT. 0 )THEN

      CALL AverageDownTo( iLevel-1, MF_uCF )
      CALL AverageDownTo( iLevel-1, MF_uGF )

    END IF

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE MakeNewLevelFromScratch


  SUBROUTINE MakeNewLevelFromCoarse( iLevel, Time, pBA, pDM ) BIND(c)

    USE FillPatchModule, ONLY: FillCoarsePatch

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

print*,'Hello and goodbye from MakeNewLevelFromCoarse'
stop 'InitializationModule.f90'

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uDF(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    CALL amrex_fluxregister_destroy( FluxRegister(iLevel) )


  END SUBROUTINE ClearLevel


  SUBROUTINE RemakeLevel( iLevel, Time, pBA, pDM ) BIND(c)

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

print*,'Hello and goodbye from RemakeLevel'
stop 'InitializationModule.f90'

  END SUBROUTINE RemakeLevel


  SUBROUTINE ErrorEstimate( iLevel, cp, Time, SetTag, ClearTag ) BIND(c)

    USE TaggingModule, ONLY: &
      TagElements_Advection1D, &
      TagElements_RiemannProblem1D, &
      TagElements_Advection2D, &
      TagElements_KelvinHelmholtz2D, &
      TagElements_Advection3D, &
      TagElements_uCF

    INTEGER,                INTENT(in), VALUE :: iLevel
    TYPE(c_ptr),            INTENT(in), VALUE :: cp
    REAL(DP),               INTENT(in), VALUE :: Time
    CHARACTER(KIND=c_char), INTENT(in), VALUE :: SetTag, ClearTag

    REAL(DP), ALLOCATABLE, SAVE :: TagCriteria(:)
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

      SELECT CASE( TRIM( ProgramName ) )

        CASE( 'Advection1D' )

          CALL TagElements_Advection1D &
                 ( iLevel, BX % lo, BX % hi, &
                   LBOUND( uCF ), UBOUND( uCF ), uCF, &
                   TagCriteria(iLevel+1), &
                   SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'RiemannProblem1D' )

          CALL TagElements_RiemannProblem1D &
                 ( iLevel, BX % lo, BX % hi, &
                   LBOUND( uCF ), UBOUND( uCF ), uCF, &
                   TagCriteria(iLevel+1), &
                   SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'Advection2D' )

          CALL TagElements_Advection2D &
                 ( iLevel, BX % lo, BX % hi, &
                   LBOUND( uCF ), UBOUND( uCF ), uCF, &
                   TagCriteria(iLevel+1), &
                   SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'KelvinHelmholtz2D' )

          CALL TagElements_KelvinHelmholtz2D &
                 ( iLevel, BX % lo, BX % hi, &
                   LBOUND( uCF ), UBOUND( uCF ), uCF, &
                   TagCriteria(iLevel+1), &
                   SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE( 'Advection3D' )

          CALL TagElements_Advection3D &
                 ( iLevel, BX % lo, BX % hi, &
                   LBOUND( uCF ), UBOUND( uCF ), uCF, &
                   TagCriteria(iLevel+1), &
                   SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

        CASE DEFAULT

          CALL TagElements_uCF &
                 ( iLevel, BX % lo, BX % hi, &
                   LBOUND( uCF ), UBOUND( uCF ), uCF, &
                   TagCriteria(iLevel+1), &
                   SetTag, ClearTag, &
                   LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

      END SELECT

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate


END MODULE InitializationModule
