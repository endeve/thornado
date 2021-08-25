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
    amrex_max_level
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

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    DescribeProgramHeaderX
  USE PolynomialBasisModuleX_Lagrange,   ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX
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
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE AMReX_BoundaryConditionsModule, ONLY: &
    lo_bc, &
    hi_bc
  USE MF_FieldsModule, ONLY: &
    CreateFields_MF, &
    MF_uGF_old, &
    MF_uGF_new, &
    MF_uCF_old, &
    MF_uCF_new, &
    MF_uPF_old, &
    MF_uPF_new, &
    MF_uAF_old, &
    MF_uAF_new, &
    MF_uDF_old, &
    MF_uDF_new, &
    FluxRegister
  USE MF_TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE InputParsingModule, ONLY: &
    InitializeParameters, &
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
    EosTableName
  USE InputOutputModuleAMReX, ONLY: &
    WriteFieldsAMReX_PlotFile
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF
  USE AverageDownModule, ONLY: &
    AverageDown
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

    CALL InitializePolynomialBasisX_Lagrange

    CALL InitializeReferenceElementX
    CALL InitializeReferenceElementX_Lagrange

    IF( EquationOfState .EQ. 'TABLE' )THEN

        CALL InitializeEquationOfState &
               ( EquationOfState_Option = EquationOfState, &
                 EquationOfStateTableName_Option &
                   = EosTableName )

    ELSE

      CALL InitializeEquationOfState &
               ( EquationOfState_Option = EquationOfState, &
                 Gamma_IDEAL_Option = Gamma_IDEAL, &
                 Verbose_Option = amrex_parallel_ioprocessor() )

    END IF

    CALL amrex_init_virtual_functions &
           ( MakeNewLevelFromScratch, &
             MakeNewLevelFromCoarse, &
             RemakeLevel, &
             ClearLevel, &
             ErrorEstimate )

    ALLOCATE( StepNo(0:amrex_max_level) )
    ALLOCATE( dt    (0:amrex_max_level) )
    ALLOCATE( t_old (0:amrex_max_level) )
    ALLOCATE( t_new (0:amrex_max_level) )

    StepNo = 0
    dt     = 0.0_DP
    t_old  = 0.0_DP
    t_new  = 0.0_DP

    CALL amrex_init_from_scratch( 0.0_DP )

    CALL AverageDown( MF_uGF_new, MF_uCF_new )

    CALL InitializeFluid_SSPRK_MF &
           ( Verbose_Option = amrex_parallel_ioprocessor() )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Initialize )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF_new, MF_uCF_new, MF_uPF_new, MF_uAF_new )

    CALL WriteFieldsAMReX_PlotFile &
           ( t_new(0), StepNo, MF_uGF_new, &
             MF_uGF_Option = MF_uGF_new, &
             MF_uCF_Option = MF_uCF_new, &
             MF_uPF_Option = MF_uPF_new, &
             MF_uAF_Option = MF_uAF_new, &
             MF_uDF_Option = MF_uDF_new )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

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

    INTEGER :: nCompCF

    BA = pBA
    DM = pDM

    t_new(iLevel) = Time
    t_old(iLevel) = Time - 1.0e200_DP

    CALL ClearLevel( iLevel )

    CALL amrex_multifab_build( MF_uGF_new(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uGF_old(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL MF_uGF_new(iLevel) % SetVal( Zero )
    CALL MF_uGF_old(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uCF_new(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uCF_old(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL MF_uCF_new(iLevel) % SetVal( Zero )
    CALL MF_uCF_old(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uPF_new(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uPF_old(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL MF_uPF_new(iLevel) % SetVal( Zero )
    CALL MF_uPF_old(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uAF_new(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uAF_old(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL MF_uAF_new(iLevel) % SetVal( Zero )
    CALL MF_uAF_old(iLevel) % SetVal( Zero )

    CALL amrex_multifab_build( MF_uDF_new(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL amrex_multifab_build( MF_uDF_old(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL MF_uDF_new(iLevel) % SetVal( Zero )
    CALL MF_uDF_old(iLevel) % SetVal( Zero )

    nCompCF = MF_uCF_new(iLevel) % nComp()

    IF( iLevel .GT. 0 .AND. do_reflux ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister(iLevel), BA, DM, &
               amrex_ref_ratio(iLevel-1), iLevel, nCompCF )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL ComputeGeometryX_MF( MF_uGF_new(iLevel) )
    CALL InitializeFields_MF( iLevel, MF_uGF_new(iLevel), MF_uCF_new(iLevel) )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE MakeNewLevelFromScratch


  SUBROUTINE MakeNewLevelFromCoarse( iLevel, Time, pBA, pDM ) BIND(c)

    USE FillPatchModule, ONLY: FillCoarsePatch_uCF, FillCoarsePatch_uGF

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM
    INTEGER :: nCompCF

print*,'Hello and goodbye from MakeNewLevelFromCoarse'
RETURN

    BA = pBA
    DM = pDM

    CALL ClearLevel( iLevel )

    t_new(iLevel) = Time
    t_old(iLevel) = Time - 1.0e200_DP

    CALL amrex_multifab_build( MF_uGF_new(iLevel), BA, DM, nDOFX * nGF, swX )
    CALL amrex_multifab_build( MF_uGF_old(iLevel), BA, DM, nDOFX * nGF, swX )

    CALL amrex_multifab_build( MF_uCF_new(iLevel), BA, DM, nDOFX * nCF, swX )
    CALL amrex_multifab_build( MF_uCF_old(iLevel), BA, DM, nDOFX * nCF, swX )

    CALL amrex_multifab_build( MF_uPF_new(iLevel), BA, DM, nDOFX * nPF, swX )
    CALL amrex_multifab_build( MF_uPF_old(iLevel), BA, DM, nDOFX * nPF, swX )

    CALL amrex_multifab_build( MF_uAF_new(iLevel), BA, DM, nDOFX * nAF, swX )
    CALL amrex_multifab_build( MF_uAF_old(iLevel), BA, DM, nDOFX * nAF, swX )

    CALL amrex_multifab_build( MF_uDF_new(iLevel), BA, DM, nDOFX * nDF, swX )
    CALL amrex_multifab_build( MF_uDF_old(iLevel), BA, DM, nDOFX * nDF, swX )

    nCompCF = MF_uCF_new(iLevel) % nComp()

    IF( iLevel .GT. 0 .AND. do_reflux ) &
      CALL amrex_fluxregister_build &
             ( FluxRegister(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
               iLevel, nCompCF )

    CALL FillCoarsePatch_uGF( iLevel, Time, MF_uGF_new(iLevel) )
    CALL FillCoarsePatch_uCF( iLevel, Time, MF_uCF_new(iLevel) )
    CALL FillCoarsePatch( iLevel, Time, MF_uPF_new(iLevel) )
    CALL FillCoarsePatch( iLevel, Time, MF_uAF_new(iLevel) )
    CALL FillCoarsePatch( iLevel, Time, MF_uDF_new(iLevel) )

  END SUBROUTINE MakeNewLevelFromCoarse


  SUBROUTINE ClearLevel( iLevel ) BIND(c)

    INTEGER, INTENT(in), VALUE :: iLevel

    CALL amrex_multifab_destroy( MF_uDF_new(iLevel) )
    CALL amrex_multifab_destroy( MF_uDF_old(iLevel) )

    CALL amrex_multifab_destroy( MF_uAF_new(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF_old(iLevel) )

    CALL amrex_multifab_destroy( MF_uPF_new(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF_old(iLevel) )

    CALL amrex_multifab_destroy( MF_uCF_new(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF_old(iLevel) )

    CALL amrex_multifab_destroy( MF_uGF_new(iLevel) )
    CALL amrex_multifab_destroy( MF_uGF_old(iLevel) )

    CALL amrex_fluxregister_destroy( FluxRegister(iLevel) )

  END SUBROUTINE ClearLevel


  SUBROUTINE RemakeLevel( iLevel, Time, pBA, pDM ) BIND(c)

    USE FillPatchModule, ONLY: FillPatch_uCF, FillPatch_uGF

    INTEGER,     INTENT(in), VALUE :: iLevel
    REAL(DP),    INTENT(in), VALUE :: Time
    TYPE(c_ptr), INTENT(in), VALUE :: pBA, pDM

    TYPE(amrex_boxarray)  :: BA
    TYPE(amrex_distromap) :: DM
    TYPE(amrex_multifab)  :: new_MF_uGF_new
    TYPE(amrex_multifab)  :: new_MF_uCF_new
    TYPE(amrex_multifab)  :: new_MF_uPF_new
    TYPE(amrex_multifab)  :: new_MF_uAF_new
    TYPE(amrex_multifab)  :: new_MF_uDF_new

    INTEGER :: nCompGF, nCompCF, nCompPF, nCompAF, nCompDF

print*,'Hello and goodbye from RemakeLevel'
RETURN

    BA = pBA
    DM = pDM

    CALL amrex_multifab_build( new_MF_uGF_new, BA, DM, nDOFX * nGF, swX )
    CALL FillPatch_uGF( iLevel, Time, new_MF_uGF_new )

    CALL amrex_multifab_build( new_MF_uCF_new, BA, DM, nDOFX * nCF, swX )
    CALL FillPatch_uCF( iLevel, Time, new_MF_uCF_new )

    CALL amrex_multifab_build( new_MF_uPF_new, BA, DM, nDOFX * nPF, swX )
    CALL FillPatch( iLevel, Time, new_MF_uPF_new )

    CALL amrex_multifab_build( new_MF_uAF_new, BA, DM, nDOFX * nAF, swX )
    CALL FillPatch( iLevel, Time, new_MF_uAF_new )

    CALL amrex_multifab_build( new_MF_uDF_new, BA, DM, nDOFX * nDF, swX )
    CALL FillPatch( iLevel, Time, new_MF_uDF_new )

    CALL ClearLevel( iLevel )

    t_new( iLevel ) = Time
    t_old( iLevel ) = Time - 1.0e200_DP

    nCompGF = new_MF_uGF_new % nComp()
    CALL amrex_multifab_build( MF_uGF_new( iLevel ), BA, DM, nCompGF, swX )
    CALL amrex_multifab_build( MF_uGF_old( iLevel ), BA, DM, nCompGF, swX )

    nCompCF = new_MF_uCF_new % nComp()
    CALL amrex_multifab_build( MF_uCF_new( iLevel ), BA, DM, nCompCF, swX )
    CALL amrex_multifab_build( MF_uCF_old( iLevel ), BA, DM, nCompCF, swX )

    nCompPF = new_MF_uPF_new % nComp()
    CALL amrex_multifab_build( MF_uPF_new( iLevel ), BA, DM, nCompPF, swX )
    CALL amrex_multifab_build( MF_uPF_old( iLevel ), BA, DM, nCompPF, swX )

    nCompAF = new_MF_uAF_new % nComp()
    CALL amrex_multifab_build( MF_uAF_new( iLevel ), BA, DM, nCompAF, swX )
    CALL amrex_multifab_build( MF_uAF_old( iLevel ), BA, DM, nCompAF, swX )

    nCompDF = new_MF_uDF_new % nComp()
    CALL amrex_multifab_build( MF_uDF_new( iLevel ), BA, DM, nCompDF, swX )
    CALL amrex_multifab_build( MF_uDF_old( iLevel ), BA, DM, nCompDF, swX )

    IF ( iLevel .GT. 0 .AND. do_reflux ) &
       CALL amrex_fluxregister_build &
              ( FluxRegister(iLevel), BA, DM, &
                amrex_ref_ratio(iLevel-1), iLevel, nCompCF )

    CALL MF_uGF_new( iLevel ) % Copy( new_MF_uGF_new, 1, 1, nCompGF, swX )
    CALL MF_uCF_new( iLevel ) % Copy( new_MF_uCF_new, 1, 1, nCompCF, swX )
    CALL MF_uPF_new( iLevel ) % Copy( new_MF_uPF_new, 1, 1, nCompPF, swX )
    CALL MF_uAF_new( iLevel ) % Copy( new_MF_uAF_new, 1, 1, nCompAF, swX )
    CALL MF_uDF_new( iLevel ) % Copy( new_MF_uDF_new, 1, 1, nCompDF, swX )

    CALL amrex_multifab_destroy( new_MF_uDF_new )
    CALL amrex_multifab_destroy( new_MF_uAF_new )
    CALL amrex_multifab_destroy( new_MF_uPF_new )
    CALL amrex_multifab_destroy( new_MF_uCF_new )
    CALL amrex_multifab_destroy( new_MF_uGF_new )

  END SUBROUTINE RemakeLevel


  SUBROUTINE ErrorEstimate( iLevel, cp, Time, SetTag, ClearTag ) BIND(c)

    USE TaggingModule, ONLY: TagElements_uCF

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
    CALL amrex_mfiter_build( MFI, MF_uCF_new( iLevel ), Tiling = UseTiling )

    DO WHILE( MFI % next() )

       BX = MFI % TileBox()

       uCF    => MF_uCF_new( iLevel ) % DataPtr( MFI )
       TagArr => Tag                  % DataPtr( MFI )

       ! TagCriteria(iLevel+1) because iLevel starts at 0 but
       ! TagCriteria starts with 1
       CALL TagElements_uCF &
              ( iLevel, BX % lo, BX % hi, LBOUND( uCF ), UBOUND( uCF ), uCF, &
                TagCriteria(iLevel+1), &
                SetTag, ClearTag, LBOUND( TagArr ), UBOUND( TagArr ), TagArr )

    END DO

    CALL amrex_mfiter_destroy( MFI )
    !$OMP END PARALLEL

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ErrorEstimate

END MODULE InitializationModule
