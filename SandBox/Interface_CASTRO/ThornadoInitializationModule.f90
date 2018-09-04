module ThornadoInitializationModule

  use KindModule, only: &
    DP
  use ProgramHeaderModule, only: &
    InitializeProgramHeader, &
    nNodesX, nNodesE, &
    iE_B0, iE_E0, iE_B1, iE_E1
  use QuadratureModule, only: &
    InitializeQuadratures
  use ReferenceElementModuleX, only: &
    InitializeReferenceElementX
  use ReferenceElementModuleE, only: &
    InitializeReferenceElementE
  use ReferenceElementModule, only: &
    InitializeReferenceElement
  use PolynomialBasisModule_Lagrange, only: &
    InitializePolynomialBasis_Lagrange
  use ReferenceElementModuleX_Lagrange, only: &
    InitializeReferenceElementX_Lagrange
  use ReferenceElementModuleE_Lagrange, only: &
    InitializeReferenceElementE_Lagrange
  use ReferenceElementModule_Lagrange, only: &
    InitializeReferenceElement_Lagrange
  use EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE
  use OpacityModule_TABLE, only: &
    InitializeOpacities_TABLE
  use MeshModule, only: &
    MeshX, MeshE, &
    CreateMesh, &
    DestroyMesh
  use GeometryFieldsModule, only: &
    CreateGeometryFields, &
    DestroyGeometryFields
  use FluidFieldsModule, only: &
    CreateFluidFields, &
    DestroyFluidFields
  use GeometryFieldsModuleE, only: &
    CreateGeometryFieldsE, &
    DestroyGeometryFieldsE, &
    uGE
  use GeometryComputationModuleE, only: &
    ComputeGeometryE
  use RadiationFieldsModule, only: &
    CreateRadiationFields, &
    DestroyRadiationFields
  use TwoMoment_ClosureModule, only: &
    InitializeClosure_TwoMoment
  use TwoMoment_PositivityLimiterModule, only: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment

  implicit none
  private

  public :: InitThornado
  public :: InitThornado_Patch
  public :: FreeThornado_Patch

contains

  subroutine InitThornado( nDimsX, nE, nSpeciesIn ) &
      bind(C, name = "InitThornado")

    use RadiationFieldsModule, only : nSpecies

    integer, INTENT(in) :: nDimsX, nE, nSpeciesIn

    integer :: nX(3), i

    nX = 1
    DO i = 1, nDimsX
      nX(i) = nX(i) + 1
    END DO

    call InitializeProgramHeader &
           ( ProgramName_Option = '', nNodes_Option = 2, &
             nX_Option = nX, nE_Option = nE )

    nSpecies = nSpeciesIn

    call InitializeQuadratures

    call InitializeReferenceElementX

    call InitializeReferenceElementE

    call InitializeReferenceElement

    call InitializePolynomialBasis_Lagrange

    call InitializeReferenceElementX_Lagrange

    call InitializeReferenceElementE_Lagrange

    call InitializeReferenceElement_Lagrange

    call InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = 'EquationOfStateTable.h5', &
             Verbose_Option = .true. )

    call InitializeOpacities_TABLE &
           ( OpacityTableName_Option &
               = 'OpacityTable.h5', &
             Verbose_Option = .true. )

  end subroutine InitThornado

  subroutine InitThornado_Patch( nX, swX, xL, xR, swE, eL_in, eR_in ) &
      bind(C, name = "InitThornado_Patch")

    use UnitsModule          , only : MeV
    use RadiationFieldsModule, only : nSpecies
    use ProgramHeaderModule  , only : nE, nNodesX, nNodesE

    integer,  INTENT(in) :: nX(3), swX(3)
    REAL(DP), INTENT(in) :: xL(3), xR(3)
    integer,  INTENT(in) :: swE
    REAL(DP), INTENT(in) :: eL_in, eR_in

    REAL(DP) :: eL, eR
    integer  :: iDim

    ! --- Convert from MeV (expected) to thornado code units ---

    eL = eL_in * MeV
    eR = eR_in * MeV

    call InitializeProgramHeader &
           ( ProgramName_Option = '', nNodes_Option = 2, &
             nX_Option = nX, swX_Option = swX, &
             xL_Option = xL, xR_Option  = xR,  &
             nE_Option = nE, swE_Option = swE, &
             eL_Option = eL, eR_Option  = eR )

    ! Note we always use 3 here even if calling from Castro 2D
    DO iDim = 1, 3

      call CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), &
               swX(iDim), xL(iDim), xR(iDim) )

    END DO

    call CreateMesh( MeshE, nE, nNodesE, swE, eL, eR )

    call CreateGeometryFields &
           ( nX, swX, CoordinateSystem_Option = 'CARTESIAN', &
             Verbose_Option = .FALSE. )

    call CreateFluidFields &
           ( nX, swX, Verbose_Option = .FALSE. )

    call CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = .FALSE. )

    call ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    call CreateRadiationFields &
           ( nX, swX, nE, swE, nSpecies_Option = nSpecies, &
             Verbose_Option = .FALSE. )

    ! --- Two-Moment Solver ---

    call InitializeClosure_TwoMoment &
           ( Verbose_Option = .FALSE. )

    call InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option = 0.0_DP, &
             Max_1_Option = 1.0_DP, &
             Min_2_Option = 0.0_DP, &
             UsePositivityLimiter_Option = .TRUE., &
             Verbose_Option = .FALSE. )

 !  call InitializePositivityLimiter_TwoMoment &
 !         ( Min_1_Option = - HUGE( 1.0_DP ), &
 !           Max_1_Option = + HUGE( 1.0_DP ), &
 !           Min_2_Option = - HUGE( 1.0_DP ), &
 !           UsePositivityLimiter_Option = .FALSE., &
 !           Verbose_Option = .FALSE. )

  end subroutine InitThornado_Patch

  subroutine FreeThornado_Patch () &
      bind(C, name = "FreeThornado_Patch")

    integer :: iDim

    DO iDim = 1, 3

      call DestroyMesh( MeshX(iDim) )

    END DO

    call DestroyMesh( MeshE )

    call DestroyGeometryFields

    call DestroyFluidFields

    call DestroyGeometryFieldsE

    call DestroyRadiationFields

    ! --- Two-Moment Solver ---

    call FinalizePositivityLimiter_TwoMoment

  end subroutine FreeThornado_Patch


end module ThornadoInitializationModule
