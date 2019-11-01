module ThornadoInitializationModule

#include "Flash.h"

  use KindModule, only: &
    DP, SqrtTiny
  use UnitsModule, only : &
    MeV
  use ProgramHeaderModule, only: &
    InitializeProgramHeader, &
    InitializeProgramHeaderX, &
    nNodesE, &
    iE_B0, iE_E0, iE_B1, iE_E1
  use DeviceModule, only: &
    InitializeDevice, &
    FinalizeDevice
  use TimersModule, only: &
    InitializeTimers, &
    FinalizeTimers
  use QuadratureModule, only: &
    InitializeQuadratures
  use ReferenceElementModuleX, only: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  use ReferenceElementModuleE, only: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  use ReferenceElementModule, only: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  use PolynomialBasisModuleX_Lagrange, only: &
    InitializePolynomialBasisX_Lagrange
  use PolynomialBasisModule_Lagrange, only: &
    InitializePolynomialBasis_Lagrange
  use ReferenceElementModuleX_Lagrange, only: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  use ReferenceElementModuleE_Lagrange, only: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  use ReferenceElementModule_Lagrange, only: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  use SubcellReconstructionModule, only: &
    InitializeSubcellReconstruction, &
    FinalizeSubcellReconstruction
#ifdef FLASH_EOS_WEAKLIB
  use EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  use wlEquationOfStateTableModule, only: &
    EquationOfStateTableType
#else
  use EquationOfStateModule_IDEAL, only: &
    InitializeEquationOfState_IDEAL, &
    FinalizeEquationOfState_IDEAL
#endif
  use OpacityModule_TABLE, only: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  use MeshModule, only: &
    MeshX, MeshE, &
    CreateMesh, &
    DestroyMesh
  use GeometryFieldsModule, only: &
    CreateGeometryFields, &
    DestroyGeometryFields
  use GeometryFieldsModuleE, only: &
    CreateGeometryFieldsE, &
    DestroyGeometryFieldsE, &
    uGE
  use GeometryComputationModuleE, only: &
    ComputeGeometryE
  use FluidFieldsModule, only: &
    CreateFluidFields, &
    DestroyFluidFields
  use RadiationFieldsModule, only: &
    CreateRadiationFields, &
    DestroyRadiationFields
  use TwoMoment_ClosureModule, only: &
    InitializeClosure_TwoMoment
  use TwoMoment_PositivityLimiterModule, only: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment
  use TwoMoment_MeshRefinementModule, only : &
    InitializeMeshRefinement_TwoMoment, &
    FinalizeMeshRefinement_TwoMoment
  use TwoMoment_DiscretizationModule_Collisions_Neutrinos, only : &
    InitializeNonlinearSolverTally, &
    FinalizeNonlinearSolverTally

  implicit none
  private

  public :: InitThornado
  public :: FreeThornado
  public :: InitThornado_Patch
  public :: FreeThornado_Patch

contains

  subroutine InitThornado &
    ( swE, eL_MeV, eR_MeV, zoomE, &
      EquationOfStateTableName_Option, External_EOS, &
      Gamma_IDEAL_Option, &
      OpacityTableName_EmAb_Option, OpacityTableName_Iso_Option, &
      OpacityTableName_NES_Option, OpacityTableName_Pair_Option )

    integer,  intent(in) :: swE
    real(dp), intent(in) :: eL_MeV, eR_MeV, zoomE

    character(len=*), intent(in), optional :: EquationOfStateTableName_Option
#ifdef FLASH_EOS_WEAKLIB
    type(EquationOfStateTableType), pointer, &
                      intent(in), optional :: External_EOS
#else
    integer,          intent(in), optional :: External_EOS
#endif

    real(dp),         intent(in), optional :: Gamma_IDEAL_Option
    character(len=*), intent(in), optional :: OpacityTableName_EmAb_Option
    character(len=*), intent(in), optional :: OpacityTableName_Iso_Option
    character(len=*), intent(in), optional :: OpacityTableName_NES_Option
    character(len=*), intent(in), optional :: OpacityTableName_Pair_Option

    integer  :: nDimsX, nE, nX(3), nNodes, i
    real(dp) :: eL, eR

    ! --- Convert from MeV (expected) to thornado code units ---

    nDimsX = NDIM
    nE = THORNADO_NE
    nNodes = THORNADO_NNODES

    eL = eL_MeV * MeV
    eR = eR_MeV * MeV

    nX = 1
    DO i = 1, nDimsX
      nX(i) = nX(i) + 1
    END DO

    call InitializeDevice

    call InitializeProgramHeader &
           ( ProgramName_Option = '', nNodes_Option = nNodes, &
             nX_Option = nX, nE_Option = nE, swE_Option = swE, &
             eL_Option = eL, eR_Option = eR, zoomE_Option = zoomE )

    call InitializeTimers

    call InitializeQuadratures

    call InitializeReferenceElementX

    call InitializeReferenceElementE

    call InitializeReferenceElement

    call InitializePolynomialBasisX_Lagrange

    call InitializePolynomialBasis_Lagrange

    call InitializeReferenceElementX_Lagrange

    call InitializeReferenceElementE_Lagrange

    call InitializeReferenceElement_Lagrange

    ! --- For Mapping Between FV and DG Representations ---

    call InitializeSubcellReconstruction

    ! --- Energy Grid ---

    call CreateMesh &
           ( MeshE, nE, nNodesE, swE, eL, eR, zoomOption = zoomE )

    call CreateGeometryFieldsE &
           ( nE, swE, Verbose_Option = .FALSE. )

    call ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    ! --- Two-Moment Solver ---

    call InitializeClosure_TwoMoment &
           ( Verbose_Option = .FALSE. )

    call InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option = 0.0_DP + SqrtTiny, &
             Max_1_Option = 1.0d0 - EPSILON(1.0d0), &
             Min_2_Option = 0.0_DP + SqrtTiny, &
             UsePositivityLimiter_Option = .FALSE., &
             Verbose_Option = .FALSE. )

    ! --- Nuclear Equation of State ---
#ifdef FLASH_EOS_WEAKLIB
    call InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = EquationOfStateTableName_Option, &
             Verbose_Option = .false. , &
             External_EOS = External_EOS )
#else
    call InitializeEquationOfState_IDEAL &
           ( Gamma_IDEAL_Option = Gamma_IDEAL_Option )
#endif
    
    ! --- Neutrino Opacities ---

    call InitializeOpacities_TABLE &
           ( OpacityTableName_EmAb_Option &
               = OpacityTableName_EmAb_Option, &
             OpacityTableName_Iso_Option &
               = OpacityTableName_Iso_Option, &
             OpacityTableName_NES_Option &
               = OpacityTableName_NES_Option, &
             OpacityTableName_Pair_Option &
               = OpacityTableName_Pair_Option, &
             Verbose_Option = .false. )

    ! --- For refinement and coarsening of DG data

    call InitializeMeshRefinement_TwoMoment

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uGE )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( uGE )
#endif

  end subroutine InitThornado

  subroutine FreeThornado(write_timers)

    logical, intent(in) :: write_timers

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: uGE )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( uGE )
#endif

    call DestroyMesh( MeshE )

    call DestroyGeometryFieldsE

    if ( write_timers ) call FinalizeTimers

    call FinalizeReferenceElementX

    call FinalizeReferenceElementX_Lagrange

    call FinalizeReferenceElementE

    call FinalizeReferenceElementE_Lagrange

    call FinalizeReferenceElement

    call FinalizeReferenceElement_Lagrange

    call FinalizeSubcellReconstruction

#ifdef FLASH_EOS_WEAKLIB
    call FinalizeEquationOfState_TABLE
#else
    call FinalizeEquationOfState_IDEAL
#endif

    call FinalizeOpacities_TABLE

    call FinalizeMeshRefinement_TwoMoment

    call FinalizePositivityLimiter_TwoMoment

    call FinalizeDevice

  end subroutine FreeThornado

  subroutine InitThornado_Patch( nX, swX, xL, xR )

    use ProgramHeaderModule, only: nE, swE, nNodesX

    integer,  intent(in) :: nX(3), swX(3)
    real(dp), intent(in) :: xL(3), xR(3)

    integer :: iDim, nSpecies

    call InitializeProgramHeaderX &
           ( nX_Option = nX, swX_Option = swX, &
             xL_Option = xL, xR_Option  = xR,  &
             reinitializeZ_Option = .TRUE. )

    DO iDim = 1, 3

      call CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), &
               swX(iDim), xL(iDim), xR(iDim) )

    END DO

    call CreateGeometryFields &
           ( nX, swX, CoordinateSystem_Option = 'CARTESIAN', &
             Verbose_Option = .FALSE. )

    call CreateFluidFields &
           ( nX, swX, Verbose_Option = .FALSE. )

    nSpecies = THORNADO_NSPECIES

    call CreateRadiationFields &
           ( nX, swX, nE, swE, nSpecies_Option = nSpecies, &
             Verbose_Option = .FALSE. )

    !call InitializeNonlinearSolverTally

  end subroutine InitThornado_Patch

  subroutine FreeThornado_Patch()

    integer :: iDim

    DO iDim = 1, 3

      call DestroyMesh( MeshX(iDim) )

    END DO

    call DestroyGeometryFields

    call DestroyFluidFields

    call DestroyRadiationFields

    !call FinalizeNonlinearSolverTally

  end subroutine FreeThornado_Patch


end module ThornadoInitializationModule
