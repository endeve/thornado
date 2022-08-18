module ThornadoInitializationModule

  use KindModule, only: &
    DP, SqrtTiny, One, Zero
  use UnitsModule, only : &
    MeV
  use ProgramHeaderModule, only: &
    InitializeProgramHeader, &
    DescribeProgramHeader, &
    InitializeProgramHeaderX, &
    bcZ, nNodesE, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iX_B0, iX_E0, iX_B1, iX_E1
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
#ifdef MICROPHYSICS_WEAKLIB
  use EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE, &
    Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y
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
    DestroyGeometryFields, &
    uGF
  use GeometryFieldsModuleE, only: &
    CreateGeometryFieldsE, &
    DestroyGeometryFieldsE, &
    uGE
  use GeometryComputationModule, only: &
    ComputeGeometryX
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
  use TwoMoment_MeshRefinementModule, only : &
    InitializeMeshRefinement_TwoMoment, &
    FinalizeMeshRefinement_TwoMoment
#ifdef TWOMOMENT_ORDER_1
  use TwoMoment_PositivityLimiterModule, only: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment
#elif TWOMOMENT_ORDER_V
  use TwoMoment_SlopeLimiterModule_OrderV, only: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment
  use TwoMoment_PositivityLimiterModule_OrderV, only: &
    InitializePositivityLimiter_TwoMoment, &
    FinalizePositivityLimiter_TwoMoment
  use Euler_PositivityLimiterModule_NonRelativistic_TABLE, only: &
    InitializePositivityLimiter_Euler_NonRelativistic_TABLE, &
    FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
#endif

  implicit none
  private

  public :: InitThornado
  public :: FreeThornado
  public :: InitThornado_Patch
  public :: FreeThornado_Patch

contains

  subroutine InitThornado &
    ( nNodes, nDimsX, nE, swE, eL_MeV, eR_MeV, zoomE, bcE, &
      EquationOfStateTableName_Option, External_EOS, &
      Gamma_IDEAL_Option, &
      PositivityLimiter_Option, UpperBry1_Option, &
      SlopeLimiter_Option, &
      OpacityTableName_EmAb_Option, OpacityTableName_Iso_Option, &
      OpacityTableName_NES_Option, OpacityTableName_Pair_Option, &
      OpacityTableName_Brem_Option, Verbose_Option )

    integer,  intent(in) :: nNodes, nDimsX, nE, swE, bcE
    real(dp), intent(in) :: eL_MeV, eR_MeV, zoomE

    character(len=*), intent(in), optional :: EquationOfStateTableName_Option

#ifdef MICROPHYSICS_WEAKLIB
    type(EquationOfStateTableType), pointer, &
                      intent(in), optional :: External_EOS
#else
    integer,          intent(in), optional :: External_EOS
#endif

    real(dp),         intent(in), optional :: Gamma_IDEAL_Option
    logical,          intent(in), optional :: PositivityLimiter_Option
    logical,          intent(in), optional :: SlopeLimiter_Option
    real(dp),         intent(in), optional :: UpperBry1_Option
    character(len=*), intent(in), optional :: OpacityTableName_EmAb_Option
    character(len=*), intent(in), optional :: OpacityTableName_Iso_Option
    character(len=*), intent(in), optional :: OpacityTableName_NES_Option
    character(len=*), intent(in), optional :: OpacityTableName_Pair_Option
    character(len=*), intent(in), optional :: OpacityTableName_Brem_Option
    logical,          intent(in), optional :: Verbose_Option

    logical  :: PositivityLimiter, SlopeLimiter, EnergyLimiter, Verbose
    integer  :: nX(3), bcX(3)
    integer  :: i
    real(dp) :: eL, eR, UpperBry1

    IF( PRESENT(PositivityLimiter_Option) )THEN
      PositivityLimiter = PositivityLimiter_Option
    ELSE
      PositivityLimiter = .FALSE.
    END IF

    IF( PRESENT(SlopeLimiter_Option) )THEN
      SlopeLimiter = SlopeLimiter_Option
    ELSE
      SlopeLimiter = .FALSE.
    END IF

    EnergyLimiter = .TRUE.
!!    IF( PRESENT(EnergyLimiter_Option) )THEN
!!      EnergyLimiter = EnergyLimiter_Option
!!    ELSE
!!      EnergyLimiter = .FALSE.
!!    END IF

    IF( PRESENT(Verbose_Option) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( PRESENT(UpperBry1_Option) )THEN
      UpperBry1 = UpperBry1_Option
    ELSE
      UpperBry1 = 1.0d0 - EPSILON(1.0d0)
    END IF

    IF(Verbose)THEN
      WRITE(*,*)
#ifdef TWOMOMENT_ORDER_V
      WRITE(*,*) 'INFO: USE TWOMOMENT_ORDER_V'
#elif TWOMOMENT_ORDER_1
      WRITE(*,*) 'INFO: USE TWOMOMENT_ORDER_1'
#else
      WRITE(*,*) 'INFO: USE Default TWOMOMENT_ORDER'
#endif
      WRITE(*,*) '---------------------------------------'
    END IF

    ! --- Convert from MeV (expected) to thornado code units ---

    IF(Verbose)THEN
      WRITE(*,'(A12,ES12.3,A4)') 'Emin = ', eL_MeV, 'MeV'
      WRITE(*,'(A12,ES12.3,A4)') 'Emax = ', eR_MeV, 'MeV'
      WRITE(*,'(A12,ES12.3)')   'ZoomE = ', zoomE
    END IF

    eL = eL_MeV * MeV
    eR = eR_MeV * MeV

    nX = 1
    DO i = 1, nDimsX
      nX(i) = nX(i) + 1
    END DO


    ! bcX is for general thornado setting
    ! flash calling thornado is handled differently
    ! check TimeSteppingModule_Flash.F90 for details
    bcX = [ 0, 0, 0 ]

    call InitializeDevice

    call InitializeProgramHeader &
           ( ProgramName_Option = '', nNodes_Option = nNodes, &
             nX_Option = nX, bcX_Option = bcX, &
             nE_Option = nE, swE_Option = swE, bcE_Option = bcE, &
             eL_Option = eL, eR_Option = eR, zoomE_Option = zoomE )

#ifdef THORNADO_DEBUG
    call DescribeProgramHeader
#endif

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
           ( nE, swE, Verbose_Option = Verbose )

    call ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    ! --- Two-Moment Solver ---

    call InitializeClosure_TwoMoment &
           ( Verbose_Option = Verbose )

#ifdef TWOMOMENT_ORDER_V
    call InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = 1.75_DP, &
             UseSlopeLimiter_Option &
               = SlopeLimiter, &
             Verbose_Option &
               = Verbose )
#endif

#ifdef TWOMOMENT_ORDER_1
    call InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option = 0.0_DP + SqrtTiny, &
             Max_1_Option = UpperBry1, &
             Min_2_Option = 0.0_DP + SqrtTiny, &
             UsePositivityLimiter_Option = PositivityLimiter, &
             Verbose_Option = Verbose )
#elif TWOMOMENT_ORDER_V
    call InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option = 0.0_DP + SqrtTiny, &
             Max_1_Option = UpperBry1, &
             Min_2_Option = 0.0_DP + SqrtTiny, &
             UsePositivityLimiter_Option = PositivityLimiter, &
             UseEnergyLimiter_Option = EnergyLimiter, &
             Verbose_Option = Verbose )
#endif

    ! --- Nuclear Equation of State ---
#ifdef MICROPHYSICS_WEAKLIB
    call InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = EquationOfStateTableName_Option, &
             Verbose_Option = Verbose , &
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
             OpacityTableName_Brem_Option &
               = OpacityTableName_Brem_Option, &
             EquationOfStateTableName_Option &
               = EquationOfStateTableName_Option, &
             Verbose_Option = Verbose )

    ! --- For refinement and coarsening of DG data

    call InitializeMeshRefinement_TwoMoment

    ! --- For applying limiter on fluid field
#if defined TWOMOMENT_ORDER_V
#if defined MICROPHYSICS_WEAKLIB
    call InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option &
               = .TRUE., &
             Verbose_Option &
               = Verbose, &
             Min_1_Option &
               = ( One + 1.0d-3 * EPSILON( One ) ) * Min_D, &
             Min_2_Option &
               = ( One + 1.0d-3 * EPSILON( One ) ) * Min_T, &
             Min_3_Option &
               = ( One + 1.0d-3 * EPSILON( One ) ) * Min_Y, &
             Max_1_Option &
               = ( One - 1.0d-3 * EPSILON( One ) ) * Max_D, &
             Max_2_Option &
               = ( One - 1.0d-3 * EPSILON( One ) ) * Max_T, &
             Max_3_Option &
               = ( One - 1.0d-3 * EPSILON( One ) ) * Max_Y )
#endif
#endif

  end subroutine InitThornado


  subroutine FreeThornado(write_timers)

    logical, intent(in) :: write_timers

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

#ifdef MICROPHYSICS_WEAKLIB
    call FinalizeEquationOfState_TABLE
#else
    call FinalizeEquationOfState_IDEAL
#endif

    call FinalizeOpacities_TABLE

    call FinalizeMeshRefinement_TwoMoment

#ifdef TWOMOMENT_ORDER_V

#ifdef MICROPHYSICS_WEAKLIB
    call FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
#endif

    call FinalizeSlopeLimiter_TwoMoment
#endif

    call FinalizePositivityLimiter_TwoMoment

    call FinalizeDevice

  end subroutine FreeThornado


  subroutine InitThornado_Patch &
    ( nX, swX, xL, xR, nSpecies, CoordinateSystem_Option )

    use ProgramHeaderModule, only: nE, swE, nNodesX

    integer,  intent(in) :: nX(3), swX(3), nSpecies
    real(dp), intent(in) :: xL(3), xR(3)
    character(len=*), intent(in), optional :: CoordinateSystem_Option

    integer :: iDim, bcX(3)
    character(24) :: CoordinateSystem

    IF( PRESENT(CoordinateSystem_Option) )THEN

      IF( TRIM(CoordinateSystem_Option) == 'spherical' )THEN
        CoordinateSystem = 'SPHERICAL'
      ELSE IF( TRIM(CoordinateSystem_Option) == 'cartesian' )THEN
        CoordinateSystem = 'CARTESIAN'
      ELSE
        print*, '[InitThornado_Patch] Invalid Coordinate System: ', &
                 CoordinateSystem_Option
      END IF

    ELSE
      CoordinateSystem = 'CARTESIAN'
    END IF

    ! bcX is for general thornado setting
    ! flash calling thornado is handled differently
    ! check TimeSteppingModule_Flash.F90 for details
    bcX = [ 0, 0, 0 ]

    call InitializeProgramHeaderX &
           ( nX_Option = nX, swX_Option = swX, bcX_Option = bcX, &
             xL_Option = xL, xR_Option  = xR,  &
             reinitializeZ_Option = .TRUE. )

    DO iDim = 1, 3

      call CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), &
               swX(iDim), xL(iDim), xR(iDim) )

    END DO

    call CreateGeometryFields &
           ( nX, swX, CoordinateSystem_Option = CoordinateSystem, &
             Verbose_Option = .FALSE. )

    CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    call CreateFluidFields &
           ( nX, swX, Verbose_Option = .FALSE. )

    call CreateRadiationFields &
           ( nX, swX, nE, swE, nSpecies_Option = nSpecies, &
             Verbose_Option = .FALSE. )

  end subroutine InitThornado_Patch


  subroutine FreeThornado_Patch()

    integer :: iDim

    DO iDim = 1, 3

      call DestroyMesh( MeshX(iDim) )

    END DO

    call DestroyGeometryFields

    call DestroyFluidFields

    call DestroyRadiationFields

  end subroutine FreeThornado_Patch

end module ThornadoInitializationModule
