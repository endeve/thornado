#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_NeutrinoMatterSolverModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, One, Two, Three, FourPi, SqrtTiny
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    AtomicMassUnit, &
    Centimeter, &
    Erg, &
    Gram, &
    MeV, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, &
    nDOFE, nDOFX
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Collisions_OuterLoop, &
    Timer_Collisions_InnerLoop, &
    Timer_Collisions_ComputeOpacity, &
    Timer_Collisions_ComputeRates, &
    Timer_Collisions_InitializeRHS, &
    Timer_Collisions_NeutrinoRHS, &
    Timer_Collisions_MatterRHS, &
    Timer_Collisions_SolveLS, &
    Timer_Collisions_UpdateFP, &
    Timer_Collisions_CheckOuter, &
    Timer_Collisions_CheckInner, &
    Timer_Opacity_D0, &
    Timer_Opacity_LimitD0, &
    Timer_Opacity_EC, &
    Timer_Opacity_ES, &
    Timer_Opacity_NES, &
    Timer_Opacity_Pair, &
    Timer_Opacity_Brem, &
    Timer_OpacityRate_NES, &
    Timer_OpacityRate_Pair, &
    Timer_OpacityRate_Brem
  USE ArrayUtilitiesModule, ONLY: &
    CreatePackIndex, &
    ArrayPack, &
    ArrayUnpack, &
    ArrayCopy
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    LinearLeastSquares, &
    LinearLeastSquares_LWORK
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    iNuE, &
    iNuE_Bar, &
    LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputePressure_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, &
    LimitEquilibriumDistributions_DG, &
    ComputeNeutrinoOpacities_EC, &
    ComputeNeutrinoOpacities_ES, &
    ComputeNeutrinoOpacities_NES, &
    ComputeNeutrinoOpacities_Pair, &
    ComputeNeutrinoOpacities_Brem, &
    ComputeNeutrinoOpacityRates_NES, &
    ComputeNeutrinoOpacityRates_Pair, &
    ComputeNeutrinoOpacityRates_Brem, &
    ComputeNeutrinoOpacityRates_LinearCorrections_NES, &
    ComputeNeutrinoOpacityRates_LinearCorrections_Pair, &
    ComputeNeutrinoOpacityRates_LinearCorrections_Brem
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeEddingtonTensorComponents_dd
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: InitializeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling_FP_Nested_AA

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_Y = One
  REAL(DP), PARAMETER :: Unit_E = Erg / Gram
  REAL(DP), PARAMETER :: Unit_T = MeV
  REAL(DP), PARAMETER :: Unit_V = SpeedOfLight

  INTEGER :: MoveLeft = 1

#if   defined( TWOMOMENT_ORDER_1 )
  INTEGER, PARAMETER :: iD  = 1
  INTEGER, PARAMETER :: iY  = 1
  INTEGER, PARAMETER :: iEf = 2
  INTEGER, PARAMETER :: iV1 = 1
  INTEGER, PARAMETER :: iV2 = 1
  INTEGER, PARAMETER :: iV3 = 1
  INTEGER, PARAMETER :: nMatterEquations = 2
#elif defined( TWOMOMENT_ORDER_V )
  INTEGER, PARAMETER :: iD  = 1
  INTEGER, PARAMETER :: iY  = 1
  INTEGER, PARAMETER :: iEf = 2
  INTEGER, PARAMETER :: iV1 = 3
  INTEGER, PARAMETER :: iV2 = 4
  INTEGER, PARAMETER :: iV3 = 5
  INTEGER, PARAMETER :: nMatterEquations = 5
#elif defined( TWOMOMENT_RELATIVISTIC )
  INTEGER, PARAMETER :: iD  = 1
  INTEGER, PARAMETER :: iY  = 2
  INTEGER, PARAMETER :: iEf = 3
  INTEGER, PARAMETER :: iV1 = 4
  INTEGER, PARAMETER :: iV2 = 5
  INTEGER, PARAMETER :: iV3 = 6
  INTEGER, PARAMETER :: nMatterEquations = 6
#endif

  INTEGER  :: nE_G, nX_G, nZ(4)
  INTEGER  :: n_FP_inner, n_FP_outer
  INTEGER  :: iE_B, iE_E

  REAL(DP), ALLOCATABLE :: E_N(:)  ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:) ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:) ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:) ! --- Integration Weights Scaled by (hc)^3
  REAL(DP), ALLOCATABLE :: W3_S(:) ! --- Integration Weights Scaled by (hc)^3
  REAL(DP), ALLOCATABLE :: FourPiEp2(:)

  ! --- Solver scratch arrays ---

  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: DnuNorm
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_Dnu, Dnu_old
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_Inu_d_1, Inu_u_1_old
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_Inu_d_2, Inu_u_2_old
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_Inu_d_3, Inu_u_3_old

  REAL(DP), DIMENSION(:), ALLOCATABLE :: Omega
  REAL(DP), DIMENSION(:), ALLOCATABLE :: E_old, Ef_old, C_Ef, S_Ef, G_Ef, U_Ef
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Y_old, C_Y, S_Y, G_Y, U_Y
  REAL(DP), DIMENSION(:), ALLOCATABLE :: T_old
  REAL(DP), DIMENSION(:), ALLOCATABLE :: V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: D_old, cD_old, C_D, S_D, G_D, U_D

  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: SqrtGm
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Dnu_0
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Sigma_Iso
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Phi_0_Iso
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Phi_1_Iso
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_EmAb, Eta_EmAb
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NES , Eta_NES
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Pair, Eta_Pair
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Brem, Eta_Brem
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__In__u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__In__u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__In__u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__Out_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__Out_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__Out_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Pro_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Pro_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Pro_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Ann_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Ann_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Ann_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_3

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_0, H_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_1, H_II_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_0, J_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_1, J_II_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: S_Sigma

  INTEGER :: LWORK_outer
  INTEGER :: LWORK_inner

  INTEGER,  DIMENSION(:)    , ALLOCATABLE :: INFO

  ! --- Solver Parameters to be initialized

  LOGICAL :: SolverParametersInitialized = .FALSE.

  LOGICAL  :: Include_NES
  LOGICAL  :: Include_Pair
  LOGICAL  :: Include_Brem
  LOGICAL  :: Include_LinCorr
  INTEGER  :: M_outer
  INTEGER  :: M_inner
  INTEGER  :: MaxIter_outer
  INTEGER  :: MaxIter_inner
  INTEGER  :: M_FP
  REAL(DP) :: Rtol_outer
  REAL(DP) :: Rtol_inner
  REAL(DP) :: wMatrRHS(nMatterEquations)
  REAL(DP) :: DnuMax
  LOGICAL  :: FreezeOpacities

  ! --- Temporary arrays for scatter/gather (packing)

  INTEGER,  DIMENSION(:)    , ALLOCATABLE, TARGET :: Error_T
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: D_T, T_T, Y_T, E_T
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: SqrtGm_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Dnu_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Inu_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Inu_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Inu_u_3_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Dnu_0_T
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Sigma_Iso_T
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Phi_0_Iso_T
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Phi_1_Iso_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_EmAb_T, Eta_EmAb_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NES_T , Eta_NES_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Pair_T, Eta_Pair_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Brem_T, Eta_Brem_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__In__u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__In__u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__In__u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__Out_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__Out_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NES__Out_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Pro_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Pro_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Pro_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Ann_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Ann_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Pair_Ann_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_3_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_0_T, H_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_1_T, H_II_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_0_T, J_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_1_T, J_II_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: S_Sigma_T

CONTAINS


  SUBROUTINE InitializeNeutrinoMatterSolver( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    INTEGER :: iE1, iE2, iN_E, iS, iN_X

    REAL(DP) :: TMP(1)
    REAL(DP), ALLOCATABLE :: AMAT(:,:), BVEC(:)

    iE_B = iZ_B(1)
    iE_E = iZ_E(1)
    nZ   = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    n_FP_outer = nMatterEquations
    n_FP_inner = nE_G * nCR * nSpecies

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )
    ALLOCATE( FourPiEp2(nE_G) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S      = W2_N / ( PlanckConstant * SpeedOfLight )**3
    W3_S      = W3_N / ( PlanckConstant * SpeedOfLight )**3
    FourPiEp2 = FourPi * E_N * E_N

    ALLOCATE( DnuNorm(nSpecies,nX_G) )

    ALLOCATE( C_Dnu    (nE_G,nSpecies,nX_G) )
    ALLOCATE( C_Inu_d_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( C_Inu_d_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( C_Inu_d_3(nE_G,nSpecies,nX_G) )

    ALLOCATE( Dnu_old    (nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_1_old(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_2_old(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_3_old(nE_G,nSpecies,nX_G) )

    ALLOCATE( Omega(nX_G) )

    ALLOCATE(    E_old(nX_G) )
    ALLOCATE(   Ef_old(nX_G) )
    ALLOCATE( C_Ef    (nX_G) )
    ALLOCATE( S_Ef    (nX_G) )
    ALLOCATE( G_Ef    (nX_G) )
    ALLOCATE( U_Ef    (nX_G) )

    ALLOCATE(   Y_old(nX_G) )
    ALLOCATE( C_Y    (nX_G) )
    ALLOCATE( S_Y    (nX_G) )
    ALLOCATE( G_Y    (nX_G) )
    ALLOCATE( U_Y    (nX_G) )

    ALLOCATE(   T_old(nX_G) )

    ALLOCATE( V_u_1_old(nX_G) )
    ALLOCATE( C_V_d_1(nX_G) )
    ALLOCATE( S_V_d_1(nX_G) )
    ALLOCATE( G_V_d_1(nX_G) )
    ALLOCATE( U_V_d_1(nX_G) )

    ALLOCATE( V_u_2_old(nX_G) )
    ALLOCATE( C_V_d_2(nX_G) )
    ALLOCATE( S_V_d_2(nX_G) )
    ALLOCATE( G_V_d_2(nX_G) )
    ALLOCATE( U_V_d_2(nX_G) )

    ALLOCATE( V_u_3_old(nX_G) )
    ALLOCATE( C_V_d_3(nX_G) )
    ALLOCATE( S_V_d_3(nX_G) )
    ALLOCATE( G_V_d_3(nX_G) )
    ALLOCATE( U_V_d_3(nX_G) )

    ALLOCATE(   D_old(nX_G) )
    ALLOCATE(  cD_old(nX_G) )
    ALLOCATE( C_D    (nX_G) )
    ALLOCATE( S_D    (nX_G) )
    ALLOCATE( G_D    (nX_G) )
    ALLOCATE( U_D    (nX_G) )

    ALLOCATE(         SqrtGm(              nX_G) )
    ALLOCATE(          Dnu_0(nE_G,nSpecies,nX_G) )
    ALLOCATE(      Sigma_Iso(nE_G,         nX_G) )
    ALLOCATE(      Phi_0_Iso(nE_G,         nX_G) )
    ALLOCATE(      Phi_1_Iso(nE_G,         nX_G) )
    ALLOCATE(       Chi_EmAb(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_EmAb(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Chi_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Eta_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Brem(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Brem(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__In__u_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__In__u_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__In__u_3(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__Out_u_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__Out_u_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__Out_u_3(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Pro_u_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Pro_u_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Pro_u_3(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Ann_u_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Ann_u_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Ann_u_3(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Pro_u_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Pro_u_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Pro_u_3(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Ann_u_1(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Ann_u_2(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Ann_u_3(nE_G,nSpecies,nX_G) )

    ALLOCATE(  H_I_0 (nE_G,nE_G,nX_G) )
    ALLOCATE(  H_I_1 (nE_G,nE_G,nX_G) )
    ALLOCATE(  H_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE(  H_II_1(nE_G,nE_G,nX_G) )
    ALLOCATE(  J_I_0 (nE_G,nE_G,nX_G) )
    ALLOCATE(  J_I_1 (nE_G,nE_G,nX_G) )
    ALLOCATE(  J_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE(  J_II_1(nE_G,nE_G,nX_G) )
    ALLOCATE( S_Sigma(nE_G,nE_G,nX_G) )

    ALLOCATE( Error_T(nX_G) )

    ALLOCATE( D_T(nX_G) )
    ALLOCATE( T_T(nX_G) )
    ALLOCATE( Y_T(nX_G) )
    ALLOCATE( E_T(nX_G) )

    ALLOCATE( SqrtGm_T(nX_G) )

    ALLOCATE(     Dnu_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_3_T(nE_G,nSpecies,nX_G) )

    ALLOCATE(             Dnu_0_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(      Sigma_Iso_T(nE_G,         nX_G) )
    ALLOCATE(      Phi_0_Iso_T(nE_G,         nX_G) )
    ALLOCATE(      Phi_1_Iso_T(nE_G,         nX_G) )
    ALLOCATE(       Chi_EmAb_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_EmAb_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Chi_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Eta_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Brem_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Brem_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__In__u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__In__u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__In__u_3_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__Out_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__Out_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_NES__Out_u_3_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Pro_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Pro_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Pro_u_3_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Ann_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Ann_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Pair_Ann_u_3_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Pro_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Pro_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Pro_u_3_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Ann_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Ann_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( L_Brem_Ann_u_3_T(nE_G,nSpecies,nX_G) )

    ALLOCATE(   H_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(   H_I_1_T(nE_G,nE_G,nX_G) )
    ALLOCATE(  H_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(  H_II_1_T(nE_G,nE_G,nX_G) )
    ALLOCATE(   J_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(   J_I_1_T(nE_G,nE_G,nX_G) )
    ALLOCATE(  J_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(  J_II_1_T(nE_G,nE_G,nX_G) )
    ALLOCATE( S_Sigma_T(nE_G,nE_G,nX_G) )

    ALLOCATE( INFO(nX_G) )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS ) &
    !$OMP MAP( alloc: INFO, &
    !$OMP             DnuNorm, &
    !$OMP             C_Dnu, &
    !$OMP             C_Inu_d_1, &
    !$OMP             C_Inu_d_2, &
    !$OMP             C_Inu_d_3, &
    !$OMP             Dnu_old, &
    !$OMP             Inu_u_1_old, &
    !$OMP             Inu_u_2_old, &
    !$OMP             Inu_u_3_old, &
    !$OMP             Omega, &
    !$OMP             E_old, Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, &
    !$OMP             D_old, cD_old, C_D, S_D, G_D, U_D, &
    !$OMP             Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$OMP             T_old, &
    !$OMP             V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$OMP             V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$OMP             V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$OMP             SqrtGm, &
    !$OMP             Dnu_0, Sigma_Iso, Phi_0_Iso, Phi_1_Iso, &
    !$OMP             Chi_EmAb, Eta_EmAb, &
    !$OMP             Chi_NES, Eta_NES, &
    !$OMP             Chi_Pair, Eta_Pair, &
    !$OMP             Chi_Brem, Eta_Brem, &
    !$OMP             L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$OMP             L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$OMP             L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$OMP             L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$OMP             L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$OMP             L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$OMP             H_I_0, H_II_0, J_I_0, J_II_0, &
    !$OMP             H_I_1, H_II_1, J_I_1, J_II_1, S_Sigma, &
    !$OMP             D_T, T_T, Y_T, E_T, Error_T, &
    !$OMP             SqrtGm_T, &
    !$OMP             Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$OMP             Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$OMP             Chi_EmAb_T, Eta_EmAb_T, &
    !$OMP             Chi_NES_T, Eta_NES_T, &
    !$OMP             Chi_Pair_T, Eta_Pair_T, &
    !$OMP             Chi_Brem_T, Eta_Brem_T, &
    !$OMP             L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$OMP             L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$OMP             L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$OMP             L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$OMP             L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$OMP             L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T, &
    !$OMP             H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$OMP             H_I_1_T, H_II_1_T, J_I_1_T, J_II_1_T, S_Sigma_T )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS ) &
    !$ACC CREATE( INFO, &
    !$ACC         DnuNorm, &
    !$ACC         C_Dnu, &
    !$ACC         C_Inu_d_1, &
    !$ACC         C_Inu_d_2, &
    !$ACC         C_Inu_d_3, &
    !$ACC         Dnu_old, &
    !$ACC         Inu_u_1_old, &
    !$ACC         Inu_u_2_old, &
    !$ACC         Inu_u_3_old, &
    !$ACC         Omega, &
    !$ACC         E_old, Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, &
    !$ACC         D_old, cD_old, C_D, S_D, G_D, U_D, &
    !$ACC         Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$ACC         T_old, &
    !$ACC         V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$ACC         V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$ACC         V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$ACC         SqrtGm, &
    !$ACC         Dnu_0, Sigma_Iso, Phi_0_Iso, Phi_1_Iso, &
    !$ACC         Chi_EmAb, Eta_EmAb, &
    !$ACC         Chi_NES, Eta_NES, &
    !$ACC         Chi_Pair, Eta_Pair, &
    !$ACC         Chi_Brem, Eta_Brem, &
    !$ACC         L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$ACC         L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$ACC         L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$ACC         L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$ACC         L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$ACC         L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$ACC         H_I_0, H_II_0, J_I_0, J_II_0, &
    !$ACC         H_I_1, H_II_1, J_I_1, J_II_1, S_Sigma, &
    !$ACC         D_T, T_T, Y_T, E_T, Error_T, &
    !$ACC         SqrtGm_T, &
    !$ACC         Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$ACC         Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$ACC         Chi_EmAb_T, Eta_EmAb_T, &
    !$ACC         Chi_NES_T, Eta_NES_T, &
    !$ACC         Chi_Pair_T, Eta_Pair_T, &
    !$ACC         Chi_Brem_T, Eta_Brem_T, &
    !$ACC         L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$ACC         L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$ACC         L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$ACC         L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$ACC         L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$ACC         L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T, &
    !$ACC         H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$ACC         H_I_1_T, H_II_1_T, J_I_1_T, J_II_1_T, S_Sigma_T )
#endif

    IF ( .NOT. Include_NES ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
      DO iE2  = 1, nE_G
      DO iE1  = 1, nE_G

        H_I_0   (iE1,iE2,iN_X) = Zero
        H_I_1   (iE1,iE2,iN_X) = Zero
        H_II_0  (iE1,iE2,iN_X) = Zero
        H_II_1  (iE1,iE2,iN_X) = Zero
        H_I_0_T (iE1,iE2,iN_X) = Zero
        H_I_1_T (iE1,iE2,iN_X) = Zero
        H_II_0_T(iE1,iE2,iN_X) = Zero
        H_II_1_T(iE1,iE2,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

    IF ( .NOT. Include_Pair ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
      DO iE2  = 1, nE_G
      DO iE1  = 1, nE_G

        J_I_0   (iE1,iE2,iN_X) = Zero
        J_I_1   (iE1,iE2,iN_X) = Zero
        J_II_0  (iE1,iE2,iN_X) = Zero
        J_II_1  (iE1,iE2,iN_X) = Zero
        J_I_0_T (iE1,iE2,iN_X) = Zero
        J_I_1_T (iE1,iE2,iN_X) = Zero
        J_II_0_T(iE1,iE2,iN_X) = Zero
        J_II_1_T(iE1,iE2,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

    IF ( .NOT. Include_Brem ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
      DO iE2  = 1, nE_G
      DO iE1  = 1, nE_G

        S_Sigma  (iE1,iE2,iN_X) = Zero
        S_Sigma_T(iE1,iE2,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

    IF ( .NOT. Include_LinCorr ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(2)
#endif
      DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        Phi_1_Iso  (iN_E,iN_X) = Zero
        Phi_1_Iso_T(iN_E,iN_X) = Zero

      END DO
      END DO

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        L_NES__In__u_1  (iN_E,iS,iN_X) = Zero
        L_NES__In__u_2  (iN_E,iS,iN_X) = Zero
        L_NES__In__u_3  (iN_E,iS,iN_X) = Zero
        L_NES__Out_u_1  (iN_E,iS,iN_X) = Zero
        L_NES__Out_u_2  (iN_E,iS,iN_X) = Zero
        L_NES__Out_u_3  (iN_E,iS,iN_X) = Zero
        L_Pair_Pro_u_1  (iN_E,iS,iN_X) = Zero
        L_Pair_Pro_u_2  (iN_E,iS,iN_X) = Zero
        L_Pair_Pro_u_3  (iN_E,iS,iN_X) = Zero
        L_Pair_Ann_u_1  (iN_E,iS,iN_X) = Zero
        L_Pair_Ann_u_2  (iN_E,iS,iN_X) = Zero
        L_Pair_Ann_u_3  (iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_1  (iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_2  (iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_3  (iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_1  (iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_2  (iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_3  (iN_E,iS,iN_X) = Zero

        L_NES__In__u_1_T(iN_E,iS,iN_X) = Zero
        L_NES__In__u_2_T(iN_E,iS,iN_X) = Zero
        L_NES__In__u_3_T(iN_E,iS,iN_X) = Zero
        L_NES__Out_u_1_T(iN_E,iS,iN_X) = Zero
        L_NES__Out_u_2_T(iN_E,iS,iN_X) = Zero
        L_NES__Out_u_3_T(iN_E,iS,iN_X) = Zero
        L_Pair_Pro_u_1_T(iN_E,iS,iN_X) = Zero
        L_Pair_Pro_u_2_T(iN_E,iS,iN_X) = Zero
        L_Pair_Pro_u_3_T(iN_E,iS,iN_X) = Zero
        L_Pair_Ann_u_1_T(iN_E,iS,iN_X) = Zero
        L_Pair_Ann_u_2_T(iN_E,iS,iN_X) = Zero
        L_Pair_Ann_u_3_T(iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_1_T(iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_2_T(iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_3_T(iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_1_T(iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_2_T(iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_3_T(iN_E,iS,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

    IF( M_outer > 3 )THEN

      ALLOCATE( AMAT(n_FP_outer,M_outer) )
      ALLOCATE( BVEC(n_FP_outer        ) )
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA &
      !$OMP MAP( alloc: AMAT, BVEC )
#elif defined(THORNADO_OACC  )
      !$ACC ENTER DATA &
      !$ACC CREATE( AMAT, BVEC )
#endif
      CALL LinearLeastSquares_LWORK &
             ( 'N', n_FP_outer, M_outer-1, 1, AMAT, n_FP_outer, &
               BVEC, n_FP_outer, TMP, LWORK_outer )
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: AMAT, BVEC )
#elif defined(THORNADO_OACC  )
      !$ACC EXIT DATA &
      !$ACC DELETE( AMAT, BVEC )
#endif
      DEALLOCATE( AMAT )
      DEALLOCATE( BVEC )

    ELSE

      LWORK_outer = 1

    END IF

    IF( M_inner > 3 )THEN

      ALLOCATE( AMAT(n_FP_inner,M_inner) )
      ALLOCATE( BVEC(n_FP_inner        ) )
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA &
      !$OMP MAP( alloc: AMAT, BVEC )
#elif defined(THORNADO_OACC  )
      !$ACC ENTER DATA &
      !$ACC CREATE( AMAT, BVEC )
#endif
      CALL LinearLeastSquares_LWORK &
             ( 'N', n_FP_inner, M_inner-1, 1, AMAT, n_FP_inner, &
               BVEC, n_FP_inner, TMP, LWORK_inner )
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: AMAT, BVEC )
#elif defined(THORNADO_OACC  )
      !$ACC EXIT DATA &
      !$ACC DELETE( AMAT, BVEC )
#endif
      DEALLOCATE( AMAT )
      DEALLOCATE( BVEC )

    ELSE

      LWORK_inner = 1

    END IF

  END SUBROUTINE InitializeNeutrinoMatterSolver


  SUBROUTINE InitializeNeutrinoMatterSolverParameters &
    ( M_outer_Option, M_inner_Option, MaxIter_outer_Option, &
      MaxIter_inner_Option, Rtol_inner_Option, Rtol_outer_Option, &
      Include_NES_Option, Include_Pair_Option, Include_Brem_Option, &
      Include_LinCorr_Option, wMatrRHS_Option, DnuMax_Option, &
      FreezeOpacities_Option, Verbose_Option )

    INTEGER , INTENT(in), OPTIONAL :: M_outer_Option
    INTEGER , INTENT(in), OPTIONAL :: M_inner_Option
    INTEGER , INTENT(in), OPTIONAL :: MaxIter_outer_Option
    INTEGER , INTENT(in), OPTIONAL :: MaxIter_inner_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rtol_inner_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rtol_outer_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_NES_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_Pair_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_Brem_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_LinCorr_Option
    REAL(DP), INTENT(in), OPTIONAL :: wMatrRHS_Option(nMatterEquations)
    REAL(DP), INTENT(in), OPTIONAL :: DnuMax_Option
    LOGICAL , INTENT(in), OPTIONAL :: FreezeOpacities_Option
    LOGICAL , INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( SolverParametersInitialized ) RETURN

    IF( PRESENT( M_outer_Option ) )THEN
      M_outer = M_outer_Option
    ELSE
      M_outer = 3
    END IF

    IF( PRESENT( M_inner_Option ) )THEN
      M_inner = M_inner_Option
    ELSE
      M_inner = 2
    END IF

    IF( PRESENT( MaxIter_outer_Option ) )THEN
      MaxIter_outer = MaxIter_outer_Option
    ELSE
      MaxIter_outer = 100
    END IF

    IF( PRESENT( MaxIter_inner_Option ) )THEN
      MaxIter_inner = MaxIter_inner_Option
    ELSE
      MaxIter_inner = 100
    END IF

    IF( PRESENT( Rtol_inner_Option ) )THEN
      Rtol_inner = Rtol_inner_Option
    ELSE
      Rtol_inner = 1.0d-08
    END IF

    IF( PRESENT( Rtol_outer_Option ) )THEN
      Rtol_outer = Rtol_outer_Option
    ELSE
      Rtol_outer = 1.0d-08
    END IF

    IF( PRESENT( Include_NES_Option ) )THEN
      Include_NES = Include_NES_Option
    ELSE
      Include_NES = .TRUE.
    END IF

    IF( PRESENT( Include_Pair_Option ) )THEN
      Include_Pair = Include_Pair_Option
    ELSE
      Include_Pair = .TRUE.
    END IF

    IF( PRESENT( Include_Brem_Option ) )THEN
      Include_Brem = Include_Brem_Option
    ELSE
      Include_Brem = .TRUE.
    END IF

    IF( PRESENT( Include_LinCorr_Option ) )THEN
      Include_LinCorr = Include_LinCorr_Option
    ELSE
      Include_LinCorr = .FALSE.
    END IF

    IF( PRESENT( wMatrRHS_Option ) )THEN
      wMatrRHS = wMatrRHS_Option
    ELSE
      wMatrRHS = One ! --- One = On, Zero = Off
    END IF

    IF( PRESENT( DnuMax_Option ) )THEN
      DnuMax = DnuMax_Option
    ELSE
      DnuMax = One - EPSILON( One )
    END IF

    IF( PRESENT( FreezeOpacities_Option ) )THEN
      FreezeOpacities = FreezeOpacities_Option
    ELSE
      FreezeOpacities = .FALSE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: Neutrino-Matter Solver Parameters:'
      WRITE(*,'(A)') '  ----------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A32,I6.6)')     '', 'M_outer: '        , M_outer
      WRITE(*,'(A4,A32,I6.6)')     '', 'MaxIter_outer: '  , MaxIter_outer
      WRITE(*,'(A4,A32,ES10.3E3)') '', 'Rtol_outer: '     , Rtol_outer
      WRITE(*,*)
      WRITE(*,'(A4,A32,I6.6)')     '', 'M_inner: '        , M_inner
      WRITE(*,'(A4,A32,I6.6)')     '', 'MaxIter_inner: '  , MaxIter_inner
      WRITE(*,'(A4,A32,ES10.3E3)') '', 'Rtol_inner: '     , Rtol_inner
      WRITE(*,*)
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NES: '    , Include_NES
      WRITE(*,'(A4,A32,L1)')       '', 'Include_Pair: '   , Include_Pair
      WRITE(*,'(A4,A32,L1)')       '', 'Include_Brem: '   , Include_Brem
      WRITE(*,'(A4,A32,L1)')       '', 'Include_LinCorr: ', Include_LinCorr
      WRITE(*,*)
!!$      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iY ): '  , INT(wMatrRHS(iY ))
!!$      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iEf): '  , INT(wMatrRHS(iEf))
!!$      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iV1): '  , INT(wMatrRHS(iV1))
!!$      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iV2): '  , INT(wMatrRHS(iV2))
!!$      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iV3): '  , INT(wMatrRHS(iV3))
!!$      WRITE(*,*)
      WRITE(*,'(A4,A32,ES10.3E3)') '', 'DnuMax: '         , DnuMax
      WRITE(*,'(A4,A32,L1)')       '', 'FreezeOpacities: ', FreezeOpacities
      WRITE(*,*)

    END IF

    M_FP = MAX( M_outer, M_inner )

    SolverParametersInitialized = .TRUE.

  END SUBROUTINE InitializeNeutrinoMatterSolverParameters


  SUBROUTINE FinalizeNeutrinoMatterSolver

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS, &
    !$OMP               INFO, &
    !$OMP               DnuNorm, &
    !$OMP               C_Dnu, &
    !$OMP               C_Inu_d_1, &
    !$OMP               C_Inu_d_2, &
    !$OMP               C_Inu_d_3, &
    !$OMP               Dnu_old, &
    !$OMP               Inu_u_1_old, &
    !$OMP               Inu_u_2_old, &
    !$OMP               Inu_u_3_old, &
    !$OMP               Omega, &
    !$OMP               E_old, Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, &
    !$OMP               D_old, cD_old, C_D, S_D, G_D, U_D, &
    !$OMP               Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$OMP               T_old, &
    !$OMP               V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$OMP               V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$OMP               V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$OMP               SqrtGm, &
    !$OMP               Dnu_0, Sigma_Iso, Phi_0_Iso, Phi_1_Iso, &
    !$OMP               Chi_EmAb, Eta_EmAb, &
    !$OMP               Chi_NES, Eta_NES, &
    !$OMP               Chi_Pair, Eta_Pair, &
    !$OMP               Chi_Brem, Eta_Brem, &
    !$OMP               L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$OMP               L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$OMP               L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$OMP               L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$OMP               L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$OMP               L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$OMP               H_I_0, H_II_0, J_I_0, J_II_0, &
    !$OMP               H_I_1, H_II_1, J_I_1, J_II_1, S_Sigma, &
    !$OMP               D_T, T_T, Y_T, E_T, Error_T, &
    !$OMP               SqrtGm_T, &
    !$OMP               Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$OMP               Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$OMP               Chi_EmAb_T, Eta_EmAb_T, &
    !$OMP               Chi_NES_T, Eta_NES_T, &
    !$OMP               Chi_Pair_T, Eta_Pair_T, &
    !$OMP               Chi_Brem_T, Eta_Brem_T, &
    !$OMP               L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$OMP               L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$OMP               L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$OMP               L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$OMP               L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$OMP               L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T, &
    !$OMP               H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$OMP               H_I_1_T, H_II_1_T, J_I_1_T, J_II_1_T, S_Sigma_T )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS, &
    !$ACC         INFO, &
    !$ACC         DnuNorm, &
    !$ACC         C_Dnu, &
    !$ACC         C_Inu_d_1, &
    !$ACC         C_Inu_d_2, &
    !$ACC         C_Inu_d_3, &
    !$ACC         Dnu_old, &
    !$ACC         Inu_u_1_old, &
    !$ACC         Inu_u_2_old, &
    !$ACC         Inu_u_3_old, &
    !$ACC         Omega, &
    !$ACC         E_old, Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, &
    !$ACC         D_old, cD_old, C_D, S_D, G_D, U_D, &
    !$ACC         Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$ACC         T_old, &
    !$ACC         V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$ACC         V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$ACC         V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$ACC         SqrtGm, &
    !$ACC         Dnu_0, Sigma_Iso, Phi_0_Iso, Phi_1_Iso, &
    !$ACC         Chi_EmAb, Eta_EmAb, &
    !$ACC         Chi_NES, Eta_NES, &
    !$ACC         Chi_Pair, Eta_Pair, &
    !$ACC         Chi_Brem, Eta_Brem, &
    !$ACC         L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$ACC         L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$ACC         L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$ACC         L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$ACC         L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$ACC         L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$ACC         H_I_0, H_II_0, J_I_0, J_II_0, &
    !$ACC         H_I_1, H_II_1, J_I_1, J_II_1, S_Sigma, &
    !$ACC         D_T, T_T, Y_T, E_T, Error_T, &
    !$ACC         SqrtGm_T, &
    !$ACC         Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$ACC         Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$ACC         Chi_EmAb_T, Eta_EmAb_T, &
    !$ACC         Chi_NES_T, Eta_NES_T, &
    !$ACC         Chi_Pair_T, Eta_Pair_T, &
    !$ACC         Chi_Brem_T, Eta_Brem_T, &
    !$ACC         L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$ACC         L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$ACC         L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$ACC         L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$ACC         L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$ACC         L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T, &
    !$ACC         H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$ACC         H_I_1_T, H_II_1_T, J_I_1_T, J_II_1_T, S_Sigma_T )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2 )
    DEALLOCATE( INFO )
    DEALLOCATE( DnuNorm )
    DEALLOCATE( C_Dnu )
    DEALLOCATE( C_Inu_d_1 )
    DEALLOCATE( C_Inu_d_2 )
    DEALLOCATE( C_Inu_d_3 )
    DEALLOCATE( Dnu_old )
    DEALLOCATE( Inu_u_1_old )
    DEALLOCATE( Inu_u_2_old )
    DEALLOCATE( Inu_u_3_old )
    DEALLOCATE( Omega )
    DEALLOCATE( E_old, Ef_old, C_Ef, S_Ef, G_Ef, U_Ef )
    DEALLOCATE( Y_old, C_Y, S_Y, G_Y, U_Y )
    DEALLOCATE( T_old )
    DEALLOCATE( V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1 )
    DEALLOCATE( V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2 )
    DEALLOCATE( V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3 )
    DEALLOCATE( cD_old, D_old, C_D, S_D, G_D, U_D )
    DEALLOCATE( SqrtGm )
    DEALLOCATE( Dnu_0, Sigma_Iso, Phi_0_Iso, Phi_1_Iso )
    DEALLOCATE( Chi_EmAb, Eta_EmAb )
    DEALLOCATE( Chi_NES, Eta_NES )
    DEALLOCATE( Chi_Pair, Eta_Pair )
    DEALLOCATE( Chi_Brem, Eta_Brem )
    DEALLOCATE( L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3 )
    DEALLOCATE( L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3 )
    DEALLOCATE( L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3 )
    DEALLOCATE( L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3 )
    DEALLOCATE( L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3 )
    DEALLOCATE( L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3 )
    DEALLOCATE( H_I_0, H_II_0, J_I_0, J_II_0 )
    DEALLOCATE( H_I_1, H_II_1, J_I_1, J_II_1, S_Sigma )
    DEALLOCATE( D_T, T_T, Y_T, E_T, SqrtGm_T, Error_T )
    DEALLOCATE( Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T )
    DEALLOCATE( Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T )
    DEALLOCATE( Chi_EmAb_T, Eta_EmAb_T )
    DEALLOCATE( Chi_NES_T, Eta_NES_T )
    DEALLOCATE( Chi_Pair_T, Eta_Pair_T )
    DEALLOCATE( Chi_Brem_T, Eta_Brem_T )
    DEALLOCATE( L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T )
    DEALLOCATE( L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T )
    DEALLOCATE( L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T )
    DEALLOCATE( L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T )
    DEALLOCATE( L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T )
    DEALLOCATE( L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T )
    DEALLOCATE( H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T )
    DEALLOCATE( H_I_1_T, H_II_1_T, J_I_1_T, J_II_1_T, S_Sigma_T )

  END SUBROUTINE FinalizeNeutrinoMatterSolver


  SUBROUTINE ComputePointsAndWeightsE( E, W2, W3 )

    REAL(DP), INTENT(out) :: E (1:nE_G)
    REAL(DP), INTENT(out) :: W2(1:nE_G)
    REAL(DP), INTENT(out) :: W3(1:nE_G)

    INTEGER :: iN_E, iE, iNodeE

    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      E (iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

      W2(iN_E) = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**2
      W3(iN_E) = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**3

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA &
    ( dt, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, D, T, Y, E, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, nIterations_Inner, nIterations_Outer )

    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: D, T, Y, E
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Alpha
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Beta_u_1, Beta_u_2, Beta_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    INTEGER,  DIMENSION(:),     INTENT(inout) :: nIterations_Inner
    INTEGER,  DIMENSION(:),     INTENT(inout) :: nIterations_Outer

    ! --- Local Variables ---

    INTEGER  :: k_outer, Mk_outer, nX_P_outer
    INTEGER  :: k_inner, Mk_inner, nX_P_inner

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: P
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: Error

    LOGICAL,  ALLOCATABLE, DIMENSION(:) :: ITERATE_outer, ITERATE_inner
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: PackIndex_outer, UnpackIndex_outer
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: PackIndex_inner, UnpackIndex_inner

    ! --- Least-squares scratch arrays ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: AMAT_outer, GVEC_outer, FVEC_outer
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: BVEC_outer, GVECm_outer, FVECm_outer
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: WORK_outer, TAU_outer, Alpha_outer

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: AMAT_inner, GVEC_inner, FVEC_inner
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: BVEC_inner, GVECm_inner, FVECm_inner
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: WORK_inner, TAU_inner, Alpha_inner

#if defined( TWOMOMENT_RELATIVISTIC )
    ALLOCATE( P(nX_G) )
#endif

    ALLOCATE( Error(nX_G) )

    ALLOCATE(     ITERATE_outer(                   nX_G) )
    ALLOCATE(   PackIndex_outer(                   nX_G) )
    ALLOCATE( UnpackIndex_outer(                   nX_G) )
    ALLOCATE(        AMAT_outer(n_FP_outer,M_outer,nX_G) )
    ALLOCATE(        GVEC_outer(n_FP_outer,M_outer,nX_G) )
    ALLOCATE(        FVEC_outer(n_FP_outer,M_outer,nX_G) )
    ALLOCATE(        BVEC_outer(n_FP_outer,        nX_G) )
    ALLOCATE(       GVECm_outer(n_FP_outer,        nX_G) )
    ALLOCATE(       FVECm_outer(n_FP_outer,        nX_G) )
    ALLOCATE(        WORK_outer(       LWORK_outer,nX_G) )
    ALLOCATE(         TAU_outer(n_FP_outer,        nX_G) )
    ALLOCATE(       Alpha_outer(           M_outer,nX_G) )

    ALLOCATE(     ITERATE_inner(                   nX_G) )
    ALLOCATE(   PackIndex_inner(                   nX_G) )
    ALLOCATE( UnpackIndex_inner(                   nX_G) )
    ALLOCATE(        AMAT_inner(n_FP_inner,M_inner,nX_G) )
    ALLOCATE(        GVEC_inner(n_FP_inner,M_inner,nX_G) )
    ALLOCATE(        FVEC_inner(n_FP_inner,M_inner,nX_G) )
    ALLOCATE(        BVEC_inner(n_FP_inner,        nX_G) )
    ALLOCATE(       GVECm_inner(n_FP_inner,        nX_G) )
    ALLOCATE(       FVECm_inner(n_FP_inner,        nX_G) )
    ALLOCATE(        WORK_inner(       LWORK_inner,nX_G) )
    ALLOCATE(         TAU_inner(n_FP_inner,        nX_G) )
    ALLOCATE(       Alpha_inner(           M_inner,nX_G) )

    ITERATE_outer = .TRUE.
    ITERATE_inner = .TRUE.

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: &
    !$OMP   ITERATE_outer, ITERATE_inner ) &
    !$OMP MAP( alloc: &
    !$OMP   PackIndex_outer, UnpackIndex_outer, &
    !$OMP   AMAT_outer, GVEC_outer, FVEC_outer, &
    !$OMP   BVEC_outer, GVECm_outer, FVECm_outer, &
    !$OMP   WORK_outer, TAU_outer, Alpha_outer, &
    !$OMP   PackIndex_inner, UnpackIndex_inner, &
    !$OMP   AMAT_inner, GVEC_inner, FVEC_inner, &
    !$OMP   BVEC_inner, GVECm_inner, FVECm_inner, &
    !$OMP   WORK_inner, TAU_inner, Alpha_inner )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( &
    !$ACC   ITERATE_outer, ITERATE_inner ) &
    !$ACC CREATE( &
    !$ACC   PackIndex_outer, UnpackIndex_outer, &
    !$ACC   AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC   BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC   WORK_outer, TAU_outer, Alpha_outer, &
    !$ACC   PackIndex_inner, UnpackIndex_inner, &
    !$ACC   AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC   BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC   WORK_inner, TAU_inner, Alpha_inner )
#endif

    SqrtGm = SQRT( Gm_dd_11 * Gm_dd_22 * Gm_dd_33 )

    ! --- Initial RHS ---

    CALL TimersStart( Timer_Collisions_InitializeRHS )

#if   defined( TWOMOMENT_ORDER_1 )

    CALL InitializeRHS_OrderOne &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

    CALL LimitNeutrinoDistribution_OrderV & ! --- Enforce Dnu < 1
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    CALL InitializeRHS_OrderV &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, T, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_RELATIVISTIC )

    CALL ComputePressure_TABLE & 
           ( D, T, Y, P )

    CALL InitializeRHS_Relativistic &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, P, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

#endif

    CALL TimersStop( Timer_Collisions_InitializeRHS )

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Collisions_ComputeOpacity )

    CALL ComputeOpacities_Packed( D, T, Y, SqrtGm )

    CALL TimersStop( Timer_Collisions_ComputeOpacity )

    ! --- Start Outer Loop ---

    CALL TimersStart( Timer_Collisions_OuterLoop )

    k_outer = 0
    DO WHILE( ANY( ITERATE_outer(:) ) .AND. k_outer < MaxIter_outer )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      Error = 0

      CALL ComputeDnuNorm( ITERATE_outer, Dnu )

      CALL CreatePackIndex &
             ( ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      IF ( k_outer > 1 .AND. .NOT. FreezeOpacities ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Collisions_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( D, T, Y, SqrtGm, ITERATE_outer, nX_P_outer, &
                 PackIndex_outer, UnpackIndex_outer )

        CALL TimersStop( Timer_Collisions_ComputeOpacity )

      END IF

      ! --- Start Inner Loop ---

      CALL TimersStart( Timer_Collisions_InnerLoop )

      k_inner = 0
      DO WHILE( ANY( ITERATE_inner(:) ) .AND. k_inner < MaxIter_inner )

        k_inner  = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_inner, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

        CALL TimersStart( Timer_Collisions_ComputeRates )

        CALL ComputeRates_Packed &
               ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, ITERATE_inner, nX_P_inner, &
                 PackIndex_inner, UnpackIndex_inner, nX_P_outer )

        CALL TimersStop( Timer_Collisions_ComputeRates )

        ! --- Right-Hand Side Vectors and Residuals (inner) ---

        CALL TimersStart( Timer_Collisions_NeutrinoRHS )

#if   defined( TWOMOMENT_ORDER_1 )

        CALL ComputeNeutrinoRHS_OrderOne &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, dt, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

        CALL ComputeNeutrinoRHS_OrderV &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, dt, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_RELATIVISTIC )

        CALL ComputeNeutrinoRHS_Relativistic &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, dt, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

#endif

        CALL TimersStop( Timer_Collisions_NeutrinoRHS )

        ! --- Anderson Acceleration (inner) ---

        CALL TimersStart( Timer_Collisions_SolveLS )

        CALL SolveLS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner, &
                 AMAT_inner, BVEC_inner, Alpha_inner, TAU_inner, &
                 LWORK_inner, WORK_inner )

        CALL TimersStop( Timer_Collisions_SolveLS )

        ! --- Update Residuals and Solution Vectors (inner) ---

        CALL TimersStart( Timer_Collisions_UpdateFP )

#if   defined( TWOMOMENT_ORDER_1 )

        CALL UpdateNeutrinoRHS_OrderOne &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

        CALL UpdateNeutrinoRHS_OrderV &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_RELATIVISTIC )

        CALL UpdateNeutrinoRHS_Relativistic &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

#endif

        ! --- Check Convergence (inner) ---

        CALL TimersStart( Timer_Collisions_CheckInner )

        CALL CheckConvergence_Inner &
               ( ITERATE_inner, n_FP_inner, k_inner, &
                 nIterations_Inner, FVECm_inner )

        CALL TimersStop( Timer_Collisions_CheckInner )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Collisions_UpdateFP )

      END DO ! --- Inner Loop ---

      CALL TimersStop( Timer_Collisions_InnerLoop )

      ! --- Right-Hand Side Vectors and Residuals (outer) ---

      CALL TimersStart( Timer_Collisions_MatterRHS )

#if   defined( TWOMOMENT_ORDER_1 )

      CALL ComputeMatterRHS_OrderOne &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

      CALL ComputeMatterRHS_OrderV &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_RELATIVISTIC )

      CALL ComputeMatterRHS_Relativistic &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, & 
               D, E, P, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

#endif

      CALL TimersStop( Timer_Collisions_MatterRHS )

      ! --- Anderson Acceleration (outer) ---

      CALL TimersStart( Timer_Collisions_SolveLS )

      CALL SolveLS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVECm_outer, GVECm_outer, FVEC_outer, GVEC_outer, &
               AMAT_outer, BVEC_outer, Alpha_outer, TAU_outer, &
               LWORK_outer, WORK_outer )

      CALL TimersStop( Timer_Collisions_SolveLS )

      ! --- Update Residuals and Solution Vectors (outer) ---

      CALL TimersStart( Timer_Collisions_UpdateFP )

#if   defined( TWOMOMENT_ORDER_1 )

      CALL UpdateMatterRHS_OrderOne &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Y, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

      CALL UpdateMatterRHS_OrderV &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Y, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_RELATIVISTIC )

      CALL UpdateMatterRHS_Relativistic &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               D, Y, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#endif

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer, Error )

#if   defined( TWOMOMENT_RELATIVISTIC )

      CALL ComputePressure_TABLE( D, T, Y, P )

#endif

      CALL CheckErrorFlag_FP &
             ( Error, k_outer, k_inner, &
               D, E, Y, T, V_u_1, V_u_2, V_u_3, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      ! --- Check Convergence (outer) ---

      CALL TimersStart( Timer_Collisions_CheckOuter )

      CALL CheckConvergence_Outer &
             ( ITERATE_outer, ITERATE_inner, n_FP_outer, k_outer, &
               nIterations_Outer, FVECm_outer )

      CALL TimersStop( Timer_Collisions_CheckOuter )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Collisions_UpdateFP )

    END DO ! --- Outer Loop ---

    CALL TimersStop( Timer_Collisions_OuterLoop )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: &
    !$OMP   ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
    !$OMP   AMAT_outer, GVEC_outer, FVEC_outer, &
    !$OMP   BVEC_outer, GVECm_outer, FVECm_outer, &
    !$OMP   WORK_outer, TAU_outer, Alpha_outer, &
    !$OMP   ITERATE_inner, PackIndex_inner, UnpackIndex_inner, &
    !$OMP   AMAT_inner, GVEC_inner, FVEC_inner, &
    !$OMP   BVEC_inner, GVECm_inner, FVECm_inner, &
    !$OMP   WORK_inner, TAU_inner, Alpha_inner )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( &
    !$ACC   ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
    !$ACC   AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC   BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC   WORK_outer, TAU_outer, Alpha_outer, &
    !$ACC   ITERATE_inner, PackIndex_inner, UnpackIndex_inner, &
    !$ACC   AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC   BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC   WORK_inner, TAU_inner, Alpha_inner )
#endif

  END SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA


  SUBROUTINE ComputeOpacities_Packed &
    ( D, T, Y, SqrtGm, MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET   :: D, T, Y, SqrtGm
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: MASK
    INTEGER,                INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                INTENT(in), OPTIONAL :: nX_P0

    INTEGER                             :: nX, nX0, iX, iE
    REAL(DP), DIMENSION(:)    , POINTER :: D_P, T_P, Y_P, SqrtGm_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Dnu_0_P
    REAL(DP), DIMENSION(:,:)  , POINTER :: Sigma_Iso_P
    REAL(DP), DIMENSION(:,:)  , POINTER :: Phi_0_Iso_P
    REAL(DP), DIMENSION(:,:)  , POINTER :: Phi_1_Iso_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_EmAb_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_0_P, H_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_1_P, H_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_0_P, J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_1_P, J_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: S_Sigma_P

    IF( PRESENT( nX_P ) )THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF( PRESENT( nX_P0 ) )THEN
      nX0 = nX_P0
    ELSE
      nX0 = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P => D_T(1:nX)
      T_P => T_T(1:nX)
      Y_P => Y_T(1:nX)

      SqrtGm_P => SqrtGm_T(1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, D, T, Y, SqrtGm, D_P, T_P, Y_P, SqrtGm_P )

      Dnu_0_P     => Dnu_0_T    (:,:,1:nX)
      Sigma_Iso_P => Sigma_Iso_T(  :,1:nX)
      Phi_0_Iso_P => Phi_0_Iso_T(  :,1:nX)
      Phi_1_Iso_P => Phi_1_Iso_T(  :,1:nX)
      Chi_EmAb_P  => Chi_EmAb_T (:,:,1:nX)

      H_I_0_P     => H_I_0_T    (:,:,1:nX)
      H_I_1_P     => H_I_1_T    (:,:,1:nX)
      H_II_0_P    => H_II_0_T   (:,:,1:nX)
      H_II_1_P    => H_II_1_T   (:,:,1:nX)

      J_I_0_P     => J_I_0_T    (:,:,1:nX)
      J_I_1_P     => J_I_1_T    (:,:,1:nX)
      J_II_0_P    => J_II_0_T   (:,:,1:nX)
      J_II_1_P    => J_II_1_T   (:,:,1:nX)

      S_Sigma_P   => S_Sigma_T  (:,:,1:nX)

    ELSE

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      SqrtGm_P => SqrtGm(:)

      Dnu_0_P     => Dnu_0    (:,:,:)
      Sigma_Iso_P => Sigma_Iso(  :,:)
      Phi_0_Iso_P => Phi_0_Iso(  :,:)
      Phi_1_Iso_P => Phi_1_Iso(  :,:)
      Chi_EmAb_P  => Chi_EmAb (:,:,:)

      H_I_0_P     => H_I_0    (:,:,:)
      H_I_1_P     => H_I_1    (:,:,:)
      H_II_0_P    => H_II_0   (:,:,:)
      H_II_1_P    => H_II_1   (:,:,:)

      J_I_0_P     => J_I_0    (:,:,:)
      J_I_1_P     => J_I_1    (:,:,:)
      J_II_0_P    => J_II_0   (:,:,:)
      J_II_1_P    => J_II_1   (:,:,:)

      S_Sigma_P   => S_Sigma  (:,:,:)

    END IF

    ! --- Equilibrium Distributions ---

    CALL TimersStart( Timer_Opacity_D0 )

    CALL ComputeEquilibriumDistributions &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Y_P, Dnu_0_P )

!!$    CALL ComputeEquilibriumDistributions_DG &
!!$           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Y_P, SqrtGm_P, Dnu_0_P )

    CALL TimersStop( Timer_Opacity_D0 )

    CALL TimersStart( Timer_Opacity_LimitD0 )

    CALL LimitEquilibriumDistributions_DG &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, Dnu_0_P )

    CALL TimersStop( Timer_Opacity_LimitD0 )

    ! --- EmAb ---

    CALL TimersStart( Timer_Opacity_EC )

    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Y_P, Chi_EmAb_P )

    CALL TimersStop( Timer_Opacity_EC )

    ! --- Isoenergetic scattering ---

    CALL TimersStart( Timer_Opacity_ES )

    CALL ComputeNeutrinoOpacities_ES &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, 1, Phi_0_Iso_P )

    IF( Include_LinCorr )THEN

      CALL ComputeNeutrinoOpacities_ES &
             ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, 2, Phi_1_Iso_P )

    END IF

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX = 1, nX
    DO iE = 1, nE_G

      Sigma_Iso_P(iE,iX) &
        = FourPiEp2(iE) * ( Phi_0_Iso_P(iE,iX) - Third * Phi_1_Iso_P(iE,iX) )

    END DO
    END DO

    CALL TimersStop( Timer_Opacity_ES )

    IF( Include_NES )THEN

      ! --- NES Scattering Functions ---

      CALL TimersStart( Timer_Opacity_NES )

      CALL ComputeNeutrinoOpacities_NES &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, 1, H_I_0_P, H_II_0_P )

      IF( Include_LinCorr )THEN

        CALL ComputeNeutrinoOpacities_NES &
               ( 1, nE_G, 1, nX, D_P, T_P, Y_P, 2, H_I_1_P, H_II_1_P )

      END IF

      CALL TimersStop( Timer_Opacity_NES )

    END IF

    IF( Include_Pair )THEN

      ! --- Pair Kernels ---

      CALL TimersStart( Timer_Opacity_Pair )

      CALL ComputeNeutrinoOpacities_Pair &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, 1, J_I_0_P, J_II_0_P )

      IF( Include_LinCorr )THEN

        CALL ComputeNeutrinoOpacities_Pair &
               ( 1, nE_G, 1, nX, D_P, T_P, Y_P, 2, J_I_1_P, J_II_1_P )

      END IF

      CALL TimersStop( Timer_Opacity_Pair )

    END IF

    IF( Include_Brem )THEN

      ! --- Brem Kernels ---

      CALL TimersStart( Timer_Opacity_Brem )

      CALL ComputeNeutrinoOpacities_Brem &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, S_Sigma_P )

      CALL TimersStop( Timer_Opacity_Brem )

    END IF

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, Dnu_0_P, Dnu_0 )
      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, Chi_EmAb_P, Chi_EmAb )
      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, Sigma_Iso_P, Sigma_Iso )

      IF ( nX < nX0 ) THEN

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, H_I_0_P, H_II_0_P, H_I_0, H_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, J_I_0_P, J_II_0_P, J_I_0, J_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, S_Sigma_P, S_Sigma )

        IF( Include_LinCorr )THEN

          CALL ArrayUnpack &
                 ( nX, MASK, PackIndex, H_I_1_P, H_II_1_P, H_I_1, H_II_1 )

          CALL ArrayUnpack &
                 ( nX, MASK, PackIndex, J_I_1_P, J_II_1_P, J_I_1, J_II_1 )

        END IF

      END IF

    END IF

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(:),     INTENT(in), TARGET   :: D
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Dnu
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_1
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_2
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_3
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: MASK
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:),     INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P0

    REAL(DP), DIMENSION(:)    , POINTER :: D_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Dnu_P, Dnu_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Inu_u_1_P, Inu_u_2_P, Inu_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_NES_P , Eta_NES_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_Pair_P, Eta_Pair_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_Brem_P, Eta_Brem_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NES__In__u_1_P, L_NES__Out_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NES__In__u_2_P, L_NES__Out_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NES__In__u_3_P, L_NES__Out_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Pair_Pro_u_1_P, L_Pair_Ann_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Pair_Pro_u_2_P, L_Pair_Ann_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Pair_Pro_u_3_P, L_Pair_Ann_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Brem_Pro_u_1_P, L_Brem_Ann_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Brem_Pro_u_2_P, L_Brem_Ann_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Brem_Pro_u_3_P, L_Brem_Ann_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_0_P, H_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_1_P, H_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_0_P, J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_1_P, J_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: S_Sigma_P

    INTEGER :: nX, nX0

    IF ( PRESENT( nX_P ) ) THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( PRESENT( nX_P0 ) ) THEN
      nX0 = nX_P0
    ELSE
      nX0 = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P => D_T(1:nX)

      CALL ArrayPack( nX, UnpackIndex, D, D_P )

      Dnu_P     => Dnu_T    (:,:,1:nX)
      Dnu_0_P   => Dnu_0_T  (:,:,1:nX)
      Inu_u_1_P => Inu_u_1_T(:,:,1:nX)
      Inu_u_2_P => Inu_u_2_T(:,:,1:nX)
      Inu_u_3_P => Inu_u_3_T(:,:,1:nX)

      CALL ArrayPack( nX, UnpackIndex, Dnu, Dnu_0, Dnu_P, Dnu_0_P )
      CALL ArrayPack( nX, UnpackIndex, &
                      Inu_u_1, Inu_u_2, Inu_u_3, Inu_u_1_P, Inu_u_2_P, Inu_u_3_P )

      Chi_NES_P        => Chi_NES_T       (:,:,1:nX)
      Eta_NES_P        => Eta_NES_T       (:,:,1:nX)
      Chi_Pair_P       => Chi_Pair_T      (:,:,1:nX)
      Eta_Pair_P       => Eta_Pair_T      (:,:,1:nX)
      Chi_Brem_P       => Chi_Brem_T      (:,:,1:nX)
      Eta_Brem_P       => Eta_Brem_T      (:,:,1:nX)

      L_NES__In__u_1_P => L_NES__In__u_1_T(:,:,1:nX)
      L_NES__In__u_2_P => L_NES__In__u_2_T(:,:,1:nX)
      L_NES__In__u_3_P => L_NES__In__u_3_T(:,:,1:nX)
      L_NES__Out_u_1_P => L_NES__Out_u_1_T(:,:,1:nX)
      L_NES__Out_u_2_P => L_NES__Out_u_2_T(:,:,1:nX)
      L_NES__Out_u_3_P => L_NES__Out_u_3_T(:,:,1:nX)
      L_Pair_Pro_u_1_P => L_Pair_Pro_u_1_T(:,:,1:nX)
      L_Pair_Pro_u_2_P => L_Pair_Pro_u_2_T(:,:,1:nX)
      L_Pair_Pro_u_3_P => L_Pair_Pro_u_3_T(:,:,1:nX)
      L_Pair_Ann_u_1_P => L_Pair_Ann_u_1_T(:,:,1:nX)
      L_Pair_Ann_u_2_P => L_Pair_Ann_u_2_T(:,:,1:nX)
      L_Pair_Ann_u_3_P => L_Pair_Ann_u_3_T(:,:,1:nX)
      L_Brem_Pro_u_1_P => L_Brem_Pro_u_1_T(:,:,1:nX)
      L_Brem_Pro_u_2_P => L_Brem_Pro_u_2_T(:,:,1:nX)
      L_Brem_Pro_u_3_P => L_Brem_Pro_u_3_T(:,:,1:nX)
      L_Brem_Ann_u_1_P => L_Brem_Ann_u_1_T(:,:,1:nX)
      L_Brem_Ann_u_2_P => L_Brem_Ann_u_2_T(:,:,1:nX)
      L_Brem_Ann_u_3_P => L_Brem_Ann_u_3_T(:,:,1:nX)

      H_I_0_P          => H_I_0_T         (:,:,1:nX)
      H_I_1_P          => H_I_1_T         (:,:,1:nX)
      H_II_0_P         => H_II_0_T        (:,:,1:nX)
      H_II_1_P         => H_II_1_T        (:,:,1:nX)

      J_I_0_P          => J_I_0_T         (:,:,1:nX)
      J_I_1_P          => J_I_1_T         (:,:,1:nX)
      J_II_0_P         => J_II_0_T        (:,:,1:nX)
      J_II_1_P         => J_II_1_T        (:,:,1:nX)

      S_Sigma_P        => S_Sigma_T       (:,:,1:nX)

      IF ( nX < nX0 .OR. FreezeOpacities ) THEN

        CALL ArrayPack &
               ( nX, UnpackIndex, H_I_0, H_II_0, H_I_0_P, H_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, J_I_0, J_II_0, J_I_0_P, J_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, S_Sigma, S_Sigma_P )

        IF( Include_LinCorr )THEN

          CALL ArrayPack &
                 ( nX, UnpackIndex, H_I_1, H_II_1, H_I_1_P, H_II_1_P )

          CALL ArrayPack &
                 ( nX, UnpackIndex, J_I_1, J_II_1, J_I_1_P, J_II_1_P )

        END IF

      END IF

    ELSE

      D_P => D(:)

      Dnu_P            => Dnu           (:,:,:)
      Dnu_0_P          => Dnu_0         (:,:,:)
      Inu_u_1_P        => Inu_u_1       (:,:,:)
      Inu_u_2_P        => Inu_u_2       (:,:,:)
      Inu_u_3_P        => Inu_u_3       (:,:,:)

      Chi_NES_P        => Chi_NES       (:,:,:)
      Eta_NES_P        => Eta_NES       (:,:,:)
      Chi_Pair_P       => Chi_Pair      (:,:,:)
      Eta_Pair_P       => Eta_Pair      (:,:,:)
      Chi_Brem_P       => Chi_Brem      (:,:,:)
      Eta_Brem_P       => Eta_Brem      (:,:,:)

      L_NES__In__u_1_P => L_NES__In__u_1(:,:,:)
      L_NES__In__u_2_P => L_NES__In__u_2(:,:,:)
      L_NES__In__u_3_P => L_NES__In__u_3(:,:,:)
      L_NES__Out_u_1_P => L_NES__Out_u_1(:,:,:)
      L_NES__Out_u_2_P => L_NES__Out_u_2(:,:,:)
      L_NES__Out_u_3_P => L_NES__Out_u_3(:,:,:)
      L_Pair_Pro_u_1_P => L_Pair_Pro_u_1(:,:,:)
      L_Pair_Pro_u_2_P => L_Pair_Pro_u_2(:,:,:)
      L_Pair_Pro_u_3_P => L_Pair_Pro_u_3(:,:,:)
      L_Pair_Ann_u_1_P => L_Pair_Ann_u_1(:,:,:)
      L_Pair_Ann_u_2_P => L_Pair_Ann_u_2(:,:,:)
      L_Pair_Ann_u_3_P => L_Pair_Ann_u_3(:,:,:)
      L_Brem_Pro_u_1_P => L_Brem_Pro_u_1(:,:,:)
      L_Brem_Pro_u_2_P => L_Brem_Pro_u_2(:,:,:)
      L_Brem_Pro_u_3_P => L_Brem_Pro_u_3(:,:,:)
      L_Brem_Ann_u_1_P => L_Brem_Ann_u_1(:,:,:)
      L_Brem_Ann_u_2_P => L_Brem_Ann_u_2(:,:,:)
      L_Brem_Ann_u_3_P => L_Brem_Ann_u_3(:,:,:)

      H_I_0_P          => H_I_0         (:,:,:)
      H_I_1_P          => H_I_1         (:,:,:)
      H_II_0_P         => H_II_0        (:,:,:)
      H_II_1_P         => H_II_1        (:,:,:)

      J_I_0_P          => J_I_0         (:,:,:)
      J_I_1_P          => J_I_1         (:,:,:)
      J_II_0_P         => J_II_0        (:,:,:)
      J_II_1_P         => J_II_1        (:,:,:)

      S_Sigma_P        => S_Sigma       (:,:,:)

    END IF

    ! --- NES Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NES )

    CALL ComputeNeutrinoOpacityRates_NES &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, H_I_0_P, H_II_0_P, &
             Eta_NES_P, Chi_NES_P )

    IF( Include_LinCorr )THEN

      ! --- Compute Linear Rate Corrections (NES) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NES &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, H_I_1_P, H_II_1_P, &
               L_NES__In__u_1_P, L_NES__In__u_2_P, L_NES__In__u_3_P, &
               L_NES__Out_u_1_P, L_NES__Out_u_2_P, L_NES__Out_u_3_P )

    END IF

    CALL TimersStop( Timer_OpacityRate_NES )

    ! --- Pair Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_Pair )

    CALL ComputeNeutrinoOpacityRates_Pair &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, J_I_0_P, J_II_0_P, &
             Eta_Pair_P, Chi_Pair_P )

    IF( Include_LinCorr )THEN

      ! --- Compute Linear Rate Corrections (Pair) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, J_I_1_P, J_II_1_P, &
               L_Pair_Pro_u_1_P, L_Pair_Pro_u_2_P, L_Pair_Pro_u_3_P, &
               L_Pair_Ann_u_1_P, L_Pair_Ann_u_2_P, L_Pair_Ann_u_3_P )

    END IF

    CALL TimersStop( Timer_OpacityRate_Pair )

    ! --- Brem Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_Brem )

    CALL ComputeNeutrinoOpacityRates_Brem &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, S_Sigma_P, &
             Eta_Brem_P, Chi_Brem_P )

    IF( Include_LinCorr )THEN

      ! --- Compute Linear Rate Corrections (Brem) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_Brem &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, S_Sigma_P, &
               L_Brem_Pro_u_1_P, L_Brem_Pro_u_2_P, L_Brem_Pro_u_3_P, &
               L_Brem_Ann_u_1_P, L_Brem_Ann_u_2_P, L_Brem_Ann_u_3_P )

    END IF

    CALL TimersStop( Timer_OpacityRate_Brem )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NES_P, Eta_NES_P, Chi_NES, Eta_NES )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_Pair_P, Eta_Pair_P, Chi_Pair, Eta_Pair )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_Brem_P, Eta_Brem_P, Chi_Brem, Eta_Brem )

      IF( Include_LinCorr )THEN

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_NES__In__u_1_P, L_NES__In__u_2_P, L_NES__In__u_3_P, &
                 L_NES__In__u_1  , L_NES__In__u_2  , L_NES__In__u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_NES__Out_u_1_P, L_NES__Out_u_2_P, L_NES__Out_u_3_P, &
                 L_NES__Out_u_1  , L_NES__Out_u_2  , L_NES__Out_u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_Pair_Pro_u_1_P, L_Pair_Pro_u_2_P, L_Pair_Pro_u_3_P, &
                 L_Pair_Pro_u_1  , L_Pair_Pro_u_2  , L_Pair_Pro_u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_Pair_Ann_u_1_P, L_Pair_Ann_u_2_P, L_Pair_Ann_u_3_P, &
                 L_Pair_Ann_u_1  , L_Pair_Ann_u_2  , L_Pair_Ann_u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_Brem_Pro_u_1_P, L_Brem_Pro_u_2_P, L_Brem_Pro_u_3_P, &
                 L_Brem_Pro_u_1  , L_Brem_Pro_u_2  , L_Brem_Pro_u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_Brem_Ann_u_1_P, L_Brem_Ann_u_2_P, L_Brem_Ann_u_3_P, &
                 L_Brem_Ann_u_1  , L_Brem_Ann_u_2  , L_Brem_Ann_u_3 )

      END IF

    END IF

  END SUBROUTINE ComputeRates_Packed


  SUBROUTINE UpdateTemperature_Packed &
    ( D, E, Y, T, MASK, nX_P, PackIndex, UnpackIndex, Error )

    REAL(DP), DIMENSION(:), INTENT(in)   , TARGET   :: D, E, Y
    REAL(DP), DIMENSION(:), INTENT(inout), TARGET   :: T
    LOGICAL,  DIMENSION(:), INTENT(in)   , OPTIONAL :: MASK
    INTEGER,                INTENT(in)   , OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in)   , OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,  DIMENSION(:), INTENT(inout), TARGET   :: Error

    REAL(DP), DIMENSION(:), POINTER :: D_P, E_P, Y_P, T_P
    INTEGER,  DIMENSION(:), POINTER :: Error_P

    INTEGER :: nX

    IF( PRESENT( nX_P ) )THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P => D_T(1:nX)
      T_P => T_T(1:nX)
      Y_P => Y_T(1:nX)
      E_P => E_T(1:nX)
      Error_P => Error_T(1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, D, Y, T, E, D_P, Y_P, T_P, E_P )
      CALL ArrayPack &
             ( nX, UnpackIndex, Error, Error_P )

    ELSE

      D_P => D(:)
      Y_P => Y(:)
      T_P => T(:)
      E_P => E(:)
      Error_P => Error(:)

    END IF

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D_P, E_P, Y_P, T_P, Error_Option = Error_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, T_P, T )
      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, Error_P, Error )

    END IF

  END SUBROUTINE UpdateTemperature_Packed


  SUBROUTINE LimitNeutrinoDistribution_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER :: iN_E, iN_X, iS

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      Dnu(iN_E,iS,iN_X) = MIN( DnuMax, Dnu(iN_E,iS,iN_X) )

    END DO
    END DO
    END DO

  END SUBROUTINE LimitNeutrinoDistribution_OrderV


  SUBROUTINE InitializeRHS_OrderOne &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Y, E
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: SUM_Y, SUM_Ef

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Y, SUM_Ef )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Y, SUM_Ef )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Y, SUM_Ef )
#endif
    DO iN_X = 1, nX_G

      ! --- Store Initial Matter State ---

      Y_old (iN_X) = Y(iN_X)
      Ef_old(iN_X) = E(iN_X)

      ! --- Scaling Factors ---

      S_Y (iN_X) = One / ( D(iN_X) * Y(iN_X) / AtomicMassUnit )
      S_Ef(iN_X) = One / ( D(iN_X) * E(iN_X) )

      ! --- Initial Guess for Matter State ---

      U_Y (iN_X) = One
      U_Ef(iN_X) = One

      SUM_Y  = Zero
      SUM_Ef = Zero

#if   defined( THORNADO_OMP_OL )
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP REDUCTION( + : SUM_Y, SUM_Ef )
#elif defined( THORNADO_OACC   )
      !$ACC LOOP VECTOR COLLAPSE(2) &
      !$ACC REDUCTION( + : SUM_Y, SUM_Ef )
#endif
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        ! --- Old States for Neutrino Number Density and Flux ---

        C_Dnu    (iN_E,iS,iN_X) = Dnu    (iN_E,iS,iN_X)
        C_Inu_d_1(iN_E,iS,iN_X) = Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)
        C_Inu_d_2(iN_E,iS,iN_X) = Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)
        C_Inu_d_3(iN_E,iS,iN_X) = Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

        IF ( iS <= iNuE_Bar ) THEN
        SUM_Y  = SUM_Y  + Dnu(iN_E,iS,iN_X) * W2_S(iN_E) * LeptonNumber(iS)
        END IF
        SUM_Ef = SUM_Ef + Dnu(iN_E,iS,iN_X) * W3_S(iN_E)

      END DO
      END DO

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_Y (iN_X) = U_Y (iN_X) + wMatrRHS(iY ) * SUM_Y  * S_Y (iN_X)
      C_Ef(iN_X) = U_Ef(iN_X) + wMatrRHS(iEf) * SUM_Ef * S_Ef(iN_X)

    END DO

  END SUBROUTINE InitializeRHS_OrderOne


  SUBROUTINE InitializeRHS_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Y, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: vDotV, vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3, Ef
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3
    REAL(DP) :: SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( vDotV, SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, Ef )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( vDotV, SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          V_d_1, V_d_2, V_d_3, Ef )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( vDotV, SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, Ef, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iN_X = 1, nX_G

      SqrtGm(iN_X) = SQRT( Gm_dd_11(iN_X) * Gm_dd_22(iN_X) * Gm_dd_33(iN_X) )

      V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
      V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
      V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

      ! --- Specific Fluid Energy ---

      vDotV =   V_u_1(iN_X) * V_d_1 &
              + V_u_2(iN_X) * V_d_2 &
              + V_u_3(iN_X) * V_d_3

      Ef = E(iN_X) + Half * vDotV

      ! --- Compute Omega (Richardson damping coeff.) based on velocity ---
      ! --- For first interation: must be consistent with initial guess ---

      Omega(iN_X) = One / ( One + SQRT( vDotV ) )

      ! --- Store Initial Matter State ---

      D_old (iN_X) = D(iN_X)
      Y_old (iN_X) = Y(iN_X)
      E_old (iN_X) = E(iN_X)
      Ef_old(iN_X) = Ef
      T_old (iN_X) = T(iN_X)

      V_u_1_old(iN_X) = V_u_1(iN_X)
      V_u_2_old(iN_X) = V_u_2(iN_X)
      V_u_3_old(iN_X) = V_u_3(iN_X)

      ! --- Scaling Factors ---

      S_Y    (iN_X) = One / ( D(iN_X) * Y (iN_X) / AtomicMassUnit )
      S_Ef   (iN_X) = One / ( D(iN_X) * Ef )
      S_V_d_1(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_2(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_3(iN_X) = One / ( D(iN_X) * SpeedOfLight )

      ! --- Initial Guess for Matter State ---

      U_Y    (iN_X) = One
      U_Ef   (iN_X) = One
      U_V_d_1(iN_X) = V_d_1 / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2 / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3 / SpeedOfLight

      SUM_Y  = Zero
      SUM_Ef = Zero
      SUM_V1 = Zero
      SUM_V2 = Zero
      SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
      !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
      !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
      !$OMP REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
      !$ACC LOOP VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
      !$ACC          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
      !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
      !$ACC REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#endif
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        ! --- Store Initial Neutrino State ---

        Dnu_old    (iN_E,iS,iN_X) = Dnu    (iN_E,iS,iN_X)
        Inu_u_1_old(iN_E,iS,iN_X) = Inu_u_1(iN_E,iS,iN_X)
        Inu_u_2_old(iN_E,iS,iN_X) = Inu_u_2(iN_E,iS,iN_X)
        Inu_u_3_old(iN_E,iS,iN_X) = Inu_u_3(iN_E,iS,iN_X)

        vDotInu =   V_u_1(iN_X) * Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
                  + V_u_2(iN_X) * Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
                  + V_u_3(iN_X) * Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

#if defined( TWOMOMENT_ORDER_V )
        CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu    (iN_E,iS,iN_X), Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )
#endif

        vDotK_d_1 &
          = ( V_u_1(iN_X) * k_dd_11 &
            + V_u_2(iN_X) * k_dd_12 &
            + V_u_3(iN_X) * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
        vDotK_d_2 &
          = ( V_u_1(iN_X) * k_dd_12 &
            + V_u_2(iN_X) * k_dd_22 &
            + V_u_3(iN_X) * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
        vDotK_d_3 &
          = ( V_u_1(iN_X) * k_dd_13 &
            + V_u_2(iN_X) * k_dd_23 &
            + V_u_3(iN_X) * k_dd_33 ) * Dnu(iN_E,iS,iN_X)

        ! --- Eulerian Neutrino Number Density ---

        Nnu = Dnu(iN_E,iS,iN_X) + vDotInu

        ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

        Enu = Dnu(iN_E,iS,iN_X) + Two * vDotInu

        ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

        Fnu_d_1 &
          = Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
              + V_d_1 * Dnu(iN_E,iS,iN_X) + vDotK_d_1

        Fnu_d_2 &
          = Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
              + V_d_2 * Dnu(iN_E,iS,iN_X) + vDotK_d_2

        Fnu_d_3 &
          = Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) &
              + V_d_3 * Dnu(iN_E,iS,iN_X) + vDotK_d_3

        ! --- Old States for Neutrino Number Density and Flux ---

        C_Dnu    (iN_E,iS,iN_X) &
          = Dnu    (iN_E,iS,iN_X) + vDotInu

        C_Inu_d_1(iN_E,iS,iN_X) &
          = Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) + vDotK_d_1

        C_Inu_d_2(iN_E,iS,iN_X) &
          = Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) + vDotK_d_2

        C_Inu_d_3(iN_E,iS,iN_X) &
          = Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) + vDotK_d_3

        IF ( iS <= iNuE_Bar ) THEN
        SUM_Y  = SUM_Y  + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
        END IF

        SUM_Ef = SUM_Ef + Enu     * W3_S(iN_E)
        SUM_V1 = SUM_V1 + Fnu_d_1 * W3_S(iN_E)
        SUM_V2 = SUM_V2 + Fnu_d_2 * W3_S(iN_E)
        SUM_V3 = SUM_V3 + Fnu_d_3 * W3_S(iN_E)

      END DO
      END DO

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_Y    (iN_X) = U_Y    (iN_X) + wMatrRHS(iY ) * SUM_Y  * S_Y    (iN_X)
      C_Ef   (iN_X) = U_Ef   (iN_X) + wMatrRHS(iEf) * SUM_Ef * S_Ef   (iN_X)
      C_V_d_1(iN_X) = U_V_d_1(iN_X) + wMatrRHS(iV1) * SUM_V1 * S_V_d_1(iN_X)
      C_V_d_2(iN_X) = U_V_d_2(iN_X) + wMatrRHS(iV2) * SUM_V2 * S_V_d_2(iN_X)
      C_V_d_3(iN_X) = U_V_d_3(iN_X) + wMatrRHS(iV3) * SUM_V3 * S_V_d_3(iN_X)

    END DO

  END SUBROUTINE InitializeRHS_OrderV


  SUBROUTINE InitializeRHS_Relativistic &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, P, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Y, E, P, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Alpha
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Beta_u_1, Beta_u_2, Beta_u_3

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotV, vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, W, vvDotK
    REAL(DP) :: DT, B_d_1, B_d_2, B_d_3, cD, SJ(3), Tau_f
    REAL(DP) :: Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3
    REAL(DP) :: Inu_d_1, Inu_d_2, Inu_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: SUM_Y, SUM_V1, SUM_V2, SUM_V3, SUM_Ef

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, vDotV, W, SJ, Tau_f, &
    !$OMP          SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, vDotV, W, SJ, Tau_f, &
    !$ACC          SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, vDotV, W, SJ, Tau_f, &
    !$OMP          SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          DT, Inu_d_1, Inu_d_2, Inu_d_3, VDotInu, &
    !$OMP          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 )
#endif
    DO iN_X = 1, nX_G

      B_d_1 = Gm_dd_11(iN_X) * Beta_u_1(iN_X)
      B_d_2 = Gm_dd_22(iN_X) * Beta_u_2(iN_X)
      B_d_3 = Gm_dd_33(iN_X) * Beta_u_3(iN_X)

      V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
      V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
      V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

      ! --- Lorentz Factor ---

      vDotV =   V_u_1(iN_X) * V_d_1 &
              + V_u_2(iN_X) * V_d_2 &
              + V_u_3(iN_X) * V_d_3

      W = 1.0_DP / SQRT( 1.0_DP - vDotV )
      
      ! --- Compute Omega (Richardson damping coeff.) based on velocity ---
      ! --- For first interation: must be consistent with initial guess ---

      Omega(iN_X) = One / ( W + SQRT( vDotV ) )
      
      ! --- Store Initial Matter State ---

      D_old (iN_X) = D(iN_X)
      Y_old (iN_X) = Y(iN_X)
      cD_old(iN_X) = W * D(iN_X)
      Ef_old(iN_X) = E(iN_X)

      ! --- Scaling Factors ---

      S_D    (iN_X) = One / ( D (iN_X) )
      S_Y    (iN_X) = One / ( Y (iN_X) )
      S_Ef   (iN_X) = One / ( E(iN_X) )
      S_V_d_1(iN_X) = One / ( SpeedOfLight )
      S_V_d_2(iN_X) = One / ( SpeedOfLight )
      S_V_d_3(iN_X) = One / ( SpeedOfLight )

      IF (MoveLeft .EQ. 1) THEN
        
        Ef_old(iN_X) = E(iN_X) &
                     + ( W - 1.0_DP ) / W

        S_Ef(iN_X)   = 1.0_DP / Ef_old(iN_X)

      END IF

      ! --- Initial Guess for Matter State ---

      U_D    (iN_X) = One
      U_Y    (iN_X) = One
      U_Ef   (iN_X) = One
      U_V_d_1(iN_X) = V_d_1 / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2 / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3 / SpeedOfLight

      SUM_Y  = Zero
      SUM_Ef = Zero
      SUM_V1 = Zero
      SUM_V2 = Zero
      SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( DT, Inu_d_1, Inu_d_2, Inu_d_3, VDotInu, &
    !$OMP          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$OMP REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( DT, Inu_d_1, Inu_d_2, Inu_d_3, VDotInu, &
    !$ACC          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$ACC REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#endif
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      DT = 1.0_DP / ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) )

      Inu_d_1 = DT * ( B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) & 
              * Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
              - DT * ( B_d_1 * V_u_2(iN_X) *Gm_dd_22(iN_X) ) * Inu_u_2(iN_E,iS,iN_X)   &
              - DT * ( B_d_1 * V_u_3(iN_X) * Gm_dd_33(iN_X) ) * Inu_u_3(iN_E,iS,iN_X) 

      Inu_d_2 = DT * ( B_d_1 * V_u_1(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) &
              * Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
              - DT * ( B_d_2 * V_u_1(iN_X) * Gm_dd_11(iN_X) ) * Inu_u_1(iN_E,iS,iN_X) &
              - DT * ( Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) * B_d_2 * V_u_3(iN_X) ) 

      Inu_d_3 = DT * ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) - Alpha(iN_X) ) &
                            * Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                            - DT * ( Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) * B_d_3 * V_u_1(iN_X) ) &
                            - DT * ( Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) * B_d_3 * V_u_2(iN_X) )

      vDotInu = V_u_1(iN_X) * Inu_d_1 &
              + V_u_2(iN_X) * Inu_d_2 &
              + V_u_3(iN_X) * Inu_d_3

#if defined( TWOMOMENT_RELATIVISTIC )
       CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu(iN_E,iS,iN_X),     Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 Alpha(iN_X), Beta_u_1(iN_X), Beta_u_2(iN_X), Beta_u_3(iN_X), & 
                 V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), k_dd  )
#endif

      vDotK_d_1 &
        = ( V_u_1(iN_X) * k_dd(1,1) &
          + V_u_2(iN_X) * k_dd(2,1) &
          + V_u_3(iN_X) * k_dd(3,1) ) * Dnu(iN_E,iS,iN_X)
      vDotK_d_2 &                                            
        = ( V_u_1(iN_X) * k_dd(1,2) &
          + V_u_2(iN_X) * k_dd(2,2) &
          + V_u_3(iN_X) * k_dd(3,2) ) * Dnu(iN_E,iS,iN_X)
      vDotK_d_3 &                                            
        = ( V_u_1(iN_X) * k_dd(1,3) &
          + V_u_2(iN_X) * k_dd(2,3) &
          + V_u_3(iN_X) * k_dd(3,3) ) * Dnu(iN_E,iS,iN_X)

       vvDotK = V_u_1(iN_X) * vDotK_d_1 &
              + V_u_2(iN_X) * vDotK_d_2 & 
              + V_u_3(iN_X) * vDotK_d_3 

      ! --- Eulerian Neutrino Number Density ---

      Nnu = W * Dnu(iN_E,iS,iN_X) + vDotInu

      ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

      Enu = W**2 * Dnu(iN_E,iS,iN_X) +  W * Two * vDotInu + vvDotK 

      ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

      Fnu_d_1 &
        = W * Inu_d_1 + W**2 * V_d_1 * Dnu(iN_E,iS,iN_X) &
        + W * vDotInu * V_d_1  + vDotK_d_1

      Fnu_d_2 &
        = W * Inu_d_2 + W**2 * V_d_2 * Dnu(iN_E,iS,iN_X) &
        + W * vDotInu * V_d_2   + vDotK_d_2

      Fnu_d_3 &
        = W * Inu_d_3 + W**2 * V_d_3 * Dnu(iN_E,iS,iN_X) &
        + W * vDotInu * V_d_3  + vDotK_d_3

      ! --- Old States for Neutrino Number Density and Flux ---

      C_Dnu    (iN_E,iS,iN_X) = W * Dnu(iN_E,iS,iN_X) + vDotInu
      C_Inu_d_1(iN_E,iS,iN_X) = W * Inu_d_1 + vDotK_d_1
      C_Inu_d_2(iN_E,iS,iN_X) = W * Inu_d_2 + vDotK_d_2
      C_Inu_d_3(iN_E,iS,iN_X) = W * Inu_d_3 + vDotK_d_3


      IF ( iS <= iNuE_Bar ) THEN
        SUM_Y  = SUM_Y  + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
      END IF

      SUM_Ef  = SUM_Ef  + Enu     * W3_S(iN_E)
      SUM_V1  = SUM_V1  + Fnu_d_1 * W3_S(iN_E)
      SUM_V2  = SUM_V2  + Fnu_d_2 * W3_S(iN_E)
      SUM_V3  = SUM_V3  + Fnu_d_3 * W3_S(iN_E)

    END DO
    END DO

      SUM_Y = ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Y

      SJ(1) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_1
      SJ(2) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_2
      SJ(3) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_3

      Tau_f = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 - P(iN_X) - cD_old(iN_X)

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_D    (iN_X) &
        = W * U_D(iN_X) 
      C_Y    (iN_X) &
        = U_Y(iN_X) + SUM_Y * S_Y(iN_X)
      C_Ef   (iN_X) &
        = ( Tau_f + SUM_Ef ) * S_Ef  (iN_X) / cD_old(iN_X)
      C_V_d_1(iN_X) &
        = ( SJ(1) + SUM_V1 ) * S_V_d_1(iN_X) / cD_old(iN_X)
      C_V_d_2(iN_X) &
        = ( SJ(2) + SUM_V2 ) * S_V_d_2(iN_X) / cD_old(iN_X)
      C_V_d_3(iN_X) &
        = ( SJ(3) + SUM_V3 ) * S_V_d_3(iN_X) / cD_old(iN_X)

    END DO

  END SUBROUTINE InitializeRHS_Relativistic


  SUBROUTINE ComputeDnuNorm( MASK, Dnu )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu

    INTEGER :: iN_X, iS

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies

      IF ( MASK(iN_X) ) THEN

        DnuNorm(iS,iN_X) = WNORM( Dnu(:,iS,iN_X), W2_S )

      END IF

    END DO
    END DO

  END SUBROUTINE ComputeDnuNorm


  SUBROUTINE ComputeMatterRHS_OrderOne &
    ( MASK, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: SUM_Y, SUM_Ef

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Y, SUM_Ef )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Y, SUM_Ef )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Y, SUM_Ef )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        SUM_Y  = Zero
        SUM_Ef = Zero

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP REDUCTION( + : SUM_Y, SUM_Ef )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR COLLAPSE(2) &
        !$ACC REDUCTION( + : SUM_Y, SUM_Ef )
#endif
        DO iS   = 1, nSpecies
        DO iN_E = 1, nE_G

          IF ( iS <= iNuE_Bar ) THEN
          SUM_Y  = SUM_Y  + Dnu(iN_E,iS,iN_X) * W2_S(iN_E) * LeptonNumber(iS)
          END IF

          SUM_Ef = SUM_Ef + Dnu(iN_E,iS,iN_X) * W3_S(iN_E)

        END DO
        END DO

        G_Y    (iN_X) = C_Y (iN_X) - wMatrRHS(iY ) * SUM_Y  * S_Y (iN_X)
        G_Ef   (iN_X) = C_Ef(iN_X) - wMatrRHS(iEf) * SUM_Ef * S_Ef(iN_X)

        Gm(iY ,iN_X) = G_Y (iN_X)
        Gm(iEf,iN_X) = G_Ef(iN_X)

        Fm(iY ,iN_X) = G_Y (iN_X) - U_Y (iN_X)
        Fm(iEf,iN_X) = G_Ef(iN_X) - U_Ef(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_OrderOne


  SUBROUTINE ComputeMatterRHS_OrderV &
    ( MASK, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3
    REAL(DP) :: SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          V_d_1, V_d_2, V_d_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

        SUM_Y  = Zero
        SUM_Ef = Zero
        SUM_V1 = Zero
        SUM_V2 = Zero
        SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
        !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
        !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
        !$OMP REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR COLLAPSE(2) &
        !$ACC PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
        !$ACC          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
        !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
        !$ACC REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#endif
        DO iS   = 1, nSpecies
        DO iN_E = 1, nE_G

          vDotInu =   V_u_1(iN_X) * Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
                    + V_u_2(iN_X) * Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
                    + V_u_3(iN_X) * Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

#if defined( TWOMOMENT_ORDER_V )
          CALL ComputeEddingtonTensorComponents_dd &
                 ( Dnu    (iN_E,iS,iN_X), Inu_u_1(iN_E,iS,iN_X), &
                   Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                   Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                   k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )
#endif

          vDotK_d_1 &
            = ( V_u_1(iN_X) * k_dd_11 &
              + V_u_2(iN_X) * k_dd_12 &
              + V_u_3(iN_X) * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
          vDotK_d_2 &
            = ( V_u_1(iN_X) * k_dd_12 &
              + V_u_2(iN_X) * k_dd_22 &
              + V_u_3(iN_X) * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
          vDotK_d_3 &
            = ( V_u_1(iN_X) * k_dd_13 &
              + V_u_2(iN_X) * k_dd_23 &
              + V_u_3(iN_X) * k_dd_33 ) * Dnu(iN_E,iS,iN_X)

          ! --- Eulerian Neutrino Number Density ---

          Nnu = Dnu(iN_E,iS,iN_X) + vDotInu

          ! --- Eulerian Neutrino Energy Density Scaled by Neutrino Energy ---

          Enu = Dnu(iN_E,iS,iN_X) + Two * vDotInu

          ! --- Eulerian Neutrino Momentum Density Scaled by Neutrino Energy ---

          Fnu_d_1 &
            = Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
                + V_d_1 * Dnu(iN_E,iS,iN_X) + vDotK_d_1

          Fnu_d_2 &
            = Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
                + V_d_2 * Dnu(iN_E,iS,iN_X) + vDotK_d_2

          Fnu_d_3 &
            = Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) &
                + V_d_3 * Dnu(iN_E,iS,iN_X) + vDotK_d_3

          IF ( iS <= iNuE_Bar ) THEN
          SUM_Y  = SUM_Y  + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
          END IF

          SUM_Ef = SUM_Ef + Enu     * W3_S(iN_E)
          SUM_V1 = SUM_V1 + Fnu_d_1 * W3_S(iN_E)
          SUM_V2 = SUM_V2 + Fnu_d_2 * W3_S(iN_E)
          SUM_V3 = SUM_V3 + Fnu_d_3 * W3_S(iN_E)

        END DO
        END DO

        G_Y    (iN_X) = C_Y    (iN_X) - wMatrRHS(iY ) * SUM_Y  * S_Y    (iN_X)
        G_Ef   (iN_X) = C_Ef   (iN_X) - wMatrRHS(iEf) * SUM_Ef * S_Ef   (iN_X)
        G_V_d_1(iN_X) = C_V_d_1(iN_X) - wMatrRHS(iV1) * SUM_V1 * S_V_d_1(iN_X)
        G_V_d_2(iN_X) = C_V_d_2(iN_X) - wMatrRHS(iV2) * SUM_V2 * S_V_d_2(iN_X)
        G_V_d_3(iN_X) = C_V_d_3(iN_X) - wMatrRHS(iV3) * SUM_V3 * S_V_d_3(iN_X)

        Gm(iY ,iN_X) = G_Y    (iN_X)
        Gm(iEf,iN_X) = G_Ef   (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

        Fm(iY ,iN_X) = G_Y    (iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = G_Ef   (iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = G_V_d_1(iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = G_V_d_2(iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = G_V_d_3(iN_X) - U_V_d_3(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_OrderV


  SUBROUTINE ComputeMatterRHS_Relativistic &
    ( MASK, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, E, P, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: D, E, P, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Alpha, Beta_u_1, Beta_u_2, Beta_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33


    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotV, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, W, h
    REAL(DP) :: DT, B_d_1, B_d_2, B_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3, Nnu
    REAL(DP) :: Inu_d_1, Inu_d_2, Inu_d_3, vDotInu
    REAL(DP) :: SUM_Ef, SUM_V1, SUM_V2, SUM_V3, SUM_Y

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotV, W, B_d_1, B_d_2, B_d_3, DT, h )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          V_d_1, V_d_2, V_d_3, vDotV, W, B_d_1, B_d_2, B_d_3, DT, h )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotV, W, B_d_1, B_d_2, B_d_3, DT, h, &
    !$OMP          Inu_d_1, Inu_d_2, Inu_d_3, k_dd, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 )
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

        vDotV = V_u_1(iN_X) * V_d_1 &
              + V_u_2(iN_X) * V_d_2 &
              + V_u_3(iN_X) * V_d_3 

        W = 1.0_DP / SQRT(1.0_DP - vDotV)

        B_d_1 = Gm_dd_11(iN_X) * Beta_u_1(iN_X)
        B_d_2 = Gm_dd_22(iN_X) * Beta_u_2(iN_X)
        B_d_3 = Gm_dd_33(iN_X) * Beta_u_3(iN_X)

        DT = 1.0_DP / ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) )
     
        SUM_Y   = Zero
        SUM_Ef  = Zero
        SUM_V1  = Zero
        SUM_V2  = Zero
        SUM_V3  = Zero

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( Inu_d_1, Inu_d_2, Inu_d_3, k_dd, &
        !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
        !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
        !$OMP REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR COLLAPSE(2) &
        !$ACC PRIVATE( Inu_d_1, Inu_d_2, Inu_d_3, k_dd, &
        !$ACC          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
        !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
        !$ACC REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#endif
        DO iS   = 1, nSpecies
        DO iN_E = 1, nE_G



          Inu_d_1 = DT * ( B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) & 
                  * Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                  - DT * ( B_d_1 * V_u_2(iN_X) * Gm_dd_22(iN_X) ) * Inu_u_2(iN_E,iS,iN_X)   &
                  - DT * ( B_d_1 * V_u_3(iN_X) * Gm_dd_33(iN_X) ) * Inu_u_3(iN_E,iS,iN_X) 

          Inu_d_2 = DT * ( B_d_1 * V_u_1(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) &
                  * Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                  - DT * ( B_d_2 * V_u_1(iN_X) * Gm_dd_11(iN_X) ) * Inu_u_1(iN_E,iS,iN_X) &
                  - DT * ( Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) * B_d_2 * V_u_3(iN_X) ) 

          Inu_d_3 = DT * ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) - Alpha(iN_X) ) &
                  * Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                  - DT * ( Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) * B_d_3 * V_u_1(iN_X) ) &
                  - DT * ( Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) * B_d_3 * V_u_2(iN_X) )

          vDotInu = V_u_1(iN_X) * Inu_d_1 &
                  + V_u_2(iN_X) * Inu_d_2 &
                  + V_u_3(iN_X) * Inu_d_3

#if defined( TWOMOMENT_RELATIVISTIC )
          CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu    (iN_E,iS,iN_X), Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 Alpha(iN_X),Beta_u_1(iN_X),Beta_u_2(iN_X),Beta_u_3(iN_X), & 
                 V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), k_dd  )
#endif

          vDotK_d_1 &
            = ( V_u_1(iN_X) * k_dd(1,1) &
            + V_u_2(iN_X) * k_dd(2,1) &
            + V_u_3(iN_X) * k_dd(3,1) ) * Dnu(iN_E,iS,iN_X)
          vDotK_d_2 &                                                  
            = ( V_u_1(iN_X) * k_dd(1,2) &                 
            + V_u_2(iN_X) * k_dd(2,2) &                 
            + V_u_3(iN_X) * k_dd(3,2) ) * Dnu(iN_E,iS,iN_X)
          vDotK_d_3 &                                                  
            = ( V_u_1(iN_X) * k_dd(1,3) &                 
            + V_u_2(iN_X) * k_dd(2,3) &                 
            + V_u_3(iN_X) * k_dd(3,3) ) * Dnu(iN_E,iS,iN_X)


          vvDotK = V_u_1(iN_X) * vDotK_d_1 &
                 + V_u_2(iN_X) * vDotK_d_2 & 
                 + V_u_3(iN_X) * vDotK_d_3 

          Nnu = W * Dnu(iN_E,iS,iN_X) + vDotInu

      ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

          Enu = W**2 * Dnu(iN_E,iS,iN_X) +  W * Two * vDotInu + vvDotK 

      ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

          Fnu_d_1 &
            = W * Inu_d_1 + W**2 * V_d_1 * Dnu(iN_E,iS,iN_X) &
            + W * vDotInu * V_d_1  + vDotK_d_1

          Fnu_d_2 &
            = W * Inu_d_2 + W**2 * V_d_2 * Dnu(iN_E,iS,iN_X) &
            + W * vDotInu * V_d_2   + vDotK_d_2

          Fnu_d_3 &
            = W * Inu_d_3 + W**2 * V_d_3 * Dnu(iN_E,iS,iN_X) &
            + W * vDotInu * V_d_3  + vDotK_d_3


          IF ( iS <= iNuE_Bar ) THEN
            SUM_Y  = SUM_Y  + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
          END IF

          SUM_Ef  = SUM_Ef  + Enu     * W3_S(iN_E)
          SUM_V1  = SUM_V1  + Fnu_d_1 * W3_S(iN_E)
          SUM_V2  = SUM_V2  + Fnu_d_2 * W3_S(iN_E)
          SUM_V3  = SUM_V3  + Fnu_d_3 * W3_S(iN_E)

        END DO
        END DO

        h = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) 

        G_D    (iN_X)  = ( 1.0_DP / W ) * C_D  (iN_X) 
        G_Y    (iN_X)  = C_Y    (iN_X) - ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Y * S_Y(iN_X)
        G_Ef   (iN_X)  = S_Ef(iN_X) * ( 1.0_DP - W ) / W**2 * ( W + ( W + 1.0_DP) * P(iN_X) / D(iN_X) ) &
                       + 1.0_DP / W * C_Ef(iN_X) - S_Ef(iN_X) / ( W * cD_old(iN_X) ) * SUM_Ef 
        G_V_d_1(iN_X)  = 1.0_DP / ( h * W ) * (  C_V_d_1(iN_X) - SUM_V1 * S_V_d_1(iN_X) / cD_old(iN_X) ) 
        G_V_d_2(iN_X)  = 1.0_DP / ( h * W ) * (  C_V_d_2(iN_X) - SUM_V2 * S_V_d_2(iN_X) / cD_old(iN_X) ) 
        G_V_d_3(iN_X)  = 1.0_DP / ( h * W ) * (  C_V_d_3(iN_X) - SUM_V3 * S_V_d_3(iN_X) / cD_old(iN_X) ) 

        IF( MoveLeft .EQ. 1) THEN

          G_Ef(iN_X) = ( 1.0_DP / W ) * C_Ef(iN_X) - S_Ef(iN_X) / ( W * cD_old(iN_X) ) * SUM_Ef &
                     + S_Ef(iN_X) * ( 1.0_DP - W ) / W**2 *  ( W + 1.0_DP) * P(iN_X) / D(iN_X) 
        END IF

        Gm(iD ,iN_X) = G_D    (iN_X)
        Gm(iY ,iN_X) = G_Y    (iN_X)
        Gm(iEf,iN_X) = G_Ef   (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

        Fm(iD ,iN_X) = G_D    (iN_X) - U_D  (iN_X)
        Fm(iY ,iN_X) = G_Y    (iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = G_Ef   (iN_X)  - U_Ef  (iN_X)
        Fm(iV1,iN_X) = G_V_d_1(iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = G_V_d_2(iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = G_V_d_3(iN_X) - U_V_d_3(iN_X)

      END IF

    END DO

  END SUBROUTINE ComputeMatterRHS_Relativistic 


  SUBROUTINE ComputeNeutrinoRHS_OrderOne &
    ( MASK, Fm, Gm, dt, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: Eta_T, Chi_T, Kappa
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: L_N, L_G1, L_G2, L_G3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Eta_T, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          Eta_T, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$ACC          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Eta_T, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS )
#endif

    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

#if defined( TWOMOMENT_ORDER_1 )
        CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu    (iN_E,iS,iN_X), Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )
#endif

        ! --- Emissivity ---

        Eta_T =   Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) &
                + Eta_NES (iN_E,iS,iN_X) &
                + Eta_Pair(iN_E,iS,iN_X) &
                + Eta_Brem(iN_E,iS,iN_X)

        ! --- Number Opacity ---

        Chi_T =   Chi_EmAb(iN_E,iS,iN_X) &
                + Chi_NES (iN_E,iS,iN_X) &
                + Chi_Pair(iN_E,iS,iN_X) &
                + Chi_Brem(iN_E,iS,iN_X)

        ! --- Number Flux Opacity ---

        Kappa = Chi_T + Sigma_Iso(iN_E,iN_X)

        ! --- Linear Correction Terms ---

        L_u_1 =   (L_NES__In__u_1(iN_E,iS,iN_X)-L_NES__Out_u_1(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_1(iN_E,iS,iN_X)-L_Pair_Ann_u_1(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_1(iN_E,iS,iN_X)-L_Brem_Ann_u_1(iN_E,iS,iN_X))

        L_u_2 =   (L_NES__In__u_2(iN_E,iS,iN_X)-L_NES__Out_u_2(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_2(iN_E,iS,iN_X)-L_Pair_Ann_u_2(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_2(iN_E,iS,iN_X)-L_Brem_Ann_u_2(iN_E,iS,iN_X))

        L_u_3 =   (L_NES__In__u_3(iN_E,iS,iN_X)-L_NES__Out_u_3(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_3(iN_E,iS,iN_X)-L_Pair_Ann_u_3(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_3(iN_E,iS,iN_X)-L_Brem_Ann_u_3(iN_E,iS,iN_X))

        L_d_1 = Gm_dd_11(iN_X) &
                  * (   L_NES__In__u_1(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_1(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_1(iN_E,iS,iN_X) )

        L_d_2 = Gm_dd_22(iN_X) &
                  * (   L_NES__In__u_2(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_2(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_2(iN_E,iS,iN_X) )

        L_d_3 = Gm_dd_33(iN_X) &
                  * (   L_NES__In__u_3(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_3(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_3(iN_E,iS,iN_X) )

        L_N  = - (   L_u_1 * Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
                   + L_u_2 * Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
                   + L_u_3 * Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) )

        L_G1 = Third * L_d_1 &
               - (   L_u_1 * k_dd_11 &
                   + L_u_2 * k_dd_12 &
                   + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)

        L_G2 = Third * L_d_2 &
               - (   L_u_1 * k_dd_12 &
                   + L_u_2 * k_dd_22 &
                   + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)

        L_G3 = Third * L_d_3 &
               - (   L_u_1 * k_dd_13 &
                   + L_u_2 * k_dd_23 &
                   + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! ! --- Number Equation ---
        !
        ! Gm(iOS+iCR_N,iN_X) &
        !   = ( One - Omega(iN_X) ) * J(iN_E,iS,iN_X) &
        !     + Omega(iN_X) &
        !         * ( C_J(iN_E,iS,iN_X) - vDotH + dt * ( Eta_T + L_N ) ) &
        !         / ( One + dt * Chi_T )
        !
        ! ! --- Number Flux 1 Equation ---
        !
        ! Gm(iOS+iCR_G1,iN_X) &
        !   = ( One - Omega(iN_X) ) * H_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
        !     + Omega(iN_X) &
        !         * ( C_H_d_1(iN_E,iS,iN_X) - vDotK_d_1 + dt * L_G1 ) &
        !         / ( One + dt * Kappa )
        !
        ! ! --- Number Flux 2 Equation ---
        !
        ! Gm(iOS+iCR_G2,iN_X) &
        !   = ( One - Omega(iN_X) ) * H_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
        !     + Omega(iN_X) &
        !         * ( C_H_d_2(iN_E,iS,iN_X) - vDotK_d_2 + dt * L_G2 ) &
        !         / ( One + dt * Kappa )
        !
        ! ! --- Number Flux 3 Equation ---
        !
        ! Gm(iOS+iCR_G3,iN_X) &
        !   = ( One - Omega(iN_X) ) * H_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) &
        !     + Omega(iN_X) &
        !         * ( C_H_d_3(iN_E,iS,iN_X) - vDotK_d_3 + dt * L_G3 ) &
        !         / ( One + dt * Kappa )

        ! --- Number Equation ---

        Gm(iOS+iCR_N,iN_X) &
          = ( C_Dnu(iN_E,iS,iN_X) + dt * ( Eta_T + L_N ) ) &
            / ( One + dt * Chi_T )

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( C_Inu_d_1(iN_E,iS,iN_X) + dt * L_G1 ) / ( One + dt * Kappa )

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( C_Inu_d_2(iN_E,iS,iN_X) + dt * L_G2 ) / ( One + dt * Kappa )

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( C_Inu_d_3(iN_E,iS,iN_X) + dt * L_G3 ) / ( One + dt * Kappa )

        Fm(iOS+iCR_N ,iN_X) &
          = Gm(iOS+iCR_N ,iN_X) - Dnu    (iN_E,iS,iN_X)

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_OrderOne


  SUBROUTINE ComputeNeutrinoRHS_OrderV &
    ( MASK, Fm, Gm, dt, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: Eta_T, Chi_T, Kappa
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: L_N, L_G1, L_G2, L_G3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Eta_T, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          Eta_T, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$ACC          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Eta_T, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS )
#endif

    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        vDotInu =   V_u_1(iN_X) * Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
                  + V_u_2(iN_X) * Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
                  + V_u_3(iN_X) * Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

#if defined( TWOMOMENT_ORDER_V )
        CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu    (iN_E,iS,iN_X), Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )
#endif

        vDotK_d_1 &
          = ( V_u_1(iN_X) * k_dd_11 &
            + V_u_2(iN_X) * k_dd_12 &
            + V_u_3(iN_X) * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
        vDotK_d_2 &
          = ( V_u_1(iN_X) * k_dd_12 &
            + V_u_2(iN_X) * k_dd_22 &
            + V_u_3(iN_X) * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
        vDotK_d_3 &
          = ( V_u_1(iN_X) * k_dd_13 &
            + V_u_2(iN_X) * k_dd_23 &
            + V_u_3(iN_X) * k_dd_33 ) * Dnu(iN_E,iS,iN_X)

        ! --- Emissivity ---

        Eta_T =   Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) &
                + Eta_NES (iN_E,iS,iN_X) &
                + Eta_Pair(iN_E,iS,iN_X) &
                + Eta_Brem(iN_E,iS,iN_X)

        ! --- Number Opacity ---

        Chi_T =   Chi_EmAb(iN_E,iS,iN_X) &
                + Chi_NES (iN_E,iS,iN_X) &
                + Chi_Pair(iN_E,iS,iN_X) &
                + Chi_Brem(iN_E,iS,iN_X)

        ! --- Number Flux Opacity ---

        Kappa = Chi_T + Sigma_Iso(iN_E,iN_X)

        ! --- Linear Correction Terms ---

        L_u_1 =   (L_NES__In__u_1(iN_E,iS,iN_X)-L_NES__Out_u_1(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_1(iN_E,iS,iN_X)-L_Pair_Ann_u_1(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_1(iN_E,iS,iN_X)-L_Brem_Ann_u_1(iN_E,iS,iN_X))

        L_u_2 =   (L_NES__In__u_2(iN_E,iS,iN_X)-L_NES__Out_u_2(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_2(iN_E,iS,iN_X)-L_Pair_Ann_u_2(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_2(iN_E,iS,iN_X)-L_Brem_Ann_u_2(iN_E,iS,iN_X))

        L_u_3 =   (L_NES__In__u_3(iN_E,iS,iN_X)-L_NES__Out_u_3(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_3(iN_E,iS,iN_X)-L_Pair_Ann_u_3(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_3(iN_E,iS,iN_X)-L_Brem_Ann_u_3(iN_E,iS,iN_X))

        L_d_1 = Gm_dd_11(iN_X) &
                  * (   L_NES__In__u_1(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_1(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_1(iN_E,iS,iN_X) )

        L_d_2 = Gm_dd_22(iN_X) &
                  * (   L_NES__In__u_2(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_2(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_2(iN_E,iS,iN_X) )

        L_d_3 = Gm_dd_33(iN_X) &
                  * (   L_NES__In__u_3(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_3(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_3(iN_E,iS,iN_X) )

        L_N  = - (   L_u_1 * Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
                   + L_u_2 * Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
                   + L_u_3 * Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) )

        L_G1 = Third * L_d_1 &
               - (   L_u_1 * k_dd_11 &
                   + L_u_2 * k_dd_12 &
                   + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)

        L_G2 = Third * L_d_2 &
               - (   L_u_1 * k_dd_12 &
                   + L_u_2 * k_dd_22 &
                   + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)

        L_G3 = Third * L_d_3 &
               - (   L_u_1 * k_dd_13 &
                   + L_u_2 * k_dd_23 &
                   + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! --- Number Equation ---

        Gm(iOS+iCR_N,iN_X) &
          = ( ( One - Omega(iN_X) ) * Dnu(iN_E,iS,iN_X) &
              + Omega(iN_X) &
                  * ( C_Dnu(iN_E,iS,iN_X) - vDotInu + dt * ( Eta_T + L_N ) ) ) &
            / ( One + Omega(iN_X) * dt * Chi_T )

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( ( One - Omega(iN_X) ) * Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X) &
              + Omega(iN_X) &
                  * ( C_Inu_d_1(iN_E,iS,iN_X) - vDotK_d_1 + dt * L_G1 ) ) &
            / ( One + Omega(iN_X) * dt * Kappa )

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( ( One - Omega(iN_X) ) * Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X) &
              + Omega(iN_X) &
                  * ( C_Inu_d_2(iN_E,iS,iN_X) - vDotK_d_2 + dt * L_G2 ) ) &
            / ( One + Omega(iN_X) * dt * Kappa )

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( ( One - Omega(iN_X) ) * Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X) &
              + Omega(iN_X) &
                  * ( C_Inu_d_3(iN_E,iS,iN_X) - vDotK_d_3 + dt * L_G3 ) ) &
            / ( One + Omega(iN_X) * dt * Kappa )

        Fm(iOS+iCR_N ,iN_X) &
          = Gm(iOS+iCR_N ,iN_X) - Dnu    (iN_E,iS,iN_X)

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_OrderV


  SUBROUTINE ComputeNeutrinoRHS_Relativistic &
    ( MASK, Fm, Gm, dt, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3  )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Alpha, Beta_u_1, Beta_u_2, Beta_u_3


    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: k_dd(3,3), vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, W, vDotV 
    REAL(DP) :: Eta_T, Chi_T, Kappa
    REAL(DP) :: Inu_d_1, Inu_d_2, Inu_d_3, B_d_1, B_d_2, B_d_3, Det


#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, Det, Inu_d_1, Inu_d_2, Inu_d_3, &
    !$OMP          vDotInu, W, k_dd, vDotV, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          Eta_T, Chi_T, Kappa, iOS )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( B_d_1, B_d_2, B_d_3, Det, Inu_d_1, Inu_d_2, Inu_d_3, &
    !$ACC          vDotInu, W, k_dd, vDotV, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          Eta_T, Chi_T, Kappa, iOS )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, Det, Inu_d_1, Inu_d_2, Inu_d_3, &
    !$OMP          vDotInu, W, k_dd, vDotV, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          Eta_T, Chi_T, Kappa, iOS )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        B_d_1 = Gm_dd_11(iN_X) * Beta_u_1(iN_X)
        B_d_2 = Gm_dd_22(iN_X) * Beta_u_2(iN_X)
        B_d_3 = Gm_dd_33(iN_X) * Beta_u_3(iN_X)

        Det = 1.0_DP / ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) )

        Inu_d_1 = Det * ( B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) & 
                * Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                - Det * ( B_d_1 * V_u_2(iN_X) *Gm_dd_22(iN_X) ) * Inu_u_2(iN_E,iS,iN_X)   &
                - Det * ( B_d_1 * V_u_3(iN_X) * Gm_dd_33(iN_X) ) * Inu_u_3(iN_E,iS,iN_X) 

        Inu_d_2 = Det * ( B_d_1 * V_u_1(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) &
                * Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                - Det * ( B_d_2 * V_u_1(iN_X) * Gm_dd_11(iN_X) ) * Inu_u_1(iN_E,iS,iN_X) &
                - Det * ( Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) * B_d_2 * V_u_3(iN_X) ) 

        Inu_d_3 = Det * ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) - Alpha(iN_X) )& 
                * Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                - Det * ( Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) * B_d_3 * V_u_1(iN_X) ) &
                - Det * ( Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) * B_d_3 * V_u_2(iN_X) )

        vDotInu = V_u_1(iN_X) * Inu_d_1 &
                + V_u_2(iN_X) * Inu_d_2 &
                + V_u_3(iN_X) * Inu_d_3

#if defined( TWOMOMENT_RELATIVISTIC )
       CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu(iN_E,iS,iN_X), Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 Alpha(iN_X), Beta_u_1(iN_X), Beta_u_2(iN_X), Beta_u_3(iN_X), & 
                 V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), k_dd  )
#endif

        vDotK_d_1 &
          = ( V_u_1(iN_X) * k_dd(1,1) &
            + V_u_2(iN_X) * k_dd(2,1) &
            + V_u_3(iN_X) * k_dd(3,1) ) * Dnu(iN_E,iS,iN_X)
        vDotK_d_2 &                                                  
          = ( V_u_1(iN_X) * k_dd(1,2) &
            + V_u_2(iN_X) * k_dd(2,2) &
            + V_u_3(iN_X) * k_dd(3,2) ) * Dnu(iN_E,iS,iN_X)
        vDotK_d_3 &                                                  
          = ( V_u_1(iN_X) * k_dd(1,3) &
            + V_u_2(iN_X) * k_dd(2,3) &
            + V_u_3(iN_X) * k_dd(3,3) ) * Dnu(iN_E,iS,iN_X)

        vDotV =   V_u_1(iN_X) * V_u_1(iN_X) * Gm_dd_11(iN_X) &
                + V_u_2(iN_X) * V_u_2(iN_X) * Gm_dd_22(iN_X) &
                + V_u_3(iN_X) * V_u_3(iN_X) * Gm_dd_33(iN_X)

        W = 1.0_DP / SQRT(1.0_DP - vDotV) 

        ! --- Emissivity ---

        Eta_T =   Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) &
                + Eta_NES (iN_E,iS,iN_X) &
                + Eta_Pair(iN_E,iS,iN_X) &
                + Eta_Brem(iN_E,iS,iN_X)

        ! --- Number Opacity ---

        Chi_T =   Chi_EmAb(iN_E,iS,iN_X) &
                + Chi_NES (iN_E,iS,iN_X) &
                + Chi_Pair(iN_E,iS,iN_X) &
                + Chi_Brem(iN_E,iS,iN_X)

        ! --- Number Flux Opacity ---

        Kappa = Chi_T + Sigma_Iso(iN_E,iN_X)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! --- Number Equation ---

        Gm(iOS+iCR_N,iN_X) &
          = ( One - Omega(iN_X) ) *     Dnu(iN_E,iS,iN_X) &
           + Omega(iN_X)   *  ( C_Dnu(iN_E,iS,iN_X) - vDotInu + dt * Eta_T ) &
              / ( W + dt * Chi_T )

        Fm(iOS+iCR_N,iN_X) &
          = Gm(iOS+iCR_N,iN_X) - Dnu(iN_E,iS,iN_X)

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( One - Omega(iN_X) ) *     Inu_d_1 &
            +       Omega(iN_X)   * ( C_Inu_d_1(iN_E,iS,iN_X) - vDotK_d_1 ) &
                                    / ( W + dt * Kappa )

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - Inu_d_1

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( One - Omega(iN_X) ) *     Inu_d_2 &
            +       Omega(iN_X)   * ( C_Inu_d_2(iN_E,iS,iN_X) - vDotK_d_2 ) &
                                    / ( W + dt * Kappa )

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - Inu_d_2

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( One - Omega(iN_X) ) *     Inu_d_3 &
            +       Omega(iN_X)   * ( C_Inu_d_3(iN_E,iS,iN_x) - vDotK_d_3 ) &
                                    / ( W + dt * Kappa )

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - Inu_d_3
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_Relativistic


  SUBROUTINE UpdateMatterRHS_OrderOne &
    ( MASK, Fm, Gm, Y, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: Y, E
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER :: iN_X

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        Fm(iY ,iN_X) = Gm(iY ,iN_X) - U_Y (iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef(iN_X)

        U_Y   (iN_X) = Gm(iY ,iN_X)
        U_Ef  (iN_X) = Gm(iEf,iN_X)

        Y     (iN_X) = U_Y (iN_X) * Y_old (iN_X)
        E     (iN_X) = U_Ef(iN_X) * Ef_old(iN_X)

      END IF

    END DO

  END SUBROUTINE UpdateMatterRHS_OrderOne


  SUBROUTINE UpdateMatterRHS_OrderV &
    ( MASK, Fm, Gm, Y, E, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_X
    REAL(DP) :: vDotV, V_d_1, V_d_2, V_d_3, Ef

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( vDotV, V_d_1, V_d_2, V_d_3, Ef )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( vDotV, V_d_1, V_d_2, V_d_3, Ef )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( vDotV, V_d_1, V_d_2, V_d_3, Ef )
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        Fm(iY ,iN_X) = Gm(iY ,iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)

        U_Y    (iN_X) = Gm(iY ,iN_X)
        U_Ef   (iN_X) = Gm(iEf,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        Y(iN_X) = U_Y    (iN_X) * Y_old (iN_X)
        Ef      = U_Ef   (iN_X) * Ef_old(iN_X)
        V_d_1   = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2   = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3   = U_V_d_3(iN_X) * SpeedOfLight

        ! --- Update Three-Velocity (Index Up) ---

        V_u_1(iN_X) = V_d_1 / Gm_dd_11(iN_X)
        V_u_2(iN_X) = V_d_2 / Gm_dd_22(iN_X)
        V_u_3(iN_X) = V_d_3 / Gm_dd_33(iN_X)

        ! --- Compute E from Ef ---

        vDotV =   V_u_1(iN_X) * V_d_1 &
                + V_u_2(iN_X) * V_d_2 &
                + V_u_3(iN_X) * V_d_3

        E(iN_X) = Ef - Half * vDotV

        ! --- Compute Omega (Richardson damping coeff.) based on velocity ---

        Omega(iN_X) = One / ( One + SQRT( vDotV ) )

      END IF

    END DO

  END SUBROUTINE UpdateMatterRHS_OrderV


  SUBROUTINE UpdateMatterRHS_Relativistic &
    ( MASK, Fm, Gm, D, Y, E, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: D, Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_X
    REAL(DP) :: vDotV, W, V_d_1, V_d_2, V_d_3, Ef_temp, P
    
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Ef_temp, V_d_1, V_d_2, V_d_3, vDotV, W )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( Ef_temp, V_d_1, V_d_2, V_d_3, vDotV, W )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Ef_temp, V_d_1, V_d_2, V_d_3, vDotV, W )
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        Fm(iD ,iN_X) = Gm(iD ,iN_X) - U_D    (iN_X)
        Fm(iY ,iN_X) = Gm(iY ,iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)

        U_D    (iN_X) = Gm(iD ,iN_X)
        U_Y    (iN_X) = Gm(iY ,iN_X)
        U_Ef   (iN_X) = Gm(iEf,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        D(iN_X) = U_D    (iN_X) * D_old (iN_X)
        Y(iN_X) = U_Y    (iN_X) * Y_old (iN_X)
        Ef_temp = U_Ef   (iN_X) * Ef_old(iN_X)
        V_d_1   = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2   = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3   = U_V_d_3(iN_X) * SpeedOfLight

        ! --- Update Three-Velocity (Index Up) ---

        V_u_1(iN_X) = V_d_1 / Gm_dd_11(iN_X)
        V_u_2(iN_X) = V_d_2 / Gm_dd_22(iN_X)
        V_u_3(iN_X) = V_d_3 / Gm_dd_33(iN_X)

        vDotV =   V_u_1(iN_X) * V_d_1 &
                + V_u_2(iN_X) * V_d_2 &
                + V_u_3(iN_X) * V_d_3

        W = 1.0_DP / SQRT( 1.0_DP - vDotV)

        ! --- Compute Omega (Richardson damping coeff.) based on velocity ---

        Omega(iN_X) = One / ( W + SQRT( vDotV ) )

        E(iN_X) = Ef_temp

        IF ( MoveLeft .EQ. 1 ) THEN

          E(iN_X) = Ef_temp - ( W - One ) / W

        END IF

      END IF

    END DO
  END SUBROUTINE UpdateMatterRHS_Relativistic 


  SUBROUTINE UpdateNeutrinoRHS_OrderOne &
    ( MASK, Fm, Gm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iOS )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iOS )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iOS )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        Fm(iOS+iCR_N ,iN_X) &
          = Gm(iOS+iCR_N ,iN_X) - Dnu    (iN_E,iS,iN_X)

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

        Dnu    (iN_E,iS,iN_X) = Gm(iOS+iCR_N ,iN_X)
        Inu_u_1(iN_E,iS,iN_X) = Gm(iOS+iCR_G1,iN_X) / Gm_dd_11(iN_X)
        Inu_u_2(iN_E,iS,iN_X) = Gm(iOS+iCR_G2,iN_X) / Gm_dd_22(iN_X)
        Inu_u_3(iN_E,iS,iN_X) = Gm(iOS+iCR_G3,iN_X) / Gm_dd_33(iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE UpdateNeutrinoRHS_OrderOne


  SUBROUTINE UpdateNeutrinoRHS_OrderV &
    ( MASK, Fm, Gm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iOS )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iOS )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iOS )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        Fm(iOS+iCR_N ,iN_X) &
          = Gm(iOS+iCR_N ,iN_X) - Dnu    (iN_E,iS,iN_X)

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

        Dnu    (iN_E,iS,iN_X) = Gm(iOS+iCR_N ,iN_X)
        Inu_u_1(iN_E,iS,iN_X) = Gm(iOS+iCR_G1,iN_X) / Gm_dd_11(iN_X)
        Inu_u_2(iN_E,iS,iN_X) = Gm(iOS+iCR_G2,iN_X) / Gm_dd_22(iN_X)
        Inu_u_3(iN_E,iS,iN_X) = Gm(iOS+iCR_G3,iN_X) / Gm_dd_33(iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE UpdateNeutrinoRHS_OrderV


  SUBROUTINE UpdateNeutrinoRHS_Relativistic &
    ( MASK, Fm, Gm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Alpha, Beta_u_1, Beta_u_2, Beta_u_3


    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: B_d_1, B_d_2, B_d_3
    REAL(DP) :: Inu_d_1, Inu_d_2, Inu_d_3, DT

 
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iOS, B_d_1, B_d_2, B_d_3, DT, Inu_d_1, Inu_d_2, Inu_d_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iOS, B_d_1, B_d_2, B_d_3, DT, Inu_d_1, Inu_d_2, Inu_d_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iOS, B_d_1, B_d_2, B_d_3, DT, Inu_d_1, Inu_d_2, Inu_d_3 )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR


        B_d_1 = Gm_dd_11(iN_X) * Beta_u_1(iN_X)
        B_d_2 = Gm_dd_22(iN_X) * Beta_u_2(iN_X)
        B_d_3 = Gm_dd_33(iN_X) * Beta_u_3(iN_X)

        DT = 1.0_DP / ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) )

        Inu_d_1 = DT * ( B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) & 
                * Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                - DT * ( B_d_1 * V_u_2(iN_X) *Gm_dd_22(iN_X) ) * Inu_u_2(iN_E,iS,iN_X)   &
                - DT * ( B_d_1 * V_u_3(iN_X) * Gm_dd_33(iN_X) ) * Inu_u_3(iN_E,iS,iN_X) 

        Inu_d_2 = DT * ( B_d_1 * V_u_1(iN_X) + B_d_3 * V_u_3(iN_X) - Alpha(iN_X) ) &
                * Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                - DT * ( B_d_2 * V_u_1(iN_X) * Gm_dd_11(iN_X) ) * Inu_u_1(iN_E,iS,iN_X) &
                - DT * ( Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) * B_d_2 * V_u_3(iN_X) ) 

        Inu_d_3 = DT * ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) - Alpha(iN_X) ) &
                * Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                - DT * ( Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) * B_d_3 * V_u_1(iN_X) ) &
                - DT * ( Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) * B_d_3 * V_u_2(iN_X) )




        Fm(iOS+iCR_N ,iN_X) = Gm(iOS+iCR_N ,iN_X) - Dnu    (iN_E,iS,iN_X)
        Fm(iOS+iCR_G1,iN_X) = Gm(iOS+iCR_G1,iN_X) - Inu_d_1
        Fm(iOS+iCR_G2,iN_X) = Gm(iOS+iCR_G2,iN_X) - Inu_d_2
        Fm(iOS+iCR_G3,iN_X) = Gm(iOS+iCR_G3,iN_X) - Inu_d_3

        Dnu(iN_E,iS,iN_X) = Gm(iOS+iCR_N ,iN_X)
        Inu_d_1           = Gm(iOS+iCR_G1,iN_X)
        Inu_d_2           = Gm(iOS+iCR_G2,iN_X)
        Inu_d_3           = Gm(iOS+iCR_G3,iN_X)


        Inu_u_1(iN_E,iS,iN_X) = ( 1.0_DP - B_d_1 * V_u_1(iN_X) / Alpha(iN_X) ) * Inu_d_1 / Gm_dd_11(iN_X)  &
                              - Inu_d_2 * B_d_1 * V_u_2(iN_X) / ( Alpha(iN_X) *Gm_dd_11(iN_X) ) &
                              - Inu_d_3 * B_d_1 * V_u_3(iN_X) / ( Gm_dd_11(iN_X) * Alpha(iN_X) )

        Inu_u_2(iN_E,iS,iN_X) =( 1.0_DP -   B_d_2 * V_u_2(iN_X) / Alpha(iN_X) ) * Inu_d_2 / Gm_dd_22(iN_X)  &
                              - Inu_d_1 * B_d_2 * V_u_1(iN_X) / ( Alpha(iN_X) *Gm_dd_22(iN_X) ) &
                              - Inu_d_3 * B_d_2 * V_u_3(iN_X) / ( Gm_dd_22(iN_X) * Alpha(iN_X) )

        Inu_u_3(iN_E,iS,iN_X) =( 1.0_DP - B_d_3 * V_u_3(iN_X) / Alpha(iN_X) ) * Inu_d_3 / Gm_dd_33(iN_X)  &
                              - Inu_d_1 * B_d_3 * V_u_1(iN_X) / ( Alpha(iN_X) *Gm_dd_33(iN_X) ) &
                              - Inu_d_2 * B_d_3 * V_u_2(iN_X) / ( Gm_dd_33(iN_X) * Alpha(iN_X) )


      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE UpdateNeutrinoRHS_Relativistic 


  SUBROUTINE SolveLS_FP &
    ( MASK, n_FP, M, Mk, Fm, Gm, F, G, A, B, Alpha, TAU, LWORK, WORK )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm, B
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F, G, A
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Alpha, TAU
    INTEGER,                    INTENT(in)    :: LWORK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: WORK

    INTEGER  :: iN_X, iFP, iM
    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
    DO iFP  = 1, n_FP
      IF( MASK(iN_X) )THEN
        F(iFP,Mk,iN_X) = Fm(iFP,iN_X)
        G(iFP,Mk,iN_X) = Gm(iFP,iN_X)
      END IF
    END DO
    END DO

    IF ( Mk > 1 ) THEN

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO COLLAPSE(2)
#endif
      DO iN_X = 1, nX_G
      DO iFP  = 1, n_FP
        IF( MASK(iN_X) )THEN
          B(iFP,iN_X) = - Fm(iFP,iN_X)
        END IF
      END DO
      END DO

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
      DO iM   = 1, Mk-1
      DO iFP  = 1, n_FP
        IF ( MASK(iN_X) ) THEN
          A(iFP,iM,iN_X) = F(iFP,iM,iN_X) - Fm(iFP,iN_X)
        END IF
      END DO
      END DO
      END DO

      IF ( Mk == 2 ) THEN

#if   defined( THORNADO_OMP_OL )
        !$OMP TARGET TEAMS DISTRIBUTE &
        !$OMP PRIVATE( AA11, AB1 )
#elif defined( THORNADO_OACC   )
        !$ACC PARALLEL LOOP GANG &
        !$ACC PRIVATE( AA11, AB1 )
#elif defined( THORNADO_OMP    )
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AB1 )
#endif
        DO iN_X = 1, nX_G
          IF( MASK(iN_X) )THEN

            AA11 = Zero
            AB1  = Zero

#if   defined(THORNADO_OMP_OL)
            !$OMP PARALLEL DO SIMD &
            !$OMP REDUCTION( +: AA11, AB1 )
#elif defined(THORNADO_OACC  )
            !$ACC LOOP VECTOR &
            !$ACC REDUCTION( +: AA11, AB1 )
#endif
            DO iFP = 1, n_FP

              AA11 = AA11 + A(iFP,1,iN_X) * A(iFP,1,iN_X)
              AB1  = AB1  + A(iFP,1,iN_X) * B(iFP  ,iN_X)

            END DO

            IF( ABS( AA11 ) < SqrtTiny )THEN

              B(1,iN_X) = Zero

            ELSE

              B(1,iN_X) = AB1 / AA11

            END IF

          END IF
        END DO

      ELSE IF ( Mk == 3 ) THEN

#if defined  ( THORNADO_OMP_OL )
        !$OMP TARGET TEAMS DISTRIBUTE &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#elif defined( THORNADO_OACC   )
        !$ACC PARALLEL LOOP GANG &
        !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#elif defined( THORNADO_OMP    )
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#endif
        DO iN_X = 1, nX_G
          IF ( MASK(iN_X) ) THEN

            AA11 = Zero
            AA12 = Zero
            AA22 = Zero
            AB1  = Zero
            AB2  = Zero

#if   defined(THORNADO_OMP_OL)
            !$OMP PARALLEL DO SIMD &
            !$OMP REDUCTION( +: AA11, AA12, AA22, AB1, AB2 )
#elif defined(THORNADO_OACC  )
            !$ACC LOOP VECTOR &
            !$ACC REDUCTION( +: AA11, AA12, AA22, AB1, AB2 )
#endif
            DO iFP = 1, n_FP

              AA11 = AA11 + A(iFP,1,iN_X) * A(iFP,1,iN_X)
              AA12 = AA12 + A(iFP,1,iN_X) * A(iFP,2,iN_X)
              AA22 = AA22 + A(iFP,2,iN_X) * A(iFP,2,iN_X)

              AB1  = AB1  + A(iFP,1,iN_X) * B(iFP,iN_X)
              AB2  = AB2  + A(iFP,2,iN_X) * B(iFP,iN_X)

            END DO

            DET_AA = AA11 * AA22 - AA12 * AA12

            IF( ABS( DET_AA ) < SqrtTiny )THEN

              B(1,iN_X) = Zero
              B(2,iN_X) = Zero

            ELSE

              B(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
              B(2,iN_X) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

            END IF

          END IF
        END DO

      ELSE IF ( Mk > 3 ) THEN

        DO iN_X = 1, nX_G
          IF( MASK(iN_X) )THEN

            CALL LinearLeastSquares &
                   ( 'N', n_FP, Mk-1, 1, A(:,:,iN_X), n_FP, B(:,iN_X), n_FP, &
                     TAU(:,iN_X), WORK(:,iN_X), LWORK, INFO(iN_X) )

          END IF
        END DO

      END IF

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( SUM1 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( SUM1 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( SUM1 )
#endif
      DO iN_X = 1, nX_G
        IF( MASK(iN_X) )THEN

          SUM1 = Zero
          DO iM = 1, Mk-1
            Alpha(iM,iN_X) = B(iM,iN_X)
            SUM1 = SUM1 + B(iM,iN_X)
          END DO
          Alpha(Mk,iN_X) = One - SUM1

        END IF
      END DO

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( SUM1 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( SUM1 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( SUM1 )
#endif
      DO iN_X = 1, nX_G
      DO iFP  = 1, n_FP
        IF( MASK(iN_X) )THEN

          SUM1 = Zero
          DO iM = 1, Mk
            SUM1 = SUM1 + G(iFP,iM,iN_X) * Alpha(iM,iN_X)
          END DO
          Gm(iFP,iN_X) = SUM1

        END IF
      END DO
      END DO

    ELSE

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G
        IF( MASK(iN_X) )THEN

          Alpha(Mk,iN_X) = One

        END IF
      END DO

    END IF

  END SUBROUTINE SolveLS_FP


  SUBROUTINE ShiftRHS_FP( MASK, n_FP, M, Mk, F, G )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F, G

    INTEGER  :: iN_X, iFP, iM
    REAL(DP) :: FTMP(1:n_FP,1:M)
    REAL(DP) :: GTMP(1:n_FP,1:M)

    IF ( Mk == M ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( FTMP, GTMP )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( FTMP, GTMP )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( FTMP, GTMP )
#endif
      DO iN_X = 1, nX_G
        IF( MASK(iN_X) )THEN

#if   defined( THORNADO_OMP_OL )
          !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM  = 1, Mk-1
          DO iFP = 1, n_FP
            FTMP(iFP,iM) = F(iFP,iM+1,iN_X)
            GTMP(iFP,iM) = G(iFP,iM+1,iN_X)
          END DO
          END DO

#if   defined( THORNADO_OMP_OL )
          !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM  = 1, Mk-1
          DO iFP = 1, n_FP
            F(iFP,iM,iN_X) = FTMP(iFP,iM)
            G(iFP,iM,iN_X) = GTMP(iFP,iM)
          END DO
          END DO

        END IF
      END DO

    END IF

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE CheckConvergence_Inner &
    ( MASK, n_FP, k_inner, nIterations_Inner, Fm )

    LOGICAL,  DIMENSION(:)    , INTENT(inout) :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, k_inner
    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Inner
    REAL(DP), DIMENSION(:,:)  , INTENT(in)    :: Fm

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X, iS, iN_E, iOS
    REAL(DP) :: Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3

    ! --- Compute Norm of F (Use Maximum Across Moments) ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 )
#endif
    DO iN_X = 1, nX_G

      IF( MASK(iN_X) )THEN

        CONVERGED = .TRUE.

        DO iS = 1, nSpecies

          Fnorm_N  = Zero
          Fnorm_G1 = Zero
          Fnorm_G2 = Zero
          Fnorm_G3 = Zero

          DO iN_E = 1, nE_G

            iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

            Fnorm_N  = Fnorm_N  + W2_S(iN_E) * Fm(iOS+iCR_N ,iN_X)**2
            Fnorm_G1 = Fnorm_G1 + W2_S(iN_E) * Fm(iOS+iCR_G1,iN_X)**2
            Fnorm_G2 = Fnorm_G2 + W2_S(iN_E) * Fm(iOS+iCR_G2,iN_X)**2
            Fnorm_G3 = Fnorm_G3 + W2_S(iN_E) * Fm(iOS+iCR_G3,iN_X)**2

          END DO

          Fnorm = SQRT( MAX( Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 ) )

          CONVERGED = CONVERGED .AND. ( Fnorm <= Rtol_inner * DnuNorm(iS,iN_X) )

        END DO

        nIterations_Inner(iN_X) &
          = nIterations_Inner(iN_X) + 1
        IF( CONVERGED )THEN
          MASK(iN_X) = .FALSE.
        END IF

      END IF

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( MASK )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( MASK )
#endif

  END SUBROUTINE CheckConvergence_Inner


  SUBROUTINE CheckConvergence_Outer &
    ( MASK_OUTER, MASK_INNER, n_FP, k_outer, nIterations_Outer, Fm )

    LOGICAL,  DIMENSION(:)    , INTENT(inout) :: MASK_OUTER, MASK_INNER
    INTEGER,                    INTENT(in)    :: n_FP, k_outer
    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Outer
    REAL(DP), DIMENSION(:,:)  , INTENT(in)    :: Fm

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X
    REAL(DP) :: Fnorm

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, Fnorm )
#endif
    DO iN_X = 1, nX_G
      IF( MASK_OUTER(iN_X) )THEN

        Fnorm = MAXVAL( ABS( Fm(:,iN_X) ) )

        CONVERGED = Fnorm <= Rtol_outer

        nIterations_Outer(iN_X) = k_outer
        IF( CONVERGED )THEN
          MASK_OUTER(iN_X) = .FALSE.
        END IF

        MASK_INNER(iN_X) = MASK_OUTER(iN_X)

      END IF
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( MASK_OUTER, MASK_INNER )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( MASK_OUTER, MASK_INNER )
#endif

  END SUBROUTINE CheckConvergence_Outer


!!$  SUBROUTINE CheckConvergence_Outer &
!!$    ( MASK_OUTER, MASK_INNER, n_FP, k_outer, nIterations_Outer, Fm )
!!$
!!$    LOGICAL,  DIMENSION(:)    , INTENT(inout) :: MASK_OUTER, MASK_INNER
!!$    INTEGER,                    INTENT(in)    :: n_FP, k_outer
!!$    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Outer
!!$    REAL(DP), DIMENSION(:,:)  , INTENT(in)    :: Fm
!!$
!!$    LOGICAL  :: CONVERGED
!!$    INTEGER  :: iX, iOS
!!$    REAL(DP) :: Fnorm_Y(nDOFX), Fnorm_Ef(nDOFX), Fnorm_V(nDOFX)
!!$
!!$#if   defined( THORNADO_OMP_OL )
!!$    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
!!$    !$OMP PRIVATE( iOS, CONVERGED, Fnorm_Y, Fnorm_Ef, Fnorm_V )
!!$#elif defined( THORNADO_OACC   )
!!$    !$ACC PARALLEL LOOP GANG VECTOR &
!!$    !$ACC PRIVATE( iOS, CONVERGED, Fnorm_Y, Fnorm_Ef, Fnorm_V )
!!$#elif defined( THORNADO_OMP    )
!!$    !$OMP PARALLEL DO &
!!$    !$OMP PRIVATE( iOS, CONVERGED, Fnorm_Y, Fnorm_Ef, Fnorm_V )
!!$#endif
!!$    DO iX = 1, nX_G / nDOFX
!!$      iOS = (iX-1) * nDOFX
!!$      IF( ANY( MASK_OUTER(iOS+1:iOS+nDOFX) ) )THEN
!!$
!!$        Fnorm_Y  =      ABS( Fm(iY ,iOS+1:iOS+nDOFX) )
!!$        Fnorm_Ef =      ABS( Fm(iEf,iOS+1:iOS+nDOFX) )
!!$        Fnorm_V  = MAX( ABS( Fm(iV1,iOS+1:iOS+nDOFX) ), &
!!$                        ABS( Fm(iV2,iOS+1:iOS+nDOFX) ), &
!!$                        ABS( Fm(iV3,iOS+1:iOS+nDOFX) ) )
!!$
!!$        CONVERGED = ALL( Fnorm_Y  <= Rtol_outer ) .AND. &
!!$                    ALL( Fnorm_Ef <= Rtol_outer ) .AND. &
!!$                    ALL( Fnorm_V  <= Rtol_outer )
!!$
!!$        IF( CONVERGED )THEN
!!$          MASK_OUTER       (iOS+1:iOS+nDOFX) = .FALSE.
!!$          nIterations_Outer(iOS+1:iOS+nDOFX) = k_outer
!!$        END IF
!!$
!!$        MASK_INNER(iOS+1:iOS+nDOFX) = MASK_OUTER(iOS+1:iOS+nDOFX)
!!$
!!$      END IF
!!$    END DO
!!$
!!$#if   defined( THORNADO_OMP_OL )
!!$    !$OMP TARGET UPDATE FROM( MASK_OUTER, MASK_INNER )
!!$#elif defined( THORNADO_OACC   )
!!$    !$ACC UPDATE HOST( MASK_OUTER, MASK_INNER )
!!$#endif
!!$
!!$  END SUBROUTINE CheckConvergence_Outer


  FUNCTION ENORM( X ) ! Delete Me?
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: ENORM
    REAL(DP), INTENT(in) :: X(:)

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM


  FUNCTION WNORM( X, W )
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: WNORM
    REAL(DP), INTENT(in) :: X(:), W(:)

    WNORM = SQRT( DOT_PRODUCT( W * X, X ) )

    RETURN
  END FUNCTION WNORM


  SUBROUTINE CheckErrorFlag_FP &
    ( Error, k_outer, k_inner, &
      D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    USE mpi

    INTEGER,  DIMENSION(:)    , INTENT(in) :: Error
    INTEGER,                    INTENT(in) :: k_outer, k_inner
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Y, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: ierr
    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: D_P, T_P, Y_P, E_P, V1_P, V2_P, V3_P
    REAL(DP) :: D0_P, T0_P, Y0_P, E0_P, V10_P, V20_P, V30_P

    IF ( ANY( Error > 0 ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM &
      !$OMP ( D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !$OMP   D_old, Y_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
      !$OMP   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST &
      !$ACC ( D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !$ACC   D_old, Y_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
      !$ACC   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#endif
      DO iN_X = 1, nX_G
        IF ( Error(iN_X) > 0 ) THEN

          D_P   = D(iN_X) / Unit_D
          Y_P   = Y(iN_X) / Unit_Y
          E_P   = E(iN_X) / Unit_E
          T_P   = T(iN_X) / Unit_T
          V1_P  = V_u_1(iN_X) / Unit_V
          V2_P  = V_u_2(iN_X) / Unit_V
          V3_P  = V_u_3(iN_X) / Unit_V

          D0_P  = D_old(iN_X) / Unit_D
          T0_P  = T_old(iN_X) / Unit_T
          Y0_P  = Y_old(iN_X) / Unit_Y
          E0_P  = E_old(iN_X) / Unit_E
          V10_P = V_u_1_old(iN_X) / Unit_V
          V20_P = V_u_2_old(iN_X) / Unit_V
          V30_P = V_u_3_old(iN_X) / Unit_V

          WRITE(*,*)                     '[SolveNeutrinoMatterCoupling_FP_Nested_AA] Error'
          WRITE(*,'(a,2i5)')             '             iN_X, Error : ', iN_X, Error(iN_X)
          WRITE(*,'(a,5x,2i23)')         '        k_outer, k_inner : ', k_outer, k_inner
          WRITE(*,'(a,5x,7es23.15)')     '   D, Y, E, T, V_u       : ', D_P, Y_P, E_P, T_P, V1_P, V2_P, V3_P
          WRITE(*,'(a,5x,7es23.15)')     '   D, Y, E, T, V_u (old) : ', D0_P, Y0_P, E0_P, T0_P, V10_P, V20_P, V30_P

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Dnu           : ', iS, ( Dnu    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Inu_u_1       : ', iS, ( Inu_u_1(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Inu_u_2       : ', iS, ( Inu_u_2(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Inu_u_3       : ', iS, ( Inu_u_3(iN_E,iS,iN_X), iN_E = 1, nE_G )

          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Dnu     (old) : ', iS, ( Dnu_old    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Inu_u_1 (old) : ', iS, ( Inu_u_1_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Inu_u_2 (old) : ', iS, ( Inu_u_2_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '      iS, Inu_u_3 (old) : ', iS, ( Inu_u_3_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

        END IF
      END DO
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    END IF

  END SUBROUTINE CheckErrorFlag_FP


END MODULE TwoMoment_NeutrinoMatterSolverModule
