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
    Kelvin, &
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
    Timer_Collisions_MatterFromRadiation, &
    Timer_Collisions_RadiationBacktrack, &
    Timer_Opacity_D0, &
    Timer_Opacity_LimitD0, &
    Timer_Opacity_EC, &
    Timer_Opacity_ES, &
    Timer_Opacity_NNS, &
    Timer_Opacity_NES, &
    Timer_Opacity_Pair, &
    Timer_Opacity_NuPair, &
    Timer_Opacity_Brem, &
    Timer_OpacityRate_NNS, &
    Timer_OpacityRate_NES, &
    Timer_OpacityRate_Pair, &
    Timer_OpacityRate_NuPair, &
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
    iNu, iNu_Bar, LeptonNumber, nChirals, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputePressure_TABLE, &
    Min_D, Min_T, Min_Y
  USE OpacityModule_TABLE, ONLY: &
    IsoIncludesNucleonScattering
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, &
    LimitEquilibriumDistributions_DG, &
    ComputeNeutrinoOpacities_EC, &
    ComputeNeutrinoOpacities_ES, &
    ComputeNeutrinoOpacities_NNS, &
    ComputeNeutrinoOpacities_NES, &
    ComputeNeutrinoOpacities_Pair, &
    ComputeNeutrinoOpacities_NuPair, &
    ComputeNeutrinoOpacities_Brem, &
    ComputeNeutrinoOpacityRates_NNS, &
    ComputeNeutrinoOpacityRates_NES, &
    ComputeNeutrinoOpacityRates_Pair, &
    ComputeNeutrinoOpacityRates_NuPair, &
    ComputeNeutrinoOpacityRates_Brem, &
    ComputeNeutrinoOpacityRates_LinearCorrections_NNS, &
    ComputeNeutrinoOpacityRates_LinearCorrections_NES, &
    ComputeNeutrinoOpacityRates_LinearCorrections_Pair, &
    ComputeNeutrinoOpacityRates_LinearCorrections_NuPair, &
    ComputeNeutrinoOpacityRates_LinearCorrections_Brem, &
    SetOpacityMask, &
    iM0, iM1, &
    nOp, iOp_EmAb, iOp_ECTable, iOp_Iso, &
    iOp_NNS, iOp_NES, iOp_Pair, iOp_NuPair, iOp_Brem
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeEddingtonTensorComponents_dd, &
    ComputeConserved_TwoMoment
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: InitializeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling_FP_Nested_AA
#if defined( TWOMOMENT_ORDER_V )
  PUBLIC :: ComputeOpacities_Packed
  PUBLIC :: ComputeRates_Packed
  PUBLIC :: ComputeTransportOpacity
  PUBLIC :: AssembleTransportOpacity
  PUBLIC :: ComputeMatterRHS_Collisions_OrderV
  PUBLIC :: ComputeMeanInverseMeanFreePaths
  PUBLIC :: Include_NNS, Include_NES, Include_Pair, Include_NuPair, Include_Brem
  PUBLIC :: nQ
  PUBLIC :: iQ_L, iQ_E, iQ_H, iQ_M1, iQ_M2, iQ_M3
  PUBLIC :: nNNSLedger, iNNSLedger_NIn, iNNSLedger_NOut, iNNSLedger_NNet, &
            iNNSLedger_EIn, iNNSLedger_EOut, iNNSLedger_QH
  PUBLIC :: nEmAbLedger, iEmAbLedger_NIn, iEmAbLedger_NOut, &
            iEmAbLedger_NNet, iEmAbLedger_EIn, iEmAbLedger_EOut, &
            iEmAbLedger_QH
  PUBLIC :: nEmAbSpectrum, iEmAbSpectrum_Dnu, iEmAbSpectrum_DnuEq, &
            iEmAbSpectrum_Chi, iEmAbSpectrum_EIn, &
            iEmAbSpectrum_EOut, iEmAbSpectrum_QH
  PUBLIC :: nC
  PUBLIC :: iC_EmAb, iC_Iso, iC_NNS, iC_NES
  PUBLIC :: iC_Pair, iC_NuPair, iC_Brem

  ! Matter collision RHS quantities.
  INTEGER, PARAMETER :: iQ_L  = 1
  INTEGER, PARAMETER :: iQ_E  = 2
  INTEGER, PARAMETER :: iQ_H  = 3
  INTEGER, PARAMETER :: iQ_M1 = 4
  INTEGER, PARAMETER :: iQ_M2 = 5
  INTEGER, PARAMETER :: iQ_M3 = 6
  INTEGER, PARAMETER :: nQ    = 6

  INTEGER, PARAMETER :: iNNSLedger_NIn  = 1
  INTEGER, PARAMETER :: iNNSLedger_NOut = 2
  INTEGER, PARAMETER :: iNNSLedger_NNet = 3
  INTEGER, PARAMETER :: iNNSLedger_EIn  = 4
  INTEGER, PARAMETER :: iNNSLedger_EOut = 5
  INTEGER, PARAMETER :: iNNSLedger_QH   = 6
  INTEGER, PARAMETER :: nNNSLedger      = 6

  INTEGER, PARAMETER :: iEmAbLedger_NIn  = 1
  INTEGER, PARAMETER :: iEmAbLedger_NOut = 2
  INTEGER, PARAMETER :: iEmAbLedger_NNet = 3
  INTEGER, PARAMETER :: iEmAbLedger_EIn  = 4
  INTEGER, PARAMETER :: iEmAbLedger_EOut = 5
  INTEGER, PARAMETER :: iEmAbLedger_QH   = 6
  INTEGER, PARAMETER :: nEmAbLedger      = 6

  INTEGER, PARAMETER :: iEmAbSpectrum_Dnu   = 1
  INTEGER, PARAMETER :: iEmAbSpectrum_DnuEq = 2
  INTEGER, PARAMETER :: iEmAbSpectrum_Chi   = 3
  INTEGER, PARAMETER :: iEmAbSpectrum_EIn   = 4
  INTEGER, PARAMETER :: iEmAbSpectrum_EOut  = 5
  INTEGER, PARAMETER :: iEmAbSpectrum_QH    = 6
  INTEGER, PARAMETER :: nEmAbSpectrum       = 6

  ! Collision-process decomposition.
  INTEGER, PARAMETER :: iC_EmAb   = 1
  INTEGER, PARAMETER :: iC_Iso    = 2
  INTEGER, PARAMETER :: iC_NNS    = 3
  INTEGER, PARAMETER :: iC_NES    = 4
  INTEGER, PARAMETER :: iC_Pair   = 5
  INTEGER, PARAMETER :: iC_NuPair = 6
  INTEGER, PARAMETER :: iC_Brem   = 7
  INTEGER, PARAMETER :: nC        = 7

  INTEGER, PARAMETER, PUBLIC :: nLeptonFamilies = 3
  INTEGER, PARAMETER, PUBLIC :: iLFA_InputConversion  = 1
  INTEGER, PARAMETER, PUBLIC :: iLFA_EntryProjection  = 2
  INTEGER, PARAMETER, PUBLIC :: iLFA_CollisionSource  = 3
  INTEGER, PARAMETER, PUBLIC :: iLFA_RawMap           = 4
  INTEGER, PARAMETER, PUBLIC :: iLFA_MapProjection    = 5
  INTEGER, PARAMETER, PUBLIC :: iLFA_FixedPointResidual = 6
  INTEGER, PARAMETER, PUBLIC :: iLFA_EquationResidual = 7
  INTEGER, PARAMETER, PUBLIC :: nLeptonFlavorAudit    = 7
#endif

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

  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: DnuNorm, Fnorm0
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
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Sigma_Iso
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_0_Iso
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_1_Iso
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_EmAb, Eta_EmAb
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NNS , Eta_NNS, R_NNS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NES , Eta_NES, R_NES
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Pair, Eta_Pair, R_Pair
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NuPair, Eta_NuPair, R_NuPair
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Brem, Eta_Brem, R_Brem
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__In__u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__In__u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__In__u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__Out_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__Out_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__Out_u_3
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
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Pro_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Pro_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Pro_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Ann_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Ann_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Ann_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_3

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_NNS_0, Phi_NbNS_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_NNS_1, Phi_NbNS_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_0, H_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_1, H_II_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_0, J_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_1, J_II_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Nu_J_I_0, Nu_J_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Nu_J_I_1, Nu_J_II_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: S_Sigma

  LOGICAL,  DIMENSION(:,:),   ALLOCATABLE, TARGET :: Opacity_Mask

  INTEGER :: LWORK_outer
  INTEGER :: LWORK_inner

  INTEGER,  DIMENSION(:)    , ALLOCATABLE :: INFO

  ! --- Solver Parameters to be initialized

  LOGICAL :: SolverParametersInitialized = .FALSE.

  LOGICAL  :: Include_NNS
  LOGICAL  :: NNS_ApplyManyBodyCorrection
  LOGICAL  :: NNS_ApplyWeakMagnetism
  LOGICAL  :: Include_NES
  LOGICAL  :: Include_Pair
  LOGICAL  :: Include_NuPair
  LOGICAL  :: Include_Brem
  LOGICAL  :: Include_LinCorr
  LOGICAL  :: UseFirstMomentRates
  LOGICAL  :: UseConservativeProjection
  LOGICAL  :: UseRadiationOnlySolve
  INTEGER  :: M_outer
  INTEGER  :: M_inner
  INTEGER  :: MaxIter_outer
  INTEGER  :: MaxIter_inner
  INTEGER  :: M_FP
  REAL(DP) :: Rtol_outer
  REAL(DP) :: Rtol_inner
  REAL(DP) :: Atol_inner
  REAL(DP) :: wMatrRHS(nMatterEquations)
  REAL(DP) :: DnuMax
  LOGICAL  :: FreezeOpacities

  ! --- Temporary arrays for scatter/gather (packing)

  INTEGER,  DIMENSION(:)    , ALLOCATABLE, TARGET :: Error_T
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: D_T, T_T, Y_T, E_T
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: T0_T, Y0_T
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: SqrtGm_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Dnu_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Inu_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Inu_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Inu_u_3_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Dnu_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Sigma_Iso_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_0_Iso_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_1_Iso_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_EmAb_T, Eta_EmAb_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NNS_T , Eta_NNS_T, R_NNS_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NES_T , Eta_NES_T, R_NES_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Pair_T, Eta_Pair_T, R_Pair_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NuPair_T, Eta_NuPair_T, R_NuPair_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Brem_T, Eta_Brem_T, R_Brem_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__In__u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__In__u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__In__u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__Out_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__Out_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NNS__Out_u_3_T
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
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Pro_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Pro_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Pro_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Ann_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Ann_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_NuPair_Ann_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Pro_u_3_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_2_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: L_Brem_Ann_u_3_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_NNS_0_T, Phi_NbNS_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Phi_NNS_1_T, Phi_NbNS_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_0_T, H_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_1_T, H_II_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_0_T, J_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_1_T, J_II_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Nu_J_I_0_T, Nu_J_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Nu_J_I_1_T, Nu_J_II_1_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: S_Sigma_T

  LOGICAL,  DIMENSION(:,:),   ALLOCATABLE, TARGET :: Opacity_Mask_T

CONTAINS


  SUBROUTINE InitializeNeutrinoMatterSolver &
    ( iZ_B, iZ_E, ComputeFirstMomentRates_Option )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)
    LOGICAL, INTENT(in), OPTIONAL :: ComputeFirstMomentRates_Option

    INTEGER :: iE1, iE2, iN_E, iS, iC, iN_X

    REAL(DP) :: TMP(1)
    REAL(DP), ALLOCATABLE :: AMAT(:,:), BVEC(:)

    iE_B = iZ_B(1)
    iE_E = iZ_E(1)
    nZ   = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    n_FP_outer = nMatterEquations
    n_FP_inner = nE_G * nCR * nSpecies

    UseFirstMomentRates = Include_LinCorr
    IF ( PRESENT( ComputeFirstMomentRates_Option ) ) &
      UseFirstMomentRates &
        = UseFirstMomentRates .OR. ComputeFirstMomentRates_Option

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
    ALLOCATE( Fnorm0 (nSpecies,nX_G) )

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
    ALLOCATE(      Sigma_Iso(nE_G,nChirals,nX_G) )
    ALLOCATE(      Phi_0_Iso(nE_G,nChirals,nX_G) )
    ALLOCATE(      Phi_1_Iso(nE_G,nChirals,nX_G) )
    ALLOCATE(       Chi_EmAb(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_EmAb(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Chi_NNS(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Eta_NNS(nE_G,nSpecies,nX_G) )
    ALLOCATE(          R_NNS(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Chi_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Eta_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(          R_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(         R_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(     Chi_NuPair(nE_G,nSpecies,nX_G) )
    ALLOCATE(     Eta_NuPair(nE_G,nSpecies,nX_G) )
    ALLOCATE(       R_NuPair(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Brem(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Brem(nE_G,nSpecies,nX_G) )
    ALLOCATE(         R_Brem(nE_G,nSpecies,nX_G) )
    ALLOCATE(         Phi_NNS_0 (nE_G,nE_G,nX_G) )
    ALLOCATE(        Phi_NbNS_0 (nE_G,nE_G,nX_G) )
    ALLOCATE(             H_I_0 (nE_G,nE_G,nX_G) )
    ALLOCATE(             H_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE(             J_I_0 (nE_G,nE_G,nX_G) )
    ALLOCATE(             J_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE(           Nu_J_I_0(nE_G,nE_G,nX_G) )
    ALLOCATE(          Nu_J_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE(   S_Sigma(nE_G,nE_G,nX_G) )

    IF ( UseFirstMomentRates ) THEN
      ALLOCATE( L_NNS__In__u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__In__u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__In__u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__Out_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__Out_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__Out_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( Phi_NNS_1 (nE_G,nE_G,nX_G) )
      ALLOCATE( Phi_NbNS_1(nE_G,nE_G,nX_G) )
      ALLOCATE( L_NES__In__u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__In__u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__In__u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__Out_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__Out_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__Out_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( H_I_1 (nE_G,nE_G,nX_G) )
      ALLOCATE( H_II_1(nE_G,nE_G,nX_G) )
      ALLOCATE( L_Pair_Pro_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Pro_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Pro_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Ann_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Ann_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Ann_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( J_I_1 (nE_G,nE_G,nX_G) )
      ALLOCATE( J_II_1(nE_G,nE_G,nX_G) )
      ALLOCATE( L_NuPair_Pro_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Pro_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Pro_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Ann_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Ann_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Ann_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( Nu_J_I_1 (nE_G,nE_G,nX_G) )
      ALLOCATE( Nu_J_II_1(nE_G,nE_G,nX_G) )
      ALLOCATE( L_Brem_Pro_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Pro_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Pro_u_3(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Ann_u_1(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Ann_u_2(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Ann_u_3(nE_G,nSpecies,nX_G) )
    END IF

    ALLOCATE(    Opacity_Mask(nOp,nX_G) )

    ALLOCATE( Error_T(nX_G) )

    ALLOCATE( D_T(nX_G) )
    ALLOCATE( T_T(nX_G) )
    ALLOCATE( Y_T(nX_G) )
    ALLOCATE( E_T(nX_G) )

    ALLOCATE( T0_T(nX_G) )
    ALLOCATE( Y0_T(nX_G) )

    ALLOCATE( SqrtGm_T(nX_G) )

    ALLOCATE(     Dnu_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_1_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_2_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Inu_u_3_T(nE_G,nSpecies,nX_G) )

    ALLOCATE(          Dnu_0_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(      Sigma_Iso_T(nE_G,nChirals,nX_G) )
    ALLOCATE(      Phi_0_Iso_T(nE_G,nChirals,nX_G) )
    ALLOCATE(      Phi_1_Iso_T(nE_G,nChirals,nX_G) )
    ALLOCATE(       Chi_EmAb_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_EmAb_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Chi_NNS_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Eta_NNS_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(          R_NNS_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Chi_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(        Eta_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(          R_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(         R_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(     Chi_NuPair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(     Eta_NuPair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       R_NuPair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Chi_Brem_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(       Eta_Brem_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(         R_Brem_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(          Phi_NNS_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(         Phi_NbNS_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(              H_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(             H_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(              J_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(             J_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(           Nu_J_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(          Nu_J_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(   S_Sigma_T(nE_G,nE_G,nX_G) )

    IF ( UseFirstMomentRates ) THEN
      ALLOCATE( L_NNS__In__u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__In__u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__In__u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__Out_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__Out_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NNS__Out_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( Phi_NNS_1_T (nE_G,nE_G,nX_G) )
      ALLOCATE( Phi_NbNS_1_T(nE_G,nE_G,nX_G) )
      ALLOCATE( L_NES__In__u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__In__u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__In__u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__Out_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__Out_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NES__Out_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( H_I_1_T (nE_G,nE_G,nX_G) )
      ALLOCATE( H_II_1_T(nE_G,nE_G,nX_G) )
      ALLOCATE( L_Pair_Pro_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Pro_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Pro_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Ann_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Ann_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Pair_Ann_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( J_I_1_T (nE_G,nE_G,nX_G) )
      ALLOCATE( J_II_1_T(nE_G,nE_G,nX_G) )
      ALLOCATE( L_NuPair_Pro_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Pro_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Pro_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Ann_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Ann_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_NuPair_Ann_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( Nu_J_I_1_T (nE_G,nE_G,nX_G) )
      ALLOCATE( Nu_J_II_1_T(nE_G,nE_G,nX_G) )
      ALLOCATE( L_Brem_Pro_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Pro_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Pro_u_3_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Ann_u_1_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Ann_u_2_T(nE_G,nSpecies,nX_G) )
      ALLOCATE( L_Brem_Ann_u_3_T(nE_G,nSpecies,nX_G) )
    END IF
 
    ALLOCATE(    Opacity_Mask_T(nOp,nX_G) )

    ALLOCATE( INFO(nX_G) )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS ) &
    !$OMP MAP( alloc: INFO, &
    !$OMP             DnuNorm, &
    !$OMP             Fnorm0, &
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
    !$OMP             Chi_NNS, Eta_NNS, R_NNS, &
    !$OMP             Chi_NES, Eta_NES, R_NES, &
    !$OMP             Chi_Pair, Eta_Pair, R_Pair, &
    !$OMP             Chi_NuPair, Eta_NuPair, R_NuPair, &
    !$OMP             Chi_Brem, Eta_Brem, R_Brem, &
    !$OMP             Phi_NNS_0, Phi_NbNS_0, &
    !$OMP             H_I_0, H_II_0, J_I_0, J_II_0, &
    !$OMP             Nu_J_I_0, Nu_J_II_0, S_Sigma, &
    !$OMP             D_T, T_T, Y_T, E_T, Error_T,  &
    !$OMP             T0_T, Y0_T, &
    !$OMP             SqrtGm_T, &
    !$OMP             Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$OMP             Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$OMP             Chi_EmAb_T, Eta_EmAb_T, &
    !$OMP             Chi_NNS_T, Eta_NNS_T, R_NNS_T, &
    !$OMP             Chi_NES_T, Eta_NES_T, R_NES_T, &
    !$OMP             Chi_Pair_T, Eta_Pair_T, R_Pair_T, &
    !$OMP             Chi_NuPair_T, Eta_NuPair_T, R_NuPair_T, &
    !$OMP             Chi_Brem_T, Eta_Brem_T, R_Brem_T, &
    !$OMP             Phi_NNS_0_T, Phi_NbNS_0_T, &
    !$OMP             H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$OMP             Nu_J_I_0_T, Nu_J_II_0_T, S_Sigma_T, &
    !$OMP             Opacity_Mask, Opacity_Mask_T )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS ) &
    !$ACC CREATE( INFO, &
    !$ACC         DnuNorm, &
    !$ACC         Fnorm0, &
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
    !$ACC         Chi_NNS, Eta_NNS, R_NNS, &
    !$ACC         Chi_NES, Eta_NES, R_NES, &
    !$ACC         Chi_Pair, Eta_Pair, R_Pair, &
    !$ACC         Chi_NuPair, Eta_NuPair, R_NuPair, &
    !$ACC         Chi_Brem, Eta_Brem, R_Brem, &
    !$ACC         Phi_NNS_0, Phi_NbNS_0, &
    !$ACC         H_I_0, H_II_0, J_I_0, J_II_0, &
    !$ACC         Nu_J_I_0, Nu_J_II_0, S_Sigma, &
    !$ACC         D_T, T_T, Y_T, E_T, Error_T, &
    !$ACC         T0_T, Y0_T, &
    !$ACC         SqrtGm_T, &
    !$ACC         Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$ACC         Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$ACC         Chi_EmAb_T, Eta_EmAb_T, &
    !$ACC         Chi_NNS_T, Eta_NNS_T, R_NNS_T, &
    !$ACC         Chi_NES_T, Eta_NES_T, R_NES_T, &
    !$ACC         Chi_Pair_T, Eta_Pair_T, R_Pair_T, &
    !$ACC         Chi_NuPair_T, Eta_NuPair_T, R_NuPair_T, &
    !$ACC         Chi_Brem_T, Eta_Brem_T, R_Brem_T, &
    !$ACC         Phi_NNS_0_T, Phi_NbNS_0_T, &
    !$ACC         H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$ACC         Nu_J_I_0_T, Nu_J_II_0_T, S_Sigma_T, &
    !$ACC         Opacity_Mask, Opacity_Mask_T )
#endif

    IF ( UseFirstMomentRates ) THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: &
    !$OMP   L_NNS__In__u_1, L_NNS__In__u_2, L_NNS__In__u_3, &
    !$OMP   L_NNS__Out_u_1, L_NNS__Out_u_2, L_NNS__Out_u_3, &
    !$OMP   Phi_NNS_1, Phi_NbNS_1, &
    !$OMP   L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$OMP   L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$OMP   H_I_1, H_II_1, &
    !$OMP   L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$OMP   L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$OMP   J_I_1, J_II_1, &
    !$OMP   L_NuPair_Pro_u_1, L_NuPair_Pro_u_2, L_NuPair_Pro_u_3, &
    !$OMP   L_NuPair_Ann_u_1, L_NuPair_Ann_u_2, L_NuPair_Ann_u_3, &
    !$OMP   Nu_J_I_1, Nu_J_II_1, &
    !$OMP   L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$OMP   L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$OMP   L_NNS__In__u_1_T, L_NNS__In__u_2_T, L_NNS__In__u_3_T, &
    !$OMP   L_NNS__Out_u_1_T, L_NNS__Out_u_2_T, L_NNS__Out_u_3_T, &
    !$OMP   Phi_NNS_1_T, Phi_NbNS_1_T, &
    !$OMP   L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$OMP   L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$OMP   H_I_1_T, H_II_1_T, &
    !$OMP   L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$OMP   L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$OMP   J_I_1_T, J_II_1_T, &
    !$OMP   L_NuPair_Pro_u_1_T, L_NuPair_Pro_u_2_T, L_NuPair_Pro_u_3_T, &
    !$OMP   L_NuPair_Ann_u_1_T, L_NuPair_Ann_u_2_T, L_NuPair_Ann_u_3_T, &
    !$OMP   Nu_J_I_1_T, Nu_J_II_1_T, &
    !$OMP   L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$OMP   L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( &
    !$ACC   L_NNS__In__u_1, L_NNS__In__u_2, L_NNS__In__u_3, &
    !$ACC   L_NNS__Out_u_1, L_NNS__Out_u_2, L_NNS__Out_u_3, &
    !$ACC   Phi_NNS_1, Phi_NbNS_1, &
    !$ACC   L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$ACC   L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$ACC   H_I_1, H_II_1, &
    !$ACC   L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$ACC   L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$ACC   J_I_1, J_II_1, &
    !$ACC   L_NuPair_Pro_u_1, L_NuPair_Pro_u_2, L_NuPair_Pro_u_3, &
    !$ACC   L_NuPair_Ann_u_1, L_NuPair_Ann_u_2, L_NuPair_Ann_u_3, &
    !$ACC   Nu_J_I_1, Nu_J_II_1, &
    !$ACC   L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$ACC   L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$ACC   L_NNS__In__u_1_T, L_NNS__In__u_2_T, L_NNS__In__u_3_T, &
    !$ACC   L_NNS__Out_u_1_T, L_NNS__Out_u_2_T, L_NNS__Out_u_3_T, &
    !$ACC   Phi_NNS_1_T, Phi_NbNS_1_T, &
    !$ACC   L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$ACC   L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$ACC   H_I_1_T, H_II_1_T, &
    !$ACC   L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$ACC   L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$ACC   J_I_1_T, J_II_1_T, &
    !$ACC   L_NuPair_Pro_u_1_T, L_NuPair_Pro_u_2_T, L_NuPair_Pro_u_3_T, &
    !$ACC   L_NuPair_Ann_u_1_T, L_NuPair_Ann_u_2_T, L_NuPair_Ann_u_3_T, &
    !$ACC   Nu_J_I_1_T, Nu_J_II_1_T, &
    !$ACC   L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$ACC   L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T )
#endif
    END IF


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
!NNS
        Phi_NNS_0    (iE1,iE2,iN_X) = Zero
        Phi_NbNS_0   (iE1,iE2,iN_X) = Zero
        Phi_NNS_0_T  (iE1,iE2,iN_X) = Zero
        Phi_NbNS_0_T (iE1,iE2,iN_X) = Zero
!NES
        H_I_0   (iE1,iE2,iN_X) = Zero
        H_II_0  (iE1,iE2,iN_X) = Zero
        H_I_0_T (iE1,iE2,iN_X) = Zero
        H_II_0_T(iE1,iE2,iN_X) = Zero
!Pair
        J_I_0   (iE1,iE2,iN_X) = Zero
        J_II_0  (iE1,iE2,iN_X) = Zero
        J_I_0_T (iE1,iE2,iN_X) = Zero
        J_II_0_T(iE1,iE2,iN_X) = Zero
!NuPair
        Nu_J_I_0   (iE1,iE2,iN_X) = Zero
        Nu_J_II_0  (iE1,iE2,iN_X) = Zero
        Nu_J_I_0_T (iE1,iE2,iN_X) = Zero
        Nu_J_II_0_T(iE1,iE2,iN_X) = Zero

!Brem
        S_Sigma  (iE1,iE2,iN_X) = Zero
        S_Sigma_T(iE1,iE2,iN_X) = Zero

      END DO
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
      DO iC   = 1, nChirals
      DO iN_E = 1, nE_G

        Phi_1_Iso  (iN_E,iC,iN_X) = Zero
        Phi_1_Iso_T(iN_E,iC,iN_X) = Zero

      END DO
      END DO
      END DO

    IF ( UseFirstMomentRates ) THEN

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

        Phi_NNS_1   (iE1,iE2,iN_X) = Zero
        Phi_NbNS_1  (iE1,iE2,iN_X) = Zero
        Phi_NNS_1_T (iE1,iE2,iN_X) = Zero
        Phi_NbNS_1_T(iE1,iE2,iN_X) = Zero
        H_I_1       (iE1,iE2,iN_X) = Zero
        H_II_1      (iE1,iE2,iN_X) = Zero
        H_I_1_T     (iE1,iE2,iN_X) = Zero
        H_II_1_T    (iE1,iE2,iN_X) = Zero
        J_I_1       (iE1,iE2,iN_X) = Zero
        J_II_1      (iE1,iE2,iN_X) = Zero
        J_I_1_T     (iE1,iE2,iN_X) = Zero
        J_II_1_T    (iE1,iE2,iN_X) = Zero
        Nu_J_I_1    (iE1,iE2,iN_X) = Zero
        Nu_J_II_1   (iE1,iE2,iN_X) = Zero
        Nu_J_I_1_T  (iE1,iE2,iN_X) = Zero
        Nu_J_II_1_T (iE1,iE2,iN_X) = Zero

      END DO
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

        L_NNS__In__u_1  (iN_E,iS,iN_X) = Zero
        L_NNS__In__u_2  (iN_E,iS,iN_X) = Zero
        L_NNS__In__u_3  (iN_E,iS,iN_X) = Zero
        L_NNS__Out_u_1  (iN_E,iS,iN_X) = Zero
        L_NNS__Out_u_2  (iN_E,iS,iN_X) = Zero
        L_NNS__Out_u_3  (iN_E,iS,iN_X) = Zero
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
        L_NuPair_Pro_u_1(iN_E,iS,iN_X) = Zero
        L_NuPair_Pro_u_2(iN_E,iS,iN_X) = Zero
        L_NuPair_Pro_u_3(iN_E,iS,iN_X) = Zero
        L_NuPair_Ann_u_1(iN_E,iS,iN_X) = Zero
        L_NuPair_Ann_u_2(iN_E,iS,iN_X) = Zero
        L_NuPair_Ann_u_3(iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_1  (iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_2  (iN_E,iS,iN_X) = Zero
        L_Brem_Pro_u_3  (iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_1  (iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_2  (iN_E,iS,iN_X) = Zero
        L_Brem_Ann_u_3  (iN_E,iS,iN_X) = Zero

        L_NNS__In__u_1_T (iN_E,iS,iN_X) = Zero
        L_NNS__In__u_2_T (iN_E,iS,iN_X) = Zero
        L_NNS__In__u_3_T (iN_E,iS,iN_X) = Zero
        L_NNS__Out_u_1_T (iN_E,iS,iN_X) = Zero
        L_NNS__Out_u_2_T (iN_E,iS,iN_X) = Zero
        L_NNS__Out_u_3_T (iN_E,iS,iN_X) = Zero
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
        L_NuPair_Pro_u_1_T(iN_E,iS,iN_X) = Zero
        L_NuPair_Pro_u_2_T(iN_E,iS,iN_X) = Zero
        L_NuPair_Pro_u_3_T(iN_E,iS,iN_X) = Zero
        L_NuPair_Ann_u_1_T(iN_E,iS,iN_X) = Zero
        L_NuPair_Ann_u_2_T(iN_E,iS,iN_X) = Zero
        L_NuPair_Ann_u_3_T(iN_E,iS,iN_X) = Zero
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
      Atol_inner_Option, Include_NNS_Option, NNS_ApplyWeakMagnetism_Option, &
      NNS_ApplyManyBodyCorrection_Option, Include_NES_Option, Include_Pair_Option, &
      Include_NuPair_Option, Include_Brem_Option, Include_LinCorr_Option, &
      wMatrRHS_Option, DnuMax_Option, FreezeOpacities_Option, &
      UseConservativeProjection_Option, Verbose_Option, &
      UseRadiationOnlySolve_Option )

    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE mpi

    INTEGER , INTENT(in), OPTIONAL :: M_outer_Option
    INTEGER , INTENT(in), OPTIONAL :: M_inner_Option
    INTEGER , INTENT(in), OPTIONAL :: MaxIter_outer_Option
    INTEGER , INTENT(in), OPTIONAL :: MaxIter_inner_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rtol_outer_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rtol_inner_Option
    REAL(DP), INTENT(in), OPTIONAL :: Atol_inner_option
    LOGICAL , INTENT(in), OPTIONAL :: Include_NNS_Option
    LOGICAL , INTENT(in), OPTIONAL :: NNS_ApplyManyBodyCorrection_Option
    LOGICAL , INTENT(in), OPTIONAL :: NNS_ApplyWeakMagnetism_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_NES_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_Pair_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_NuPair_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_Brem_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_LinCorr_Option
    REAL(DP), INTENT(in), OPTIONAL :: wMatrRHS_Option(nMatterEquations)
    REAL(DP), INTENT(in), OPTIONAL :: DnuMax_Option
    LOGICAL , INTENT(in), OPTIONAL :: FreezeOpacities_Option
    LOGICAL , INTENT(in), OPTIONAL :: UseConservativeProjection_Option
    LOGICAL , INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL , INTENT(in), OPTIONAL :: UseRadiationOnlySolve_Option

    LOGICAL :: Verbose
    INTEGER :: ierr, myid

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

    IF( PRESENT( Rtol_outer_Option ) )THEN
      Rtol_outer = Rtol_outer_Option
    ELSE
      Rtol_outer = 1.0d-08
    END IF

    IF( PRESENT( Rtol_inner_Option ) )THEN
      Rtol_inner = Rtol_inner_Option
    ELSE
      Rtol_inner = 1.0d-08
    END IF

    IF( PRESENT( Atol_inner_Option ) )THEN
      Atol_inner = Atol_inner_Option
    ELSE
      Atol_inner = 1.0d-02 * Rtol_inner
    END IF

    IF( PRESENT( Include_NNS_Option ) )THEN
      Include_NNS = Include_NNS_Option
    ELSE
      Include_NNS = .TRUE.
    END IF

    IF ( Include_NNS .AND. IsoIncludesNucleonScattering ) THEN

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

      IF ( myid == 0 ) THEN

        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: Include_NNS is enabled, but the Iso table already'
        WRITE(ERROR_UNIT,'(A)') &
          '       contains isoenergetic scattering on free nucleons.'
        WRITE(ERROR_UNIT,'(A)') &
          '       Use an Iso table without scattering on n/p.'

        FLUSH(ERROR_UNIT)

      END IF

      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

    END IF

    IF( PRESENT( NNS_ApplyManyBodyCorrection_Option ) )THEN
      NNS_ApplyManyBodyCorrection = NNS_ApplyManyBodyCorrection_Option
    ELSE
      NNS_ApplyManyBodyCorrection = .FALSE.
    END IF

    IF( PRESENT( NNS_ApplyWeakMagnetism_Option ) )THEN
      NNS_ApplyWeakMagnetism = NNS_ApplyWeakMagnetism_Option
    ELSE
      NNS_ApplyWeakMagnetism = .TRUE.
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

    IF( PRESENT( Include_NuPair_Option ) )THEN
      Include_NuPair = Include_NuPair_Option
    ELSE
      Include_NuPair = .FALSE.
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

    IF( PRESENT( UseConservativeProjection_Option ) )THEN
      UseConservativeProjection = UseConservativeProjection_Option
    ELSE
      UseConservativeProjection = .FALSE.
    END IF

    IF( PRESENT( UseRadiationOnlySolve_Option ) )THEN
      UseRadiationOnlySolve = UseRadiationOnlySolve_Option
    ELSE
      UseRadiationOnlySolve = .FALSE.
    END IF

    IF( UseRadiationOnlySolve .AND. UseConservativeProjection )THEN

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

      IF( myid == 0 )THEN
        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: The invariant-constrained solve already enforces the'
        WRITE(ERROR_UNIT,'(A)') &
          '       conserved matter state. Disable UseConservativeProjection.'
        FLUSH(ERROR_UNIT)
      END IF

      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

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
      WRITE(*,'(A4,A32,ES10.3E3)') '', 'Atol_inner: '     , Atol_inner
      WRITE(*,'(A4,A32,L1)')       '', 'Conservative projection: ', &
        UseConservativeProjection
      WRITE(*,'(A4,A32,L1)')       '', 'Invariant-constrained solve: ', &
        UseRadiationOnlySolve
      WRITE(*,*)
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NNS: '    , Include_NNS
      WRITE(*,'(A4,A32,L1)')       '', 'NNS_manybody: '   , NNS_ApplyManyBodyCorrection
      WRITE(*,'(A4,A32,L1)')       '', 'NNS_weak_mag: '   , NNS_ApplyWeakMagnetism
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NES: '    , Include_NES
      WRITE(*,'(A4,A32,L1)')       '', 'Include_Pair: '   , Include_Pair
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NuPair: ' , Include_NuPair
      WRITE(*,'(A4,A32,L1)')       '', 'Include_Brem: '   , Include_Brem
      WRITE(*,'(A4,A32,L1)')       '', 'Include_LinCorr: ', Include_LinCorr
      WRITE(*,*)
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iY ): '  , INT(wMatrRHS(iY ))
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iEf): '  , INT(wMatrRHS(iEf))
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iV1): '  , INT(wMatrRHS(iV1))
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iV2): '  , INT(wMatrRHS(iV2))
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iV3): '  , INT(wMatrRHS(iV3))
      WRITE(*,*)
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
    !$OMP               Fnorm0, &
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
    !$OMP               Chi_NNS, Eta_NNS, R_NNS, &
    !$OMP               Chi_NES, Eta_NES, R_NES, &
    !$OMP               Chi_Pair, Eta_Pair, R_Pair, &
    !$OMP               Chi_NuPair, Eta_NuPair, R_NuPair, &
    !$OMP               Chi_Brem, Eta_Brem, R_Brem, &
    !$OMP               Phi_NNS_0, Phi_NbNS_0, &
    !$OMP               H_I_0, H_II_0, J_I_0, J_II_0, &
    !$OMP               Nu_J_I_0, Nu_J_II_0, S_Sigma, &
    !$OMP               D_T, T_T, Y_T, E_T, Error_T,  &
    !$OMP               T0_T, Y0_T, &
    !$OMP               SqrtGm_T, &
    !$OMP               Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$OMP               Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$OMP               Chi_EmAb_T, Eta_EmAb_T, &
    !$OMP               Chi_NNS_T, Eta_NNS_T, R_NNS_T, &
    !$OMP               Chi_NES_T, Eta_NES_T, R_NES_T, &
    !$OMP               Chi_Pair_T, Eta_Pair_T, R_Pair_T, &
    !$OMP               Chi_NuPair_T, Eta_NuPair_T, R_NuPair_T, &
    !$OMP               Chi_Brem_T, Eta_Brem_T, R_Brem_T, &
    !$OMP               Phi_NNS_0_T, Phi_NbNS_0_T, &
    !$OMP               H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$OMP               Nu_J_I_0_T, Nu_J_II_0_T, S_Sigma_T, &
    !$OMP               Opacity_Mask, Opacity_Mask_T )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatrRHS, &
    !$ACC         INFO, &
    !$ACC         DnuNorm, &
    !$ACC         Fnorm0, &
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
    !$ACC         Chi_NNS, Eta_NNS, R_NNS, &
    !$ACC         Chi_NES, Eta_NES, R_NES, &
    !$ACC         Chi_Pair, Eta_Pair, R_Pair, &
    !$ACC         Chi_NuPair, Eta_NuPair, R_NuPair, &
    !$ACC         Chi_Brem, Eta_Brem, R_Brem, &
    !$ACC         Phi_NNS_0, Phi_NbNS_0, &
    !$ACC         H_I_0, H_II_0, J_I_0, J_II_0, &
    !$ACC         Nu_J_I_0, Nu_J_II_0, S_Sigma, &
    !$ACC         D_T, T_T, Y_T, E_T, Error_T,  &
    !$ACC         T0_T, Y0_T, &
    !$ACC         SqrtGm_T, &
    !$ACC         Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T, &
    !$ACC         Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T, &
    !$ACC         Chi_EmAb_T, Eta_EmAb_T, &
    !$ACC         Chi_NNS_T, Eta_NNS_T, R_NNS_T, &
    !$ACC         Chi_NES_T, Eta_NES_T, R_NES_T, &
    !$ACC         Chi_Pair_T, Eta_Pair_T, R_Pair_T, &
    !$ACC         Chi_NuPair_T, Eta_NuPair_T, R_NuPair_T, &
    !$ACC         Chi_Brem_T, Eta_Brem_T, R_Brem_T, &
    !$ACC         Phi_NNS_0_T, Phi_NbNS_0_T, &
    !$ACC         H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, &
    !$ACC         Nu_J_I_0_T, Nu_J_II_0_T, S_Sigma_T, &
    !$ACC         Opacity_Mask, Opacity_Mask_T )
#endif

    IF ( UseFirstMomentRates ) THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: &
    !$OMP   L_NNS__In__u_1, L_NNS__In__u_2, L_NNS__In__u_3, &
    !$OMP   L_NNS__Out_u_1, L_NNS__Out_u_2, L_NNS__Out_u_3, &
    !$OMP   Phi_NNS_1, Phi_NbNS_1, &
    !$OMP   L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$OMP   L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$OMP   H_I_1, H_II_1, &
    !$OMP   L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$OMP   L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$OMP   J_I_1, J_II_1, &
    !$OMP   L_NuPair_Pro_u_1, L_NuPair_Pro_u_2, L_NuPair_Pro_u_3, &
    !$OMP   L_NuPair_Ann_u_1, L_NuPair_Ann_u_2, L_NuPair_Ann_u_3, &
    !$OMP   Nu_J_I_1, Nu_J_II_1, &
    !$OMP   L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$OMP   L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$OMP   L_NNS__In__u_1_T, L_NNS__In__u_2_T, L_NNS__In__u_3_T, &
    !$OMP   L_NNS__Out_u_1_T, L_NNS__Out_u_2_T, L_NNS__Out_u_3_T, &
    !$OMP   Phi_NNS_1_T, Phi_NbNS_1_T, &
    !$OMP   L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$OMP   L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$OMP   H_I_1_T, H_II_1_T, &
    !$OMP   L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$OMP   L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$OMP   J_I_1_T, J_II_1_T, &
    !$OMP   L_NuPair_Pro_u_1_T, L_NuPair_Pro_u_2_T, L_NuPair_Pro_u_3_T, &
    !$OMP   L_NuPair_Ann_u_1_T, L_NuPair_Ann_u_2_T, L_NuPair_Ann_u_3_T, &
    !$OMP   Nu_J_I_1_T, Nu_J_II_1_T, &
    !$OMP   L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$OMP   L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( &
    !$ACC   L_NNS__In__u_1, L_NNS__In__u_2, L_NNS__In__u_3, &
    !$ACC   L_NNS__Out_u_1, L_NNS__Out_u_2, L_NNS__Out_u_3, &
    !$ACC   Phi_NNS_1, Phi_NbNS_1, &
    !$ACC   L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3, &
    !$ACC   L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$ACC   H_I_1, H_II_1, &
    !$ACC   L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3, &
    !$ACC   L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$ACC   J_I_1, J_II_1, &
    !$ACC   L_NuPair_Pro_u_1, L_NuPair_Pro_u_2, L_NuPair_Pro_u_3, &
    !$ACC   L_NuPair_Ann_u_1, L_NuPair_Ann_u_2, L_NuPair_Ann_u_3, &
    !$ACC   Nu_J_I_1, Nu_J_II_1, &
    !$ACC   L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3, &
    !$ACC   L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$ACC   L_NNS__In__u_1_T, L_NNS__In__u_2_T, L_NNS__In__u_3_T, &
    !$ACC   L_NNS__Out_u_1_T, L_NNS__Out_u_2_T, L_NNS__Out_u_3_T, &
    !$ACC   Phi_NNS_1_T, Phi_NbNS_1_T, &
    !$ACC   L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T, &
    !$ACC   L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T, &
    !$ACC   H_I_1_T, H_II_1_T, &
    !$ACC   L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T, &
    !$ACC   L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T, &
    !$ACC   J_I_1_T, J_II_1_T, &
    !$ACC   L_NuPair_Pro_u_1_T, L_NuPair_Pro_u_2_T, L_NuPair_Pro_u_3_T, &
    !$ACC   L_NuPair_Ann_u_1_T, L_NuPair_Ann_u_2_T, L_NuPair_Ann_u_3_T, &
    !$ACC   Nu_J_I_1_T, Nu_J_II_1_T, &
    !$ACC   L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T, &
    !$ACC   L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T )
#endif
    END IF

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2 )
    DEALLOCATE( INFO )
    DEALLOCATE( DnuNorm )
    DEALLOCATE( Fnorm0 )
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
    DEALLOCATE( Chi_NNS, Eta_NNS, R_NNS )
    DEALLOCATE( Chi_NES, Eta_NES, R_NES )
    DEALLOCATE( Chi_Pair, Eta_Pair, R_Pair )
    DEALLOCATE( Chi_NuPair, Eta_NuPair, R_NuPair )
    DEALLOCATE( Chi_Brem, Eta_Brem, R_Brem )
    DEALLOCATE( Phi_NNS_0, Phi_NbNS_0 )
    DEALLOCATE( H_I_0, H_II_0, J_I_0, J_II_0 )
    DEALLOCATE( S_Sigma )
    DEALLOCATE( Nu_J_I_0, Nu_J_II_0 )
    DEALLOCATE( D_T, T_T, Y_T, E_T, SqrtGm_T, Error_T )
    DEALLOCATE( T0_T, Y0_T )
    DEALLOCATE( Dnu_T, Inu_u_1_T, Inu_u_2_T, Inu_u_3_T )
    DEALLOCATE( Dnu_0_T, Sigma_Iso_T, Phi_0_Iso_T, Phi_1_Iso_T )
    DEALLOCATE( Chi_EmAb_T, Eta_EmAb_T )
    DEALLOCATE( Chi_NNS_T, Eta_NNS_T, R_NNS_T )
    DEALLOCATE( Chi_NES_T, Eta_NES_T, R_NES_T )
    DEALLOCATE( Chi_Pair_T, Eta_Pair_T, R_Pair_T )
    DEALLOCATE( Chi_NuPair_T, Eta_NuPair_T, R_NuPair_T )
    DEALLOCATE( Chi_Brem_T, Eta_Brem_T, R_Brem_T )
    DEALLOCATE( Phi_NNS_0_T, Phi_NbNS_0_T )
    DEALLOCATE( H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T )
    DEALLOCATE( S_Sigma_T )
    DEALLOCATE( Nu_J_I_0_T, Nu_J_II_0_T )

    IF ( UseFirstMomentRates ) THEN
      DEALLOCATE( L_NNS__In__u_1, L_NNS__In__u_2, L_NNS__In__u_3 )
      DEALLOCATE( L_NNS__Out_u_1, L_NNS__Out_u_2, L_NNS__Out_u_3 )
      DEALLOCATE( Phi_NNS_1, Phi_NbNS_1 )
      DEALLOCATE( L_NES__In__u_1, L_NES__In__u_2, L_NES__In__u_3 )
      DEALLOCATE( L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3 )
      DEALLOCATE( H_I_1, H_II_1 )
      DEALLOCATE( L_Pair_Pro_u_1, L_Pair_Pro_u_2, L_Pair_Pro_u_3 )
      DEALLOCATE( L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3 )
      DEALLOCATE( J_I_1, J_II_1 )
      DEALLOCATE( L_NuPair_Pro_u_1, L_NuPair_Pro_u_2, L_NuPair_Pro_u_3 )
      DEALLOCATE( L_NuPair_Ann_u_1, L_NuPair_Ann_u_2, L_NuPair_Ann_u_3 )
      DEALLOCATE( Nu_J_I_1, Nu_J_II_1 )
      DEALLOCATE( L_Brem_Pro_u_1, L_Brem_Pro_u_2, L_Brem_Pro_u_3 )
      DEALLOCATE( L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3 )
      DEALLOCATE( L_NNS__In__u_1_T, L_NNS__In__u_2_T, L_NNS__In__u_3_T )
      DEALLOCATE( L_NNS__Out_u_1_T, L_NNS__Out_u_2_T, L_NNS__Out_u_3_T )
      DEALLOCATE( Phi_NNS_1_T, Phi_NbNS_1_T )
      DEALLOCATE( L_NES__In__u_1_T, L_NES__In__u_2_T, L_NES__In__u_3_T )
      DEALLOCATE( L_NES__Out_u_1_T, L_NES__Out_u_2_T, L_NES__Out_u_3_T )
      DEALLOCATE( H_I_1_T, H_II_1_T )
      DEALLOCATE( L_Pair_Pro_u_1_T, L_Pair_Pro_u_2_T, L_Pair_Pro_u_3_T )
      DEALLOCATE( L_Pair_Ann_u_1_T, L_Pair_Ann_u_2_T, L_Pair_Ann_u_3_T )
      DEALLOCATE( J_I_1_T, J_II_1_T )
      DEALLOCATE( L_NuPair_Pro_u_1_T, L_NuPair_Pro_u_2_T, L_NuPair_Pro_u_3_T )
      DEALLOCATE( L_NuPair_Ann_u_1_T, L_NuPair_Ann_u_2_T, L_NuPair_Ann_u_3_T )
      DEALLOCATE( Nu_J_I_1_T, Nu_J_II_1_T )
      DEALLOCATE( L_Brem_Pro_u_1_T, L_Brem_Pro_u_2_T, L_Brem_Pro_u_3_T )
      DEALLOCATE( L_Brem_Ann_u_1_T, L_Brem_Ann_u_2_T, L_Brem_Ann_u_3_T )
    END IF
    DEALLOCATE( Opacity_Mask, Opacity_Mask_T )

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


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE ComputeTransportOpacity &
    ( D, T, Y, SqrtGm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      TransportOpacity, PatchHasCollisions_Option )

    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:), SqrtGm(:)
    REAL(DP), INTENT(in)  :: Dnu(:,:,:)
    REAL(DP), INTENT(in)  :: Inu_u_1(:,:,:), Inu_u_2(:,:,:), Inu_u_3(:,:,:)
    REAL(DP), INTENT(out) :: TransportOpacity(:,:,:)
    LOGICAL,  INTENT(out), OPTIONAL :: PatchHasCollisions_Option

    INTEGER :: ierr, myid
    LOGICAL :: PatchHasCollisions

    IF( Include_LinCorr )THEN

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

      IF( myid == 0 )THEN
        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: AP LLF scaling is incompatible with first-moment rates.'
        FLUSH(ERROR_UNIT)
      END IF

      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

    END IF

    CALL SetOpacityMask( nX_G, nDOFX, D, T, Y, Opacity_Mask )

    PatchHasCollisions &
      =     ANY( Opacity_Mask(iOp_EmAb,   :) ) &
        .OR. ANY( Opacity_Mask(iOp_ECTable,:) ) &
        .OR. ANY( Opacity_Mask(iOp_Iso,    :) ) &
        .OR. ( Include_NNS    .AND. ANY( Opacity_Mask(iOp_NNS,   :) ) ) &
        .OR. ( Include_NES    .AND. ANY( Opacity_Mask(iOp_NES,   :) ) ) &
        .OR. ( Include_Pair   .AND. ANY( Opacity_Mask(iOp_Pair,  :) ) ) &
        .OR. ( Include_NuPair .AND. ANY( Opacity_Mask(iOp_NuPair,:) ) ) &
        .OR. ( Include_Brem   .AND. ANY( Opacity_Mask(iOp_Brem,  :) ) )

    IF( .NOT. PatchHasCollisions )THEN

      IF( PRESENT( PatchHasCollisions_Option ) ) &
        PatchHasCollisions_Option = .FALSE.
      CALL ZeroTransportOpacity( TransportOpacity )
      RETURN

    END IF

    CALL ComputeOpacities_Packed &
           ( D, T, Y, SqrtGm, nX_P0 = nX_G, T0 = T, Y0 = Y, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL ComputeRates_Packed &
           ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, nX_P0 = nX_G, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL AssembleTransportOpacity( TransportOpacity )

    IF( PRESENT( PatchHasCollisions_Option ) ) &
      PatchHasCollisions_Option = .TRUE.

  END SUBROUTINE ComputeTransportOpacity


  SUBROUTINE ZeroTransportOpacity( TransportOpacity )

    REAL(DP), INTENT(out) :: TransportOpacity(:,:,:)

    INTEGER :: iN_E, iN_X, iS

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( TransportOpacity )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, SIZE( TransportOpacity, 3 )
    DO iS   = 1, SIZE( TransportOpacity, 2 )
    DO iN_E = 1, SIZE( TransportOpacity, 1 )

      TransportOpacity(iN_E,iS,iN_X) = Zero

    END DO
    END DO
    END DO

  END SUBROUTINE ZeroTransportOpacity


  SUBROUTINE AssembleTransportOpacity( TransportOpacity )

    REAL(DP), INTENT(out) :: TransportOpacity(:,:,:)

    INTEGER :: iChiral, iN_E, iN_X, iS
    LOGICAL :: UseNNS, UseNES, UsePair, UseNuPair, UseBrem

    UseNNS    = Include_NNS
    UseNES    = Include_NES
    UsePair   = Include_Pair
    UseNuPair = Include_NuPair
    UseBrem   = Include_Brem

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iChiral )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iChiral ) &
    !$ACC FIRSTPRIVATE( UseNNS, UseNES, UsePair, UseNuPair, UseBrem ) &
    !$ACC PRESENT( TransportOpacity, Sigma_Iso, Chi_EmAb, Chi_NNS, &
    !$ACC          Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE( iChiral )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      iChiral = nChirals - ( LeptonNumber(iS) + 1 ) / nChirals

      TransportOpacity(iN_E,iS,iN_X) &
        =   Sigma_Iso(iN_E,iChiral,iN_X) &
          + Chi_EmAb  (iN_E,iS,iN_X)

      IF( UseNNS ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) + Chi_NNS(iN_E,iS,iN_X)

      IF( UseNES ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) + Chi_NES(iN_E,iS,iN_X)

      IF( UsePair ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) + Chi_Pair(iN_E,iS,iN_X)

      IF( UseNuPair ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) + Chi_NuPair(iN_E,iS,iN_X)

      IF( UseBrem ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) + Chi_Brem(iN_E,iS,iN_X)

    END DO
    END DO
    END DO

  END SUBROUTINE AssembleTransportOpacity


  SUBROUTINE ComputeMeanInverseMeanFreePaths &
    ( D, T, Y, SqrtGm, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, MeanInvMFP )

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:), SqrtGm(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP), INTENT(in)  :: Dnu(:,:,:)
    REAL(DP), INTENT(in)  :: Inu_u_1(:,:,:), Inu_u_2(:,:,:), Inu_u_3(:,:,:)
    REAL(DP), INTENT(out) :: MeanInvMFP(:,:,:)

    INTEGER :: iChiral, iN_E, iN_X, iProcess, iS
    REAL(DP) :: Denominator, IdotI, IdotL, InvW2_Ref, WeightE
    REAL(DP) :: Numerator(nC)
    CALL SetOpacityMask( nX_G, nDOFX, D, T, Y, Opacity_Mask )

    ! This is an instantaneous diagnostic.  The final state is therefore
    ! both the evaluation state and the EC-table reference state.
    CALL ComputeOpacities_Packed &
           ( D, T, Y, SqrtGm, nX_P0 = nX_G, T0 = T, Y0 = Y, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL ComputeRates_Packed &
           ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, nX_P0 = nX_G, &
             Opacity_Mask_Input = Opacity_Mask )

    ! This diagnostic is a ratio of energy integrals.  Normalize the
    ! quadrature measure before accumulating it so its scale cancels exactly
    ! while the intermediate sums remain well away from underflow.
    InvW2_Ref = One / MAXVAL( W2_N )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iChiral, iN_E, iProcess, IdotI, IdotL, WeightE, &
    !$OMP          Denominator, Numerator )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iChiral, iN_E, iProcess, IdotI, IdotL, WeightE, &
    !$ACC          Denominator, Numerator ) &
    !$ACC PRESENT( MeanInvMFP, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33, Dnu_0, &
    !$ACC          Chi_EmAb, Phi_0_Iso, Phi_1_Iso, &
    !$ACC          Chi_NNS, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, &
    !$ACC          Eta_NNS, Eta_NES, Eta_Pair, Eta_NuPair, Eta_Brem, &
    !$ACC          L_NNS__Out_u_1, L_NNS__Out_u_2, L_NNS__Out_u_3, &
    !$ACC          L_NES__Out_u_1, L_NES__Out_u_2, L_NES__Out_u_3, &
    !$ACC          L_Pair_Ann_u_1, L_Pair_Ann_u_2, L_Pair_Ann_u_3, &
    !$ACC          L_NuPair_Ann_u_1, L_NuPair_Ann_u_2, L_NuPair_Ann_u_3, &
    !$ACC          L_Brem_Ann_u_1, L_Brem_Ann_u_2, L_Brem_Ann_u_3, &
    !$ACC          W2_N, FourPiEp2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iChiral, iN_E, iProcess, IdotI, IdotL, WeightE, &
    !$OMP          Denominator, Numerator )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies

      Denominator = Zero
      DO iProcess = 1, nC
        Numerator(iProcess) = Zero
      END DO

      iChiral = nChirals - ( LeptonNumber(iS) + 1 ) / nChirals

      DO iN_E = 1, nE_G

        WeightE = W2_N(iN_E) * InvW2_Ref

        Denominator = Denominator &
          + WeightE * Dnu(iN_E,iS,iN_X)

        IdotI = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)**2 &
              + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)**2 &
              + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)**2

        Numerator(iC_EmAb) = Numerator(iC_EmAb) &
          + WeightE * Dnu(iN_E,iS,iN_X) &
            * Chi_EmAb(iN_E,iS,iN_X) &
            * ( One - Dnu_0(iN_E,iS,iN_X) )

        Numerator(iC_Iso) = Numerator(iC_Iso) + WeightE &
          * FourPiEp2(iN_E) &
          * ( Dnu(iN_E,iS,iN_X) * Phi_0_Iso(iN_E,iChiral,iN_X) &
                * ( One - Dnu(iN_E,iS,iN_X) ) &
              - Phi_1_Iso(iN_E,iChiral,iN_X) * IdotI )

        IF ( Include_NNS ) THEN
          IdotL = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                    * L_NNS__Out_u_1(iN_E,iS,iN_X) &
                + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                    * L_NNS__Out_u_2(iN_E,iS,iN_X) &
                + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                    * L_NNS__Out_u_3(iN_E,iS,iN_X)
          Numerator(iC_NNS) = Numerator(iC_NNS) + WeightE &
            * ( Dnu(iN_E,iS,iN_X) &
                  * ( Chi_NNS(iN_E,iS,iN_X) - Eta_NNS(iN_E,iS,iN_X) ) &
                - IdotL )
        END IF

        IF ( Include_NES ) THEN
          IdotL = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                    * L_NES__Out_u_1(iN_E,iS,iN_X) &
                + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                    * L_NES__Out_u_2(iN_E,iS,iN_X) &
                + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                    * L_NES__Out_u_3(iN_E,iS,iN_X)
          Numerator(iC_NES) = Numerator(iC_NES) + WeightE &
            * ( Dnu(iN_E,iS,iN_X) &
                  * ( Chi_NES(iN_E,iS,iN_X) - Eta_NES(iN_E,iS,iN_X) ) &
                - IdotL )
        END IF

        IF ( Include_Pair ) THEN
          IdotL = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                    * L_Pair_Ann_u_1(iN_E,iS,iN_X) &
                + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                    * L_Pair_Ann_u_2(iN_E,iS,iN_X) &
                + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                    * L_Pair_Ann_u_3(iN_E,iS,iN_X)
          Numerator(iC_Pair) = Numerator(iC_Pair) + WeightE &
            * ( Dnu(iN_E,iS,iN_X) &
                  * ( Chi_Pair(iN_E,iS,iN_X) - Eta_Pair(iN_E,iS,iN_X) ) &
                + IdotL )
        END IF

        IF ( Include_NuPair ) THEN
          IdotL = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                    * L_NuPair_Ann_u_1(iN_E,iS,iN_X) &
                + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                    * L_NuPair_Ann_u_2(iN_E,iS,iN_X) &
                + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                    * L_NuPair_Ann_u_3(iN_E,iS,iN_X)
          Numerator(iC_NuPair) = Numerator(iC_NuPair) + WeightE &
            * ( Dnu(iN_E,iS,iN_X) &
                  * ( Chi_NuPair(iN_E,iS,iN_X) &
                      - Eta_NuPair(iN_E,iS,iN_X) ) &
                + IdotL )
        END IF

        IF ( Include_Brem ) THEN
          IdotL = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
                    * L_Brem_Ann_u_1(iN_E,iS,iN_X) &
                + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
                    * L_Brem_Ann_u_2(iN_E,iS,iN_X) &
                + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X) &
                    * L_Brem_Ann_u_3(iN_E,iS,iN_X)
          Numerator(iC_Brem) = Numerator(iC_Brem) + WeightE &
            * ( Dnu(iN_E,iS,iN_X) &
                  * ( Chi_Brem(iN_E,iS,iN_X) - Eta_Brem(iN_E,iS,iN_X) ) &
                + IdotL )
        END IF

      END DO

      IF ( Denominator > Zero ) THEN
        DO iProcess = 1, nC
          MeanInvMFP(iProcess,iS,iN_X) &
            = Numerator(iProcess) / Denominator
        END DO
      ELSE
        DO iProcess = 1, nC
          MeanInvMFP(iProcess,iS,iN_X) = Zero
        END DO
      END IF

    END DO
    END DO

  END SUBROUTINE ComputeMeanInverseMeanFreePaths
#endif


  SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA &
    ( dt, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, D, T, Y, E, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, nIterations_Inner, nIterations_Outer, &
      nIterations_Inner_Final, nIterations_Outer_Final, &
      RelativeResidual_Inner_First, &
      MatterCoupling_Mask_Option &
#if defined( TWOMOMENT_ORDER_V )
      , ComputeMatterEquationResidual, MatterEquationResidual_Rate_Final, &
        InitialRadiationLimiter_Rate_Final, &
        RadiationResidual_Relative_Final, &
        ConservativeProjectionRadiationResidual_Relative_Final, &
        ConservativeProjectionMatterHeating_Rate_Final, &
        ConservativeProjectionStatus_Final, &
        LeptonFlavorAudit_Rate_Final &
#endif
      )

    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE mpi

    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: D, T, Y, E
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Alpha
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Beta_u_1, Beta_u_2, Beta_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    INTEGER,  DIMENSION(:),     INTENT(inout) :: nIterations_Inner
    INTEGER,  DIMENSION(:),     INTENT(inout) :: nIterations_Outer
    INTEGER,  DIMENSION(:),     INTENT(out)   :: nIterations_Inner_Final
    INTEGER,  DIMENSION(:),     INTENT(out)   :: nIterations_Outer_Final
    REAL(DP), DIMENSION(:),     INTENT(out)   :: RelativeResidual_Inner_First
    LOGICAL,  DIMENSION(:),     INTENT(out), OPTIONAL :: MatterCoupling_Mask_Option
#if defined( TWOMOMENT_ORDER_V )
    LOGICAL,                    INTENT(in)    :: ComputeMatterEquationResidual
    REAL(DP), DIMENSION(:,:),   INTENT(out)   :: MatterEquationResidual_Rate_Final
    REAL(DP), DIMENSION(:,:),   INTENT(out)   :: InitialRadiationLimiter_Rate_Final
    REAL(DP), DIMENSION(:),     INTENT(out)   :: RadiationResidual_Relative_Final
    REAL(DP), DIMENSION(:),     INTENT(out)   :: &
      ConservativeProjectionRadiationResidual_Relative_Final
    REAL(DP), DIMENSION(:),     INTENT(out)   :: &
      ConservativeProjectionMatterHeating_Rate_Final
    INTEGER,  DIMENSION(:),     INTENT(out)   :: ConservativeProjectionStatus_Final
    REAL(DP), DIMENSION(:,:,:), INTENT(out)   :: LeptonFlavorAudit_Rate_Final
#endif

    ! --- Local Variables ---

    INTEGER  :: ierr, myid
    INTEGER  :: k_outer, Mk_outer, nX_P_outer
    INTEGER  :: k_inner, Mk_inner, nX_P_inner
    INTEGER  :: iN_X
#if defined( TWOMOMENT_ORDER_V )
    INTEGER  :: k_backtrack, k_terminal
#endif

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: P
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: Error

    LOGICAL,  ALLOCATABLE, DIMENSION(:) :: ITERATE_outer, ITERATE_inner
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: PackIndex_outer, UnpackIndex_outer
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: PackIndex_inner, UnpackIndex_inner

#if defined( TWOMOMENT_ORDER_V )
    LOGICAL :: TerminalProjectionMode
    LOGICAL :: TerminalOpacityCacheCurrent, TerminalRatesCacheCurrent
    LOGICAL,  ALLOCATABLE, DIMENSION(:) :: &
      ProjectionMask, RadiationBacktrackMask
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: &
      T_Conservative, Y_Conservative, E_Conservative, &
      V_u_1_Conservative, V_u_2_Conservative, V_u_3_Conservative
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: &
      FVECm_projection, GVECm_projection, &
      InitialRadiationLimiter_ConservationCorrection, &
      RadiationLeptonInvariant
#endif

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
    Error = 0

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

#if defined( TWOMOMENT_ORDER_V )
    ALLOCATE( ProjectionMask(nX_G) )
    ALLOCATE( RadiationBacktrackMask(nX_G) )
    ALLOCATE( T_Conservative(nX_G), Y_Conservative(nX_G) )
    ALLOCATE( E_Conservative(nX_G) )
    ALLOCATE( V_u_1_Conservative(nX_G) )
    ALLOCATE( V_u_2_Conservative(nX_G) )
    ALLOCATE( V_u_3_Conservative(nX_G) )
    ALLOCATE( FVECm_projection(n_FP_outer,nX_G) )
    ALLOCATE( GVECm_projection(n_FP_outer,nX_G) )
    ALLOCATE &
      ( InitialRadiationLimiter_ConservationCorrection &
          (nMatterEquations,nX_G) )
    ALLOCATE( RadiationLeptonInvariant(nLeptonFamilies,nX_G) )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      nIterations_Inner(iN_X) = 0
      nIterations_Outer(iN_X) = 0
      nIterations_Inner_Final(iN_X) = 0
      nIterations_Outer_Final(iN_X) = 0
      RelativeResidual_Inner_First(iN_X) = Zero

#if defined( TWOMOMENT_ORDER_V )
      MatterEquationResidual_Rate_Final(:,iN_X) = Zero
      InitialRadiationLimiter_Rate_Final(:,iN_X) = Zero
      RadiationResidual_Relative_Final(iN_X) = Zero
      ConservativeProjectionRadiationResidual_Relative_Final(iN_X) = Zero
      ConservativeProjectionMatterHeating_Rate_Final(iN_X) = Zero
      ConservativeProjectionStatus_Final(iN_X) = -2
      LeptonFlavorAudit_Rate_Final(:,:,iN_X) = Zero
#endif

    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Error ) &
    !$OMP MAP( alloc: &
    !$OMP   ITERATE_outer, ITERATE_inner, &
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
    !$ACC COPYIN( Error ) &
    !$ACC CREATE( &
    !$ACC   ITERATE_outer, ITERATE_inner, &
    !$ACC   PackIndex_outer, UnpackIndex_outer, &
    !$ACC   AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC   BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC   WORK_outer, TAU_outer, Alpha_outer, &
    !$ACC   PackIndex_inner, UnpackIndex_inner, &
    !$ACC   AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC   BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC   WORK_inner, TAU_inner, Alpha_inner )
#endif

#if defined( TWOMOMENT_ORDER_V )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: ProjectionMask, RadiationBacktrackMask, &
    !$OMP             T_Conservative, Y_Conservative, E_Conservative, &
    !$OMP             V_u_1_Conservative, V_u_2_Conservative, &
    !$OMP             V_u_3_Conservative, &
    !$OMP             FVECm_projection, GVECm_projection, &
    !$OMP             InitialRadiationLimiter_ConservationCorrection, &
    !$OMP             RadiationLeptonInvariant )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( ProjectionMask, RadiationBacktrackMask, &
    !$ACC         T_Conservative, Y_Conservative, E_Conservative, &
    !$ACC         V_u_1_Conservative, V_u_2_Conservative, &
    !$ACC         V_u_3_Conservative, &
    !$ACC         FVECm_projection, GVECm_projection, &
    !$ACC         InitialRadiationLimiter_ConservationCorrection, &
    !$ACC         RadiationLeptonInvariant )
#endif
#endif


    CALL SetOpacityMask &
         ( nX_G, nDOFX, D, T, Y, Opacity_Mask )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( Opacity_Mask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( Opacity_Mask )
#endif

    IF( PRESENT( MatterCoupling_Mask_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MatterCoupling_Mask_Option )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        ITERATE_outer(iN_X) = &
            Opacity_Mask(iOp_EmAb,   iN_X) .OR. &
            Opacity_Mask(iOp_ECTable,iN_X) .OR. &
            Opacity_Mask(iOp_Iso,    iN_X) .OR. &
            ( Include_NNS    .AND. Opacity_Mask(iOp_NNS,   iN_X) ) .OR. &
            ( Include_NES    .AND. Opacity_Mask(iOp_NES,   iN_X) ) .OR. &
            ( Include_Pair   .AND. Opacity_Mask(iOp_Pair,  iN_X) ) .OR. &
            ( Include_NuPair .AND. Opacity_Mask(iOp_NuPair,iN_X) ) .OR. &
            ( Include_Brem   .AND. Opacity_Mask(iOp_Brem,  iN_X) )

        MatterCoupling_Mask_Option(iN_X) = ITERATE_outer(iN_X)

      END DO

    ELSE

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        ITERATE_outer(iN_X) = &
            Opacity_Mask(iOp_EmAb,   iN_X) .OR. &
            Opacity_Mask(iOp_ECTable,iN_X) .OR. &
            Opacity_Mask(iOp_Iso,    iN_X) .OR. &
            ( Include_NNS    .AND. Opacity_Mask(iOp_NNS,   iN_X) ) .OR. &
            ( Include_NES    .AND. Opacity_Mask(iOp_NES,   iN_X) ) .OR. &
            ( Include_Pair   .AND. Opacity_Mask(iOp_Pair,  iN_X) ) .OR. &
            ( Include_NuPair .AND. Opacity_Mask(iOp_NuPair,iN_X) ) .OR. &
            ( Include_Brem   .AND. Opacity_Mask(iOp_Brem,  iN_X) )

      END DO

    END IF

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( ITERATE_outer )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( ITERATE_outer )
#endif

#if defined( TWOMOMENT_ORDER_V )
    ProjectionMask = ITERATE_outer
    TerminalProjectionMode = .FALSE.
    TerminalOpacityCacheCurrent = .FALSE.
    TerminalRatesCacheCurrent = .FALSE.
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( ProjectionMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( ProjectionMask )
#endif
#endif

    SqrtGm = SQRT( Gm_dd_11 * Gm_dd_22 * Gm_dd_33 )

    ! --- Initial RHS ---

    CALL TimersStart( Timer_Collisions_InitializeRHS )

#if   defined( TWOMOMENT_ORDER_1 )

    CALL InitializeRHS_OrderOne &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

    CALL ApplyNeutrinoDistributionLimit_OrderV &
           ( dt, ITERATE_outer, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Nnu, &
             V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             InitialRadiationLimiter_Rate_Final, &
             InitialRadiationLimiter_ConservationCorrection, &
             LeptonFlavorAudit_Rate_Final(:,iLFA_InputConversion,:), &
             LeptonFlavorAudit_Rate_Final(:,iLFA_EntryProjection,:) )

    CALL InitializeRHS_OrderV &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, T, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
             InitialRadiationLimiter_ConservationCorrection )

    IF( UseRadiationOnlySolve ) &
      CALL InitializeRadiationLeptonInvariants_OrderV &
             ( ITERATE_outer, C_Dnu, RadiationLeptonInvariant )

#elif defined( TWOMOMENT_RELATIVISTIC )

    CALL ComputePressure_TABLE & 
           ( D, T, Y, P )

    CALL InitializeRHS_Relativistic &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Y, E, P, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

#endif

    CALL TimersStop( Timer_Collisions_InitializeRHS )

    CALL CreatePackIndex &
           ( ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

    ! --- Compute Opacity Kernels ---

#if defined( TWOMOMENT_ORDER_V )
    IF( .NOT. UseRadiationOnlySolve .OR. FreezeOpacities )THEN
#endif

    CALL TimersStart( Timer_Collisions_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( D, T, Y, SqrtGm, ITERATE_outer, nX_P_outer, &
             PackIndex_outer, UnpackIndex_outer, &
             T0 = T_old, Y0 = Y_old, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL TimersStop( Timer_Collisions_ComputeOpacity )

#if defined( TWOMOMENT_ORDER_V )
    END IF
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( UseRadiationOnlySolve )THEN

      CALL TimersStart( Timer_Collisions_OuterLoop )
      CALL TimersStart( Timer_Collisions_InnerLoop )

      ITERATE_inner = ITERATE_outer
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( ITERATE_inner )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( ITERATE_inner )
#endif

      k_outer = 0
      k_inner = 0

      DO WHILE( ANY( ITERATE_inner(:) ) .AND. k_inner < MaxIter_inner )

        k_inner = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_inner, nX_P_inner, &
                 PackIndex_inner, UnpackIndex_inner )

        CALL TimersStart( Timer_Collisions_MatterFromRadiation )

        CALL RecoverMatterFromRadiationInvariants_OrderV &
               ( ITERATE_inner, ProjectionMask, &
                 FVECm_outer, GVECm_outer, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 D, T, Y, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 InitialRadiationLimiter_ConservationCorrection, &
                 RadiationLeptonInvariant, &
                 nX_P_inner, PackIndex_inner, UnpackIndex_inner, Error )

        CALL TimersStop( Timer_Collisions_MatterFromRadiation )

        RadiationBacktrackMask = ITERATE_inner .AND. Error /= 0

        IF( ANY( RadiationBacktrackMask ) .AND. k_inner > 1 )THEN

          CALL TimersStart( Timer_Collisions_RadiationBacktrack )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( RadiationBacktrackMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( RadiationBacktrackMask )
#endif

          k_backtrack = 0

          DO WHILE &
            ( ANY( RadiationBacktrackMask ) .AND. k_backtrack < 20 )

            k_backtrack = k_backtrack + 1

            CALL BacktrackRadiationOnlyIterate_OrderV &
                   ( RadiationBacktrackMask, MAX( 1, Mk_inner - 1 ), &
                     FVEC_inner, GVEC_inner, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33 )

            CALL ProjectRadiationOnlyAndersonUpdate_OrderV &
                   ( RadiationBacktrackMask, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     V_u_1, V_u_2, V_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                     RadiationLeptonInvariant )

            CALL CreatePackIndex &
                   ( RadiationBacktrackMask, nX_P_inner, &
                     PackIndex_inner, UnpackIndex_inner )

            CALL TimersStart( Timer_Collisions_MatterFromRadiation )

            CALL RecoverMatterFromRadiationInvariants_OrderV &
                   ( RadiationBacktrackMask, ProjectionMask, &
                     FVECm_outer, GVECm_outer, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     D, T, Y, E, V_u_1, V_u_2, V_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                     InitialRadiationLimiter_ConservationCorrection, &
                     RadiationLeptonInvariant, &
                     nX_P_inner, PackIndex_inner, UnpackIndex_inner, Error )

            CALL TimersStop( Timer_Collisions_MatterFromRadiation )

            RadiationBacktrackMask &
              = RadiationBacktrackMask .AND. Error /= 0

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( RadiationBacktrackMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( RadiationBacktrackMask )
#endif

          END DO

          CALL CreatePackIndex &
                 ( ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner )

          CALL TimersStop( Timer_Collisions_RadiationBacktrack )

        END IF

        CALL CheckError_FP &
               ( ITERATE_outer, ITERATE_inner, 0, k_inner, Error, &
                 D, Y, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 IgnoreIterationLimit_Option = .TRUE. )

        IF( .NOT. FreezeOpacities )THEN

          CALL TimersStart( Timer_Collisions_ComputeOpacity )

          CALL ComputeOpacities_Packed &
                 ( D, T, Y, SqrtGm, &
                   ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner, &
                   T0 = T_old, Y0 = Y_old, &
                   Opacity_Mask_Input = Opacity_Mask )

          CALL TimersStop( Timer_Collisions_ComputeOpacity )

        END IF

        CALL ComputeDnuNorm( ITERATE_inner, Dnu )

        CALL TimersStart( Timer_Collisions_ComputeRates )

        IF( FreezeOpacities )THEN

          CALL ComputeRates_Packed &
                 ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner, &
                   nX_P_outer, Opacity_Mask )

        ELSE

          CALL ComputeRates_Packed &
                 ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner, &
                   nX_P_inner, Opacity_Mask )

        END IF

        CALL TimersStop( Timer_Collisions_ComputeRates )

        CALL TimersStart( Timer_Collisions_NeutrinoRHS )

        CALL ComputeNeutrinoRHS_OrderV &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, dt, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        CALL ProjectRadiationOnlyFixedPointMap_OrderV &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        CALL TimersStop( Timer_Collisions_NeutrinoRHS )

        CALL TimersStart( Timer_Collisions_CheckInner )

        CALL CheckConvergence_RadiationOnly &
               ( ITERATE_inner, k_inner, &
                 nIterations_Inner, nIterations_Outer, &
                 nIterations_Inner_Final, nIterations_Outer_Final, &
                 RelativeResidual_Inner_First, &
                 RadiationResidual_Relative_Final, FVECm_inner )

        CALL TimersStop( Timer_Collisions_CheckInner )

        IF( ANY( ITERATE_inner(:) ) )THEN

          CALL TimersStart( Timer_Collisions_SolveLS )

          CALL SolveLS_FP &
                 ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                   FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner, &
                   AMAT_inner, BVEC_inner, Alpha_inner, TAU_inner, &
                   LWORK_inner, WORK_inner )

          CALL TimersStop( Timer_Collisions_SolveLS )

          CALL TimersStart( Timer_Collisions_UpdateFP )

          CALL UpdateNeutrinoRHS_OrderV &
                 ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33 )

          CALL ProjectRadiationOnlyAndersonUpdate_OrderV &
                 ( ITERATE_inner, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   RadiationLeptonInvariant )

          CALL ShiftRHS_FP &
                 ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                   FVEC_inner, GVEC_inner )

          CALL TimersStop( Timer_Collisions_UpdateFP )

        END IF

      END DO

      IF( ANY( ITERATE_inner(:) ) )THEN

        CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

        IF( myid == 0 )THEN
          WRITE(ERROR_UNIT,'(A,I0,A)') &
            'ERROR: Invariant-constrained radiation solve failed after ', &
            k_inner, ' iterations.'
          WRITE(ERROR_UNIT,'(A)') &
            '       The nested solver fallback is intentionally disabled.'
          FLUSH(ERROR_UNIT)
        END IF

        CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

      END IF

      CALL TimersStop( Timer_Collisions_InnerLoop )
      CALL TimersStop( Timer_Collisions_OuterLoop )

    END IF

    IF( .NOT. UseRadiationOnlySolve )THEN
#endif

    ! --- Start Outer Loop ---

    CALL TimersStart( Timer_Collisions_OuterLoop )

    k_outer = 0
#if defined( TWOMOMENT_ORDER_V )
    k_terminal = 0
    DO WHILE &
      ( ANY( ITERATE_outer(:) ) &
        .AND. (       ( .NOT. TerminalProjectionMode &
                        .AND. k_outer < MaxIter_outer ) &
                .OR. ( TerminalProjectionMode &
                        .AND. k_terminal < MaxIter_outer ) ) )

      IF( TerminalProjectionMode )THEN
        k_terminal = k_terminal + 1
      ELSE
        k_outer = k_outer + 1
        Mk_outer = MIN( M_outer, k_outer )
      END IF
#else
    DO WHILE( ANY( ITERATE_outer(:) ) .AND. k_outer < MaxIter_outer )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )
#endif

      Error = 0

      CALL ComputeDnuNorm( ITERATE_outer, Dnu )

      IF ( k_outer > 1 ) THEN

        CALL CreatePackIndex &
               ( ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        IF ( .NOT. FreezeOpacities ) THEN

#if defined( TWOMOMENT_ORDER_V )
          IF( .NOT. TerminalOpacityCacheCurrent )THEN
#endif

          ! --- Recompute Opacity Kernels ---

          CALL TimersStart( Timer_Collisions_ComputeOpacity )

          CALL ComputeOpacities_Packed &
                 ( D, T, Y, SqrtGm, ITERATE_outer, nX_P_outer, &
                   PackIndex_outer, UnpackIndex_outer, &
                   T0 = T_old, Y0 = Y_old, &
                   Opacity_Mask_Input = Opacity_Mask )

          CALL TimersStop( Timer_Collisions_ComputeOpacity )

#if defined( TWOMOMENT_ORDER_V )
          END IF
#endif

        END IF

#if defined( TWOMOMENT_ORDER_V )
        TerminalOpacityCacheCurrent = .FALSE.
#endif

      END IF

      ! --- Start Inner Loop ---

      CALL TimersStart( Timer_Collisions_InnerLoop )

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        IF( ITERATE_outer(iN_X) )THEN
          nIterations_Inner_Final(iN_X) = 0
#if defined( TWOMOMENT_ORDER_V )
          IF( TerminalProjectionMode )THEN
            nIterations_Outer(iN_X) = nIterations_Outer(iN_X) + 1
            nIterations_Outer_Final(iN_X) = k_outer + k_terminal
          ELSE
            nIterations_Outer_Final(iN_X) = k_outer
          END IF
#else
          nIterations_Outer_Final(iN_X) = k_outer
#endif
          RelativeResidual_Inner_First(iN_X) = Zero
        END IF

      END DO

      ITERATE_inner = ITERATE_outer
#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET UPDATE TO( ITERATE_inner )
#elif defined( THORNADO_OACC   )
      !$ACC UPDATE DEVICE( ITERATE_inner )
#endif

      k_inner = 0
      DO WHILE( ANY( ITERATE_inner(:) ) .AND. k_inner < MaxIter_inner )

        k_inner  = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_inner, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

#if defined( TWOMOMENT_ORDER_V )
        IF( .NOT. (       TerminalProjectionMode &
                   .AND. TerminalRatesCacheCurrent &
                   .AND. k_inner == 1 ) )THEN
#endif
        CALL TimersStart( Timer_Collisions_ComputeRates )

        CALL ComputeRates_Packed &
               ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, ITERATE_inner, nX_P_inner, &
                 PackIndex_inner, UnpackIndex_inner, nX_P_outer, Opacity_Mask )

        CALL TimersStop( Timer_Collisions_ComputeRates )
#if defined( TWOMOMENT_ORDER_V )
        END IF

        IF( k_inner == 1 ) TerminalRatesCacheCurrent = .FALSE.
#endif

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

        ! --- Update Solution Vectors (inner) ---

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
                 nIterations_Inner, nIterations_Inner_Final, &
                 RelativeResidual_Inner_First, FVECm_inner )

        CALL CheckError_Inner &
               ( ITERATE_inner, k_outer, k_inner, &
                 D, Y, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        CALL TimersStop( Timer_Collisions_CheckInner )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Collisions_UpdateFP )

      END DO ! --- Inner Loop ---

      CALL TimersStop( Timer_Collisions_InnerLoop )

#if defined( TWOMOMENT_ORDER_V )
      IF( UseConservativeProjection .AND. TerminalProjectionMode )THEN

        ProjectionMask = ITERATE_outer
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( ProjectionMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( ProjectionMask )
#endif

        CALL ApplyConservativeProjection_OrderV &
               ( dt, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 D, T, Y, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                 ProjectionMask, ITERATE_outer, &
                 nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
                 FVECm_inner, GVECm_inner, &
                 FVECm_projection, GVECm_projection, &
                 T_Conservative, Y_Conservative, E_Conservative, &
                 V_u_1_Conservative, V_u_2_Conservative, &
                 V_u_3_Conservative, Error, &
                 RadiationResidual_Relative_Final, &
                 ConservativeProjectionRadiationResidual_Relative_Final, &
                 ConservativeProjectionMatterHeating_Rate_Final, &
                 ConservativeProjectionStatus_Final )

        TerminalOpacityCacheCurrent = ANY( ITERATE_outer )
        TerminalRatesCacheCurrent = ANY( ITERATE_outer )

        CALL CheckError_FP &
               ( ITERATE_inner, ITERATE_inner, k_outer, k_inner, Error, &
                 D, Y, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        IF(      ANY( ITERATE_outer ) &
            .AND. k_terminal >= MaxIter_outer )THEN

          CALL CheckError_FP &
                 ( ITERATE_outer, ITERATE_inner, &
                   MaxIter_outer, k_inner, Error, &
                   D, Y, E, T, V_u_1, V_u_2, V_u_3, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        END IF

        CYCLE

      END IF
#endif

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

      ! --- Update Solution Vectors (outer) ---

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

      ! --- Check Convergence (outer) ---

      CALL TimersStart( Timer_Collisions_CheckOuter )

      CALL CheckConvergence_Outer &
             ( ITERATE_outer, ITERATE_inner, n_FP_outer, k_outer, &
               nIterations_Outer, FVECm_outer )

      CALL CheckError_FP &
             ( ITERATE_outer, ITERATE_inner, k_outer, k_inner, Error, &
               D, Y, E, T, V_u_1, V_u_2, V_u_3, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL TimersStop( Timer_Collisions_CheckOuter )

      IF ( .NOT. ANY( ITERATE_outer ) )THEN

#if defined( TWOMOMENT_ORDER_V )
        IF( UseConservativeProjection )THEN

          TerminalProjectionMode = .TRUE.

          CALL ApplyConservativeProjection_OrderV &
               ( dt, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 D, T, Y, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                 ProjectionMask, ITERATE_outer, &
                 nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
                 FVECm_inner, GVECm_inner, &
                 FVECm_projection, GVECm_projection, &
                 T_Conservative, Y_Conservative, E_Conservative, &
                 V_u_1_Conservative, V_u_2_Conservative, &
                 V_u_3_Conservative, Error, &
                 RadiationResidual_Relative_Final, &
                 ConservativeProjectionRadiationResidual_Relative_Final, &
                 ConservativeProjectionMatterHeating_Rate_Final, &
                 ConservativeProjectionStatus_Final )

          TerminalOpacityCacheCurrent = ANY( ITERATE_outer )
          TerminalRatesCacheCurrent = ANY( ITERATE_outer )

          CALL CheckError_FP &
               ( ITERATE_inner, ITERATE_inner, k_outer, k_inner, Error, &
                 D, Y, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        END IF
#endif

        CALL TimersStop( Timer_Collisions_UpdateFP )

        CYCLE

      END IF

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Collisions_UpdateFP )

    END DO ! --- Outer Loop ---

    CALL TimersStop( Timer_Collisions_OuterLoop )

#if defined( TWOMOMENT_ORDER_V )
    END IF
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( ComputeMatterEquationResidual )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G
        ProjectionMask(iN_X) = nIterations_Outer_Final(iN_X) > 0
      END DO

      CALL ComputeMatterRHS_OrderV &
             ( ProjectionMask, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL CreatePackIndex &
             ( ProjectionMask, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

      CALL ComputeOpacities_Packed &
             ( D, T, Y, SqrtGm, ProjectionMask, nX_P_inner, &
               PackIndex_inner, UnpackIndex_inner, &
               T0 = T_old, Y0 = Y_old, Opacity_Mask_Input = Opacity_Mask )

      CALL ComputeRates_Packed &
             ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               ProjectionMask, nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
               nX_P_inner, Opacity_Mask )

      CALL ComputeNeutrinoRHS_OrderV &
             ( ProjectionMask, FVECm_inner, GVECm_inner, dt, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL ComputeLeptonFlavorAudit_OrderV &
             ( ProjectionMask, dt, &
               UseRadiationOnlySolve, &
               FVECm_inner, GVECm_inner, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               LeptonFlavorAudit_Rate_Final )

    END IF

    ! Convert the final matter-equation residual after the conservative
    ! projection into the physical rate variables used by the diagnostics.

    IF( ComputeMatterEquationResidual )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      MatterEquationResidual_Rate_Final(:,iN_X) = Zero

      IF( nIterations_Outer_Final(iN_X) > 0 )THEN

        MatterEquationResidual_Rate_Final(iQ_L,iN_X) &
          = D_old(iN_X) * Y_old(iN_X) / AtomicMassUnit &
              * FVECm_outer(iY,iN_X) / dt

        MatterEquationResidual_Rate_Final(iQ_E,iN_X) &
          = D_old(iN_X) * Ef_old(iN_X) &
              * FVECm_outer(iEf,iN_X) / dt

        MatterEquationResidual_Rate_Final(iQ_M1,iN_X) &
          = D_old(iN_X) * SpeedOfLight &
              * FVECm_outer(iV1,iN_X) / dt
        MatterEquationResidual_Rate_Final(iQ_M2,iN_X) &
          = D_old(iN_X) * SpeedOfLight &
              * FVECm_outer(iV2,iN_X) / dt
        MatterEquationResidual_Rate_Final(iQ_M3,iN_X) &
          = D_old(iN_X) * SpeedOfLight &
              * FVECm_outer(iV3,iN_X) / dt

        MatterEquationResidual_Rate_Final(iQ_H,iN_X) &
          = MatterEquationResidual_Rate_Final(iQ_E,iN_X) &
              - V_u_1(iN_X) &
                  * MatterEquationResidual_Rate_Final(iQ_M1,iN_X) &
              - V_u_2(iN_X) &
                  * MatterEquationResidual_Rate_Final(iQ_M2,iN_X) &
              - V_u_3(iN_X) &
                  * MatterEquationResidual_Rate_Final(iQ_M3,iN_X)

      END IF

    END DO

    END IF
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: &
    !$OMP   Error, ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
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
    !$ACC   Error, ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
    !$ACC   AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC   BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC   WORK_outer, TAU_outer, Alpha_outer, &
    !$ACC   ITERATE_inner, PackIndex_inner, UnpackIndex_inner, &
    !$ACC   AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC   BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC   WORK_inner, TAU_inner, Alpha_inner )
#endif

#if defined( TWOMOMENT_ORDER_V )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ProjectionMask, RadiationBacktrackMask, &
    !$OMP               T_Conservative, Y_Conservative, E_Conservative, &
    !$OMP               V_u_1_Conservative, V_u_2_Conservative, &
    !$OMP               V_u_3_Conservative, &
    !$OMP               FVECm_projection, GVECm_projection, &
    !$OMP               InitialRadiationLimiter_ConservationCorrection, &
    !$OMP               RadiationLeptonInvariant )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ProjectionMask, RadiationBacktrackMask, &
    !$ACC         T_Conservative, Y_Conservative, E_Conservative, &
    !$ACC         V_u_1_Conservative, V_u_2_Conservative, &
    !$ACC         V_u_3_Conservative, &
    !$ACC         FVECm_projection, GVECm_projection, &
    !$ACC         InitialRadiationLimiter_ConservationCorrection, &
    !$ACC         RadiationLeptonInvariant )
#endif
#endif

  END SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE InitializeRadiationLeptonInvariants_OrderV &
    ( MASK, Nnu, RadiationLeptonInvariant )

    LOGICAL,  DIMENSION(:),     INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Nnu
    REAL(DP), DIMENSION(:,:),   INTENT(out) :: RadiationLeptonInvariant

    INTEGER  :: iFamily, iN_E, iN_X, iS
    REAL(DP) :: LeptonIntegral

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_E, iS, LeptonIntegral )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG COLLAPSE(2) &
    !$ACC PRIVATE( iN_E, iS, LeptonIntegral )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_E, iS, LeptonIntegral )
#endif
    DO iN_X = 1, nX_G
    DO iFamily = 1, nLeptonFamilies

      LeptonIntegral = Zero

      IF( MASK(iN_X) )THEN
        DO iS = 1, nSpecies
          IF( ( iS + 1 ) / 2 /= iFamily ) CYCLE
          DO iN_E = 1, nE_G
            LeptonIntegral = LeptonIntegral &
              + W2_S(iN_E) * LeptonNumber(iS) * Nnu(iN_E,iS,iN_X)
          END DO
        END DO
      END IF

      RadiationLeptonInvariant(iFamily,iN_X) = LeptonIntegral

    END DO
    END DO

  END SUBROUTINE InitializeRadiationLeptonInvariants_OrderV


  SUBROUTINE FindRadiationLeptonInvariantViolations_OrderV &
    ( MASK, ViolationMask, RadiationLeptonInvariant, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:),     INTENT(in)  :: MASK
    LOGICAL,  DIMENSION(:),     INTENT(out) :: ViolationMask
    REAL(DP), DIMENSION(:,:),   INTENT(in)  :: RadiationLeptonInvariant
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)  :: &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP), PARAMETER :: InvariantTolerance = 64.0_DP * EPSILON(One)

    INTEGER  :: iFamily, iN_E, iN_X, iS
    REAL(DP) :: Current, DeltaLepton, LeptonScale, Residual, vDotI

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE( iFamily, iN_E, iS, Current, DeltaLepton, &
    !$OMP          LeptonScale, Residual, vDotI )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( iFamily, iN_E, iS, Current, DeltaLepton, &
    !$ACC          LeptonScale, Residual, vDotI )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iFamily, iN_E, iS, Current, DeltaLepton, &
    !$OMP          LeptonScale, Residual, vDotI )
#endif
    DO iN_X = 1, nX_G

      ViolationMask(iN_X) = .FALSE.
      IF( .NOT. MASK(iN_X) ) CYCLE

      DO iFamily = 2, MIN( nLeptonFamilies, ( nSpecies + 1 ) / 2 )

        Current = Zero
        LeptonScale = Zero

        DO iS = 1, nSpecies
          IF( ( iS + 1 ) / 2 /= iFamily ) CYCLE
          DO iN_E = 1, nE_G
            vDotI = V_u_1(iN_X) * Gm_dd_11(iN_X) &
                      * Inu_u_1(iN_E,iS,iN_X) &
                  + V_u_2(iN_X) * Gm_dd_22(iN_X) &
                      * Inu_u_2(iN_E,iS,iN_X) &
                  + V_u_3(iN_X) * Gm_dd_33(iN_X) &
                      * Inu_u_3(iN_E,iS,iN_X)
            DeltaLepton = W2_S(iN_E) * LeptonNumber(iS) &
              * ( Dnu(iN_E,iS,iN_X) + vDotI )
            Current = Current + DeltaLepton
            LeptonScale = LeptonScale + ABS( DeltaLepton )
          END DO
        END DO

        Residual = RadiationLeptonInvariant(iFamily,iN_X) - Current
        IF( ABS( Residual ) > InvariantTolerance &
              * MAX( LeptonScale, &
                     ABS( RadiationLeptonInvariant(iFamily,iN_X) ) ) )THEN
          ViolationMask(iN_X) = .TRUE.
          EXIT
        END IF

      END DO

    END DO

  END SUBROUTINE FindRadiationLeptonInvariantViolations_OrderV


  SUBROUTINE EnforceRadiationLeptonInvariants_OrderV &
    ( MASK, RadiationLeptonInvariant, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: RadiationLeptonInvariant
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER, PARAMETER :: MaxInvariantCorrections = 4
    REAL(DP), PARAMETER :: InvariantTolerance = 64.0_DP * EPSILON(One)

    INTEGER  :: iFamily, iN_E, iN_X, iS, k
    REAL(DP) :: Capacity, Current, DeltaLepton, D_End, FluxMagnitude
    REAL(DP) :: FluxScale, Fraction, I1_End, I2_End, I3_End
    REAL(DP) :: LeptonScale, N_Current, N_End, Residual, SignedChange
    REAL(DP) :: vDotI, vDotI_End

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_E, iS, k, Capacity, Current, DeltaLepton, D_End, &
    !$OMP          FluxMagnitude, FluxScale, Fraction, &
    !$OMP          I1_End, I2_End, I3_End, LeptonScale, &
    !$OMP          N_Current, N_End, Residual, SignedChange, &
    !$OMP          vDotI, vDotI_End )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG COLLAPSE(2) &
    !$ACC PRIVATE( iN_E, iS, k, Capacity, Current, DeltaLepton, D_End, &
    !$ACC          FluxMagnitude, FluxScale, Fraction, &
    !$ACC          I1_End, I2_End, I3_End, LeptonScale, &
    !$ACC          N_Current, N_End, Residual, SignedChange, &
    !$ACC          vDotI, vDotI_End )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_E, iS, k, Capacity, Current, DeltaLepton, D_End, &
    !$OMP          FluxMagnitude, FluxScale, Fraction, &
    !$OMP          I1_End, I2_End, I3_End, LeptonScale, &
    !$OMP          N_Current, N_End, Residual, SignedChange, &
    !$OMP          vDotI, vDotI_End )
#endif
    DO iN_X = 1, nX_G
    DO iFamily = 2, MIN( nLeptonFamilies, ( nSpecies + 1 ) / 2 )

      IF( .NOT. MASK(iN_X) ) CYCLE

      DO k = 1, MaxInvariantCorrections

        Current = Zero
        LeptonScale = Zero

        DO iS = 1, nSpecies
          IF( ( iS + 1 ) / 2 /= iFamily ) CYCLE
          DO iN_E = 1, nE_G
            vDotI = V_u_1(iN_X) * Gm_dd_11(iN_X) &
                      * Inu_u_1(iN_E,iS,iN_X) &
                  + V_u_2(iN_X) * Gm_dd_22(iN_X) &
                      * Inu_u_2(iN_E,iS,iN_X) &
                  + V_u_3(iN_X) * Gm_dd_33(iN_X) &
                      * Inu_u_3(iN_E,iS,iN_X)
            N_Current = Dnu(iN_E,iS,iN_X) + vDotI
            DeltaLepton = W2_S(iN_E) * LeptonNumber(iS) * N_Current
            Current = Current + DeltaLepton
            LeptonScale = LeptonScale + ABS( DeltaLepton )
          END DO
        END DO

        Residual = RadiationLeptonInvariant(iFamily,iN_X) - Current
        IF( ABS( Residual ) <= InvariantTolerance &
              * MAX( LeptonScale, &
                     ABS( RadiationLeptonInvariant(iFamily,iN_X) ) ) ) EXIT

        Capacity = Zero

        DO iS = 1, nSpecies
          IF( ( iS + 1 ) / 2 /= iFamily ) CYCLE
          DO iN_E = 1, nE_G
            vDotI = V_u_1(iN_X) * Gm_dd_11(iN_X) &
                      * Inu_u_1(iN_E,iS,iN_X) &
                  + V_u_2(iN_X) * Gm_dd_22(iN_X) &
                      * Inu_u_2(iN_E,iS,iN_X) &
                  + V_u_3(iN_X) * Gm_dd_33(iN_X) &
                      * Inu_u_3(iN_E,iS,iN_X)
            N_Current = Dnu(iN_E,iS,iN_X) + vDotI

            ! Reduce only the species whose removal moves the signed family
            ! residual toward zero; do not populate empty energy groups.
            IF( Residual * LeptonNumber(iS) >= Zero ) CYCLE

            FluxMagnitude = SQRT( MAX( Zero, &
                Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)**2 &
              + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)**2 &
              + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)**2 ) )
            IF( FluxMagnitude > Zero )THEN
              FluxScale = MIN( One, SqrtTiny / FluxMagnitude )
            ELSE
              FluxScale = Zero
            END IF
            D_End = SqrtTiny
            I1_End = FluxScale * Inu_u_1(iN_E,iS,iN_X)
            I2_End = FluxScale * Inu_u_2(iN_E,iS,iN_X)
            I3_End = FluxScale * Inu_u_3(iN_E,iS,iN_X)
            vDotI_End = V_u_1(iN_X) * Gm_dd_11(iN_X) * I1_End &
                      + V_u_2(iN_X) * Gm_dd_22(iN_X) * I2_End &
                      + V_u_3(iN_X) * Gm_dd_33(iN_X) * I3_End
            N_End = D_End + vDotI_End
            SignedChange = W2_S(iN_E) * LeptonNumber(iS) &
              * ( N_End - N_Current )

            IF( Residual * SignedChange > Zero ) &
              Capacity = Capacity + SignedChange
          END DO
        END DO

        IF( Residual * Capacity <= Zero ) EXIT
        Fraction = MIN( One, ABS( Residual / Capacity ) )

        DO iS = 1, nSpecies
          IF( ( iS + 1 ) / 2 /= iFamily ) CYCLE
          DO iN_E = 1, nE_G
            vDotI = V_u_1(iN_X) * Gm_dd_11(iN_X) &
                      * Inu_u_1(iN_E,iS,iN_X) &
                  + V_u_2(iN_X) * Gm_dd_22(iN_X) &
                      * Inu_u_2(iN_E,iS,iN_X) &
                  + V_u_3(iN_X) * Gm_dd_33(iN_X) &
                      * Inu_u_3(iN_E,iS,iN_X)
            N_Current = Dnu(iN_E,iS,iN_X) + vDotI

            IF( Residual * LeptonNumber(iS) >= Zero ) CYCLE

            FluxMagnitude = SQRT( MAX( Zero, &
                Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)**2 &
              + Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)**2 &
              + Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)**2 ) )
            IF( FluxMagnitude > Zero )THEN
              FluxScale = MIN( One, SqrtTiny / FluxMagnitude )
            ELSE
              FluxScale = Zero
            END IF
            D_End = SqrtTiny
            I1_End = FluxScale * Inu_u_1(iN_E,iS,iN_X)
            I2_End = FluxScale * Inu_u_2(iN_E,iS,iN_X)
            I3_End = FluxScale * Inu_u_3(iN_E,iS,iN_X)
            vDotI_End = V_u_1(iN_X) * Gm_dd_11(iN_X) * I1_End &
                      + V_u_2(iN_X) * Gm_dd_22(iN_X) * I2_End &
                      + V_u_3(iN_X) * Gm_dd_33(iN_X) * I3_End
            N_End = D_End + vDotI_End
            SignedChange = W2_S(iN_E) * LeptonNumber(iS) &
              * ( N_End - N_Current )

            IF( Residual * SignedChange > Zero )THEN
              Dnu(iN_E,iS,iN_X) = Dnu(iN_E,iS,iN_X) &
                + Fraction * ( D_End - Dnu(iN_E,iS,iN_X) )
              Inu_u_1(iN_E,iS,iN_X) = Inu_u_1(iN_E,iS,iN_X) &
                + Fraction * ( I1_End - Inu_u_1(iN_E,iS,iN_X) )
              Inu_u_2(iN_E,iS,iN_X) = Inu_u_2(iN_E,iS,iN_X) &
                + Fraction * ( I2_End - Inu_u_2(iN_E,iS,iN_X) )
              Inu_u_3(iN_E,iS,iN_X) = Inu_u_3(iN_E,iS,iN_X) &
                + Fraction * ( I3_End - Inu_u_3(iN_E,iS,iN_X) )
            END IF
          END DO
        END DO

      END DO

    END DO
    END DO

  END SUBROUTINE EnforceRadiationLeptonInvariants_OrderV


  SUBROUTINE ProjectRadiationOnlyAndersonUpdate_OrderV &
    ( MASK, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      RadiationLeptonInvariant )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: RadiationLeptonInvariant

    INTEGER :: iN_E, iN_X, iS

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP FIRSTPRIVATE( DnuMax )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC FIRSTPRIVATE( DnuMax )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) FIRSTPRIVATE( DnuMax )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) ) &
        CALL ProjectPrimitiveRadiationMoment_OrderV &
               ( Dnu    (iN_E,iS,iN_X), &
                 Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), &
                 Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 DnuMax )

    END DO
    END DO
    END DO

    CALL EnforceRadiationLeptonInvariants_OrderV &
           ( MASK, RadiationLeptonInvariant, &
             Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
             V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

  END SUBROUTINE ProjectRadiationOnlyAndersonUpdate_OrderV


  SUBROUTINE ProjectRadiationOnlyFixedPointMap_OrderV &
    ( MASK, Fm, Gm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iOS, iS
    REAL(DP) :: Dnu_Map, Inu_Map_u_1, Inu_Map_u_2, Inu_Map_u_3

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP FIRSTPRIVATE( DnuMax ) &
    !$OMP PRIVATE( iOS, Dnu_Map, Inu_Map_u_1, Inu_Map_u_2, Inu_Map_u_3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC FIRSTPRIVATE( DnuMax ) &
    !$ACC PRIVATE( iOS, Dnu_Map, Inu_Map_u_1, Inu_Map_u_2, Inu_Map_u_3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) FIRSTPRIVATE( DnuMax ) &
    !$OMP PRIVATE( iOS, Dnu_Map, Inu_Map_u_1, Inu_Map_u_2, Inu_Map_u_3 )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        Dnu_Map = Gm(iOS+iCR_N,iN_X)
        Inu_Map_u_1 = Gm(iOS+iCR_G1,iN_X) / Gm_dd_11(iN_X)
        Inu_Map_u_2 = Gm(iOS+iCR_G2,iN_X) / Gm_dd_22(iN_X)
        Inu_Map_u_3 = Gm(iOS+iCR_G3,iN_X) / Gm_dd_33(iN_X)

        CALL ProjectPrimitiveRadiationMoment_OrderV &
               ( Dnu_Map, Inu_Map_u_1, Inu_Map_u_2, Inu_Map_u_3, &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 DnuMax )

        Gm(iOS+iCR_N ,iN_X) = Dnu_Map
        Gm(iOS+iCR_G1,iN_X) = Inu_Map_u_1 * Gm_dd_11(iN_X)
        Gm(iOS+iCR_G2,iN_X) = Inu_Map_u_2 * Gm_dd_22(iN_X)
        Gm(iOS+iCR_G3,iN_X) = Inu_Map_u_3 * Gm_dd_33(iN_X)

        Fm(iOS+iCR_N ,iN_X) &
          = Gm(iOS+iCR_N ,iN_X) - Dnu(iN_E,iS,iN_X)
        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) &
              - Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)
        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) &
              - Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)
        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) &
              - Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ProjectRadiationOnlyFixedPointMap_OrderV


  SUBROUTINE ComputeLeptonFlavorAudit_OrderV &
    ( MASK, dt, ApplyMapProjection, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, AuditRate )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP),                   INTENT(in)    :: dt
    LOGICAL,                    INTENT(in)    :: ApplyMapProjection
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: AuditRate

    INTEGER  :: iFamily, iN_E, iN_X, iOS, iS
    REAL(DP) :: Chi_T, CollisionRate, SignedWeight
    REAL(DP) :: D_Map, I_Map_u_1, I_Map_u_2, I_Map_u_3
    REAL(DP) :: N_Current, N_Map, N_Old, N_Projected
    REAL(DP) :: vDotI_Current, vDotI_Old

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE( iFamily, iN_E, iOS, iS, Chi_T, CollisionRate, &
    !$OMP          SignedWeight, D_Map, I_Map_u_1, I_Map_u_2, I_Map_u_3, &
    !$OMP          N_Current, N_Map, N_Old, N_Projected, &
    !$OMP          vDotI_Current, vDotI_Old ) FIRSTPRIVATE( DnuMax )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( iFamily, iN_E, iOS, iS, Chi_T, CollisionRate, &
    !$ACC          SignedWeight, D_Map, I_Map_u_1, I_Map_u_2, I_Map_u_3, &
    !$ACC          N_Current, N_Map, N_Old, N_Projected, &
    !$ACC          vDotI_Current, vDotI_Old ) FIRSTPRIVATE( DnuMax )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iFamily, iN_E, iOS, iS, Chi_T, CollisionRate, &
    !$OMP          SignedWeight, D_Map, I_Map_u_1, I_Map_u_2, I_Map_u_3, &
    !$OMP          N_Current, N_Map, N_Old, N_Projected, &
    !$OMP          vDotI_Current, vDotI_Old ) FIRSTPRIVATE( DnuMax )
#endif
    DO iN_X = 1, nX_G

      AuditRate(:,iLFA_CollisionSource:nLeptonFlavorAudit,iN_X) = Zero
      IF( .NOT. MASK(iN_X) ) CYCLE

      DO iS = 1, nSpecies
      DO iN_E = 1, nE_G

        iFamily = ( iS + 1 ) / 2
        IF( iFamily > nLeptonFamilies ) CYCLE

        SignedWeight = W2_S(iN_E) * LeptonNumber(iS)
        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        vDotI_Current &
          = V_u_1(iN_X) * Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X) &
          + V_u_2(iN_X) * Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X) &
          + V_u_3(iN_X) * Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)
        vDotI_Old &
          = V_u_1_old(iN_X) * Gm_dd_11(iN_X) * Inu_u_1_old(iN_E,iS,iN_X) &
          + V_u_2_old(iN_X) * Gm_dd_22(iN_X) * Inu_u_2_old(iN_E,iS,iN_X) &
          + V_u_3_old(iN_X) * Gm_dd_33(iN_X) * Inu_u_3_old(iN_E,iS,iN_X)

        N_Current = Dnu(iN_E,iS,iN_X) + vDotI_Current
        N_Old = Dnu_old(iN_E,iS,iN_X) + vDotI_Old

        D_Map = Gm(iOS+iCR_N,iN_X)
        I_Map_u_1 = Gm(iOS+iCR_G1,iN_X) / Gm_dd_11(iN_X)
        I_Map_u_2 = Gm(iOS+iCR_G2,iN_X) / Gm_dd_22(iN_X)
        I_Map_u_3 = Gm(iOS+iCR_G3,iN_X) / Gm_dd_33(iN_X)
        N_Map = D_Map &
          + V_u_1(iN_X) * Gm(iOS+iCR_G1,iN_X) &
          + V_u_2(iN_X) * Gm(iOS+iCR_G2,iN_X) &
          + V_u_3(iN_X) * Gm(iOS+iCR_G3,iN_X)

        Chi_T = Chi_EmAb(iN_E,iS,iN_X) &
              + Chi_NNS(iN_E,iS,iN_X) + Chi_NES(iN_E,iS,iN_X) &
              + Chi_Pair(iN_E,iS,iN_X) + Chi_NuPair(iN_E,iS,iN_X) &
              + Chi_Brem(iN_E,iS,iN_X)
        CollisionRate &
          = ( Fm(iOS+iCR_N,iN_X) &
                * ( One + Omega(iN_X) * dt * Chi_T ) / Omega(iN_X) &
              - ( C_Dnu(iN_E,iS,iN_X) - vDotI_Current &
                    - Dnu(iN_E,iS,iN_X) ) ) / dt

        N_Projected = N_Map
        IF( ApplyMapProjection )THEN
          CALL ProjectPrimitiveRadiationMoment_OrderV &
                 ( D_Map, I_Map_u_1, I_Map_u_2, I_Map_u_3, &
                   Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                   DnuMax )
          N_Projected = D_Map &
            + V_u_1(iN_X) * Gm_dd_11(iN_X) * I_Map_u_1 &
            + V_u_2(iN_X) * Gm_dd_22(iN_X) * I_Map_u_2 &
            + V_u_3(iN_X) * Gm_dd_33(iN_X) * I_Map_u_3
        END IF

        AuditRate(iFamily,iLFA_CollisionSource,iN_X) &
          = AuditRate(iFamily,iLFA_CollisionSource,iN_X) &
              + SignedWeight * CollisionRate
        AuditRate(iFamily,iLFA_RawMap,iN_X) &
          = AuditRate(iFamily,iLFA_RawMap,iN_X) &
              + SignedWeight * ( N_Map - N_Old ) / dt
        AuditRate(iFamily,iLFA_MapProjection,iN_X) &
          = AuditRate(iFamily,iLFA_MapProjection,iN_X) &
              + SignedWeight * ( N_Projected - N_Map ) / dt
        AuditRate(iFamily,iLFA_FixedPointResidual,iN_X) &
          = AuditRate(iFamily,iLFA_FixedPointResidual,iN_X) &
              + SignedWeight * ( N_Projected - N_Current ) / dt
        AuditRate(iFamily,iLFA_EquationResidual,iN_X) &
          = AuditRate(iFamily,iLFA_EquationResidual,iN_X) &
              + SignedWeight &
                  * ( ( N_Current - N_Old ) / dt - CollisionRate )

      END DO
      END DO

    END DO

  END SUBROUTINE ComputeLeptonFlavorAudit_OrderV


  SUBROUTINE ProjectPrimitiveRadiationMoment_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, DnuUpperBound )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, DnuUpperBound

    REAL(DP) :: Dnu_Projected, FluxMagnitude
    REAL(DP) :: FluxMagnitude_Projected, FluxScale

    FluxMagnitude &
      = SQRT( MAX( Zero, &
            Gm_dd_11 * Inu_u_1**2 &
          + Gm_dd_22 * Inu_u_2**2 &
          + Gm_dd_33 * Inu_u_3**2 ) )

    ! Orthogonal projection in the local orthonormal moment variables onto
    ! the second-order cone |I| <= Dnu, followed by the occupation bound.
    ! Each energy/species tuple is projected independently.

    IF( FluxMagnitude <= Dnu )THEN

      Dnu_Projected = Dnu
      FluxMagnitude_Projected = FluxMagnitude

    ELSE IF( FluxMagnitude <= - Dnu )THEN

      Dnu_Projected = Zero
      FluxMagnitude_Projected = Zero

    ELSE

      Dnu_Projected = Half * ( Dnu + FluxMagnitude )
      FluxMagnitude_Projected = Dnu_Projected

    END IF

    IF( Dnu_Projected > DnuUpperBound )THEN
      Dnu_Projected = DnuUpperBound
      FluxMagnitude_Projected = MIN( FluxMagnitude, DnuUpperBound )
    ELSE IF( Dnu_Projected < SqrtTiny )THEN
      Dnu_Projected = SqrtTiny
      FluxMagnitude_Projected = MIN( FluxMagnitude, SqrtTiny )
    END IF

    IF( FluxMagnitude > Zero )THEN
      FluxScale = FluxMagnitude_Projected / FluxMagnitude
    ELSE
      FluxScale = Zero
    END IF

    Dnu = Dnu_Projected
    Inu_u_1 = FluxScale * Inu_u_1
    Inu_u_2 = FluxScale * Inu_u_2
    Inu_u_3 = FluxScale * Inu_u_3

  END SUBROUTINE ProjectPrimitiveRadiationMoment_OrderV


  SUBROUTINE BacktrackRadiationOnlyIterate_OrderV &
    ( MASK, iHistory, F, G, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: iHistory
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: F, G
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iOS, iS
    REAL(DP) :: Previous

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iOS, Previous )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iOS, Previous )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE( iOS, Previous )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        Previous &
          = G(iOS+iCR_N,iHistory,iN_X) &
              - F(iOS+iCR_N,iHistory,iN_X)
        Dnu(iN_E,iS,iN_X) = Half * ( Dnu(iN_E,iS,iN_X) + Previous )

        Previous &
          = G(iOS+iCR_G1,iHistory,iN_X) &
              - F(iOS+iCR_G1,iHistory,iN_X)
        Inu_u_1(iN_E,iS,iN_X) &
          = Half &
              * ( Inu_u_1(iN_E,iS,iN_X) &
                    + Previous / Gm_dd_11(iN_X) )

        Previous &
          = G(iOS+iCR_G2,iHistory,iN_X) &
              - F(iOS+iCR_G2,iHistory,iN_X)
        Inu_u_2(iN_E,iS,iN_X) &
          = Half &
              * ( Inu_u_2(iN_E,iS,iN_X) &
                    + Previous / Gm_dd_22(iN_X) )

        Previous &
          = G(iOS+iCR_G3,iHistory,iN_X) &
              - F(iOS+iCR_G3,iHistory,iN_X)
        Inu_u_3(iN_E,iS,iN_X) &
          = Half &
              * ( Inu_u_3(iN_E,iS,iN_X) &
                    + Previous / Gm_dd_33(iN_X) )

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE BacktrackRadiationOnlyIterate_OrderV


  SUBROUTINE RecoverMatterFromRadiationState_OrderV &
    ( MASK, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, Y, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      RadiationLimiterCorrection, Error )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: RadiationLimiterCorrection
    INTEGER,  DIMENSION(:),     INTENT(inout) :: Error

    INTEGER  :: iN_E, iN_X, iS, SolveInfo
    REAL(DP) :: A(3,3), B(3), VHat(3)
    REAL(DP) :: Ef, SqrtG1, SqrtG2, SqrtG3, vDotV
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13
    REAL(DP) :: k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: SUM_D, SUM_I1, SUM_I2, SUM_I3
    REAL(DP) :: SUM_K11, SUM_K12, SUM_K13
    REAL(DP) :: SUM_K22, SUM_K23, SUM_K33
    REAL(DP) :: SUM_L, SUM_LI1, SUM_LI2, SUM_LI3
    REAL(DP) :: SUM_Ef_Old, SUM_L_Old
    REAL(DP) :: SUM_P1_Old, SUM_P2_Old, SUM_P3_Old
    REAL(DP) :: SUM_Y, SUM_Ef
    REAL(DP) :: U_Y_New, U_Ef_New
    REAL(DP) :: U_V1_New, U_V2_New, U_V3_New

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( A, B, VHat, SolveInfo, Ef, SqrtG1, SqrtG2, SqrtG3, &
    !$OMP          vDotV, SUM_D, SUM_I1, SUM_I2, SUM_I3, &
    !$OMP          SUM_K11, SUM_K12, SUM_K13, &
    !$OMP          SUM_K22, SUM_K23, SUM_K33, &
    !$OMP          SUM_L, SUM_LI1, SUM_LI2, SUM_LI3, &
    !$OMP          SUM_Ef_Old, SUM_L_Old, &
    !$OMP          SUM_P1_Old, SUM_P2_Old, SUM_P3_Old, SUM_Y, SUM_Ef, &
    !$OMP          U_Y_New, U_Ef_New, U_V1_New, U_V2_New, U_V3_New )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( A, B, VHat, SolveInfo, Ef, SqrtG1, SqrtG2, SqrtG3, &
    !$ACC          vDotV, SUM_D, SUM_I1, SUM_I2, SUM_I3, &
    !$ACC          SUM_K11, SUM_K12, SUM_K13, &
    !$ACC          SUM_K22, SUM_K23, SUM_K33, &
    !$ACC          SUM_L, SUM_LI1, SUM_LI2, SUM_LI3, &
    !$ACC          SUM_Ef_Old, SUM_L_Old, &
    !$ACC          SUM_P1_Old, SUM_P2_Old, SUM_P3_Old, SUM_Y, SUM_Ef, &
    !$ACC          U_Y_New, U_Ef_New, U_V1_New, U_V2_New, U_V3_New )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( A, B, VHat, SolveInfo, Ef, SqrtG1, SqrtG2, SqrtG3, &
    !$OMP          vDotV, k_dd_11, k_dd_12, k_dd_13, &
    !$OMP          k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          SUM_D, SUM_I1, SUM_I2, SUM_I3, &
    !$OMP          SUM_K11, SUM_K12, SUM_K13, &
    !$OMP          SUM_K22, SUM_K23, SUM_K33, &
    !$OMP          SUM_L, SUM_LI1, SUM_LI2, SUM_LI3, &
    !$OMP          SUM_Ef_Old, SUM_L_Old, &
    !$OMP          SUM_P1_Old, SUM_P2_Old, SUM_P3_Old, SUM_Y, SUM_Ef, &
    !$OMP          U_Y_New, U_Ef_New, U_V1_New, U_V2_New, U_V3_New )
#endif
    DO iN_X = 1, nX_G

      IF( .NOT. MASK(iN_X) ) CYCLE

      Error(iN_X) = 0
      INFO(iN_X) = 0

      SqrtG1 = SQRT( Gm_dd_11(iN_X) )
      SqrtG2 = SQRT( Gm_dd_22(iN_X) )
      SqrtG3 = SQRT( Gm_dd_33(iN_X) )

      SUM_D   = Zero
      SUM_I1  = Zero
      SUM_I2  = Zero
      SUM_I3  = Zero
      SUM_K11 = Zero
      SUM_K12 = Zero
      SUM_K13 = Zero
      SUM_K22 = Zero
      SUM_K23 = Zero
      SUM_K33 = Zero
      SUM_L   = Zero
      SUM_LI1 = Zero
      SUM_LI2 = Zero
      SUM_LI3 = Zero
      SUM_Ef_Old = Zero
      SUM_L_Old  = Zero
      SUM_P1_Old = Zero
      SUM_P2_Old = Zero
      SUM_P3_Old = Zero

#if   defined(THORNADO_OMP_OL)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, &
    !$OMP          k_dd_22, k_dd_23, k_dd_33 ) &
    !$OMP REDUCTION( + : SUM_D, SUM_I1, SUM_I2, SUM_I3, &
    !$OMP                 SUM_K11, SUM_K12, SUM_K13, &
    !$OMP                 SUM_K22, SUM_K23, SUM_K33, &
    !$OMP                 SUM_L, SUM_LI1, SUM_LI2, SUM_LI3, &
    !$OMP                 SUM_Ef_Old, SUM_L_Old, &
    !$OMP                 SUM_P1_Old, SUM_P2_Old, SUM_P3_Old )
#elif defined(THORNADO_OACC)
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, &
    !$ACC          k_dd_22, k_dd_23, k_dd_33 ) &
    !$ACC REDUCTION( + : SUM_D, SUM_I1, SUM_I2, SUM_I3, &
    !$ACC                 SUM_K11, SUM_K12, SUM_K13, &
    !$ACC                 SUM_K22, SUM_K23, SUM_K33, &
    !$ACC                 SUM_L, SUM_LI1, SUM_LI2, SUM_LI3, &
    !$ACC                 SUM_Ef_Old, SUM_L_Old, &
    !$ACC                 SUM_P1_Old, SUM_P2_Old, SUM_P3_Old )
#endif
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        CALL ComputeEddingtonTensorComponents_dd &
               ( Dnu(iN_E,iS,iN_X), &
                 Inu_u_1(iN_E,iS,iN_X), &
                 Inu_u_2(iN_E,iS,iN_X), &
                 Inu_u_3(iN_E,iS,iN_X), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 k_dd_11, k_dd_12, k_dd_13, &
                 k_dd_22, k_dd_23, k_dd_33 )

        SUM_D = SUM_D + W3_S(iN_E) * Dnu(iN_E,iS,iN_X)
        SUM_I1 = SUM_I1 + W3_S(iN_E) * SqrtG1 &
          * Inu_u_1(iN_E,iS,iN_X)
        SUM_I2 = SUM_I2 + W3_S(iN_E) * SqrtG2 &
          * Inu_u_2(iN_E,iS,iN_X)
        SUM_I3 = SUM_I3 + W3_S(iN_E) * SqrtG3 &
          * Inu_u_3(iN_E,iS,iN_X)

        SUM_K11 = SUM_K11 + W3_S(iN_E) * Dnu(iN_E,iS,iN_X) &
          * k_dd_11 / ( SqrtG1 * SqrtG1 )
        SUM_K12 = SUM_K12 + W3_S(iN_E) * Dnu(iN_E,iS,iN_X) &
          * k_dd_12 / ( SqrtG1 * SqrtG2 )
        SUM_K13 = SUM_K13 + W3_S(iN_E) * Dnu(iN_E,iS,iN_X) &
          * k_dd_13 / ( SqrtG1 * SqrtG3 )
        SUM_K22 = SUM_K22 + W3_S(iN_E) * Dnu(iN_E,iS,iN_X) &
          * k_dd_22 / ( SqrtG2 * SqrtG2 )
        SUM_K23 = SUM_K23 + W3_S(iN_E) * Dnu(iN_E,iS,iN_X) &
          * k_dd_23 / ( SqrtG2 * SqrtG3 )
        SUM_K33 = SUM_K33 + W3_S(iN_E) * Dnu(iN_E,iS,iN_X) &
          * k_dd_33 / ( SqrtG3 * SqrtG3 )

        SUM_Ef_Old = SUM_Ef_Old + W3_S(iN_E) &
          * ( Two * C_Dnu(iN_E,iS,iN_X) - Dnu_old(iN_E,iS,iN_X) )
        SUM_P1_Old = SUM_P1_Old + W3_S(iN_E) &
          * ( C_Inu_d_1(iN_E,iS,iN_X) / SqrtG1 &
              + SqrtG1 * V_u_1_old(iN_X) * Dnu_old(iN_E,iS,iN_X) )
        SUM_P2_Old = SUM_P2_Old + W3_S(iN_E) &
          * ( C_Inu_d_2(iN_E,iS,iN_X) / SqrtG2 &
              + SqrtG2 * V_u_2_old(iN_X) * Dnu_old(iN_E,iS,iN_X) )
        SUM_P3_Old = SUM_P3_Old + W3_S(iN_E) &
          * ( C_Inu_d_3(iN_E,iS,iN_X) / SqrtG3 &
              + SqrtG3 * V_u_3_old(iN_X) * Dnu_old(iN_E,iS,iN_X) )

        IF( iS <= iNuE_Bar )THEN
          SUM_L_Old = SUM_L_Old + W2_S(iN_E) * LeptonNumber(iS) &
            * C_Dnu(iN_E,iS,iN_X)
          SUM_L = SUM_L + W2_S(iN_E) * LeptonNumber(iS) &
            * Dnu(iN_E,iS,iN_X)
          SUM_LI1 = SUM_LI1 + W2_S(iN_E) * LeptonNumber(iS) &
            * SqrtG1 * Inu_u_1(iN_E,iS,iN_X)
          SUM_LI2 = SUM_LI2 + W2_S(iN_E) * LeptonNumber(iS) &
            * SqrtG2 * Inu_u_2(iN_E,iS,iN_X)
          SUM_LI3 = SUM_LI3 + W2_S(iN_E) * LeptonNumber(iS) &
            * SqrtG3 * Inu_u_3(iN_E,iS,iN_X)
        END IF

      END DO
      END DO

      A = Zero
      A(1,1) = One + wMatrRHS(iV1) &
        * ( SUM_D + SUM_K11 ) / D(iN_X)
      A(1,2) = wMatrRHS(iV1) * SUM_K12 / D(iN_X)
      A(1,3) = wMatrRHS(iV1) * SUM_K13 / D(iN_X)
      A(2,1) = wMatrRHS(iV2) * SUM_K12 / D(iN_X)
      A(2,2) = One + wMatrRHS(iV2) &
        * ( SUM_D + SUM_K22 ) / D(iN_X)
      A(2,3) = wMatrRHS(iV2) * SUM_K23 / D(iN_X)
      A(3,1) = wMatrRHS(iV3) * SUM_K13 / D(iN_X)
      A(3,2) = wMatrRHS(iV3) * SUM_K23 / D(iN_X)
      A(3,3) = One + wMatrRHS(iV3) &
        * ( SUM_D + SUM_K33 ) / D(iN_X)

      B(1) = SqrtG1 * V_u_1_old(iN_X) &
        + wMatrRHS(iV1) / D(iN_X) &
          * ( SUM_P1_Old &
              + RadiationLimiterCorrection(iV1,iN_X) / SqrtG1 &
              - SUM_I1 )
      B(2) = SqrtG2 * V_u_2_old(iN_X) &
        + wMatrRHS(iV2) / D(iN_X) &
          * ( SUM_P2_Old &
              + RadiationLimiterCorrection(iV2,iN_X) / SqrtG2 &
              - SUM_I2 )
      B(3) = SqrtG3 * V_u_3_old(iN_X) &
        + wMatrRHS(iV3) / D(iN_X) &
          * ( SUM_P3_Old &
              + RadiationLimiterCorrection(iV3,iN_X) / SqrtG3 &
              - SUM_I3 )

      CALL SolveLinearSystem3x3_OrderV( A, B, VHat, SolveInfo )
      INFO(iN_X) = SolveInfo

      IF( SolveInfo /= 0 )THEN
        Error(iN_X) = 1
        Fm(:,iN_X) = Zero
        Gm(iY ,iN_X) = U_Y    (iN_X)
        Gm(iEf,iN_X) = U_Ef   (iN_X)
        Gm(iV1,iN_X) = U_V_d_1(iN_X)
        Gm(iV2,iN_X) = U_V_d_2(iN_X)
        Gm(iV3,iN_X) = U_V_d_3(iN_X)
        CYCLE
      END IF

      V_u_1(iN_X) = VHat(1) / SqrtG1
      V_u_2(iN_X) = VHat(2) / SqrtG2
      V_u_3(iN_X) = VHat(3) / SqrtG3

      SUM_Y = SUM_L &
        + VHat(1) * SUM_LI1 &
        + VHat(2) * SUM_LI2 &
        + VHat(3) * SUM_LI3
      SUM_Ef = SUM_D + Two &
        * (   VHat(1) * SUM_I1 &
            + VHat(2) * SUM_I2 &
            + VHat(3) * SUM_I3 )

      Y(iN_X) = Y_old(iN_X) + AtomicMassUnit / D(iN_X) &
        * wMatrRHS(iY) &
        * (   SUM_L_Old + RadiationLimiterCorrection(iY,iN_X) &
            - SUM_Y )
      Ef = Ef_old(iN_X) + wMatrRHS(iEf) / D(iN_X) &
        * (   SUM_Ef_Old + RadiationLimiterCorrection(iEf,iN_X) &
            - SUM_Ef )

      U_Y_New = Y(iN_X) / Y_old(iN_X)
      U_Ef_New = Ef / Ef_old(iN_X)
      U_V1_New = SqrtG1 * VHat(1) / SpeedOfLight
      U_V2_New = SqrtG2 * VHat(2) / SpeedOfLight
      U_V3_New = SqrtG3 * VHat(3) / SpeedOfLight

      Fm(iY ,iN_X) = U_Y_New  - U_Y    (iN_X)
      Fm(iEf,iN_X) = U_Ef_New - U_Ef   (iN_X)
      Fm(iV1,iN_X) = U_V1_New - U_V_d_1(iN_X)
      Fm(iV2,iN_X) = U_V2_New - U_V_d_2(iN_X)
      Fm(iV3,iN_X) = U_V3_New - U_V_d_3(iN_X)

      U_Y    (iN_X) = U_Y_New
      U_Ef   (iN_X) = U_Ef_New
      U_V_d_1(iN_X) = U_V1_New
      U_V_d_2(iN_X) = U_V2_New
      U_V_d_3(iN_X) = U_V3_New

      G_Y    (iN_X) = U_Y_New
      G_Ef   (iN_X) = U_Ef_New
      G_V_d_1(iN_X) = U_V1_New
      G_V_d_2(iN_X) = U_V2_New
      G_V_d_3(iN_X) = U_V3_New

      Gm(iY ,iN_X) = U_Y_New
      Gm(iEf,iN_X) = U_Ef_New
      Gm(iV1,iN_X) = U_V1_New
      Gm(iV2,iN_X) = U_V2_New
      Gm(iV3,iN_X) = U_V3_New

      vDotV = VHat(1)**2 + VHat(2)**2 + VHat(3)**2
      E(iN_X) = Ef - Half * vDotV
      Omega(iN_X) = One / ( One + SQRT( vDotV ) )

    END DO

  END SUBROUTINE RecoverMatterFromRadiationState_OrderV


  SUBROUTINE SolveLinearSystem3x3_OrderV( A_Input, B_Input, X, Info )
#if   defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: A_Input(3,3), B_Input(3)
    REAL(DP), INTENT(out) :: X(3)
    INTEGER,  INTENT(out) :: Info

    INTEGER  :: i, j, k, Pivot
    REAL(DP) :: A(3,3), B(3), Factor, PivotMagnitude, Scale, Swap

    Scale = MAX( MAXVAL( ABS( A_Input ) ), One )
    A = A_Input / Scale
    B = B_Input / Scale
    Info = 0

    DO k = 1, 2

      Pivot = k
      PivotMagnitude = ABS( A(k,k) )
      DO i = k + 1, 3
        IF( ABS( A(i,k) ) > PivotMagnitude )THEN
          Pivot = i
          PivotMagnitude = ABS( A(i,k) )
        END IF
      END DO

      IF( PivotMagnitude <= 128.0_DP * EPSILON(One) )THEN
        Info = k
        X = Zero
        RETURN
      END IF

      IF( Pivot /= k )THEN
        DO j = k, 3
          Swap = A(k,j)
          A(k,j) = A(Pivot,j)
          A(Pivot,j) = Swap
        END DO
        Swap = B(k)
        B(k) = B(Pivot)
        B(Pivot) = Swap
      END IF

      DO i = k + 1, 3
        Factor = A(i,k) / A(k,k)
        A(i,k) = Zero
        DO j = k + 1, 3
          A(i,j) = A(i,j) - Factor * A(k,j)
        END DO
        B(i) = B(i) - Factor * B(k)
      END DO

    END DO

    IF( ABS( A(3,3) ) <= 128.0_DP * EPSILON(One) )THEN
      Info = 3
      X = Zero
      RETURN
    END IF

    X(3) = B(3) / A(3,3)
    X(2) = ( B(2) - A(2,3) * X(3) ) / A(2,2)
    X(1) = ( B(1) - A(1,2) * X(2) - A(1,3) * X(3) ) / A(1,1)

  END SUBROUTINE SolveLinearSystem3x3_OrderV


  SUBROUTINE RecoverMatterFromRadiationInvariants_OrderV &
    ( MASK, WorkingMask, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, T, Y, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      RadiationLimiterCorrection, &
      RadiationLeptonInvariant, &
      nX_P, PackIndex, UnpackIndex, Error )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: WorkingMask
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      T, Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: RadiationLimiterCorrection
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: RadiationLeptonInvariant
    INTEGER,                    INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(:),     INTENT(in)    :: PackIndex, UnpackIndex
    INTEGER,  DIMENSION(:),     INTENT(inout) :: Error

    INTEGER, PARAMETER :: MaxIter_InvariantCoupling = 8

    INTEGER :: iN_X, kInvariant

    WorkingMask = MASK

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( WorkingMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( WorkingMask )
#endif

    kInvariant = 0
    DO WHILE &
      ( ANY( WorkingMask(:) ) &
        .AND. kInvariant < MaxIter_InvariantCoupling )

      kInvariant = kInvariant + 1

      CALL TimersStart( Timer_Collisions_MatterRHS )

      CALL RecoverMatterFromRadiationState_OrderV &
             ( WorkingMask, Fm, Gm, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, Y, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               RadiationLimiterCorrection, Error )

      CALL TimersStop( Timer_Collisions_MatterRHS )

      CALL FindRadiationLeptonInvariantViolations_OrderV &
             ( MASK, WorkingMask, RadiationLeptonInvariant, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( WorkingMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( WorkingMask )
#endif

      IF( .NOT. ANY( WorkingMask(:) ) ) EXIT

      CALL EnforceRadiationLeptonInvariants_OrderV &
             ( WorkingMask, RadiationLeptonInvariant, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      IF( WorkingMask(iN_X) ) INFO(iN_X) = 1

    END DO

    CALL UpdateTemperature_Packed &
           ( D, E, Y, T, MASK, nX_P, PackIndex, UnpackIndex, Error )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      IF( MASK(iN_X) .AND. INFO(iN_X) /= 0 ) Error(iN_X) = 1

    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( Error )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( Error )
#endif

  END SUBROUTINE RecoverMatterFromRadiationInvariants_OrderV


  SUBROUTINE CheckConvergence_RadiationOnly &
    ( MASK, k, nIterations_Inner, nIterations_Outer, &
      nIterations_Inner_Final, nIterations_Outer_Final, &
      RelativeResidual_First, RelativeResidual_Final, Fm )

    LOGICAL,  DIMENSION(:),   INTENT(inout) :: MASK
    INTEGER,                  INTENT(in)    :: k
    INTEGER,  DIMENSION(:),   INTENT(inout) :: &
      nIterations_Inner, nIterations_Outer
    INTEGER,  DIMENSION(:),   INTENT(inout) :: &
      nIterations_Inner_Final, nIterations_Outer_Final
    REAL(DP), DIMENSION(:),   INTENT(inout) :: &
      RelativeResidual_First, RelativeResidual_Final
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Fm

    LOGICAL  :: Converged
    INTEGER  :: iN_E, iN_X, iOS, iS
    REAL(DP) :: Fnorm, Fnorm_G1, Fnorm_G2, Fnorm_G3, Fnorm_N
    REAL(DP) :: RelativeResidual

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Converged, iN_E, iOS, iS, Fnorm, Fnorm_G1, Fnorm_G2, &
    !$OMP          Fnorm_G3, Fnorm_N, RelativeResidual )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( Converged, iN_E, iOS, iS, Fnorm, Fnorm_G1, Fnorm_G2, &
    !$ACC          Fnorm_G3, Fnorm_N, RelativeResidual )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Converged, iN_E, iOS, iS, Fnorm, Fnorm_G1, Fnorm_G2, &
    !$OMP          Fnorm_G3, Fnorm_N, RelativeResidual )
#endif
    DO iN_X = 1, nX_G

      IF( MASK(iN_X) )THEN

        Converged = .TRUE.
        RelativeResidual_Final(iN_X) = Zero

        DO iS = 1, nSpecies

          Fnorm_N  = Zero
          Fnorm_G1 = Zero
          Fnorm_G2 = Zero
          Fnorm_G3 = Zero

          DO iN_E = 1, nE_G

            iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

            Fnorm_N &
              = Fnorm_N + W2_S(iN_E) * Fm(iOS+iCR_N,iN_X)**2
            Fnorm_G1 &
              = Fnorm_G1 + W2_S(iN_E) * Fm(iOS+iCR_G1,iN_X)**2
            Fnorm_G2 &
              = Fnorm_G2 + W2_S(iN_E) * Fm(iOS+iCR_G2,iN_X)**2
            Fnorm_G3 &
              = Fnorm_G3 + W2_S(iN_E) * Fm(iOS+iCR_G3,iN_X)**2

          END DO

          Fnorm = SQRT( MAX( Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 ) )

          IF( k == 1 ) Fnorm0(iS,iN_X) = Fnorm

          IF( DnuNorm(iS,iN_X) > Zero )THEN
            RelativeResidual = Fnorm / DnuNorm(iS,iN_X)
          ELSE IF( Fnorm == Zero )THEN
            RelativeResidual = Zero
          ELSE
            RelativeResidual = HUGE(One)
          END IF

          RelativeResidual_Final(iN_X) &
            = MAX( RelativeResidual_Final(iN_X), RelativeResidual )

          Converged = Converged .AND. &
            ( ( Fnorm <= Rtol_inner * DnuNorm(iS,iN_X) ) .OR. &
              ( Fnorm <= Atol_inner * Fnorm0(iS,iN_X) ) )

        END DO

        IF( k == 1 ) &
          RelativeResidual_First(iN_X) = RelativeResidual_Final(iN_X)

        nIterations_Inner(iN_X) = nIterations_Inner(iN_X) + 1
        nIterations_Outer(iN_X) = k
        nIterations_Inner_Final(iN_X) = 1
        nIterations_Outer_Final(iN_X) = k

        IF( Converged ) MASK(iN_X) = .FALSE.

      END IF

    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( MASK, nIterations_Inner, nIterations_Outer )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( MASK, nIterations_Inner, nIterations_Outer )
#endif

  END SUBROUTINE CheckConvergence_RadiationOnly


  SUBROUTINE ApplyConservativeProjection_OrderV &
    ( dt, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, T, Y, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
      ProjectionMask, ContinuationMask, &
      nX_P, PackIndex, UnpackIndex, &
      FVECm_inner, GVECm_inner, &
      FVECm_projection, GVECm_projection, &
      T_Conservative, Y_Conservative, E_Conservative, &
      V_u_1_Conservative, V_u_2_Conservative, V_u_3_Conservative, &
      ProjectionError, &
      RadiationResidual_Relative_Final, &
      ConservativeProjectionRadiationResidual_Relative_Final, &
      ConservativeProjectionMatterHeating_Rate_Final, &
      ConservativeProjectionStatus_Final )

    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: T, Y, E
    REAL(DP), DIMENSION(:),     INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:),     INTENT(in)    :: SqrtGm
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: ProjectionMask
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: ContinuationMask
    INTEGER,                    INTENT(inout) :: nX_P
    INTEGER,  DIMENSION(:),     INTENT(inout) :: PackIndex, UnpackIndex
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: FVECm_inner, GVECm_inner
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: &
      FVECm_projection, GVECm_projection
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      T_Conservative, Y_Conservative, E_Conservative
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      V_u_1_Conservative, V_u_2_Conservative, V_u_3_Conservative
    INTEGER,  DIMENSION(:),     INTENT(inout) :: ProjectionError
    REAL(DP), DIMENSION(:),     INTENT(inout) :: RadiationResidual_Relative_Final
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      ConservativeProjectionRadiationResidual_Relative_Final
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      ConservativeProjectionMatterHeating_Rate_Final
    INTEGER,  DIMENSION(:),     INTENT(inout) :: ConservativeProjectionStatus_Final

    INTEGER, PARAMETER :: MaxIter_ConservativeProjection = 32
    REAL(DP), PARAMETER :: Rtol_ConservativeProjection &
      = 100.0_DP * EPSILON(One)

    INTEGER  :: iN_X, k_projection
    REAL(DP) :: Ef, ProjectionResidual, Q_E, Q_M1, Q_M2, Q_M3
    REAL(DP) :: V_d_1, V_d_2, V_d_3, vDotV
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      ContinuationMask(iN_X) = .FALSE.

      T_Conservative(iN_X) = T(iN_X)
      Y_Conservative(iN_X) = Y(iN_X)
      E_Conservative(iN_X) = E(iN_X)
      V_u_1_Conservative(iN_X) = V_u_1(iN_X)
      V_u_2_Conservative(iN_X) = V_u_2(iN_X)
      V_u_3_Conservative(iN_X) = V_u_3(iN_X)

      IF( ProjectionMask(iN_X) )THEN
        ProjectionError(iN_X) = 0
        RadiationResidual_Relative_Final(iN_X) = Zero
        ConservativeProjectionRadiationResidual_Relative_Final(iN_X) = Zero
        ConservativeProjectionMatterHeating_Rate_Final(iN_X) = Zero
        ConservativeProjectionStatus_Final(iN_X) = 0
      END IF

    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( ProjectionMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( ProjectionMask )
#endif

    CALL ComputeDnuNorm( ProjectionMask, Dnu )

    CALL ComputeRadiationResidualRelative &
           ( ProjectionMask, FVECm_inner, &
             RadiationResidual_Relative_Final )

    ! The Order-v Eulerian radiation energy and momentum depend on the matter
    ! velocity.  Consequently, one application of the outer fixed-point map
    ! is not yet an exactly conservative matter state when it changes the
    ! velocity.  With the radiation moments fixed, iterate only this cheap
    ! matter-conservation map before evaluating the projected opacities.

    k_projection = 0
    DO WHILE &
      ( ANY( ProjectionMask(:) ) &
        .AND. k_projection < MaxIter_ConservativeProjection )

      k_projection = k_projection + 1

      CALL TimersStart( Timer_Collisions_MatterRHS )

      CALL ComputeMatterRHS_OrderV &
             ( ProjectionMask, &
               FVECm_projection, GVECm_projection, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1_Conservative, V_u_2_Conservative, &
               V_u_3_Conservative, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL TimersStop( Timer_Collisions_MatterRHS )

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( Ef, ProjectionResidual, V_d_1, V_d_2, V_d_3, vDotV )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( Ef, ProjectionResidual, V_d_1, V_d_2, V_d_3, vDotV )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( Ef, ProjectionResidual, V_d_1, V_d_2, V_d_3, vDotV )
#endif
      DO iN_X = 1, nX_G

        IF( ProjectionMask(iN_X) )THEN

          V_d_1 = Gm_dd_11(iN_X) * V_u_1_Conservative(iN_X)
          V_d_2 = Gm_dd_22(iN_X) * V_u_2_Conservative(iN_X)
          V_d_3 = Gm_dd_33(iN_X) * V_u_3_Conservative(iN_X)
          vDotV =   V_u_1_Conservative(iN_X) * V_d_1 &
                  + V_u_2_Conservative(iN_X) * V_d_2 &
                  + V_u_3_Conservative(iN_X) * V_d_3
          Ef = E_Conservative(iN_X) + Half * vDotV

          ProjectionResidual &
            = MAX &
                ( ABS &
                    ( GVECm_projection(iY,iN_X) &
                      - Y_Conservative(iN_X) / Y_old(iN_X) ), &
                  ABS &
                    ( GVECm_projection(iEf,iN_X) &
                      - Ef / Ef_old(iN_X) ), &
                  ABS &
                    ( GVECm_projection(iV1,iN_X) &
                      - V_d_1 / SpeedOfLight ), &
                  ABS &
                    ( GVECm_projection(iV2,iN_X) &
                      - V_d_2 / SpeedOfLight ), &
                  ABS &
                    ( GVECm_projection(iV3,iN_X) &
                      - V_d_3 / SpeedOfLight ) )

          Y_Conservative(iN_X) &
            = GVECm_projection(iY,iN_X) * Y_old(iN_X)
          Ef = GVECm_projection(iEf,iN_X) * Ef_old(iN_X)

          V_d_1 = GVECm_projection(iV1,iN_X) * SpeedOfLight
          V_d_2 = GVECm_projection(iV2,iN_X) * SpeedOfLight
          V_d_3 = GVECm_projection(iV3,iN_X) * SpeedOfLight

          V_u_1_Conservative(iN_X) = V_d_1 / Gm_dd_11(iN_X)
          V_u_2_Conservative(iN_X) = V_d_2 / Gm_dd_22(iN_X)
          V_u_3_Conservative(iN_X) = V_d_3 / Gm_dd_33(iN_X)

          vDotV =   V_u_1_Conservative(iN_X) * V_d_1 &
                  + V_u_2_Conservative(iN_X) * V_d_2 &
                  + V_u_3_Conservative(iN_X) * V_d_3

          E_Conservative(iN_X) = Ef - Half * vDotV

          IF( ProjectionResidual <= Rtol_ConservativeProjection ) &
            ProjectionMask(iN_X) = .FALSE.

        END IF

      END DO

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( ProjectionMask )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( ProjectionMask )
#endif

    END DO

    IF( ANY( ProjectionMask(:) ) )THEN

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        IF( ProjectionMask(iN_X) )THEN
          ConservativeProjectionRadiationResidual_Relative_Final(iN_X) &
            = HUGE(One)
          ConservativeProjectionStatus_Final(iN_X) = -3
          ProjectionError(iN_X) = 1
        END IF

      END DO

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      ProjectionMask(iN_X) &
        = ConservativeProjectionStatus_Final(iN_X) == 0

    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( ProjectionMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( ProjectionMask )
#endif

    CALL CreatePackIndex &
           ( ProjectionMask, nX_P, PackIndex, UnpackIndex )

    IF( nX_P > 0 )THEN

      CALL UpdateTemperature_Packed &
             ( D, E_Conservative, Y_Conservative, T_Conservative, &
               ProjectionMask, nX_P, PackIndex, UnpackIndex, ProjectionError )

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        IF( ProjectionMask(iN_X) .AND. ProjectionError(iN_X) /= 0 )THEN
          ProjectionMask(iN_X) = .FALSE.
          ConservativeProjectionRadiationResidual_Relative_Final(iN_X) &
            = HUGE(One)
          ConservativeProjectionStatus_Final(iN_X) = -1
        END IF

      END DO

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( ProjectionMask )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( ProjectionMask )
#endif

      CALL CreatePackIndex &
             ( ProjectionMask, nX_P, PackIndex, UnpackIndex )

    END IF

    ! Heating implied by the exactly conservative fixed-radiation matter
    ! state.  This uses the same total-energy-minus-velocity-times-momentum
    ! definition as MatterIncrement_Rate_Option.

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Ef, Q_E, Q_M1, Q_M2, Q_M3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotV )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( Ef, Q_E, Q_M1, Q_M2, Q_M3, &
    !$ACC          V_d_1, V_d_2, V_d_3, vDotV )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Ef, Q_E, Q_M1, Q_M2, Q_M3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotV )
#endif
    DO iN_X = 1, nX_G

      IF( ProjectionMask(iN_X) )THEN

        V_d_1 = Gm_dd_11(iN_X) * V_u_1_Conservative(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2_Conservative(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3_Conservative(iN_X)

        vDotV =   V_u_1_Conservative(iN_X) * V_d_1 &
                + V_u_2_Conservative(iN_X) * V_d_2 &
                + V_u_3_Conservative(iN_X) * V_d_3
        Ef = E_Conservative(iN_X) + Half * vDotV

        Q_E = D_old(iN_X) * ( Ef - Ef_old(iN_X) ) / dt
        Q_M1 = D_old(iN_X) &
          * ( V_d_1 - Gm_dd_11(iN_X) * V_u_1_old(iN_X) ) / dt
        Q_M2 = D_old(iN_X) &
          * ( V_d_2 - Gm_dd_22(iN_X) * V_u_2_old(iN_X) ) / dt
        Q_M3 = D_old(iN_X) &
          * ( V_d_3 - Gm_dd_33(iN_X) * V_u_3_old(iN_X) ) / dt

        ConservativeProjectionMatterHeating_Rate_Final(iN_X) &
          = Q_E &
              - V_u_1_Conservative(iN_X) * Q_M1 &
              - V_u_2_Conservative(iN_X) * Q_M2 &
              - V_u_3_Conservative(iN_X) * Q_M3

        T(iN_X) = T_Conservative(iN_X)
        Y(iN_X) = Y_Conservative(iN_X)
        E(iN_X) = E_Conservative(iN_X)
        V_u_1(iN_X) = V_u_1_Conservative(iN_X)
        V_u_2(iN_X) = V_u_2_Conservative(iN_X)
        V_u_3(iN_X) = V_u_3_Conservative(iN_X)

        U_Y(iN_X) = Y_Conservative(iN_X) / Y_old(iN_X)
        U_Ef(iN_X) = Ef / Ef_old(iN_X)
        U_V_d_1(iN_X) = V_d_1 / SpeedOfLight
        U_V_d_2(iN_X) = V_d_2 / SpeedOfLight
        U_V_d_3(iN_X) = V_d_3 / SpeedOfLight
        Omega(iN_X) = One / ( One + SQRT( vDotV ) )

      END IF

    END DO

    IF( nX_P > 0 )THEN

      IF( .NOT. FreezeOpacities )THEN

        CALL TimersStart( Timer_Collisions_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( D, T_Conservative, Y_Conservative, SqrtGm, &
                 ProjectionMask, nX_P, PackIndex, UnpackIndex, &
                 T0 = T_old, Y0 = Y_old, &
                 Opacity_Mask_Input = Opacity_Mask )

        CALL TimersStop( Timer_Collisions_ComputeOpacity )

      END IF

      CALL TimersStart( Timer_Collisions_ComputeRates )

      CALL ComputeRates_Packed &
             ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               ProjectionMask, nX_P, PackIndex, UnpackIndex, nX_P, Opacity_Mask )

      CALL TimersStop( Timer_Collisions_ComputeRates )

      CALL TimersStart( Timer_Collisions_NeutrinoRHS )

      CALL ComputeNeutrinoRHS_OrderV &
             ( ProjectionMask, FVECm_inner, GVECm_inner, dt, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1_Conservative, V_u_2_Conservative, V_u_3_Conservative, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL TimersStop( Timer_Collisions_NeutrinoRHS )

      CALL ComputeRadiationResidualRelative &
             ( ProjectionMask, FVECm_inner, &
               ConservativeProjectionRadiationResidual_Relative_Final, &
               ContinueIteration_Option = ContinuationMask )

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        IF( ProjectionMask(iN_X) )THEN
          ConservativeProjectionStatus_Final(iN_X) &
            = MERGE( 0, 1, ContinuationMask(iN_X) )
        END IF

      END DO

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( ContinuationMask, ProjectionError )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( ContinuationMask, ProjectionError )
#endif

  END SUBROUTINE ApplyConservativeProjection_OrderV


  SUBROUTINE ComputeRadiationResidualRelative &
    ( MASK, Fm, RelativeResidual, ContinueIteration_Option )

    LOGICAL,  DIMENSION(:),   INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: Fm
    REAL(DP), DIMENSION(:),   INTENT(inout) :: RelativeResidual
    LOGICAL,  DIMENSION(:),   INTENT(inout), OPTIONAL :: &
      ContinueIteration_Option

    LOGICAL  :: Converged
    INTEGER  :: iN_E, iN_X, iOS, iS
    REAL(DP) :: Fnorm, Fnorm_G1, Fnorm_G2, Fnorm_G3, Fnorm_N
    REAL(DP) :: SpeciesRelativeResidual

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Converged, iN_E, iOS, iS, Fnorm, Fnorm_G1, Fnorm_G2, &
    !$OMP          Fnorm_G3, Fnorm_N, SpeciesRelativeResidual )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( Converged, iN_E, iOS, iS, Fnorm, Fnorm_G1, Fnorm_G2, &
    !$ACC          Fnorm_G3, Fnorm_N, SpeciesRelativeResidual )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Converged, iN_E, iOS, iS, Fnorm, Fnorm_G1, Fnorm_G2, &
    !$OMP          Fnorm_G3, Fnorm_N, SpeciesRelativeResidual )
#endif
    DO iN_X = 1, nX_G

      IF( MASK(iN_X) )THEN

        RelativeResidual(iN_X) = Zero
        Converged = .TRUE.

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

          IF( DnuNorm(iS,iN_X) > Zero )THEN
            SpeciesRelativeResidual = Fnorm / DnuNorm(iS,iN_X)
          ELSE IF( Fnorm == Zero )THEN
            SpeciesRelativeResidual = Zero
          ELSE
            SpeciesRelativeResidual = HUGE(One)
          END IF

          RelativeResidual(iN_X) &
            = MAX( RelativeResidual(iN_X), SpeciesRelativeResidual )

          Converged = Converged .AND. &
            ( ( Fnorm <= Rtol_inner * DnuNorm(iS,iN_X) ) .OR. &
              ( Fnorm <= Atol_inner * Fnorm0(iS,iN_X) ) )

        END DO

        IF( PRESENT( ContinueIteration_Option ) ) &
          ContinueIteration_Option(iN_X) = .NOT. Converged

      END IF

    END DO

  END SUBROUTINE ComputeRadiationResidualRelative
#endif


  SUBROUTINE ComputeOpacities_Packed &
    ( D, T, Y, SqrtGm, MASK, nX_P, PackIndex, &
      UnpackIndex, nX_P0, T0, Y0, Opacity_Mask_Input )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET   :: D, T, Y, SqrtGm
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: MASK
    INTEGER,                INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                INTENT(in), OPTIONAL :: nX_P0
    REAL(DP), DIMENSION(:), INTENT(in), TARGET, OPTIONAL :: T0, Y0
    LOGICAL, DIMENSION(:,:), INTENT(in), TARGET :: Opacity_Mask_Input

    INTEGER                             :: nX, nX0, iX, iE, iC
    REAL(DP), DIMENSION(:)    , POINTER :: D_P, T_P, Y_P, SqrtGm_P
    REAL(DP), DIMENSION(:)    , POINTER :: T0_P, Y0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Dnu_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Sigma_Iso_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_Iso_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_1_Iso_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_EmAb_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_NNS_0_P, Phi_NbNS_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_NNS_1_P, Phi_NbNS_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_0_P, H_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_1_P, H_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_0_P, J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_1_P, J_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Nu_J_I_0_P, Nu_J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Nu_J_I_1_P, Nu_J_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: S_Sigma_P

    LOGICAL,  DIMENSION(:,:),   POINTER :: Opacity_Mask_P

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

      Opacity_Mask_P => Opacity_Mask_T(:,1:nX)
  
      CALL ArrayPack &
             ( nX, UnpackIndex, Opacity_Mask_Input, Opacity_Mask_P )

      D_P => D_T(1:nX)
      T_P => T_T(1:nX)
      Y_P => Y_T(1:nX)

      IF ( PRESENT ( T0 ) ) THEN
        T0_P => T0_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, T0, T0_P ) 
      ELSE
        T0_P => T_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, T, T0_P ) 
      ENDIF

      IF ( PRESENT ( Y0 ) ) THEN
        Y0_P => Y0_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, Y0, Y0_P ) 
      ELSE
        Y0_P => Y_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, Y, Y0_P ) 
      ENDIF

      SqrtGm_P => SqrtGm_T(1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, D, T, Y, SqrtGm, D_P, T_P, Y_P, SqrtGm_P )

      Dnu_0_P     => Dnu_0_T    (:,:,1:nX)
      Sigma_Iso_P => Sigma_Iso_T(:,:,1:nX)
      Phi_0_Iso_P => Phi_0_Iso_T(:,:,1:nX)
      Phi_1_Iso_P => Phi_1_Iso_T(:,:,1:nX)
      Chi_EmAb_P  => Chi_EmAb_T (:,:,1:nX)

      Phi_NNS_0_P   => Phi_NNS_0_T  (:,:,1:nX)
      Phi_NbNS_0_P  => Phi_NbNS_0_T (:,:,1:nX)

      H_I_0_P     => H_I_0_T    (:,:,1:nX)
      H_II_0_P    => H_II_0_T   (:,:,1:nX)

      J_I_0_P     => J_I_0_T    (:,:,1:nX)
      J_II_0_P    => J_II_0_T   (:,:,1:nX)

      Nu_J_I_0_P  => Nu_J_I_0_T (:,:,1:nX)
      Nu_J_II_0_P => Nu_J_II_0_T(:,:,1:nX)

      IF ( UseFirstMomentRates ) THEN
        Phi_NNS_1_P  => Phi_NNS_1_T (:,:,1:nX)
        Phi_NbNS_1_P => Phi_NbNS_1_T(:,:,1:nX)
        H_I_1_P      => H_I_1_T     (:,:,1:nX)
        H_II_1_P     => H_II_1_T    (:,:,1:nX)
        J_I_1_P      => J_I_1_T     (:,:,1:nX)
        J_II_1_P     => J_II_1_T    (:,:,1:nX)
        Nu_J_I_1_P   => Nu_J_I_1_T  (:,:,1:nX)
        Nu_J_II_1_P  => Nu_J_II_1_T (:,:,1:nX)
      END IF

      S_Sigma_P   => S_Sigma_T  (:,:,1:nX)

    ELSE

      Opacity_Mask_P => Opacity_Mask_Input(:,:)

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      IF ( PRESENT ( T0 ) ) THEN
        T0_P => T0(:)
      ELSE
        T0_P => T(:)
      ENDIF

      IF ( PRESENT ( Y0 ) ) THEN
        Y0_P => Y0(:)
      ELSE
        Y0_P => Y(:)
      ENDIF

      SqrtGm_P => SqrtGm(:)

      Dnu_0_P     => Dnu_0    (:,:,:)
      Sigma_Iso_P => Sigma_Iso(:,:,:)
      Phi_0_Iso_P => Phi_0_Iso(:,:,:)
      Phi_1_Iso_P => Phi_1_Iso(:,:,:)
      Chi_EmAb_P  => Chi_EmAb (:,:,:)

      Phi_NNS_0_P   => Phi_NNS_0  (:,:,:)
      Phi_NbNS_0_P  => Phi_NbNS_0 (:,:,:)

      H_I_0_P     => H_I_0    (:,:,:)
      H_II_0_P    => H_II_0   (:,:,:)

      J_I_0_P     => J_I_0    (:,:,:)
      J_II_0_P    => J_II_0   (:,:,:)

      Nu_J_I_0_P  => Nu_J_I_0 (:,:,:)
      Nu_J_II_0_P => Nu_J_II_0(:,:,:)

      IF ( UseFirstMomentRates ) THEN
        Phi_NNS_1_P  => Phi_NNS_1 (:,:,:)
        Phi_NbNS_1_P => Phi_NbNS_1(:,:,:)
        H_I_1_P      => H_I_1     (:,:,:)
        H_II_1_P     => H_II_1    (:,:,:)
        J_I_1_P      => J_I_1     (:,:,:)
        J_II_1_P     => J_II_1    (:,:,:)
        Nu_J_I_1_P   => Nu_J_I_1  (:,:,:)
        Nu_J_II_1_P  => Nu_J_II_1 (:,:,:)
      END IF

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
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, Dnu_0_P, DnuMax )

    CALL TimersStop( Timer_Opacity_LimitD0 )

    ! --- EmAb ---

    CALL TimersStart( Timer_Opacity_EC )

    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Y_P, Dnu_0_P, Chi_EmAb_P, &
             Opacity_Mask_P, T_old = T0_P, Y_old = Y0_P )

    !Chi_EmAb_P = Zero

    CALL TimersStop( Timer_Opacity_EC )

    ! --- Isoenergetic scattering ---

    CALL TimersStart( Timer_Opacity_ES )

    CALL ComputeNeutrinoOpacities_ES &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, &
             iM0, Phi_0_Iso_P, Opacity_Mask_P )

    CALL ComputeNeutrinoOpacities_ES &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, &
             iM1, Phi_1_Iso_P, Opacity_Mask_P )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = 1, nX
    DO iC = 1, nChirals
    DO iE = 1, nE_G

      Sigma_Iso_P(iE,iC,iX) &
        = FourPiEp2(iE) * ( Phi_0_Iso_P(iE,iC,iX) - Third * Phi_1_Iso_P(iE,iC,iX) )

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_Opacity_ES )

    IF( Include_NNS )THEN

      ! --- NNS Scattering Functions ---

      CALL TimersStart( Timer_Opacity_NNS )

      CALL ComputeNeutrinoOpacities_NNS &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM0, &
               Phi_NNS_0_P, Phi_NbNS_0_P, &
               ApplyWeakMagnetism = NNS_ApplyWeakMagnetism, &
               ApplyManyBodyCorrection = NNS_ApplyManyBodyCorrection )

      IF ( UseFirstMomentRates ) &
        CALL ComputeNeutrinoOpacities_NNS &
               ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM1, &
                 Phi_NNS_1_P, Phi_NbNS_1_P, &
                 ApplyWeakMagnetism = NNS_ApplyWeakMagnetism, &
                 ApplyManyBodyCorrection = NNS_ApplyManyBodyCorrection )

      CALL TimersStop( Timer_Opacity_NNS )

    END IF

    IF( Include_NES )THEN

      ! --- NES Scattering Functions ---

      CALL TimersStart( Timer_Opacity_NES )

      CALL ComputeNeutrinoOpacities_NES &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM0, H_I_0_P, H_II_0_P )

      IF( UseFirstMomentRates )THEN

        CALL ComputeNeutrinoOpacities_NES &
               ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM1, H_I_1_P, H_II_1_P )

      END IF

      CALL TimersStop( Timer_Opacity_NES )

    END IF

    IF( Include_Pair )THEN

      ! --- Pair Kernels ---

      CALL TimersStart( Timer_Opacity_Pair )

      CALL ComputeNeutrinoOpacities_Pair &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM0, J_I_0_P, J_II_0_P )

      IF( UseFirstMomentRates )THEN

        CALL ComputeNeutrinoOpacities_Pair &
               ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM1, J_I_1_P, J_II_1_P )

      END IF

      CALL TimersStop( Timer_Opacity_Pair )

    END IF

    IF( Include_NuPair )THEN

      ! --- NuPair Kernels ---

      IF( ANY( Opacity_Mask_Input(iOp_NuPair,:) ) ) THEN

        CALL TimersStart( Timer_Opacity_NuPair )

        CALL ComputeNeutrinoOpacities_NuPair &
               ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM0, &
                 Nu_J_I_0_P, Nu_J_II_0_P )

        IF ( UseFirstMomentRates ) &
          CALL ComputeNeutrinoOpacities_NuPair &
                 ( 1, nE_G, 1, nX, D_P, T_P, Y_P, iM1, &
                   Nu_J_I_1_P, Nu_J_II_1_P )

        CALL TimersStop( Timer_Opacity_NuPair )

      END IF

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
               ( nX, MASK, PackIndex, Phi_NNS_0_P, Phi_NNS_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, Phi_NbNS_0_P, Phi_NbNS_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, H_I_0_P, H_II_0_P, H_I_0, H_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, J_I_0_P, J_II_0_P, J_I_0, J_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, Nu_J_I_0_P, Nu_J_II_0_P, Nu_J_I_0, Nu_J_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, S_Sigma_P, S_Sigma )

        IF( UseFirstMomentRates )THEN

          CALL ArrayUnpack &
                 ( nX, MASK, PackIndex, &
                   Phi_NNS_1_P, Phi_NbNS_1_P, Phi_NNS_1, Phi_NbNS_1 )

          CALL ArrayUnpack &
                 ( nX, MASK, PackIndex, H_I_1_P, H_II_1_P, H_I_1, H_II_1 )

          CALL ArrayUnpack &
                 ( nX, MASK, PackIndex, J_I_1_P, J_II_1_P, J_I_1, J_II_1 )

          CALL ArrayUnpack &
                 ( nX, MASK, PackIndex, &
                   Nu_J_I_1_P, Nu_J_II_1_P, Nu_J_I_1, Nu_J_II_1 )

        END IF

      END IF

    END IF

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( D, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, MASK, nX_P, PackIndex, UnpackIndex, &
      nX_P0, Opacity_Mask_Input )

    REAL(DP), DIMENSION(:),     INTENT(in), TARGET   :: D
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Dnu
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_1
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_2
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_3
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: MASK
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:),     INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P0
    LOGICAL,  DIMENSION(:,:),   INTENT(in), TARGET   :: Opacity_Mask_Input

    REAL(DP), DIMENSION(:)    , POINTER :: D_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Dnu_P, Dnu_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Inu_u_1_P, Inu_u_2_P, Inu_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_NNS_P , Eta_NNS_P, R_NNS_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_NES_P , Eta_NES_P, R_NES_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_Pair_P, Eta_Pair_P, R_Pair_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_NuPair_P, Eta_NuPair_P, R_NuPair_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_Brem_P, Eta_Brem_P, R_Brem_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NNS__In__u_1_P, L_NNS__Out_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NNS__In__u_2_P, L_NNS__Out_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NNS__In__u_3_P, L_NNS__Out_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NES__In__u_1_P, L_NES__Out_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NES__In__u_2_P, L_NES__Out_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NES__In__u_3_P, L_NES__Out_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Pair_Pro_u_1_P, L_Pair_Ann_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Pair_Pro_u_2_P, L_Pair_Ann_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Pair_Pro_u_3_P, L_Pair_Ann_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NuPair_Pro_u_1_P, L_NuPair_Ann_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NuPair_Pro_u_2_P, L_NuPair_Ann_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_NuPair_Pro_u_3_P, L_NuPair_Ann_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Brem_Pro_u_1_P, L_Brem_Ann_u_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Brem_Pro_u_2_P, L_Brem_Ann_u_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: L_Brem_Pro_u_3_P, L_Brem_Ann_u_3_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_NNS_0_P, Phi_NbNS_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_NNS_1_P, Phi_NbNS_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_0_P, H_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_1_P, H_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_0_P, J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_1_P, J_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Nu_J_I_0_P, Nu_J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Nu_J_I_1_P, Nu_J_II_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: S_Sigma_P

    LOGICAL,  DIMENSION(:,:),   POINTER :: Opacity_Mask_P

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

      Opacity_Mask_P => Opacity_Mask_T(:,1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, Opacity_Mask_Input, Opacity_Mask_P )

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

      Chi_NNS_P        => Chi_NNS_T       (:,:,1:nX)
      Eta_NNS_P        => Eta_NNS_T       (:,:,1:nX)
      R_NNS_P          => R_NNS_T         (:,:,1:nX)
      Chi_NES_P        => Chi_NES_T       (:,:,1:nX)
      Eta_NES_P        => Eta_NES_T       (:,:,1:nX)
      R_NES_P          => R_NES_T         (:,:,1:nX)
      Chi_Pair_P       => Chi_Pair_T      (:,:,1:nX)
      Eta_Pair_P       => Eta_Pair_T      (:,:,1:nX)
      R_Pair_P         => R_Pair_T        (:,:,1:nX)
      Chi_NuPair_P     => Chi_NuPair_T    (:,:,1:nX)
      Eta_NuPair_P     => Eta_NuPair_T    (:,:,1:nX)
      R_NuPair_P       => R_NuPair_T      (:,:,1:nX)
      Chi_Brem_P       => Chi_Brem_T      (:,:,1:nX)
      Eta_Brem_P       => Eta_Brem_T      (:,:,1:nX)
      R_Brem_P         => R_Brem_T        (:,:,1:nX)

      Phi_NNS_0_P        => Phi_NNS_0_T       (:,:,1:nX)
      Phi_NbNS_0_P       => Phi_NbNS_0_T      (:,:,1:nX)

      H_I_0_P          => H_I_0_T         (:,:,1:nX)
      H_II_0_P         => H_II_0_T        (:,:,1:nX)

      J_I_0_P          => J_I_0_T         (:,:,1:nX)
      J_II_0_P         => J_II_0_T        (:,:,1:nX)

      Nu_J_I_0_P       => Nu_J_I_0_T      (:,:,1:nX)
      Nu_J_II_0_P      => Nu_J_II_0_T     (:,:,1:nX)

      IF ( UseFirstMomentRates ) THEN
        L_NNS__In__u_1_P  => L_NNS__In__u_1_T (:,:,1:nX)
        L_NNS__In__u_2_P  => L_NNS__In__u_2_T (:,:,1:nX)
        L_NNS__In__u_3_P  => L_NNS__In__u_3_T (:,:,1:nX)
        L_NNS__Out_u_1_P  => L_NNS__Out_u_1_T (:,:,1:nX)
        L_NNS__Out_u_2_P  => L_NNS__Out_u_2_T (:,:,1:nX)
        L_NNS__Out_u_3_P  => L_NNS__Out_u_3_T (:,:,1:nX)
        L_NES__In__u_1_P  => L_NES__In__u_1_T (:,:,1:nX)
        L_NES__In__u_2_P  => L_NES__In__u_2_T (:,:,1:nX)
        L_NES__In__u_3_P  => L_NES__In__u_3_T (:,:,1:nX)
        L_NES__Out_u_1_P  => L_NES__Out_u_1_T (:,:,1:nX)
        L_NES__Out_u_2_P  => L_NES__Out_u_2_T (:,:,1:nX)
        L_NES__Out_u_3_P  => L_NES__Out_u_3_T (:,:,1:nX)
        L_Pair_Pro_u_1_P  => L_Pair_Pro_u_1_T (:,:,1:nX)
        L_Pair_Pro_u_2_P  => L_Pair_Pro_u_2_T (:,:,1:nX)
        L_Pair_Pro_u_3_P  => L_Pair_Pro_u_3_T (:,:,1:nX)
        L_Pair_Ann_u_1_P  => L_Pair_Ann_u_1_T (:,:,1:nX)
        L_Pair_Ann_u_2_P  => L_Pair_Ann_u_2_T (:,:,1:nX)
        L_Pair_Ann_u_3_P  => L_Pair_Ann_u_3_T (:,:,1:nX)
        L_NuPair_Pro_u_1_P => L_NuPair_Pro_u_1_T(:,:,1:nX)
        L_NuPair_Pro_u_2_P => L_NuPair_Pro_u_2_T(:,:,1:nX)
        L_NuPair_Pro_u_3_P => L_NuPair_Pro_u_3_T(:,:,1:nX)
        L_NuPair_Ann_u_1_P => L_NuPair_Ann_u_1_T(:,:,1:nX)
        L_NuPair_Ann_u_2_P => L_NuPair_Ann_u_2_T(:,:,1:nX)
        L_NuPair_Ann_u_3_P => L_NuPair_Ann_u_3_T(:,:,1:nX)
        L_Brem_Pro_u_1_P  => L_Brem_Pro_u_1_T (:,:,1:nX)
        L_Brem_Pro_u_2_P  => L_Brem_Pro_u_2_T (:,:,1:nX)
        L_Brem_Pro_u_3_P  => L_Brem_Pro_u_3_T (:,:,1:nX)
        L_Brem_Ann_u_1_P  => L_Brem_Ann_u_1_T (:,:,1:nX)
        L_Brem_Ann_u_2_P  => L_Brem_Ann_u_2_T (:,:,1:nX)
        L_Brem_Ann_u_3_P  => L_Brem_Ann_u_3_T (:,:,1:nX)
        Phi_NNS_1_P       => Phi_NNS_1_T      (:,:,1:nX)
        Phi_NbNS_1_P      => Phi_NbNS_1_T     (:,:,1:nX)
        H_I_1_P           => H_I_1_T          (:,:,1:nX)
        H_II_1_P          => H_II_1_T         (:,:,1:nX)
        J_I_1_P           => J_I_1_T          (:,:,1:nX)
        J_II_1_P          => J_II_1_T         (:,:,1:nX)
        Nu_J_I_1_P        => Nu_J_I_1_T       (:,:,1:nX)
        Nu_J_II_1_P       => Nu_J_II_1_T      (:,:,1:nX)
      END IF

      S_Sigma_P        => S_Sigma_T       (:,:,1:nX)

      IF ( nX < nX0 .OR. FreezeOpacities ) THEN

        CALL ArrayPack &
               ( nX, UnpackIndex, Phi_NNS_0, Phi_NNS_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, Phi_NbNS_0, Phi_NbNS_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, H_I_0, H_II_0, H_I_0_P, H_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, J_I_0, J_II_0, J_I_0_P, J_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, Nu_J_I_0, Nu_J_II_0, Nu_J_I_0_P, Nu_J_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, S_Sigma, S_Sigma_P )

        IF( UseFirstMomentRates )THEN

          CALL ArrayPack &
                 ( nX, UnpackIndex, &
                   Phi_NNS_1, Phi_NbNS_1, Phi_NNS_1_P, Phi_NbNS_1_P )

          CALL ArrayPack &
                 ( nX, UnpackIndex, H_I_1, H_II_1, H_I_1_P, H_II_1_P )

          CALL ArrayPack &
                 ( nX, UnpackIndex, J_I_1, J_II_1, J_I_1_P, J_II_1_P )

          CALL ArrayPack &
                 ( nX, UnpackIndex, &
                   Nu_J_I_1, Nu_J_II_1, Nu_J_I_1_P, Nu_J_II_1_P )

        END IF

      END IF

    ELSE

      Opacity_Mask_P => Opacity_Mask_Input(:,1:nX)

      D_P => D(:)

      Dnu_P            => Dnu           (:,:,:)
      Dnu_0_P          => Dnu_0         (:,:,:)
      Inu_u_1_P        => Inu_u_1       (:,:,:)
      Inu_u_2_P        => Inu_u_2       (:,:,:)
      Inu_u_3_P        => Inu_u_3       (:,:,:)

      Chi_NNS_P        => Chi_NNS       (:,:,:)
      Eta_NNS_P        => Eta_NNS       (:,:,:)
      R_NNS_P          => R_NNS         (:,:,:)
      Chi_NES_P        => Chi_NES       (:,:,:)
      Eta_NES_P        => Eta_NES       (:,:,:)
      R_NES_P          => R_NES         (:,:,:)
      Chi_Pair_P       => Chi_Pair      (:,:,:)
      Eta_Pair_P       => Eta_Pair      (:,:,:)
      R_Pair_P         => R_Pair        (:,:,:)
      Chi_NuPair_P     => Chi_NuPair    (:,:,:)
      Eta_NuPair_P     => Eta_NuPair    (:,:,:)
      R_NuPair_P       => R_NuPair      (:,:,:)
      Chi_Brem_P       => Chi_Brem      (:,:,:)
      Eta_Brem_P       => Eta_Brem      (:,:,:)
      R_Brem_P         => R_Brem        (:,:,:)

      IF ( UseFirstMomentRates ) THEN
      L_NNS__In__u_1_P => L_NNS__In__u_1(:,:,:)
      L_NNS__In__u_2_P => L_NNS__In__u_2(:,:,:)
      L_NNS__In__u_3_P => L_NNS__In__u_3(:,:,:)
      L_NNS__Out_u_1_P => L_NNS__Out_u_1(:,:,:)
      L_NNS__Out_u_2_P => L_NNS__Out_u_2(:,:,:)
      L_NNS__Out_u_3_P => L_NNS__Out_u_3(:,:,:)
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
      L_NuPair_Pro_u_1_P => L_NuPair_Pro_u_1(:,:,:)
      L_NuPair_Pro_u_2_P => L_NuPair_Pro_u_2(:,:,:)
      L_NuPair_Pro_u_3_P => L_NuPair_Pro_u_3(:,:,:)
      L_NuPair_Ann_u_1_P => L_NuPair_Ann_u_1(:,:,:)
      L_NuPair_Ann_u_2_P => L_NuPair_Ann_u_2(:,:,:)
      L_NuPair_Ann_u_3_P => L_NuPair_Ann_u_3(:,:,:)
      L_Brem_Pro_u_1_P => L_Brem_Pro_u_1(:,:,:)
      L_Brem_Pro_u_2_P => L_Brem_Pro_u_2(:,:,:)
      L_Brem_Pro_u_3_P => L_Brem_Pro_u_3(:,:,:)
      L_Brem_Ann_u_1_P => L_Brem_Ann_u_1(:,:,:)
      L_Brem_Ann_u_2_P => L_Brem_Ann_u_2(:,:,:)
      L_Brem_Ann_u_3_P => L_Brem_Ann_u_3(:,:,:)
      END IF

      Phi_NNS_0_P        => Phi_NNS_0       (:,:,:)
      Phi_NbNS_0_P       => Phi_NbNS_0      (:,:,:)

      H_I_0_P          => H_I_0         (:,:,:)
      H_II_0_P         => H_II_0        (:,:,:)

      J_I_0_P          => J_I_0         (:,:,:)
      J_II_0_P         => J_II_0        (:,:,:)

      Nu_J_I_0_P       => Nu_J_I_0      (:,:,:)
      Nu_J_II_0_P      => Nu_J_II_0     (:,:,:)

      IF ( UseFirstMomentRates ) THEN
        Phi_NNS_1_P  => Phi_NNS_1 (:,:,:)
        Phi_NbNS_1_P => Phi_NbNS_1(:,:,:)
        H_I_1_P      => H_I_1     (:,:,:)
        H_II_1_P     => H_II_1    (:,:,:)
        J_I_1_P      => J_I_1     (:,:,:)
        J_II_1_P     => J_II_1    (:,:,:)
        Nu_J_I_1_P   => Nu_J_I_1  (:,:,:)
        Nu_J_II_1_P  => Nu_J_II_1 (:,:,:)
      END IF

      S_Sigma_P        => S_Sigma       (:,:,:)

    END IF

    ! --- NNS Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NNS )

    CALL ComputeNeutrinoOpacityRates_NNS &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, &
             Phi_NNS_0_P, Phi_NbNS_0_P, Eta_NNS_P, Chi_NNS_P, R_NNS_P, &
             Opacity_Mask_P )

    IF ( Include_NNS .AND. UseFirstMomentRates ) &
      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NNS &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, &
               Phi_NNS_1_P, Phi_NbNS_1_P, &
               L_NNS__In__u_1_P, L_NNS__In__u_2_P, L_NNS__In__u_3_P, &
               L_NNS__Out_u_1_P, L_NNS__Out_u_2_P, L_NNS__Out_u_3_P, &
               Opacity_Mask_P )

    CALL TimersStop( Timer_OpacityRate_NNS )

    ! --- NES Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NES )

    CALL ComputeNeutrinoOpacityRates_NES &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, &
             H_I_0_P, H_II_0_P, Eta_NES_P, Chi_NES_P, R_NES_P, &
             Opacity_Mask_P )

    IF( Include_NES .AND. UseFirstMomentRates )THEN

      ! --- Compute Linear Rate Corrections (NES) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NES &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, H_I_1_P, H_II_1_P, &
               L_NES__In__u_1_P, L_NES__In__u_2_P, L_NES__In__u_3_P, &
               L_NES__Out_u_1_P, L_NES__Out_u_2_P, L_NES__Out_u_3_P, &
               Opacity_Mask_P )

    END IF

    CALL TimersStop( Timer_OpacityRate_NES )

    ! --- Pair Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_Pair )

    CALL ComputeNeutrinoOpacityRates_Pair &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, &
             J_I_0_P, J_II_0_P, Eta_Pair_P, Chi_Pair_P, R_Pair_P, &
             Opacity_Mask_P )

    IF( Include_Pair .AND. UseFirstMomentRates )THEN

      ! --- Compute Linear Rate Corrections (Pair) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, J_I_1_P, J_II_1_P, &
               L_Pair_Pro_u_1_P, L_Pair_Pro_u_2_P, L_Pair_Pro_u_3_P, &
               L_Pair_Ann_u_1_P, L_Pair_Ann_u_2_P, L_Pair_Ann_u_3_P, &
               Opacity_Mask_P )

    END IF

    CALL TimersStop( Timer_OpacityRate_Pair )

    ! --- NuPair Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NuPair )

    CALL ComputeNeutrinoOpacityRates_NuPair &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, &
             Nu_J_I_0_P, Nu_J_II_0_P, Eta_NuPair_P, Chi_NuPair_P, &
             R_NuPair_P, Opacity_Mask_P )

    IF ( Include_NuPair .AND. UseFirstMomentRates ) &
      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NuPair &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, &
               Nu_J_I_1_P, Nu_J_II_1_P, &
               L_NuPair_Pro_u_1_P, L_NuPair_Pro_u_2_P, L_NuPair_Pro_u_3_P, &
               L_NuPair_Ann_u_1_P, L_NuPair_Ann_u_2_P, L_NuPair_Ann_u_3_P, &
               Opacity_Mask_P )

    CALL TimersStop( Timer_OpacityRate_NuPair )

    ! --- Brem Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_Brem )

    CALL ComputeNeutrinoOpacityRates_Brem &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, Dnu_P, Dnu_0_P, &
             S_Sigma_P, Eta_Brem_P, Chi_Brem_P, R_Brem_P, &
             Opacity_Mask_P )

    IF( Include_Brem .AND. UseFirstMomentRates )THEN

      ! --- Compute Linear Rate Corrections (Brem) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_Brem &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, Dnu_0_P, S_Sigma_P, &
               L_Brem_Pro_u_1_P, L_Brem_Pro_u_2_P, L_Brem_Pro_u_3_P, &
               L_Brem_Ann_u_1_P, L_Brem_Ann_u_2_P, L_Brem_Ann_u_3_P, &
               Opacity_Mask_P )

    END IF

    CALL TimersStop( Timer_OpacityRate_Brem )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NNS_P, Eta_NNS_P, Chi_NNS, Eta_NNS )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, R_NNS_P, R_NNS )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NES_P, Eta_NES_P, Chi_NES, Eta_NES )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, R_NES_P, R_NES )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_Pair_P, Eta_Pair_P, Chi_Pair, Eta_Pair )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, R_Pair_P, R_Pair )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NuPair_P, Eta_NuPair_P, Chi_NuPair, Eta_NuPair )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, R_NuPair_P, R_NuPair )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_Brem_P, Eta_Brem_P, Chi_Brem, Eta_Brem )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, R_Brem_P, R_Brem )

      IF( UseFirstMomentRates )THEN

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_NNS__In__u_1_P, L_NNS__In__u_2_P, L_NNS__In__u_3_P, &
                 L_NNS__In__u_1  , L_NNS__In__u_2  , L_NNS__In__u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_NNS__Out_u_1_P, L_NNS__Out_u_2_P, L_NNS__Out_u_3_P, &
                 L_NNS__Out_u_1  , L_NNS__Out_u_2  , L_NNS__Out_u_3 )

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
                 L_NuPair_Pro_u_1_P, L_NuPair_Pro_u_2_P, L_NuPair_Pro_u_3_P, &
                 L_NuPair_Pro_u_1  , L_NuPair_Pro_u_2  , L_NuPair_Pro_u_3 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 L_NuPair_Ann_u_1_P, L_NuPair_Ann_u_2_P, L_NuPair_Ann_u_3_P, &
                 L_NuPair_Ann_u_1  , L_NuPair_Ann_u_2  , L_NuPair_Ann_u_3 )

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


  SUBROUTINE ApplyNeutrinoDistributionLimit_OrderV &
    ( dt, MASK, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Nnu_Conserved, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      LimiterRate, ConservationCorrection, &
      LeptonFlavorInputConversionRate, LeptonFlavorLimiterRate )

    REAL(DP),                   INTENT(in)    :: dt
    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Nnu_Conserved
    REAL(DP), DIMENSION(:),     INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:),   INTENT(out)   :: LimiterRate
    REAL(DP), DIMENSION(:,:),   INTENT(out)   :: ConservationCorrection
    REAL(DP), DIMENSION(:,:),   INTENT(out)   :: &
      LeptonFlavorInputConversionRate
    REAL(DP), DIMENSION(:,:),   INTENT(out)   :: LeptonFlavorLimiterRate

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: Dnu_Input, DeltaDnu
    REAL(DP) :: Inu_u_1_Input, Inu_u_2_Input, Inu_u_3_Input
    REAL(DP) :: DeltaInu_u_1, DeltaInu_u_2, DeltaInu_u_3
    REAL(DP) :: DeltaLepton, DeltaLeptonElectron
    REAL(DP) :: DeltaLeptonMuon, DeltaLeptonTau, DeltaNumber
    REAL(DP) :: InputLeptonElectron, InputLeptonMuon, InputLeptonTau
    REAL(DP) :: DeltaEnergy, DeltaMomentum1, DeltaMomentum2, DeltaMomentum3
    REAL(DP) :: DeltaFnu_d_1, DeltaFnu_d_2, DeltaFnu_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: Nnu_Input, Nnu_Limited
    REAL(DP) :: Gnu_Input_d_1, Gnu_Input_d_2, Gnu_Input_d_3
    REAL(DP) :: Gnu_Limited_d_1, Gnu_Limited_d_2, Gnu_Limited_d_3

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE FIRSTPRIVATE( DnuMax ) &
    !$OMP PRIVATE( V_d_1, V_d_2, V_d_3, &
    !$OMP          DeltaLepton, DeltaLeptonElectron, &
    !$OMP          DeltaLeptonMuon, DeltaLeptonTau, DeltaNumber, DeltaEnergy, &
    !$OMP          InputLeptonElectron, InputLeptonMuon, InputLeptonTau, &
    !$OMP          DeltaMomentum1, DeltaMomentum2, DeltaMomentum3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG FIRSTPRIVATE( DnuMax ) &
    !$ACC PRIVATE( V_d_1, V_d_2, V_d_3, &
    !$ACC          DeltaLepton, DeltaLeptonElectron, &
    !$ACC          DeltaLeptonMuon, DeltaLeptonTau, DeltaNumber, DeltaEnergy, &
    !$ACC          InputLeptonElectron, InputLeptonMuon, InputLeptonTau, &
    !$ACC          DeltaMomentum1, DeltaMomentum2, DeltaMomentum3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO FIRSTPRIVATE( DnuMax ) &
    !$OMP PRIVATE( V_d_1, V_d_2, V_d_3, &
    !$OMP          DeltaLepton, DeltaLeptonElectron, &
    !$OMP          DeltaLeptonMuon, DeltaLeptonTau, DeltaNumber, DeltaEnergy, &
    !$OMP          InputLeptonElectron, InputLeptonMuon, InputLeptonTau, &
    !$OMP          DeltaMomentum1, DeltaMomentum2, DeltaMomentum3, &
    !$OMP          Dnu_Input, DeltaDnu, &
    !$OMP          Inu_u_1_Input, Inu_u_2_Input, Inu_u_3_Input, &
    !$OMP          DeltaInu_u_1, DeltaInu_u_2, DeltaInu_u_3, &
    !$OMP          DeltaFnu_d_1, DeltaFnu_d_2, DeltaFnu_d_3, &
    !$OMP          Nnu_Input, Nnu_Limited, &
    !$OMP          Gnu_Input_d_1, Gnu_Input_d_2, Gnu_Input_d_3, &
    !$OMP          Gnu_Limited_d_1, Gnu_Limited_d_2, Gnu_Limited_d_3 )
#endif
    DO iN_X = 1, nX_G

      DeltaLepton = Zero
      DeltaLeptonElectron = Zero
      DeltaLeptonMuon = Zero
      DeltaLeptonTau = Zero
      InputLeptonElectron = Zero
      InputLeptonMuon = Zero
      InputLeptonTau = Zero
      DeltaEnergy = Zero
      DeltaMomentum1 = Zero
      DeltaMomentum2 = Zero
      DeltaMomentum3 = Zero

      IF( MASK(iN_X) )THEN

        V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

#if   defined(THORNADO_OMP_OL)
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( Dnu_Input, DeltaDnu, &
        !$OMP          Inu_u_1_Input, Inu_u_2_Input, Inu_u_3_Input, &
        !$OMP          DeltaInu_u_1, DeltaInu_u_2, DeltaInu_u_3, &
        !$OMP          DeltaFnu_d_1, DeltaFnu_d_2, DeltaFnu_d_3, &
        !$OMP          Nnu_Input, Nnu_Limited, &
        !$OMP          Gnu_Input_d_1, Gnu_Input_d_2, Gnu_Input_d_3, &
        !$OMP          Gnu_Limited_d_1, Gnu_Limited_d_2, &
        !$OMP          Gnu_Limited_d_3 ) &
        !$OMP REDUCTION( + : DeltaLepton, DeltaLeptonElectron, &
        !$OMP                 DeltaLeptonMuon, DeltaLeptonTau, DeltaEnergy, &
        !$OMP                 InputLeptonElectron, InputLeptonMuon, &
        !$OMP                 InputLeptonTau, &
        !$OMP                 DeltaMomentum1, DeltaMomentum2, DeltaMomentum3 )
#elif defined(THORNADO_OACC)
        !$ACC LOOP VECTOR COLLAPSE(2) &
        !$ACC PRIVATE( Dnu_Input, DeltaDnu, &
        !$ACC          Inu_u_1_Input, Inu_u_2_Input, Inu_u_3_Input, &
        !$ACC          DeltaInu_u_1, DeltaInu_u_2, DeltaInu_u_3, &
        !$ACC          DeltaFnu_d_1, DeltaFnu_d_2, DeltaFnu_d_3, &
        !$ACC          Nnu_Input, Nnu_Limited, &
        !$ACC          Gnu_Input_d_1, Gnu_Input_d_2, Gnu_Input_d_3, &
        !$ACC          Gnu_Limited_d_1, Gnu_Limited_d_2, &
        !$ACC          Gnu_Limited_d_3 ) &
        !$ACC REDUCTION( + : DeltaLepton, DeltaLeptonElectron, &
        !$ACC                 DeltaLeptonMuon, DeltaLeptonTau, DeltaEnergy, &
        !$ACC                 InputLeptonElectron, InputLeptonMuon, &
        !$ACC                 InputLeptonTau, &
        !$ACC                 DeltaMomentum1, DeltaMomentum2, DeltaMomentum3 )
#endif
        DO iS   = 1, nSpecies
        DO iN_E = 1, nE_G

          Dnu_Input = Dnu(iN_E,iS,iN_X)
          Inu_u_1_Input = Inu_u_1(iN_E,iS,iN_X)
          Inu_u_2_Input = Inu_u_2(iN_E,iS,iN_X)
          Inu_u_3_Input = Inu_u_3(iN_E,iS,iN_X)

          Nnu_Input = Dnu_Input &
            + V_u_1(iN_X) * Gm_dd_11(iN_X) * Inu_u_1_Input &
            + V_u_2(iN_X) * Gm_dd_22(iN_X) * Inu_u_2_Input &
            + V_u_3(iN_X) * Gm_dd_33(iN_X) * Inu_u_3_Input
          DeltaNumber = ( Nnu_Input - Nnu_Conserved(iN_E,iS,iN_X) ) &
                          * W2_S(iN_E) * LeptonNumber(iS)
          SELECT CASE( ( iS + 1 ) / 2 )
            CASE( 1 )
              InputLeptonElectron = InputLeptonElectron + DeltaNumber
            CASE( 2 )
              InputLeptonMuon = InputLeptonMuon + DeltaNumber
            CASE( 3 )
              InputLeptonTau = InputLeptonTau + DeltaNumber
          END SELECT

          CALL ProjectPrimitiveRadiationMoment_OrderV &
                 ( Dnu    (iN_E,iS,iN_X), &
                   Inu_u_1(iN_E,iS,iN_X), &
                   Inu_u_2(iN_E,iS,iN_X), &
                   Inu_u_3(iN_E,iS,iN_X), &
                   Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                   DnuMax )

          DeltaDnu = Dnu(iN_E,iS,iN_X) - Dnu_Input
          DeltaInu_u_1 = Inu_u_1(iN_E,iS,iN_X) - Inu_u_1_Input
          DeltaInu_u_2 = Inu_u_2(iN_E,iS,iN_X) - Inu_u_2_Input
          DeltaInu_u_3 = Inu_u_3(iN_E,iS,iN_X) - Inu_u_3_Input

          IF(      DeltaDnu     /= Zero &
              .OR. DeltaInu_u_1 /= Zero &
              .OR. DeltaInu_u_2 /= Zero &
              .OR. DeltaInu_u_3 /= Zero )THEN

            CALL ComputeConserved_TwoMoment &
                   ( Dnu_Input, &
                     Inu_u_1_Input, Inu_u_2_Input, Inu_u_3_Input, &
                     Nnu_Input, Gnu_Input_d_1, Gnu_Input_d_2, Gnu_Input_d_3, &
                     V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

            CALL ComputeConserved_TwoMoment &
                   ( Dnu(iN_E,iS,iN_X), &
                     Inu_u_1(iN_E,iS,iN_X), &
                     Inu_u_2(iN_E,iS,iN_X), &
                     Inu_u_3(iN_E,iS,iN_X), &
                     Nnu_Limited, &
                     Gnu_Limited_d_1, Gnu_Limited_d_2, Gnu_Limited_d_3, &
                     V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

            DeltaFnu_d_1 &
              = V_d_1 * DeltaDnu &
                  + Gnu_Limited_d_1 - Gnu_Input_d_1
            DeltaFnu_d_2 &
              = V_d_2 * DeltaDnu &
                  + Gnu_Limited_d_2 - Gnu_Input_d_2
            DeltaFnu_d_3 &
              = V_d_3 * DeltaDnu &
                  + Gnu_Limited_d_3 - Gnu_Input_d_3

            DeltaNumber = ( Nnu_Limited - Nnu_Input ) &
                            * W2_S(iN_E) * LeptonNumber(iS)
            SELECT CASE( ( iS + 1 ) / 2 )
              CASE( 1 )
                DeltaLeptonElectron = DeltaLeptonElectron + DeltaNumber
              CASE( 2 )
                DeltaLeptonMuon = DeltaLeptonMuon + DeltaNumber
              CASE( 3 )
                DeltaLeptonTau = DeltaLeptonTau + DeltaNumber
            END SELECT

            DeltaEnergy = DeltaEnergy &
              + ( DeltaDnu &
                  + Two &
                      * (   V_d_1 * DeltaInu_u_1 &
                          + V_d_2 * DeltaInu_u_2 &
                          + V_d_3 * DeltaInu_u_3 ) ) * W3_S(iN_E)
            DeltaMomentum1 &
              = DeltaMomentum1 + DeltaFnu_d_1 * W3_S(iN_E)
            DeltaMomentum2 &
              = DeltaMomentum2 + DeltaFnu_d_2 * W3_S(iN_E)
            DeltaMomentum3 &
              = DeltaMomentum3 + DeltaFnu_d_3 * W3_S(iN_E)

          END IF

        END DO
        END DO

      END IF

      LimiterRate(:,iN_X) = Zero
      DeltaLepton = DeltaLeptonElectron
      LimiterRate(iQ_L,iN_X) = DeltaLepton / dt
      LimiterRate(iQ_E,iN_X) = DeltaEnergy / dt
      LimiterRate(iQ_M1,iN_X) = DeltaMomentum1 / dt
      LimiterRate(iQ_M2,iN_X) = DeltaMomentum2 / dt
      LimiterRate(iQ_M3,iN_X) = DeltaMomentum3 / dt
      LimiterRate(iQ_H,iN_X) &
        = LimiterRate(iQ_E,iN_X) &
            - V_u_1(iN_X) * LimiterRate(iQ_M1,iN_X) &
            - V_u_2(iN_X) * LimiterRate(iQ_M2,iN_X) &
            - V_u_3(iN_X) * LimiterRate(iQ_M3,iN_X)

      ConservationCorrection(:,iN_X) = Zero
      ConservationCorrection(iY ,iN_X) = -DeltaLepton
      ConservationCorrection(iEf,iN_X) = -DeltaEnergy
      ConservationCorrection(iV1,iN_X) = -DeltaMomentum1
      ConservationCorrection(iV2,iN_X) = -DeltaMomentum2
      ConservationCorrection(iV3,iN_X) = -DeltaMomentum3

      LeptonFlavorLimiterRate(1,iN_X) = DeltaLeptonElectron / dt
      LeptonFlavorLimiterRate(2,iN_X) = DeltaLeptonMuon / dt
      LeptonFlavorLimiterRate(3,iN_X) = DeltaLeptonTau / dt
      LeptonFlavorInputConversionRate(1,iN_X) = InputLeptonElectron / dt
      LeptonFlavorInputConversionRate(2,iN_X) = InputLeptonMuon / dt
      LeptonFlavorInputConversionRate(3,iN_X) = InputLeptonTau / dt

    END DO

  END SUBROUTINE ApplyNeutrinoDistributionLimit_OrderV


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
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Nnu_Clipped, &
      Gnu_d_1_Clipped, Gnu_d_2_Clipped, Gnu_d_3_Clipped, &
      RadiationLimiterCorrection )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Y, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: &
      Nnu_Clipped, Gnu_d_1_Clipped, Gnu_d_2_Clipped, Gnu_d_3_Clipped
    REAL(DP), DIMENSION(:,:)  , INTENT(in)  :: RadiationLimiterCorrection

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

        Nnu_Clipped(iN_E,iS,iN_X) = C_Dnu(iN_E,iS,iN_X)
        Gnu_d_1_Clipped(iN_E,iS,iN_X) = C_Inu_d_1(iN_E,iS,iN_X)
        Gnu_d_2_Clipped(iN_E,iS,iN_X) = C_Inu_d_2(iN_E,iS,iN_X)
        Gnu_d_3_Clipped(iN_E,iS,iN_X) = C_Inu_d_3(iN_E,iS,iN_X)

        IF ( iS <= iNuE_Bar ) THEN
        SUM_Y  = SUM_Y  + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
        END IF

        SUM_Ef = SUM_Ef + Enu     * W3_S(iN_E)
        SUM_V1 = SUM_V1 + Fnu_d_1 * W3_S(iN_E)
        SUM_V2 = SUM_V2 + Fnu_d_2 * W3_S(iN_E)
        SUM_V3 = SUM_V3 + Fnu_d_3 * W3_S(iN_E)

      END DO
      END DO

      ! --- Store the Old Radiation Contribution to Each Conserved Total ---
      ! --- Include the radiation removed by the occupation-number bound. ---

      C_Y(iN_X) = wMatrRHS(iY) &
        * ( SUM_Y + RadiationLimiterCorrection(iY,iN_X) ) * S_Y(iN_X)
      C_Ef(iN_X) = wMatrRHS(iEf) &
        * ( SUM_Ef + RadiationLimiterCorrection(iEf,iN_X) ) * S_Ef(iN_X)
      C_V_d_1(iN_X) = wMatrRHS(iV1) &
        * ( SUM_V1 + RadiationLimiterCorrection(iV1,iN_X) ) &
        * S_V_d_1(iN_X)
      C_V_d_2(iN_X) = wMatrRHS(iV2) &
        * ( SUM_V2 + RadiationLimiterCorrection(iV2,iN_X) ) &
        * S_V_d_2(iN_X)
      C_V_d_3(iN_X) = wMatrRHS(iV3) &
        * ( SUM_V3 + RadiationLimiterCorrection(iV3,iN_X) ) &
        * S_V_d_3(iN_X)

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


      IF ( wMatrRHS(iD) .EQ. 1.0_DP) THEN
        C_D    (iN_X) &
          = W * U_D(iN_X)
      ELSE  
        C_D    (iN_X) &
          = U_D(iN_X)
      END IF

      IF ( wMatrRHS(iY) .EQ. 1.0_DP) THEN
        C_Y    (iN_X) &
          = U_Y(iN_X) + SUM_Y * S_Y(iN_X)
      ELSE
        C_Y    (iN_X) &
          = U_Y(iN_X)
      END IF
       
      IF ( wMatrRHS(iEF) .EQ. 1.0_DP) THEN
        C_Ef   (iN_X) &
          = ( Tau_f + SUM_Ef ) * S_Ef  (iN_X) / cD_old(iN_X)
      ELSE
        C_Ef   (iN_X) &
          = U_Ef(iN_X)
      END IF

      IF ( wMatrRHS(iV1) .EQ. 1.0_DP) THEN
        C_V_d_1(iN_X) &
          = ( SJ(1) + SUM_V1 ) * S_V_d_1(iN_X) / cD_old(iN_X)
      ELSE
        C_V_d_1(iN_X) &
          = U_V_d_1(iN_X)
      END IF


      IF ( wMatrRHS(iV2) .EQ. 1.0_DP) THEN
        C_V_d_2(iN_X) &
          = ( SJ(2) + SUM_V2 ) * S_V_d_2(iN_X) / cD_old(iN_X)
      ELSE
        C_V_d_2(iN_X) &
          = U_V_d_2(iN_X)
      END IF


      IF ( wMatrRHS(iV3) .EQ. 1.0_DP) THEN
        C_V_d_3(iN_X) &
          = ( SJ(3) + SUM_V3 ) * S_V_d_3(iN_X) / cD_old(iN_X)
      ELSE
        C_V_d_3(iN_X) &
          = U_V_d_3(iN_X)
      END IF





    END DO

  END SUBROUTINE InitializeRHS_Relativistic


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE ComputeMatterRHS_Collisions_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, MASK, Q_Matter, &
      NNS_Ledger_Option, QH_ByProcessSpecies_Option, EmAb_Ledger_Option, &
      EmAb_Spectrum_Option, RadiationCollisionRHS_Option )

    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)  :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    LOGICAL,  DIMENSION(:),     INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: Q_Matter
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: NNS_Ledger_Option
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: QH_ByProcessSpecies_Option
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: EmAb_Ledger_Option
    REAL(DP), DIMENSION(:,:,:,:), INTENT(out), OPTIONAL :: EmAb_Spectrum_Option
    REAL(DP), DIMENSION(:,:,:,:), INTENT(out), OPTIONAL :: &
      RadiationCollisionRHS_Option

    REAL(DP), ALLOCATABLE :: Q_SP(:,:,:,:)
    REAL(DP), ALLOCATABLE :: NNS_Ledger(:,:,:)
    REAL(DP), ALLOCATABLE :: EmAb_Ledger(:,:,:)
    REAL(DP), ALLOCATABLE :: EmAb_Spectrum(:,:,:,:)
    REAL(DP), ALLOCATABLE :: RadiationCollisionRHS(:,:,:,:)
    REAL(DP), ALLOCATABLE :: QH_ByProcessSpecies(:,:,:)

    INTEGER  :: ierr, myid
    LOGICAL  :: ComputeEmAbLedger, InvalidShape
    INTEGER  :: iG, iL, iN_E, iN_X, iP, iQ, iS
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: RN(nC)
    REAL(DP) :: RG1(nC)
    REAL(DP) :: RG2(nC)
    REAL(DP) :: RG3(nC)
    REAL(DP) :: QL, QE, QM1, QM2, QM3
    REAL(DP) :: RN_Total, RG1_Total, RG2_Total, RG3_Total
    REAL(DP) :: V2

    ComputeEmAbLedger = PRESENT( EmAb_Ledger_Option )

    InvalidShape = SIZE(MASK) /= nX_G .OR. &
                   SIZE(Q_Matter,1) /= nQ .OR. &
                   SIZE(Q_Matter,2) /= 1 + nSpecies + nC .OR. &
                   SIZE(Q_Matter,3) /= nX_G

    IF( PRESENT( QH_ByProcessSpecies_Option ) )THEN
      InvalidShape = InvalidShape .OR. &
                     SIZE(QH_ByProcessSpecies_Option,1) /= nC .OR. &
                     SIZE(QH_ByProcessSpecies_Option,2) /= nSpecies .OR. &
                     SIZE(QH_ByProcessSpecies_Option,3) /= nX_G
    END IF

    IF( PRESENT( EmAb_Ledger_Option ) )THEN
      InvalidShape = InvalidShape .OR. &
                     SIZE(EmAb_Ledger_Option,1) /= nEmAbLedger .OR. &
                     SIZE(EmAb_Ledger_Option,2) /= nSpecies .OR. &
                     SIZE(EmAb_Ledger_Option,3) /= nX_G
    END IF

    IF( PRESENT( EmAb_Spectrum_Option ) )THEN
      InvalidShape = InvalidShape .OR. &
                     SIZE(EmAb_Spectrum_Option,1) /= nEmAbSpectrum .OR. &
                     SIZE(EmAb_Spectrum_Option,2) /= nE_G .OR. &
                     SIZE(EmAb_Spectrum_Option,3) /= nSpecies .OR. &
                     SIZE(EmAb_Spectrum_Option,4) /= nX_G
    END IF

    IF( PRESENT( RadiationCollisionRHS_Option ) )THEN
      InvalidShape = InvalidShape .OR. &
                     SIZE(RadiationCollisionRHS_Option,1) /= nCR .OR. &
                     SIZE(RadiationCollisionRHS_Option,2) /= nE_G .OR. &
                     SIZE(RadiationCollisionRHS_Option,3) /= nSpecies .OR. &
                     SIZE(RadiationCollisionRHS_Option,4) /= nX_G
    END IF

    IF( InvalidShape )THEN

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

      IF ( myid == 0 ) THEN
        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: Invalid output shape in ComputeMatterRHS_Collisions_OrderV.'
        FLUSH(ERROR_UNIT)
      END IF

      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

    END IF

    ALLOCATE( Q_SP(5,nC,nSpecies,nX_G) )
    ALLOCATE( NNS_Ledger(nNNSLedger,nSpecies,nX_G) )
    ALLOCATE( EmAb_Ledger(nEmAbLedger,nSpecies,nX_G) )
    ALLOCATE( RadiationCollisionRHS(nCR,nE_G,nSpecies,nX_G) )
    RadiationCollisionRHS = Zero

    IF( PRESENT( QH_ByProcessSpecies_Option ) )THEN
      ALLOCATE( QH_ByProcessSpecies(nC,nSpecies,nX_G) )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: QH_ByProcessSpecies )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( QH_ByProcessSpecies )
#endif
    END IF


    IF( PRESENT( EmAb_Spectrum_Option ) )THEN
      ALLOCATE( EmAb_Spectrum(nEmAbSpectrum,nE_G,nSpecies,nX_G) )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: EmAb_Spectrum )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( EmAb_Spectrum )
#endif
    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( to: RadiationCollisionRHS )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN( RadiationCollisionRHS )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: Q_SP, NNS_Ledger, EmAb_Ledger )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( Q_SP, NNS_Ledger, EmAb_Ledger )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( Q_SP, Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iP   = 1, nC
    DO iQ   = 1, 5
      Q_SP(iQ,iP,iS,iN_X) = Zero
    END DO
    END DO
    END DO
    END DO

    IF( ComputeEmAbLedger )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( EmAb_Ledger )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iL   = 1, nEmAbLedger
      EmAb_Ledger(iL,iS,iN_X) = Zero
    END DO
    END DO
    END DO

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( NNS_Ledger )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iL   = 1, nNNSLedger
      NNS_Ledger(iL,iS,iN_X) = Zero
    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iG   = 1, 1 + nSpecies + nC
    DO iQ   = 1, nQ
      Q_Matter(iQ,iG,iN_X) = Zero
    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP FIRSTPRIVATE( ComputeEmAbLedger ) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$OMP          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$OMP          RG3_Total, QL, QE, QM1, QM2, QM3, iN_E, iP, iG )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC FIRSTPRIVATE( ComputeEmAbLedger ) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$ACC          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$ACC          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$ACC          RG3_Total, QL, QE, QM1, QM2, QM3, iN_E, iP, iG ) &
    !$ACC PRESENT( Q_SP, NNS_Ledger, EmAb_Ledger, RadiationCollisionRHS )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP FIRSTPRIVATE( ComputeEmAbLedger ) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$OMP          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$OMP          RG3_Total, QL, QE, QM1, QM2, QM3, iN_E, iP, iG )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies

      IF ( MASK(iN_X) ) THEN

        V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

        DO iN_E = 1, nE_G

          I_d_1 = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)
          I_d_2 = Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)
          I_d_3 = Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)

          CALL ComputeEddingtonTensorComponents_dd &
                 ( Dnu(iN_E,iS,iN_X), &
                   Inu_u_1(iN_E,iS,iN_X), Inu_u_2(iN_E,iS,iN_X), &
                   Inu_u_3(iN_E,iS,iN_X), &
                   Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                   k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33 )

          RN  = Zero
          RG1 = Zero
          RG2 = Zero
          RG3 = Zero

          RN(iC_EmAb) = Chi_EmAb(iN_E,iS,iN_X) &
                         * ( Dnu_0(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) )
          RG1(iC_EmAb) = -Chi_EmAb(iN_E,iS,iN_X) * I_d_1
          RG2(iC_EmAb) = -Chi_EmAb(iN_E,iS,iN_X) * I_d_2
          RG3(iC_EmAb) = -Chi_EmAb(iN_E,iS,iN_X) * I_d_3

          IF( ComputeEmAbLedger )THEN
            EmAb_Ledger(iEmAbLedger_NIn,iS,iN_X) = &
              EmAb_Ledger(iEmAbLedger_NIn,iS,iN_X) &
                + Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) * W2_S(iN_E)
            EmAb_Ledger(iEmAbLedger_NOut,iS,iN_X) = &
              EmAb_Ledger(iEmAbLedger_NOut,iS,iN_X) &
                + Chi_EmAb(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W2_S(iN_E)
            EmAb_Ledger(iEmAbLedger_NNet,iS,iN_X) = &
              EmAb_Ledger(iEmAbLedger_NNet,iS,iN_X) &
                + RN(iC_EmAb) * W2_S(iN_E)
            EmAb_Ledger(iEmAbLedger_EIn,iS,iN_X) = &
              EmAb_Ledger(iEmAbLedger_EIn,iS,iN_X) &
                + Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) * W3_S(iN_E)
            EmAb_Ledger(iEmAbLedger_EOut,iS,iN_X) = &
              EmAb_Ledger(iEmAbLedger_EOut,iS,iN_X) &
                + Chi_EmAb(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W3_S(iN_E)
          END IF

          iG = nChirals - ( LeptonNumber(iS) + 1 ) / nChirals
          RG1(iC_Iso) = -Sigma_Iso(iN_E,iG,iN_X) * I_d_1
          RG2(iC_Iso) = -Sigma_Iso(iN_E,iG,iN_X) * I_d_2
          RG3(iC_Iso) = -Sigma_Iso(iN_E,iG,iN_X) * I_d_3

          IF( Include_NNS )THEN
            RN(iC_NNS) = R_NNS(iN_E,iS,iN_X)
            RG1(iC_NNS) = -Chi_NNS(iN_E,iS,iN_X) * I_d_1
            RG2(iC_NNS) = -Chi_NNS(iN_E,iS,iN_X) * I_d_2
            RG3(iC_NNS) = -Chi_NNS(iN_E,iS,iN_X) * I_d_3

            IF( Include_LinCorr )THEN
              L_u_1 = L_NNS__In__u_1(iN_E,iS,iN_X) &
                      - L_NNS__Out_u_1(iN_E,iS,iN_X)
              L_u_2 = L_NNS__In__u_2(iN_E,iS,iN_X) &
                      - L_NNS__Out_u_2(iN_E,iS,iN_X)
              L_u_3 = L_NNS__In__u_3(iN_E,iS,iN_X) &
                      - L_NNS__Out_u_3(iN_E,iS,iN_X)
              L_d_1 = Gm_dd_11(iN_X) * L_NNS__In__u_1(iN_E,iS,iN_X)
              L_d_2 = Gm_dd_22(iN_X) * L_NNS__In__u_2(iN_E,iS,iN_X)
              L_d_3 = Gm_dd_33(iN_X) * L_NNS__In__u_3(iN_E,iS,iN_X)
              RN(iC_NNS) = RN(iC_NNS) &
                           - ( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
              RG1(iC_NNS) = RG1(iC_NNS) + Third * L_d_1 &
                            - ( L_u_1 * k_dd_11 + L_u_2 * k_dd_12 &
                              + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
              RG2(iC_NNS) = RG2(iC_NNS) + Third * L_d_2 &
                            - ( L_u_1 * k_dd_12 + L_u_2 * k_dd_22 &
                              + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
              RG3(iC_NNS) = RG3(iC_NNS) + Third * L_d_3 &
                            - ( L_u_1 * k_dd_13 + L_u_2 * k_dd_23 &
                              + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)
            END IF
          END IF

          IF( Include_NNS )THEN
            NNS_Ledger(iNNSLedger_NIn ,iS,iN_X) = &
              NNS_Ledger(iNNSLedger_NIn ,iS,iN_X) &
                + Eta_NNS(iN_E,iS,iN_X) * W2_S(iN_E)
            NNS_Ledger(iNNSLedger_NOut,iS,iN_X) = &
              NNS_Ledger(iNNSLedger_NOut,iS,iN_X) &
                + Chi_NNS(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W2_S(iN_E)
            NNS_Ledger(iNNSLedger_NNet,iS,iN_X) = &
              NNS_Ledger(iNNSLedger_NNet,iS,iN_X) &
                + RN(iC_NNS) * W2_S(iN_E)
            NNS_Ledger(iNNSLedger_EIn ,iS,iN_X) = &
              NNS_Ledger(iNNSLedger_EIn ,iS,iN_X) &
                + Eta_NNS(iN_E,iS,iN_X) * W3_S(iN_E)
            NNS_Ledger(iNNSLedger_EOut,iS,iN_X) = &
              NNS_Ledger(iNNSLedger_EOut,iS,iN_X) &
                + Chi_NNS(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W3_S(iN_E)
          END IF

          IF( Include_NES )THEN
            RN(iC_NES) = R_NES(iN_E,iS,iN_X)
            RG1(iC_NES) = -Chi_NES(iN_E,iS,iN_X) * I_d_1
            RG2(iC_NES) = -Chi_NES(iN_E,iS,iN_X) * I_d_2
            RG3(iC_NES) = -Chi_NES(iN_E,iS,iN_X) * I_d_3

            IF( Include_LinCorr )THEN
              L_u_1 = L_NES__In__u_1(iN_E,iS,iN_X) &
                      - L_NES__Out_u_1(iN_E,iS,iN_X)
              L_u_2 = L_NES__In__u_2(iN_E,iS,iN_X) &
                      - L_NES__Out_u_2(iN_E,iS,iN_X)
              L_u_3 = L_NES__In__u_3(iN_E,iS,iN_X) &
                      - L_NES__Out_u_3(iN_E,iS,iN_X)
              L_d_1 = Gm_dd_11(iN_X) * L_NES__In__u_1(iN_E,iS,iN_X)
              L_d_2 = Gm_dd_22(iN_X) * L_NES__In__u_2(iN_E,iS,iN_X)
              L_d_3 = Gm_dd_33(iN_X) * L_NES__In__u_3(iN_E,iS,iN_X)
              RN(iC_NES) = RN(iC_NES) &
                           - ( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
              RG1(iC_NES) = RG1(iC_NES) + Third * L_d_1 &
                            - ( L_u_1 * k_dd_11 + L_u_2 * k_dd_12 &
                              + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
              RG2(iC_NES) = RG2(iC_NES) + Third * L_d_2 &
                            - ( L_u_1 * k_dd_12 + L_u_2 * k_dd_22 &
                              + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
              RG3(iC_NES) = RG3(iC_NES) + Third * L_d_3 &
                            - ( L_u_1 * k_dd_13 + L_u_2 * k_dd_23 &
                              + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)
            END IF
          END IF

          IF( Include_Pair )THEN
            RN(iC_Pair) = R_Pair(iN_E,iS,iN_X)
            RG1(iC_Pair) = -Chi_Pair(iN_E,iS,iN_X) * I_d_1
            RG2(iC_Pair) = -Chi_Pair(iN_E,iS,iN_X) * I_d_2
            RG3(iC_Pair) = -Chi_Pair(iN_E,iS,iN_X) * I_d_3

            IF( Include_LinCorr )THEN
              L_u_1 = L_Pair_Ann_u_1(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_1(iN_E,iS,iN_X)
              L_u_2 = L_Pair_Ann_u_2(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_2(iN_E,iS,iN_X)
              L_u_3 = L_Pair_Ann_u_3(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_3(iN_E,iS,iN_X)
              L_d_1 = -Gm_dd_11(iN_X) * L_Pair_Pro_u_1(iN_E,iS,iN_X)
              L_d_2 = -Gm_dd_22(iN_X) * L_Pair_Pro_u_2(iN_E,iS,iN_X)
              L_d_3 = -Gm_dd_33(iN_X) * L_Pair_Pro_u_3(iN_E,iS,iN_X)
              RN(iC_Pair) = RN(iC_Pair) &
                            - ( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
              RG1(iC_Pair) = RG1(iC_Pair) + Third * L_d_1 &
                             - ( L_u_1 * k_dd_11 + L_u_2 * k_dd_12 &
                               + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
              RG2(iC_Pair) = RG2(iC_Pair) + Third * L_d_2 &
                             - ( L_u_1 * k_dd_12 + L_u_2 * k_dd_22 &
                               + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
              RG3(iC_Pair) = RG3(iC_Pair) + Third * L_d_3 &
                             - ( L_u_1 * k_dd_13 + L_u_2 * k_dd_23 &
                               + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)
            END IF
          END IF

          IF( Include_NuPair )THEN
            RN(iC_NuPair) = R_NuPair(iN_E,iS,iN_X)
            RG1(iC_NuPair) = -Chi_NuPair(iN_E,iS,iN_X) * I_d_1
            RG2(iC_NuPair) = -Chi_NuPair(iN_E,iS,iN_X) * I_d_2
            RG3(iC_NuPair) = -Chi_NuPair(iN_E,iS,iN_X) * I_d_3

            IF( Include_LinCorr )THEN
              L_u_1 = L_NuPair_Ann_u_1(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_1(iN_E,iS,iN_X)
              L_u_2 = L_NuPair_Ann_u_2(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_2(iN_E,iS,iN_X)
              L_u_3 = L_NuPair_Ann_u_3(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_3(iN_E,iS,iN_X)
              L_d_1 = -Gm_dd_11(iN_X) * L_NuPair_Pro_u_1(iN_E,iS,iN_X)
              L_d_2 = -Gm_dd_22(iN_X) * L_NuPair_Pro_u_2(iN_E,iS,iN_X)
              L_d_3 = -Gm_dd_33(iN_X) * L_NuPair_Pro_u_3(iN_E,iS,iN_X)
              RN(iC_NuPair) = RN(iC_NuPair) &
                              - ( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
              RG1(iC_NuPair) = RG1(iC_NuPair) + Third * L_d_1 &
                               - ( L_u_1 * k_dd_11 + L_u_2 * k_dd_12 &
                                 + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
              RG2(iC_NuPair) = RG2(iC_NuPair) + Third * L_d_2 &
                               - ( L_u_1 * k_dd_12 + L_u_2 * k_dd_22 &
                                 + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
              RG3(iC_NuPair) = RG3(iC_NuPair) + Third * L_d_3 &
                               - ( L_u_1 * k_dd_13 + L_u_2 * k_dd_23 &
                                 + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)
            END IF
          END IF

          IF( Include_Brem )THEN
            RN(iC_Brem) = R_Brem(iN_E,iS,iN_X)
            RG1(iC_Brem) = -Chi_Brem(iN_E,iS,iN_X) * I_d_1
            RG2(iC_Brem) = -Chi_Brem(iN_E,iS,iN_X) * I_d_2
            RG3(iC_Brem) = -Chi_Brem(iN_E,iS,iN_X) * I_d_3

            IF( Include_LinCorr )THEN
              L_u_1 = L_Brem_Ann_u_1(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_1(iN_E,iS,iN_X)
              L_u_2 = L_Brem_Ann_u_2(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_2(iN_E,iS,iN_X)
              L_u_3 = L_Brem_Ann_u_3(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_3(iN_E,iS,iN_X)
              L_d_1 = -Gm_dd_11(iN_X) * L_Brem_Pro_u_1(iN_E,iS,iN_X)
              L_d_2 = -Gm_dd_22(iN_X) * L_Brem_Pro_u_2(iN_E,iS,iN_X)
              L_d_3 = -Gm_dd_33(iN_X) * L_Brem_Pro_u_3(iN_E,iS,iN_X)
              RN(iC_Brem) = RN(iC_Brem) &
                            - ( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
              RG1(iC_Brem) = RG1(iC_Brem) + Third * L_d_1 &
                             - ( L_u_1 * k_dd_11 + L_u_2 * k_dd_12 &
                               + L_u_3 * k_dd_13 ) * Dnu(iN_E,iS,iN_X)
              RG2(iC_Brem) = RG2(iC_Brem) + Third * L_d_2 &
                             - ( L_u_1 * k_dd_12 + L_u_2 * k_dd_22 &
                               + L_u_3 * k_dd_23 ) * Dnu(iN_E,iS,iN_X)
              RG3(iC_Brem) = RG3(iC_Brem) + Third * L_d_3 &
                             - ( L_u_1 * k_dd_13 + L_u_2 * k_dd_23 &
                               + L_u_3 * k_dd_33 ) * Dnu(iN_E,iS,iN_X)
            END IF
          END IF

          RN_Total  = SUM( RN  )
          RG1_Total = SUM( RG1 )
          RG2_Total = SUM( RG2 )
          RG3_Total = SUM( RG3 )

          RadiationCollisionRHS(iCR_N,iN_E,iS,iN_X) &
            = RN_Total + V_u_1(iN_X) * RG1_Total &
                       + V_u_2(iN_X) * RG2_Total &
                       + V_u_3(iN_X) * RG3_Total
          RadiationCollisionRHS(iCR_G1,iN_E,iS,iN_X) &
            = RG1_Total + V_d_1 * RN_Total
          RadiationCollisionRHS(iCR_G2,iN_E,iS,iN_X) &
            = RG2_Total + V_d_2 * RN_Total
          RadiationCollisionRHS(iCR_G3,iN_E,iS,iN_X) &
            = RG3_Total + V_d_3 * RN_Total

          DO iP = 1, nC

            QL = Zero
            IF ( iS <= iNuE_Bar ) &
              QL = -LeptonNumber(iS) * RN(iP) * W2_S(iN_E)

            QE  = -( RN(iP) + V_u_1(iN_X) * RG1(iP) &
                              + V_u_2(iN_X) * RG2(iP) &
                              + V_u_3(iN_X) * RG3(iP) ) * W3_S(iN_E)
            QM1 = -( RG1(iP) + V_d_1 * RN(iP) ) * W3_S(iN_E)
            QM2 = -( RG2(iP) + V_d_2 * RN(iP) ) * W3_S(iN_E)
            QM3 = -( RG3(iP) + V_d_3 * RN(iP) ) * W3_S(iN_E)

            IF( iP == iC_NNS ) THEN
              NNS_Ledger(iNNSLedger_QH,iS,iN_X) = &
                NNS_Ledger(iNNSLedger_QH,iS,iN_X) &
                  + QE - V_u_1(iN_X) * QM1 &
                       - V_u_2(iN_X) * QM2 &
                       - V_u_3(iN_X) * QM3
            END IF

            IF( ComputeEmAbLedger .AND. iP == iC_EmAb ) THEN
              EmAb_Ledger(iEmAbLedger_QH,iS,iN_X) = &
                EmAb_Ledger(iEmAbLedger_QH,iS,iN_X) &
                  + QE - V_u_1(iN_X) * QM1 &
                       - V_u_2(iN_X) * QM2 &
                       - V_u_3(iN_X) * QM3
            END IF

            Q_SP(1,iP,iS,iN_X) = Q_SP(1,iP,iS,iN_X) + QL
            Q_SP(2,iP,iS,iN_X) = Q_SP(2,iP,iS,iN_X) + QE
            Q_SP(3,iP,iS,iN_X) = Q_SP(3,iP,iS,iN_X) + QM1
            Q_SP(4,iP,iS,iN_X) = Q_SP(4,iP,iS,iN_X) + QM2
            Q_SP(5,iP,iS,iN_X) = Q_SP(5,iP,iS,iN_X) + QM3

          END DO

        END DO

      END IF

    END DO
    END DO

    IF( PRESENT( EmAb_Spectrum_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( V2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( V2 ) &
    !$ACC PRESENT( EmAb_Spectrum )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE( V2 )
#endif
      DO iN_X = 1, nX_G
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        IF( MASK(iN_X) )THEN
          V2 = Gm_dd_11(iN_X) * V_u_1(iN_X)**2 &
             + Gm_dd_22(iN_X) * V_u_2(iN_X)**2 &
             + Gm_dd_33(iN_X) * V_u_3(iN_X)**2

          EmAb_Spectrum(iEmAbSpectrum_Dnu,iN_E,iS,iN_X) &
            = Dnu(iN_E,iS,iN_X)
          EmAb_Spectrum(iEmAbSpectrum_DnuEq,iN_E,iS,iN_X) &
            = Dnu_0(iN_E,iS,iN_X)
          EmAb_Spectrum(iEmAbSpectrum_Chi,iN_E,iS,iN_X) &
            = Chi_EmAb(iN_E,iS,iN_X)
          EmAb_Spectrum(iEmAbSpectrum_EIn,iN_E,iS,iN_X) &
            = Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) &
                * W3_S(iN_E)
          EmAb_Spectrum(iEmAbSpectrum_EOut,iN_E,iS,iN_X) &
            = Chi_EmAb(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) &
                * W3_S(iN_E)
          EmAb_Spectrum(iEmAbSpectrum_QH,iN_E,iS,iN_X) &
            = ( One - V2 ) * Chi_EmAb(iN_E,iS,iN_X) &
                * ( Dnu(iN_E,iS,iN_X) - Dnu_0(iN_E,iS,iN_X) ) &
                * W3_S(iN_E)
        ELSE
          EmAb_Spectrum(:,iN_E,iS,iN_X) = Zero
        END IF

      END DO
      END DO
      END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( from: EmAb_Spectrum )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA COPYOUT( EmAb_Spectrum )
#endif
      EmAb_Spectrum_Option = EmAb_Spectrum
      DEALLOCATE( EmAb_Spectrum )

    END IF

    IF( PRESENT( RadiationCollisionRHS_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( from: RadiationCollisionRHS )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA COPYOUT( RadiationCollisionRHS )
#endif
      RadiationCollisionRHS_Option = RadiationCollisionRHS
    ELSE
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: RadiationCollisionRHS )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( RadiationCollisionRHS )
#endif
    END IF
    DEALLOCATE( RadiationCollisionRHS )

    ! Species sums.
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) PRIVATE( iP )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE( iP ) &
    !$ACC PRESENT( Q_SP, Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iP )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
      DO iP = 1, nC
        Q_Matter(iQ_L ,1+iS,iN_X) = Q_Matter(iQ_L ,1+iS,iN_X) &
                                             + Q_SP(1,iP,iS,iN_X)
        Q_Matter(iQ_E ,1+iS,iN_X) = Q_Matter(iQ_E ,1+iS,iN_X) &
                                             + Q_SP(2,iP,iS,iN_X)
        Q_Matter(iQ_M1,1+iS,iN_X) = Q_Matter(iQ_M1,1+iS,iN_X) &
                                             + Q_SP(3,iP,iS,iN_X)
        Q_Matter(iQ_M2,1+iS,iN_X) = Q_Matter(iQ_M2,1+iS,iN_X) &
                                             + Q_SP(4,iP,iS,iN_X)
        Q_Matter(iQ_M3,1+iS,iN_X) = Q_Matter(iQ_M3,1+iS,iN_X) &
                                             + Q_SP(5,iP,iS,iN_X)
      END DO
    END DO
    END DO

    ! Process sums.
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) PRIVATE( iS )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE( iS ) &
    !$ACC PRESENT( Q_SP, Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iS )
#endif
    DO iN_X = 1, nX_G
    DO iP   = 1, nC
      DO iS = 1, nSpecies
        Q_Matter(iQ_L ,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_L ,1+nSpecies+iP,iN_X) + Q_SP(1,iP,iS,iN_X)
        Q_Matter(iQ_E ,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_E ,1+nSpecies+iP,iN_X) + Q_SP(2,iP,iS,iN_X)
        Q_Matter(iQ_M1,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_M1,1+nSpecies+iP,iN_X) + Q_SP(3,iP,iS,iN_X)
        Q_Matter(iQ_M2,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_M2,1+nSpecies+iP,iN_X) + Q_SP(4,iP,iS,iN_X)
        Q_Matter(iQ_M3,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_M3,1+nSpecies+iP,iN_X) + Q_SP(5,iP,iS,iN_X)
      END DO
    END DO
    END DO

    ! Totals and internal-energy source.
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) PRIVATE( iS )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE( iS ) &
    !$ACC PRESENT( Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iS )
#endif
    DO iN_X = 1, nX_G
    DO iG   = 1, 1 + nSpecies + nC

      IF ( iG == 1 ) THEN
        DO iS = 1, nSpecies
          Q_Matter(iQ_L ,iG,iN_X) = Q_Matter(iQ_L ,iG,iN_X) &
                                            + Q_Matter(iQ_L ,1+iS,iN_X)
          Q_Matter(iQ_E ,iG,iN_X) = Q_Matter(iQ_E ,iG,iN_X) &
                                            + Q_Matter(iQ_E ,1+iS,iN_X)
          Q_Matter(iQ_M1,iG,iN_X) = Q_Matter(iQ_M1,iG,iN_X) &
                                            + Q_Matter(iQ_M1,1+iS,iN_X)
          Q_Matter(iQ_M2,iG,iN_X) = Q_Matter(iQ_M2,iG,iN_X) &
                                            + Q_Matter(iQ_M2,1+iS,iN_X)
          Q_Matter(iQ_M3,iG,iN_X) = Q_Matter(iQ_M3,iG,iN_X) &
                                            + Q_Matter(iQ_M3,1+iS,iN_X)
        END DO
      END IF

      Q_Matter(iQ_H,iG,iN_X) = Q_Matter(iQ_E,iG,iN_X) &
        - V_u_1(iN_X) * Q_Matter(iQ_M1,iG,iN_X) &
        - V_u_2(iN_X) * Q_Matter(iQ_M2,iG,iN_X) &
        - V_u_3(iN_X) * Q_Matter(iQ_M3,iG,iN_X)

    END DO
    END DO

    IF( PRESENT( QH_ByProcessSpecies_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( Q_SP, QH_ByProcessSpecies )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
      DO iS   = 1, nSpecies
      DO iP   = 1, nC
        QH_ByProcessSpecies(iP,iS,iN_X) = Q_SP(2,iP,iS,iN_X) &
          - V_u_1(iN_X) * Q_SP(3,iP,iS,iN_X) &
          - V_u_2(iN_X) * Q_SP(4,iP,iS,iN_X) &
          - V_u_3(iN_X) * Q_SP(5,iP,iS,iN_X)
      END DO
      END DO
      END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( from: QH_ByProcessSpecies )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA COPYOUT( QH_ByProcessSpecies )
#endif
      QH_ByProcessSpecies_Option = QH_ByProcessSpecies

      DEALLOCATE( QH_ByProcessSpecies )
    END IF

    IF( PRESENT( NNS_Ledger_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP( from: NNS_Ledger )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA COPYOUT( NNS_Ledger )
#endif
      NNS_Ledger_Option = NNS_Ledger
    ELSE
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP( release: NNS_Ledger )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA DELETE( NNS_Ledger )
#endif
    END IF

    IF( PRESENT( EmAb_Ledger_Option ) )THEN
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP( from: EmAb_Ledger )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA COPYOUT( EmAb_Ledger )
#endif
      EmAb_Ledger_Option = EmAb_Ledger
    ELSE
#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP( release: EmAb_Ledger )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA DELETE( EmAb_Ledger )
#endif
    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: Q_SP )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( Q_SP )
#endif

    DEALLOCATE( Q_SP, NNS_Ledger, EmAb_Ledger )

  END SUBROUTINE ComputeMatterRHS_Collisions_OrderV
#endif


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

        ! Form the conservation residual from matter and radiation
        ! increments.  This is algebraically equivalent to subtracting the
        ! new radiation state from the absolute conserved total, but it does
        ! not first combine the small radiation contribution with a matter
        ! state of order unity.

        Fm(iY,iN_X) = C_Y(iN_X) &
          - wMatrRHS(iY) * SUM_Y * S_Y(iN_X) &
          - ( U_Y(iN_X) - One )

        Fm(iEf,iN_X) = C_Ef(iN_X) &
          - wMatrRHS(iEf) * SUM_Ef * S_Ef(iN_X) &
          - ( U_Ef(iN_X) - One )

        Fm(iV1,iN_X) = C_V_d_1(iN_X) &
          - wMatrRHS(iV1) * SUM_V1 * S_V_d_1(iN_X) &
          - ( U_V_d_1(iN_X) &
              - Gm_dd_11(iN_X) * V_u_1_old(iN_X) / SpeedOfLight )

        Fm(iV2,iN_X) = C_V_d_2(iN_X) &
          - wMatrRHS(iV2) * SUM_V2 * S_V_d_2(iN_X) &
          - ( U_V_d_2(iN_X) &
              - Gm_dd_22(iN_X) * V_u_2_old(iN_X) / SpeedOfLight )

        Fm(iV3,iN_X) = C_V_d_3(iN_X) &
          - wMatrRHS(iV3) * SUM_V3 * S_V_d_3(iN_X) &
          - ( U_V_d_3(iN_X) &
              - Gm_dd_33(iN_X) * V_u_3_old(iN_X) / SpeedOfLight )

        G_Y    (iN_X) = U_Y    (iN_X) + Fm(iY ,iN_X)
        G_Ef   (iN_X) = U_Ef   (iN_X) + Fm(iEf,iN_X)
        G_V_d_1(iN_X) = U_V_d_1(iN_X) + Fm(iV1,iN_X)
        G_V_d_2(iN_X) = U_V_d_2(iN_X) + Fm(iV2,iN_X)
        G_V_d_3(iN_X) = U_V_d_3(iN_X) + Fm(iV3,iN_X)

        Gm(iY ,iN_X) = G_Y    (iN_X)
        Gm(iEf,iN_X) = G_Ef   (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

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

        IF ( wMatrRHS(iD) .EQ. 1.0_DP) THEN
          G_D    (iN_X)  = ( 1.0_DP / W ) * C_D  (iN_X) 
        ELSE  
          G_D    (iN_X) &
            = C_D(iN_X)
        END IF

        IF ( wMatrRHS(iY) .EQ. 1.0_DP) THEN
          G_Y    (iN_X) &
            = C_Y    (iN_X) - ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Y * S_Y(iN_X)
        ELSE
          G_Y    (iN_X) &
            = C_Y(iN_X)
        END IF
         
        IF ( wMatrRHS(iEF) .EQ. 1.0_DP) THEN
          G_Ef   (iN_X)  = S_Ef(iN_X) * ( 1.0_DP - W ) / W**2 * ( W + ( W + 1.0_DP) * P(iN_X) / D(iN_X) ) &
                         + 1.0_DP / W * C_Ef(iN_X) - S_Ef(iN_X) / ( W * cD_old(iN_X) ) * SUM_Ef 
        ELSE
          G_Ef   (iN_X) &
            = C_Ef(iN_X)
        END IF

        IF ( wMatrRHS(iV1) .EQ. 1.0_DP) THEN
          G_V_d_1(iN_X)  &
            = 1.0_DP / ( h * W ) * (  C_V_d_1(iN_X) - SUM_V1 * S_V_d_1(iN_X) / cD_old(iN_X) ) 
        ELSE
          G_V_d_1(iN_X) &
            = C_V_d_1(iN_X)
        END IF


        IF ( wMatrRHS(iV2) .EQ. 1.0_DP) THEN
          G_V_d_2(iN_X) &
            = 1.0_DP / ( h * W ) * (  C_V_d_2(iN_X) - SUM_V2 * S_V_d_2(iN_X) / cD_old(iN_X) ) 
        ELSE
          G_V_d_2(iN_X) &
            = C_V_d_2(iN_X)
        END IF


        IF ( wMatrRHS(iV3) .EQ. 1.0_DP) THEN
          G_V_d_3(iN_X) &
            = 1.0_DP / ( h * W ) * (  C_V_d_3(iN_X) - SUM_V3 * S_V_d_3(iN_X) / cD_old(iN_X) ) 
        ELSE
          G_V_d_3(iN_X) &
            = C_V_d_3(iN_X)
        END IF


        

        IF( MoveLeft .EQ. 1) THEN
          IF ( wMatrRHS(iEF) .EQ. 1.0_DP) THEN
            G_Ef(iN_X) = ( 1.0_DP / W ) * C_Ef(iN_X) - S_Ef(iN_X) / ( W * cD_old(iN_X) ) * SUM_Ef &
                       + S_Ef(iN_X) * ( 1.0_DP - W ) / W**2 *  ( W + 1.0_DP) * P(iN_X) / D(iN_X) 
          ELSE
            G_Ef   (iN_X) &
              = C_Ef(iN_X)
          END IF

        END IF

        Gm(iD ,iN_X) = G_D    (iN_X)
        Gm(iY ,iN_X) = G_Y    (iN_X)
        Gm(iEf,iN_X) = G_Ef   (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

        Fm(iD ,iN_X) = G_D    (iN_X) - U_D  (iN_X)
        Fm(iY ,iN_X) = G_Y    (iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = G_Ef   (iN_X) - U_Ef  (iN_X)
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

    INTEGER  :: iN_E, iN_X, iS, iOS, iC
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: R_N, Chi_T, Kappa
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: L_N, L_G1, L_G2, L_G3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          R_N, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS, iC )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          R_N, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$ACC          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS, iC )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          R_N, Chi_T, Kappa, L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS, iC )
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

        ! --- Number Opacity ---

        Chi_T =   Chi_EmAb  (iN_E,iS,iN_X) &
                + Chi_NNS   (iN_E,iS,iN_X) &
                + Chi_NES   (iN_E,iS,iN_X) &
                + Chi_Pair  (iN_E,iS,iN_X) &
                + Chi_NuPair(iN_E,iS,iN_X) &
                + Chi_Brem  (iN_E,iS,iN_X)

        ! --- Number Flux Opacity ---

        iC = nChirals - (LeptonNumber(iS) + 1) / nChirals

        Kappa = Chi_T + Sigma_Iso(iN_E,iC,iN_X)

        ! --- Linear Correction Terms ---

        IF ( Include_LinCorr ) THEN

        L_u_1 =   (L_NNS__In__u_1(iN_E,iS,iN_X)-L_NNS__Out_u_1(iN_E,iS,iN_X)) &
                + (L_NES__In__u_1(iN_E,iS,iN_X)-L_NES__Out_u_1(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_1(iN_E,iS,iN_X)-L_Pair_Ann_u_1(iN_E,iS,iN_X)) &
                - (L_NuPair_Pro_u_1(iN_E,iS,iN_X)-L_NuPair_Ann_u_1(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_1(iN_E,iS,iN_X)-L_Brem_Ann_u_1(iN_E,iS,iN_X))

        L_u_2 =   (L_NNS__In__u_2(iN_E,iS,iN_X)-L_NNS__Out_u_2(iN_E,iS,iN_X)) &
                + (L_NES__In__u_2(iN_E,iS,iN_X)-L_NES__Out_u_2(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_2(iN_E,iS,iN_X)-L_Pair_Ann_u_2(iN_E,iS,iN_X)) &
                - (L_NuPair_Pro_u_2(iN_E,iS,iN_X)-L_NuPair_Ann_u_2(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_2(iN_E,iS,iN_X)-L_Brem_Ann_u_2(iN_E,iS,iN_X))

        L_u_3 =   (L_NNS__In__u_3(iN_E,iS,iN_X)-L_NNS__Out_u_3(iN_E,iS,iN_X)) &
                + (L_NES__In__u_3(iN_E,iS,iN_X)-L_NES__Out_u_3(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_3(iN_E,iS,iN_X)-L_Pair_Ann_u_3(iN_E,iS,iN_X)) &
                - (L_NuPair_Pro_u_3(iN_E,iS,iN_X)-L_NuPair_Ann_u_3(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_3(iN_E,iS,iN_X)-L_Brem_Ann_u_3(iN_E,iS,iN_X))

        L_d_1 = Gm_dd_11(iN_X) &
                  * (   L_NNS__In__u_1(iN_E,iS,iN_X) &
                      + L_NES__In__u_1(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_1(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_1(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_1(iN_E,iS,iN_X) )

        L_d_2 = Gm_dd_22(iN_X) &
                  * (   L_NNS__In__u_2(iN_E,iS,iN_X) &
                      + L_NES__In__u_2(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_2(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_2(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_2(iN_E,iS,iN_X) )

        L_d_3 = Gm_dd_33(iN_X) &
                  * (   L_NNS__In__u_3(iN_E,iS,iN_X) &
                      + L_NES__In__u_3(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_3(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_3(iN_E,iS,iN_X) &
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

        ELSE

          L_N  = Zero
          L_G1 = Zero
          L_G2 = Zero
          L_G3 = Zero

        END IF

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        R_N = Chi_EmAb(iN_E,iS,iN_X) &
                * ( Dnu_0(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) ) &
              + R_NNS(iN_E,iS,iN_X) &
              + R_NES(iN_E,iS,iN_X) &
              + R_Pair(iN_E,iS,iN_X) &
              + R_NuPair(iN_E,iS,iN_X) &
              + R_Brem(iN_E,iS,iN_X) &
              + L_N

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

        Fm(iOS+iCR_N,iN_X) &
          = (   C_Dnu(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) &
              + dt * R_N ) / ( One + dt * Chi_T )

        Gm(iOS+iCR_N,iN_X) &
          = Dnu(iN_E,iS,iN_X) + Fm(iOS+iCR_N,iN_X)

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( C_Inu_d_1(iN_E,iS,iN_X) + dt * L_G1 ) / ( One + dt * Kappa )

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( C_Inu_d_2(iN_E,iS,iN_X) + dt * L_G2 ) / ( One + dt * Kappa )

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( C_Inu_d_3(iN_E,iS,iN_X) + dt * L_G3 ) / ( One + dt * Kappa )

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

    INTEGER  :: iN_E, iN_X, iS, iOS, iC
    REAL(DP) :: vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: R_N, Chi_NonEmAb, Chi_T, Kappa
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: L_N, L_G1, L_G2, L_G3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          R_N, Chi_NonEmAb, Chi_T, Kappa, &
    !$OMP          I_d_1, I_d_2, I_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS, iC )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          R_N, Chi_NonEmAb, Chi_T, Kappa, &
    !$ACC          I_d_1, I_d_2, I_d_3, &
    !$ACC          L_u_1, L_u_2, L_u_3, &
    !$ACC          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS, iC )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          R_N, Chi_NonEmAb, Chi_T, Kappa, &
    !$OMP          I_d_1, I_d_2, I_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, &
    !$OMP          L_d_1, L_d_2, L_d_3, L_N, L_G1, L_G2, L_G3, iOS, iC )
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

        ! --- Number Opacity ---

        Chi_NonEmAb =   Chi_NNS   (iN_E,iS,iN_X) &
                       + Chi_NES   (iN_E,iS,iN_X) &
                       + Chi_Pair  (iN_E,iS,iN_X) &
                       + Chi_NuPair(iN_E,iS,iN_X) &
                       + Chi_Brem  (iN_E,iS,iN_X)

        Chi_T = Chi_EmAb(iN_E,iS,iN_X) + Chi_NonEmAb

        ! --- Number Flux Opacity ---

        iC = nChirals - (LeptonNumber(iS) + 1) / nChirals

        Kappa = Chi_T + Sigma_Iso(iN_E,iC,iN_X)

        ! --- Linear Correction Terms ---

        IF ( Include_LinCorr ) THEN

        L_u_1 =   (L_NNS__In__u_1(iN_E,iS,iN_X)-L_NNS__Out_u_1(iN_E,iS,iN_X)) &
                + (L_NES__In__u_1(iN_E,iS,iN_X)-L_NES__Out_u_1(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_1(iN_E,iS,iN_X)-L_Pair_Ann_u_1(iN_E,iS,iN_X)) &
                - (L_NuPair_Pro_u_1(iN_E,iS,iN_X)-L_NuPair_Ann_u_1(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_1(iN_E,iS,iN_X)-L_Brem_Ann_u_1(iN_E,iS,iN_X))

        L_u_2 =   (L_NNS__In__u_2(iN_E,iS,iN_X)-L_NNS__Out_u_2(iN_E,iS,iN_X)) &
                + (L_NES__In__u_2(iN_E,iS,iN_X)-L_NES__Out_u_2(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_2(iN_E,iS,iN_X)-L_Pair_Ann_u_2(iN_E,iS,iN_X)) &
                - (L_NuPair_Pro_u_2(iN_E,iS,iN_X)-L_NuPair_Ann_u_2(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_2(iN_E,iS,iN_X)-L_Brem_Ann_u_2(iN_E,iS,iN_X))

        L_u_3 =   (L_NNS__In__u_3(iN_E,iS,iN_X)-L_NNS__Out_u_3(iN_E,iS,iN_X)) &
                + (L_NES__In__u_3(iN_E,iS,iN_X)-L_NES__Out_u_3(iN_E,iS,iN_X)) &
                - (L_Pair_Pro_u_3(iN_E,iS,iN_X)-L_Pair_Ann_u_3(iN_E,iS,iN_X)) &
                - (L_NuPair_Pro_u_3(iN_E,iS,iN_X)-L_NuPair_Ann_u_3(iN_E,iS,iN_X)) &
                - (L_Brem_Pro_u_3(iN_E,iS,iN_X)-L_Brem_Ann_u_3(iN_E,iS,iN_X))

        L_d_1 = Gm_dd_11(iN_X) &
                  * (   L_NNS__In__u_1(iN_E,iS,iN_X) &
                      + L_NES__In__u_1(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_1(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_1(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_1(iN_E,iS,iN_X) )

        L_d_2 = Gm_dd_22(iN_X) &
                  * (   L_NNS__In__u_2(iN_E,iS,iN_X) &
                      + L_NES__In__u_2(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_2(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_2(iN_E,iS,iN_X) &
                      - L_Brem_Pro_u_2(iN_E,iS,iN_X) )

        L_d_3 = Gm_dd_33(iN_X) &
                  * (   L_NNS__In__u_3(iN_E,iS,iN_X) &
                      + L_NES__In__u_3(iN_E,iS,iN_X) &
                      - L_Pair_Pro_u_3(iN_E,iS,iN_X) &
                      - L_NuPair_Pro_u_3(iN_E,iS,iN_X) &
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

        ELSE

          L_N  = Zero
          L_G1 = Zero
          L_G2 = Zero
          L_G3 = Zero

        END IF

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! Form the collision source as a sum of process-wise net rates.
        ! The well-balanced rates vanish exactly at their
        ! represented discrete equilibria.

        R_N = Chi_EmAb(iN_E,iS,iN_X) &
                * ( Dnu_0(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) ) &
              + R_NNS(iN_E,iS,iN_X) &
              + R_NES(iN_E,iS,iN_X) &
              + R_Pair(iN_E,iS,iN_X) &
              + R_NuPair(iN_E,iS,iN_X) &
              + R_Brem(iN_E,iS,iN_X) &
              + L_N

        I_d_1 = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)
        I_d_2 = Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)
        I_d_3 = Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)

        ! --- Number Equation ---

        Fm(iOS+iCR_N,iN_X) &
          = Omega(iN_X) &
              * ( C_Dnu(iN_E,iS,iN_X) - vDotInu &
                    - Dnu(iN_E,iS,iN_X) + dt * R_N ) &
            / ( One + Omega(iN_X) * dt * Chi_T )

        Gm(iOS+iCR_N,iN_X) &
          = Dnu(iN_E,iS,iN_X) + Fm(iOS+iCR_N,iN_X)

        ! --- Number Flux 1 Equation ---

        Fm(iOS+iCR_G1,iN_X) &
          = Omega(iN_X) &
              * ( C_Inu_d_1(iN_E,iS,iN_X) - vDotK_d_1 - I_d_1 &
                    + dt * ( L_G1 - Kappa * I_d_1 ) ) &
            / ( One + Omega(iN_X) * dt * Kappa )

        Gm(iOS+iCR_G1,iN_X) = I_d_1 + Fm(iOS+iCR_G1,iN_X)

        ! --- Number Flux 2 Equation ---

        Fm(iOS+iCR_G2,iN_X) &
          = Omega(iN_X) &
              * ( C_Inu_d_2(iN_E,iS,iN_X) - vDotK_d_2 - I_d_2 &
                    + dt * ( L_G2 - Kappa * I_d_2 ) ) &
            / ( One + Omega(iN_X) * dt * Kappa )

        Gm(iOS+iCR_G2,iN_X) = I_d_2 + Fm(iOS+iCR_G2,iN_X)

        ! --- Number Flux 3 Equation ---

        Fm(iOS+iCR_G3,iN_X) &
          = Omega(iN_X) &
              * ( C_Inu_d_3(iN_E,iS,iN_X) - vDotK_d_3 - I_d_3 &
                    + dt * ( L_G3 - Kappa * I_d_3 ) ) &
            / ( One + Omega(iN_X) * dt * Kappa )

        Gm(iOS+iCR_G3,iN_X) = I_d_3 + Fm(iOS+iCR_G3,iN_X)

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


    INTEGER  :: iN_E, iN_X, iS, iOS, iC
    REAL(DP) :: k_dd(3,3), vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, W, vDotV 
    REAL(DP) :: R_N, Chi_T, Kappa
    REAL(DP) :: Inu_d_1, Inu_d_2, Inu_d_3, B_d_1, B_d_2, B_d_3, Det


#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, Det, Inu_d_1, Inu_d_2, Inu_d_3, &
    !$OMP          vDotInu, W, k_dd, vDotV, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          R_N, Chi_T, Kappa, iOS, iC )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( B_d_1, B_d_2, B_d_3, Det, Inu_d_1, Inu_d_2, Inu_d_3, &
    !$ACC          vDotInu, W, k_dd, vDotV, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          R_N, Chi_T, Kappa, iOS, iC )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, Det, Inu_d_1, Inu_d_2, Inu_d_3, &
    !$OMP          vDotInu, W, k_dd, vDotV, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          R_N, Chi_T, Kappa, iOS, iC )
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

        ! --- Number Opacity ---

        Chi_T =   Chi_EmAb  (iN_E,iS,iN_X) &
                + Chi_NNS   (iN_E,iS,iN_X) &
                + Chi_NES   (iN_E,iS,iN_X) &
                + Chi_Pair  (iN_E,iS,iN_X) &
                + Chi_NuPair(iN_E,iS,iN_X) &
                + Chi_Brem  (iN_E,iS,iN_X)

        ! --- Number Flux Opacity ---

        iC = nChirals - (LeptonNumber(iS) + 1) / nChirals

        Kappa = Chi_T + Sigma_Iso(iN_E,iC,iN_X)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        R_N = Chi_EmAb(iN_E,iS,iN_X) &
                * ( Dnu_0(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) ) &
              + R_NNS(iN_E,iS,iN_X) &
              + R_NES(iN_E,iS,iN_X) &
              + R_Pair(iN_E,iS,iN_X) &
              + R_NuPair(iN_E,iS,iN_X) &
              + R_Brem(iN_E,iS,iN_X)

        ! --- Number Equation ---

        Fm(iOS+iCR_N,iN_X) &
          = Omega(iN_X) &
              * (   C_Dnu(iN_E,iS,iN_X) - vDotInu &
                  - W * Dnu(iN_E,iS,iN_X) &
                  + dt * Alpha(iN_X) * R_N ) &
              / ( W + dt * Alpha(iN_X) * Chi_T )

        Gm(iOS+iCR_N,iN_X) &
          = Dnu(iN_E,iS,iN_X) + Fm(iOS+iCR_N,iN_X)

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( One - Omega(iN_X) ) *     Inu_d_1 &
            +       Omega(iN_X)   * ( C_Inu_d_1(iN_E,iS,iN_X) - vDotK_d_1 ) &
                                    / ( W + dt * Alpha(iN_X) * Kappa )

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - Inu_d_1

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( One - Omega(iN_X) ) *     Inu_d_2 &
            +       Omega(iN_X)   * ( C_Inu_d_2(iN_E,iS,iN_X) - vDotK_d_2 ) &
                                    / ( W + dt * Alpha(iN_X) * Kappa )

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - Inu_d_2

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( One - Omega(iN_X) ) *     Inu_d_3 &
            +       Omega(iN_X)   * ( C_Inu_d_3(iN_E,iS,iN_x) - vDotK_d_3 ) &
                                    / ( W + dt * Alpha(iN_X) * Kappa )

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

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FTMP, GTMP

    IF ( Mk == M ) THEN
 
      ALLOCATE( FTMP(1:n_FP,1:M,1:nX_G))
      ALLOCATE( GTMP(1:n_FP,1:M,1:nX_G))

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP MAP( alloc: FTMP, GTMP )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG &
      !$ACC CREATE( FTMP, GTMP )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO 
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
            FTMP(iFP,iM,iN_X) = F(iFP,iM+1,iN_X)
            GTMP(iFP,iM,iN_X) = G(iFP,iM+1,iN_X)
          END DO
          END DO

#if   defined( THORNADO_OMP_OL )
          !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM  = 1, Mk-1
          DO iFP = 1, n_FP
            F(iFP,iM,iN_X) = FTMP(iFP,iM,iN_X)
            G(iFP,iM,iN_X) = GTMP(iFP,iM,iN_X)
          END DO
          END DO

        END IF
      END DO

    END IF

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE CheckConvergence_Inner &
    ( MASK, n_FP, k_inner, nIterations_Inner, &
      nIterations_Inner_Final, RelativeResidual_Inner_First, Fm )

    LOGICAL,  DIMENSION(:)    , INTENT(inout) :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, k_inner
    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Inner
    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Inner_Final
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: RelativeResidual_Inner_First
    REAL(DP), DIMENSION(:,:)  , INTENT(in)    :: Fm

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X, iS, iN_E, iOS
    REAL(DP) :: Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3
    REAL(DP) :: RelativeResidual

    ! --- Compute Norm of F (Use Maximum Across Moments) ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, &
    !$OMP          Fnorm_G3, RelativeResidual )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, &
    !$ACC          Fnorm_G3, RelativeResidual )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, &
    !$OMP          Fnorm_G3, RelativeResidual )
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

          IF( k_inner == 1 )THEN

            Fnorm0(iS,iN_X) = Fnorm

            IF( DnuNorm(iS,iN_X) > Zero )THEN
              RelativeResidual = Fnorm / DnuNorm(iS,iN_X)
            ELSE IF( Fnorm == Zero )THEN
              RelativeResidual = Zero
            ELSE
              RelativeResidual = HUGE(One)
            END IF

            RelativeResidual_Inner_First(iN_X) &
              = MAX( RelativeResidual_Inner_First(iN_X), RelativeResidual )

          END IF

          CONVERGED = CONVERGED .AND. &
                      ( ( Fnorm <= Rtol_inner * DnuNorm(iS,iN_X) ) .OR. &
                        ( Fnorm <= Atol_inner * Fnorm0 (iS,iN_X) ) )

        END DO

        nIterations_Inner(iN_X) = nIterations_Inner(iN_X) + 1
        nIterations_Inner_Final(iN_X) = k_inner

        IF( CONVERGED )THEN
          MASK(iN_X) = .FALSE.
        END IF

      END IF

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( MASK, nIterations_Inner )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( MASK, nIterations_Inner )
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
    !$OMP TARGET UPDATE FROM( MASK_OUTER, MASK_INNER, nIterations_Outer )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( MASK_OUTER, MASK_INNER, nIterations_Outer )
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


  SUBROUTINE CheckError_Inner &
    ( MASK_inner, k_outer, k_inner, &
      D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    USE mpi

    LOGICAL,  DIMENSION(:)    , INTENT(in) :: MASK_inner
    INTEGER,                    INTENT(in) :: k_outer, k_inner
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Y, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: ierr
    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: D_P, T_P, Y_P, E_P, V1_P, V2_P, V3_P
    REAL(DP) :: D0_P, T0_P, Y0_P, E0_P, V10_P, V20_P, V30_P
    REAL(DP) :: Min_E, Min_E_0

    CHARACTER(len=100) :: old_state_filename

    IF ( ( ANY( MASK_inner ) .and. k_inner >= MaxIter_inner ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM &
      !$OMP ( D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !$OMP   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !$OMP   D_old, Y_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
      !$OMP   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST &
      !$ACC ( D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !$ACC   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !$ACC   D_old, Y_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
      !$ACC   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#endif
      DO iN_X = 1, nX_G
        IF ( ( MASK_inner(iN_X) .and. k_inner >= MaxIter_inner ) ) THEN

          WRITE(old_state_filename,*) iN_X
          old_state_filename = "old_state_iN_X_"//TRIM(ADJUSTL(old_state_filename))//".dat"

          OPEN (UNIT=17, FILE=old_state_filename,ACTION='write')

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

          CALL ComputeSpecificInternalEnergy_TABLE(D(iN_X),     Min_T, Y(iN_X),     Min_E)
          CALL ComputeSpecificInternalEnergy_TABLE(D_old(iN_X), Min_T, Y_old(iN_X), Min_E_0)

          WRITE(17,'(7es23.15)') D0_P, T0_P, Y0_P, E0_P, V10_P, V20_P, V30_P
          WRITE(17,'(i5)')       nE_G

          WRITE(*,*)                      '[SolveNeutrinoMatterCoupling_FP_Nested_AA] Non-coverged inner loop'
          WRITE(*,'(a,2i5)')              '              iN_X        : ', iN_X
          WRITE(*,'(a,5x,2i23)')          '         k_outer, k_inner : ', k_outer, k_inner
          WRITE(*,'(a,5x,7es23.15)')      '    D, Y, E, T, V_u       : ', D_P, Y_P, E_P, T_P, V1_P, V2_P, V3_P
          WRITE(*,'(a,5x,7es23.15)')      '    D, Y, E, T, V_u (old) : ', D0_P, Y0_P, E0_P, T0_P, V10_P, V20_P, V30_P
          WRITE(*,'(a,5x,2es23.15)')      '         Min E, Min E old : ', Min_E / Unit_e, Min_E_0 / Unit_E

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Dnu           : ', iS, ( Dnu    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_1       : ', iS, ( Inu_u_1(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_2       : ', iS, ( Inu_u_2(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_3       : ', iS, ( Inu_u_3(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Dnu     (old) : ', iS, ( Dnu_old    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_1 (old) : ', iS, ( Inu_u_1_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_2 (old) : ', iS, ( Inu_u_2_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_3 (old) : ', iS, ( Inu_u_3_old(iN_E,iS,iN_X), iN_E = 1, nE_G )

          WRITE(17,'(i5,100es23.15)') iS,( Dnu_old    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(17,'(i5,100es23.15)') iS,( Inu_u_1_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(17,'(i5,100es23.15)') iS,( Inu_u_2_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(17,'(i5,100es23.15)') iS,( Inu_u_3_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

          CLOSE(17)

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Nnu           : ', iS, ( Nnu    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Gnu_d_1       : ', iS, ( Gnu_d_1(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Gnu_d_2       : ', iS, ( Gnu_d_2(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Gnu_d_3       : ', iS, ( Gnu_d_3(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

        END IF
      END DO
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    END IF

  END SUBROUTINE CheckError_Inner


  SUBROUTINE CheckError_FP &
    ( MASK_outer, MASK_inner, k_outer, k_inner, Error, &
      D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, IgnoreIterationLimit_Option )

    USE mpi

    LOGICAL,  DIMENSION(:)    , INTENT(in) :: MASK_outer, MASK_inner
    INTEGER,                    INTENT(in) :: k_outer, k_inner
    INTEGER,  DIMENSION(:)    , INTENT(in) :: Error
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Y, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    LOGICAL, INTENT(in), OPTIONAL :: IgnoreIterationLimit_Option

    INTEGER  :: ierr
    INTEGER  :: iN_E, iN_X, iS
    LOGICAL  :: IgnoreIterationLimit
    REAL(DP) :: D_P, T_P, Y_P, E_P, V1_P, V2_P, V3_P
    REAL(DP) :: D0_P, T0_P, Y0_P, E0_P, V10_P, V20_P, V30_P
    REAL(DP) :: Min_E, Min_E_0

    CHARACTER(len=100) :: old_state_filename

    IgnoreIterationLimit = .FALSE.
    IF( PRESENT( IgnoreIterationLimit_Option ) ) &
      IgnoreIterationLimit = IgnoreIterationLimit_Option

    IF (        ANY( Error > 0 ) &
         .or. ( .NOT. IgnoreIterationLimit &
                  .AND. ANY( MASK_inner ) .AND. k_inner >= MaxIter_inner ) &
         .or. ( .NOT. IgnoreIterationLimit &
                  .AND. ANY( MASK_outer ) .AND. k_outer >= MaxIter_outer ) ) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM &
      !$OMP ( D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !$OMP   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !$OMP   D_old, Y_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
      !$OMP   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST &
      !$ACC ( D, Y, E, T, V_u_1, V_u_2, V_u_3, &
      !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !$ACC   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      !$ACC   D_old, Y_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
      !$ACC   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#endif
      DO iN_X = 1, nX_G
        IF (        Error(iN_X) > 0 &
             .or. ( .NOT. IgnoreIterationLimit &
                      .AND. MASK_inner(iN_X) .AND. k_inner >= MaxIter_inner ) &
             .or. ( .NOT. IgnoreIterationLimit &
                      .AND. MASK_outer(iN_X) .AND. k_outer >= MaxIter_outer ) ) THEN


          WRITE(old_state_filename,*) iN_X
          old_state_filename = "old_state_iN_X_"//TRIM(ADJUSTL(old_state_filename))//".dat"

          OPEN (UNIT=17, FILE=old_state_filename,ACTION='write')

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

          CALL ComputeSpecificInternalEnergy_TABLE(D(iN_X),     Min_T, Y(iN_X),     Min_E)
          CALL ComputeSpecificInternalEnergy_TABLE(D_old(iN_X), Min_T, Y_old(iN_X), Min_E_0)

          WRITE(17,'(7es23.15)') D0_P, T0_P, Y0_P, E0_P, V10_P, V20_P, V30_P
          WRITE(17,'(i5)')       nE_G

          WRITE(*,*)                      '[SolveNeutrinoMatterCoupling_FP_Nested_AA] Error'
          WRITE(*,'(a,2i5)')              '              iN_X, Error : ', iN_X, Error(iN_X)
          WRITE(*,'(a,5x,2i23)')          '         k_outer, k_inner : ', k_outer, k_inner
          WRITE(*,'(a,5x,7es23.15)')      '    D, Y, E, T, V_u       : ', D_P, Y_P, E_P, T_P, V1_P, V2_P, V3_P
          WRITE(*,'(a,5x,7es23.15)')      '    D, Y, E, T, V_u (old) : ', D0_P, Y0_P, E0_P, T0_P, V10_P, V20_P, V30_P
          WRITE(*,'(a,5x,2es23.15)')      '         Min E, Min E old : ', Min_E / Unit_e, Min_E_0 / Unit_E

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Dnu           : ', iS, ( Dnu    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_1       : ', iS, ( Inu_u_1(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_2       : ', iS, ( Inu_u_2(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_3       : ', iS, ( Inu_u_3(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Dnu     (old) : ', iS, ( Dnu_old    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_1 (old) : ', iS, ( Inu_u_1_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_2 (old) : ', iS, ( Inu_u_2_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Inu_u_3 (old) : ', iS, ( Inu_u_3_old(iN_E,iS,iN_X), iN_E = 1, nE_G )

          WRITE(17,'(i5,100es23.15)') iS,( Dnu_old    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(17,'(i5,100es23.15)') iS,( Inu_u_1_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(17,'(i5,100es23.15)') iS,( Inu_u_2_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(17,'(i5,100es23.15)') iS,( Inu_u_3_old(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

          CLOSE(17)

          DO iS = 1, nSpecies
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Nnu           : ', iS, ( Nnu    (iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Gnu_d_1       : ', iS, ( Gnu_d_1(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Gnu_d_2       : ', iS, ( Gnu_d_2(iN_E,iS,iN_X), iN_E = 1, nE_G )
          WRITE(*,'(a,5x,i5,100es23.15)') '       iS, Gnu_d_3       : ', iS, ( Gnu_d_3(iN_E,iS,iN_X), iN_E = 1, nE_G )
          END DO

        END IF
      END DO
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    END IF

  END SUBROUTINE CheckError_FP


END MODULE TwoMoment_NeutrinoMatterSolverModule
