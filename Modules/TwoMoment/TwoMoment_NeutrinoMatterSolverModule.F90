#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif

MODULE TwoMoment_NeutrinoMatterSolverModule

  USE KindModule, ONLY: &
    DP, Zero, Third, Half, One, Two, Three, FourPi, SqrtTiny
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    BoltzmannConstant, &
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
    Timer_Opacity_EC, &
    Timer_Opacity_ES, &
    Timer_Opacity_NNS, &
    Timer_Opacity_NES, &
    Timer_Opacity_Pair, &
    Timer_Opacity_NuPair, &
    Timer_Opacity_Brem, &
    Timer_Opacity_GeneralLepton, &
    Timer_OpacityRate_NNS, &
    Timer_OpacityRate_NES, &
    Timer_OpacityRate_Pair, &
    Timer_OpacityRate_NuPair, &
    Timer_OpacityRate_Brem, &
    Timer_OpacityRate_GeneralLepton
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
    iNuE, iNuE_Bar, iNuM, iNuM_Bar, &
    iNu, iNu_Bar, LeptonNumber, nChirals, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeMuonChemicalPotential_TABLE, &
    ComputeMuonFractionFromDegeneracy_TABLE, &
    ComputePressure_TABLE, &
    Min_D, Min_T, Min_Y, Max_T
  USE OpacityModule_TABLE, ONLY: &
    IsoIncludesNucleonScattering
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, &
    ComputeNeutrinoOpacities_EC, &
    WriteMuonChargedCurrentOpacityDiagnostic, &
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
    ComputeDetailedBalanceFactors, &
    SetOpacityMask, &
    iM0, iM1, &
    nOp, iOp_EmAb_Electron, iOp_EmAb_Muon, &
    iOp_ECTable, iOp_Iso, &
    iOp_NNS, iOp_NES, iOp_Pair, iOp_NuPair, iOp_Brem, &
    iOp_NMS, iOp_LFE, iOp_LFC, iOp_IMD
  USE NeutrinoLeptonInteractionsModule, ONLY: &
    InitializeNeutrinoLeptonInteractions, &
    FinalizeNeutrinoLeptonInteractions, &
    FinalizeNeutrinoLeptonInteractionKernels, &
    ComputeNeutrinoLeptonKernels, ComputeNeutrinoLeptonRates, &
    Chi_GL, R_GL, nGL, iGL_NMS, iGL_LFE, iGL_LFC, iGL_IMD
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
  PUBLIC :: FinalizeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling_FP_Nested_AA

  INTEGER, PARAMETER :: iReducedSolver_Anderson    = 1
  INTEGER, PARAMETER :: iReducedSolver_SchurNewton = 2
  INTEGER, PARAMETER :: iReducedSolver_CoupledEtaSchur = 3

  ! Matter collision RHS quantities used by the shared limiter accounting.
  INTEGER, PARAMETER :: iQ_Le = 1
  INTEGER, PARAMETER :: iQ_E  = 2
  INTEGER, PARAMETER :: iQ_H  = 3
  INTEGER, PARAMETER :: iQ_M1 = 4
  INTEGER, PARAMETER :: iQ_M2 = 5
  INTEGER, PARAMETER :: iQ_M3 = 6
  INTEGER, PARAMETER :: iQ_Lm = 7
  INTEGER, PARAMETER :: nQ    = 7

#if defined( TWOMOMENT_ORDER_V )
  PUBLIC :: ComputeOpacities_Packed
  PUBLIC :: ComputeRates_Packed
  PUBLIC :: ComputeTransportOpacity
  PUBLIC :: AssembleTransportOpacity
  PUBLIC :: ComputeMatterRHS_Collisions_OrderV
  PUBLIC :: ComputeMeanInverseMeanFreePaths
  PUBLIC :: Include_NNS, Include_NES, Include_Pair, Include_NuPair, Include_Brem
  PUBLIC :: Include_NMS, Include_LFE, Include_LFC, Include_IMD
  PUBLIC :: nQ
  PUBLIC :: iQ_Le, iQ_Lm, iQ_E, iQ_H, iQ_M1, iQ_M2, iQ_M3
  PUBLIC :: nNNSRateDecomposition, iNNSRateDecomposition_NIn, iNNSRateDecomposition_NOut, iNNSRateDecomposition_NNet, &
            iNNSRateDecomposition_EIn, iNNSRateDecomposition_EOut, iNNSRateDecomposition_QH
  PUBLIC :: nEmAbRateDecomposition, iEmAbRateDecomposition_NIn, iEmAbRateDecomposition_NOut, &
            iEmAbRateDecomposition_NNet, iEmAbRateDecomposition_EIn, iEmAbRateDecomposition_EOut, &
            iEmAbRateDecomposition_QH
  PUBLIC :: nEmAbSpectrum, iEmAbSpectrum_Dnu, iEmAbSpectrum_DnuEq, &
            iEmAbSpectrum_Chi, iEmAbSpectrum_EIn, &
            iEmAbSpectrum_EOut, iEmAbSpectrum_QH
  PUBLIC :: nC
  PUBLIC :: iC_EmAb, iC_Iso, iC_NNS, iC_NES
  PUBLIC :: iC_Pair, iC_NuPair, iC_Brem
  PUBLIC :: iC_NMS, iC_LFE, iC_LFC, iC_IMD

  INTEGER, PARAMETER :: iNNSRateDecomposition_NIn  = 1
  INTEGER, PARAMETER :: iNNSRateDecomposition_NOut = 2
  INTEGER, PARAMETER :: iNNSRateDecomposition_NNet = 3
  INTEGER, PARAMETER :: iNNSRateDecomposition_EIn  = 4
  INTEGER, PARAMETER :: iNNSRateDecomposition_EOut = 5
  INTEGER, PARAMETER :: iNNSRateDecomposition_QH   = 6
  INTEGER, PARAMETER :: nNNSRateDecomposition      = 6

  INTEGER, PARAMETER :: iEmAbRateDecomposition_NIn  = 1
  INTEGER, PARAMETER :: iEmAbRateDecomposition_NOut = 2
  INTEGER, PARAMETER :: iEmAbRateDecomposition_NNet = 3
  INTEGER, PARAMETER :: iEmAbRateDecomposition_EIn  = 4
  INTEGER, PARAMETER :: iEmAbRateDecomposition_EOut = 5
  INTEGER, PARAMETER :: iEmAbRateDecomposition_QH   = 6
  INTEGER, PARAMETER :: nEmAbRateDecomposition      = 6

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
  INTEGER, PARAMETER :: iC_NMS    = 8
  INTEGER, PARAMETER :: iC_LFE    = 9
  INTEGER, PARAMETER :: iC_LFC    = 10
  INTEGER, PARAMETER :: iC_IMD    = 11
  INTEGER, PARAMETER :: nC        = 11

  INTEGER, PARAMETER, PUBLIC :: nLeptonFamilies = 3
  INTEGER, PARAMETER, PUBLIC :: iLFC_InputConversion  = 1
  INTEGER, PARAMETER, PUBLIC :: iLFC_EntryProjection  = 2
  INTEGER, PARAMETER, PUBLIC :: iLFC_CollisionSource  = 3
  INTEGER, PARAMETER, PUBLIC :: iLFC_RawMap           = 4
  INTEGER, PARAMETER, PUBLIC :: iLFC_MapProjection    = 5
  INTEGER, PARAMETER, PUBLIC :: iLFC_FixedPointResidual = 6
  INTEGER, PARAMETER, PUBLIC :: iLFC_EquationResidual = 7
  INTEGER, PARAMETER, PUBLIC :: nLeptonFlavorCheck    = 7
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
  INTEGER, PARAMETER :: iYe = 1
  INTEGER, PARAMETER :: iYm = 2
  INTEGER, PARAMETER :: iEf = 3
  INTEGER, PARAMETER :: iV1 = 1
  INTEGER, PARAMETER :: iV2 = 1
  INTEGER, PARAMETER :: iV3 = 1
  INTEGER, PARAMETER :: nMatterEquations = 3
#elif defined( TWOMOMENT_ORDER_V )
  INTEGER, PARAMETER :: iD  = 1
  INTEGER, PARAMETER :: iYe = 1
  INTEGER, PARAMETER :: iYm = 2
  INTEGER, PARAMETER :: iEf = 3
  INTEGER, PARAMETER :: iV1 = 4
  INTEGER, PARAMETER :: iV2 = 5
  INTEGER, PARAMETER :: iV3 = 6
  INTEGER, PARAMETER :: nMatterEquations = 6
#elif defined( TWOMOMENT_RELATIVISTIC )
  INTEGER, PARAMETER :: iD  = 1
  INTEGER, PARAMETER :: iYe = 2
  INTEGER, PARAMETER :: iYm = 3
  INTEGER, PARAMETER :: iEf = 4
  INTEGER, PARAMETER :: iV1 = 5
  INTEGER, PARAMETER :: iV2 = 6
  INTEGER, PARAMETER :: iV3 = 7
  INTEGER, PARAMETER :: nMatterEquations = 7
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
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Ye_old, C_Ye, S_Ye, G_Ye, U_Ye
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Ym_old, C_Ym, S_Ym, G_Ym, U_Ym
  REAL(DP), DIMENSION(:), ALLOCATABLE :: T_old
  REAL(DP), DIMENSION(:), ALLOCATABLE :: V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: D_old, cD_old, C_D, S_D, G_D, U_D

  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: SqrtGm, T_Opacity
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
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: &
    DetailedBalanceScattering_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: &
    DetailedBalancePair_T

  LOGICAL,  DIMENSION(:,:),   ALLOCATABLE, TARGET :: Opacity_Mask
  LOGICAL,  DIMENSION(:,:),   ALLOCATABLE, TARGET :: MatterSource_Mask

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
  LOGICAL  :: Include_NMS
  LOGICAL  :: Include_LFE
  LOGICAL  :: Include_LFC
  LOGICAL  :: Include_IMD
  INTEGER  :: GeneralLeptonQuadratureOrder
  LOGICAL  :: Include_LinCorr
  LOGICAL  :: UseFirstMomentRates
  LOGICAL  :: UseConservativeProjection
  LOGICAL  :: UseRadiationOnlySolve
  LOGICAL  :: UseMuonLeptonSource
  INTEGER  :: ReducedSolverMethod
  LOGICAL  :: DebugMuonCoupling
  LOGICAL  :: TraceRadiationOnlySolve
  LOGICAL  :: AnalyzeRadiationOnlyJacobian
  CHARACTER(512) :: RadiationOnlyTraceFile
  INTEGER :: RadiationOnlySolveCallCount = 0
  INTEGER :: RadiationOnlyTraceUnit = 0
  INTEGER :: RadiationOnlyTraceSpectrumUnit = 0
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
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: D_T, T_T, Ye_T, Ym_T, E_T
  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: T0_T, Ye0_T, Ym0_T
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
    CALL InitializeNeutrinoLeptonInteractions &
           ( E_N, nX_G, nTheta_Option = GeneralLeptonQuadratureOrder, &
             EnableGeneralInteractions_Option &
               = Include_NMS .OR. Include_LFE .OR. &
                 Include_LFC .OR. Include_IMD )

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

    ALLOCATE(   Ye_old(nX_G) )
    ALLOCATE( C_Ye    (nX_G) )
    ALLOCATE( S_Ye    (nX_G) )
    ALLOCATE( G_Ye    (nX_G) )
    ALLOCATE( U_Ye    (nX_G) )

    ALLOCATE(   Ym_old(nX_G) )
    ALLOCATE( C_Ym    (nX_G) )
    ALLOCATE( S_Ym    (nX_G) )
    ALLOCATE( G_Ym    (nX_G) )
    ALLOCATE( U_Ym    (nX_G) )

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
    ALLOCATE(       T_Opacity(              nX_G) )
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
    ALLOCATE( MatterSource_Mask(nMatterEquations,nX_G) )

    ALLOCATE( Error_T(nX_G) )

    ALLOCATE( D_T (nX_G) )
    ALLOCATE( T_T (nX_G) )
    ALLOCATE( Ye_T(nX_G) )
    ALLOCATE( Ym_T(nX_G) )
    ALLOCATE( E_T (nX_G) )

    ALLOCATE( T0_T (nX_G) )
    ALLOCATE( Ye0_T(nX_G) )
    ALLOCATE( Ym0_T(nX_G) )

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
    ALLOCATE( DetailedBalanceScattering_T(nE_G,nE_G,nX_G) )
    ALLOCATE( DetailedBalancePair_T(nE_G,nE_G,nX_G) )

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
    !$OMP             Ye_old, C_Ye, S_Ye, G_Ye, U_Ye, &
    !$OMP             Ym_old, C_Ym, S_Ym, G_Ym, U_Ym, &
    !$OMP             T_old, &
    !$OMP             V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$OMP             V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$OMP             V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$OMP             SqrtGm, T_Opacity, &
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
    !$OMP             D_T, T_T, Ye_T, Ym_T, E_T, Error_T,  &
    !$OMP             T0_T, Ye0_T, Ym0_T, &
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
    !$OMP             DetailedBalanceScattering_T, &
    !$OMP             DetailedBalancePair_T, &
    !$OMP             Opacity_Mask, Opacity_Mask_T, MatterSource_Mask )
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
    !$ACC         Ye_old, C_Ye, S_Ye, G_Ye, U_Ye, &
    !$ACC         Ym_old, C_Ym, S_Ym, G_Ym, U_Ym, &
    !$ACC         T_old, &
    !$ACC         V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$ACC         V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$ACC         V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$ACC         SqrtGm, T_Opacity, &
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
    !$ACC         D_T, T_T, Ye_T, Ym_T, E_T, Error_T, &
    !$ACC         T0_T, Ye0_T, Ym0_T, &
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
    !$ACC         DetailedBalanceScattering_T, &
    !$ACC         DetailedBalancePair_T, &
    !$ACC         Opacity_Mask, Opacity_Mask_T, MatterSource_Mask )
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
      Include_NuPair_Option, Include_Brem_Option, &
      Include_NMS_Option, Include_LFE_Option, Include_LFC_Option, &
      Include_IMD_Option, GeneralLeptonQuadratureOrder_Option, &
      Include_LinCorr_Option, &
      wMatrRHS_Option, DnuMax_Option, FreezeOpacities_Option, &
      UseConservativeProjection_Option, Verbose_Option, &
      UseRadiationOnlySolve_Option, DebugMuonCoupling_Option, &
      TraceRadiationOnlySolve_Option, RadiationOnlyTraceFile_Option, &
      AnalyzeRadiationOnlyJacobian_Option, ReducedSolverMethod_Option, &
      UseMuonLeptonSource_Option )

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
    LOGICAL , INTENT(in), OPTIONAL :: Include_NMS_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_LFE_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_LFC_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_IMD_Option
    INTEGER , INTENT(in), OPTIONAL :: GeneralLeptonQuadratureOrder_Option
    LOGICAL , INTENT(in), OPTIONAL :: Include_LinCorr_Option
    REAL(DP), INTENT(in), OPTIONAL :: wMatrRHS_Option(nMatterEquations)
    REAL(DP), INTENT(in), OPTIONAL :: DnuMax_Option
    LOGICAL , INTENT(in), OPTIONAL :: FreezeOpacities_Option
    LOGICAL , INTENT(in), OPTIONAL :: UseConservativeProjection_Option
    LOGICAL , INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL , INTENT(in), OPTIONAL :: UseRadiationOnlySolve_Option
    LOGICAL , INTENT(in), OPTIONAL :: DebugMuonCoupling_Option
    LOGICAL , INTENT(in), OPTIONAL :: TraceRadiationOnlySolve_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: RadiationOnlyTraceFile_Option
    LOGICAL , INTENT(in), OPTIONAL :: AnalyzeRadiationOnlyJacobian_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: ReducedSolverMethod_Option
    LOGICAL, INTENT(in), OPTIONAL :: UseMuonLeptonSource_Option

    LOGICAL :: Verbose, MuonSourceSupported
    INTEGER :: ierr, myid

    IF( SolverParametersInitialized ) RETURN

    RadiationOnlySolveCallCount = 0

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

    Include_NMS = .FALSE.
    Include_LFE = .FALSE.
    Include_LFC = .FALSE.
    Include_IMD = .FALSE.
    IF( PRESENT(Include_NMS_Option) ) Include_NMS = Include_NMS_Option
    IF( PRESENT(Include_LFE_Option) ) Include_LFE = Include_LFE_Option
    IF( PRESENT(Include_LFC_Option) ) Include_LFC = Include_LFC_Option
    IF( PRESENT(Include_IMD_Option) ) Include_IMD = Include_IMD_Option

    GeneralLeptonQuadratureOrder = 16
    IF( PRESENT(GeneralLeptonQuadratureOrder_Option) ) &
      GeneralLeptonQuadratureOrder = GeneralLeptonQuadratureOrder_Option
    IF( GeneralLeptonQuadratureOrder < 1 )THEN
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      IF( myid == 0 )THEN
        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: General-lepton quadrature order must be positive.'
        FLUSH(ERROR_UNIT)
      END IF
      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
    END IF

    IF( PRESENT( Include_LinCorr_Option ) )THEN
      Include_LinCorr = Include_LinCorr_Option
    ELSE
      Include_LinCorr = .FALSE.
    END IF

    IF( Include_LinCorr .AND. &
        ( Include_NMS .OR. Include_LFE .OR. Include_LFC .OR. Include_IMD ) )THEN
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      IF( myid == 0 )THEN
        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: General lepton reactions currently provide only zeroth-moment kernels.'
        WRITE(ERROR_UNIT,'(A)') &
          '       Disable Include_LinCorr or implement their first angular moments.'
        FLUSH(ERROR_UNIT)
      END IF
      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
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

    ReducedSolverMethod = iReducedSolver_Anderson
    IF( PRESENT(ReducedSolverMethod_Option) )THEN
      SELECT CASE( TRIM(ADJUSTL(ReducedSolverMethod_Option)) )
        CASE( 'anderson', 'ANDERSON' )
          ReducedSolverMethod = iReducedSolver_Anderson
        CASE( 'schur-newton', 'SCHUR-NEWTON', &
              'schur_newton', 'SCHUR_NEWTON' )
          ReducedSolverMethod = iReducedSolver_SchurNewton
        CASE( 'coupled-eta-schur', 'COUPLED-ETA-SCHUR' )
          ReducedSolverMethod = iReducedSolver_CoupledEtaSchur
#if !defined(EOSMODE_COMPOSE) || !defined(INTERPOLATION_SPLIT_TABLE_COMBINED)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
          IF( myid == 0 ) WRITE(ERROR_UNIT,'(A)') &
            'coupled-eta-schur requires EOSMODE=COMPOSE, INTERPOLATION_SPLIT_TABLE=COMBINED.'
          CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
#endif
        CASE DEFAULT
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          IF( myid == 0 )THEN
            WRITE(ERROR_UNIT,'(A,A)') &
              'ERROR: Unknown reduced nonlinear method: ', &
              TRIM(ReducedSolverMethod_Option)
            WRITE(ERROR_UNIT,'(A)') &
              '       Expected anderson, schur-newton, or coupled-eta-schur.'
            FLUSH(ERROR_UNIT)
          END IF
          CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
      END SELECT
    END IF

    UseMuonLeptonSource = .FALSE.
    IF( PRESENT( UseMuonLeptonSource_Option ) ) &
      UseMuonLeptonSource = UseMuonLeptonSource_Option

    MuonSourceSupported = .FALSE.
#if defined( TWOMOMENT_ORDER_V )
    MuonSourceSupported = .NOT. UseConservativeProjection .AND. &
      ( .NOT. UseRadiationOnlySolve .OR. &
        ReducedSolverMethod == iReducedSolver_CoupledEtaSchur )
#endif
    IF( UseMuonLeptonSource )THEN
      IF( .NOT. MuonSourceSupported )THEN
        CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
        IF( myid == 0 ) WRITE(ERROR_UNIT,'(A)') &
          'Muon source update requires Order-v nested or coupled-eta-schur, without terminal projection.'
        CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
      END IF
    END IF

    IF( ReducedSolverMethod == iReducedSolver_CoupledEtaSchur )THEN
      IF( .NOT. UseRadiationOnlySolve )THEN
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        IF( myid == 0 ) WRITE(ERROR_UNIT,'(A)') &
          'coupled-eta-schur requires UseRadiationOnlySolve=.TRUE.'
        CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      END IF
    END IF

    IF( PRESENT( DebugMuonCoupling_Option ) )THEN
      DebugMuonCoupling = DebugMuonCoupling_Option
    ELSE
      DebugMuonCoupling = .FALSE.
    END IF

    IF( PRESENT( TraceRadiationOnlySolve_Option ) )THEN
      TraceRadiationOnlySolve = TraceRadiationOnlySolve_Option
    ELSE
      TraceRadiationOnlySolve = .FALSE.
    END IF

    RadiationOnlyTraceFile = 'reduced_solver_trace'
    IF( PRESENT( RadiationOnlyTraceFile_Option ) ) &
      RadiationOnlyTraceFile = TRIM( RadiationOnlyTraceFile_Option )

    IF( PRESENT( AnalyzeRadiationOnlyJacobian_Option ) )THEN
      AnalyzeRadiationOnlyJacobian = AnalyzeRadiationOnlyJacobian_Option
    ELSE
      AnalyzeRadiationOnlyJacobian = .FALSE.
    END IF

    IF( AnalyzeRadiationOnlyJacobian .AND. .NOT. UseRadiationOnlySolve )THEN

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

      IF( myid == 0 )THEN
        WRITE(ERROR_UNIT,'(A)') &
          'ERROR: Reduced-solver Jacobian analysis requires the'
        WRITE(ERROR_UNIT,'(A)') &
          '       invariant-constrained radiation solve.'
        FLUSH(ERROR_UNIT)
      END IF

      CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

    END IF

    IF( AnalyzeRadiationOnlyJacobian .AND. &
        ReducedSolverMethod == iReducedSolver_CoupledEtaSchur )THEN
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      IF( myid == 0 ) WRITE(ERROR_UNIT,'(A)') &
        'The eliminated-Ymu Jacobian diagnostic does not apply to coupled-eta-schur.'
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
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
      WRITE(*,'(A4,A32,L1)')       '', 'Direct muon-lepton source: ', &
        UseMuonLeptonSource
      IF( ReducedSolverMethod == iReducedSolver_Anderson )THEN
        WRITE(*,'(A4,A32,A)') '', 'Reduced nonlinear method: ', 'anderson'
      ELSE IF( ReducedSolverMethod == iReducedSolver_SchurNewton )THEN
        WRITE(*,'(A4,A32,A)') '', 'Reduced nonlinear method: ', &
          'schur-newton'
      ELSE
        WRITE(*,'(A4,A32,A)') '', 'Reduced nonlinear method: ', &
          'coupled-eta-schur'
      END IF
      WRITE(*,'(A4,A32,L1)')       '', 'Debug muon coupling: ', &
        DebugMuonCoupling
      WRITE(*,'(A4,A32,L1)')       '', 'Trace invariant solve: ', &
        TraceRadiationOnlySolve
      WRITE(*,'(A4,A32,L1)')       '', 'Analyze reduced Jacobian: ', &
        AnalyzeRadiationOnlyJacobian
      IF( TraceRadiationOnlySolve ) &
        WRITE(*,'(A4,A32,A)') '', 'Trace file base: ', &
          TRIM(RadiationOnlyTraceFile)
      WRITE(*,*)
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NNS: '    , Include_NNS
      WRITE(*,'(A4,A32,L1)')       '', 'NNS_manybody: '   , NNS_ApplyManyBodyCorrection
      WRITE(*,'(A4,A32,L1)')       '', 'NNS_weak_mag: '   , NNS_ApplyWeakMagnetism
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NES: '    , Include_NES
      WRITE(*,'(A4,A32,L1)')       '', 'Include_Pair: '   , Include_Pair
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NuPair: ' , Include_NuPair
      WRITE(*,'(A4,A32,L1)')       '', 'Include_Brem: '   , Include_Brem
      WRITE(*,'(A4,A32,L1)')       '', 'Include_NMS: '    , Include_NMS
      WRITE(*,'(A4,A32,L1)')       '', 'Include_LFE: '    , Include_LFE
      WRITE(*,'(A4,A32,L1)')       '', 'Include_LFC: '    , Include_LFC
      WRITE(*,'(A4,A32,L1)')       '', 'Include_IMD: '    , Include_IMD
      WRITE(*,'(A4,A32,I0)')       '', 'General-lepton quadrature: ', &
        GeneralLeptonQuadratureOrder
      WRITE(*,'(A4,A32,L1)')       '', 'Include_LinCorr: ', Include_LinCorr
      WRITE(*,*)
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iYe): '  , INT(wMatrRHS(iYe))
      WRITE(*,'(A4,A32,I1.1)')     '', 'wMatrRHS(iYm): '  , INT(wMatrRHS(iYm))
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


  SUBROUTINE FinalizeNeutrinoMatterSolverParameters

    CALL FinalizeNeutrinoLeptonInteractionKernels
    SolverParametersInitialized = .FALSE.

  END SUBROUTINE FinalizeNeutrinoMatterSolverParameters


  SUBROUTINE FinalizeNeutrinoMatterSolver

    CALL FinalizeNeutrinoLeptonInteractions

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
    !$OMP               Ye_old, C_Ye, S_Ye, G_Ye, U_Ye, &
    !$OMP               Ym_old, C_Ym, S_Ym, G_Ym, U_Ym, &
    !$OMP               T_old, &
    !$OMP               V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$OMP               V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$OMP               V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$OMP               SqrtGm, T_Opacity, &
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
    !$OMP               D_T, T_T, Ye_T, Ym_T, E_T, Error_T,  &
    !$OMP               T0_T, Ye0_T, Ym0_T, &
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
    !$OMP               DetailedBalanceScattering_T, &
    !$OMP               DetailedBalancePair_T, &
    !$OMP               Opacity_Mask, Opacity_Mask_T, MatterSource_Mask )
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
    !$ACC         Ye_old, C_Ye, S_Ye, G_Ye, U_Ye, &
    !$ACC         Ym_old, C_Ym, S_Ym, G_Ym, U_Ym, &
    !$ACC         T_old, &
    !$ACC         V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, &
    !$ACC         V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, &
    !$ACC         V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, &
    !$ACC         SqrtGm, T_Opacity, &
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
    !$ACC         D_T, T_T, Ye_T, Ym_T, E_T, Error_T,  &
    !$ACC         T0_T, Ye0_T, Ym0_T, &
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
    !$ACC         DetailedBalanceScattering_T, &
    !$ACC         DetailedBalancePair_T, &
    !$ACC         Opacity_Mask, Opacity_Mask_T, MatterSource_Mask )
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
    DEALLOCATE( Ye_old, C_Ye, S_Ye, G_Ye, U_Ye )
    DEALLOCATE( Ym_old, C_Ym, S_Ym, G_Ym, U_Ym )
    DEALLOCATE( T_old )
    DEALLOCATE( V_u_1_old, C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1 )
    DEALLOCATE( V_u_2_old, C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2 )
    DEALLOCATE( V_u_3_old, C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3 )
    DEALLOCATE( cD_old, D_old, C_D, S_D, G_D, U_D )
    DEALLOCATE( SqrtGm, T_Opacity )
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
    DEALLOCATE( D_T, T_T, Ye_T, Ym_T, E_T, SqrtGm_T, Error_T )
    DEALLOCATE( T0_T, Ye0_T, Ym0_T )
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
    DEALLOCATE( DetailedBalanceScattering_T )
    DEALLOCATE( DetailedBalancePair_T )
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
    DEALLOCATE( Opacity_Mask, Opacity_Mask_T, MatterSource_Mask )

    IF( RadiationOnlyTraceUnit /= 0 )THEN
      CLOSE( RadiationOnlyTraceUnit )
      RadiationOnlyTraceUnit = 0
    END IF
    IF( RadiationOnlyTraceSpectrumUnit /= 0 )THEN
      CLOSE( RadiationOnlyTraceSpectrumUnit )
      RadiationOnlyTraceSpectrumUnit = 0
    END IF

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
    ( D, T, Ye, Ym, SqrtGm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      TransportOpacity, PatchHasCollisions_Option )

    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

    REAL(DP), INTENT(in)  :: D(:), T(:), Ye(:), Ym(:), SqrtGm(:)
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

    CALL SetOpacityMask( nX_G, nDOFX, D, T, Ye, Ym, Opacity_Mask )

    PatchHasCollisions &
      =     ( nSpecies >= iNuE_Bar .AND. &
              ANY( Opacity_Mask(iOp_EmAb_Electron,:) ) ) &
        .OR. ( nSpecies >= iNuM_Bar .AND. &
              ANY( Opacity_Mask(iOp_EmAb_Muon,:) ) ) &
        .OR. ( nSpecies >= iNuE_Bar .AND. &
              ANY( Opacity_Mask(iOp_ECTable,:) ) ) &
        .OR. ANY( Opacity_Mask(iOp_Iso,    :) ) &
        .OR. ( Include_NNS    .AND. ANY( Opacity_Mask(iOp_NNS,   :) ) ) &
        .OR. ( Include_NES    .AND. ANY( Opacity_Mask(iOp_NES,   :) ) ) &
        .OR. ( Include_Pair   .AND. ANY( Opacity_Mask(iOp_Pair,  :) ) ) &
        .OR. ( Include_NuPair .AND. ANY( Opacity_Mask(iOp_NuPair,:) ) ) &
        .OR. ( Include_Brem   .AND. ANY( Opacity_Mask(iOp_Brem,  :) ) ) &
        .OR. ( Include_NMS    .AND. ANY( Opacity_Mask(iOp_NMS,   :) ) ) &
        .OR. ( Include_LFE    .AND. ANY( Opacity_Mask(iOp_LFE,   :) ) ) &
        .OR. ( Include_LFC    .AND. ANY( Opacity_Mask(iOp_LFC,   :) ) ) &
        .OR. ( Include_IMD    .AND. ANY( Opacity_Mask(iOp_IMD,   :) ) )

    IF( .NOT. PatchHasCollisions )THEN

      IF( PRESENT( PatchHasCollisions_Option ) ) &
        PatchHasCollisions_Option = .FALSE.
      CALL ZeroTransportOpacity( TransportOpacity )
      RETURN

    END IF

    CALL ComputeOpacities_Packed &
           ( D, T, Ye, Ym, SqrtGm, nX_P0 = nX_G, &
             T0 = T, Ye0 = Ye, Ym0 = Ym, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL ComputeRates_Packed &
           ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, nX_P0 = nX_G, &
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
    LOGICAL :: UseNMS, UseLFE, UseLFC, UseIMD

    UseNNS    = Include_NNS
    UseNES    = Include_NES
    UsePair   = Include_Pair
    UseNuPair = Include_NuPair
    UseBrem   = Include_Brem
    UseNMS    = Include_NMS
    UseLFE    = Include_LFE
    UseLFC    = Include_LFC
    UseIMD    = Include_IMD

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iChiral )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iChiral ) &
    !$ACC FIRSTPRIVATE( UseNNS, UseNES, UsePair, UseNuPair, UseBrem, &
    !$ACC               UseNMS, UseLFE, UseLFC, UseIMD ) &
    !$ACC PRESENT( TransportOpacity, Sigma_Iso, Chi_EmAb, Chi_NNS, &
    !$ACC          Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, Chi_GL )
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

      IF( UseNMS ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) &
              + Chi_GL(iN_E,iS,iN_X,iGL_NMS)

      IF( UseLFE ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) &
              + Chi_GL(iN_E,iS,iN_X,iGL_LFE)

      IF( UseLFC ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) &
              + Chi_GL(iN_E,iS,iN_X,iGL_LFC)

      IF( UseIMD ) &
        TransportOpacity(iN_E,iS,iN_X) &
          = TransportOpacity(iN_E,iS,iN_X) &
              + Chi_GL(iN_E,iS,iN_X,iGL_IMD)

    END DO
    END DO
    END DO

  END SUBROUTINE AssembleTransportOpacity


  SUBROUTINE ComputeMeanInverseMeanFreePaths &
    ( D, T, Ye, Ym, SqrtGm, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, MeanInvMFP )

    REAL(DP), INTENT(in)  :: D(:), T(:), Ye(:), Ym(:), SqrtGm(:)
    REAL(DP), INTENT(in)  :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP), INTENT(in)  :: Dnu(:,:,:)
    REAL(DP), INTENT(in)  :: Inu_u_1(:,:,:), Inu_u_2(:,:,:), Inu_u_3(:,:,:)
    REAL(DP), INTENT(out) :: MeanInvMFP(:,:,:)

    INTEGER :: iChiral, iN_E, iN_X, iProcess, iS
    REAL(DP) :: Denominator, IdotI, IdotL, InvW2_Ref, WeightE
    REAL(DP) :: Numerator(nC)
    CALL SetOpacityMask( nX_G, nDOFX, D, T, Ye, Ym, Opacity_Mask )

    ! This is an instantaneous diagnostic.  The final state is therefore
    ! both the evaluation state and the EC-table reference state.
    CALL ComputeOpacities_Packed &
           ( D, T, Ye, Ym, SqrtGm, nX_P0 = nX_G, &
             T0 = T, Ye0 = Ye, Ym0 = Ym, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL ComputeRates_Packed &
           ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, nX_P0 = nX_G, &
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

        IF( Include_NMS ) &
          Numerator(iC_NMS) = Numerator(iC_NMS) + WeightE &
            * Dnu(iN_E,iS,iN_X) * Chi_GL(iN_E,iS,iN_X,iGL_NMS)

        IF( Include_LFE ) &
          Numerator(iC_LFE) = Numerator(iC_LFE) + WeightE &
            * Dnu(iN_E,iS,iN_X) * Chi_GL(iN_E,iS,iN_X,iGL_LFE)

        IF( Include_LFC ) &
          Numerator(iC_LFC) = Numerator(iC_LFC) + WeightE &
            * Dnu(iN_E,iS,iN_X) * Chi_GL(iN_E,iS,iN_X,iGL_LFC)

        IF( Include_IMD ) &
          Numerator(iC_IMD) = Numerator(iC_IMD) + WeightE &
            * Dnu(iN_E,iS,iN_X) * Chi_GL(iN_E,iS,iN_X,iGL_IMD)

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
    ( dt, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, D, T, Ye, Ym, E, &
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
        LeptonFlavorCheck_Rate_Final &
#endif
      )

    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE mpi

    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: D, T, Ye, Ym, E
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
    REAL(DP), DIMENSION(:,:,:), INTENT(out)   :: LeptonFlavorCheck_Rate_Final
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
      ProjectionMask, RadiationBacktrackMask, &
      TraceEvaluationMask, TraceHasPrevious, TraceHasPrevious2, &
      SchurHasPrevious
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TraceBacktrackCount
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: &
      T_Conservative, Ye_Conservative, Ym_Conservative, E_Conservative, &
      V_u_1_Conservative, V_u_2_Conservative, V_u_3_Conservative
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: &
      FVECm_projection, GVECm_projection, &
      InitialRadiationLimiter_ConservationCorrection, &
      RadiationLeptonInvariant, &
      TracePreviousState, TracePrevious2State, TracePreviousResidual, &
      TraceAlpha
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: &
      SchurPreviousInput, SchurPreviousMap, SchurPreviousResidual, &
      SchurGain, SchurDamping, &
      T_RadiationRecoveryBase, Ye_RadiationRecoveryBase, &
      Ym_RadiationRecoveryBase, E_RadiationRecoveryBase, &
      V1_RadiationRecoveryBase, V2_RadiationRecoveryBase, &
      V3_RadiationRecoveryBase
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: &
      TraceCandidateStep, TraceProjectedStep, TraceProjectionCorrection
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
    ALLOCATE( T_Conservative(nX_G), Ye_Conservative(nX_G) )
    ALLOCATE( Ym_Conservative(nX_G) )
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
    ALLOCATE( SchurHasPrevious(nX_G) )
    ALLOCATE( SchurPreviousInput(nX_G), SchurPreviousMap(nX_G) )
    ALLOCATE( SchurPreviousResidual(nX_G) )
    ALLOCATE( SchurGain(nX_G), SchurDamping(nX_G) )
    ALLOCATE &
      ( T_RadiationRecoveryBase(nX_G), &
        Ye_RadiationRecoveryBase(nX_G), &
        Ym_RadiationRecoveryBase(nX_G), &
        E_RadiationRecoveryBase(nX_G), &
        V1_RadiationRecoveryBase(nX_G), &
        V2_RadiationRecoveryBase(nX_G), &
        V3_RadiationRecoveryBase(nX_G) )
    SchurHasPrevious = .FALSE.
    SchurPreviousInput = Zero
    SchurPreviousMap = Zero
    SchurPreviousResidual = Zero
    SchurGain = Zero
    SchurDamping = One
    T_RadiationRecoveryBase = Zero
    Ye_RadiationRecoveryBase = Zero
    Ym_RadiationRecoveryBase = Zero
    E_RadiationRecoveryBase = Zero
    V1_RadiationRecoveryBase = Zero
    V2_RadiationRecoveryBase = Zero
    V3_RadiationRecoveryBase = Zero
    IF( TraceRadiationOnlySolve .AND. UseRadiationOnlySolve )THEN
      ALLOCATE( TraceEvaluationMask(nX_G) )
      ALLOCATE( TraceHasPrevious(nX_G), TraceHasPrevious2(nX_G) )
      ALLOCATE( TraceBacktrackCount(nX_G) )
      ALLOCATE( TracePreviousState(n_FP_inner,nX_G) )
      ALLOCATE( TracePrevious2State(n_FP_inner,nX_G) )
      ALLOCATE( TracePreviousResidual(n_FP_inner,nX_G) )
      ALLOCATE( TraceAlpha(M_inner,nX_G) )
      ALLOCATE( TraceCandidateStep(nX_G), TraceProjectedStep(nX_G) )
      ALLOCATE( TraceProjectionCorrection(nX_G) )
      TraceEvaluationMask = .FALSE.
      TraceHasPrevious = .FALSE.
      TraceHasPrevious2 = .FALSE.
      TraceBacktrackCount = 0
      TracePreviousState = Zero
      TracePrevious2State = Zero
      TracePreviousResidual = Zero
      TraceAlpha = Zero
      TraceCandidateStep = Zero
      TraceProjectedStep = Zero
      TraceProjectionCorrection = Zero
    END IF
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
      LeptonFlavorCheck_Rate_Final(:,:,iN_X) = Zero
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
    !$OMP             T_Conservative, Ye_Conservative, Ym_Conservative, &
    !$OMP             E_Conservative, &
    !$OMP             V_u_1_Conservative, V_u_2_Conservative, &
    !$OMP             V_u_3_Conservative, &
    !$OMP             FVECm_projection, GVECm_projection, &
    !$OMP             InitialRadiationLimiter_ConservationCorrection, &
    !$OMP             RadiationLeptonInvariant, &
    !$OMP             SchurHasPrevious, SchurPreviousInput, &
    !$OMP             SchurPreviousMap, SchurPreviousResidual, &
    !$OMP             SchurGain, SchurDamping, &
    !$OMP             T_RadiationRecoveryBase, &
    !$OMP             Ye_RadiationRecoveryBase, &
    !$OMP             Ym_RadiationRecoveryBase, &
    !$OMP             E_RadiationRecoveryBase, &
    !$OMP             V1_RadiationRecoveryBase, &
    !$OMP             V2_RadiationRecoveryBase, &
    !$OMP             V3_RadiationRecoveryBase )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( ProjectionMask, RadiationBacktrackMask, &
    !$ACC         T_Conservative, Ye_Conservative, Ym_Conservative, &
    !$ACC         E_Conservative, &
    !$ACC         V_u_1_Conservative, V_u_2_Conservative, &
    !$ACC         V_u_3_Conservative, &
    !$ACC         FVECm_projection, GVECm_projection, &
    !$ACC         InitialRadiationLimiter_ConservationCorrection, &
    !$ACC         RadiationLeptonInvariant, &
    !$ACC         SchurHasPrevious, SchurPreviousInput, &
    !$ACC         SchurPreviousMap, SchurPreviousResidual, &
    !$ACC         SchurGain, SchurDamping, &
    !$ACC         T_RadiationRecoveryBase, &
    !$ACC         Ye_RadiationRecoveryBase, &
    !$ACC         Ym_RadiationRecoveryBase, &
    !$ACC         E_RadiationRecoveryBase, &
    !$ACC         V1_RadiationRecoveryBase, &
    !$ACC         V2_RadiationRecoveryBase, &
    !$ACC         V3_RadiationRecoveryBase )
#endif
#endif

#if defined( TWOMOMENT_ORDER_V )
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP ( SchurHasPrevious, SchurPreviousInput, SchurPreviousMap, &
    !$OMP   SchurPreviousResidual, SchurGain, SchurDamping, &
    !$OMP   T_RadiationRecoveryBase, Ye_RadiationRecoveryBase, &
    !$OMP   Ym_RadiationRecoveryBase, E_RadiationRecoveryBase, &
    !$OMP   V1_RadiationRecoveryBase, V2_RadiationRecoveryBase, &
    !$OMP   V3_RadiationRecoveryBase )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( SchurHasPrevious, SchurPreviousInput, SchurPreviousMap, &
    !$ACC   SchurPreviousResidual, SchurGain, SchurDamping, &
    !$ACC   T_RadiationRecoveryBase, Ye_RadiationRecoveryBase, &
    !$ACC   Ym_RadiationRecoveryBase, E_RadiationRecoveryBase, &
    !$ACC   V1_RadiationRecoveryBase, V2_RadiationRecoveryBase, &
    !$ACC   V3_RadiationRecoveryBase )
#endif
#endif


    CALL SetOpacityMask &
         ( nX_G, nDOFX, D, T, Ye, Ym, Opacity_Mask )

    CALL SetMatterSourceMask &
         ( Opacity_Mask, MatterSource_Mask )

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

        ITERATE_outer(iN_X) = MatterSource_Mask(iEf,iN_X)

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

        ITERATE_outer(iN_X) = MatterSource_Mask(iEf,iN_X)

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
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Ye, Ym, E, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

    CALL ApplyNeutrinoDistributionLimit_OrderV &
           ( dt, ITERATE_outer, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Nnu, &
             V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             InitialRadiationLimiter_Rate_Final, &
             InitialRadiationLimiter_ConservationCorrection, &
             LeptonFlavorCheck_Rate_Final(:,iLFC_InputConversion,:), &
             LeptonFlavorCheck_Rate_Final(:,iLFC_EntryProjection,:) )

    CALL InitializeRHS_OrderV &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, &
             Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
             InitialRadiationLimiter_ConservationCorrection )

    IF( UseRadiationOnlySolve ) &
      CALL InitializeRadiationLeptonInvariants_OrderV &
             ( ITERATE_outer, C_Dnu, RadiationLeptonInvariant )

#elif defined( TWOMOMENT_RELATIVISTIC )

    CALL ComputePressure_TABLE & 
           ( D, T, Ye, Ym, P )

    CALL InitializeRHS_Relativistic &
           ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Ye, Ym, E, P, V_u_1, V_u_2, V_u_3, &
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
           ( D, T, Ye, Ym, SqrtGm, ITERATE_outer, nX_P_outer, &
             PackIndex_outer, UnpackIndex_outer, &
             T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, &
             Opacity_Mask_Input = Opacity_Mask )

    CALL TimersStop( Timer_Collisions_ComputeOpacity )

#if defined( TWOMOMENT_ORDER_V )
    IF( DebugMuonCoupling ) &
      CALL CheckMuonCouplingIteration_OrderV &
             ( 'initial opacity construction', dt, 0, 0, ITERATE_outer, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               CheckOpacity_Option = .TRUE. )
#endif

#if defined( TWOMOMENT_ORDER_V )
    END IF
#endif

#if defined( TWOMOMENT_ORDER_V )
    IF( UseRadiationOnlySolve )THEN

      IF( TraceRadiationOnlySolve .OR. AnalyzeRadiationOnlyJacobian ) &
        RadiationOnlySolveCallCount = RadiationOnlySolveCallCount + 1

      IF( TraceRadiationOnlySolve )THEN
        CALL InitializeRadiationOnlyTrace_OrderV
      END IF

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

      IF( ReducedSolverMethod == iReducedSolver_CoupledEtaSchur )THEN
        CALL SolveRadiationMuonEtaSchur_OrderV &
               ( dt, ITERATE_inner, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                 FVECm_outer, GVECm_outer, FVECm_inner, GVECm_inner, &
                 InitialRadiationLimiter_ConservationCorrection, &
                 RadiationLeptonInvariant, nIterations_Inner, &
                 nIterations_Outer, nIterations_Inner_Final, &
                 nIterations_Outer_Final, RelativeResidual_Inner_First, &
                 RadiationResidual_Relative_Final )
      END IF

      DO WHILE( ANY( ITERATE_inner(:) ) .AND. k_inner < MaxIter_inner )

        k_inner = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_inner, nX_P_inner, &
                 PackIndex_inner, UnpackIndex_inner )

        CALL TimersStart( Timer_Collisions_MatterFromRadiation )

        CALL SaveRadiationRecoveryMatterState_OrderV &
               ( ITERATE_inner, T, Ye, Ym, E, &
                 V_u_1, V_u_2, V_u_3, &
                 T_RadiationRecoveryBase, Ye_RadiationRecoveryBase, &
                 Ym_RadiationRecoveryBase, E_RadiationRecoveryBase, &
                 V1_RadiationRecoveryBase, V2_RadiationRecoveryBase, &
                 V3_RadiationRecoveryBase )

        CALL RecoverMatterFromRadiationInvariants_OrderV &
               ( ITERATE_inner, ProjectionMask, &
                 FVECm_outer, GVECm_outer, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
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

            IF( TraceRadiationOnlySolve )THEN
              DO iN_X = 1, nX_G
                IF( RadiationBacktrackMask(iN_X) ) &
                  TraceBacktrackCount(iN_X) = TraceBacktrackCount(iN_X) + 1
              END DO
            END IF

            CALL RestoreRadiationRecoveryMatterState_OrderV &
                   ( RadiationBacktrackMask, &
                     T_RadiationRecoveryBase, Ye_RadiationRecoveryBase, &
                     Ym_RadiationRecoveryBase, E_RadiationRecoveryBase, &
                     V1_RadiationRecoveryBase, V2_RadiationRecoveryBase, &
                     V3_RadiationRecoveryBase, &
                     T, Ye, Ym, E, V_u_1, V_u_2, V_u_3 )

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
                     D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
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

        ! An inadmissible Anderson iterate is recoverable here: the
        ! backtracking loop above first retries progressively shorter
        ! updates.  Report it as a fatal matter-recovery failure only if
        ! that recovery mechanism has been exhausted.
        IF( ( DebugMuonCoupling .OR. TraceRadiationOnlySolve ) &
            .AND. ANY( ITERATE_inner .AND. Error /= 0 ) ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'matter recovery after backtracking', &
                   dt, 0, k_inner, ITERATE_inner, &
                   D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   Error_Option = Error )

        CALL CheckError_FP &
               ( ITERATE_outer, ITERATE_inner, 0, k_inner, Error, &
                 D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 IgnoreIterationLimit_Option = .TRUE. )

        IF( .NOT. FreezeOpacities )THEN

          CALL TimersStart( Timer_Collisions_ComputeOpacity )

          CALL ComputeOpacities_Packed &
                 ( D, T, Ye, Ym, SqrtGm, &
                   ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner, &
                   T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, &
                   Opacity_Mask_Input = Opacity_Mask )

          CALL TimersStop( Timer_Collisions_ComputeOpacity )

          IF( DebugMuonCoupling ) &
            CALL CheckMuonCouplingIteration_OrderV &
                   ( 'opacity construction', dt, 0, k_inner, ITERATE_inner, &
                     D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     CheckOpacity_Option = .TRUE. )

        END IF

        CALL ComputeDnuNorm( ITERATE_inner, Dnu )

        CALL TimersStart( Timer_Collisions_ComputeRates )

        IF( FreezeOpacities )THEN

          CALL ComputeRates_Packed &
                 ( D, T_Opacity, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner, &
                   nX_P_outer, Opacity_Mask )

        ELSE

          CALL ComputeRates_Packed &
                 ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   ITERATE_inner, nX_P_inner, &
                   PackIndex_inner, UnpackIndex_inner, &
                   nX_P_inner, Opacity_Mask )

        END IF

        CALL TimersStop( Timer_Collisions_ComputeRates )

        IF( DebugMuonCoupling ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'rate construction', dt, 0, k_inner, ITERATE_inner, &
                   D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE. )

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

        IF( DebugMuonCoupling ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'fixed-point map', dt, 0, k_inner, ITERATE_inner, &
                   D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE., &
                   Fm_Option = FVECm_inner, Gm_Option = GVECm_inner )

        IF( ReducedSolverMethod == iReducedSolver_SchurNewton ) &
          CALL PrepareRadiationOnlySchurNewton_OrderV &
                 ( ITERATE_inner, GVECm_inner, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   SchurHasPrevious, SchurPreviousInput, &
                   SchurPreviousMap, SchurPreviousResidual, &
                   SchurGain, SchurDamping )

        CALL TimersStop( Timer_Collisions_NeutrinoRHS )

        CALL TimersStart( Timer_Collisions_CheckInner )

        IF( TraceRadiationOnlySolve ) &
          TraceEvaluationMask = ITERATE_inner

        CALL CheckConvergence_RadiationOnly &
               ( ITERATE_inner, k_inner, &
                 nIterations_Inner, nIterations_Outer, &
                 nIterations_Inner_Final, nIterations_Outer_Final, &
                 RelativeResidual_Inner_First, &
                 RadiationResidual_Relative_Final, FVECm_inner )

        CALL TimersStop( Timer_Collisions_CheckInner )

        IF( TraceRadiationOnlySolve ) &
          CALL WriteRadiationOnlyIterationTrace_OrderV &
                 ( dt, k_inner, TraceEvaluationMask, ITERATE_inner, Error, &
                   D, T, Ye, Ym, E, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   FVECm_inner, GVECm_inner, RadiationLeptonInvariant, &
                   InitialRadiationLimiter_ConservationCorrection, &
                   TraceHasPrevious, TraceHasPrevious2, &
                   TracePreviousState, TracePrevious2State, &
                   TracePreviousResidual, TraceBacktrackCount, &
                   TraceCandidateStep, TraceProjectedStep, &
                   TraceProjectionCorrection, TraceAlpha, &
                   SchurGain, SchurDamping )

        ! Capture the iteration-limit state before the Anderson update.  At
        ! this point the matter state, opacities, rates, fixed-point map, and
        ! residual all belong to the same evaluated radiation iterate.
        IF( ( DebugMuonCoupling .OR. TraceRadiationOnlySolve ) &
            .AND. ANY( ITERATE_inner(:) ) &
            .AND. k_inner >= MaxIter_inner ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'radiation-solve iteration limit', dt, 0, k_inner, &
                   ITERATE_inner, D, T, Ye, Ym, E, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE., &
                   Fm_Option = FVECm_inner, Gm_Option = GVECm_inner, &
                   ForceFailure_Option = .TRUE. )

        IF( ANY( ITERATE_inner(:) ) )THEN

          IF( ReducedSolverMethod == iReducedSolver_SchurNewton ) &
            CALL ApplyRadiationOnlySchurNewtonPreconditioner_OrderV &
                   ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                     SchurDamping )

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

          IF( TraceRadiationOnlySolve ) &
            CALL CaptureRadiationOnlyUpdateTrace_OrderV &
                   ( ITERATE_inner, Mk_inner, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, GVECm_inner, &
                     Alpha_inner, TracePreviousState, &
                     TraceCandidateStep, TraceProjectedStep, &
                     TraceProjectionCorrection, TraceAlpha )

          IF( DebugMuonCoupling ) &
            CALL CheckMuonCouplingIteration_OrderV &
                   ( 'Anderson-updated iterate', dt, 0, k_inner, ITERATE_inner, &
                     D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     CheckOpacity_Option = .TRUE., &
                     CheckRates_Option = .TRUE., &
                     Fm_Option = FVECm_inner, Gm_Option = GVECm_inner )

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

      IF( AnalyzeRadiationOnlyJacobian &
            .AND. ReducedSolverMethod /= iReducedSolver_CoupledEtaSchur &
            .AND. RadiationOnlySolveCallCount == 1 ) &
        CALL AnalyzeRadiationOnlyJacobian_OrderV &
               ( dt, nIterations_Outer_Final, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                 FVECm_outer, GVECm_outer, &
                 FVECm_inner, GVECm_inner, &
                 ITERATE_inner, RadiationBacktrackMask, Error, &
                 PackIndex_inner, UnpackIndex_inner, &
                 InitialRadiationLimiter_ConservationCorrection, &
                 RadiationLeptonInvariant )

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
                 ( D, T, Ye, Ym, SqrtGm, ITERATE_outer, nX_P_outer, &
                   PackIndex_outer, UnpackIndex_outer, &
                   T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, &
                   Opacity_Mask_Input = Opacity_Mask )

          CALL TimersStop( Timer_Collisions_ComputeOpacity )

          IF( DebugMuonCoupling ) &
            CALL CheckMuonCouplingIteration_OrderV &
                   ( 'nested opacity construction', dt, k_outer, 0, &
                     ITERATE_outer, D, T, Ye, Ym, E, &
                     V_u_1, V_u_2, V_u_3, &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                     Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                     CheckOpacity_Option = .TRUE. )

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

        IF( FreezeOpacities )THEN
          CALL ComputeRates_Packed &
                 ( D, T_Opacity, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   ITERATE_inner, nX_P_inner, PackIndex_inner, &
                   UnpackIndex_inner, nX_P_outer, Opacity_Mask )
        ELSE
          CALL ComputeRates_Packed &
                 ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   ITERATE_inner, nX_P_inner, PackIndex_inner, &
                   UnpackIndex_inner, nX_P_outer, Opacity_Mask )
        END IF

        CALL TimersStop( Timer_Collisions_ComputeRates )

#if defined( TWOMOMENT_ORDER_V )
        IF( DebugMuonCoupling ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'nested rate construction', dt, k_outer, k_inner, &
                   ITERATE_inner, D, T, Ye, Ym, E, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE. )
#endif
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

#if defined( TWOMOMENT_ORDER_V )
        IF( DebugMuonCoupling ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'nested fixed-point map', dt, k_outer, k_inner, &
                   ITERATE_inner, D, T, Ye, Ym, E, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE., &
                   Fm_Option = FVECm_inner, Gm_Option = GVECm_inner )
#endif

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

#if defined( TWOMOMENT_ORDER_V )
        IF( DebugMuonCoupling .AND. k_inner >= MaxIter_inner &
              .AND. ANY(ITERATE_inner) ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'nested inner iteration limit', dt, k_outer, k_inner, &
                   ITERATE_inner, D, T, Ye, Ym, E, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE., &
                   Fm_Option = FVECm_inner, Gm_Option = GVECm_inner, &
                   ForceFailure_Option = .TRUE. )
#endif

        CALL CheckError_Inner &
               ( ITERATE_inner, k_outer, k_inner, &
                 D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        CALL TimersStop( Timer_Collisions_CheckInner )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Collisions_UpdateFP )

#if defined( TWOMOMENT_ORDER_V )
        IF( DebugMuonCoupling ) &
          CALL CheckMuonCouplingIteration_OrderV &
                 ( 'nested Anderson-updated iterate', dt, k_outer, k_inner, &
                   ITERATE_inner, D, T, Ye, Ym, E, &
                   V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   CheckOpacity_Option = .TRUE., &
                   CheckRates_Option = .TRUE., &
                   Fm_Option = FVECm_inner, Gm_Option = GVECm_inner )
#endif

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
                 D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                 ProjectionMask, ITERATE_outer, &
                 nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
                 FVECm_inner, GVECm_inner, &
                 FVECm_projection, GVECm_projection, &
                 T_Conservative, Ye_Conservative, Ym_Conservative, &
                 E_Conservative, &
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
                 D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        IF(      ANY( ITERATE_outer ) &
            .AND. k_terminal >= MaxIter_outer )THEN

          CALL CheckError_FP &
                 ( ITERATE_outer, ITERATE_inner, &
                   MaxIter_outer, k_inner, Error, &
                   D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
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

      ! The inner map's cached transfer rates can precede its last radiation
      ! update. Refresh only the general-lepton rates for the accepted state;
      ! no opacity/kernel evaluation and no NNS/NES/pair rate pass is needed.
      IF( UseMuonLeptonSource .AND. (Include_LFE .OR. Include_LFC .OR. Include_IMD) )THEN
        CALL TimersStart(Timer_OpacityRate_GeneralLepton)
        IF( FreezeOpacities )THEN
          CALL ComputeNeutrinoLeptonRates &
            (T_Opacity,W2_N,Dnu,Opacity_Mask(iOp_NMS,:),Opacity_Mask(iOp_LFE,:), &
             Opacity_Mask(iOp_LFC,:),Opacity_Mask(iOp_IMD,:), &
             Include_NMS,Include_LFE,Include_LFC,Include_IMD,nX_P_outer,UnpackIndex_outer)
        ELSE
          CALL ComputeNeutrinoLeptonRates &
            (T,W2_N,Dnu,Opacity_Mask(iOp_NMS,:),Opacity_Mask(iOp_LFE,:), &
             Opacity_Mask(iOp_LFC,:),Opacity_Mask(iOp_IMD,:), &
             Include_NMS,Include_LFE,Include_LFC,Include_IMD,nX_P_outer,UnpackIndex_outer)
        END IF
        CALL TimersStop(Timer_OpacityRate_GeneralLepton)
      END IF

      CALL ComputeMatterRHS_OrderV &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               dt, InitialRadiationLimiter_ConservationCorrection(iYm,:) )

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
               Ye, Ym, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_ORDER_V )

      CALL UpdateMatterRHS_OrderV &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#elif defined( TWOMOMENT_RELATIVISTIC )

      CALL UpdateMatterRHS_Relativistic &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               D, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#endif

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Ye, Ym, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer, Error )

#if defined( TWOMOMENT_ORDER_V )
      IF( DebugMuonCoupling ) &
        CALL CheckMuonCouplingIteration_OrderV &
               ( 'nested matter update', dt, k_outer, k_inner, &
                 ITERATE_outer, D, T, Ye, Ym, E, &
                 V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Error_Option = Error, &
                 CheckOpacity_Option = .TRUE., &
                 CheckRates_Option = .TRUE. )
#endif

#if   defined( TWOMOMENT_RELATIVISTIC )

      CALL ComputePressure_TABLE( D, T, Ye, Ym, P )

#endif

      ! --- Check Convergence (outer) ---

      CALL TimersStart( Timer_Collisions_CheckOuter )

      CALL CheckConvergence_Outer &
             ( ITERATE_outer, ITERATE_inner, n_FP_outer, k_outer, &
               nIterations_Outer, FVECm_outer )

      ! CALL PrintStatus_FP &
      !       ( ITERATE_outer, ITERATE_inner, k_outer, k_inner, Error, &
      !         FVECm_outer, FVECm_inner, GVECm_outer, GVECm_inner, &
      !         D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
      !         Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      !         Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      !         Gm_dd_11, Gm_dd_22, Gm_dd_33, 1 )

#if defined( TWOMOMENT_ORDER_V )
      IF( DebugMuonCoupling .AND. k_outer >= MaxIter_outer &
            .AND. ANY(ITERATE_outer) ) &
        CALL CheckMuonCouplingIteration_OrderV &
               ( 'nested outer iteration limit', dt, k_outer, k_inner, &
                 ITERATE_outer, D, T, Ye, Ym, E, &
                 V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 Error_Option = Error, &
                 CheckOpacity_Option = .TRUE., &
                 CheckRates_Option = .TRUE., &
                 ForceFailure_Option = .TRUE. )
#endif

      CALL CheckError_FP &
             ( ITERATE_outer, ITERATE_inner, k_outer, k_inner, Error, &
               D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
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
                 D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                 ProjectionMask, ITERATE_outer, &
                 nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
                 FVECm_inner, GVECm_inner, &
                 FVECm_projection, GVECm_projection, &
                 T_Conservative, Ye_Conservative, Ym_Conservative, &
                 E_Conservative, &
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
                 D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
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

      CALL CreatePackIndex &
             ( ProjectionMask, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

      CALL ComputeOpacities_Packed &
             ( D, T, Ye, Ym, SqrtGm, ProjectionMask, nX_P_inner, &
               PackIndex_inner, UnpackIndex_inner, &
               T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, Opacity_Mask_Input = Opacity_Mask )

      CALL ComputeRates_Packed &
             ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               ProjectionMask, nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
               nX_P_inner, Opacity_Mask )

      CALL ComputeMatterRHS_OrderV &
             ( ProjectionMask, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               dt, InitialRadiationLimiter_ConservationCorrection(iYm,:) )

      CALL ComputeNeutrinoRHS_OrderV &
             ( ProjectionMask, FVECm_inner, GVECm_inner, dt, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL ComputeLeptonFlavorCheck_OrderV &
             ( ProjectionMask, dt, &
               UseRadiationOnlySolve, &
               FVECm_inner, GVECm_inner, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               LeptonFlavorCheck_Rate_Final )

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

        MatterEquationResidual_Rate_Final(iQ_Le,iN_X) &
          = D_old(iN_X) * Ye_old(iN_X) / AtomicMassUnit &
              * FVECm_outer(iYe,iN_X) / dt

        MatterEquationResidual_Rate_Final(iQ_Lm,iN_X) &
          = D_old(iN_X) * Ym_old(iN_X) / AtomicMassUnit &
              * FVECm_outer(iYm,iN_X) / dt

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
    !$OMP               T_Conservative, Ye_Conservative, Ym_Conservative, &
    !$OMP               E_Conservative, &
    !$OMP               V_u_1_Conservative, V_u_2_Conservative, &
    !$OMP               V_u_3_Conservative, &
    !$OMP               FVECm_projection, GVECm_projection, &
    !$OMP               InitialRadiationLimiter_ConservationCorrection, &
    !$OMP               RadiationLeptonInvariant, &
    !$OMP               SchurHasPrevious, SchurPreviousInput, &
    !$OMP               SchurPreviousMap, SchurPreviousResidual, &
    !$OMP               SchurGain, SchurDamping, &
    !$OMP               T_RadiationRecoveryBase, &
    !$OMP               Ye_RadiationRecoveryBase, &
    !$OMP               Ym_RadiationRecoveryBase, &
    !$OMP               E_RadiationRecoveryBase, &
    !$OMP               V1_RadiationRecoveryBase, &
    !$OMP               V2_RadiationRecoveryBase, &
    !$OMP               V3_RadiationRecoveryBase )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ProjectionMask, RadiationBacktrackMask, &
    !$ACC         T_Conservative, Ye_Conservative, Ym_Conservative, &
    !$ACC         E_Conservative, &
    !$ACC         V_u_1_Conservative, V_u_2_Conservative, &
    !$ACC         V_u_3_Conservative, &
    !$ACC         FVECm_projection, GVECm_projection, &
    !$ACC         InitialRadiationLimiter_ConservationCorrection, &
    !$ACC         RadiationLeptonInvariant, &
    !$ACC         SchurHasPrevious, SchurPreviousInput, &
    !$ACC         SchurPreviousMap, SchurPreviousResidual, &
    !$ACC         SchurGain, SchurDamping, &
    !$ACC         T_RadiationRecoveryBase, &
    !$ACC         Ye_RadiationRecoveryBase, &
    !$ACC         Ym_RadiationRecoveryBase, &
    !$ACC         E_RadiationRecoveryBase, &
    !$ACC         V1_RadiationRecoveryBase, &
    !$ACC         V2_RadiationRecoveryBase, &
    !$ACC         V3_RadiationRecoveryBase )
#endif
#endif

    IF( ALLOCATED(TraceEvaluationMask) )THEN
      DEALLOCATE( TraceEvaluationMask, TraceHasPrevious, TraceHasPrevious2 )
      DEALLOCATE( TraceBacktrackCount )
      DEALLOCATE( TracePreviousState, TracePrevious2State )
      DEALLOCATE( TracePreviousResidual, TraceAlpha )
      DEALLOCATE( TraceCandidateStep, TraceProjectedStep )
      DEALLOCATE( TraceProjectionCorrection )
    END IF

    DEALLOCATE( SchurHasPrevious )
    DEALLOCATE( SchurPreviousInput, SchurPreviousMap )
    DEALLOCATE( SchurPreviousResidual, SchurGain, SchurDamping )
    DEALLOCATE &
      ( T_RadiationRecoveryBase, Ye_RadiationRecoveryBase, &
        Ym_RadiationRecoveryBase, E_RadiationRecoveryBase, &
        V1_RadiationRecoveryBase, V2_RadiationRecoveryBase, &
        V3_RadiationRecoveryBase )

  END SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA


  SUBROUTINE SetMatterSourceMask &
    ( OpacityMask, MatterSourceMask )

    LOGICAL,  INTENT(in)  :: OpacityMask(:,:)
    LOGICAL,  INTENT(out) :: MatterSourceMask(:,:)

    INTEGER :: iN_X, iNodeX, iX
    LOGICAL :: CollisionActive
    LOGICAL :: ElectronEmAbActive, ElectronLeptonActive
    LOGICAL :: MuonEmAbActive, FlavorExchangeActive
    LOGICAL :: UseNNS, UseNES, UsePair, UseNuPair, UseBrem
    LOGICAL :: UseNMS, UseLFE, UseLFC, UseIMD

    UseNNS    = Include_NNS
    UseNES    = Include_NES
    UsePair   = Include_Pair
    UseNuPair = Include_NuPair
    UseBrem   = Include_Brem
    UseNMS    = Include_NMS
    UseLFE    = Include_LFE
    UseLFC    = Include_LFC
    UseIMD    = Include_IMD

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP FIRSTPRIVATE( UseNNS, UseNES, UsePair, UseNuPair, UseBrem, &
    !$OMP               UseNMS, UseLFE, UseLFC, UseIMD ) &
    !$OMP PRIVATE( iN_X, iNodeX, CollisionActive, ElectronEmAbActive, &
    !$OMP          ElectronLeptonActive, MuonEmAbActive, FlavorExchangeActive )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC FIRSTPRIVATE( UseNNS, UseNES, UsePair, UseNuPair, UseBrem, &
    !$ACC               UseNMS, UseLFE, UseLFC, UseIMD ) &
    !$ACC PRIVATE( iN_X, iNodeX, CollisionActive, ElectronEmAbActive, &
    !$ACC          ElectronLeptonActive, MuonEmAbActive, FlavorExchangeActive )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP FIRSTPRIVATE( UseNNS, UseNES, UsePair, UseNuPair, UseBrem, &
    !$OMP               UseNMS, UseLFE, UseLFC, UseIMD ) &
    !$OMP PRIVATE( iN_X, iNodeX, CollisionActive, ElectronEmAbActive, &
    !$OMP          ElectronLeptonActive, MuonEmAbActive, FlavorExchangeActive )
#endif
    DO iX = 1, nX_G / nDOFX

      iN_X = 1 + ( iX - 1 ) * nDOFX

      ElectronEmAbActive = nSpecies >= iNuE_Bar .AND. &
        OpacityMask(iOp_EmAb_Electron,iN_X)
      MuonEmAbActive = nSpecies >= iNuM_Bar .AND. &
        OpacityMask(iOp_EmAb_Muon,iN_X)

      ElectronLeptonActive = nSpecies >= iNuE_Bar .AND. &
        ( ElectronEmAbActive .OR. OpacityMask(iOp_ECTable,iN_X) )

      FlavorExchangeActive = nSpecies >= iNuM_Bar .AND. &
        ( ( UseLFE .AND. OpacityMask(iOp_LFE,iN_X) ) .OR. &
          ( UseLFC .AND. OpacityMask(iOp_LFC,iN_X) ) .OR. &
          ( UseIMD .AND. OpacityMask(iOp_IMD,iN_X) ) )

      ElectronLeptonActive = ElectronLeptonActive .OR. FlavorExchangeActive
      MuonEmAbActive = MuonEmAbActive .OR. FlavorExchangeActive

      CollisionActive = &
          ElectronLeptonActive .OR. &
          MuonEmAbActive .OR. &
          OpacityMask(iOp_Iso,    iN_X) .OR. &
          ( UseNNS    .AND. OpacityMask(iOp_NNS,   iN_X) ) .OR. &
          ( UseNES    .AND. OpacityMask(iOp_NES,   iN_X) ) .OR. &
          ( UsePair   .AND. OpacityMask(iOp_Pair,  iN_X) ) .OR. &
          ( UseNuPair .AND. OpacityMask(iOp_NuPair,iN_X) ) .OR. &
          ( UseBrem   .AND. OpacityMask(iOp_Brem,  iN_X) ) .OR. &
          ( UseNMS    .AND. OpacityMask(iOp_NMS,   iN_X) ) .OR. &
          FlavorExchangeActive

      DO iNodeX = 1, nDOFX
        iN_X = iNodeX + ( iX - 1 ) * nDOFX

#if defined( TWOMOMENT_ORDER_1 )
        MatterSourceMask(iYe,iN_X) = ElectronLeptonActive
        MatterSourceMask(iYm,iN_X) = MuonEmAbActive
        MatterSourceMask(iEf,iN_X) = CollisionActive
#elif defined( TWOMOMENT_ORDER_V )
        MatterSourceMask(iYe,iN_X) = ElectronLeptonActive
        MatterSourceMask(iYm,iN_X) = MuonEmAbActive
        MatterSourceMask(iEf,iN_X) = CollisionActive
        MatterSourceMask(iV1,iN_X) = CollisionActive
        MatterSourceMask(iV2,iN_X) = CollisionActive
        MatterSourceMask(iV3,iN_X) = CollisionActive
#elif defined( TWOMOMENT_RELATIVISTIC )
        MatterSourceMask(iD ,iN_X) = CollisionActive
        MatterSourceMask(iYe,iN_X) = ElectronLeptonActive
        MatterSourceMask(iYm,iN_X) = MuonEmAbActive
        MatterSourceMask(iEf,iN_X) = CollisionActive
        MatterSourceMask(iV1,iN_X) = CollisionActive
        MatterSourceMask(iV2,iN_X) = CollisionActive
        MatterSourceMask(iV3,iN_X) = CollisionActive
#endif
      END DO

    END DO

  END SUBROUTINE SetMatterSourceMask


  REAL(DP) FUNCTION MatterSourceWeight( iEquation, iN_X )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iEquation, iN_X

    IF( MatterSource_Mask(iEquation,iN_X) )THEN
      MatterSourceWeight = wMatrRHS(iEquation)
    ELSE
      MatterSourceWeight = Zero
    END IF

  END FUNCTION MatterSourceWeight


  LOGICAL FUNCTION LeptonFamilyCouplesToMatter( iFamily, iN_X )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iFamily, iN_X

    SELECT CASE( iFamily )
      CASE( 1 )
        LeptonFamilyCouplesToMatter = MatterSource_Mask(iYe,iN_X)
      CASE( 2 )
        LeptonFamilyCouplesToMatter = MatterSource_Mask(iYm,iN_X)
      CASE DEFAULT
        LeptonFamilyCouplesToMatter = .FALSE.
    END SELECT

    ! A zero material weight is a prescribed-reservoir experiment: an
    ! active charged/flavor reaction can still change radiation number.
    ! Only absence of the physical reaction permits a radiation-only
    ! signed-number constraint. The eta unknown is separately disabled
    ! when its material weight is zero.

  END FUNCTION LeptonFamilyCouplesToMatter


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

      DO iFamily = 1, MIN( nLeptonFamilies, ( nSpecies + 1 ) / 2 )

        IF( LeptonFamilyCouplesToMatter( iFamily, iN_X ) ) CYCLE

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
    DO iFamily = 1, MIN( nLeptonFamilies, ( nSpecies + 1 ) / 2 )

      IF( .NOT. MASK(iN_X) ) CYCLE
      IF( LeptonFamilyCouplesToMatter( iFamily, iN_X ) ) CYCLE

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


  SUBROUTINE SolveRadiationMuonEtaSchur_OrderV &
    ( dt, Active, Dnu, I1, I2, I3, D, T, Ye, Ym, E, V1, V2, V3, &
      g1, g2, g3, SqrtGm, Fmatter, Gmatter, Frad, Grad, &
      LimiterCorrection, LeptonInvariant, nInner, nOuter, &
      nInnerFinal, nOuterFinal, FirstResidual, FinalResidual )

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE mpi

    REAL(DP), INTENT(in) :: dt
    LOGICAL, INTENT(inout) :: Active(:)
    REAL(DP), INTENT(inout) :: Dnu(:,:,:), I1(:,:,:), I2(:,:,:), I3(:,:,:)
    REAL(DP), INTENT(in) :: D(:), g1(:), g2(:), g3(:), SqrtGm(:)
    REAL(DP), INTENT(inout) :: T(:), Ye(:), Ym(:), E(:), V1(:), V2(:), V3(:)
    REAL(DP), INTENT(inout) :: Fmatter(:,:), Gmatter(:,:), Frad(:,:), Grad(:,:)
    REAL(DP), INTENT(in) :: LimiterCorrection(:,:), LeptonInvariant(:,:)
    INTEGER, INTENT(inout) :: nInner(:), nOuter(:), nInnerFinal(:), nOuterFinal(:)
    REAL(DP), INTENT(inout) :: FirstResidual(:), FinalResidual(:)

    ! The radiation residual is x-P G(x,eta), with the existing independent
    ! spectral bound projection. In inventory mode the muon residual is
    ! Y_EOS-Y_conserved and the Schur step uses the equivalent row H=q-L F.
    ! In direct-source mode both the accepted row and the Newton row are
    ! Y_EOS-Y_old-dt*S_mu, plus the existing initial-limiter compensation.
    ! Number-conserving reaction residuals never enter that material row.
    ! No inverse muon susceptibility appears in this residual evaluation.
    INTEGER, PARAMETER :: KrylovDimension = 24
    REAL(DP), PARAMETER :: FDStep = 1.0e-6_DP
    REAL(DP), PARAMETER :: Armijo = 1.0e-4_DP
    INTEGER :: k, p, s, eN, os, bt, ierr, rank, unit, io, nInitialPacked
    REAL(DP) :: mu, densityScale, leptonScale, denom, dqz, dqy
    REAL(DP), ALLOCATABLE :: X(:,:), Eta(:), F(:,:), Q(:), Scale(:,:), QScale(:)
    REAL(DP), ALLOCATABLE :: B(:,:), EtaJac(:), Z(:,:), Y(:,:), StepX(:,:), StepEta(:)
    REAL(DP), ALLOCATABLE :: TrialX(:,:), TrialEta(:), TrialF(:,:), TrialQ(:)
    REAL(DP), ALLOCATABLE :: EtaPerturb(:), QZ(:), QY(:), MuonConserved(:)
    REAL(DP), ALLOCATABLE :: Merit(:), TrialMerit(:), Lambda(:), QAcceptance(:)
    REAL(DP), ALLOCATABLE :: WarmT(:), Fbase(:,:), Qbase(:)
    REAL(DP), ALLOCATABLE :: Hvalue(:), Hbase(:), Hsource(:), RawFm(:,:)
    REAL(DP), ALLOCATABLE :: NumberSource(:,:,:), NumberFactor(:,:,:), NumberFraction(:,:,:)
    REAL(DP), ALLOCATABLE :: MaterialSource(:), MaterialAbsoluteSource(:)
    LOGICAL, ALLOCATABLE :: EtaActive(:), EvalMask(:), Pending(:), Bad(:), WorkMask(:)
    INTEGER, ALLOCATABLE :: Err(:), Pack(:), Unpack(:), InvariantPivot(:,:), EvalReason(:)
    CHARACTER(512) :: TraceName

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM(Dnu,I1,I2,I3,D,T,Ye,Ym,E,V1,V2,V3, &
    !$OMP                    C_Dnu,Ym_old,W2_S,MatterSource_Mask,T_Opacity, &
    !$OMP                    LimiterCorrection,LeptonInvariant,g1,g2,g3)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST(Dnu,I1,I2,I3,D,T,Ye,Ym,E,V1,V2,V3, &
    !$ACC             C_Dnu,Ym_old,W2_S,MatterSource_Mask,T_Opacity, &
    !$ACC             LimiterCorrection,LeptonInvariant,g1,g2,g3)
#endif
    ALLOCATE(X(n_FP_inner,nX_G),Eta(nX_G),F(n_FP_inner,nX_G),Q(nX_G))
    ALLOCATE(Scale(n_FP_inner,nX_G),QScale(nX_G),QAcceptance(nX_G))
    ALLOCATE(B(n_FP_inner,nX_G),EtaJac(nX_G),Z(n_FP_inner,nX_G),Y(n_FP_inner,nX_G))
    ALLOCATE(StepX(n_FP_inner,nX_G),StepEta(nX_G),TrialX(n_FP_inner,nX_G))
    ALLOCATE(TrialEta(nX_G),TrialF(n_FP_inner,nX_G),TrialQ(nX_G))
    ALLOCATE(EtaPerturb(nX_G),QZ(nX_G),QY(nX_G),MuonConserved(nX_G))
    ALLOCATE(Merit(nX_G),TrialMerit(nX_G),Lambda(nX_G),WarmT(nX_G))
    ALLOCATE(Fbase(n_FP_inner,nX_G),Qbase(nX_G))
    ALLOCATE(Hvalue(nX_G),Hbase(nX_G),Hsource(nX_G),RawFm(n_FP_inner,nX_G))
    ALLOCATE(NumberSource(nE_G,nSpecies,nX_G),NumberFactor(nE_G,nSpecies,nX_G))
    ALLOCATE(NumberFraction(nE_G,nSpecies,nX_G))
    IF( UseMuonLeptonSource ) ALLOCATE(MaterialSource(nX_G),MaterialAbsoluteSource(nX_G))
    ALLOCATE(EtaActive(nX_G),EvalMask(nX_G),Pending(nX_G),Bad(nX_G),WorkMask(nX_G))
    ALLOCATE(Err(nX_G),Pack(nX_G),Unpack(nX_G))
    ALLOCATE(EvalReason(nX_G))
    EvalReason = 0
    ALLOCATE(InvariantPivot(nLeptonFamilies,nX_G))
    EtaActive = MatterSource_Mask(iYm,:) .AND. Active .AND. wMatrRHS(iYm) /= Zero
    IF( nSpecies < iNuM_Bar ) EtaActive = .FALSE.
    Eta = Zero
    nInitialPacked = COUNT(Active)
    WarmT = T
    Scale = One
    QScale = One
    Lambda = One
    StepEta = Zero
    StepX = Zero
    Merit = Zero
    TrialMerit = Zero
    Err = 0
    F = Zero
    Q = Zero
    DO p = 1, nX_G
      IF( EtaActive(p) )THEN
        CALL ComputeMuonChemicalPotential_TABLE(D(p),T(p),Ye(p),Ym(p),mu)
        Eta(p) = mu / (BoltzmannConstant*T(p))
      END IF
      DO s = 1, nSpecies
      DO eN = 1, nE_G
        os = ((eN-1)+(s-1)*nE_G)*nCR
        X(os+iCR_N,p) = Dnu(eN,s,p)
        X(os+iCR_G1,p) = g1(p)*I1(eN,s,p)
        X(os+iCR_G2,p) = g2(p)*I2(eN,s,p)
        X(os+iCR_G3,p) = g3(p)*I3(eN,s,p)
      END DO
      END DO
    END DO
    ! Replace one redundant density equation in each uncoupled family
    ! with its signed-number invariant. Choose a populated pivot once per
    ! solve; do not project states during Jacobian or residual evaluation.
    InvariantPivot = 0
    DO p = 1, nX_G
      DO s = 1, MIN(nLeptonFamilies,nSpecies/2)
        IF( LeptonFamilyCouplesToMatter(s,p) ) CYCLE
        eN = MAXLOC(W2_S*MAX(Dnu(:,2*s-1,p),Dnu(:,2*s,p)),DIM=1)
        InvariantPivot(s,p) = ((eN-1)+(2*s-2)*nE_G)*nCR+iCR_N
      END DO
    END DO

    ! Work masks and indices are needed by the existing packed GPU physics
    ! routines. Newton/Krylov control and small Hessenberg solves are on host.
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(alloc: EvalMask,WorkMask,Err,Pack,Unpack,MuonConserved,EtaActive)
    !$OMP TARGET ENTER DATA MAP(alloc: NumberSource,NumberFactor,NumberFraction)
    IF( UseMuonLeptonSource )THEN
      !$OMP TARGET ENTER DATA MAP(alloc: MaterialSource,MaterialAbsoluteSource)
    END IF
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE(EvalMask,WorkMask,Err,Pack,Unpack,MuonConserved,EtaActive)
    !$ACC ENTER DATA CREATE(NumberSource,NumberFactor,NumberFraction)
    IF( UseMuonLeptonSource )THEN
      !$ACC ENTER DATA CREATE(MaterialSource,MaterialAbsoluteSource)
    END IF
#endif
    unit = 0
    IF( TraceRadiationOnlySolve )THEN
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      WRITE(TraceName,'(A,A,I6.6,A)') TRIM(RadiationOnlyTraceFile), &
        '.coupled_eta.rank',rank,'.dat'
      IF( RadiationOnlySolveCallCount == 1 )THEN
        OPEN(NEWUNIT=unit,FILE=TRIM(TraceName),STATUS='REPLACE',IOSTAT=io)
      ELSE
        OPEN(NEWUNIT=unit,FILE=TRIM(TraceName),POSITION='APPEND',IOSTAT=io)
      END IF
      IF( io /= 0 ) CALL Fail('cannot open coupled-eta trace',0)
      WRITE(unit,'(A)') &
        '# call iteration point eta Ym Ycons Qmu radiation_relative merit lambda Hmu rho_g_cm3 T_MeV'
    END IF

    DO k = 1, MaxIter_inner
      IF( .NOT. ANY(Active) ) EXIT
      CALL Evaluate(Active,X,Eta,F,Q,Bad,.TRUE.)
      IF( ANY(Bad .AND. Active) ) CALL Fail('invalid accepted thermodynamic state', &
        FINDLOC(Bad .AND. Active,.TRUE.,DIM=1))
      WarmT = T

      ! Freeze all row/coordinate scales during this Newton step and its
      ! line search. Recompute them only at the next accepted iterate.
      DO p = 1, nX_G
        IF( .NOT. Active(p) ) CYCLE
        leptonScale = ABS(Ym(p))+ABS(Ym_old(p))
        DO s = 1, nSpecies
        DO eN = 1, nE_G
          os = ((eN-1)+(s-1)*nE_G)*nCR
          densityScale = MAX(ABS(X(os+iCR_N,p)),ABS(Grad(os+iCR_N,p)),SqrtTiny)
          Scale(os+iCR_N,p) = densityScale
          Scale(os+iCR_G1,p) = densityScale*SQRT(g1(p))
          Scale(os+iCR_G2,p) = densityScale*SQRT(g2(p))
          Scale(os+iCR_G3,p) = densityScale*SQRT(g3(p))
          IF( .NOT. UseMuonLeptonSource .AND. &
              (s == iNuM .OR. s == iNuM_Bar) ) leptonScale = leptonScale &
            + AtomicMassUnit/D(p)*W2_S(eN) &
              * (ABS(C_Dnu(eN,s,p))+ABS(X(os+iCR_N,p)) &
                 + ABS(V1(p)*X(os+iCR_G1,p)) &
                 + ABS(V2(p)*X(os+iCR_G2,p))+ABS(V3(p)*X(os+iCR_G3,p)))
        END DO
        END DO
        leptonScale = leptonScale &
          + ABS(AtomicMassUnit/D(p)*LimiterCorrection(iYm,p))
        IF( UseMuonLeptonSource ) leptonScale = leptonScale &
          + ABS(dt*AtomicMassUnit/D(p)*MatterSourceWeight(iYm,p))*MaterialAbsoluteSource(p)
        QScale(p) = MAX(leptonScale,SqrtTiny)
        ! Inventory mode uses an inventory-roundoff bound. Source mode
        ! instead scales this row to the material abundance and reaction
        ! magnitudes, without a radiation-inventory tolerance floor.
        QAcceptance(p) = 64.0_DP*EPSILON(One)*QScale(p)
        ! A direct source is a nonlinear equation, not an algebraically
        ! eliminated inventory. Use the configured nonlinear tolerance on
        ! its OWN scale, without accepting an absolute 1e-8 fraction error.
        IF( UseMuonLeptonSource ) &
          QAcceptance(p) = MAX(Rtol_inner,64.0_DP*EPSILON(One))*QScale(p)
      END DO
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(Active)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(Active)
#endif
      CALL ComputeDnuNorm(Active,Dnu)
      Pending = Active
      CALL CheckConvergence_RadiationOnly &
        (Active,k,nInner,nOuter,nInnerFinal,nOuterFinal,FirstResidual,FinalResidual,Frad)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Frad,Grad,nInnerFinal,nOuterFinal,FirstResidual,FinalResidual)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST(Frad,Grad,nInnerFinal,nOuterFinal,FirstResidual,FinalResidual)
#endif
      DO p = 1, nX_G
        IF( Pending(p) .AND. EtaActive(p) .AND. ABS(Q(p)) > QAcceptance(p) ) &
          Active(p) = .TRUE.
      END DO
      CALL FindRadiationLeptonInvariantViolations_OrderV &
        (EvalMask,WorkMask,LeptonInvariant,Dnu,I1,I2,I3,V1,V2,V3,g1,g2,g3)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(WorkMask)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST(WorkMask)
#endif
      Active = Active .OR. WorkMask
      CALL Norm(Active,F,Hvalue,Merit)
      IF( unit /= 0 )THEN
        DO p = 1, nX_G
          IF( Pending(p) ) WRITE(unit,'(3I9,10ES24.15E3)') &
            RadiationOnlySolveCallCount,k,p,Eta(p),Ym(p),MuonConserved(p), &
            Q(p),FinalResidual(p),Merit(p),Lambda(p),Hvalue(p), &
            D(p)/(Gram/Centimeter**3),BoltzmannConstant*T(p)/MeV
        END DO
        FLUSH(unit)
      END IF
      IF( .NOT. ANY(Active) ) EXIT
      Fbase = F
      Qbase = Q
      Hbase = Hvalue

      IF( k == 1 )THEN
        ! M1 has no unique flux derivative at radiation vacuum. Try the
        ! existing frozen-matter map as a first guess, keeping the implicit
        ! input invariants fixed. This initializes Newton rather than
        ! globalizing a Newton correction: only joint admissibility is
        ! required here, not descent from the potentially vacuum input.
        Pending = Active
        StepX = Zero
        DO p = 1, nX_G
          IF( Active(p) ) StepX(:,p) = Grad(:,p)-X(:,p)
        END DO
        Lambda = One
        DO bt = 1, 12
          DO p = 1, nX_G
            TrialX(:,p) = X(:,p)+Lambda(p)*StepX(:,p)
          END DO
          CALL ProjectTrial(Pending,TrialX)
          CALL Evaluate(Pending,TrialX,Eta,TrialF,TrialQ,Bad,.TRUE.)
          CALL Norm(Pending,TrialF,Hvalue,TrialMerit)
          DO p = 1, nX_G
            IF( .NOT. Pending(p) ) CYCLE
            IF( .NOT. Bad(p) )THEN
              X(:,p) = TrialX(:,p)
              Pending(p) = .FALSE.
              CYCLE
            END IF
            Lambda(p) = Half*Lambda(p)
          END DO
          IF( .NOT. ANY(Pending) ) EXIT
        END DO
        ! Re-evaluate every point at its own accepted state before Newton;
        ! rejected predictor trials must not supply coefficients to a step.
        CYCLE
      END IF

      ! b=F_eta and d=Q_eta. The radiation block A includes the
      ! temperature, EOS, rates, velocity, and energy-group response.
      EtaPerturb = One
      CALL Differentiate(Active .AND. EtaActive,Zero*X,EtaPerturb,B,EtaJac,.TRUE.)

      B = B/Scale
      EtaJac = EtaJac/QScale
      CALL RadiationLinearSolve(-Fbase/Scale,Z)
      CALL RadiationLinearSolve(B,Y)
      CALL Differentiate(Active,Scale*Z,Zero*Eta,TrialF,QZ,.TRUE.)
      CALL Differentiate(Active,Scale*Y,Zero*Eta,TrialF,QY,.TRUE.)
      QZ = QZ/QScale
      QY = QY/QScale
      StepEta = Zero
      DO p = 1, nX_G
        IF( .NOT. Active(p) .OR. .NOT. EtaActive(p) ) CYCLE
        denom = EtaJac(p)-QY(p)
        IF( ABS(Hbase(p)+QZ(p)*QScale(p)) <= QAcceptance(p) ) CYCLE
        IF( ABS(denom) <= 128.0_DP*EPSILON(One)*MAX(ABS(EtaJac(p)),ABS(QY(p)),TINY(One)) )THEN
          CALL Fail('singular muon Schur complement',p)
        ELSE
          StepEta(p) = -(Hbase(p)/QScale(p)+QZ(p))/denom
        END IF
      END DO
      DO p = 1, nX_G
        ! Scalar trust modification: a susceptibility-limited correction
        ! can be enormous at pair symmetry. Compress only that scalar
        ! displacement, then reconstruct the radiation block response.
        ! This retains the Schur-Newton step to cubic order near a root.
        StepEta(p) = ASINH(StepEta(p))
        StepX(:,p) = Scale(:,p)*(Z(:,p)-Y(:,p)*StepEta(p))
        Lambda(p) = One
      END DO

      Pending = Active
      DO bt = 1, 24
        DO p = 1, nX_G
          TrialEta(p) = Eta(p)+Lambda(p)*StepEta(p)
          TrialX(:,p) = X(:,p)+Lambda(p)*StepX(:,p)
        END DO
        CALL ProjectTrial(Pending,TrialX)
        CALL Evaluate(Pending,TrialX,TrialEta,TrialF,TrialQ,Bad,.TRUE.)
        CALL Norm(Pending,TrialF,Hvalue,TrialMerit)
        DO p = 1, nX_G
          IF( .NOT. Pending(p) ) CYCLE
          IF( .NOT. Bad(p) )THEN
            IF( TrialMerit(p) <= (One-Armijo*Lambda(p))*Merit(p) )THEN
              X(:,p) = TrialX(:,p)
              Eta(p) = TrialEta(p)
              Pending(p) = .FALSE.
              CYCLE
            END IF
          END IF
          Lambda(p) = Half*Lambda(p)
        END DO
        IF( .NOT. ANY(Pending) ) EXIT
      END DO
      IF( ANY(Pending) ) CALL Fail('joint radiation/eta line search exhausted', &
        FINDLOC(Pending,.TRUE.,DIM=1))
    END DO
    IF( ANY(Active) ) CALL Fail('coupled radiation/eta Newton iteration limit', &
      FINDLOC(Active,.TRUE.,DIM=1))
    IF( unit /= 0 ) CLOSE(unit)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO(Active)
    !$OMP TARGET EXIT DATA MAP(delete: EvalMask,WorkMask,Err,Pack,Unpack,MuonConserved,EtaActive)
    !$OMP TARGET EXIT DATA MAP(delete: NumberSource,NumberFactor,NumberFraction)
    IF( UseMuonLeptonSource )THEN
      !$OMP TARGET EXIT DATA MAP(delete: MaterialSource,MaterialAbsoluteSource)
    END IF
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE(Active)
    !$ACC EXIT DATA DELETE(EvalMask,WorkMask,Err,Pack,Unpack,MuonConserved,EtaActive)
    !$ACC EXIT DATA DELETE(NumberSource,NumberFactor,NumberFraction)
    IF( UseMuonLeptonSource )THEN
      !$ACC EXIT DATA DELETE(MaterialSource,MaterialAbsoluteSource)
    END IF
#endif

  CONTAINS

    SUBROUTINE Fail(reason,point)
      CHARACTER(*), INTENT(in) :: reason
      INTEGER, INTENT(in) :: point
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      WRITE(ERROR_UNIT,'(A,I0,3A,I0,A,I0)') 'Coupled radiation/eta failure, rank ', &
        rank,': ',reason,'; point ',point,'; Newton iteration ',k
      IF( point > 0 )THEN
        WRITE(ERROR_UNIT,'(A,4ES24.15E3)') &
          ' rho [g/cm3], accepted T [MeV], eta, Qmu: ', &
          D(point)/(Gram/Centimeter**3),BoltzmannConstant*WarmT(point)/MeV,Eta(point),Q(point)
        WRITE(ERROR_UNIT,'(A,I0)') ' evaluation failure code: ',EvalReason(point)
        IF( EvalReason(point) == 1 )THEN
          DO s=1,nSpecies
          DO eN=1,nE_G
            densityScale=SQRT(g1(point)*I1(eN,s,point)**2+g2(point)*I2(eN,s,point)**2 &
                             +g3(point)*I3(eN,s,point)**2)
            IF( Dnu(eN,s,point) < Zero .OR. Dnu(eN,s,point) > DnuMax .OR. &
                densityScale > Dnu(eN,s,point)*(One+64*EPSILON(One)) )THEN
              WRITE(ERROR_UNIT,'(A,2I6,3ES24.15E3)') ' inadmissible species/energy, N, |I|, Nmax: ', &
                s,eN,Dnu(eN,s,point),densityScale,DnuMax
            END IF
          END DO
          END DO
        END IF
        IF( INDEX(reason,'line search') > 0 )THEN
          WRITE(ERROR_UNIT,'(A,4ES24.15E3)') ' base/trial merit, delta_eta, lambda: ', &
            Merit(point),TrialMerit(point),StepEta(point),Lambda(point)
          WRITE(ERROR_UNIT,'(A,4ES24.15E3)') ' Schur d, cy, S, rhs: ', &
            EtaJac(point),QY(point),EtaJac(point)-QY(point), &
            -Hbase(point)/QScale(point)-QZ(point)
        END IF
      END IF
      FLUSH(ERROR_UNIT)
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    END SUBROUTINE Fail

    SUBROUTINE Evaluate(mask,state,etaState,residual,qResidual,failed,withRadiation)
      LOGICAL, INTENT(in) :: mask(:),withRadiation
      REAL(DP), INTENT(in) :: state(:,:),etaState(:)
      REAL(DP), INTENT(out) :: residual(:,:),qResidual(:)
      LOGICAL, INTENT(out) :: failed(:)
      INTEGER :: j,a,tup,off,nPacked,localError,family,pivot,energyPivot
      REAL(DP) :: fluxnorm,susceptibility,invariantChange,numberChange,sourceSum
      EvalMask = mask
      Err = 0
      WHERE(mask) EvalReason = 0
      failed = .FALSE.
      residual = Zero
      qResidual = Zero
      WHERE(mask) Hvalue = Zero
      WHERE(mask) Hsource = Zero
      WHERE(mask) T = WarmT
      DO j = 1, nX_G
        IF( .NOT. mask(j) ) CYCLE
        DO a = 1, nSpecies
        DO tup = 1, nE_G
          off = ((tup-1)+(a-1)*nE_G)*nCR
          IF( .NOT. ALL(IEEE_IS_FINITE(state(off+1:off+nCR,j))) )THEN
            Err(j) = 1
            CYCLE
          END IF
          Dnu(tup,a,j) = state(off+iCR_N,j)
          I1(tup,a,j) = state(off+iCR_G1,j)/g1(j)
          I2(tup,a,j) = state(off+iCR_G2,j)/g2(j)
          I3(tup,a,j) = state(off+iCR_G3,j)/g3(j)
          fluxnorm = HYPOT(HYPOT(state(off+iCR_G1,j)/SQRT(g1(j)), &
            state(off+iCR_G2,j)/SQRT(g2(j))),state(off+iCR_G3,j)/SQRT(g3(j)))
          IF( Dnu(tup,a,j) < Zero .OR. Dnu(tup,a,j) > DnuMax &
                .OR. fluxnorm > Dnu(tup,a,j)*(One+64*EPSILON(One)) ) Err(j) = 1
        END DO
        END DO
      END DO
      WHERE(mask .AND. Err /= 0) EvalReason = 1
      EvalMask = EvalMask .AND. Err == 0
      CALL TimersStart(Timer_Collisions_MatterFromRadiation)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO(EvalMask,Err,Dnu,I1,I2,I3,T,MuonConserved,EtaActive)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE(EvalMask,Err,Dnu,I1,I2,I3,T,MuonConserved,EtaActive)
#endif
      CALL RecoverMatterFromRadiationState_OrderV &
        (EvalMask,Fmatter,Gmatter,Dnu,I1,I2,I3,D,Ye,Ym,E,V1,V2,V3,g1,g2,g3, &
         LimiterCorrection,Err,RetainMuonFraction_Option=EtaActive, &
         MuonFractionConserved_Option=MuonConserved)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM(Ye,Ym,E,V1,V2,V3,Err,MuonConserved,U_Ym,G_Ym,Omega)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST(Ye,Ym,E,V1,V2,V3,Err,MuonConserved,U_Ym,G_Ym,Omega)
#endif
      WHERE(EvalMask .AND. Err /= 0) EvalReason = 2
      WorkMask = EvalMask .AND. .NOT. EtaActive .AND. Err == 0
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(WorkMask)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(WorkMask)
#endif
      CALL CreatePackIndex(WorkMask,nPacked,Pack,Unpack)
      IF( nPacked > 0 ) CALL UpdateTemperature_Packed &
        (D,E,Ye,Ym,T,WorkMask,nPacked,Pack,Unpack,Err,DisplayErrors_Option=.FALSE.)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(T,Err)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST(T,Err)
#endif
      WHERE(WorkMask .AND. Err /= 0) EvalReason = 4
      DO j = 1, nX_G
        IF( .NOT. EvalMask(j) .OR. Err(j) /= 0 ) CYCLE
        IF( EtaActive(j) )THEN
          CALL RecoverMuonEtaTemperature_OrderV &
            (D(j),E(j),Ye(j),etaState(j),T(j),Ym(j),susceptibility,localError)
          Err(j) = localError
          IF( localError /= 0 ) EvalReason(j) = 3
          qResidual(j) = Ym(j)-MuonConserved(j)
          IF( localError == 0 )THEN
            U_Ym(j) = Ym(j)
            G_Ym(j) = Ym(j)
          END IF
        END IF
      END DO
      failed = mask .AND. Err /= 0
      EvalMask = mask .AND. .NOT. failed
      CALL TimersStop(Timer_Collisions_MatterFromRadiation)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO(EvalMask,Err,T,Ye,Ym,E,V1,V2,V3,U_Ym,G_Ym)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE(EvalMask,Err,T,Ye,Ym,E,V1,V2,V3,U_Ym,G_Ym)
#endif
      CALL CreatePackIndex(EvalMask,nPacked,Pack,Unpack)
      IF( withRadiation .AND. nPacked > 0 )THEN
        CALL TimersStart(Timer_Collisions_ComputeOpacity)
        IF( .NOT. FreezeOpacities ) CALL ComputeOpacities_Packed &
          (D,T,Ye,Ym,SqrtGm,EvalMask,nPacked,Pack,Unpack, &
           T0=T_old,Ye0=Ye_old,Ym0=Ym_old,Opacity_Mask_Input=Opacity_Mask)
        CALL TimersStop(Timer_Collisions_ComputeOpacity)
        CALL TimersStart(Timer_Collisions_ComputeRates)
        IF( FreezeOpacities )THEN
          CALL ComputeRates_Packed &
            (D,T_Opacity,Dnu,I1,I2,I3,EvalMask,nPacked,Pack,Unpack,nInitialPacked,Opacity_Mask)
        ELSE
          CALL ComputeRates_Packed &
            (D,T,Dnu,I1,I2,I3,EvalMask,nPacked,Pack,Unpack,nPacked,Opacity_Mask)
        END IF
        CALL TimersStop(Timer_Collisions_ComputeRates)
        CALL TimersStart(Timer_Collisions_NeutrinoRHS)
        CALL ComputeNeutrinoRHS_OrderV &
          (EvalMask,Frad,Grad,dt,Dnu,I1,I2,I3,V1,V2,V3,g1,g2,g3, &
           NumberSource_Option=NumberSource,NumberFactor_Option=NumberFactor, &
           NumberFraction_Option=NumberFraction)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE FROM(Frad,NumberSource,NumberFactor,NumberFraction)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE HOST(Frad,NumberSource,NumberFactor,NumberFraction)
#endif
        RawFm = Frad
        CALL ProjectRadiationOnlyFixedPointMap_OrderV &
          (EvalMask,Frad,Grad,Dnu,I1,I2,I3,g1,g2,g3)
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE FROM(Frad,Grad)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE HOST(Frad,Grad)
#endif
        IF( UseMuonLeptonSource )THEN
          CALL ComputeMuonMaterialSource_OrderV(EvalMask,Dnu,MaterialSource,MaterialAbsoluteSource)
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM(MaterialSource,MaterialAbsoluteSource)
#elif defined(THORNADO_OACC)
          !$ACC UPDATE HOST(MaterialSource,MaterialAbsoluteSource)
#endif
        END IF
        DO j = 1, nX_G
          IF( EvalMask(j) )THEN
            residual(:,j) = -Frad(:,j)
            IF( EtaActive(j) )THEN
              sourceSum = Zero
              IF( .NOT. UseMuonLeptonSource )THEN
              DO a = iNuM,iNuM_Bar
              DO tup = 1,nE_G
                off = ((tup-1)+(a-1)*nE_G)*nCR
                numberChange = state(off+iCR_N,j)-C_Dnu(tup,a,j) &
                  +V1(j)*state(off+iCR_G1,j)+V2(j)*state(off+iCR_G2,j) &
                  +V3(j)*state(off+iCR_G3,j)
                sourceSum = sourceSum+W2_S(tup)*LeptonNumber(a) &
                  * (NumberFraction(tup,a,j)*numberChange &
                    +dt*NumberSource(tup,a,j)/(Omega(j)*NumberFactor(tup,a,j)) &
                    +(Frad(off+iCR_N,j)-RawFm(off+iCR_N,j))/Omega(j))
              END DO
              END DO
              END IF
              IF( UseMuonLeptonSource )THEN
                Hsource(j) = -AtomicMassUnit/D(j)*MatterSourceWeight(iYm,j) &
                  * (dt*MaterialSource(j)+LimiterCorrection(iYm,j))
              ELSE
                Hsource(j) = AtomicMassUnit/D(j)*MatterSourceWeight(iYm,j) &
                  * (sourceSum-LimiterCorrection(iYm,j))
              END IF
              Hvalue(j) = (Ym(j)-Ym_old(j))+Hsource(j)
              IF( UseMuonLeptonSource ) qResidual(j) = Hvalue(j)
            END IF
            DO family = 1, MIN(nLeptonFamilies,nSpecies/2)
              pivot = InvariantPivot(family,j)
              IF( pivot == 0 ) CYCLE
              invariantChange = Zero
              DO a = 2*family-1, 2*family
              DO tup = 1, nE_G
                off = ((tup-1)+(a-1)*nE_G)*nCR
                numberChange = state(off+iCR_N,j)-C_Dnu(tup,a,j) &
                  + V1(j)*state(off+iCR_G1,j)+V2(j)*state(off+iCR_G2,j) &
                  + V3(j)*state(off+iCR_G3,j)
                invariantChange = invariantChange+W2_S(tup)*LeptonNumber(a)*numberChange
              END DO
              END DO
              energyPivot = MOD((pivot-iCR_N)/nCR,nE_G)+1
              residual(pivot,j) = invariantChange/W2_S(energyPivot)
            END DO
            Frad(:,j) = -residual(:,j)
            IF( .NOT. ALL(IEEE_IS_FINITE(residual(:,j))) )THEN
              failed(j) = .TRUE.
              EvalReason(j) = 5
            END IF
          END IF
        END DO
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET UPDATE TO(Frad)
#elif defined(THORNADO_OACC)
        !$ACC UPDATE DEVICE(Frad)
#endif
        CALL TimersStop(Timer_Collisions_NeutrinoRHS)
      END IF
    END SUBROUTINE Evaluate

    SUBROUTINE ProjectTrial(mask,state)
      LOGICAL, INTENT(in) :: mask(:)
      REAL(DP), INTENT(inout) :: state(:,:)
      INTEGER :: j,a,tup,off
      REAL(DP) :: occupation,f1,f2,f3
      DO j = 1,nX_G
        IF( .NOT. mask(j) ) CYCLE
        DO a = 1,nSpecies
        DO tup = 1,nE_G
          off = ((tup-1)+(a-1)*nE_G)*nCR
          occupation = state(off+iCR_N,j)
          f1 = state(off+iCR_G1,j)/g1(j)
          f2 = state(off+iCR_G2,j)/g2(j)
          f3 = state(off+iCR_G3,j)/g3(j)
          CALL ProjectPrimitiveRadiationMoment_OrderV(occupation,f1,f2,f3, &
            g1(j),g2(j),g3(j),DnuMax)
          state(off+iCR_N,j) = occupation
          state(off+iCR_G1,j) = g1(j)*f1
          state(off+iCR_G2,j) = g2(j)*f2
          state(off+iCR_G3,j) = g3(j)*f3
        END DO
        END DO
      END DO
    END SUBROUTINE ProjectTrial

    SUBROUTINE Norm(mask,residual,qResidual,value)
      LOGICAL, INTENT(in) :: mask(:)
      REAL(DP), INTENT(in) :: residual(:,:),qResidual(:)
      REAL(DP), INTENT(out) :: value(:)
      INTEGER :: j
      value = Zero
      DO j = 1, nX_G
        IF( mask(j) ) value(j) = SQRT(SUM((residual(:,j)/Scale(:,j) &
          /MAX(Rtol_inner,EPSILON(One)))**2)/n_FP_inner &
          +(qResidual(j)/QAcceptance(j))**2)
      END DO
    END SUBROUTINE Norm

    SUBROUTINE Differentiate(mask,dx,de,df,dq,withRadiation)
      LOGICAL, INTENT(in) :: mask(:),withRadiation
      REAL(DP), INTENT(in) :: dx(:,:),de(:)
      REAL(DP), INTENT(out) :: df(:,:),dq(:)
      REAL(DP), ALLOCATABLE :: h(:),xp(:,:),ep(:),fp(:,:),qp(:),fmProbe(:,:),qmProbe(:)
      REAL(DP), ALLOCATABLE :: inward(:,:),forwardFirst(:),sourceFirst(:)
      LOGICAL, ALLOCATABLE :: todo(:),badProbe(:),badFirst(:)
      INTEGER :: j,attempt,a,tup,off
      REAL(DP) :: fluxnorm,directionnorm,inwardDensity
      ALLOCATE(h(nX_G),xp(n_FP_inner,nX_G),ep(nX_G),fp(n_FP_inner,nX_G),qp(nX_G))
      ALLOCATE(fmProbe(n_FP_inner,nX_G),qmProbe(nX_G),inward(n_FP_inner,nX_G))
      ALLOCATE(forwardFirst(nX_G),sourceFirst(nX_G))
      ALLOCATE(todo(nX_G),badProbe(nX_G),badFirst(nX_G))
      todo = mask
      h = FDStep
      df = Zero
      dq = Zero
      inward = Zero
      DO j = 1, nX_G
        IF( MAXVAL(ABS(dx(:,j)/Scale(:,j)))+ABS(de(j)) == Zero ) todo(j) = .FALSE.
        IF( .NOT. todo(j) ) CYCLE
        h(j) = FDStep/MAX(One,MAXVAL(ABS(dx(:,j)/Scale(:,j))),ABS(de(j)))
        DO a = 1, nSpecies
        DO tup = 1, nE_G
          off = ((tup-1)+(a-1)*nE_G)*nCR
          inwardDensity = Zero
          directionnorm = HYPOT(HYPOT(dx(off+iCR_G1,j)/SQRT(g1(j)), &
            dx(off+iCR_G2,j)/SQRT(g2(j))),dx(off+iCR_G3,j)/SQRT(g3(j)))
          IF( X(off+iCR_N,j) >= DnuMax-FDStep*Scale(off+iCR_N,j) ) &
            inwardDensity = -ABS(dx(off+iCR_N,j))
          IF( X(off+iCR_N,j) <= FDStep*Scale(off+iCR_N,j) ) &
            inwardDensity = Two*(ABS(dx(off+iCR_N,j))+directionnorm)
          inward(off+iCR_N,j) = inwardDensity
          fluxnorm = HYPOT(HYPOT(X(off+iCR_G1,j)/SQRT(g1(j)), &
            X(off+iCR_G2,j)/SQRT(g2(j))),X(off+iCR_G3,j)/SQRT(g3(j)))
          IF( X(off+iCR_N,j) > FDStep*Scale(off+iCR_N,j) .AND. fluxnorm > Zero .AND. &
              fluxnorm >= X(off+iCR_N,j)*(One-10*FDStep) )THEN
            inward(off+iCR_G1:off+iCR_G3,j) = -X(off+iCR_G1:off+iCR_G3,j)/fluxnorm &
              * (Two*ABS(dx(off+iCR_N,j))+Two*directionnorm)
          END IF
        END DO
        END DO
      END DO
      DO attempt = 1, 24
        DO j = 1, nX_G
          xp(:,j) = X(:,j)+Half*h(j)*(inward(:,j)+dx(:,j))
          ep(j) = Eta(j)+Half*h(j)*de(j)
        END DO
        CALL Evaluate(todo,xp,ep,fp,qp,badFirst,withRadiation)
        forwardFirst = Ym
        sourceFirst = Hsource
        DO j = 1, nX_G
          xp(:,j) = X(:,j)+Half*h(j)*(inward(:,j)-dx(:,j))
          ep(j) = Eta(j)-Half*h(j)*de(j)
        END DO
        CALL Evaluate(todo,xp,ep,fmProbe,qmProbe,badProbe,withRadiation)
        DO j = 1, nX_G
          IF( .NOT. todo(j) ) CYCLE
          IF( badProbe(j) .OR. badFirst(j) )THEN
            h(j) = Half*h(j)
          ELSE
            IF( withRadiation ) df(:,j) = (fp(:,j)-fmProbe(:,j))/h(j)
            IF( EtaActive(j) ) dq(j) = &
              ((forwardFirst(j)-Ym(j))+(sourceFirst(j)-Hsource(j)))/h(j)
            todo(j) = .FALSE.
          END IF
        END DO
        IF( .NOT. ANY(todo) ) EXIT
      END DO
      IF( ANY(todo) )THEN
        j=FINDLOC(todo,.TRUE.,DIM=1)
        DO a=1,nSpecies
        DO tup=1,nE_G
          off=((tup-1)+(a-1)*nE_G)*nCR
          IF( Dnu(tup,a,j) < Zero .OR. SQRT(g1(j)*I1(tup,a,j)**2 &
              +g2(j)*I2(tup,a,j)**2+g3(j)*I3(tup,a,j)**2) &
              > Dnu(tup,a,j)*(One+64*EPSILON(One)) )THEN
            WRITE(ERROR_UNIT,'(A,3I6,5ES24.15E3)') ' probe point/species/energy, h, baseN, scaleN, dxN, inwardN: ', &
              j,a,tup,h(j),X(off+iCR_N,j),Scale(off+iCR_N,j),dx(off+iCR_N,j),inward(off+iCR_N,j)
            WRITE(ERROR_UNIT,'(A,3ES24.15E3)') ' base I1, dxI1, inwardI1: ', &
              X(off+iCR_G1,j),dx(off+iCR_G1,j),inward(off+iCR_G1,j)
          END IF
        END DO
        END DO
        CALL Fail('no admissible Jacobian probe',j)
      END IF
    END SUBROUTINE Differentiate

    SUBROUTINE RadiationLinearSolve(rhs,solution)
      REAL(DP), INTENT(in) :: rhs(:,:)
      REAL(DP), INTENT(out) :: solution(:,:)
      REAL(DP), ALLOCATABLE :: basis(:,:,:),hess(:,:,:),cs(:,:),sn(:,:),r(:,:)
      REAL(DP), ALLOCATABLE :: w(:,:),jq(:),beta(:),coeff(:,:),av(:,:)
      LOGICAL, ALLOCATABLE :: todo(:)
      INTEGER, ALLOCATABLE :: used(:)
      INTEGER :: j,m,a,pass
      REAL(DP) :: dot,normw,t1,t2,rotation
      ALLOCATE(basis(n_FP_inner,KrylovDimension+1,nX_G))
      ALLOCATE(hess(KrylovDimension+1,KrylovDimension,nX_G))
      ALLOCATE(cs(KrylovDimension,nX_G),sn(KrylovDimension,nX_G))
      ALLOCATE(r(KrylovDimension+1,nX_G),coeff(KrylovDimension,nX_G),used(nX_G))
      ALLOCATE(w(n_FP_inner,nX_G),av(n_FP_inner,nX_G),jq(nX_G),beta(nX_G),todo(nX_G))
      basis = Zero
      hess = Zero
      r = Zero
      coeff = Zero
      used = 0
      solution = Zero
      todo = Active
      DO j = 1, nX_G
        beta(j) = SQRT(SUM(rhs(:,j)**2))
        IF( beta(j) <= TINY(One) ) todo(j) = .FALSE.
        IF( todo(j) ) basis(:,1,j) = rhs(:,j)/beta(j)
        r(1,j) = beta(j)
      END DO
      DO m = 1, KrylovDimension
        IF( .NOT. ANY(todo) ) EXIT
        CALL Differentiate(todo,Scale*basis(:,m,:),Zero*Eta,av,jq,.TRUE.)
        w = av/Scale
        DO j = 1, nX_G
          IF( .NOT. todo(j) ) CYCLE
          ! Twice-reorthogonalized Arnoldi protects a tiny Schur response
          ! from loss of orthogonality in strongly relaxing spectra.
          DO pass = 1, 2
            DO a = 1, m
              dot = SUM(w(:,j)*basis(:,a,j))
              hess(a,m,j) = hess(a,m,j)+dot
              w(:,j) = w(:,j)-dot*basis(:,a,j)
            END DO
          END DO
          normw = SQRT(SUM(w(:,j)**2))
          hess(m+1,m,j) = normw
          IF( normw > TINY(One) ) basis(:,m+1,j) = w(:,j)/normw
          DO a = 1, m-1
            t1 = cs(a,j)*hess(a,m,j)+sn(a,j)*hess(a+1,m,j)
            t2 = -sn(a,j)*hess(a,m,j)+cs(a,j)*hess(a+1,m,j)
            hess(a,m,j) = t1
            hess(a+1,m,j) = t2
          END DO
          rotation = SQRT(hess(m,m,j)**2+hess(m+1,m,j)**2)
          IF( rotation <= TINY(One) ) CALL Fail('singular radiation Newton block',j)
          cs(m,j) = hess(m,m,j)/rotation
          sn(m,j) = hess(m+1,m,j)/rotation
          hess(m,m,j) = rotation
          hess(m+1,m,j) = Zero
          r(m+1,j) = -sn(m,j)*r(m,j)
          r(m,j) = cs(m,j)*r(m,j)
          used(j) = m
          IF( ABS(r(m+1,j)) <= 1.0e-7_DP*beta(j) .OR. normw <= TINY(One) ) todo(j) = .FALSE.
        END DO
      END DO
      DO j = 1, nX_G
        IF( .NOT. Active(j) .OR. used(j) == 0 ) CYCLE
        IF( ABS(r(used(j)+1,j)) > 1.0e-5_DP*beta(j) ) &
          CALL Fail('radiation Krylov block did not converge',j)
        DO a = used(j), 1, -1
          coeff(a,j) = (r(a,j)-SUM(hess(a,a+1:used(j),j)*coeff(a+1:used(j),j)))/hess(a,a,j)
        END DO
        DO a = 1, used(j)
          solution(:,j) = solution(:,j)+coeff(a,j)*basis(:,a,j)
        END DO
      END DO
    END SUBROUTINE RadiationLinearSolve

  END SUBROUTINE SolveRadiationMuonEtaSchur_OrderV


  SUBROUTINE RecoverMuonEtaTemperature_OrderV(D,Etarget,Ye,Eta,T,Ym,Chi,Error)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(in) :: D,Etarget,Ye,Eta
    REAL(DP), INTENT(inout) :: T,Ym
    REAL(DP), INTENT(out) :: Chi
    INTEGER, INTENT(out) :: Error
    REAL(DP) :: tlo,thi,tc,tp,ec,ep,yc,yp,chp,loResidual,hiResidual,residual
    REAL(DP) :: slope,delta,tNext,initialT
    INTEGER :: iteration,probeError,attempt
    LOGICAL :: haveLo,haveHi

    ! Recover temperature at fixed eta using the forward EOS, including
    ! the T dependence of the net muon fraction. This is an EOS scalar
    ! recovery, not a second muon-conservation iteration.
    Error = 1
    Chi = Zero
    initialT = MAX(Min_T,MIN(T,Max_T))
    tc = initialT
    tlo = Min_T
    thi = Max_T
    haveLo = .FALSE.
    haveHi = .FALSE.
    loResidual = Zero
    hiResidual = Zero
    DO iteration = 1, 64
      yc = Zero
      CALL ComputeMuonFractionFromDegeneracy_TABLE(D,tc,Ye,Eta,yc,Chi,probeError)
      IF( probeError /= 0 ) RETURN
      CALL ComputeSpecificInternalEnergy_TABLE(D,tc,Ye,yc,ec)
      residual = ec-Etarget
      IF( ABS(residual) <= 100*EPSILON(One)*MAX(ABS(ec),ABS(Etarget),SqrtTiny) )THEN
        T = tc
        Ym = yc
        Error = 0
        RETURN
      END IF
      IF( residual < Zero )THEN
        tlo = tc
        loResidual = residual
        haveLo = .TRUE.
      ELSE
        thi = tc
        hiResidual = residual
        haveHi = .TRUE.
      END IF
      ! Numerical heat capacity along eta=constant; avoid the legacy EOS
      ! derivative path, whose muon-composition derivatives are incomplete.
      tp = MIN(Max_T,tc*(One+1.0e-4_DP))
      IF( tp == tc ) tp = MAX(Min_T,tc*(One-1.0e-4_DP))
      yp = yc
      CALL ComputeMuonFractionFromDegeneracy_TABLE(D,tp,Ye,Eta,yp,chp,probeError)
      slope = Zero
      IF( probeError == 0 .AND. tp /= tc )THEN
        CALL ComputeSpecificInternalEnergy_TABLE(D,tp,Ye,yp,ep)
        slope = (ep-ec)/(tp-tc)
      END IF
      IF( slope > Zero )THEN
        delta = MAX(-Half*tc,MIN(Half*tc,-residual/slope))
        tNext = tc+delta
      ELSE IF( residual < Zero )THEN
        tNext = MIN(Max_T,1.5_DP*tc)
      ELSE
        tNext = MAX(Min_T,Half*tc)
      END IF
      IF( haveLo .AND. haveHi )THEN
        IF( tNext <= tlo .OR. tNext >= thi ) &
          tNext = tlo-loResidual*(thi-tlo)/(hiResidual-loResidual)
        tNext = MAX(tlo+0.01_DP*(thi-tlo),MIN(tNext,thi-0.01_DP*(thi-tlo)))
      END IF
      tNext = MAX(Min_T,MIN(tNext,Max_T))
      ! Composition can restrict the valid temperature interval at fixed
      ! eta. Shorten this EOS proposal rather than accepting a clipped Ymu.
      DO attempt = 1, 24
        yp = yc
        CALL ComputeMuonFractionFromDegeneracy_TABLE(D,tNext,Ye,Eta,yp,chp,probeError)
        IF( probeError == 0 ) EXIT
        tNext = Half*(tc+tNext)
      END DO
      IF( probeError /= 0 .OR. tNext == tc ) RETURN
      tc = tNext
    END DO
  END SUBROUTINE RecoverMuonEtaTemperature_OrderV


  SUBROUTINE PrepareRadiationOnlySchurNewton_OrderV &
    ( MASK, Gm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      HasPrevious, PreviousInput, PreviousMap, PreviousResidual, &
      Gain, Damping )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: HasPrevious
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      PreviousInput, PreviousMap, PreviousResidual, Gain, Damping

    REAL(DP), PARAMETER :: InitialDamping = 0.20_DP
    REAL(DP), PARAMETER :: MinimumDamping = 0.02_DP
    REAL(DP), PARAMETER :: MaximumDamping = One
    REAL(DP), PARAMETER :: ResidualGrowthLimit = 1.25_DP
    REAL(DP), PARAMETER :: SecantFloor = 128.0_DP * EPSILON(One)

    INTEGER  :: iN_E, iN_X, iOS, iS
    REAL(DP) :: CurrentInput, CurrentMap, CurrentResidual
    REAL(DP) :: DeltaInput, DeltaMap, Denominator, LocalDamping
    REAL(DP) :: LocalGain, Scale, vDotI
    LOGICAL  :: ReliableSecant

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE( iN_E, iOS, iS, CurrentInput, CurrentMap, &
    !$OMP          CurrentResidual, DeltaInput, DeltaMap, Denominator, &
    !$OMP          LocalDamping, LocalGain, Scale, vDotI, &
    !$OMP          ReliableSecant )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( iN_E, iOS, iS, CurrentInput, CurrentMap, &
    !$ACC          CurrentResidual, DeltaInput, DeltaMap, Denominator, &
    !$ACC          LocalDamping, LocalGain, Scale, vDotI, &
    !$ACC          ReliableSecant )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iN_E, iOS, iS, CurrentInput, CurrentMap, &
    !$OMP          CurrentResidual, DeltaInput, DeltaMap, Denominator, &
    !$OMP          LocalDamping, LocalGain, Scale, vDotI, &
    !$OMP          ReliableSecant )
#endif
    DO iN_X = 1, nX_G

      IF( .NOT. MASK(iN_X) ) CYCLE

      Gain(iN_X) = Zero
      Damping(iN_X) = One

      IF( nSpecies < iNuM_Bar ) CYCLE
      IF( .NOT. MatterSource_Mask(iYm,iN_X) ) CYCLE

      CurrentInput = Zero
      CurrentMap = Zero
      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        iOS = ( (iN_E-1) + (iS-1)*nE_G ) * nCR
        vDotI = V_u_1(iN_X) * Gm_dd_11(iN_X) &
                  * Inu_u_1(iN_E,iS,iN_X) &
              + V_u_2(iN_X) * Gm_dd_22(iN_X) &
                  * Inu_u_2(iN_E,iS,iN_X) &
              + V_u_3(iN_X) * Gm_dd_33(iN_X) &
                  * Inu_u_3(iN_E,iS,iN_X)
        CurrentInput = CurrentInput + W2_S(iN_E) * LeptonNumber(iS) &
          * ( Dnu(iN_E,iS,iN_X) + vDotI )
        CurrentMap = CurrentMap + W2_S(iN_E) * LeptonNumber(iS) &
          * ( Gm(iOS+iCR_N,iN_X) &
              + V_u_1(iN_X) * Gm(iOS+iCR_G1,iN_X) &
              + V_u_2(iN_X) * Gm(iOS+iCR_G2,iN_X) &
              + V_u_3(iN_X) * Gm(iOS+iCR_G3,iN_X) )
      END DO
      END DO

      CurrentResidual = CurrentMap - CurrentInput
      LocalDamping = InitialDamping
      LocalGain = Zero

      IF( HasPrevious(iN_X) )THEN
        DeltaInput = CurrentInput - PreviousInput(iN_X)
        DeltaMap = CurrentMap - PreviousMap(iN_X)
        Scale = MAX( ABS(CurrentInput), ABS(CurrentMap), &
                     ABS(PreviousInput(iN_X)), &
                     ABS(PreviousMap(iN_X)), SqrtTiny )
        ReliableSecant = ABS(DeltaInput) > SecantFloor * Scale
        IF( ReliableSecant )THEN
          LocalGain = DeltaMap / DeltaInput
          Denominator = One - LocalGain
          ReliableSecant = LocalGain == LocalGain &
            .AND. ABS(LocalGain) < HUGE(One) &
            .AND. Denominator > SecantFloor
        END IF
        IF( ReliableSecant ) &
          LocalDamping = MIN( MaximumDamping, &
            MAX( MinimumDamping, One / Denominator ) )
        IF( ABS(CurrentResidual) > ResidualGrowthLimit &
              * MAX(ABS(PreviousResidual(iN_X)),SqrtTiny) ) &
          LocalDamping = MIN( LocalDamping, &
            MAX(MinimumDamping,Half*Damping(iN_X)) )
      END IF

      Gain(iN_X) = LocalGain
      Damping(iN_X) = LocalDamping
      PreviousInput(iN_X) = CurrentInput
      PreviousMap(iN_X) = CurrentMap
      PreviousResidual(iN_X) = CurrentResidual
      HasPrevious(iN_X) = .TRUE.

    END DO

  END SUBROUTINE PrepareRadiationOnlySchurNewton_OrderV


  SUBROUTINE ApplyRadiationOnlySchurNewtonPreconditioner_OrderV &
    ( MASK, Fm, Gm, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Damping )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Damping

    INTEGER  :: iN_E, iN_X, iOS, iS
    REAL(DP) :: StateMoment

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iOS, StateMoment )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iOS, StateMoment )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE( iOS, StateMoment )
#endif
    DO iN_X = 1, nX_G
    DO iS = iNuM, MIN(iNuM_Bar,nSpecies)
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) .AND. MatterSource_Mask(iYm,iN_X) )THEN
        iOS = ( (iN_E-1) + (iS-1)*nE_G ) * nCR

        StateMoment = Dnu(iN_E,iS,iN_X)
        Gm(iOS+iCR_N,iN_X) = StateMoment &
          + Damping(iN_X) * Fm(iOS+iCR_N,iN_X)
        Fm(iOS+iCR_N,iN_X) = Damping(iN_X) * Fm(iOS+iCR_N,iN_X)

        StateMoment = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)
        Gm(iOS+iCR_G1,iN_X) = StateMoment &
          + Damping(iN_X) * Fm(iOS+iCR_G1,iN_X)
        Fm(iOS+iCR_G1,iN_X) = Damping(iN_X) * Fm(iOS+iCR_G1,iN_X)

        StateMoment = Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)
        Gm(iOS+iCR_G2,iN_X) = StateMoment &
          + Damping(iN_X) * Fm(iOS+iCR_G2,iN_X)
        Fm(iOS+iCR_G2,iN_X) = Damping(iN_X) * Fm(iOS+iCR_G2,iN_X)

        StateMoment = Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)
        Gm(iOS+iCR_G3,iN_X) = StateMoment &
          + Damping(iN_X) * Fm(iOS+iCR_G3,iN_X)
        Fm(iOS+iCR_G3,iN_X) = Damping(iN_X) * Fm(iOS+iCR_G3,iN_X)
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyRadiationOnlySchurNewtonPreconditioner_OrderV


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


  SUBROUTINE ComputeLeptonFlavorCheck_OrderV &
    ( MASK, dt, ApplyMapProjection, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, CheckRate )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP),                   INTENT(in)    :: dt
    LOGICAL,                    INTENT(in)    :: ApplyMapProjection
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: CheckRate

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

      CheckRate(:,iLFC_CollisionSource:nLeptonFlavorCheck,iN_X) = Zero
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
              + Chi_Brem(iN_E,iS,iN_X) &
              + Chi_GL(iN_E,iS,iN_X,iGL_NMS) &
              + Chi_GL(iN_E,iS,iN_X,iGL_LFE) &
              + Chi_GL(iN_E,iS,iN_X,iGL_LFC) &
              + Chi_GL(iN_E,iS,iN_X,iGL_IMD)
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

        CheckRate(iFamily,iLFC_CollisionSource,iN_X) &
          = CheckRate(iFamily,iLFC_CollisionSource,iN_X) &
              + SignedWeight * CollisionRate
        CheckRate(iFamily,iLFC_RawMap,iN_X) &
          = CheckRate(iFamily,iLFC_RawMap,iN_X) &
              + SignedWeight * ( N_Map - N_Old ) / dt
        CheckRate(iFamily,iLFC_MapProjection,iN_X) &
          = CheckRate(iFamily,iLFC_MapProjection,iN_X) &
              + SignedWeight * ( N_Projected - N_Map ) / dt
        CheckRate(iFamily,iLFC_FixedPointResidual,iN_X) &
          = CheckRate(iFamily,iLFC_FixedPointResidual,iN_X) &
              + SignedWeight * ( N_Projected - N_Current ) / dt
        CheckRate(iFamily,iLFC_EquationResidual,iN_X) &
          = CheckRate(iFamily,iLFC_EquationResidual,iN_X) &
              + SignedWeight &
                  * ( ( N_Current - N_Old ) / dt - CollisionRate )

      END DO
      END DO

    END DO

  END SUBROUTINE ComputeLeptonFlavorCheck_OrderV


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
      D, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      RadiationLimiterCorrection, Error, &
      RetainMuonFraction_Option, MuonFractionConserved_Option )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      Ye, Ym, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: RadiationLimiterCorrection
    INTEGER,  DIMENSION(:),     INTENT(inout) :: Error
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: RetainMuonFraction_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: MuonFractionConserved_Option

    INTEGER  :: iN_E, iN_X, iS, SolveInfo
    REAL(DP) :: A(3,3), B(3), VHat(3)
    REAL(DP) :: Ef, SqrtG1, SqrtG2, SqrtG3, vDotV
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13
    REAL(DP) :: k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: SUM_D
    REAL(DP) :: SUM_K11, SUM_K12, SUM_K13
    REAL(DP) :: SUM_K22, SUM_K23, SUM_K33
    REAL(DP) :: SUM_DeltaLe, SUM_DeltaLm
    REAL(DP) :: SUM_DeltaEf
    REAL(DP) :: SUM_DeltaP1, SUM_DeltaP2, SUM_DeltaP3
    REAL(DP) :: DeltaEf, DeltaLe, DeltaLm, DeltaN
    REAL(DP) :: U_Ye_New, U_Ym_New, U_Ef_New
    REAL(DP) :: U_V1_New, U_V2_New, U_V3_New
    REAL(DP) :: wYe, wYm, wEf, wV1, wV2, wV3

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( A, B, VHat, SolveInfo, Ef, SqrtG1, SqrtG2, SqrtG3, &
    !$OMP          vDotV, SUM_D, &
    !$OMP          SUM_K11, SUM_K12, SUM_K13, &
    !$OMP          SUM_K22, SUM_K23, SUM_K33, &
    !$OMP          SUM_DeltaLe, SUM_DeltaLm, &
    !$OMP          SUM_DeltaEf, &
    !$OMP          SUM_DeltaP1, SUM_DeltaP2, SUM_DeltaP3, &
    !$OMP          DeltaEf, DeltaLe, DeltaLm, DeltaN, &
    !$OMP          U_Ye_New, U_Ym_New, U_Ef_New, &
    !$OMP          U_V1_New, U_V2_New, U_V3_New, &
    !$OMP          wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( A, B, VHat, SolveInfo, Ef, SqrtG1, SqrtG2, SqrtG3, &
    !$ACC          vDotV, SUM_D, &
    !$ACC          SUM_K11, SUM_K12, SUM_K13, &
    !$ACC          SUM_K22, SUM_K23, SUM_K33, &
    !$ACC          SUM_DeltaLe, SUM_DeltaLm, &
    !$ACC          SUM_DeltaEf, &
    !$ACC          SUM_DeltaP1, SUM_DeltaP2, SUM_DeltaP3, &
    !$ACC          DeltaEf, DeltaLe, DeltaLm, DeltaN, &
    !$ACC          U_Ye_New, U_Ym_New, U_Ef_New, &
    !$ACC          U_V1_New, U_V2_New, U_V3_New, &
    !$ACC          wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( A, B, VHat, SolveInfo, Ef, SqrtG1, SqrtG2, SqrtG3, &
    !$OMP          vDotV, k_dd_11, k_dd_12, k_dd_13, &
    !$OMP          k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          SUM_D, &
    !$OMP          SUM_K11, SUM_K12, SUM_K13, &
    !$OMP          SUM_K22, SUM_K23, SUM_K33, &
    !$OMP          SUM_DeltaLe, SUM_DeltaLm, &
    !$OMP          SUM_DeltaEf, &
    !$OMP          SUM_DeltaP1, SUM_DeltaP2, SUM_DeltaP3, &
    !$OMP          DeltaEf, DeltaLe, DeltaLm, DeltaN, &
    !$OMP          U_Ye_New, U_Ym_New, U_Ef_New, &
    !$OMP          U_V1_New, U_V2_New, U_V3_New, &
    !$OMP          wYe, wYm, wEf, wV1, wV2, wV3 )
#endif
    DO iN_X = 1, nX_G

      IF( .NOT. MASK(iN_X) ) CYCLE

      Error(iN_X) = 0
      INFO(iN_X) = 0

      wYe = MatterSourceWeight( iYe, iN_X )
      wYm = MatterSourceWeight( iYm, iN_X )
      wEf = MatterSourceWeight( iEf, iN_X )
      wV1 = MatterSourceWeight( iV1, iN_X )
      wV2 = MatterSourceWeight( iV2, iN_X )
      wV3 = MatterSourceWeight( iV3, iN_X )

      SqrtG1 = SQRT( Gm_dd_11(iN_X) )
      SqrtG2 = SQRT( Gm_dd_22(iN_X) )
      SqrtG3 = SQRT( Gm_dd_33(iN_X) )

      SUM_D   = Zero
      SUM_K11 = Zero
      SUM_K12 = Zero
      SUM_K13 = Zero
      SUM_K22 = Zero
      SUM_K23 = Zero
      SUM_K33 = Zero
      SUM_DeltaP1 = Zero
      SUM_DeltaP2 = Zero
      SUM_DeltaP3 = Zero

#if   defined(THORNADO_OMP_OL)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, &
    !$OMP          k_dd_22, k_dd_23, k_dd_33 ) &
    !$OMP REDUCTION( + : SUM_D, &
    !$OMP                 SUM_K11, SUM_K12, SUM_K13, &
    !$OMP                 SUM_K22, SUM_K23, SUM_K33, &
    !$OMP                 SUM_DeltaP1, SUM_DeltaP2, SUM_DeltaP3 )
#elif defined(THORNADO_OACC)
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, &
    !$ACC          k_dd_22, k_dd_23, k_dd_33 ) &
    !$ACC REDUCTION( + : SUM_D, &
    !$ACC                 SUM_K11, SUM_K12, SUM_K13, &
    !$ACC                 SUM_K22, SUM_K23, SUM_K33, &
    !$ACC                 SUM_DeltaP1, SUM_DeltaP2, SUM_DeltaP3 )
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

        SUM_DeltaP1 = SUM_DeltaP1 + W3_S(iN_E) &
          * ( C_Inu_d_1(iN_E,iS,iN_X) / SqrtG1 &
              + SqrtG1 * V_u_1_old(iN_X) * Dnu_old(iN_E,iS,iN_X) &
              - SqrtG1 * Inu_u_1(iN_E,iS,iN_X) )
        SUM_DeltaP2 = SUM_DeltaP2 + W3_S(iN_E) &
          * ( C_Inu_d_2(iN_E,iS,iN_X) / SqrtG2 &
              + SqrtG2 * V_u_2_old(iN_X) * Dnu_old(iN_E,iS,iN_X) &
              - SqrtG2 * Inu_u_2(iN_E,iS,iN_X) )
        SUM_DeltaP3 = SUM_DeltaP3 + W3_S(iN_E) &
          * ( C_Inu_d_3(iN_E,iS,iN_X) / SqrtG3 &
              + SqrtG3 * V_u_3_old(iN_X) * Dnu_old(iN_E,iS,iN_X) &
              - SqrtG3 * Inu_u_3(iN_E,iS,iN_X) )

      END DO
      END DO

      A = Zero
      A(1,1) = One + wV1 &
        * ( SUM_D + SUM_K11 ) / D(iN_X)
      A(1,2) = wV1 * SUM_K12 / D(iN_X)
      A(1,3) = wV1 * SUM_K13 / D(iN_X)
      A(2,1) = wV2 * SUM_K12 / D(iN_X)
      A(2,2) = One + wV2 &
        * ( SUM_D + SUM_K22 ) / D(iN_X)
      A(2,3) = wV2 * SUM_K23 / D(iN_X)
      A(3,1) = wV3 * SUM_K13 / D(iN_X)
      A(3,2) = wV3 * SUM_K23 / D(iN_X)
      A(3,3) = One + wV3 &
        * ( SUM_D + SUM_K33 ) / D(iN_X)

      B(1) = SqrtG1 * V_u_1_old(iN_X) &
        + wV1 / D(iN_X) &
          * ( SUM_DeltaP1 &
              + RadiationLimiterCorrection(iV1,iN_X) / SqrtG1 )
      B(2) = SqrtG2 * V_u_2_old(iN_X) &
        + wV2 / D(iN_X) &
          * ( SUM_DeltaP2 &
              + RadiationLimiterCorrection(iV2,iN_X) / SqrtG2 )
      B(3) = SqrtG3 * V_u_3_old(iN_X) &
        + wV3 / D(iN_X) &
          * ( SUM_DeltaP3 &
              + RadiationLimiterCorrection(iV3,iN_X) / SqrtG3 )

      CALL SolveLinearSystem3x3_OrderV( A, B, VHat, SolveInfo )
      INFO(iN_X) = SolveInfo

      IF( SolveInfo /= 0 )THEN
        Error(iN_X) = 1
        Fm(:,iN_X) = Zero
        Gm(iYe,iN_X) = U_Ye(iN_X)
        Gm(iYm,iN_X) = U_Ym(iN_X)
        Gm(iEf,iN_X) = U_Ef   (iN_X)
        Gm(iV1,iN_X) = U_V_d_1(iN_X)
        Gm(iV2,iN_X) = U_V_d_2(iN_X)
        Gm(iV3,iN_X) = U_V_d_3(iN_X)
        CYCLE
      END IF

      V_u_1(iN_X) = VHat(1) / SqrtG1
      V_u_2(iN_X) = VHat(2) / SqrtG2
      V_u_3(iN_X) = VHat(3) / SqrtG3

      ! Form the accepted old-to-new Eulerian number and energy changes
      ! before the spectral reduction.  This avoids subtracting separately
      ! accumulated inventories, including their order-v flux contributions.

      SUM_DeltaLe = Zero
      SUM_DeltaLm = Zero
      SUM_DeltaEf = Zero

#if   defined(THORNADO_OMP_OL)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) PRIVATE( DeltaEf, DeltaN ) &
    !$OMP REDUCTION( + : SUM_DeltaLe, SUM_DeltaLm, SUM_DeltaEf )
#elif defined(THORNADO_OACC)
    !$ACC LOOP VECTOR COLLAPSE(2) PRIVATE( DeltaEf, DeltaN ) &
    !$ACC REDUCTION( + : SUM_DeltaLe, SUM_DeltaLm, SUM_DeltaEf )
#endif
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        DeltaN = C_Dnu(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) &
          - VHat(1) * SqrtG1 * Inu_u_1(iN_E,iS,iN_X) &
          - VHat(2) * SqrtG2 * Inu_u_2(iN_E,iS,iN_X) &
          - VHat(3) * SqrtG3 * Inu_u_3(iN_E,iS,iN_X)

        DeltaEf = Two * C_Dnu(iN_E,iS,iN_X) &
          - Dnu_old(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) &
          - Two * VHat(1) * SqrtG1 * Inu_u_1(iN_E,iS,iN_X) &
          - Two * VHat(2) * SqrtG2 * Inu_u_2(iN_E,iS,iN_X) &
          - Two * VHat(3) * SqrtG3 * Inu_u_3(iN_E,iS,iN_X)

        SUM_DeltaEf = SUM_DeltaEf + W3_S(iN_E) * DeltaEf

        IF( iS <= iNuE_Bar )THEN
          SUM_DeltaLe = SUM_DeltaLe &
            + W2_S(iN_E) * LeptonNumber(iS) * DeltaN
        END IF

        IF( iS == iNuM .OR. iS == iNuM_Bar )THEN
          SUM_DeltaLm = SUM_DeltaLm &
            + W2_S(iN_E) * LeptonNumber(iS) * DeltaN
        END IF

      END DO
      END DO

      DeltaLe = SUM_DeltaLe + RadiationLimiterCorrection(iYe,iN_X)
      DeltaLm = SUM_DeltaLm + RadiationLimiterCorrection(iYm,iN_X)

      Ye(iN_X) = Ye_old(iN_X) + AtomicMassUnit / D(iN_X) &
        * wYe * DeltaLe
      U_Ym_New = Ym_old(iN_X) + AtomicMassUnit / D(iN_X) * wYm * DeltaLm
      IF( PRESENT(MuonFractionConserved_Option) ) &
        MuonFractionConserved_Option(iN_X) = U_Ym_New
      IF( PRESENT(RetainMuonFraction_Option) )THEN
        IF( .NOT. RetainMuonFraction_Option(iN_X) ) Ym(iN_X) = U_Ym_New
      ELSE
        Ym(iN_X) = U_Ym_New
      END IF
      Ef = Ef_old(iN_X) + wEf / D(iN_X) &
        * ( SUM_DeltaEf + RadiationLimiterCorrection(iEf,iN_X) )

      U_Ye_New = Ye(iN_X) / Ye_old(iN_X)
      U_Ym_New = Ym(iN_X)
      U_Ef_New = Ef / Ef_old(iN_X)
      U_V1_New = SqrtG1 * VHat(1) / SpeedOfLight
      U_V2_New = SqrtG2 * VHat(2) / SpeedOfLight
      U_V3_New = SqrtG3 * VHat(3) / SpeedOfLight

      Fm(iYe,iN_X) = U_Ye_New - U_Ye(iN_X)
      Fm(iYm,iN_X) = U_Ym_New - U_Ym(iN_X)
      Fm(iEf,iN_X) = U_Ef_New - U_Ef   (iN_X)
      Fm(iV1,iN_X) = U_V1_New - U_V_d_1(iN_X)
      Fm(iV2,iN_X) = U_V2_New - U_V_d_2(iN_X)
      Fm(iV3,iN_X) = U_V3_New - U_V_d_3(iN_X)

      U_Ye(iN_X) = U_Ye_New
      U_Ym(iN_X) = U_Ym_New
      U_Ef   (iN_X) = U_Ef_New
      U_V_d_1(iN_X) = U_V1_New
      U_V_d_2(iN_X) = U_V2_New
      U_V_d_3(iN_X) = U_V3_New

      G_Ye(iN_X) = U_Ye_New
      G_Ym(iN_X) = U_Ym_New
      G_Ef   (iN_X) = U_Ef_New
      G_V_d_1(iN_X) = U_V1_New
      G_V_d_2(iN_X) = U_V2_New
      G_V_d_3(iN_X) = U_V3_New

      Gm(iYe,iN_X) = U_Ye_New
      Gm(iYm,iN_X) = U_Ym_New
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
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
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
      T, Ye, Ym, E, V_u_1, V_u_2, V_u_3
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
               D, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
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
           ( D, E, Ye, Ym, T, MASK, nX_P, PackIndex, UnpackIndex, Error, &
             DisplayErrors_Option = .FALSE. )

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


  SUBROUTINE EvaluateRadiationOnlyMapForJacobian_OrderV &
    ( dt, iPoint, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
      FVECm_outer, GVECm_outer, FVECm, GVECm, &
      MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
      RadiationLimiterCorrection, RadiationLeptonInvariant, &
      MapDensity, EtaMu, EvaluationFailed, RecoverMatter_Option )

    REAL(DP),                   INTENT(in)    :: dt
    INTEGER,                    INTENT(in)    :: iPoint
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      T, Ye, Ym, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: &
      FVECm_outer, GVECm_outer, FVECm, GVECm
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: MASK, WorkingMask
    INTEGER,  DIMENSION(:),     INTENT(inout) :: &
      Error, PackIndex, UnpackIndex
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: &
      RadiationLimiterCorrection, RadiationLeptonInvariant
    REAL(DP), DIMENSION(:),     INTENT(out)   :: MapDensity
    REAL(DP),                   INTENT(out)   :: EtaMu
    LOGICAL,                    INTENT(out)   :: EvaluationFailed
    LOGICAL,                    INTENT(in), OPTIONAL :: RecoverMatter_Option

    INTEGER  :: iN_E, iOS, iS, nX_P
    LOGICAL  :: RecoverMatter
    REAL(DP) :: Mumu

    MASK = .FALSE.
    MASK(iPoint) = .TRUE.
    WorkingMask = .FALSE.
    Error = 0

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP ( MASK, WorkingMask, Error, &
    !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   T, Ye, Ym, E, V_u_1, V_u_2, V_u_3 )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( MASK, WorkingMask, Error, &
    !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   T, Ye, Ym, E, V_u_1, V_u_2, V_u_3 )
#endif

    CALL CreatePackIndex &
           ( MASK, nX_P, PackIndex, UnpackIndex )

    RecoverMatter = .TRUE.
    IF( PRESENT(RecoverMatter_Option) ) RecoverMatter = RecoverMatter_Option

    IF( RecoverMatter ) &
      CALL RecoverMatterFromRadiationInvariants_OrderV &
             ( MASK, WorkingMask, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               RadiationLimiterCorrection, RadiationLeptonInvariant, &
               nX_P, PackIndex, UnpackIndex, Error )

    EvaluationFailed = Error(iPoint) /= 0
    IF( EvaluationFailed )THEN
      MapDensity = Zero
      EtaMu = Zero
      RETURN
    END IF

    IF( .NOT. FreezeOpacities ) &
      CALL ComputeOpacities_Packed &
             ( D, T, Ye, Ym, SqrtGm, MASK, nX_P, &
               PackIndex, UnpackIndex, &
               T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, &
               Opacity_Mask_Input = Opacity_Mask )

    CALL ComputeDnuNorm( MASK, Dnu )

    IF( FreezeOpacities )THEN
      CALL ComputeRates_Packed &
             ( D, T_Opacity, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               MASK, nX_P, PackIndex, UnpackIndex, &
               nX_P, Opacity_Mask )
    ELSE
      CALL ComputeRates_Packed &
             ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               MASK, nX_P, PackIndex, UnpackIndex, &
               nX_P, Opacity_Mask )
    END IF

    CALL ComputeNeutrinoRHS_OrderV &
           ( MASK, FVECm, GVECm, dt, &
             Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
             V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    CALL ProjectRadiationOnlyFixedPointMap_OrderV &
           ( MASK, FVECm, GVECm, &
             Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( GVECm, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, Error )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( GVECm, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, Error )
#endif

    EvaluationFailed = Error(iPoint) /= 0
    MapDensity = Zero
    DO iS = iNuM, iNuM_Bar
    DO iN_E = 1, nE_G
      iOS = ( (iN_E-1) + (iS-1)*nE_G ) * nCR
      MapDensity(iN_E+(iS-iNuM)*nE_G) = GVECm(iOS+iCR_N,iPoint)
    END DO
    END DO

    CALL ComputeMuonChemicalPotential_TABLE &
           ( D(iPoint), T(iPoint), Ye(iPoint), Ym(iPoint), Mumu )
    EtaMu = Mumu / MAX( T(iPoint), SqrtTiny )

  END SUBROUTINE EvaluateRadiationOnlyMapForJacobian_OrderV


  SUBROUTINE ApplyRadiationOnlyDensityJacobian_OrderV &
    ( dt, Step, iPoint, Direction, Scale, &
      Dnu_Base, Inu_u_1_Base, Inu_u_2_Base, Inu_u_3_Base, &
      Matter_Base, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
      FVECm_outer, GVECm_outer, FVECm, GVECm, &
      MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
      RadiationLimiterCorrection, RadiationLeptonInvariant, &
      JacobianDirection, dLInput, dLMap, dEtaMu, dYm, &
      dTemperature, dEInput, &
      EvaluationFailed, JacobianDirectionFixedMatter_Option, &
      dLMapFixedMatter_Option )

    REAL(DP),                   INTENT(in)    :: dt, Step
    INTEGER,                    INTENT(in)    :: iPoint
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Direction, Scale
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: &
      Dnu_Base, Inu_u_1_Base, Inu_u_2_Base, Inu_u_3_Base
    REAL(DP), DIMENSION(7),     INTENT(in)    :: Matter_Base
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      T, Ye, Ym, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: &
      FVECm_outer, GVECm_outer, FVECm, GVECm
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: MASK, WorkingMask
    INTEGER,  DIMENSION(:),     INTENT(inout) :: &
      Error, PackIndex, UnpackIndex
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: &
      RadiationLimiterCorrection, RadiationLeptonInvariant
    REAL(DP), DIMENSION(:),     INTENT(out)   :: JacobianDirection
    REAL(DP),                   INTENT(out)   :: &
      dLInput, dLMap, dEtaMu, dYm, dTemperature, dEInput
    LOGICAL,                    INTENT(out)   :: EvaluationFailed
    REAL(DP), DIMENSION(:),     INTENT(out), OPTIONAL :: &
      JacobianDirectionFixedMatter_Option
    REAL(DP),                   INTENT(out), OPTIONAL :: &
      dLMapFixedMatter_Option

    INTEGER  :: iN_E, iS, iV
    LOGICAL  :: ComputeFixedMatter
    REAL(DP) :: EtaFixed, EtaMinus, EtaPlus
    REAL(DP) :: TMinus, TPlus, YmMinus, YmPlus
    REAL(DP), ALLOCATABLE :: MapMinus(:), MapPlus(:)

    ALLOCATE( MapMinus(SIZE(Direction)), MapPlus(SIZE(Direction)) )

    DO iS = 1, nSpecies
    DO iN_E = 1, nE_G
      Dnu    (iN_E,iS,iPoint) = Dnu_Base    (iN_E,iS)
      Inu_u_1(iN_E,iS,iPoint) = Inu_u_1_Base(iN_E,iS)
      Inu_u_2(iN_E,iS,iPoint) = Inu_u_2_Base(iN_E,iS)
      Inu_u_3(iN_E,iS,iPoint) = Inu_u_3_Base(iN_E,iS)
    END DO
    END DO
    T(iPoint) = Matter_Base(1)
    Ye(iPoint) = Matter_Base(2)
    Ym(iPoint) = Matter_Base(3)
    E(iPoint) = Matter_Base(4)
    V_u_1(iPoint) = Matter_Base(5)
    V_u_2(iPoint) = Matter_Base(6)
    V_u_3(iPoint) = Matter_Base(7)
    DO iS = iNuM, iNuM_Bar
    DO iN_E = 1, nE_G
      iV = iN_E + (iS-iNuM)*nE_G
      Dnu(iN_E,iS,iPoint) = Dnu(iN_E,iS,iPoint) &
        + Step * Scale(iV) * Direction(iV)
    END DO
    END DO

    CALL EvaluateRadiationOnlyMapForJacobian_OrderV &
           ( dt, iPoint, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
             D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
             FVECm_outer, GVECm_outer, FVECm, GVECm, &
             MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
             RadiationLimiterCorrection, RadiationLeptonInvariant, &
             MapPlus, EtaPlus, EvaluationFailed )
    IF( EvaluationFailed )THEN
      DEALLOCATE( MapMinus, MapPlus )
      RETURN
    END IF
    TPlus = T(iPoint)
    YmPlus = Ym(iPoint)

    DO iS = 1, nSpecies
    DO iN_E = 1, nE_G
      Dnu    (iN_E,iS,iPoint) = Dnu_Base    (iN_E,iS)
      Inu_u_1(iN_E,iS,iPoint) = Inu_u_1_Base(iN_E,iS)
      Inu_u_2(iN_E,iS,iPoint) = Inu_u_2_Base(iN_E,iS)
      Inu_u_3(iN_E,iS,iPoint) = Inu_u_3_Base(iN_E,iS)
    END DO
    END DO
    T(iPoint) = Matter_Base(1)
    Ye(iPoint) = Matter_Base(2)
    Ym(iPoint) = Matter_Base(3)
    E(iPoint) = Matter_Base(4)
    V_u_1(iPoint) = Matter_Base(5)
    V_u_2(iPoint) = Matter_Base(6)
    V_u_3(iPoint) = Matter_Base(7)
    DO iS = iNuM, iNuM_Bar
    DO iN_E = 1, nE_G
      iV = iN_E + (iS-iNuM)*nE_G
      Dnu(iN_E,iS,iPoint) = Dnu(iN_E,iS,iPoint) &
        - Step * Scale(iV) * Direction(iV)
    END DO
    END DO

    CALL EvaluateRadiationOnlyMapForJacobian_OrderV &
           ( dt, iPoint, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
             D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
             FVECm_outer, GVECm_outer, FVECm, GVECm, &
             MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
             RadiationLimiterCorrection, RadiationLeptonInvariant, &
             MapMinus, EtaMinus, EvaluationFailed )
    IF( EvaluationFailed )THEN
      DEALLOCATE( MapMinus, MapPlus )
      RETURN
    END IF
    TMinus = T(iPoint)
    YmMinus = Ym(iPoint)

    JacobianDirection = Zero
    dLInput = Zero
    dLMap = Zero
    dEInput = Zero
    DO iS = iNuM, iNuM_Bar
    DO iN_E = 1, nE_G
      iV = iN_E + (iS-iNuM)*nE_G
      IF( Scale(iV) > Zero ) &
        JacobianDirection(iV) &
          = ( MapPlus(iV) - MapMinus(iV) ) &
              / ( Two * Step * Scale(iV) )
      dLInput = dLInput + W2_S(iN_E) * LeptonNumber(iS) &
        * Scale(iV) * Direction(iV)
      dLMap = dLMap + W2_S(iN_E) * LeptonNumber(iS) &
        * ( MapPlus(iV) - MapMinus(iV) ) / ( Two * Step )
      dEInput = dEInput + W3_S(iN_E) &
        * Scale(iV) * Direction(iV)
    END DO
    END DO
    dEtaMu = ( EtaPlus - EtaMinus ) / ( Two * Step )
    dYm = ( YmPlus - YmMinus ) / ( Two * Step )
    dTemperature = ( TPlus - TMinus ) / ( Two * Step )

    ComputeFixedMatter &
      = PRESENT(JacobianDirectionFixedMatter_Option) &
          .AND. PRESENT(dLMapFixedMatter_Option)
    IF( PRESENT(JacobianDirectionFixedMatter_Option) .NEQV. &
        PRESENT(dLMapFixedMatter_Option) ) &
      ERROR STOP 'Incomplete fixed-matter Jacobian output request.'

    IF( ComputeFixedMatter )THEN

      DO iS = 1, nSpecies
      DO iN_E = 1, nE_G
        Dnu    (iN_E,iS,iPoint) = Dnu_Base    (iN_E,iS)
        Inu_u_1(iN_E,iS,iPoint) = Inu_u_1_Base(iN_E,iS)
        Inu_u_2(iN_E,iS,iPoint) = Inu_u_2_Base(iN_E,iS)
        Inu_u_3(iN_E,iS,iPoint) = Inu_u_3_Base(iN_E,iS)
      END DO
      END DO
      T(iPoint) = Matter_Base(1)
      Ye(iPoint) = Matter_Base(2)
      Ym(iPoint) = Matter_Base(3)
      E(iPoint) = Matter_Base(4)
      V_u_1(iPoint) = Matter_Base(5)
      V_u_2(iPoint) = Matter_Base(6)
      V_u_3(iPoint) = Matter_Base(7)
      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        iV = iN_E + (iS-iNuM)*nE_G
        Dnu(iN_E,iS,iPoint) = Dnu(iN_E,iS,iPoint) &
          + Step * Scale(iV) * Direction(iV)
      END DO
      END DO

      CALL EvaluateRadiationOnlyMapForJacobian_OrderV &
             ( dt, iPoint, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
               FVECm_outer, GVECm_outer, FVECm, GVECm, &
               MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
               RadiationLimiterCorrection, RadiationLeptonInvariant, &
               MapPlus, EtaFixed, EvaluationFailed, &
               RecoverMatter_Option = .FALSE. )
      IF( EvaluationFailed )THEN
        DEALLOCATE( MapMinus, MapPlus )
        RETURN
      END IF

      DO iS = 1, nSpecies
      DO iN_E = 1, nE_G
        Dnu    (iN_E,iS,iPoint) = Dnu_Base    (iN_E,iS)
        Inu_u_1(iN_E,iS,iPoint) = Inu_u_1_Base(iN_E,iS)
        Inu_u_2(iN_E,iS,iPoint) = Inu_u_2_Base(iN_E,iS)
        Inu_u_3(iN_E,iS,iPoint) = Inu_u_3_Base(iN_E,iS)
      END DO
      END DO
      T(iPoint) = Matter_Base(1)
      Ye(iPoint) = Matter_Base(2)
      Ym(iPoint) = Matter_Base(3)
      E(iPoint) = Matter_Base(4)
      V_u_1(iPoint) = Matter_Base(5)
      V_u_2(iPoint) = Matter_Base(6)
      V_u_3(iPoint) = Matter_Base(7)
      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        iV = iN_E + (iS-iNuM)*nE_G
        Dnu(iN_E,iS,iPoint) = Dnu(iN_E,iS,iPoint) &
          - Step * Scale(iV) * Direction(iV)
      END DO
      END DO

      CALL EvaluateRadiationOnlyMapForJacobian_OrderV &
             ( dt, iPoint, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
               FVECm_outer, GVECm_outer, FVECm, GVECm, &
               MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
               RadiationLimiterCorrection, RadiationLeptonInvariant, &
               MapMinus, EtaFixed, EvaluationFailed, &
               RecoverMatter_Option = .FALSE. )
      IF( EvaluationFailed )THEN
        DEALLOCATE( MapMinus, MapPlus )
        RETURN
      END IF

      JacobianDirectionFixedMatter_Option = Zero
      dLMapFixedMatter_Option = Zero
      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        iV = iN_E + (iS-iNuM)*nE_G
        IF( Scale(iV) > Zero ) &
          JacobianDirectionFixedMatter_Option(iV) &
            = ( MapPlus(iV) - MapMinus(iV) ) &
                / ( Two * Step * Scale(iV) )
        dLMapFixedMatter_Option &
          = dLMapFixedMatter_Option &
              + W2_S(iN_E) * LeptonNumber(iS) &
                  * ( MapPlus(iV) - MapMinus(iV) ) / ( Two * Step )
      END DO
      END DO

    END IF

    DEALLOCATE( MapMinus, MapPlus )

  END SUBROUTINE ApplyRadiationOnlyDensityJacobian_OrderV


  SUBROUTINE AnalyzeRadiationOnlyJacobian_OrderV &
    ( dt, Iterations, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
      FVECm_outer, GVECm_outer, FVECm, GVECm, &
      MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
      RadiationLimiterCorrection, RadiationLeptonInvariant )

    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE mpi

    INTEGER, PARAMETER :: MaxAnalysisPoints = 8
    INTEGER, PARAMETER :: RequestedArnoldiDimension = 12
    REAL(DP), PARAMETER :: DifferenceStep = 1.0e-5_DP
    REAL(DP), PARAMETER :: OccupationFloor = 1.0e-12_DP

    REAL(DP),                   INTENT(in)    :: dt
    INTEGER,  DIMENSION(:),     INTENT(in)    :: Iterations
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: D
    REAL(DP), DIMENSION(:),     INTENT(inout) :: &
      T, Ye, Ym, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)    :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: &
      FVECm_outer, GVECm_outer, FVECm, GVECm
    LOGICAL,  DIMENSION(:),     INTENT(inout) :: MASK, WorkingMask
    INTEGER,  DIMENSION(:),     INTENT(inout) :: &
      Error, PackIndex, UnpackIndex
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: &
      RadiationLimiterCorrection, RadiationLeptonInvariant

    CHARACTER(640) :: FileName
    INTEGER :: AnalysisPoints(MaxAnalysisPoints)
    INTEGER :: BestIterations, BestPoint, i, iA, iB, iN_E, iPoint
    INTEGER :: iS, iV, ierr, iUnit, j, mEffective, myid, nActive, nMu
    LOGICAL :: EvaluationFailed
    REAL(DP) :: Beta, Correction, dEInput, dEInputHalf
    REAL(DP) :: dEtaMu, dEtaMuHalf, dEtaMudL, dEtaMudYm
    REAL(DP) :: dLMapdEta, dTemperature, dTemperatureHalf
    REAL(DP) :: dTemperaturedL, dYm, dYmHalf, dYmdL
    REAL(DP) :: dLInput, dLInputHalf, dLMap, dLMapHalf
    REAL(DP) :: dLMapFixedMatter, DirectGain, FeedbackGain
    REAL(DP) :: FDConsistency, MuChemicalPotential, NormDirection
    REAL(DP) :: DirectNorm, DirectRayleigh, FeedbackNorm, FeedbackRayleigh
    REAL(DP) :: MuLeptonGain, Rayleigh, SchurDerivative, UpperOccupation
    REAL(DP) :: Matter_Base(7)
    REAL(DP), ALLOCATABLE :: Dnu_Base(:,:), Inu_u_1_Base(:,:)
    REAL(DP), ALLOCATABLE :: Inu_u_2_Base(:,:), Inu_u_3_Base(:,:)
    REAL(DP), ALLOCATABLE :: Direction(:), DirectionFixedMatter(:)
    REAL(DP), ALLOCATABLE :: DirectionHalf(:), EnergyDirection(:)
    REAL(DP), ALLOCATABLE :: H(:,:), Q(:,:), Scale(:), Work(:)

    IF( nSpecies < iNuM_Bar ) RETURN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$OMP   Opacity_Mask, Dnu_0, Chi_EmAb, Iterations )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$ACC   Opacity_Mask, Dnu_0, Chi_EmAb, Iterations )
#endif

    AnalysisPoints = 0
    DO iA = 1, MaxAnalysisPoints
      BestPoint = 0
      BestIterations = -1
      DO iPoint = 1, nX_G
        IF( .NOT. Opacity_Mask(iOp_EmAb_Muon,iPoint) ) CYCLE
        IF( ANY(AnalysisPoints == iPoint) ) CYCLE
        IF( Iterations(iPoint) > BestIterations )THEN
          BestIterations = Iterations(iPoint)
          BestPoint = iPoint
        END IF
      END DO
      IF( BestPoint == 0 ) EXIT
      AnalysisPoints(iA) = BestPoint
    END DO

    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    WRITE(FileName,'(A,".jacobian.rank",I6.6,".dat")') &
      TRIM(RadiationOnlyTraceFile), myid
    OPEN( NEWUNIT = iUnit, FILE = TRIM(FileName), &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(iUnit,'(A)') '# Local Jacobian analysis of the reduced collision map'
    WRITE(iUnit,'(A)') &
      '# P: point metadata; S: scalar muon-lepton response; ' // &
      'H: Arnoldi matrix; Q: Arnoldi basis; E: equilibrium spectrum'
    WRITE(iUnit,'(A,ES26.16E3)') '# centered_difference_step: ', &
      DifferenceStep

    nMu = 2 * nE_G
    ALLOCATE( Dnu_Base(nE_G,nSpecies) )
    ALLOCATE( Inu_u_1_Base(nE_G,nSpecies) )
    ALLOCATE( Inu_u_2_Base(nE_G,nSpecies) )
    ALLOCATE( Inu_u_3_Base(nE_G,nSpecies) )
    ALLOCATE( Direction(nMu), DirectionFixedMatter(nMu) )
    ALLOCATE( DirectionHalf(nMu), EnergyDirection(nMu) )
    ALLOCATE( Scale(nMu), Work(nMu) )
    ALLOCATE( Q(nMu,RequestedArnoldiDimension+1) )
    ALLOCATE( H(RequestedArnoldiDimension+1,RequestedArnoldiDimension) )

    DO iA = 1, MaxAnalysisPoints

      iPoint = AnalysisPoints(iA)
      IF( iPoint == 0 ) EXIT

      Dnu_Base = Dnu(:,:,iPoint)
      Inu_u_1_Base = Inu_u_1(:,:,iPoint)
      Inu_u_2_Base = Inu_u_2(:,:,iPoint)
      Inu_u_3_Base = Inu_u_3(:,:,iPoint)
      Matter_Base = [ T(iPoint), Ye(iPoint), Ym(iPoint), E(iPoint), &
                      V_u_1(iPoint), V_u_2(iPoint), V_u_3(iPoint) ]

      UpperOccupation = MIN( DnuMax, One )
      Scale = Zero
      Direction = Zero
      EnergyDirection = Zero
      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        iV = iN_E + (iS-iNuM)*nE_G
        IF( Dnu_Base(iN_E,iS) > OccupationFloor * UpperOccupation &
              .AND. Dnu_Base(iN_E,iS) &
                < (One-OccupationFloor)*UpperOccupation )THEN
          Scale(iV) = Dnu_Base(iN_E,iS) &
            * ( One - Dnu_Base(iN_E,iS) / UpperOccupation )
          Direction(iV) = W2_S(iN_E) * LeptonNumber(iS) * Scale(iV)
          EnergyDirection(iV) = W3_S(iN_E) * Scale(iV)
        END IF
      END DO
      END DO

      nActive = COUNT( Scale > Zero )
      IF( nActive == 0 ) CYCLE

      Correction = DOT_PRODUCT(Direction,EnergyDirection) &
        / MAX( DOT_PRODUCT(EnergyDirection,EnergyDirection), SqrtTiny )
      Direction = Direction - Correction * EnergyDirection
      NormDirection = SQRT( DOT_PRODUCT(Direction,Direction) )
      IF( NormDirection <= SqrtTiny )THEN
        WRITE(ERROR_UNIT,'(A,I0)') &
          'Reduced-solver Jacobian analysis has no muon mode at point ', &
          iPoint
        CYCLE
      END IF
      Q = Zero
      H = Zero
      Q(:,1) = Direction / NormDirection

      CALL ApplyRadiationOnlyDensityJacobian_OrderV &
             ( dt, DifferenceStep, iPoint, Q(:,1), Scale, &
               Dnu_Base, Inu_u_1_Base, Inu_u_2_Base, Inu_u_3_Base, &
               Matter_Base, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
               FVECm_outer, GVECm_outer, FVECm, GVECm, &
               MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
               RadiationLimiterCorrection, RadiationLeptonInvariant, &
               Direction, dLInput, dLMap, dEtaMu, dYm, &
               dTemperature, dEInput, &
               EvaluationFailed )
      IF( EvaluationFailed ) CYCLE

      CALL ApplyRadiationOnlyDensityJacobian_OrderV &
             ( dt, Half*DifferenceStep, iPoint, Q(:,1), Scale, &
               Dnu_Base, Inu_u_1_Base, Inu_u_2_Base, Inu_u_3_Base, &
               Matter_Base, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
               FVECm_outer, GVECm_outer, FVECm, GVECm, &
               MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
               RadiationLimiterCorrection, RadiationLeptonInvariant, &
               DirectionHalf, dLInputHalf, dLMapHalf, dEtaMuHalf, &
               dYmHalf, dTemperatureHalf, dEInputHalf, EvaluationFailed, &
               JacobianDirectionFixedMatter_Option = DirectionFixedMatter, &
               dLMapFixedMatter_Option = dLMapFixedMatter )
      IF( EvaluationFailed ) CYCLE

      FDConsistency = SQRT &
        ( DOT_PRODUCT(Direction-DirectionHalf,Direction-DirectionHalf) ) &
          / MAX( SQRT(DOT_PRODUCT(DirectionHalf,DirectionHalf)), SqrtTiny )
      Rayleigh = DOT_PRODUCT( Q(:,1), DirectionHalf )
      MuLeptonGain = Zero
      DirectGain = Zero
      FeedbackGain = Zero
      dEtaMudL = Zero
      dYmdL = Zero
      dEtaMudYm = Zero
      dTemperaturedL = Zero
      dLMapdEta = Zero
      SchurDerivative = Zero
      IF( ABS(dLInputHalf) > SqrtTiny )THEN
        MuLeptonGain = dLMapHalf / dLInputHalf
        DirectGain = dLMapFixedMatter / dLInputHalf
        FeedbackGain = MuLeptonGain - DirectGain
        dEtaMudL = dEtaMuHalf / dLInputHalf
        dYmdL = dYmHalf / dLInputHalf
        dTemperaturedL = dTemperatureHalf / dLInputHalf
        SchurDerivative = One - MuLeptonGain
      END IF
      IF( ABS(dYmHalf) > SqrtTiny ) &
        dEtaMudYm = dEtaMuHalf / dYmHalf
      IF( ABS(dEtaMuHalf) > SqrtTiny ) &
        dLMapdEta = dLMapHalf / dEtaMuHalf
      DirectRayleigh = DOT_PRODUCT( Q(:,1), DirectionFixedMatter )
      FeedbackRayleigh = Rayleigh - DirectRayleigh
      DirectNorm = SQRT &
        ( DOT_PRODUCT(DirectionFixedMatter,DirectionFixedMatter) )
      DirectionFixedMatter = DirectionHalf - DirectionFixedMatter
      FeedbackNorm = SQRT &
        ( DOT_PRODUCT(DirectionFixedMatter,DirectionFixedMatter) )

      Work = DirectionHalf
      mEffective = RequestedArnoldiDimension
      DO j = 1, RequestedArnoldiDimension

        IF( j > 1 )THEN
          CALL ApplyRadiationOnlyDensityJacobian_OrderV &
                 ( dt, Half*DifferenceStep, iPoint, Q(:,j), Scale, &
                   Dnu_Base, Inu_u_1_Base, Inu_u_2_Base, &
                   Inu_u_3_Base, Matter_Base, &
                   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
                   FVECm_outer, GVECm_outer, FVECm, GVECm, &
                   MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
                   RadiationLimiterCorrection, RadiationLeptonInvariant, &
                   Work, dLInput, dLMap, dEtaMu, dYm, &
                   dTemperature, dEInput, &
                   EvaluationFailed )
          IF( EvaluationFailed )THEN
            mEffective = j - 1
            EXIT
          END IF
        END IF

        DO i = 1, j
          H(i,j) = DOT_PRODUCT( Q(:,i), Work )
          Work = Work - H(i,j) * Q(:,i)
        END DO
        DO i = 1, j
          Correction = DOT_PRODUCT( Q(:,i), Work )
          H(i,j) = H(i,j) + Correction
          Work = Work - Correction * Q(:,i)
        END DO

        Beta = SQRT( DOT_PRODUCT(Work,Work) )
        H(j+1,j) = Beta
        IF( j == RequestedArnoldiDimension ) EXIT
        IF( Beta <= 1.0e3_DP*EPSILON(One) )THEN
          mEffective = j
          EXIT
        END IF
        Q(:,j+1) = Work / Beta

      END DO

      Dnu(:,:,iPoint) = Dnu_Base
      Inu_u_1(:,:,iPoint) = Inu_u_1_Base
      Inu_u_2(:,:,iPoint) = Inu_u_2_Base
      Inu_u_3(:,:,iPoint) = Inu_u_3_Base
      T(iPoint) = Matter_Base(1)
      Ye(iPoint) = Matter_Base(2)
      Ym(iPoint) = Matter_Base(3)
      E(iPoint) = Matter_Base(4)
      V_u_1(iPoint) = Matter_Base(5)
      V_u_2(iPoint) = Matter_Base(6)
      V_u_3(iPoint) = Matter_Base(7)
      CALL EvaluateRadiationOnlyMapForJacobian_OrderV &
             ( dt, iPoint, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
               FVECm_outer, GVECm_outer, FVECm, GVECm, &
               MASK, WorkingMask, Error, PackIndex, UnpackIndex, &
               RadiationLimiterCorrection, RadiationLeptonInvariant, &
               Direction, dEtaMu, EvaluationFailed )
      IF( EvaluationFailed ) CYCLE

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( Dnu_0, Chi_EmAb )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( Dnu_0, Chi_EmAb )
#endif

      CALL ComputeMuonChemicalPotential_TABLE &
             ( D(iPoint), Matter_Base(1), Matter_Base(2), Matter_Base(3), &
               MuChemicalPotential )
      WRITE(iUnit,'(A,4(1X,I0),6(1X,ES26.16E3))') &
        'P', iPoint, Iterations(iPoint), nActive, mEffective, &
        D(iPoint)/Unit_D, Matter_Base(1)/MeV, Matter_Base(2), &
        Matter_Base(3), MuChemicalPotential/MeV, &
        MuChemicalPotential/MAX(Matter_Base(1),SqrtTiny)
      WRITE(iUnit,'(A,1X,I0,19(1X,ES26.16E3))') &
        'S', iPoint, dLInputHalf, dLMapHalf, &
        MuLeptonGain, DirectGain, FeedbackGain, &
        dEtaMudL, dYmdL, dEtaMudYm, dTemperaturedL, &
        dLMapdEta, SchurDerivative, &
        Rayleigh, DirectRayleigh, FeedbackRayleigh, &
        SQRT(DOT_PRODUCT(DirectionHalf,DirectionHalf)), &
        DirectNorm, FeedbackNorm, &
        FDConsistency, dEInputHalf

      DO j = 1, mEffective
      DO i = 1, mEffective + 1
        WRITE(iUnit,'(A,3(1X,I0),1X,ES26.16E3)') &
          'H', iPoint, i, j, H(i,j)
      END DO
      END DO
      DO j = 1, mEffective
      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        iV = iN_E + (iS-iNuM)*nE_G
        WRITE(iUnit,'(A,4(1X,I0),3(1X,ES26.16E3))') &
          'Q', iPoint, j, iS, iN_E, E_N(iN_E)/MeV, &
          Q(iV,j), Scale(iV)
      END DO
      END DO
      END DO

      DO iS = iNuM, iNuM_Bar
      DO iN_E = 1, nE_G
        WRITE(iUnit,'(A,3(1X,I0),5(1X,ES26.16E3))') &
          'E', iPoint, iS, iN_E, E_N(iN_E)/MeV, &
          Dnu_Base(iN_E,iS), Dnu_0(iN_E,iS,iPoint), &
          Chi_EmAb(iN_E,iS,iPoint), Scale(iN_E+(iS-iNuM)*nE_G)
      END DO
      END DO

      FLUSH(iUnit)

      Dnu(:,:,iPoint) = Dnu_Base
      Inu_u_1(:,:,iPoint) = Inu_u_1_Base
      Inu_u_2(:,:,iPoint) = Inu_u_2_Base
      Inu_u_3(:,:,iPoint) = Inu_u_3_Base
      T(iPoint) = Matter_Base(1)
      Ye(iPoint) = Matter_Base(2)
      Ym(iPoint) = Matter_Base(3)
      E(iPoint) = Matter_Base(4)
      V_u_1(iPoint) = Matter_Base(5)
      V_u_2(iPoint) = Matter_Base(6)
      V_u_3(iPoint) = Matter_Base(7)

    END DO

    CLOSE(iUnit)

    ! The work masks and error array belong to the completed nonlinear solve.
    ! Restore their successful post-solve state after the diagnostic map calls.
    MASK = .FALSE.
    WorkingMask = .FALSE.
    Error = 0

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$OMP   MASK, WorkingMask, Error )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$ACC   MASK, WorkingMask, Error )
#endif

    DEALLOCATE( Dnu_Base, Inu_u_1_Base, Inu_u_2_Base, Inu_u_3_Base )
    DEALLOCATE( Direction, DirectionFixedMatter )
    DEALLOCATE( DirectionHalf, EnergyDirection )
    DEALLOCATE( Scale, Work, Q, H )

  END SUBROUTINE AnalyzeRadiationOnlyJacobian_OrderV


  SUBROUTINE InitializeRadiationOnlyTrace_OrderV

    USE mpi

    INTEGER :: ierr, iM, myid
    CHARACTER(640) :: FileName, SpectrumFileName

    IF( RadiationOnlyTraceUnit /= 0 ) RETURN

    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    WRITE(FileName,'(A,".rank",I6.6,".dat")') &
      TRIM(RadiationOnlyTraceFile), myid
    WRITE(SpectrumFileName,'(A,".spectra.rank",I6.6,".dat")') &
      TRIM(RadiationOnlyTraceFile), myid
    ! InitializeCollisions/FinalizeCollisions bracket each implicit stage, so
    ! the trace units are closed between stage solves.  Replace the files for
    ! the first solve in the process and append all subsequent solves; using
    ! STATUS='REPLACE' unconditionally silently discarded the earlier stages.
    IF( RadiationOnlySolveCallCount == 1 )THEN
      OPEN( NEWUNIT = RadiationOnlyTraceUnit, FILE = TRIM(FileName), &
            STATUS = 'REPLACE', ACTION = 'WRITE' )
      OPEN( NEWUNIT = RadiationOnlyTraceSpectrumUnit, &
            FILE = TRIM(SpectrumFileName), &
            STATUS = 'REPLACE', ACTION = 'WRITE' )
    ELSE
      OPEN( NEWUNIT = RadiationOnlyTraceUnit, FILE = TRIM(FileName), &
            STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE' )
      OPEN( NEWUNIT = RadiationOnlyTraceSpectrumUnit, &
            FILE = TRIM(SpectrumFileName), &
            STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE' )
      RETURN
    END IF

    WRITE(RadiationOnlyTraceUnit,'(A)') &
      '# Invariant-constrained radiation solver iteration trace'
    WRITE(RadiationOnlyTraceUnit,'(A,I0)') '# n_species: ', nSpecies
    WRITE(RadiationOnlyTraceUnit,'(A,I0)') '# n_energy_nodes: ', nE_G
    WRITE(RadiationOnlyTraceUnit,'(A,I0)') '# anderson_depth: ', M_inner
    WRITE(RadiationOnlyTraceUnit,'(A)',ADVANCE='NO') &
      '# columns: call iteration point active_after converged error backtracks ' // &
      'dt omega_dt schur_gain schur_damping ' // &
      'rho T_MeV Ye Ym delta_T_rel delta_Ye delta_Ym E mu_mu_MeV eta_mu ' // &
      'mu_cc_active mu_total_lepton_invariant_error ' // &
      'mu_total_lepton_invariant_relative residual_relative ' // &
      'residual_density_relative residual_momentum_relative ' // &
      'convergence_measure residual_absolute residual_ratio ' // &
      'residual_change residual_cosine step_relative step2_relative ' // &
      'candidate_step projected_step projection_correction ' // &
      'state_max map_max Dnu_min Dnu_max mu_equilibrium_deviation ' // &
      'mu_chi_max mu_stiffness_max collision_stiffness_max ' // &
      'worst_species worst_energy worst_moment'
    DO iM = 1, nSpecies
      WRITE(RadiationOnlyTraceUnit,'(A,I0)',ADVANCE='NO') &
        ' residual_species_', iM
    END DO
    DO iM = 1, M_inner
      WRITE(RadiationOnlyTraceUnit,'(A,I0)',ADVANCE='NO') ' alpha_', iM
    END DO
    WRITE(RadiationOnlyTraceUnit,*)
    WRITE(RadiationOnlyTraceSpectrumUnit,'(A)') &
      '# Terminal spectra for the invariant-constrained radiation solve'
    WRITE(RadiationOnlyTraceSpectrumUnit,'(A)') &
      '# columns: call iteration point converged species energy_index ' // &
      'energy_MeV D D0 I1 I2 I3 F_N F_G1 F_G2 F_G3 ' // &
      'G_N G_G1 G_G2 G_G3 Chi_EmAb Sigma_Iso Chi_NNS Chi_NES ' // &
      'Chi_Pair Chi_NuPair Chi_Brem'
    FLUSH(RadiationOnlyTraceUnit)
    FLUSH(RadiationOnlyTraceSpectrumUnit)

  END SUBROUTINE InitializeRadiationOnlyTrace_OrderV


  SUBROUTINE PackRadiationStateForTrace_OrderV &
    ( iN_X, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, State )

    INTEGER, INTENT(in) :: iN_X
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:), INTENT(out) :: State

    INTEGER :: iN_E, iOS, iS

    State = Zero

    DO iS = 1, nSpecies
    DO iN_E = 1, nE_G
      iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR
      State(iOS+iCR_N ) = Dnu(iN_E,iS,iN_X)
      State(iOS+iCR_G1) = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)
      State(iOS+iCR_G2) = Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)
      State(iOS+iCR_G3) = Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)
    END DO
    END DO

  END SUBROUTINE PackRadiationStateForTrace_OrderV


  SUBROUTINE WriteRadiationOnlyIterationTrace_OrderV &
    ( dt, k, EvaluationMask, ActiveMask, Error, D, T, Ye, Ym, E, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Fm, Gm, &
      RadiationLeptonInvariant, RadiationLimiterCorrection, &
      HasPrevious, HasPrevious2, &
      PreviousState, Previous2State, PreviousResidual, BacktrackCount, &
      CandidateStep, ProjectedStep, ProjectionCorrection, AlphaTrace, &
      SchurGain, SchurDamping )

    REAL(DP), INTENT(in) :: dt
    INTEGER, INTENT(in) :: k
    LOGICAL, DIMENSION(:), INTENT(in) :: EvaluationMask, ActiveMask
    INTEGER, DIMENSION(:), INTENT(in) :: Error
    REAL(DP), DIMENSION(:), INTENT(in) :: D, T, Ye, Ym, E
    REAL(DP), DIMENSION(:), INTENT(in) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:), INTENT(in) :: Fm, Gm
    REAL(DP), DIMENSION(:,:), INTENT(in) :: RadiationLeptonInvariant
    REAL(DP), DIMENSION(:,:), INTENT(in) :: RadiationLimiterCorrection
    LOGICAL, DIMENSION(:), INTENT(inout) :: HasPrevious, HasPrevious2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: &
      PreviousState, Previous2State, PreviousResidual
    INTEGER, DIMENSION(:), INTENT(inout) :: BacktrackCount
    REAL(DP), DIMENSION(:), INTENT(inout) :: &
      CandidateStep, ProjectedStep, ProjectionCorrection
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: AlphaTrace
    REAL(DP), DIMENSION(:), INTENT(in) :: SchurGain, SchurDamping

    INTEGER :: iChiral, iN_E, iN_X, iOS, iS, iV
    INTEGER :: WorstEnergy, WorstMoment, WorstSpecies, WorstVectorIndex
    INTEGER :: MuCCActive
    REAL(DP) :: ConvergenceMeasure, CurrentNorm, Denominator
    REAL(DP) :: DnuMinimum, DnuMaximum, DnuNormLocal
    REAL(DP) :: EtaMu
    REAL(DP) :: Fnorm, Fnorm_G1, Fnorm_G2, Fnorm_G3, Fnorm_N
    REAL(DP) :: InitialCriterion, MapMaximum, MuChiMaximum
    REAL(DP) :: MuEquilibriumDeviation, MuEquilibriumNorm
    REAL(DP) :: CollisionOpacity, CollisionStiffnessMaximum
    REAL(DP) :: MuStiffnessMaximum
    REAL(DP) :: MuLeptonCurrent, MuLeptonInvariantError
    REAL(DP) :: MuLeptonInvariantRelative, MuLeptonScale, Mumu
    REAL(DP) :: ResidualAbsolute, ResidualChange, ResidualCosine
    REAL(DP) :: ResidualDensityRelative, ResidualMomentumRelative
    REAL(DP) :: ResidualDot, ResidualNorm2, ResidualPreviousNorm2
    REAL(DP) :: ResidualRatio, ResidualRelative
    REAL(DP) :: StepRelative, Step2Relative, vDotI
    REAL(DP), ALLOCATABLE :: CurrentState(:), SpeciesResidual(:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Dnu_0, Chi_EmAb, Sigma_Iso, &
    !$OMP   Chi_NNS, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, &
    !$OMP   Fm, Gm, Error, ActiveMask, Fnorm0, Omega, &
    !$OMP   SchurGain, SchurDamping, &
    !$OMP   RadiationLimiterCorrection, T_old, Ye_old, Ym_old )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Dnu_0, Chi_EmAb, Sigma_Iso, &
    !$ACC   Chi_NNS, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, &
    !$ACC   Fm, Gm, Error, ActiveMask, Fnorm0, Omega, &
    !$ACC   SchurGain, SchurDamping, &
    !$ACC   RadiationLimiterCorrection, T_old, Ye_old, Ym_old )
#endif

    ALLOCATE( CurrentState(n_FP_inner), SpeciesResidual(nSpecies) )

    DO iN_X = 1, nX_G

      IF( .NOT. EvaluationMask(iN_X) ) CYCLE

      CALL PackRadiationStateForTrace_OrderV &
             ( iN_X, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, CurrentState )

      CurrentNorm = MAX( MAXVAL(ABS(CurrentState)), SqrtTiny )
      StepRelative = Zero
      Step2Relative = Zero
      IF( HasPrevious(iN_X) ) &
        StepRelative = MAXVAL(ABS(CurrentState-PreviousState(:,iN_X))) &
                         / CurrentNorm
      IF( HasPrevious2(iN_X) ) &
        Step2Relative = MAXVAL(ABS(CurrentState-Previous2State(:,iN_X))) &
                          / CurrentNorm

      ResidualAbsolute = MAXVAL(ABS(Fm(:,iN_X)))
      WorstVectorIndex = MAXLOC(ABS(Fm(:,iN_X)),DIM=1)
      WorstMoment = MOD(WorstVectorIndex-1,nCR) + 1
      iV = (WorstVectorIndex-1) / nCR
      WorstEnergy = MOD(iV,nE_G) + 1
      WorstSpecies = iV / nE_G + 1
      MapMaximum = MAXVAL(ABS(Gm(:,iN_X)))

      ResidualNorm2 = DOT_PRODUCT(Fm(:,iN_X),Fm(:,iN_X))
      ResidualPreviousNorm2 = DOT_PRODUCT &
        ( PreviousResidual(:,iN_X), PreviousResidual(:,iN_X) )
      ResidualDot = DOT_PRODUCT(Fm(:,iN_X),PreviousResidual(:,iN_X))
      ResidualChange = Zero
      ResidualCosine = Zero
      IF( HasPrevious(iN_X) )THEN
        ResidualChange = MAXVAL &
          ( ABS(Fm(:,iN_X)-PreviousResidual(:,iN_X)) ) &
            / MAX(ResidualAbsolute,SqrtTiny)
        IF( ResidualNorm2 > Zero .AND. ResidualPreviousNorm2 > Zero ) &
          ResidualCosine = ResidualDot &
            / SQRT(ResidualNorm2*ResidualPreviousNorm2)
      END IF

      ResidualRelative = Zero
      ResidualDensityRelative = Zero
      ResidualMomentumRelative = Zero
      ConvergenceMeasure = Zero
      DO iS = 1, nSpecies
        Fnorm_N  = Zero
        Fnorm_G1 = Zero
        Fnorm_G2 = Zero
        Fnorm_G3 = Zero
        DnuNormLocal = Zero
        DO iN_E = 1, nE_G
          iOS = ( (iN_E-1) + (iS-1)*nE_G ) * nCR
          Fnorm_N = Fnorm_N + W2_S(iN_E)*Fm(iOS+iCR_N,iN_X)**2
          Fnorm_G1 = Fnorm_G1 + W2_S(iN_E)*Fm(iOS+iCR_G1,iN_X)**2
          Fnorm_G2 = Fnorm_G2 + W2_S(iN_E)*Fm(iOS+iCR_G2,iN_X)**2
          Fnorm_G3 = Fnorm_G3 + W2_S(iN_E)*Fm(iOS+iCR_G3,iN_X)**2
          DnuNormLocal = DnuNormLocal &
            + W2_S(iN_E)*Dnu(iN_E,iS,iN_X)**2
        END DO
        Fnorm = SQRT(MAX(Fnorm_N,Fnorm_G1,Fnorm_G2,Fnorm_G3))
        DnuNormLocal = SQRT(DnuNormLocal)
        IF( DnuNormLocal > Zero )THEN
          SpeciesResidual(iS) = Fnorm / DnuNormLocal
        ELSE IF( Fnorm == Zero )THEN
          SpeciesResidual(iS) = Zero
        ELSE
          SpeciesResidual(iS) = HUGE(One)
        END IF
        ResidualRelative = MAX(ResidualRelative,SpeciesResidual(iS))
        IF( DnuNormLocal > Zero )THEN
          ResidualDensityRelative = MAX &
            ( ResidualDensityRelative, SQRT(Fnorm_N)/DnuNormLocal )
          ResidualMomentumRelative = MAX &
            ( ResidualMomentumRelative, &
              SQRT(MAX(Fnorm_G1,Fnorm_G2,Fnorm_G3))/DnuNormLocal )
        END IF
        InitialCriterion = Atol_inner*Fnorm0(iS,iN_X)
        Denominator = MAX(Rtol_inner*DnuNormLocal,InitialCriterion)
        IF( Denominator > Zero )THEN
          ConvergenceMeasure = MAX(ConvergenceMeasure,Fnorm/Denominator)
        ELSE IF( Fnorm > Zero )THEN
          ConvergenceMeasure = HUGE(One)
        END IF
      END DO

      ResidualRatio = Zero
      IF( HasPrevious(iN_X) ) &
        ResidualRatio = ResidualAbsolute &
          / MAX(MAXVAL(ABS(PreviousResidual(:,iN_X))),SqrtTiny)

      DnuMinimum = MINVAL(Dnu(:,:,iN_X))
      DnuMaximum = MAXVAL(Dnu(:,:,iN_X))
      MuEquilibriumDeviation = Zero
      MuEquilibriumNorm = Zero
      MuChiMaximum = Zero
      MuStiffnessMaximum = Zero
      CollisionStiffnessMaximum = Zero
      MuLeptonCurrent = Zero
      MuLeptonInvariantError = Zero
      MuLeptonInvariantRelative = Zero
      MuLeptonScale = Zero
      Mumu = Zero
      EtaMu = Zero
      MuCCActive = 0
      IF( nSpecies >= iNuM_Bar )THEN
        MuCCActive = MERGE(1,0,Opacity_Mask(iOp_EmAb_Muon,iN_X))
        CALL ComputeMuonChemicalPotential_TABLE &
               ( D(iN_X), T(iN_X), Ye(iN_X), Ym(iN_X), Mumu )
        IF( T(iN_X) > Zero ) &
          EtaMu = Mumu / ( BoltzmannConstant * T(iN_X) )
        MuChiMaximum = MAXVAL(Chi_EmAb(:,iNuM:iNuM_Bar,iN_X))
        MuStiffnessMaximum = Omega(iN_X)*dt*MuChiMaximum
        DO iS = iNuM, iNuM_Bar
        DO iN_E = 1, nE_G
          MuEquilibriumDeviation = MuEquilibriumDeviation &
            + W2_S(iN_E) &
                * (Dnu(iN_E,iS,iN_X)-Dnu_0(iN_E,iS,iN_X))**2
          MuEquilibriumNorm = MuEquilibriumNorm &
            + W2_S(iN_E)*Dnu_0(iN_E,iS,iN_X)**2
          vDotI = V_u_1(iN_X)*Gm_dd_11(iN_X)*Inu_u_1(iN_E,iS,iN_X) &
            + V_u_2(iN_X)*Gm_dd_22(iN_X)*Inu_u_2(iN_E,iS,iN_X) &
            + V_u_3(iN_X)*Gm_dd_33(iN_X)*Inu_u_3(iN_E,iS,iN_X)
          MuLeptonCurrent = MuLeptonCurrent + W2_S(iN_E)*LeptonNumber(iS) &
            * (Dnu(iN_E,iS,iN_X)+vDotI)
          MuLeptonScale = MuLeptonScale &
            + ABS( W2_S(iN_E) * LeptonNumber(iS) &
                * (Dnu(iN_E,iS,iN_X)+vDotI) )
        END DO
        END DO
        MuEquilibriumDeviation = SQRT(MuEquilibriumDeviation) &
          / MAX(SQRT(MuEquilibriumNorm),SqrtTiny)
        MuLeptonInvariantError = MuLeptonCurrent &
          - RadiationLeptonInvariant(2,iN_X)
        IF( MuCCActive == 1 ) &
          MuLeptonInvariantError = MuLeptonInvariantError &
            + D(iN_X) / AtomicMassUnit &
                * ( Ym(iN_X) - Ym_old(iN_X) ) &
            - RadiationLimiterCorrection(iYm,iN_X)
        IF( MuCCActive == 1 ) &
          MuLeptonScale = MuLeptonScale &
            + ABS( D(iN_X) / AtomicMassUnit &
                * ( Ym(iN_X) - Ym_old(iN_X) ) ) &
            + ABS( RadiationLimiterCorrection(iYm,iN_X) )
        MuLeptonScale = MAX &
          ( MuLeptonScale, ABS( RadiationLeptonInvariant(2,iN_X) ) )
        IF( MuLeptonScale > Zero ) &
          MuLeptonInvariantRelative = ABS( MuLeptonInvariantError ) &
            / MuLeptonScale
      END IF
      DO iS = 1, nSpecies
        iChiral = nChirals - ( LeptonNumber(iS) + 1 ) / nChirals
        DO iN_E = 1, nE_G
          CollisionOpacity = Chi_EmAb(iN_E,iS,iN_X) &
            + Sigma_Iso(iN_E,iChiral,iN_X) &
            + Chi_NNS(iN_E,iS,iN_X) + Chi_NES(iN_E,iS,iN_X) &
            + Chi_Pair(iN_E,iS,iN_X) + Chi_NuPair(iN_E,iS,iN_X) &
            + Chi_Brem(iN_E,iS,iN_X)
          CollisionStiffnessMaximum = MAX &
            ( CollisionStiffnessMaximum, Omega(iN_X)*dt*CollisionOpacity )
        END DO
      END DO

      WRITE(RadiationOnlyTraceUnit,*) &
        RadiationOnlySolveCallCount, k, iN_X, MERGE(1,0,ActiveMask(iN_X)), &
        MERGE(0,1,ActiveMask(iN_X)), Error(iN_X), BacktrackCount(iN_X), &
        dt, Omega(iN_X)*dt, SchurGain(iN_X), SchurDamping(iN_X), &
        D(iN_X)/Unit_D, T(iN_X)/MeV, Ye(iN_X), Ym(iN_X), &
        (T(iN_X)-T_old(iN_X))/MAX(ABS(T_old(iN_X)),SqrtTiny), &
        Ye(iN_X)-Ye_old(iN_X), Ym(iN_X)-Ym_old(iN_X), E(iN_X)/Unit_E, &
        Mumu/MeV, EtaMu, &
        MuCCActive, MuLeptonInvariantError, &
        MuLeptonInvariantRelative, &
        ResidualRelative, ResidualDensityRelative, &
        ResidualMomentumRelative, ConvergenceMeasure, ResidualAbsolute, &
        ResidualRatio, ResidualChange, ResidualCosine, &
        StepRelative, Step2Relative, CandidateStep(iN_X), &
        ProjectedStep(iN_X), ProjectionCorrection(iN_X), &
        CurrentNorm, MapMaximum, DnuMinimum, DnuMaximum, &
        MuEquilibriumDeviation, MuChiMaximum, MuStiffnessMaximum, &
        CollisionStiffnessMaximum, &
        WorstSpecies, WorstEnergy, WorstMoment, &
        SpeciesResidual, AlphaTrace(:,iN_X)

      IF( .NOT. ActiveMask(iN_X) .OR. k >= MaxIter_inner )THEN
        DO iS = 1, nSpecies
        iChiral = nChirals - ( LeptonNumber(iS) + 1 ) / nChirals
        DO iN_E = 1, nE_G
          iOS = ( (iN_E-1) + (iS-1)*nE_G ) * nCR
          WRITE(RadiationOnlyTraceSpectrumUnit,*) &
            RadiationOnlySolveCallCount, k, iN_X, &
            MERGE(0,1,ActiveMask(iN_X)), iS, iN_E, E_N(iN_E)/MeV, &
            Dnu(iN_E,iS,iN_X), Dnu_0(iN_E,iS,iN_X), &
            Inu_u_1(iN_E,iS,iN_X), Inu_u_2(iN_E,iS,iN_X), &
            Inu_u_3(iN_E,iS,iN_X), &
            Fm(iOS+iCR_N,iN_X), Fm(iOS+iCR_G1,iN_X), &
            Fm(iOS+iCR_G2,iN_X), Fm(iOS+iCR_G3,iN_X), &
            Gm(iOS+iCR_N,iN_X), Gm(iOS+iCR_G1,iN_X), &
            Gm(iOS+iCR_G2,iN_X), Gm(iOS+iCR_G3,iN_X), &
            Chi_EmAb(iN_E,iS,iN_X), &
            Sigma_Iso(iN_E,iChiral,iN_X), &
            Chi_NNS(iN_E,iS,iN_X), Chi_NES(iN_E,iS,iN_X), &
            Chi_Pair(iN_E,iS,iN_X), Chi_NuPair(iN_E,iS,iN_X), &
            Chi_Brem(iN_E,iS,iN_X)
        END DO
        END DO
      END IF

      Previous2State(:,iN_X) = PreviousState(:,iN_X)
      PreviousState(:,iN_X) = CurrentState
      PreviousResidual(:,iN_X) = Fm(:,iN_X)
      HasPrevious2(iN_X) = HasPrevious(iN_X)
      HasPrevious(iN_X) = .TRUE.
      BacktrackCount(iN_X) = 0

    END DO

    FLUSH(RadiationOnlyTraceUnit)
    FLUSH(RadiationOnlyTraceSpectrumUnit)
    DEALLOCATE( CurrentState, SpeciesResidual )

  END SUBROUTINE WriteRadiationOnlyIterationTrace_OrderV


  SUBROUTINE CaptureRadiationOnlyUpdateTrace_OrderV &
    ( MASK, Mk, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm, Alpha, PreviousState, &
      CandidateStep, ProjectedStep, ProjectionCorrection, AlphaTrace )

    LOGICAL, DIMENSION(:), INTENT(in) :: MASK
    INTEGER, INTENT(in) :: Mk
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:), INTENT(in) :: Gm, Alpha, PreviousState
    REAL(DP), DIMENSION(:), INTENT(out) :: &
      CandidateStep, ProjectedStep, ProjectionCorrection
    REAL(DP), DIMENSION(:,:), INTENT(out) :: AlphaTrace

    INTEGER :: iN_X
    REAL(DP) :: Scale
    REAL(DP), ALLOCATABLE :: ProjectedState(:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, Gm, Alpha )
#endif

    ALLOCATE( ProjectedState(n_FP_inner) )
    DO iN_X = 1, nX_G
      IF( .NOT. MASK(iN_X) ) CYCLE
      CALL PackRadiationStateForTrace_OrderV &
             ( iN_X, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, ProjectedState )
      Scale = MAX(MAXVAL(ABS(PreviousState(:,iN_X))),SqrtTiny)
      CandidateStep(iN_X) = MAXVAL(ABS(Gm(:,iN_X)-PreviousState(:,iN_X))) &
                              / Scale
      ProjectedStep(iN_X) = MAXVAL &
        (ABS(ProjectedState-PreviousState(:,iN_X))) / Scale
      ProjectionCorrection(iN_X) = MAXVAL(ABS(ProjectedState-Gm(:,iN_X))) &
                                      / Scale
      AlphaTrace(:,iN_X) = Zero
      AlphaTrace(1:Mk,iN_X) = Alpha(1:Mk,iN_X)
    END DO
    DEALLOCATE( ProjectedState )

  END SUBROUTINE CaptureRadiationOnlyUpdateTrace_OrderV




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

          IF( Fnorm /= Fnorm .OR. Fnorm > HUGE(One) )THEN
            RelativeResidual = HUGE(One)
            Converged = .FALSE.
          ELSE IF( DnuNorm(iS,iN_X) > Zero )THEN
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
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, SqrtGm, &
      ProjectionMask, ContinuationMask, &
      nX_P, PackIndex, UnpackIndex, &
      FVECm_inner, GVECm_inner, &
      FVECm_projection, GVECm_projection, &
      T_Conservative, Ye_Conservative, Ym_Conservative, E_Conservative, &
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
    REAL(DP), DIMENSION(:),     INTENT(inout) :: T, Ye, Ym, E
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
      T_Conservative, Ye_Conservative, Ym_Conservative, E_Conservative
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
      Ye_Conservative(iN_X) = Ye(iN_X)
      Ym_Conservative(iN_X) = Ym(iN_X)
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
                    ( GVECm_projection(iYe,iN_X) &
                      - Ye_Conservative(iN_X) / Ye_old(iN_X) ), &
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

          ProjectionResidual &
            = MAX( ProjectionResidual, &
                ABS( GVECm_projection(iYm,iN_X) &
                       - Ym_Conservative(iN_X) ) )

          Ye_Conservative(iN_X) &
            = GVECm_projection(iYe,iN_X) * Ye_old(iN_X)
          Ym_Conservative(iN_X) = GVECm_projection(iYm,iN_X)
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
             ( D, E_Conservative, Ye_Conservative, Ym_Conservative, &
               T_Conservative, &
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
    ! definition as MatterUpdate_Rate_Option.

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
        Ye(iN_X) = Ye_Conservative(iN_X)
        Ym(iN_X) = Ym_Conservative(iN_X)
        E(iN_X) = E_Conservative(iN_X)
        V_u_1(iN_X) = V_u_1_Conservative(iN_X)
        V_u_2(iN_X) = V_u_2_Conservative(iN_X)
        V_u_3(iN_X) = V_u_3_Conservative(iN_X)

        U_Ye(iN_X) = Ye_Conservative(iN_X) / Ye_old(iN_X)
        U_Ym(iN_X) = Ym_Conservative(iN_X)
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
               ( D, T_Conservative, Ye_Conservative, Ym_Conservative, SqrtGm, &
                 ProjectionMask, nX_P, PackIndex, UnpackIndex, &
                 T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, &
                 Opacity_Mask_Input = Opacity_Mask )

        CALL TimersStop( Timer_Collisions_ComputeOpacity )

      END IF

      CALL TimersStart( Timer_Collisions_ComputeRates )

      IF( FreezeOpacities )THEN
        CALL ComputeRates_Packed &
               ( D, T_Opacity, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 ProjectionMask, nX_P, PackIndex, UnpackIndex, &
                 nX_P, Opacity_Mask )
      ELSE
        CALL ComputeRates_Packed &
               ( D, T_Conservative, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                 ProjectionMask, nX_P, PackIndex, UnpackIndex, &
                 nX_P, Opacity_Mask )
      END IF

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
    ( D, T, Ye, Ym, SqrtGm, MASK, nX_P, PackIndex, &
      UnpackIndex, nX_P0, T0, Ye0, Ym0, Opacity_Mask_Input )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET   :: D, T, Ye, Ym, SqrtGm
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: MASK
    INTEGER,                INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                INTENT(in), OPTIONAL :: nX_P0
    REAL(DP), DIMENSION(:), INTENT(in), TARGET, OPTIONAL :: T0, Ye0, Ym0
    LOGICAL, DIMENSION(:,:), INTENT(in), TARGET :: Opacity_Mask_Input

    INTEGER                             :: nX, nX0, iX, iE, iC
    REAL(DP), DIMENSION(:)    , POINTER :: D_P, T_P, Ye_P, Ym_P, SqrtGm_P
    REAL(DP), DIMENSION(:)    , POINTER :: T0_P, Ye0_P, Ym0_P
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
    IF ( nX <= 0 ) RETURN

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
      Ye_P => Ye_T(1:nX)
      Ym_P => Ym_T(1:nX)

      IF ( PRESENT ( T0 ) ) THEN
        T0_P => T0_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, T0, T0_P ) 
      ELSE
        T0_P => T_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, T, T0_P ) 
      ENDIF

      IF ( PRESENT ( Ye0 ) ) THEN
        Ye0_P => Ye0_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, Ye0, Ye0_P )
      ELSE
        Ye0_P => Ye_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, Ye, Ye0_P )
      ENDIF

      IF ( PRESENT ( Ym0 ) ) THEN
        Ym0_P => Ym0_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, Ym0, Ym0_P )
      ELSE
        Ym0_P => Ym_T(1:nX)

        CALL ArrayPack &
               (nX, UnpackIndex, Ym, Ym0_P )
      ENDIF

      SqrtGm_P => SqrtGm_T(1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, D, T, Ye, D_P, T_P, Ye_P )
      CALL ArrayPack &
             ( nX, UnpackIndex, Ym, SqrtGm, Ym_P, SqrtGm_P )

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
      Ye_P => Ye(:)
      Ym_P => Ym(:)

      IF ( PRESENT ( T0 ) ) THEN
        T0_P => T0(:)
      ELSE
        T0_P => T(:)
      ENDIF

      IF ( PRESENT ( Ye0 ) ) THEN
        Ye0_P => Ye0(:)
      ELSE
        Ye0_P => Ye(:)
      ENDIF

      IF ( PRESENT ( Ym0 ) ) THEN
        Ym0_P => Ym0(:)
      ELSE
        Ym0_P => Ym(:)
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
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Ye_P, Ym_P, &
             Dnu_0_P, OpacityMask_Option = Opacity_Mask_P )

    CALL TimersStop( Timer_Opacity_D0 )

    ! Dnu_0 is a thermodynamic quantity entering charged-current detailed
    ! balance.  Limiting it to the bound imposed on the evolved radiation
    ! polynomial changes the microscopic inverse rate and is therefore not
    ! permissible here.  Bounds on the evolved state are enforced separately.

    ! --- EmAb ---

    CALL TimersStart( Timer_Opacity_EC )

    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Ye_P, Ym_P, &
             Dnu_0_P, Chi_EmAb_P, Opacity_Mask_P, &
             T_old = T0_P, Ye_old = Ye0_P, &
             Ym_old = Ym0_P )

    !Chi_EmAb_P = Zero

    CALL TimersStop( Timer_Opacity_EC )

    ! --- Isoenergetic scattering ---

    CALL TimersStart( Timer_Opacity_ES )

    CALL ComputeNeutrinoOpacities_ES &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Ye_P, Ym_P, &
             iM0, Phi_0_Iso_P, Opacity_Mask_P )

    CALL ComputeNeutrinoOpacities_ES &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Ye_P, Ym_P, &
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
             ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM0, &
               Phi_NNS_0_P, Phi_NbNS_0_P, &
               ApplyWeakMagnetism = NNS_ApplyWeakMagnetism, &
               ApplyManyBodyCorrection = NNS_ApplyManyBodyCorrection )

      IF ( UseFirstMomentRates ) &
        CALL ComputeNeutrinoOpacities_NNS &
               ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM1, &
                 Phi_NNS_1_P, Phi_NbNS_1_P, &
                 ApplyWeakMagnetism = NNS_ApplyWeakMagnetism, &
                 ApplyManyBodyCorrection = NNS_ApplyManyBodyCorrection )

      CALL TimersStop( Timer_Opacity_NNS )

    END IF

    IF( Include_NES )THEN

      ! --- NES Scattering Functions ---

      CALL TimersStart( Timer_Opacity_NES )

      CALL ComputeNeutrinoOpacities_NES &
             ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM0, H_I_0_P, H_II_0_P )

      IF( UseFirstMomentRates )THEN

        CALL ComputeNeutrinoOpacities_NES &
               ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM1, H_I_1_P, H_II_1_P )

      END IF

      CALL TimersStop( Timer_Opacity_NES )

    END IF

    IF( Include_Pair )THEN

      ! --- Pair Kernels ---

      CALL TimersStart( Timer_Opacity_Pair )

      CALL ComputeNeutrinoOpacities_Pair &
             ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM0, J_I_0_P, J_II_0_P )

      IF( UseFirstMomentRates )THEN

        CALL ComputeNeutrinoOpacities_Pair &
               ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM1, J_I_1_P, J_II_1_P )

      END IF

      CALL TimersStop( Timer_Opacity_Pair )

    END IF

    IF( Include_NuPair )THEN

      ! --- NuPair Kernels ---

      IF( ANY( Opacity_Mask_Input(iOp_NuPair,:) ) ) THEN

        CALL TimersStart( Timer_Opacity_NuPair )

        CALL ComputeNeutrinoOpacities_NuPair &
               ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM0, &
                 Nu_J_I_0_P, Nu_J_II_0_P )

        IF ( UseFirstMomentRates ) &
          CALL ComputeNeutrinoOpacities_NuPair &
                 ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, iM1, &
                   Nu_J_I_1_P, Nu_J_II_1_P )

        CALL TimersStop( Timer_Opacity_NuPair )

      END IF

    END IF

    IF( Include_Brem )THEN

      ! --- Brem Kernels ---

      CALL TimersStart( Timer_Opacity_Brem )

      CALL ComputeNeutrinoOpacities_Brem &
             ( 1, nE_G, 1, nX, D_P, T_P, Ye_P, Ym_P, S_Sigma_P )

      CALL TimersStop( Timer_Opacity_Brem )

    END IF

    IF( Include_NMS .OR. Include_LFE .OR. Include_LFC .OR. Include_IMD )THEN

      CALL TimersStart( Timer_Opacity_GeneralLepton )
      CALL ComputeNeutrinoLeptonKernels &
             ( D, T, Ye, Ym, &
               Opacity_Mask_Input(iOp_NMS,:), &
               Opacity_Mask_Input(iOp_LFE,:), &
               Opacity_Mask_Input(iOp_LFC,:), &
               Opacity_Mask_Input(iOp_IMD,:), &
               Include_NMS, Include_LFE, Include_LFC, Include_IMD, &
               nX, UnpackIndex )
      CALL TimersStop( Timer_Opacity_GeneralLepton )

    END IF

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, T_P, T_Opacity )
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

    ELSE

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iX = 1, nX_G
        T_Opacity(iX) = T(iX)
      END DO

    END IF

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( D, T_Rates, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      MASK, nX_P, PackIndex, UnpackIndex, &
      nX_P0, Opacity_Mask_Input )

    REAL(DP), DIMENSION(:),     INTENT(in), TARGET   :: D, T_Rates
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Dnu
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_1
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_2
    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: Inu_u_3
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: MASK
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:),     INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P0
    LOGICAL,  DIMENSION(:,:),   INTENT(in), TARGET   :: Opacity_Mask_Input

    ! T_Rates must be the temperature used to construct the opacity kernels
    ! consumed by this call.  Passing it explicitly avoids coupling detailed
    ! balance to the lifetime or packing order of hidden module state.

    REAL(DP), DIMENSION(:)    , POINTER :: D_P, T_P
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
      T_P => T_T(1:nX)

      CALL ArrayPack( nX, UnpackIndex, D, D_P )
      CALL ArrayPack( nX, UnpackIndex, T_Rates, T_P )

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
      T_P => T_Rates(:)

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

    CALL ComputeDetailedBalanceFactors &
           ( 1, nE_G, 1, nX, E_N, T_P, &
             Include_NNS .OR. Include_NES, &
             Include_Pair .OR. Include_NuPair .OR. Include_Brem, &
             DetailedBalanceScattering_T(:,:,1:nX), &
             DetailedBalancePair_T(:,:,1:nX) )

    ! --- NNS Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NNS )

    CALL ComputeNeutrinoOpacityRates_NNS &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
             Dnu_P, &
             Phi_NNS_0_P, Phi_NbNS_0_P, Eta_NNS_P, Chi_NNS_P, R_NNS_P, &
             Opacity_Mask_P, DetailedBalanceScattering_T(:,:,1:nX) )

    IF ( Include_NNS .AND. UseFirstMomentRates ) &
      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NNS &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, &
               Phi_NNS_1_P, Phi_NbNS_1_P, &
               L_NNS__In__u_1_P, L_NNS__In__u_2_P, L_NNS__In__u_3_P, &
               L_NNS__Out_u_1_P, L_NNS__Out_u_2_P, L_NNS__Out_u_3_P, &
               Opacity_Mask_P, DetailedBalanceScattering_T(:,:,1:nX) )

    CALL TimersStop( Timer_OpacityRate_NNS )

    ! --- NES Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NES )

    CALL ComputeNeutrinoOpacityRates_NES &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
             Dnu_P, &
             H_I_0_P, H_II_0_P, Eta_NES_P, Chi_NES_P, R_NES_P, &
             Opacity_Mask_P, DetailedBalanceScattering_T(:,:,1:nX) )

    IF( Include_NES .AND. UseFirstMomentRates )THEN

      ! --- Compute Linear Rate Corrections (NES) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NES &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, H_I_1_P, H_II_1_P, &
               L_NES__In__u_1_P, L_NES__In__u_2_P, L_NES__In__u_3_P, &
               L_NES__Out_u_1_P, L_NES__Out_u_2_P, L_NES__Out_u_3_P, &
               Opacity_Mask_P, DetailedBalanceScattering_T(:,:,1:nX) )

    END IF

    CALL TimersStop( Timer_OpacityRate_NES )

    ! --- Pair Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_Pair )

    CALL ComputeNeutrinoOpacityRates_Pair &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
             Dnu_P, &
             J_I_0_P, J_II_0_P, Eta_Pair_P, Chi_Pair_P, R_Pair_P, &
             Opacity_Mask_P, DetailedBalancePair_T(:,:,1:nX) )

    IF( Include_Pair .AND. UseFirstMomentRates )THEN

      ! --- Compute Linear Rate Corrections (Pair) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, J_I_1_P, J_II_1_P, &
               L_Pair_Pro_u_1_P, L_Pair_Pro_u_2_P, L_Pair_Pro_u_3_P, &
               L_Pair_Ann_u_1_P, L_Pair_Ann_u_2_P, L_Pair_Ann_u_3_P, &
               Opacity_Mask_P, DetailedBalancePair_T(:,:,1:nX) )

    END IF

    CALL TimersStop( Timer_OpacityRate_Pair )

    ! --- NuPair Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_NuPair )

    CALL ComputeNeutrinoOpacityRates_NuPair &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
             Dnu_P, &
             Nu_J_I_0_P, Nu_J_II_0_P, Eta_NuPair_P, Chi_NuPair_P, &
             R_NuPair_P, Opacity_Mask_P, DetailedBalancePair_T(:,:,1:nX) )

    IF ( Include_NuPair .AND. UseFirstMomentRates ) &
      CALL ComputeNeutrinoOpacityRates_LinearCorrections_NuPair &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, &
               Nu_J_I_1_P, Nu_J_II_1_P, &
               L_NuPair_Pro_u_1_P, L_NuPair_Pro_u_2_P, L_NuPair_Pro_u_3_P, &
               L_NuPair_Ann_u_1_P, L_NuPair_Ann_u_2_P, L_NuPair_Ann_u_3_P, &
               Opacity_Mask_P, DetailedBalancePair_T(:,:,1:nX) )

    CALL TimersStop( Timer_OpacityRate_NuPair )

    ! --- Brem Emissivities and Opacities ---

    CALL TimersStart( Timer_OpacityRate_Brem )

    CALL ComputeNeutrinoOpacityRates_Brem &
           ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
             Dnu_P, &
             S_Sigma_P, Eta_Brem_P, Chi_Brem_P, R_Brem_P, &
             Opacity_Mask_P, DetailedBalancePair_T(:,:,1:nX) )

    IF( Include_Brem .AND. UseFirstMomentRates )THEN

      ! --- Compute Linear Rate Corrections (Brem) ---

      CALL ComputeNeutrinoOpacityRates_LinearCorrections_Brem &
             ( 1, nE_G, 1, nSpecies, 1, nX, D_P, T_P, E_N, W2_N, &
               Inu_u_1_P, Inu_u_2_P, Inu_u_3_P, S_Sigma_P, &
               L_Brem_Pro_u_1_P, L_Brem_Pro_u_2_P, L_Brem_Pro_u_3_P, &
               L_Brem_Ann_u_1_P, L_Brem_Ann_u_2_P, L_Brem_Ann_u_3_P, &
               Opacity_Mask_P, DetailedBalancePair_T(:,:,1:nX) )

    END IF

    CALL TimersStop( Timer_OpacityRate_Brem )

    IF( Include_NMS .OR. Include_LFE .OR. Include_LFC .OR. Include_IMD )THEN

      CALL TimersStart( Timer_OpacityRate_GeneralLepton )
      CALL ComputeNeutrinoLeptonRates &
             ( T_Rates, W2_N, Dnu, &
               Opacity_Mask_Input(iOp_NMS,:), &
               Opacity_Mask_Input(iOp_LFE,:), &
               Opacity_Mask_Input(iOp_LFC,:), &
               Opacity_Mask_Input(iOp_IMD,:), &
               Include_NMS, Include_LFE, Include_LFC, Include_IMD, &
               nX, UnpackIndex )
      CALL TimersStop( Timer_OpacityRate_GeneralLepton )

    END IF

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
    ( D, E, Ye, Ym, T, MASK, nX_P, PackIndex, UnpackIndex, Error, &
      DisplayErrors_Option )

    REAL(DP), DIMENSION(:), INTENT(in)   , TARGET   :: D, E, Ye, Ym
    REAL(DP), DIMENSION(:), INTENT(inout), TARGET   :: T
    LOGICAL,  DIMENSION(:), INTENT(in)   , OPTIONAL :: MASK
    INTEGER,                INTENT(in)   , OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in)   , OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,  DIMENSION(:), INTENT(inout), TARGET   :: Error
    LOGICAL,                INTENT(in)   , OPTIONAL :: DisplayErrors_Option

    REAL(DP), DIMENSION(:), POINTER :: D_P, E_P, Ye_P, Ym_P, T_P
    INTEGER,  DIMENSION(:), POINTER :: Error_P

    INTEGER :: nX

    IF( PRESENT( nX_P ) )THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P  => D_T (1:nX)
      T_P  => T_T (1:nX)
      Ye_P => Ye_T(1:nX)
      Ym_P => Ym_T(1:nX)
      E_P  => E_T (1:nX)
      Error_P => Error_T(1:nX)

        CALL ArrayPack &
             ( nX, UnpackIndex, D, Ye, Ym, D_P, Ye_P, Ym_P )
        CALL ArrayPack &
             ( nX, UnpackIndex, T, E, T_P, E_P )
      CALL ArrayPack &
             ( nX, UnpackIndex, Error, Error_P )

    ELSE

      D_P  => D (:)
      Ye_P => Ye(:)
      Ym_P => Ym(:)
      T_P  => T (:)
      E_P  => E (:)
      Error_P => Error(:)

    END IF

    IF( PRESENT(DisplayErrors_Option) )THEN
      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D_P, E_P, Ye_P, Ym_P, T_P, Error_Option = Error_P, &
               DisplayErrors_Option = DisplayErrors_Option )
    ELSE
      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D_P, E_P, Ye_P, Ym_P, T_P, Error_Option = Error_P )
    END IF

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
      LimiterRate(iQ_Le,iN_X) = DeltaLepton / dt
      LimiterRate(iQ_Lm,iN_X) = DeltaLeptonMuon / dt
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
      ConservationCorrection(iYe ,iN_X) = -DeltaLepton
      ConservationCorrection(iYm ,iN_X) = -DeltaLeptonMuon
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
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Ye, Ym, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Ye, Ym, E
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: SUM_Ye, SUM_Ym, SUM_Ef
    REAL(DP) :: wYe, wYm, wEf

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, wYe, wYm, wEf )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, wYe, wYm, wEf )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, wYe, wYm, wEf )
#endif
    DO iN_X = 1, nX_G

      wYe = MatterSourceWeight( iYe, iN_X )
      wYm = MatterSourceWeight( iYm, iN_X )
      wEf = MatterSourceWeight( iEf, iN_X )

      ! --- Store Initial Matter State ---

      Ye_old (iN_X) = Ye(iN_X)
      Ym_old (iN_X) = Ym(iN_X)
      Ef_old (iN_X) = E(iN_X)

      ! --- Scaling Factors ---

      S_Ye(iN_X) = One / ( D(iN_X) * Ye(iN_X) / AtomicMassUnit )
      IF (Ym(iN_X) <= 1d-100) THEN
        S_Ym (iN_X) = Zero
      ELSE
        S_Ym(iN_X)  = One / ( D(iN_X) * Ym(iN_X) / AtomicMassUnit )
      END IF
      S_Ef(iN_X) = One / ( D(iN_X) * E(iN_X) )

      ! --- Initial Guess for Matter State ---

      U_Ye(iN_X) = One
      U_Ym(iN_X) = One
      U_Ef(iN_X) = One

      SUM_Ye = Zero
      SUM_Ym = Zero
      SUM_Ef = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef )
#endif
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        ! --- Old States for Neutrino Number Density and Flux ---

        C_Dnu    (iN_E,iS,iN_X) = Dnu    (iN_E,iS,iN_X)
        C_Inu_d_1(iN_E,iS,iN_X) = Inu_u_1(iN_E,iS,iN_X) * Gm_dd_11(iN_X)
        C_Inu_d_2(iN_E,iS,iN_X) = Inu_u_2(iN_E,iS,iN_X) * Gm_dd_22(iN_X)
        C_Inu_d_3(iN_E,iS,iN_X) = Inu_u_3(iN_E,iS,iN_X) * Gm_dd_33(iN_X)

        IF ( iS <= iNuE_Bar ) THEN
          SUM_Ye  = SUM_Ye  + Dnu(iN_E,iS,iN_X) * W2_S(iN_E) * LeptonNumber(iS)
        END IF
        IF ( iS == iNuM .OR. iS == iNuM_Bar) THEN
          SUM_Ym  = SUM_Ym  + Dnu(iN_E,iS,iN_X) * W2_S(iN_E) * LeptonNumber(iS)
        END IF
        SUM_Ef = SUM_Ef + Dnu(iN_E,iS,iN_X) * W3_S(iN_E)

      END DO
      END DO

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_Ye(iN_X) = U_Ye(iN_X) + wYe * SUM_Ye * S_Ye(iN_X)
      C_Ym(iN_X) = U_Ym(iN_X) + wYm * SUM_Ym * S_Ym(iN_X)
      C_Ef(iN_X) = U_Ef(iN_X) + wEf * SUM_Ef * S_Ef(iN_X)

    END DO

  END SUBROUTINE InitializeRHS_OrderOne


  SUBROUTINE InitializeRHS_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Nnu_Clipped, &
      Gnu_d_1_Clipped, Gnu_d_2_Clipped, Gnu_d_3_Clipped, &
      RadiationLimiterCorrection )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: &
      Nnu_Clipped, Gnu_d_1_Clipped, Gnu_d_2_Clipped, Gnu_d_3_Clipped
    REAL(DP), DIMENSION(:,:)  , INTENT(in)  :: RadiationLimiterCorrection

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: vDotV, vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3, Ef
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3
    REAL(DP) :: SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3
    REAL(DP) :: wYe, wYm, wEf, wV1, wV2, wV3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( vDotV, SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, Ef, &
    !$OMP          wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( vDotV, SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          V_d_1, V_d_2, V_d_3, Ef, &
    !$ACC          wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( vDotV, SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, Ef, &
    !$OMP          wYe, wYm, wEf, wV1, wV2, wV3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iN_X = 1, nX_G

      wYe = MatterSourceWeight( iYe, iN_X )
      wYm = MatterSourceWeight( iYm, iN_X )
      wEf = MatterSourceWeight( iEf, iN_X )
      wV1 = MatterSourceWeight( iV1, iN_X )
      wV2 = MatterSourceWeight( iV2, iN_X )
      wV3 = MatterSourceWeight( iV3, iN_X )

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
      Ye_old(iN_X) = Ye(iN_X)
      Ym_old(iN_X) = Ym(iN_X)
      E_old (iN_X) = E(iN_X)
      Ef_old(iN_X) = Ef
      T_old (iN_X) = T(iN_X)

      V_u_1_old(iN_X) = V_u_1(iN_X)
      V_u_2_old(iN_X) = V_u_2(iN_X)
      V_u_3_old(iN_X) = V_u_3(iN_X)

      ! --- Scaling Factors ---

      S_Ye   (iN_X) = One / ( D(iN_X) * Ye(iN_X) / AtomicMassUnit )
      ! Use an additive muon fraction as the nonlinear matter variable.
      ! Unlike a ratio to Ym_old, this remains regular when Ym_old is zero.
      S_Ym   (iN_X) = AtomicMassUnit / D(iN_X)
      S_Ef   (iN_X) = One / ( D(iN_X) * Ef )
      S_V_d_1(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_2(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_3(iN_X) = One / ( D(iN_X) * SpeedOfLight )

      ! --- Initial Guess for Matter State ---

      U_Ye   (iN_X) = One
      U_Ym   (iN_X) = Ym(iN_X)
      U_Ef   (iN_X) = One
      U_V_d_1(iN_X) = V_d_1 / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2 / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3 / SpeedOfLight

      SUM_Ye = Zero
      SUM_Ym = Zero
      SUM_Ef = Zero
      SUM_V1 = Zero
      SUM_V2 = Zero
      SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$OMP REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$ACC REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
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
          SUM_Ye = SUM_Ye + Nnu * W2_S(iN_E) * LeptonNumber(iS)
        END IF
        IF ( iS == iNuM .OR. iS == iNuM_Bar ) THEN
          SUM_Ym = SUM_Ym + Nnu * W2_S(iN_E) * LeptonNumber(iS)
        END IF

        SUM_Ef = SUM_Ef + Enu     * W3_S(iN_E)
        SUM_V1 = SUM_V1 + Fnu_d_1 * W3_S(iN_E)
        SUM_V2 = SUM_V2 + Fnu_d_2 * W3_S(iN_E)
        SUM_V3 = SUM_V3 + Fnu_d_3 * W3_S(iN_E)

      END DO
      END DO

      ! --- Store the Old Radiation Contribution to Each Conserved Total ---
      ! --- Include the radiation removed by the occupation-number bound. ---

      C_Ye(iN_X) = wYe &
        * ( SUM_Ye + RadiationLimiterCorrection(iYe,iN_X) ) * S_Ye(iN_X)
      C_Ym(iN_X) = wYm &
        * ( SUM_Ym + RadiationLimiterCorrection(iYm,iN_X) ) * S_Ym(iN_X)
      C_Ef(iN_X) = wEf &
        * ( SUM_Ef + RadiationLimiterCorrection(iEf,iN_X) ) * S_Ef(iN_X)
      C_V_d_1(iN_X) = wV1 &
        * ( SUM_V1 + RadiationLimiterCorrection(iV1,iN_X) ) &
        * S_V_d_1(iN_X)
      C_V_d_2(iN_X) = wV2 &
        * ( SUM_V2 + RadiationLimiterCorrection(iV2,iN_X) ) &
        * S_V_d_2(iN_X)
      C_V_d_3(iN_X) = wV3 &
        * ( SUM_V3 + RadiationLimiterCorrection(iV3,iN_X) ) &
        * S_V_d_3(iN_X)

    END DO

  END SUBROUTINE InitializeRHS_OrderV


  SUBROUTINE InitializeRHS_Relativistic &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, D, Ye, Ym, E, P, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Alpha, Beta_u_1, Beta_u_2, Beta_u_3 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Ye, Ym, E, P, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Alpha
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Beta_u_1, Beta_u_2, Beta_u_3

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotV, vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, W, vvDotK
    REAL(DP) :: DT, B_d_1, B_d_2, B_d_3, cD, SJ(3), Tau_f
    REAL(DP) :: Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3
    REAL(DP) :: Inu_d_1, Inu_d_2, Inu_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: SUM_Ye, SUM_Ym, SUM_V1, SUM_V2, SUM_V3, SUM_Ef
    REAL(DP) :: wD, wYe, wYm, wEf, wV1, wV2, wV3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, vDotV, W, SJ, Tau_f, &
    !$OMP          SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          wD, wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, vDotV, W, SJ, Tau_f, &
    !$ACC          SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          wD, wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( B_d_1, B_d_2, B_d_3, V_d_1, V_d_2, V_d_3, vDotV, W, SJ, Tau_f, &
    !$OMP          SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          wD, wYe, wYm, wEf, wV1, wV2, wV3, &
    !$OMP          DT, Inu_d_1, Inu_d_2, Inu_d_3, VDotInu, &
    !$OMP          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 )
#endif
    DO iN_X = 1, nX_G

      wD  = MatterSourceWeight( iD,  iN_X )
      wYe = MatterSourceWeight( iYe, iN_X )
      wYm = MatterSourceWeight( iYm, iN_X )
      wEf = MatterSourceWeight( iEf, iN_X )
      wV1 = MatterSourceWeight( iV1, iN_X )
      wV2 = MatterSourceWeight( iV2, iN_X )
      wV3 = MatterSourceWeight( iV3, iN_X )

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

      D_old (iN_X) = D (iN_X)
      Ye_old(iN_X) = Ye(iN_X)
      Ym_old(iN_X) = Ym(iN_X)
      cD_old(iN_X) = D (iN_X) * W
      Ef_old(iN_X) = E (iN_X)

      ! --- Scaling Factors ---

      S_D    (iN_X) = One / ( D (iN_X) )
      S_Ye   (iN_X) = One / ( Ye(iN_X) )
      IF (Ym(iN_X) <= 1d-100) THEN
        S_Ym (iN_X) = Zero
      ELSE
        S_Ym (iN_X) = One / ( Ym(iN_X) )
      END IF
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
      U_Ye   (iN_X) = One
      U_Ym   (iN_X) = One
      U_Ef   (iN_X) = One
      U_V_d_1(iN_X) = V_d_1 / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2 / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3 / SpeedOfLight

      SUM_Ye = Zero
      SUM_Ym = Zero
      SUM_Ef = Zero
      SUM_V1 = Zero
      SUM_V2 = Zero
      SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( DT, Inu_d_1, Inu_d_2, Inu_d_3, VDotInu, &
    !$OMP          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$OMP REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( DT, Inu_d_1, Inu_d_2, Inu_d_3, VDotInu, &
    !$ACC          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$ACC REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
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
        SUM_Ye = SUM_Ye + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
      END IF
      IF ( iS == iNuM .OR. iS == iNuM_Bar ) THEN
        SUM_Ym = SUM_Ym + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
      END IF

      SUM_Ef  = SUM_Ef  + Enu     * W3_S(iN_E)
      SUM_V1  = SUM_V1  + Fnu_d_1 * W3_S(iN_E)
      SUM_V2  = SUM_V2  + Fnu_d_2 * W3_S(iN_E)
      SUM_V3  = SUM_V3  + Fnu_d_3 * W3_S(iN_E)

    END DO
    END DO

      SUM_Ye = ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Ye
      SUM_Ym = ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Ym

      SJ(1) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_1
      SJ(2) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_2
      SJ(3) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_3

      Tau_f = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 - P(iN_X) - cD_old(iN_X)

      ! --- Include Old Matter State in Constant (C) Terms ---


      IF ( wD .EQ. One ) THEN
        C_D    (iN_X) &
          = W * U_D(iN_X)
      ELSE  
        C_D    (iN_X) &
          = U_D(iN_X)
      END IF

      IF ( wYe .EQ. One ) THEN
        C_Ye    (iN_X) &
          = U_Ye(iN_X) + SUM_Ye * S_Ye(iN_X)
      ELSE
        C_Ye    (iN_X) &
          = U_Ye(iN_X)
      END IF
      IF ( wYm .EQ. One ) THEN
        C_Ym    (iN_X) &
          = U_Ym(iN_X) + SUM_Ym * S_Ym(iN_X)
      ELSE
        C_Ym    (iN_X) &
          = U_Ym(iN_X)
      END IF
       
      IF ( wEf .EQ. One ) THEN
        C_Ef   (iN_X) &
          = ( Tau_f + SUM_Ef ) * S_Ef  (iN_X) / cD_old(iN_X)
      ELSE
        C_Ef   (iN_X) &
          = U_Ef(iN_X)
      END IF

      IF ( wV1 .EQ. One ) THEN
        C_V_d_1(iN_X) &
          = ( SJ(1) + SUM_V1 ) * S_V_d_1(iN_X) / cD_old(iN_X)
      ELSE
        C_V_d_1(iN_X) &
          = U_V_d_1(iN_X)
      END IF


      IF ( wV2 .EQ. One ) THEN
        C_V_d_2(iN_X) &
          = ( SJ(2) + SUM_V2 ) * S_V_d_2(iN_X) / cD_old(iN_X)
      ELSE
        C_V_d_2(iN_X) &
          = U_V_d_2(iN_X)
      END IF


      IF ( wV3 .EQ. One ) THEN
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
      NNS_RateDecomposition_Option, QH_ByProcessSpecies_Option, EmAb_RateDecomposition_Option, &
      EmAb_Spectrum_Option, RadiationCollisionRHS_Option )

    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)  :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    LOGICAL,  DIMENSION(:),     INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: Q_Matter
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: NNS_RateDecomposition_Option
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: QH_ByProcessSpecies_Option
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: EmAb_RateDecomposition_Option
    REAL(DP), DIMENSION(:,:,:,:), INTENT(out), OPTIONAL :: EmAb_Spectrum_Option
    REAL(DP), DIMENSION(:,:,:,:), INTENT(out), OPTIONAL :: &
      RadiationCollisionRHS_Option

    REAL(DP), ALLOCATABLE :: Q_SP(:,:,:,:), QLm_SP(:,:,:)
    REAL(DP), ALLOCATABLE :: EmAb_Spectrum(:,:,:,:)
    REAL(DP), ALLOCATABLE :: RadiationCollisionRHS(:,:,:,:)
    REAL(DP), ALLOCATABLE :: QH_ByProcessSpecies(:,:,:)

    INTEGER  :: ierr, myid
    LOGICAL  :: ComputeRadiationCollisionRHS, InvalidShape
    INTEGER  :: iG, iN_E, iN_X, iP, iQ, iS
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: RN(nC)
    REAL(DP) :: RG1(nC)
    REAL(DP) :: RG2(nC)
    REAL(DP) :: RG3(nC)
    REAL(DP) :: QL, QLm, QE, QM1, QM2, QM3
    REAL(DP) :: RN_Total, RG1_Total, RG2_Total, RG3_Total
    REAL(DP) :: V2

    ComputeRadiationCollisionRHS = PRESENT( RadiationCollisionRHS_Option )

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

    IF( PRESENT( NNS_RateDecomposition_Option ) )THEN
      InvalidShape = InvalidShape .OR. &
                     SIZE(NNS_RateDecomposition_Option,1) &
                       /= nNNSRateDecomposition .OR. &
                     SIZE(NNS_RateDecomposition_Option,2) /= nSpecies .OR. &
                     SIZE(NNS_RateDecomposition_Option,3) /= nX_G
    END IF

    IF( PRESENT( EmAb_RateDecomposition_Option ) )THEN
      InvalidShape = InvalidShape .OR. &
                     SIZE(EmAb_RateDecomposition_Option,1) /= nEmAbRateDecomposition .OR. &
                     SIZE(EmAb_RateDecomposition_Option,2) /= nSpecies .OR. &
                     SIZE(EmAb_RateDecomposition_Option,3) /= nX_G
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
    ALLOCATE( QLm_SP(nC,nSpecies,nX_G) )
    IF( ComputeRadiationCollisionRHS )THEN
      ALLOCATE( RadiationCollisionRHS(nCR,nE_G,nSpecies,nX_G) )
    ELSE
      ALLOCATE( RadiationCollisionRHS(1,1,1,1) )
    END IF

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
    !$OMP TARGET ENTER DATA MAP( alloc: RadiationCollisionRHS )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( RadiationCollisionRHS )
#endif

    IF( ComputeRadiationCollisionRHS )THEN
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( RadiationCollisionRHS )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iN_X = 1, nX_G
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G
      DO iQ   = 1, nCR
        RadiationCollisionRHS(iQ,iN_E,iS,iN_X) = Zero
      END DO
      END DO
      END DO
      END DO
    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: Q_SP, QLm_SP )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( Q_SP, QLm_SP )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) PRESENT( QLm_SP )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iP   = 1, nC
      QLm_SP(iP,iS,iN_X) = Zero
    END DO
    END DO
    END DO

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
    !$OMP FIRSTPRIVATE( ComputeRadiationCollisionRHS ) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$OMP          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$OMP          RG3_Total, QL, QLm, QE, QM1, QM2, QM3, iN_E, iP, iG )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC FIRSTPRIVATE( ComputeRadiationCollisionRHS ) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$ACC          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$ACC          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$ACC          RG3_Total, QL, QLm, QE, QM1, QM2, QM3, iN_E, iP, iG ) &
    !$ACC PRESENT( Q_SP, QLm_SP, RadiationCollisionRHS )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP FIRSTPRIVATE( ComputeRadiationCollisionRHS ) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$OMP          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$OMP          RG3_Total, QL, QLm, QE, QM1, QM2, QM3, iN_E, iP, iG )
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

          IF( Include_NMS )THEN
            RN(iC_NMS) = R_GL(iN_E,iS,iN_X,iGL_NMS)
            RG1(iC_NMS) = -Chi_GL(iN_E,iS,iN_X,iGL_NMS) * I_d_1
            RG2(iC_NMS) = -Chi_GL(iN_E,iS,iN_X,iGL_NMS) * I_d_2
            RG3(iC_NMS) = -Chi_GL(iN_E,iS,iN_X,iGL_NMS) * I_d_3
          END IF

          IF( Include_LFE )THEN
            RN(iC_LFE) = R_GL(iN_E,iS,iN_X,iGL_LFE)
            RG1(iC_LFE) = -Chi_GL(iN_E,iS,iN_X,iGL_LFE) * I_d_1
            RG2(iC_LFE) = -Chi_GL(iN_E,iS,iN_X,iGL_LFE) * I_d_2
            RG3(iC_LFE) = -Chi_GL(iN_E,iS,iN_X,iGL_LFE) * I_d_3
          END IF

          IF( Include_LFC )THEN
            RN(iC_LFC) = R_GL(iN_E,iS,iN_X,iGL_LFC)
            RG1(iC_LFC) = -Chi_GL(iN_E,iS,iN_X,iGL_LFC) * I_d_1
            RG2(iC_LFC) = -Chi_GL(iN_E,iS,iN_X,iGL_LFC) * I_d_2
            RG3(iC_LFC) = -Chi_GL(iN_E,iS,iN_X,iGL_LFC) * I_d_3
          END IF

          IF( Include_IMD )THEN
            RN(iC_IMD) = R_GL(iN_E,iS,iN_X,iGL_IMD)
            RG1(iC_IMD) = -Chi_GL(iN_E,iS,iN_X,iGL_IMD) * I_d_1
            RG2(iC_IMD) = -Chi_GL(iN_E,iS,iN_X,iGL_IMD) * I_d_2
            RG3(iC_IMD) = -Chi_GL(iN_E,iS,iN_X,iGL_IMD) * I_d_3
          END IF

          IF( ComputeRadiationCollisionRHS )THEN
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
          END IF

          DO iP = 1, nC

            QL = Zero
            IF ( iS <= iNuE_Bar ) &
              QL = -LeptonNumber(iS) * RN(iP) * W2_S(iN_E)

            QLm = Zero
            IF ( iS == iNuM .OR. iS == iNuM_Bar ) &
              QLm = -LeptonNumber(iS) * RN(iP) * W2_S(iN_E)

            QE  = -( RN(iP) + V_u_1(iN_X) * RG1(iP) &
                              + V_u_2(iN_X) * RG2(iP) &
                              + V_u_3(iN_X) * RG3(iP) ) * W3_S(iN_E)
            QM1 = -( RG1(iP) + V_d_1 * RN(iP) ) * W3_S(iN_E)
            QM2 = -( RG2(iP) + V_d_2 * RN(iP) ) * W3_S(iN_E)
            QM3 = -( RG3(iP) + V_d_3 * RN(iP) ) * W3_S(iN_E)

            Q_SP(1,iP,iS,iN_X) = Q_SP(1,iP,iS,iN_X) + QL
            QLm_SP(iP,iS,iN_X) = QLm_SP(iP,iS,iN_X) + QLm
            Q_SP(2,iP,iS,iN_X) = Q_SP(2,iP,iS,iN_X) + QE
            Q_SP(3,iP,iS,iN_X) = Q_SP(3,iP,iS,iN_X) + QM1
            Q_SP(4,iP,iS,iN_X) = Q_SP(4,iP,iS,iN_X) + QM2
            Q_SP(5,iP,iS,iN_X) = Q_SP(5,iP,iS,iN_X) + QM3

          END DO

        END DO

      END IF

    END DO
    END DO

    IF( PRESENT( NNS_RateDecomposition_Option ) ) &
      CALL ComputeNNSRateDecomposition_OrderV &
             ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, MASK, Q_SP, &
               V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               NNS_RateDecomposition_Option )

    IF( PRESENT( EmAb_RateDecomposition_Option ) ) &
      CALL ComputeEmAbRateDecomposition_OrderV &
             ( Dnu, MASK, Q_SP, V_u_1, V_u_2, V_u_3, &
               EmAb_RateDecomposition_Option )

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
    !$ACC PRESENT( Q_SP, QLm_SP, Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iP )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
      DO iP = 1, nC
        Q_Matter(iQ_Le,1+iS,iN_X) = Q_Matter(iQ_Le,1+iS,iN_X) &
                                             + Q_SP(1,iP,iS,iN_X)
        Q_Matter(iQ_Lm,1+iS,iN_X) = Q_Matter(iQ_Lm,1+iS,iN_X) &
                                             + QLm_SP(iP,iS,iN_X)
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
    !$ACC PRESENT( Q_SP, QLm_SP, Q_Matter )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iS )
#endif
    DO iN_X = 1, nX_G
    DO iP   = 1, nC
      DO iS = 1, nSpecies
        Q_Matter(iQ_Le,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_Le,1+nSpecies+iP,iN_X) + Q_SP(1,iP,iS,iN_X)
        Q_Matter(iQ_Lm,1+nSpecies+iP,iN_X) &
          = Q_Matter(iQ_Lm,1+nSpecies+iP,iN_X) + QLm_SP(iP,iS,iN_X)
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
          Q_Matter(iQ_Le,iG,iN_X) = Q_Matter(iQ_Le,iG,iN_X) &
                                            + Q_Matter(iQ_Le,1+iS,iN_X)
          Q_Matter(iQ_Lm,iG,iN_X) = Q_Matter(iQ_Lm,iG,iN_X) &
                                            + Q_Matter(iQ_Lm,1+iS,iN_X)
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

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( release: Q_SP, QLm_SP )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE( Q_SP, QLm_SP )
#endif

    DEALLOCATE( Q_SP, QLm_SP )

  END SUBROUTINE ComputeMatterRHS_Collisions_OrderV


  SUBROUTINE ComputeNNSRateDecomposition_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, MASK, Q_SP, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      NNS_RateDecomposition )

    REAL(DP), DIMENSION(:,:,:),   INTENT(in)  :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    LOGICAL,  DIMENSION(:),       INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: Q_SP
    REAL(DP), DIMENSION(:),       INTENT(in)  :: &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:),   INTENT(out) :: NNS_RateDecomposition

    INTEGER :: iL, iN_E, iN_X, iS
    REAL(DP) :: I_d_1, I_d_2, I_d_3, L_u_1, L_u_2, L_u_3, RN_NNS

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: NNS_RateDecomposition )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( NNS_RateDecomposition )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( NNS_RateDecomposition )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iL   = 1, nNNSRateDecomposition
      NNS_RateDecomposition(iL,iS,iN_X) = Zero
    END DO
    END DO
    END DO

    IF( Include_NNS )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_E, I_d_1, I_d_2, I_d_3, L_u_1, L_u_2, L_u_3, RN_NNS )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iN_E, I_d_1, I_d_2, I_d_3, &
    !$ACC          L_u_1, L_u_2, L_u_3, RN_NNS ) &
    !$ACC PRESENT( NNS_RateDecomposition, Q_SP )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_E, I_d_1, I_d_2, I_d_3, L_u_1, L_u_2, L_u_3, RN_NNS )
#endif
      DO iN_X = 1, nX_G
      DO iS   = 1, nSpecies

        IF( MASK(iN_X) )THEN
          DO iN_E = 1, nE_G
            RN_NNS = R_NNS(iN_E,iS,iN_X)
            IF( Include_LinCorr )THEN
              I_d_1 = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)
              I_d_2 = Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)
              I_d_3 = Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)
              L_u_1 = L_NNS__In__u_1(iN_E,iS,iN_X) &
                      - L_NNS__Out_u_1(iN_E,iS,iN_X)
              L_u_2 = L_NNS__In__u_2(iN_E,iS,iN_X) &
                      - L_NNS__Out_u_2(iN_E,iS,iN_X)
              L_u_3 = L_NNS__In__u_3(iN_E,iS,iN_X) &
                      - L_NNS__Out_u_3(iN_E,iS,iN_X)
              RN_NNS = RN_NNS &
                         - ( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
            END IF
            NNS_RateDecomposition(iNNSRateDecomposition_NIn,iS,iN_X) = &
              NNS_RateDecomposition(iNNSRateDecomposition_NIn,iS,iN_X) &
                + Eta_NNS(iN_E,iS,iN_X) * W2_S(iN_E)
            NNS_RateDecomposition(iNNSRateDecomposition_NOut,iS,iN_X) = &
              NNS_RateDecomposition(iNNSRateDecomposition_NOut,iS,iN_X) &
                + Chi_NNS(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W2_S(iN_E)
            NNS_RateDecomposition(iNNSRateDecomposition_NNet,iS,iN_X) = &
              NNS_RateDecomposition(iNNSRateDecomposition_NNet,iS,iN_X) &
                + RN_NNS * W2_S(iN_E)
            NNS_RateDecomposition(iNNSRateDecomposition_EIn,iS,iN_X) = &
              NNS_RateDecomposition(iNNSRateDecomposition_EIn,iS,iN_X) &
                + Eta_NNS(iN_E,iS,iN_X) * W3_S(iN_E)
            NNS_RateDecomposition(iNNSRateDecomposition_EOut,iS,iN_X) = &
              NNS_RateDecomposition(iNNSRateDecomposition_EOut,iS,iN_X) &
                + Chi_NNS(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W3_S(iN_E)
          END DO
          NNS_RateDecomposition(iNNSRateDecomposition_QH,iS,iN_X) &
            = Q_SP(2,iC_NNS,iS,iN_X) &
                - V_u_1(iN_X) * Q_SP(3,iC_NNS,iS,iN_X) &
                - V_u_2(iN_X) * Q_SP(4,iC_NNS,iS,iN_X) &
                - V_u_3(iN_X) * Q_SP(5,iC_NNS,iS,iN_X)
        END IF

      END DO
      END DO

    END IF

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( from: NNS_RateDecomposition )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA COPYOUT( NNS_RateDecomposition )
#endif

  END SUBROUTINE ComputeNNSRateDecomposition_OrderV


  SUBROUTINE ComputeEmAbRateDecomposition_OrderV &
    ( Dnu, MASK, Q_SP, V_u_1, V_u_2, V_u_3, EmAb_RateDecomposition )

    REAL(DP), DIMENSION(:,:,:),   INTENT(in)  :: Dnu
    LOGICAL,  DIMENSION(:),       INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: Q_SP
    REAL(DP), DIMENSION(:),       INTENT(in)  :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:),   INTENT(out) :: EmAb_RateDecomposition

    INTEGER :: iL, iN_E, iN_X, iS

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( alloc: EmAb_RateDecomposition )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE( EmAb_RateDecomposition )
#endif

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( EmAb_RateDecomposition )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iL   = 1, nEmAbRateDecomposition
      EmAb_RateDecomposition(iL,iS,iN_X) = Zero
    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE( iN_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE( iN_E ) &
    !$ACC PRESENT( EmAb_RateDecomposition, Q_SP )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iN_E )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies

      IF( MASK(iN_X) )THEN
        DO iN_E = 1, nE_G
          EmAb_RateDecomposition(iEmAbRateDecomposition_NIn,iS,iN_X) = &
            EmAb_RateDecomposition(iEmAbRateDecomposition_NIn,iS,iN_X) &
              + Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) * W2_S(iN_E)
          EmAb_RateDecomposition(iEmAbRateDecomposition_NOut,iS,iN_X) = &
            EmAb_RateDecomposition(iEmAbRateDecomposition_NOut,iS,iN_X) &
              + Chi_EmAb(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W2_S(iN_E)
          EmAb_RateDecomposition(iEmAbRateDecomposition_NNet,iS,iN_X) = &
            EmAb_RateDecomposition(iEmAbRateDecomposition_NNet,iS,iN_X) &
              + Chi_EmAb(iN_E,iS,iN_X) &
                  * ( Dnu_0(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) ) &
                  * W2_S(iN_E)
          EmAb_RateDecomposition(iEmAbRateDecomposition_EIn,iS,iN_X) = &
            EmAb_RateDecomposition(iEmAbRateDecomposition_EIn,iS,iN_X) &
              + Chi_EmAb(iN_E,iS,iN_X) * Dnu_0(iN_E,iS,iN_X) * W3_S(iN_E)
          EmAb_RateDecomposition(iEmAbRateDecomposition_EOut,iS,iN_X) = &
            EmAb_RateDecomposition(iEmAbRateDecomposition_EOut,iS,iN_X) &
              + Chi_EmAb(iN_E,iS,iN_X) * Dnu(iN_E,iS,iN_X) * W3_S(iN_E)
        END DO
        EmAb_RateDecomposition(iEmAbRateDecomposition_QH,iS,iN_X) &
          = Q_SP(2,iC_EmAb,iS,iN_X) &
              - V_u_1(iN_X) * Q_SP(3,iC_EmAb,iS,iN_X) &
              - V_u_2(iN_X) * Q_SP(4,iC_EmAb,iS,iN_X) &
              - V_u_3(iN_X) * Q_SP(5,iC_EmAb,iS,iN_X)
      END IF

    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( from: EmAb_RateDecomposition )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA COPYOUT( EmAb_RateDecomposition )
#endif

  END SUBROUTINE ComputeEmAbRateDecomposition_OrderV
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
    REAL(DP) :: SUM_Ye, SUM_Ym, SUM_Ef
    REAL(DP) :: wYe, wYm, wEf

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, wYe, wYm, wEf )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, wYe, wYm, wEf )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, wYe, wYm, wEf )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        wYe = MatterSourceWeight( iYe, iN_X )
        wYm = MatterSourceWeight( iYm, iN_X )
        wEf = MatterSourceWeight( iEf, iN_X )

        SUM_Ye = Zero
        SUM_Ym = Zero
        SUM_Ef = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef )
#endif
        DO iS   = 1, nSpecies
        DO iN_E = 1, nE_G

          IF ( iS <= iNuE_Bar ) THEN
          SUM_Ye = SUM_Ye + Dnu(iN_E,iS,iN_X) * W2_S(iN_E) * LeptonNumber(iS)
          END IF
          IF ( iS == iNuM .OR. iS == iNuM_Bar ) THEN
            SUM_Ym = SUM_Ym + Dnu(iN_E,iS,iN_X) * W2_S(iN_E) * LeptonNumber(iS)
          END IF

          SUM_Ef = SUM_Ef + Dnu(iN_E,iS,iN_X) * W3_S(iN_E)

        END DO
        END DO

        G_Ye(iN_X) = C_Ye(iN_X) - wYe * SUM_Ye * S_Ye(iN_X)
        G_Ym(iN_X) = C_Ym(iN_X) - wYm * SUM_Ym * S_Ym(iN_X)
        G_Ef(iN_X) = C_Ef(iN_X) - wEf * SUM_Ef * S_Ef(iN_X)

        Gm(iYe,iN_X) = G_Ye(iN_X)
        Gm(iYm,iN_X) = G_Ym(iN_X)
        Gm(iEf,iN_X) = G_Ef(iN_X)

        Fm(iYe,iN_X) = G_Ye(iN_X) - U_Ye(iN_X)
        Fm(iYm,iN_X) = G_Ym(iN_X) - U_Ym(iN_X)
        Fm(iEf,iN_X) = G_Ef(iN_X) - U_Ef(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_OrderOne


  PURE SUBROUTINE MuonMaterialNumberSource_Point &
    ( Weights, Chi, Equilibrium, Distribution, Transfer, &
      ActiveCC, ActiveTransfer, Source, AbsoluteSource )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), DIMENSION(:),     INTENT(in)  :: Weights
    REAL(DP), DIMENSION(:,:),   INTENT(in)  :: Chi, Equilibrium, Distribution
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Transfer
    LOGICAL,                   INTENT(in)  :: ActiveCC
    LOGICAL,  DIMENSION(3),     INTENT(in)  :: ActiveTransfer
    REAL(DP),                  INTENT(out) :: Source, AbsoluteSource

    INTEGER  :: iN_E, iS, iP
    REAL(DP) :: Total, Correction, Rate, SignedWeight

    ! --- Material Muon-Number Source ---
    ! Species 1 is nu_mu, species 2 anti-nu_mu. Number-conserving
    ! reactions are excluded from this material source by construction.
    ! Compensated summation requires no unsafe floating-point reassociation.

    Total          = Zero
    Correction     = Zero
    AbsoluteSource = Zero

    DO iN_E = 1, SIZE( Weights )
      DO iS = 1, 2

        SignedWeight = REAL( 2 * iS - 3, DP ) * Weights(iN_E)

        IF( ActiveCC )THEN
          Rate = ( Chi(iN_E,iS) &
                     * ( Equilibrium(iN_E,iS) - Distribution(iN_E,iS) ) ) &
                   * SignedWeight
          CALL AddMuonSourceCompensated( Rate, Total, Correction )
          AbsoluteSource = AbsoluteSource + ABS( Rate )
        END IF

        DO iP = 1, 3
          IF( .NOT. ActiveTransfer(iP) ) CYCLE
          Rate = SignedWeight * Transfer(iN_E,iS,iP)
          CALL AddMuonSourceCompensated( Rate, Total, Correction )
          AbsoluteSource = AbsoluteSource + ABS( Rate )
        END DO

      END DO
    END DO

    Source = Total + Correction

  END SUBROUTINE MuonMaterialNumberSource_Point


  PURE SUBROUTINE AddMuonSourceCompensated( Value, Total, Correction )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: Value
    REAL(DP), INTENT(inout) :: Total, Correction

    REAL(DP) :: Next

    Next = Total + Value
    IF( ABS( Total ) >= ABS( Value ) )THEN
      Correction = Correction + ( ( Total - Next ) + Value )
    ELSE
      Correction = Correction + ( ( Value - Next ) + Total )
    END IF
    Total = Next

  END SUBROUTINE AddMuonSourceCompensated


  SUBROUTINE MuonMaterialSourceAtPoint_OrderV &
    ( iN_X, Dnu, TransferRequested, Source, AbsoluteSource )

    ! The implemented Order-v number equation is
    ! C_Dnu - Dnu - v_dot_I + dt*R_N = 0. Its MATERIAL number source
    ! is therefore -sum(W2*lepton_sign*R_N), as in the Q_Lm diagnostic;
    ! no additional v_dot_R_G transformation belongs in this update.
#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,                   INTENT(in)  :: iN_X
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu
    LOGICAL,  DIMENSION(3),     INTENT(in)  :: TransferRequested
    REAL(DP),                  INTENT(out) :: Source, AbsoluteSource

    LOGICAL, DIMENSION(3) :: TransferActive

    Source         = Zero
    AbsoluteSource = Zero
    IF( nSpecies < iNuM_Bar ) RETURN
    IF( .NOT. MatterSource_Mask(iYm,iN_X) ) RETURN

    TransferActive(1) = TransferRequested(1) .AND. Opacity_Mask(iOp_LFE,iN_X)
    TransferActive(2) = TransferRequested(2) .AND. Opacity_Mask(iOp_LFC,iN_X)
    TransferActive(3) = TransferRequested(3) .AND. Opacity_Mask(iOp_IMD,iN_X)

    CALL MuonMaterialNumberSource_Point &
           ( W2_S, Chi_EmAb(:,iNuM:iNuM_Bar,iN_X), &
             Dnu_0(:,iNuM:iNuM_Bar,iN_X), Dnu(:,iNuM:iNuM_Bar,iN_X), &
             R_GL(:,iNuM:iNuM_Bar,iN_X,iGL_LFE:iGL_IMD), &
             Opacity_Mask(iOp_EmAb_Muon,iN_X), TransferActive, Source, AbsoluteSource )

  END SUBROUTINE MuonMaterialSourceAtPoint_OrderV


  SUBROUTINE ComputeMuonMaterialSource_OrderV &
    ( MASK, Dnu, Source, AbsoluteSource )

    LOGICAL,  DIMENSION(:),     INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: Dnu
    REAL(DP), DIMENSION(:),     INTENT(out) :: Source, AbsoluteSource

    INTEGER :: iN_X
    LOGICAL, DIMENSION(3) :: TransferRequested

    TransferRequested = [ Include_LFE, Include_LFC, Include_IMD ]

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO FIRSTPRIVATE( TransferRequested )
#elif defined( THORNADO_OACC )
    !$ACC PARALLEL LOOP GANG VECTOR FIRSTPRIVATE( TransferRequested ) &
    !$ACC PRESENT( MASK, Dnu, Source, AbsoluteSource )
#elif defined( THORNADO_OMP )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      Source(iN_X)         = Zero
      AbsoluteSource(iN_X) = Zero
      IF( MASK(iN_X) ) &
        CALL MuonMaterialSourceAtPoint_OrderV &
               ( iN_X, Dnu, TransferRequested, Source(iN_X), AbsoluteSource(iN_X) )
    END DO

  END SUBROUTINE ComputeMuonMaterialSource_OrderV

  SUBROUTINE ComputeMatterRHS_OrderV &
    ( MASK, Fm, Gm, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, dt_Option, MuonLimiterCorrection_Option )

    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP),                   INTENT(in), OPTIONAL :: dt_Option
    REAL(DP), DIMENSION(:),      INTENT(in), OPTIONAL :: MuonLimiterCorrection_Option

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3
    REAL(DP) :: SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3
    REAL(DP) :: wYe, wYm, wEf, wV1, wV2, wV3
    REAL(DP) :: SourceDt, MuonSource, MuonAbsoluteSource, MuonIncrement
    LOGICAL, DIMENSION(3) :: TransferRequested
    INTEGER :: ierr

    SourceDt = Zero
    TransferRequested = [ Include_LFE, Include_LFC, Include_IMD ]
    IF( UseMuonLeptonSource )THEN
      IF( .NOT. PRESENT( dt_Option ) .OR. .NOT. PRESENT( MuonLimiterCorrection_Option ) )THEN
        WRITE(ERROR_UNIT,'(A)') 'Missing stage duration/limiter correction for direct muon source.'
        CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
        RETURN
      END IF
      SourceDt = dt_Option
    END IF

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP FIRSTPRIVATE( UseMuonLeptonSource, SourceDt, TransferRequested ) &
    !$OMP PRIVATE( MuonSource, MuonAbsoluteSource, MuonIncrement, &
    !$OMP          SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, &
    !$OMP          wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC FIRSTPRIVATE( UseMuonLeptonSource, SourceDt, TransferRequested ) &
    !$ACC PRIVATE( MuonSource, MuonAbsoluteSource, MuonIncrement, &
    !$ACC          SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          V_d_1, V_d_2, V_d_3, &
    !$ACC          wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( MuonSource, MuonAbsoluteSource, MuonIncrement, &
    !$OMP          SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, &
    !$OMP          wYe, wYm, wEf, wV1, wV2, wV3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        wYe = MatterSourceWeight( iYe, iN_X )
        wYm = MatterSourceWeight( iYm, iN_X )
        wEf = MatterSourceWeight( iEf, iN_X )
        wV1 = MatterSourceWeight( iV1, iN_X )
        wV2 = MatterSourceWeight( iV2, iN_X )
        wV3 = MatterSourceWeight( iV3, iN_X )

        V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)

        SUM_Ye = Zero
        SUM_Ym = Zero
        SUM_Ef = Zero
        SUM_V1 = Zero
        SUM_V2 = Zero
        SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$OMP REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$ACC REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
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
            SUM_Ye = SUM_Ye + Nnu   * W2_S(iN_E) * LeptonNumber(iS)
          END IF
          IF ( iS == iNuM .OR. iS == iNuM_Bar ) THEN
            SUM_Ym = SUM_Ym + Nnu   * W2_S(iN_E) * LeptonNumber(iS)
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

        Fm(iYe,iN_X) = C_Ye(iN_X) &
          - wYe * SUM_Ye * S_Ye(iN_X) &
          - ( U_Ye(iN_X) - One )

        Fm(iYm,iN_X) = C_Ym(iN_X) &
          - wYm * SUM_Ym * S_Ym(iN_X) &
          - ( U_Ym(iN_X) - Ym_old(iN_X) )

        IF( UseMuonLeptonSource )THEN
          MuonIncrement = Zero
          IF( wYm /= Zero )THEN
            CALL MuonMaterialSourceAtPoint_OrderV &
                   ( iN_X, Dnu, TransferRequested, MuonSource, MuonAbsoluteSource )
            MuonIncrement = wYm * S_Ym(iN_X) &
              * ( SourceDt * MuonSource + MuonLimiterCorrection_Option(iN_X) )
          END IF
          Fm(iYm,iN_X) = MuonIncrement - ( U_Ym(iN_X) - Ym_old(iN_X) )
        END IF

        Fm(iEf,iN_X) = C_Ef(iN_X) &
          - wEf * SUM_Ef * S_Ef(iN_X) &
          - ( U_Ef(iN_X) - One )

        Fm(iV1,iN_X) = C_V_d_1(iN_X) &
          - wV1 * SUM_V1 * S_V_d_1(iN_X) &
          - ( U_V_d_1(iN_X) &
              - Gm_dd_11(iN_X) * V_u_1_old(iN_X) / SpeedOfLight )

        Fm(iV2,iN_X) = C_V_d_2(iN_X) &
          - wV2 * SUM_V2 * S_V_d_2(iN_X) &
          - ( U_V_d_2(iN_X) &
              - Gm_dd_22(iN_X) * V_u_2_old(iN_X) / SpeedOfLight )

        Fm(iV3,iN_X) = C_V_d_3(iN_X) &
          - wV3 * SUM_V3 * S_V_d_3(iN_X) &
          - ( U_V_d_3(iN_X) &
              - Gm_dd_33(iN_X) * V_u_3_old(iN_X) / SpeedOfLight )

        G_Ye   (iN_X) = U_Ye(iN_X) + Fm(iYe,iN_X)
        G_Ym   (iN_X) = U_Ym(iN_X) + Fm(iYm,iN_X)
        G_Ef   (iN_X) = U_Ef   (iN_X) + Fm(iEf,iN_X)
        G_V_d_1(iN_X) = U_V_d_1(iN_X) + Fm(iV1,iN_X)
        G_V_d_2(iN_X) = U_V_d_2(iN_X) + Fm(iV2,iN_X)
        G_V_d_3(iN_X) = U_V_d_3(iN_X) + Fm(iV3,iN_X)

        Gm(iYe,iN_X) = G_Ye   (iN_X)
        Gm(iYm,iN_X) = G_Ym   (iN_X)
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
    REAL(DP) :: SUM_Ef, SUM_V1, SUM_V2, SUM_V3, SUM_Ye, SUM_Ym
    REAL(DP) :: wD, wYe, wYm, wEf, wV1, wV2, wV3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotV, W, B_d_1, B_d_2, B_d_3, DT, h, &
    !$OMP          wD, wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$ACC          V_d_1, V_d_2, V_d_3, vDotV, W, B_d_1, B_d_2, B_d_3, DT, h, &
    !$ACC          wD, wYe, wYm, wEf, wV1, wV2, wV3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          V_d_1, V_d_2, V_d_3, vDotV, W, B_d_1, B_d_2, B_d_3, DT, h, &
    !$OMP          wD, wYe, wYm, wEf, wV1, wV2, wV3, &
    !$OMP          Inu_d_1, Inu_d_2, Inu_d_3, k_dd, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 )
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        wD  = MatterSourceWeight( iD,  iN_X )
        wYe = MatterSourceWeight( iYe, iN_X )
        wYm = MatterSourceWeight( iYm, iN_X )
        wEf = MatterSourceWeight( iEf, iN_X )
        wV1 = MatterSourceWeight( iV1, iN_X )
        wV2 = MatterSourceWeight( iV2, iN_X )
        wV3 = MatterSourceWeight( iV3, iN_X )

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
     
        SUM_Ye  = Zero
        SUM_Ym  = Zero
        SUM_Ef  = Zero
        SUM_V1  = Zero
        SUM_V2  = Zero
        SUM_V3  = Zero

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Inu_d_1, Inu_d_2, Inu_d_3, k_dd, &
    !$OMP          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$OMP          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$OMP REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Inu_d_1, Inu_d_2, Inu_d_3, k_dd, &
    !$ACC          vDotInu, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, &
    !$ACC          Nnu, Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3 ) &
    !$ACC REDUCTION( + : SUM_Ye, SUM_Ym, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
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
            SUM_Ye = SUM_Ye + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
          END IF
          IF ( iS == iNuM .OR. iS == iNuM_Bar ) THEN
            SUM_Ym = SUM_Ym + Nnu     * W2_S(iN_E) * LeptonNumber(iS)
          END IF

          SUM_Ef  = SUM_Ef  + Enu     * W3_S(iN_E)
          SUM_V1  = SUM_V1  + Fnu_d_1 * W3_S(iN_E)
          SUM_V2  = SUM_V2  + Fnu_d_2 * W3_S(iN_E)
          SUM_V3  = SUM_V3  + Fnu_d_3 * W3_S(iN_E)

        END DO
        END DO

        h = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) 

        IF ( wD .EQ. One ) THEN
          G_D    (iN_X)  = ( 1.0_DP / W ) * C_D  (iN_X) 
        ELSE  
          G_D    (iN_X) &
            = C_D(iN_X)
        END IF

        IF ( wYe .EQ. One ) THEN
          G_Ye   (iN_X) &
            = C_Ye   (iN_X) - ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Ye * S_Ye(iN_X)
        ELSE
          G_Ye    (iN_X) &
            = C_Ye(iN_X)
        END IF

        IF ( wYm .EQ. One ) THEN
          G_Ym    (iN_X) &
            = C_Ym(iN_X) - ( AtomicMassUnit / cD_old(iN_X) ) * SUM_Ym * S_Ym(iN_X)
        ELSE
          G_Ym    (iN_X) &
            = C_Ym(iN_X)
        END IF
         
        IF ( wEf .EQ. One ) THEN
          G_Ef   (iN_X)  = S_Ef(iN_X) * ( 1.0_DP - W ) / W**2 * ( W + ( W + 1.0_DP) * P(iN_X) / D(iN_X) ) &
                         + 1.0_DP / W * C_Ef(iN_X) - S_Ef(iN_X) / ( W * cD_old(iN_X) ) * SUM_Ef 
        ELSE
          G_Ef   (iN_X) &
            = C_Ef(iN_X)
        END IF

        IF ( wV1 .EQ. One ) THEN
          G_V_d_1(iN_X)  &
            = 1.0_DP / ( h * W ) * (  C_V_d_1(iN_X) - SUM_V1 * S_V_d_1(iN_X) / cD_old(iN_X) ) 
        ELSE
          G_V_d_1(iN_X) &
            = C_V_d_1(iN_X)
        END IF


        IF ( wV2 .EQ. One ) THEN
          G_V_d_2(iN_X) &
            = 1.0_DP / ( h * W ) * (  C_V_d_2(iN_X) - SUM_V2 * S_V_d_2(iN_X) / cD_old(iN_X) ) 
        ELSE
          G_V_d_2(iN_X) &
            = C_V_d_2(iN_X)
        END IF


        IF ( wV3 .EQ. One ) THEN
          G_V_d_3(iN_X) &
            = 1.0_DP / ( h * W ) * (  C_V_d_3(iN_X) - SUM_V3 * S_V_d_3(iN_X) / cD_old(iN_X) ) 
        ELSE
          G_V_d_3(iN_X) &
            = C_V_d_3(iN_X)
        END IF


        

        IF( MoveLeft .EQ. 1) THEN
          IF ( wEf .EQ. One ) THEN
            G_Ef(iN_X) = ( 1.0_DP / W ) * C_Ef(iN_X) - S_Ef(iN_X) / ( W * cD_old(iN_X) ) * SUM_Ef &
                       + S_Ef(iN_X) * ( 1.0_DP - W ) / W**2 *  ( W + 1.0_DP) * P(iN_X) / D(iN_X) 
          ELSE
            G_Ef   (iN_X) &
              = C_Ef(iN_X)
          END IF

        END IF

        Gm(iD ,iN_X) = G_D    (iN_X)
        Gm(iYe,iN_X) = G_Ye   (iN_X)
        Gm(iYm,iN_X) = G_Ym   (iN_X)
        Gm(iEf,iN_X) = G_Ef   (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

        Fm(iD ,iN_X) = G_D    (iN_X) - U_D    (iN_X)
        Fm(iYe,iN_X) = G_Ye   (iN_X) - U_Ye   (iN_X)
        Fm(iYm,iN_X) = G_Ym   (iN_X) - U_Ym   (iN_X)
        Fm(iEf,iN_X) = G_Ef   (iN_X) - U_Ef   (iN_X)
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
      Gm_dd_11, Gm_dd_22, Gm_dd_33, NumberSource_Option, NumberFactor_Option, &
      NumberFraction_Option )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: NumberSource_Option, NumberFactor_Option
    REAL(DP), DIMENSION(:,:,:), INTENT(out), OPTIONAL :: NumberFraction_Option

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
                       + Chi_Brem  (iN_E,iS,iN_X) &
                       + Chi_GL(iN_E,iS,iN_X,iGL_NMS) &
                       + Chi_GL(iN_E,iS,iN_X,iGL_LFE) &
                       + Chi_GL(iN_E,iS,iN_X,iGL_LFC) &
                       + Chi_GL(iN_E,iS,iN_X,iGL_IMD)

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
        ! EmAb vanishes identically when Dnu == Dnu_0.  The inelastic and
        ! pair rates satisfy physical detailed balance and vanish at the
        ! thermodynamic Fermi--Dirac state in exact arithmetic; direct
        ! floating-point gain--loss evaluation can leave roundoff residuals.

        R_N = Chi_EmAb(iN_E,iS,iN_X) &
                * ( Dnu_0(iN_E,iS,iN_X) - Dnu(iN_E,iS,iN_X) ) &
              + R_NNS(iN_E,iS,iN_X) &
              + R_NES(iN_E,iS,iN_X) &
              + R_Pair(iN_E,iS,iN_X) &
              + R_NuPair(iN_E,iS,iN_X) &
              + R_Brem(iN_E,iS,iN_X) &
              + R_GL(iN_E,iS,iN_X,iGL_NMS) &
              + R_GL(iN_E,iS,iN_X,iGL_LFE) &
              + R_GL(iN_E,iS,iN_X,iGL_LFC) &
              + R_GL(iN_E,iS,iN_X,iGL_IMD) &
              + L_N

        IF( PRESENT(NumberSource_Option) ) NumberSource_Option(iN_E,iS,iN_X) = R_N
        IF( PRESENT(NumberFactor_Option) ) NumberFactor_Option(iN_E,iS,iN_X) &
          = One/Omega(iN_X)+dt*Chi_T
        IF( PRESENT(NumberFraction_Option) ) NumberFraction_Option(iN_E,iS,iN_X) &
          = Omega(iN_X)*dt*Chi_T/(One+Omega(iN_X)*dt*Chi_T)

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
    ( MASK, Fm, Gm, Ye, Ym, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: Ye, Ym, E
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

        Fm(iYe,iN_X) = Gm(iYe,iN_X) - U_Ye(iN_X)
        Fm(iYm,iN_X) = Gm(iYm,iN_X) - U_Ym(iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef(iN_X)

        U_Ye  (iN_X) = Gm(iYe,iN_X)
        U_Ym  (iN_X) = Gm(iYm,iN_X)
        U_Ef  (iN_X) = Gm(iEf,iN_X)

        Ye    (iN_X) = U_Ye(iN_X) * Ye_old(iN_X)
        Ym    (iN_X) = U_Ym(iN_X) * Ym_old(iN_X)
        E     (iN_X) = U_Ef(iN_X) * Ef_old(iN_X)

      END IF

    END DO

  END SUBROUTINE UpdateMatterRHS_OrderOne


  SUBROUTINE UpdateMatterRHS_OrderV &
    ( MASK, Fm, Gm, Ye, Ym, E, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: Ye, Ym, E, V_u_1, V_u_2, V_u_3
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

        Fm(iYe,iN_X) = Gm(iYe,iN_X) - U_Ye   (iN_X)
        Fm(iYm,iN_X) = Gm(iYm,iN_X) - U_Ym   (iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)

        U_Ye   (iN_X) = Gm(iYe,iN_X)
        U_Ym   (iN_X) = Gm(iYm,iN_X)
        U_Ef   (iN_X) = Gm(iEf,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        Ye(iN_X) = U_Ye   (iN_X) * Ye_old(iN_X)
        Ym(iN_X) = U_Ym(iN_X)
        Ef       = U_Ef   (iN_X) * Ef_old(iN_X)
        V_d_1    = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2    = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3    = U_V_d_3(iN_X) * SpeedOfLight

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
    ( MASK, Fm, Gm, D, Ye, Ym, E, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: D, Ye, Ym, E, V_u_1, V_u_2, V_u_3
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
        Fm(iYe,iN_X) = Gm(iYe,iN_X) - U_Ye   (iN_X)
        Fm(iYm,iN_X) = Gm(iYm,iN_X) - U_Ym   (iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)

        U_D    (iN_X) = Gm(iD ,iN_X)
        U_Ye   (iN_X) = Gm(iYe,iN_X)
        U_Ym   (iN_X) = Gm(iYm,iN_X)
        U_Ef   (iN_X) = Gm(iEf,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        D (iN_X) = U_D    (iN_X) * D_old (iN_X)
        Ye(iN_X) = U_Ye   (iN_X) * Ye_old(iN_X)
        Ym(iN_X) = U_Ym   (iN_X) * Ym_old(iN_X)
        Ef_temp  = U_Ef   (iN_X) * Ef_old(iN_X)
        V_d_1    = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2    = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3    = U_V_d_3(iN_X) * SpeedOfLight

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


  SUBROUTINE SaveRadiationRecoveryMatterState_OrderV &
    ( MASK, T, Ye, Ym, E, V1, V2, V3, &
      T_Base, Ye_Base, Ym_Base, E_Base, V1_Base, V2_Base, V3_Base )

    LOGICAL,  DIMENSION(:), INTENT(in)  :: MASK
    REAL(DP), DIMENSION(:), INTENT(in)  :: T, Ye, Ym, E, V1, V2, V3
    REAL(DP), DIMENSION(:), INTENT(out) :: &
      T_Base, Ye_Base, Ym_Base, E_Base, V1_Base, V2_Base, V3_Base

    INTEGER :: iN_X

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF( MASK(iN_X) )THEN
        T_Base (iN_X) = T (iN_X)
        Ye_Base(iN_X) = Ye(iN_X)
        Ym_Base(iN_X) = Ym(iN_X)
        E_Base (iN_X) = E (iN_X)
        V1_Base(iN_X) = V1(iN_X)
        V2_Base(iN_X) = V2(iN_X)
        V3_Base(iN_X) = V3(iN_X)
      END IF
    END DO

  END SUBROUTINE SaveRadiationRecoveryMatterState_OrderV


  SUBROUTINE RestoreRadiationRecoveryMatterState_OrderV &
    ( MASK, T_Base, Ye_Base, Ym_Base, E_Base, V1_Base, V2_Base, V3_Base, &
      T, Ye, Ym, E, V1, V2, V3 )

    LOGICAL,  DIMENSION(:), INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:), INTENT(in)    :: &
      T_Base, Ye_Base, Ym_Base, E_Base, V1_Base, V2_Base, V3_Base
    REAL(DP), DIMENSION(:), INTENT(inout) :: T, Ye, Ym, E, V1, V2, V3

    INTEGER :: iN_X

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF( MASK(iN_X) )THEN
        T (iN_X) = T_Base (iN_X)
        Ye(iN_X) = Ye_Base(iN_X)
        Ym(iN_X) = Ym_Base(iN_X)
        E (iN_X) = E_Base (iN_X)
        V1(iN_X) = V1_Base(iN_X)
        V2(iN_X) = V2_Base(iN_X)
        V3(iN_X) = V3_Base(iN_X)
      END IF
    END DO

  END SUBROUTINE RestoreRadiationRecoveryMatterState_OrderV


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
    !$OMP FIRSTPRIVATE(UseMuonLeptonSource) &
    !$OMP PRIVATE( CONVERGED, Fnorm )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC FIRSTPRIVATE(UseMuonLeptonSource) &
    !$ACC PRIVATE( CONVERGED, Fnorm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, Fnorm )
#endif
    DO iN_X = 1, nX_G
      IF( MASK_OUTER(iN_X) )THEN

        Fnorm = MAXVAL( ABS( Fm(:,iN_X) ) )

#if defined(TWOMOMENT_ORDER_V)
        ! An absolute 1e-8 test on an additive Ym would accept a tiny but
        ! order-unity relative composition step. Only the opt-in source
        ! formulation changes this row's scaling; other fields are unchanged.
        IF( UseMuonLeptonSource ) Fnorm = MAX(Fnorm,ABS(Fm(iYm,iN_X)) &
          /MAX(ABS(U_Ym(iN_X))+ABS(G_Ym(iN_X))+ABS(Ym_old(iN_X)),SqrtTiny))
#endif

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


#if defined( TWOMOMENT_ORDER_V )
  SUBROUTINE CheckMuonCouplingIteration_OrderV &
    ( Stage, dt, k_outer, k_inner, MASK, &
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Error_Option, CheckOpacity_Option, CheckRates_Option, &
      Fm_Option, Gm_Option, ForceFailure_Option )

    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE mpi

    CHARACTER(*),                  INTENT(in) :: Stage
    REAL(DP),                      INTENT(in) :: dt
    INTEGER,                       INTENT(in) :: k_outer, k_inner
    LOGICAL,  DIMENSION(:),        INTENT(in) :: MASK
    REAL(DP), DIMENSION(:),        INTENT(in) :: &
      D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),        INTENT(in) :: &
      Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:,:,:),    INTENT(in) :: &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    INTEGER,  DIMENSION(:),        INTENT(in), OPTIONAL :: Error_Option
    LOGICAL,                       INTENT(in), OPTIONAL :: &
      CheckOpacity_Option, CheckRates_Option
    REAL(DP), DIMENSION(:,:),      INTENT(in), OPTIONAL :: &
      Fm_Option, Gm_Option
    LOGICAL,                       INTENT(in), OPTIONAL :: &
      ForceFailure_Option

    INTEGER  :: ierr, myid, FailureRank, CandidateRank
    INTEGER  :: iN_E, iN_X, iS, iV, iFile
    INTEGER  :: FailureEnergy, FailurePoint, FailureSpecies
    INTEGER  :: iOS
    LOGICAL  :: CheckOpacity, CheckRates, Failed
    REAL(DP) :: Mumu, MuonLeptonCurrent, MuonLeptonOld, vDotI
    REAL(DP) :: Dnu0Value, ChiValue, Fm_N, Gm_N
    REAL(DP) :: RNNSValue, RNESValue, RPairValue, RNuPairValue, RBremValue
    REAL(DP) :: ResidualMax, MapMax, RelativeResidualMax
    REAL(DP) :: Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3
    REAL(DP) :: Dnorm, SpeciesRelativeResidual
    CHARACTER(64)  :: FailureQuantity
    CHARACTER(640) :: FileName

    CheckOpacity = .FALSE.
    CheckRates   = .FALSE.
    IF( PRESENT( CheckOpacity_Option ) ) &
      CheckOpacity = CheckOpacity_Option
    IF( PRESENT( CheckRates_Option ) ) &
      CheckRates = CheckRates_Option

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Dnu_old, &
    !$OMP   Inu_u_1_old, Inu_u_2_old, Inu_u_3_old, &
    !$OMP   V_u_1_old, V_u_2_old, V_u_3_old, Ym_old, Opacity_Mask )
    IF( CheckOpacity )THEN
      !$OMP TARGET UPDATE FROM( Dnu_0, Chi_EmAb )
    END IF
    IF( CheckRates )THEN
      !$OMP TARGET UPDATE FROM &
      !$OMP ( R_NNS, R_NES, R_Pair, R_NuPair, R_Brem )
    END IF
    IF( PRESENT( Fm_Option ) )THEN
      !$OMP TARGET UPDATE FROM( Fm_Option )
    END IF
    IF( PRESENT( Gm_Option ) )THEN
      !$OMP TARGET UPDATE FROM( Gm_Option )
    END IF
    IF( PRESENT( Error_Option ) )THEN
      !$OMP TARGET UPDATE FROM( Error_Option )
    END IF
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( D, T, Ye, Ym, E, V_u_1, V_u_2, V_u_3, &
    !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, Dnu_old, &
    !$ACC   Inu_u_1_old, Inu_u_2_old, Inu_u_3_old, &
    !$ACC   V_u_1_old, V_u_2_old, V_u_3_old, Ym_old, Opacity_Mask )
    IF( CheckOpacity )THEN
      !$ACC UPDATE HOST( Dnu_0, Chi_EmAb )
    END IF
    IF( CheckRates )THEN
      !$ACC UPDATE HOST &
      !$ACC ( R_NNS, R_NES, R_Pair, R_NuPair, R_Brem )
    END IF
    IF( PRESENT( Fm_Option ) )THEN
      !$ACC UPDATE HOST( Fm_Option )
    END IF
    IF( PRESENT( Gm_Option ) )THEN
      !$ACC UPDATE HOST( Gm_Option )
    END IF
    IF( PRESENT( Error_Option ) )THEN
      !$ACC UPDATE HOST( Error_Option )
    END IF
#endif

    Failed          = .FALSE.
    FailurePoint    = 0
    FailureEnergy   = 0
    FailureSpecies  = 0
    FailureQuantity = ''

    PointLoop: DO iN_X = 1, nX_G

      IF( .NOT. MASK(iN_X) ) CYCLE

      IF( PRESENT( ForceFailure_Option ) )THEN
        IF( ForceFailure_Option )THEN
          Failed = .TRUE.
          FailureQuantity = 'iteration limit reached'
        END IF
      END IF

      IF( PRESENT( Error_Option ) )THEN
        IF( Error_Option(iN_X) /= 0 )THEN
          Failed = .TRUE.
          FailureQuantity = 'matter-recovery error flag'
        END IF
      END IF

      IF( .NOT. Failed )THEN
        IF( .NOT. IsFiniteDebugValue(D(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'D'
        ELSE IF( .NOT. IsFiniteDebugValue(T(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'T'
        ELSE IF( .NOT. IsFiniteDebugValue(Ye(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'Ye'
        ELSE IF( .NOT. IsFiniteDebugValue(Ym(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'Ym'
        ELSE IF( .NOT. IsFiniteDebugValue(E(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'E'
        ELSE IF( .NOT. IsFiniteDebugValue(V_u_1(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'V_u_1'
        ELSE IF( .NOT. IsFiniteDebugValue(V_u_2(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'V_u_2'
        ELSE IF( .NOT. IsFiniteDebugValue(V_u_3(iN_X)) )THEN
          Failed = .TRUE.; FailureQuantity = 'V_u_3'
        END IF
      END IF

      IF( .NOT. Failed )THEN
        SpeciesLoop: DO iS = 1, nSpecies
          DO iN_E = 1, nE_G
            IF( .NOT. IsFiniteDebugValue(Dnu(iN_E,iS,iN_X)) )THEN
              Failed = .TRUE.; FailureQuantity = 'Dnu'
            ELSE IF( .NOT. IsFiniteDebugValue(Inu_u_1(iN_E,iS,iN_X)) )THEN
              Failed = .TRUE.; FailureQuantity = 'Inu_u_1'
            ELSE IF( .NOT. IsFiniteDebugValue(Inu_u_2(iN_E,iS,iN_X)) )THEN
              Failed = .TRUE.; FailureQuantity = 'Inu_u_2'
            ELSE IF( .NOT. IsFiniteDebugValue(Inu_u_3(iN_E,iS,iN_X)) )THEN
              Failed = .TRUE.; FailureQuantity = 'Inu_u_3'
            ELSE IF( CheckOpacity )THEN
              IF( .NOT. IsFiniteDebugValue(Dnu_0(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'Dnu_0'
              ELSE IF( .NOT. IsFiniteDebugValue(Chi_EmAb(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'Chi_EmAb'
              ELSE IF( Chi_EmAb(iN_E,iS,iN_X) < Zero )THEN
                Failed = .TRUE.; FailureQuantity = 'negative Chi_EmAb'
              END IF
            END IF
            IF( .NOT. Failed .AND. CheckRates )THEN
              IF( .NOT. IsFiniteDebugValue(R_NNS(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'R_NNS'
              ELSE IF( .NOT. IsFiniteDebugValue(R_NES(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'R_NES'
              ELSE IF( .NOT. IsFiniteDebugValue(R_Pair(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'R_Pair'
              ELSE IF( .NOT. IsFiniteDebugValue(R_NuPair(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'R_NuPair'
              ELSE IF( .NOT. IsFiniteDebugValue(R_Brem(iN_E,iS,iN_X)) )THEN
                Failed = .TRUE.; FailureQuantity = 'R_Brem'
              END IF
            END IF
            IF( Failed )THEN
              FailureEnergy  = iN_E
              FailureSpecies = iS
              EXIT SpeciesLoop
            END IF
          END DO
        END DO SpeciesLoop
      END IF

      IF( .NOT. Failed .AND. PRESENT(Fm_Option) )THEN
        DO iV = 1, SIZE(Fm_Option,1)
          IF( .NOT. IsFiniteDebugValue(Fm_Option(iV,iN_X)) )THEN
            Failed = .TRUE.
            FailureQuantity = 'fixed-point residual'
            EXIT
          END IF
        END DO
      END IF

      IF( .NOT. Failed .AND. PRESENT(Gm_Option) )THEN
        DO iV = 1, SIZE(Gm_Option,1)
          IF( .NOT. IsFiniteDebugValue(Gm_Option(iV,iN_X)) )THEN
            Failed = .TRUE.
            FailureQuantity = 'fixed-point map'
            EXIT
          END IF
        END DO
      END IF

      IF( Failed )THEN
        FailurePoint = iN_X
        EXIT PointLoop
      END IF

    END DO PointLoop

    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    CandidateRank = HUGE(CandidateRank)
    IF( Failed ) CandidateRank = myid
    CALL MPI_ALLREDUCE &
           ( CandidateRank, FailureRank, 1, MPI_INTEGER, MPI_MIN, &
             MPI_COMM_WORLD, ierr )

    IF( FailureRank == HUGE(FailureRank) ) RETURN

    IF( myid == FailureRank )THEN

      Mumu = Zero
      IF(      IsFiniteDebugValue(D(FailurePoint)) &
          .AND. IsFiniteDebugValue(T(FailurePoint)) &
          .AND. IsFiniteDebugValue(Ye(FailurePoint)) &
          .AND. IsFiniteDebugValue(Ym(FailurePoint)) &
          .AND. D(FailurePoint) > Zero &
          .AND. T(FailurePoint) > Zero ) &
        CALL ComputeMuonChemicalPotential_TABLE &
               ( D(FailurePoint), T(FailurePoint), &
                 Ye(FailurePoint), Ym(FailurePoint), Mumu )

      MuonLeptonCurrent = Zero
      MuonLeptonOld     = Zero
      ResidualMax       = Zero
      MapMax            = Zero
      RelativeResidualMax = Zero
      IF( PRESENT(Fm_Option) ) &
        ResidualMax = MAXVAL(ABS(Fm_Option(:,FailurePoint)))
      IF( PRESENT(Gm_Option) ) &
        MapMax = MAXVAL(ABS(Gm_Option(:,FailurePoint)))
      IF( PRESENT(Fm_Option) )THEN
        DO iS = 1, nSpecies
          Fnorm_N  = Zero
          Fnorm_G1 = Zero
          Fnorm_G2 = Zero
          Fnorm_G3 = Zero
          DO iN_E = 1, nE_G
            iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR
            Fnorm_N = Fnorm_N &
              + W2_S(iN_E) * Fm_Option(iOS+iCR_N,FailurePoint)**2
            Fnorm_G1 = Fnorm_G1 &
              + W2_S(iN_E) * Fm_Option(iOS+iCR_G1,FailurePoint)**2
            Fnorm_G2 = Fnorm_G2 &
              + W2_S(iN_E) * Fm_Option(iOS+iCR_G2,FailurePoint)**2
            Fnorm_G3 = Fnorm_G3 &
              + W2_S(iN_E) * Fm_Option(iOS+iCR_G3,FailurePoint)**2
          END DO
          Fnorm = SQRT( MAX( Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 ) )
          Dnorm = WNORM( Dnu(:,iS,FailurePoint), W2_S )
          IF( Dnorm > Zero )THEN
            SpeciesRelativeResidual = Fnorm / Dnorm
          ELSE IF( Fnorm == Zero )THEN
            SpeciesRelativeResidual = Zero
          ELSE
            SpeciesRelativeResidual = HUGE(One)
          END IF
          RelativeResidualMax &
            = MAX( RelativeResidualMax, SpeciesRelativeResidual )
        END DO
      END IF
      IF( nSpecies >= iNuM_Bar )THEN
        DO iS = iNuM, iNuM_Bar
          DO iN_E = 1, nE_G
            vDotI = V_u_1(FailurePoint) * Gm_dd_11(FailurePoint) &
                      * Inu_u_1(iN_E,iS,FailurePoint) &
                    + V_u_2(FailurePoint) * Gm_dd_22(FailurePoint) &
                      * Inu_u_2(iN_E,iS,FailurePoint) &
                    + V_u_3(FailurePoint) * Gm_dd_33(FailurePoint) &
                      * Inu_u_3(iN_E,iS,FailurePoint)
            MuonLeptonCurrent = MuonLeptonCurrent &
              + W2_S(iN_E) * LeptonNumber(iS) &
                  * ( Dnu(iN_E,iS,FailurePoint) + vDotI )
            vDotI = V_u_1_old(FailurePoint) * Gm_dd_11(FailurePoint) &
                      * Inu_u_1_old(iN_E,iS,FailurePoint) &
                    + V_u_2_old(FailurePoint) * Gm_dd_22(FailurePoint) &
                      * Inu_u_2_old(iN_E,iS,FailurePoint) &
                    + V_u_3_old(FailurePoint) * Gm_dd_33(FailurePoint) &
                      * Inu_u_3_old(iN_E,iS,FailurePoint)
            MuonLeptonOld = MuonLeptonOld &
              + W2_S(iN_E) * LeptonNumber(iS) &
                  * ( Dnu_old(iN_E,iS,FailurePoint) + vDotI )
          END DO
        END DO
      END IF

      WRITE(ERROR_UNIT,'(/,A)') 'Muon-coupling iteration diagnostic'
      WRITE(ERROR_UNIT,'(A,A)')       '  failure stage             = ', TRIM(Stage)
      WRITE(ERROR_UNIT,'(A,A)')       '  failure condition         = ', TRIM(FailureQuantity)
      WRITE(ERROR_UNIT,'(A,I0)')      '  MPI rank                  = ', myid
      WRITE(ERROR_UNIT,'(A,I0)')      '  spatial point             = ', FailurePoint
      WRITE(ERROR_UNIT,'(A,2(I0,1X))')'  outer/inner iteration     = ', k_outer, k_inner
      IF( FailureEnergy > 0 .AND. FailureSpecies > 0 ) &
        WRITE(ERROR_UNIT,'(A,2(I0,1X))')'  energy/species index      = ', &
          FailureEnergy, FailureSpecies
      WRITE(ERROR_UNIT,'(A,ES26.16E3)') '  rho [g cm^-3]             = ', &
        D(FailurePoint) / Unit_D
      WRITE(ERROR_UNIT,'(A,ES26.16E3)') '  T [MeV]                   = ', &
        T(FailurePoint) / Unit_T
      WRITE(ERROR_UNIT,'(A,3ES26.16E3)')'  Ye, Ym, Ym_old            = ', &
        Ye(FailurePoint), Ym(FailurePoint), Ym_old(FailurePoint)
      WRITE(ERROR_UNIT,'(A,ES26.16E3)') '  Ym - Ym_old               = ', &
        Ym(FailurePoint) - Ym_old(FailurePoint)
      WRITE(ERROR_UNIT,'(A,ES26.16E3)') '  mu_mu [MeV]               = ', Mumu / MeV
      WRITE(ERROR_UNIT,'(A,ES26.16E3)') '  mu_mu/T                   = ', &
        Mumu / MAX(T(FailurePoint), TINY(One))
      WRITE(ERROR_UNIT,'(A,2ES26.16E3)')'  mu-family L current/old   = ', &
        MuonLeptonCurrent, MuonLeptonOld
      WRITE(ERROR_UNIT,'(A,2ES26.16E3)')'  max |residual|, |map|    = ', &
        ResidualMax, MapMax
      WRITE(ERROR_UNIT,'(A,ES26.16E3)') '  maximum relative residual = ', &
        RelativeResidualMax
      WRITE(ERROR_UNIT,'(A,L1)')      '  muon charged current mask = ', &
        Opacity_Mask(iOp_EmAb_Muon,FailurePoint)
      FLUSH(ERROR_UNIT)

      IF( TraceRadiationOnlySolve )THEN
        WRITE(FileName,'(A,".failure.rank",I6.6,".dat")') &
          TRIM(RadiationOnlyTraceFile), myid
      ELSE
        WRITE(FileName,'(A,I0,A)') &
          'muon_coupling_failure_rank', myid, '.dat'
      END IF
      OPEN( NEWUNIT = iFile, FILE = TRIM(FileName), STATUS = 'REPLACE', &
            ACTION = 'WRITE' )
      WRITE(iFile,'(A,A)') '# stage: ', TRIM(Stage)
      WRITE(iFile,'(A,A)') '# quantity: ', TRIM(FailureQuantity)
      WRITE(iFile,'(A,3(I0,1X))') '# rank point outer_iter: ', &
        myid, FailurePoint, k_outer
      WRITE(iFile,'(A,I0)') '# inner_iter: ', k_inner
      WRITE(iFile,'(A,6ES26.16E3)') '# rho T Ye Ym Ym_old mu_mu: ', &
        D(FailurePoint) / Unit_D, T(FailurePoint) / Unit_T, &
        Ye(FailurePoint), Ym(FailurePoint), Ym_old(FailurePoint), Mumu / MeV
      WRITE(iFile,'(A,4ES26.16E3)') '# Lmu_current Lmu_old residual_max map_max: ', &
        MuonLeptonCurrent, MuonLeptonOld, ResidualMax, MapMax
      WRITE(iFile,'(A,ES26.16E3)') '# relative_residual_max: ', &
        RelativeResidualMax
      WRITE(iFile,'(A)') &
        '# iE iS E_MeV Dnu I1 I2 I3 Dnu0 ChiEmAb dtChi EmAbNet ' // &
        'RNNS RNES RPair RNuPair RBrem Fmap Gmap'
      IF( nSpecies >= iNuM_Bar )THEN
        DO iS = iNuM, iNuM_Bar
          DO iN_E = 1, nE_G
            Fm_N = Zero
            Gm_N = Zero
            Dnu0Value = Zero
            ChiValue  = Zero
            RNNSValue   = Zero
            RNESValue   = Zero
            RPairValue  = Zero
            RNuPairValue = Zero
            RBremValue  = Zero
            iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR + iCR_N
            IF( PRESENT(Fm_Option) ) Fm_N = Fm_Option(iOS,FailurePoint)
            IF( PRESENT(Gm_Option) ) Gm_N = Gm_Option(iOS,FailurePoint)
            IF( CheckOpacity )THEN
              Dnu0Value = Dnu_0(iN_E,iS,FailurePoint)
              ChiValue  = Chi_EmAb(iN_E,iS,FailurePoint)
            END IF
            IF( CheckRates )THEN
              RNNSValue    = R_NNS(iN_E,iS,FailurePoint)
              RNESValue    = R_NES(iN_E,iS,FailurePoint)
              RPairValue   = R_Pair(iN_E,iS,FailurePoint)
              RNuPairValue = R_NuPair(iN_E,iS,FailurePoint)
              RBremValue   = R_Brem(iN_E,iS,FailurePoint)
            END IF
            WRITE(iFile,'(2I6,16ES26.16E3)') &
              iN_E, iS, E_N(iN_E) / MeV, &
              Dnu(iN_E,iS,FailurePoint), &
              Inu_u_1(iN_E,iS,FailurePoint), &
              Inu_u_2(iN_E,iS,FailurePoint), &
              Inu_u_3(iN_E,iS,FailurePoint), &
              Dnu0Value, ChiValue, dt * ChiValue, &
              ChiValue * ( Dnu0Value - Dnu(iN_E,iS,FailurePoint) ), &
              RNNSValue, RNESValue, RPairValue, RNuPairValue, RBremValue, &
              Fm_N, Gm_N
          END DO
        END DO
      END IF

      ! Re-enter the opacity stack only when this diagnostic was called
      ! after a successful opacity construction.  In particular, an
      ! inadmissible matter-recovery trial may contain T=0 or an
      ! EOS-inadmissible composition and must not be queried again here.
      IF( nSpecies >= iNuM_Bar .AND. CheckOpacity )THEN
        IF( IsFiniteDebugValue(D(FailurePoint)) &
            .AND. IsFiniteDebugValue(T(FailurePoint)) &
            .AND. IsFiniteDebugValue(Ye(FailurePoint)) &
            .AND. IsFiniteDebugValue(Ym(FailurePoint)) &
            .AND. D(FailurePoint) > Zero &
            .AND. T(FailurePoint) > Zero )THEN
          IF( PRESENT( Error_Option ) )THEN
            IF( Error_Option(FailurePoint) == 0 ) &
              CALL WriteMuonChargedCurrentOpacityDiagnostic &
                     ( iFile, E_N, D(FailurePoint), T(FailurePoint), &
                       Ye(FailurePoint), Ym(FailurePoint), &
                       Dnu_0(:,:,FailurePoint) )
          ELSE
            CALL WriteMuonChargedCurrentOpacityDiagnostic &
                   ( iFile, E_N, D(FailurePoint), T(FailurePoint), &
                     Ye(FailurePoint), Ym(FailurePoint), &
                     Dnu_0(:,:,FailurePoint) )
          END IF
        END IF
      END IF
      CLOSE(iFile)

    END IF

    CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )

  END SUBROUTINE CheckMuonCouplingIteration_OrderV


  LOGICAL FUNCTION IsFiniteDebugValue( X )

    REAL(DP), INTENT(in) :: X

    IsFiniteDebugValue = X == X .AND. ABS(X) <= HUGE(X)

  END FUNCTION IsFiniteDebugValue
#endif


  SUBROUTINE CheckError_Inner &
    ( MASK_inner, k_outer, k_inner, &
      D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    USE mpi

    LOGICAL,  DIMENSION(:)    , INTENT(in) :: MASK_inner
    INTEGER,                    INTENT(in) :: k_outer, k_inner
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: ierr
    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, E_P, V1_P, V2_P, V3_P
    REAL(DP) :: D0_P, T0_P, Ye0_P, Ym0_P, E0_P, V10_P, V20_P, V30_P
    REAL(DP) :: Min_E, Min_E_0

    CHARACTER(len=100) :: old_state_filename

    IF ( ( ANY( MASK_inner ) .and. k_inner >= MaxIter_inner ) ) THEN
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
    !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
    !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP   D_old, Ye_old, Ym_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
    !$OMP   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
    !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
    !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC   D_old, Ye_old, Ym_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
    !$ACC   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#endif
      DO iN_X = 1, nX_G
        IF ( ( MASK_inner(iN_X) .and. k_inner >= MaxIter_inner ) ) THEN

          WRITE(old_state_filename,*) iN_X
          old_state_filename = "old_state_iN_X_"//TRIM(ADJUSTL(old_state_filename))//".dat"

          OPEN (UNIT=17, FILE=old_state_filename,ACTION='write')

          D_P   = D (iN_X) / Unit_D
          Ye_P  = Ye(iN_X) / Unit_Y
          Ym_P  = Ym(iN_X) / Unit_Y
          E_P   = E (iN_X) / Unit_E
          T_P   = T (iN_X) / Unit_T
          V1_P  = V_u_1(iN_X) / Unit_V
          V2_P  = V_u_2(iN_X) / Unit_V
          V3_P  = V_u_3(iN_X) / Unit_V

          D0_P   = D_old (iN_X) / Unit_D
          T0_P   = T_old (iN_X) / Unit_T
          Ye0_P  = Ye_old(iN_X) / Unit_Y
          Ym0_P  = Ym_old(iN_X) / Unit_Y
          E0_P   = E_old (iN_X) / Unit_E
          V10_P  = V_u_1_old(iN_X) / Unit_V
          V20_P  = V_u_2_old(iN_X) / Unit_V
          V30_P  = V_u_3_old(iN_X) / Unit_V

          CALL ComputeSpecificInternalEnergy_TABLE(D(iN_X),     Min_T, Ye(iN_X),    Ym(iN_X),     Min_E)
          CALL ComputeSpecificInternalEnergy_TABLE(D_old(iN_X), Min_T, Ye_old(iN_X),Ym_old(iN_X), Min_E_0)

          WRITE(17,'(8es23.15)') D0_P, T0_P, Ye0_P, Ym0_P, E0_P, V10_P, V20_P, V30_P
          WRITE(17,'(i5)')       nE_G

          WRITE(*,*)                      '[SolveNeutrinoMatterCoupling_FP_Nested_AA] Non-coverged inner loop'
          WRITE(*,'(a,2i5)')              '              iN_X        : ', iN_X
          WRITE(*,'(a,5x,2i23)')          '         k_outer, k_inner : ', k_outer, k_inner
          WRITE(*,'(a,5x,8es23.15)')      '    D, Ye, Ym, E, T, V_u       : ', D_P, Ye_P, Ym_P, E_P, T_P, V1_P, V2_P, V3_P
          WRITE(*,'(a,5x,8es23.15)')      '    D, Ye, Ym, E, T, V_u (old) : ', D0_P, Ye0_P, Ym0_P, E0_P, T0_P, V10_P, V20_P, V30_P
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
      D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, IgnoreIterationLimit_Option )

    USE mpi

    LOGICAL,  DIMENSION(:)    , INTENT(in) :: MASK_outer, MASK_inner
    INTEGER,                    INTENT(in) :: k_outer, k_inner
    INTEGER,  DIMENSION(:)    , INTENT(in) :: Error
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    LOGICAL, INTENT(in), OPTIONAL :: IgnoreIterationLimit_Option

    INTEGER  :: ierr
    INTEGER  :: iN_E, iN_X, iS
    LOGICAL  :: IgnoreIterationLimit
    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, E_P, V1_P, V2_P, V3_P
    REAL(DP) :: D0_P, T0_P, Ye0_P, Ym0_P, E0_P, V10_P, V20_P, V30_P
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
    !$OMP ( D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
    !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
    !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP   D_old, Ye_old, Ym_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
    !$OMP   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
    !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
    !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC   D_old, Ye_old, Ym_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
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

          D_P   = D (iN_X) / Unit_D
          Ye_P  = Ye(iN_X) / Unit_Y
          Ym_P  = Ym(iN_X) / Unit_Y
          E_P   = E (iN_X) / Unit_E
          T_P   = T (iN_X) / Unit_T
          V1_P  = V_u_1(iN_X) / Unit_V
          V2_P  = V_u_2(iN_X) / Unit_V
          V3_P  = V_u_3(iN_X) / Unit_V

          D0_P   = D_old (iN_X) / Unit_D
          T0_P   = T_old (iN_X) / Unit_T
          Ye0_P  = Ye_old(iN_X) / Unit_Y
          Ym0_P  = Ym_old(iN_X) / Unit_Y
          E0_P   = E_old (iN_X) / Unit_E
          V10_P  = V_u_1_old(iN_X) / Unit_V
          V20_P  = V_u_2_old(iN_X) / Unit_V
          V30_P  = V_u_3_old(iN_X) / Unit_V

          CALL ComputeSpecificInternalEnergy_TABLE(D(iN_X),     Min_T, Ye(iN_X),     Ym(iN_X),     Min_E)
          CALL ComputeSpecificInternalEnergy_TABLE(D_old(iN_X), Min_T, Ye_old(iN_X), Ym_old(iN_X), Min_E_0)

          WRITE(17,'(8es23.15)') D0_P, T0_P, Ye0_P, Ym0_P, E0_P, V10_P, V20_P, V30_P
          WRITE(17,'(i5)')       nE_G

          WRITE(*,*)                      '[SolveNeutrinoMatterCoupling_FP_Nested_AA] Error'
          WRITE(*,'(a,2i5)')              '              iN_X, Error : ', iN_X, Error(iN_X)
          WRITE(*,'(a,5x,2i23)')          '         k_outer, k_inner : ', k_outer, k_inner
          WRITE(*,'(a,5x,8es23.15)')      '    D, Ye, Ym, E, T, V_u       : ', D_P, Ye_P, Ym_P, E_P, T_P, V1_P, V2_P, V3_P
          WRITE(*,'(a,5x,8es23.15)')      '    D, Ye, Ym, E, T, V_u (old) : ', D0_P, Ye0_P, Ym0_P, E0_P, T0_P, V10_P, V20_P, V30_P
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


  SUBROUTINE PrintStatus_FP &
    ( MASK_outer, MASK_inner, k_outer, k_inner, Error, &
      Fm_outer, Fm_inner, Gm_outer, Gm_inner, &
      D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
      Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, iN_X )

    LOGICAL,  DIMENSION(:)    , INTENT(in) :: MASK_outer, MASK_inner
    INTEGER,                    INTENT(in) :: k_outer, k_inner
    INTEGER,  DIMENSION(:)    , INTENT(in) :: Error
    REAL(DP), DIMENSION(:,:)  , INTENT(in) :: Fm_outer, Fm_inner, Gm_outer, Gm_inner
    REAL(DP), DIMENSION(:)    , INTENT(in) :: D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dnu, Inu_u_1, Inu_u_2, Inu_u_3
    REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3
    REAL(DP), DIMENSION(:)    , INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER,                    INTENT(in) :: iN_X

    INTEGER  :: iN_E, iS, i
    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, E_P, V1_P, V2_P, V3_P
    REAL(DP) :: D0_P, T0_P, Ye0_P, Ym0_P, E0_P, V10_P, V20_P, V30_P
    REAL(DP) :: Min_E, Min_E_0

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( Fm_outer, Fm_inner, &
    !$OMP   Gm_outer, Gm_inner, &
    !$OMP   U_Ye, U_Ym, U_Ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$OMP   G_Ye, G_Ym, G_Ef, G_V_d_1, G_V_d_2, G_V_d_3, &
    !$OMP   D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
    !$OMP   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$OMP   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
    !$OMP   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$OMP   D_old, Ye_old, Ym_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
    !$OMP   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( Fm_outer, Fm_inner, &
    !$ACC   Gm_outer, Gm_inner, &
    !$ACC   U_Ye, U_Ym, U_Ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$ACC   D, Ye, Ym, E, T, V_u_1, V_u_2, V_u_3, &
    !$ACC   Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
    !$ACC   Nnu, Gnu_d_1, Gnu_d_2, Gnu_d_3, &
    !$ACC   Gm_dd_11, Gm_dd_22, Gm_dd_33, &
    !$ACC   D_old, Ye_old, Ym_old, E_old, T_old, V_u_1_old, V_u_2_old, V_u_3_old, &
    !$ACC   Dnu_old, Inu_u_1_old, Inu_u_2_old, Inu_u_3_old )
#endif

    D_P   = D (iN_X) / Unit_D
    Ye_P  = Ye(iN_X) / Unit_Y
    Ym_P  = Ym(iN_X) / Unit_Y
    E_P   = E (iN_X) / Unit_E
    T_P   = T (iN_X) / Unit_T
    V1_P  = V_u_1(iN_X) / Unit_V
    V2_P  = V_u_2(iN_X) / Unit_V
    V3_P  = V_u_3(iN_X) / Unit_V

    D0_P   = D_old (iN_X) / Unit_D
    T0_P   = T_old (iN_X) / Unit_T
    Ye0_P  = Ye_old(iN_X) / Unit_Y
    Ym0_P  = Ym_old(iN_X) / Unit_Y
    E0_P   = E_old (iN_X) / Unit_E
    V10_P  = V_u_1_old(iN_X) / Unit_V
    V20_P  = V_u_2_old(iN_X) / Unit_V
    V30_P  = V_u_3_old(iN_X) / Unit_V

    WRITE(*,*)                      '[SolveNeutrinoMatterCoupling_FP_Nested_AA] Convergence Status'
    WRITE(*,'(1x,4a8,5a23)')        'iQ', 'iN_X', 'k_outer', 'k_inner', 'Q_old', 'Q', 'U', 'F(iQ)', 'G(iQ)'
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iYe', iN_X, k_outer, k_inner, Ye0_P, Ye_P, U_Ye   (iN_X), Fm_outer(iYe,iN_X), G_Ye   (iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iYm', iN_X, k_outer, k_inner, Ym0_P, Ym_P, U_Ym   (iN_X), Fm_outer(iYm,iN_X), G_Ym   (iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iEf', iN_X, k_outer, k_inner,  E0_P,  E_P, U_Ef   (iN_X), Fm_outer(iEf,iN_X), G_Ef   (iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iV1', iN_X, k_outer, k_inner, V10_P, V1_P, U_V_d_1(iN_X), Fm_outer(iV1,iN_X), G_V_d_1(iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iV2', iN_X, k_outer, k_inner, V20_P, V2_P, U_V_d_2(iN_X), Fm_outer(iV2,iN_X), G_V_d_2(iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iV3', iN_X, k_outer, k_inner, V30_P, V3_P, U_V_d_3(iN_X), Fm_outer(iV3,iN_X), G_V_d_3(iN_X)

  END SUBROUTINE PrintStatus_FP


END MODULE TwoMoment_NeutrinoMatterSolverModule
