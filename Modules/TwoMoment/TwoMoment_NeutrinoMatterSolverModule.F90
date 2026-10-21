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
  USE ReferenceElementModuleE_Lagrange, ONLY: InterpMat_E
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
    ComputeElectronNeutrinoChemicalPotential_TABLE, &
    ComputeMuonNeutrinoChemicalPotential_TABLE, &
    ComputePressure_TABLE, &
    Min_D, Min_T, Min_Y, Max_T, Max_Y
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
    GetCollisionKernelScale, &
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
    ComputeNeutrinoLeptonKernels, ComputeNeutrinoLeptonRates, GetGeneralLeptonCoefficient, &
    Chi_GL, R_GL, Source_GL, nGL, iGL_NMS, iGL_LFE, iGL_LFC, iGL_IMD
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeEddingtonTensorComponents_dd, &
    ComputeConserved_TwoMoment
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
#if defined(TWOMOMENT_ORDER_V)
  USE TwoMoment_RetainedArithmeticModule
  USE TwoMoment_NewtonSolverModule_OrderV, ONLY: SolveNeutrinoMatterCoupling_Newton_OrderV_Batch
#endif
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: InitializeNeutrinoMatterSolverParameters
  PUBLIC :: FinalizeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling
  LOGICAL, PUBLIC :: UseNewtonSolve = .FALSE.


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
  REAL(DP), ALLOCATABLE, TARGET, PUBLIC :: NewtonRadiationIncrement(:,:,:,:)
  REAL(DP), ALLOCATABLE, TARGET, PUBLIC :: NewtonFluidIncrement(:,:)
  REAL(DP), ALLOCATABLE, TARGET :: AcceptedNewtonEmAbRate(:,:,:)
#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET(UseNewtonSolve,NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE(UseNewtonSolve,NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#endif
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
  LOGICAL :: UseMuonLeptonSource
  LOGICAL :: DebugMuonCoupling
  LOGICAL :: TraceNewtonSolve
  CHARACTER(512) :: NewtonTraceFile
  REAL(DP) :: NewtonFluidSourceTolerance
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
      Verbose_Option, SolverMethod_Option, DebugMuonCoupling_Option, &
      TraceNewtonSolve_Option, NewtonTraceFile_Option, &
      UseMuonLeptonSource_Option, NewtonFluidSourceTolerance_Option )


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
    LOGICAL , INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL , INTENT(in), OPTIONAL :: DebugMuonCoupling_Option
    LOGICAL, INTENT(in), OPTIONAL :: UseMuonLeptonSource_Option
    REAL(DP), INTENT(in), OPTIONAL :: NewtonFluidSourceTolerance_Option

    LOGICAL :: Verbose
    CHARACTER(*), INTENT(in), OPTIONAL :: SolverMethod_Option
    LOGICAL, INTENT(in), OPTIONAL :: TraceNewtonSolve_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: NewtonTraceFile_Option
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
      DnuMax = One
    END IF

    IF( PRESENT( FreezeOpacities_Option ) )THEN
      FreezeOpacities = FreezeOpacities_Option
    ELSE
      FreezeOpacities = .FALSE.
    END IF

    UseNewtonSolve = .FALSE.
    IF( PRESENT(SolverMethod_Option) )THEN
      SELECT CASE( TRIM(ADJUSTL(SolverMethod_Option)) )
        CASE( 'nested-aa' )
          UseNewtonSolve = .FALSE.
        CASE( 'balanced-newton' )
          UseNewtonSolve = .TRUE.
        CASE DEFAULT
          WRITE(ERROR_UNIT,'(A)') 'Expected SolverMethod=nested-aa or balanced-newton.'
          CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      END SELECT
    END IF
#if !defined(TWOMOMENT_ORDER_V) || !defined(MOMENT_CLOSURE_MINERBO)
    IF( UseNewtonSolve )THEN
      WRITE(ERROR_UNIT,'(A)') 'Balanced Newton currently requires Order-v and Minerbo closure.'
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    END IF
#endif
    UseMuonLeptonSource = .FALSE.
    IF( PRESENT(UseMuonLeptonSource_Option) ) UseMuonLeptonSource = UseMuonLeptonSource_Option
    IF( UseNewtonSolve )THEN
      IF( DnuMax /= One .OR. Include_LinCorr .OR. FreezeOpacities .OR. &
          ANY(wMatrRHS /= One) .OR. UseMuonLeptonSource )THEN
        WRITE(ERROR_UNIT,'(A)') 'Balanced Newton requires DnuMax=1, full fluid coupling, current opacities, '// &
          'no angular corrections, and no nested-source override.'
        CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      END IF
    END IF
    NewtonFluidSourceTolerance = Rtol_inner
    IF( PRESENT(NewtonFluidSourceTolerance_Option) )THEN
      IF( .NOT. IEEE_IS_FINITE(NewtonFluidSourceTolerance_Option) .OR. NewtonFluidSourceTolerance_Option<Zero )THEN
        WRITE(ERROR_UNIT,'(A)') 'Newton source tolerance must be finite and nonnegative (zero inherits Rtol_inner).'
        CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      END IF
      IF( NewtonFluidSourceTolerance_Option>Zero ) NewtonFluidSourceTolerance=NewtonFluidSourceTolerance_Option
    END IF
    DebugMuonCoupling = .FALSE.
    IF( PRESENT(DebugMuonCoupling_Option) ) DebugMuonCoupling = DebugMuonCoupling_Option
    TraceNewtonSolve = .FALSE.
    IF( PRESENT(TraceNewtonSolve_Option) ) TraceNewtonSolve = TraceNewtonSolve_Option
    NewtonTraceFile = 'newton_trace'
    IF( PRESENT(NewtonTraceFile_Option) ) NewtonTraceFile = TRIM(NewtonTraceFile_Option)

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
      WRITE(*,'(A4,A32,L1)') '', 'Balanced Newton: ', UseNewtonSolve
      WRITE(*,'(A4,A32,L1)') '', 'Direct muon-lepton source: ', UseMuonLeptonSource
      WRITE(*,'(A4,A32,L1)') '', 'Debug muon coupling: ', DebugMuonCoupling
      WRITE(*,'(A4,A32,L1)') '', 'Trace Newton: ', TraceNewtonSolve
      IF( UseNewtonSolve ) &
        WRITE(*,'(A4,A32,ES10.3E3)') '', 'Fluid source tolerance: ', NewtonFluidSourceTolerance
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

#if defined(TWOMOMENT_ORDER_V)
    IF( ALLOCATED(NewtonRadiationIncrement) )THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA MAP(release:NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA DELETE(NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#endif
      DEALLOCATE(NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
    END IF
#endif

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


  SUBROUTINE SolveNeutrinoMatterCoupling &
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
        CollisionRecoveryRadiationResidual_Relative_Final, &
        CollisionRecoveryMatterHeating_Rate_Final, &
        CollisionRecoveryStatus_Final, &
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
      CollisionRecoveryRadiationResidual_Relative_Final
    REAL(DP), DIMENSION(:),     INTENT(out)   :: &
      CollisionRecoveryMatterHeating_Rate_Final
    INTEGER,  DIMENSION(:),     INTENT(out)   :: CollisionRecoveryStatus_Final
    REAL(DP), DIMENSION(:,:,:), INTENT(out)   :: LeptonFlavorCheck_Rate_Final
#endif

    ! --- Local Variables ---

    INTEGER  :: ierr, myid
    INTEGER  :: k_outer, Mk_outer, nX_P_outer
    INTEGER  :: k_inner, Mk_inner, nX_P_inner
    INTEGER  :: iN_X
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: P
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: Error

    LOGICAL,  ALLOCATABLE, DIMENSION(:) :: ITERATE_outer, ITERATE_inner
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: PackIndex_outer, UnpackIndex_outer
    INTEGER,  ALLOCATABLE, DIMENSION(:) :: PackIndex_inner, UnpackIndex_inner

#if defined( TWOMOMENT_ORDER_V )
    LOGICAL, ALLOCATABLE :: CollisionCheckMask(:)
    REAL(DP), ALLOCATABLE :: InitialRadiationLimiter_ConservationCorrection(:,:)
#endif

    ! --- Least-squares scratch arrays ---

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: AMAT_outer, GVEC_outer, FVEC_outer
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: BVEC_outer, GVECm_outer, FVECm_outer
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: WORK_outer, TAU_outer, Alpha_outer

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: AMAT_inner, GVEC_inner, FVEC_inner
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: BVEC_inner, GVECm_inner, FVECm_inner
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)   :: WORK_inner, TAU_inner, Alpha_inner

#if defined(TWOMOMENT_ORDER_V)
    IF( UseNewtonSolve )THEN
      CALL SolveNeutrinoMatterCoupling_Newton_OrderV &
        ( dt,Dnu,Inu_u_1,Inu_u_2,Inu_u_3,V_u_1,V_u_2,V_u_3,D,T,Ye,Ym,E, &
          Gm_dd_11,Gm_dd_22,Gm_dd_33,Nnu,Gnu_d_1,Gnu_d_2,Gnu_d_3, &
          nIterations_Inner,nIterations_Outer,nIterations_Inner_Final,nIterations_Outer_Final, &
          RelativeResidual_Inner_First,MatterCoupling_Mask_Option, &
          MatterEquationResidual_Rate_Final,InitialRadiationLimiter_Rate_Final,RadiationResidual_Relative_Final, &
          CollisionRecoveryRadiationResidual_Relative_Final,CollisionRecoveryMatterHeating_Rate_Final, &
          CollisionRecoveryStatus_Final,LeptonFlavorCheck_Rate_Final,ComputeMatterEquationResidual )
      RETURN
    END IF
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO(UseNewtonSolve)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE(UseNewtonSolve)
#endif
#endif

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
    ALLOCATE( CollisionCheckMask(nX_G) )
    ALLOCATE( InitialRadiationLimiter_ConservationCorrection(nMatterEquations,nX_G) )
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
      CollisionRecoveryRadiationResidual_Relative_Final(iN_X) = Zero
      CollisionRecoveryMatterHeating_Rate_Final(iN_X) = Zero
      CollisionRecoveryStatus_Final(iN_X) = -2
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
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(alloc:CollisionCheckMask,InitialRadiationLimiter_ConservationCorrection)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA CREATE(CollisionCheckMask,InitialRadiationLimiter_ConservationCorrection)
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
    CollisionCheckMask = ITERATE_outer
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( CollisionCheckMask )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( CollisionCheckMask )
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




    ! --- Start Outer Loop ---

    CALL TimersStart( Timer_Collisions_OuterLoop )

    k_outer = 0
    DO WHILE( ANY( ITERATE_outer(:) ) .AND. k_outer < MaxIter_outer )
      k_outer = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      Error = 0

      CALL ComputeDnuNorm( ITERATE_outer, Dnu )

      IF ( k_outer > 1 ) THEN

        CALL CreatePackIndex &
               ( ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        IF ( .NOT. FreezeOpacities ) THEN


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


        END IF



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
            nIterations_Outer_Final(iN_X) = k_outer
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
      IF( UseMuonLeptonSource .AND. .NOT. .FALSE. &
            .AND. (Include_LFE .OR. Include_LFC .OR. Include_IMD) )THEN
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
    IF( ComputeMatterEquationResidual )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G
        CollisionCheckMask(iN_X) = nIterations_Outer_Final(iN_X) > 0
      END DO

      CALL CreatePackIndex &
             ( CollisionCheckMask, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

      CALL ComputeOpacities_Packed &
             ( D, T, Ye, Ym, SqrtGm, CollisionCheckMask, nX_P_inner, &
               PackIndex_inner, UnpackIndex_inner, &
               T0 = T_old, Ye0 = Ye_old, Ym0 = Ym_old, Opacity_Mask_Input = Opacity_Mask )

      CALL ComputeRates_Packed &
             ( D, T, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               CollisionCheckMask, nX_P_inner, PackIndex_inner, UnpackIndex_inner, &
               nX_P_inner, Opacity_Mask )

      CALL ComputeMatterRHS_OrderV &
             ( CollisionCheckMask, FVECm_outer, GVECm_outer, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, &
               dt, InitialRadiationLimiter_ConservationCorrection(iYm,:) )

      CALL ComputeNeutrinoRHS_OrderV &
             ( CollisionCheckMask, FVECm_inner, GVECm_inner, dt, &
               Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
               V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL ComputeLeptonFlavorCheck_OrderV &
             ( CollisionCheckMask, dt, &
               .FALSE., &
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
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP(release:CollisionCheckMask,InitialRadiationLimiter_ConservationCorrection)
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE(CollisionCheckMask,InitialRadiationLimiter_ConservationCorrection)
#endif
#endif

  END SUBROUTINE SolveNeutrinoMatterCoupling


#if defined(TWOMOMENT_ORDER_V)
  SUBROUTINE SolveNeutrinoMatterCoupling_Newton_OrderV &
    ( dt,Dnu,I1,I2,I3,V1,V2,V3,D,T,Ye,Ym,E,G1,G2,G3,Nnu,F1,F2,F3, &
      IterInner,IterOuter,FinalInner,FinalOuter,FirstResidual,CouplingMask, &
      MatterResidual,ClipRate,FinalResidual,ProjectionResidual,ProjectionHeating,ProjectionStatus,FlavorCheck, &
      ComputeCollisionChecks )

    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE FluidFieldsModule, ONLY: nCF,iCF_D,iCF_S1,iCF_S2,iCF_S3,iCF_E,iCF_Ne,iCF_Nm
    USE OpacityModule_TABLE, ONLY: C1,C2,C1_NuPair,C2_NuPair

    REAL(DP), INTENT(in) :: dt,G1(:),G2(:),G3(:)
    REAL(DP), INTENT(inout) :: Dnu(:,:,:),I1(:,:,:),I2(:,:,:),I3(:,:,:),V1(:),V2(:),V3(:)
    REAL(DP), INTENT(inout) :: D(:),T(:),Ye(:),Ym(:),E(:),Nnu(:,:,:),F1(:,:,:),F2(:,:,:),F3(:,:,:)
    INTEGER, INTENT(out) :: IterInner(:),IterOuter(:),FinalInner(:),FinalOuter(:),ProjectionStatus(:)
    REAL(DP), INTENT(out) :: FirstResidual(:),MatterResidual(:,:),ClipRate(:,:),FinalResidual(:)
    REAL(DP), INTENT(out) :: ProjectionResidual(:),ProjectionHeating(:),FlavorCheck(:,:,:)
    LOGICAL, INTENT(out), OPTIONAL :: CouplingMask(:)

    INTEGER :: NF,K,M,NA,NX,NB,Limit,Ierr,Rank,P,Q,R,S,J,B,IX,Start,Count,Point,TraceUnit,IOS
    INTEGER :: MaskCode
    INTEGER, ALLOCATABLE, TARGET :: Indices(:),Active(:),BatchIterations(:),BatchQueries(:),BatchFailure(:)
    INTEGER, ALLOCATABLE, TARGET :: QueryPack(:),QueryUnpack(:),LocalError(:)
    INTEGER, ALLOCATABLE :: Codes(:),AllIndices(:)
    LOGICAL, INTENT(in) :: ComputeCollisionChecks
    LOGICAL :: Enabled(9),ElectronActive,MuonActive,RecordChecks
    LOGICAL, ALLOCATABLE, TARGET :: StageMask(:),QueryMask(:),LocalElectron(:),LocalMuon(:),LocalOpacityMask(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: D0(:,:),I0(:,:,:),V0(:,:),Fluid0(:,:),W(:,:),U(:,:),Roots(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: BatchRadiationError(:),BatchConservationError(:)
    REAL(DP), ALLOCATABLE, TARGET :: BatchEmAb(:,:),BatchCheck(:,:,:)
    REAL(DP), ALLOCATABLE, TARGET :: WorkT(:),WorkYe(:),WorkYm(:),WorkE(:),WorkSqrtG(:)
    REAL(DP), ALLOCATABLE, TARGET :: OldE(:),OldYe(:),OldYm(:),OldV(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: LocalD(:),LocalT(:),LocalYe(:),LocalYm(:),LocalE(:),MuE(:),MuM(:)
    REAL(DP), ALLOCATABLE, TARGET :: DBScatter(:,:,:),DBPair(:,:,:)
    REAL(DP), TARGET :: KernelScale(5)
    TYPE(RetainedReal), ALLOCATABLE, TARGET :: BatchX(:,:),BatchIncrement(:,:,:)
    TYPE(RetainedReal) :: TF(7),DV(3),DI(3),NewD,NewI,DeltaEnergy,DeltaLepton
    REAL(DP) :: NClip,GClip(3),NBefore,GBefore(3),DBefore,DeltaDClip,DeltaNClip
    REAL(DP) :: V(3),Fwd,Back,Diag,Scale,ChiValue,PairSusceptibility,ProbeY

    RecordChecks = ComputeCollisionChecks .OR. TraceNewtonSolve
    NF = nSpecies/2
    K = NF+4
    M = nE_G*nSpecies
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,Rank,Ierr)
    IF( nSpecies < 2 .OR. MOD(nSpecies,2) /= 0 .OR. nSpecies > 6 .OR. dt <= Zero ) &
      CALL Fail('invalid species count or stage interval',0,0)
    TraceUnit = 0
    IF( TraceNewtonSolve )THEN
      OPEN(NEWUNIT=TraceUnit,FILE=TRIM(NewtonTraceFile)//'.newton_rank'//TRIM(IntegerString(Rank))//'.dat', &
        STATUS='UNKNOWN',POSITION='APPEND',IOSTAT=IOS)
      IF( IOS /= 0 ) CALL Fail('cannot open Newton solver trace',0,IOS)
    END IF

    IF( ALLOCATED(NewtonRadiationIncrement) )THEN
      IF( SIZE(NewtonRadiationIncrement,1) /= nE_G .OR. SIZE(NewtonRadiationIncrement,2) /= nSpecies &
          .OR. SIZE(NewtonRadiationIncrement,3) /= nX_G )THEN
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET EXIT DATA MAP(release:NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#elif defined(THORNADO_OACC)
        !$ACC EXIT DATA DELETE(NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#endif
        DEALLOCATE(NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
      END IF
    END IF
    IF( .NOT. ALLOCATED(NewtonRadiationIncrement) )THEN
      ALLOCATE(NewtonRadiationIncrement(nE_G,nSpecies,nX_G,4),NewtonFluidIncrement(nX_G,nCF))
      ALLOCATE(AcceptedNewtonEmAbRate(nE_G,nSpecies,nX_G))
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET ENTER DATA MAP(alloc:NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#elif defined(THORNADO_OACC)
      !$ACC ENTER DATA CREATE(NewtonRadiationIncrement,NewtonFluidIncrement,AcceptedNewtonEmAbRate)
#endif
    END IF
    ALLOCATE(StageMask(nX_G),QueryMask(nX_G),QueryPack(nX_G),QueryUnpack(nX_G))
    ALLOCATE(WorkT(nX_G),WorkYe(nX_G),WorkYm(nX_G),WorkE(nX_G),WorkSqrtG(nX_G))
    ALLOCATE(OldE(nX_G),OldYe(nX_G),OldYm(nX_G),OldV(3,nX_G),Codes(nX_G),AllIndices(nX_G))
    CALL SetOpacityMask(nX_G,nDOFX,D,T,Ye,Ym,Opacity_Mask)
    CALL SetMatterSourceMask(Opacity_Mask,MatterSource_Mask)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO(UseNewtonSolve)
    !$OMP TARGET UPDATE FROM(Opacity_Mask,MatterSource_Mask)
    !$OMP TARGET ENTER DATA MAP(alloc:StageMask,QueryMask,QueryPack,QueryUnpack,&
    !$OMP   WorkT,WorkYe,WorkYm,WorkE,WorkSqrtG,OldE,OldYe,OldYm,OldV)
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE(UseNewtonSolve)
    !$ACC UPDATE SELF(Opacity_Mask,MatterSource_Mask)
    !$ACC ENTER DATA CREATE(StageMask,QueryMask,QueryPack,QueryUnpack,&
    !$ACC   WorkT,WorkYe,WorkYm,WorkE,WorkSqrtG,OldE,OldYe,OldYm,OldV)
#endif
    DO P = 1, 5
      KernelScale(P) = GetCollisionKernelScale(P)
    END DO
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP(to:KernelScale)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE(Q,S,NClip,GClip,NBefore,GBefore,DBefore,DeltaDClip,DeltaNClip,J,V)
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA COPYIN(KernelScale)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE(Q,S,NClip,GClip,NBefore,GBefore,DBefore,DeltaDClip,DeltaNClip,J,V)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO PRIVATE(Q,S,NClip,GClip,NBefore,GBefore,DBefore,DeltaDClip,DeltaNClip,J,V)
#endif
    DO IX = 1, nX_G
      StageMask(IX) = MatterSource_Mask(iEf,IX)
      IF( PRESENT(CouplingMask) ) CouplingMask(IX) = StageMask(IX)
      WorkT(IX) = T(IX)
      WorkYe(IX) = Ye(IX)
      WorkYm(IX) = Ym(IX)
      WorkE(IX) = E(IX)
      WorkSqrtG(IX) = SQRT(G1(IX)*G2(IX)*G3(IX))
      OldE(IX) = E(IX)
      OldYe(IX) = Ye(IX)
      OldYm(IX) = Ym(IX)
      OldV(:,IX) = [SQRT(G1(IX))*V1(IX),SQRT(G2(IX))*V2(IX),SQRT(G3(IX))*V3(IX)]
      IterInner(IX) = 0
      IterOuter(IX) = 0
      FinalInner(IX) = 0
      FinalOuter(IX) = 0
      FirstResidual(IX) = Zero
      FinalResidual(IX) = Zero
      MatterResidual(:,IX) = Zero
      ProjectionResidual(IX) = Zero
      ProjectionHeating(IX) = Zero
      ProjectionStatus(IX) = -2
      ClipRate(:,IX) = Zero
      FlavorCheck(:,:,IX) = Zero
      NewtonRadiationIncrement(:,:,IX,:) = Zero
      NewtonFluidIncrement(IX,:) = Zero
      AcceptedNewtonEmAbRate(:,:,IX) = Zero
      IF( .NOT. StageMask(IX) ) CYCLE
      DO S = 1, nSpecies
      DO Q = 1, nE_G
        DBefore = Dnu(Q,S,IX)
        CALL ComputeConserved_TwoMoment(Dnu(Q,S,IX),I1(Q,S,IX),I2(Q,S,IX),I3(Q,S,IX), &
          NBefore,GBefore(1),GBefore(2),GBefore(3),V1(IX),V2(IX),V3(IX),G1(IX),G2(IX),G3(IX))
        CALL ProjectPrimitiveRadiationMoment_OrderV(Dnu(Q,S,IX),I1(Q,S,IX),I2(Q,S,IX),I3(Q,S,IX), &
          G1(IX),G2(IX),G3(IX),DnuMax)
        CALL ComputeConserved_TwoMoment(Dnu(Q,S,IX),I1(Q,S,IX),I2(Q,S,IX),I3(Q,S,IX), &
          NClip,GClip(1),GClip(2),GClip(3),V1(IX),V2(IX),V3(IX),G1(IX),G2(IX),G3(IX))
        NewtonRadiationIncrement(Q,S,IX,1) = NClip-Nnu(Q,S,IX)
        NewtonRadiationIncrement(Q,S,IX,2) = GClip(1)-F1(Q,S,IX)
        NewtonRadiationIncrement(Q,S,IX,3) = GClip(2)-F2(Q,S,IX)
        NewtonRadiationIncrement(Q,S,IX,4) = GClip(3)-F3(Q,S,IX)
        DeltaDClip = Dnu(Q,S,IX)-DBefore
        DeltaNClip = NClip-NBefore
        ClipRate(iQ_E,IX) = ClipRate(iQ_E,IX)+W3_S(Q)*(Two*DeltaNClip-DeltaDClip)/dt
        ClipRate(iQ_M1,IX) = ClipRate(iQ_M1,IX)+W3_S(Q)*(GClip(1)-GBefore(1)+G1(IX)*V1(IX)*DeltaDClip)/dt
        ClipRate(iQ_M2,IX) = ClipRate(iQ_M2,IX)+W3_S(Q)*(GClip(2)-GBefore(2)+G2(IX)*V2(IX)*DeltaDClip)/dt
        ClipRate(iQ_M3,IX) = ClipRate(iQ_M3,IX)+W3_S(Q)*(GClip(3)-GBefore(3)+G3(IX)*V3(IX)*DeltaDClip)/dt
        IF( S <= 2 ) ClipRate(iQ_Le,IX) = ClipRate(iQ_Le,IX) &
          +LeptonNumber(S)*W2_S(Q)*DeltaNClip/dt
        IF( S == 3 .OR. S == 4 ) ClipRate(iQ_Lm,IX) = ClipRate(iQ_Lm,IX) &
          +LeptonNumber(S)*W2_S(Q)*DeltaNClip/dt
        FlavorCheck((S+1)/2,iLFC_EntryProjection,IX) = FlavorCheck((S+1)/2,iLFC_EntryProjection,IX) &
          +LeptonNumber(S)*W2_S(Q)*DeltaNClip/dt
        FlavorCheck((S+1)/2,iLFC_InputConversion,IX) = FlavorCheck((S+1)/2,iLFC_InputConversion,IX) &
          +LeptonNumber(S)*W2_S(Q)*(NBefore-Nnu(Q,S,IX))/dt
        Nnu(Q,S,IX) = NClip
        F1(Q,S,IX) = GClip(1)
        F2(Q,S,IX) = GClip(2)
        F3(Q,S,IX) = GClip(3)
      END DO
      END DO
      ClipRate(iQ_H,IX) = ClipRate(iQ_E,IX)-V1(IX)*ClipRate(iQ_M1,IX) &
        -V2(IX)*ClipRate(iQ_M2,IX)-V3(IX)*ClipRate(iQ_M3,IX)
    END DO

    ! Cohorts have identical frozen ELEMENT masks; batch boundaries are
    ! multiples of nDOFX. No density/physics boundary is moved per Newton step.
    Limit = nX_G
    Codes = 0
    DO IX = 1, nX_G
      DO P = 1, nOp
        IF( Opacity_Mask(P,IX) ) Codes(IX) = IBSET(Codes(IX),P-1)
      END DO
    END DO
    DO WHILE( ANY(Codes /= 0) )
      Point = FINDLOC(Codes /= 0,.TRUE.,DIM=1)
      MaskCode = Codes(Point)
      ElectronActive = MatterSource_Mask(iYe,Point)
      MuonActive = MatterSource_Mask(iYm,Point)
      IF( .NOT. MatterSource_Mask(iEf,Point) )THEN
        WHERE(Codes == MaskCode) Codes = 0
        CYCLE
      END IF
      Count = 0
      DO IX = 1, nX_G
        IF( Codes(IX) /= MaskCode ) CYCLE
        Count = Count+1
        AllIndices(Count) = IX
      END DO
      NA = 4
      IF( ElectronActive ) NA = NA+1
      IF( MuonActive ) NA = NA+1
      NX = 4*M+NA+K
      ALLOCATE(Active(NA))
      J = 0
      IF( ElectronActive )THEN
        J = J+1
        Active(J) = 1
      END IF
      IF( MuonActive )THEN
        J = J+1
        Active(J) = 2
      END IF
      Active(J+1:NA) = [(P,P=NF+1,K)]
      Enabled = [Include_NNS .AND. Opacity_Mask(iOp_NNS,Point), &
                 Include_NES .AND. Opacity_Mask(iOp_NES,Point), &
                 Include_Pair .AND. Opacity_Mask(iOp_Pair,Point), &
                 Include_NuPair .AND. Opacity_Mask(iOp_NuPair,Point), &
                 Include_Brem .AND. Opacity_Mask(iOp_Brem,Point), &
                 Include_NMS .AND. Opacity_Mask(iOp_NMS,Point), &
                 Include_LFE .AND. Opacity_Mask(iOp_LFE,Point), &
                 Include_LFC .AND. Opacity_Mask(iOp_LFC,Point), &
                 Include_IMD .AND. Opacity_Mask(iOp_IMD,Point)]
      IF( TraceUnit /= 0 )THEN
        WRITE(TraceUnit,'(A,3I10,9L3)') 'newton_cohort ',MaskCode,Count,NA,Enabled
        WRITE(TraceUnit,'(A,*(I6))') 'newton_active ',Active
      END IF
      DO Start = 1, Count, Limit
        NB = MIN(Limit,Count-Start+1)
        ALLOCATE(Indices(NB),D0(M,NB),I0(M,3,NB),V0(3,NB),Fluid0(3,NB),W(M,NB),U(M,NB),Roots(3,NB))
        ALLOCATE(BatchX(NX,NB),BatchIncrement(M,4,NB),BatchIterations(NB),BatchQueries(NB),BatchFailure(NB))
        ALLOCATE(BatchRadiationError(NB),BatchConservationError(NB))
        ALLOCATE(BatchEmAb(M,NB))
        IF( RecordChecks )THEN
          ALLOCATE(BatchCheck(K,3,NB))
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET ENTER DATA MAP(alloc:BatchCheck)
#elif defined(THORNADO_OACC)
          !$ACC ENTER DATA CREATE(BatchCheck)
#endif
        END IF
        ALLOCATE(LocalD(NB),LocalT(NB),LocalYe(NB),LocalYm(NB),LocalE(NB),MuE(NB),MuM(NB),LocalError(NB))
        ALLOCATE(LocalElectron(NB),LocalMuon(NB),LocalOpacityMask(nOp,NB))
        ALLOCATE(DBScatter(nE_G,nE_G,NB),DBPair(nE_G,nE_G,NB))
        Indices = AllIndices(Start:Start+NB-1)
        IF( TraceUnit /= 0 )THEN
          WRITE(TraceUnit,'(A,3I10,ES24.15E3)') 'newton_batch ',Indices(1),Indices(NB),NB,dt
          FLUSH(TraceUnit)
        END IF
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET ENTER DATA MAP(to:Indices,Active) &
        !$OMP MAP(alloc:D0,I0,V0,Fluid0,W,U,Roots,BatchX,BatchIncrement,BatchIterations,BatchQueries, &
        !$OMP   BatchRadiationError,BatchConservationError,BatchEmAb,&
        !$OMP   LocalD,LocalT,LocalYe,LocalYm,LocalE,MuE,MuM,LocalError, &
        !$OMP   LocalElectron,LocalMuon,LocalOpacityMask,DBScatter,DBPair)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(IX,S,Q,J)
#elif defined(THORNADO_OACC)
        !$ACC ENTER DATA COPYIN(Indices,Active)
        !$ACC ENTER DATA CREATE(D0,I0,V0,Fluid0,W,U,Roots,BatchX,BatchIncrement,BatchIterations,BatchQueries, &
        !$ACC   BatchRadiationError,BatchConservationError,BatchEmAb,&
        !$ACC   LocalD,LocalT,LocalYe,LocalYm,LocalE,MuE,MuM,LocalError, &
        !$ACC   LocalElectron,LocalMuon,LocalOpacityMask,DBScatter,DBPair)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(IX,S,Q,J)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(IX,S,Q,J)
#endif
        DO B = 1, NB
          IX = Indices(B)
          Roots(:,B) = SQRT([G1(IX),G2(IX),G3(IX)])
          V0(:,B) = OldV(:,IX)
          Fluid0(:,B) = [OldE(IX),OldYe(IX),OldYm(IX)]
          DO S = 1, nSpecies
          DO Q = 1, nE_G
            J = (S-1)*nE_G+Q
            D0(J,B) = Dnu(Q,S,IX)
            I0(J,:,B) = Roots(:,B)*[I1(Q,S,IX),I2(Q,S,IX),I3(Q,S,IX)]
            W(J,B) = AtomicMassUnit*W2_S(Q)/D(IX)
            U(J,B) = W3_S(Q)/D(IX)
          END DO
          END DO
        END DO
        CALL SolveNeutrinoMatterCoupling_Newton_OrderV_Batch(nE_G,nSpecies,nDOFE,Indices,Active,Enabled,E_N,WeightsE,InterpMat_E, &
          W,U,D0,I0,V0,Fluid0,dt,DnuMax,Rtol_inner,MaxIter_inner,QueryCoefficients,BatchX,BatchIncrement, &
          BatchIterations,BatchRadiationError,BatchConservationError,BatchQueries,BatchFailure,TraceUnit,BatchEmAb,BatchCheck, &
          FluidSourceTolerance_Option=NewtonFluidSourceTolerance)
        IF( ANY(BatchFailure /= 0) )THEN
          B = FINDLOC(BatchFailure /= 0,.TRUE.,DIM=1)
          CALL Fail('nonlinear batch failure',Indices(B),BatchFailure(B))
        END IF
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(IX,J,Q,S,TF,DV,NewD,NewI)
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(IX,J,Q,S,TF,DV,NewD,NewI)
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO PRIVATE(IX,J,Q,S,TF,DV,NewD,NewI)
#endif
        DO B = 1, NB
          IX = Indices(B)
          TF = 0.0_DP
          DO J = 1, NA
            TF(Active(J)) = BatchX(4*M+J,B)
          END DO
          DV = -TF(NF+2:K)
          Ye(IX) = RetainedValue(OldYe(IX)-TF(1))
          Ym(IX) = OldYm(IX)
          IF( NF >= 2 ) Ym(IX) = RetainedValue(OldYm(IX)-TF(2))
          E(IX) = RetainedValue(OldE(IX)-TF(NF+1) &
            -OldV(1,IX)*DV(1)-OldV(2,IX)*DV(2)-OldV(3,IX)*DV(3) &
            -0.5_DP*(DV(1)*DV(1)+DV(2)*DV(2)+DV(3)*DV(3)))
          T(IX) = WorkT(IX)
          V1(IX) = RetainedValue(OldV(1,IX)+DV(1))/Roots(1,B)
          V2(IX) = RetainedValue(OldV(2,IX)+DV(2))/Roots(2,B)
          V3(IX) = RetainedValue(OldV(3,IX)+DV(3))/Roots(3,B)
          NewtonFluidIncrement(IX,iCF_Ne) = -D(IX)/AtomicMassUnit*RetainedValue(TF(1))
          IF( NF >= 2 ) NewtonFluidIncrement(IX,iCF_Nm) = -D(IX)/AtomicMassUnit*RetainedValue(TF(2))
          NewtonFluidIncrement(IX,iCF_E) = -D(IX)*RetainedValue(TF(NF+1))
          NewtonFluidIncrement(IX,iCF_S1) = -D(IX)*Roots(1,B)*RetainedValue(TF(NF+2))
          NewtonFluidIncrement(IX,iCF_S2) = -D(IX)*Roots(2,B)*RetainedValue(TF(NF+3))
          NewtonFluidIncrement(IX,iCF_S3) = -D(IX)*Roots(3,B)*RetainedValue(TF(NF+4))
          IterInner(IX) = BatchIterations(B)
          IterOuter(IX) = 1
          FinalInner(IX) = BatchIterations(B)
          FinalOuter(IX) = 1
          FinalResidual(IX) = BatchRadiationError(B)
          IF( RecordChecks )THEN
            MatterResidual(iQ_Le,IX) = -D(IX)/AtomicMassUnit*BatchCheck(1,2,B)/dt
            IF( NF >= 2 ) MatterResidual(iQ_Lm,IX) = -D(IX)/AtomicMassUnit*BatchCheck(2,2,B)/dt
            MatterResidual(iQ_E,IX) = -D(IX)*BatchCheck(NF+1,2,B)/dt
            MatterResidual(iQ_M1,IX) = -D(IX)*Roots(1,B)*BatchCheck(NF+2,2,B)/dt
            MatterResidual(iQ_M2,IX) = -D(IX)*Roots(2,B)*BatchCheck(NF+3,2,B)/dt
            MatterResidual(iQ_M3,IX) = -D(IX)*Roots(3,B)*BatchCheck(NF+4,2,B)/dt
            MatterResidual(iQ_H,IX) = MatterResidual(iQ_E,IX)-V1(IX)*MatterResidual(iQ_M1,IX) &
              -V2(IX)*MatterResidual(iQ_M2,IX)-V3(IX)*MatterResidual(iQ_M3,IX)
            DO J = 1,NF
              FlavorCheck(J,iLFC_CollisionSource,IX) = D(IX)/AtomicMassUnit*BatchCheck(J,3,B)
              ! The process-source map is transfer-g, not the accepted
              ! radiation change transfer+H; their difference is H+g.
              FlavorCheck(J,iLFC_RawMap,IX) = D(IX)/AtomicMassUnit &
                *RetainedValue(TF(J)-BatchCheck(J,2,B))/dt
              FlavorCheck(J,iLFC_EquationResidual,IX) = D(IX)/AtomicMassUnit &
                *(BatchCheck(J,1,B)+BatchCheck(J,2,B))/dt
              FlavorCheck(J,iLFC_FixedPointResidual,IX) = -FlavorCheck(J,iLFC_EquationResidual,IX)
            END DO
          END IF
          DO S = 1, nSpecies
          DO Q = 1, nE_G
            J = (S-1)*nE_G+Q
            AcceptedNewtonEmAbRate(Q,S,IX) = BatchEmAb(J,B)
            Dnu(Q,S,IX) = RetainedValue(D0(J,B)+BatchX(J,B))
            I1(Q,S,IX) = RetainedValue(I0(J,1,B)+BatchX(M+J,B))/Roots(1,B)
            I2(Q,S,IX) = RetainedValue(I0(J,2,B)+BatchX(2*M+J,B))/Roots(2,B)
            I3(Q,S,IX) = RetainedValue(I0(J,3,B)+BatchX(3*M+J,B))/Roots(3,B)
            NewtonRadiationIncrement(Q,S,IX,1) = RetainedValue( &
              Retain(NewtonRadiationIncrement(Q,S,IX,1))+BatchIncrement(J,1,B))
            NewtonRadiationIncrement(Q,S,IX,2) = RetainedValue( &
              Retain(NewtonRadiationIncrement(Q,S,IX,2))+Roots(1,B)*BatchIncrement(J,2,B))
            NewtonRadiationIncrement(Q,S,IX,3) = RetainedValue( &
              Retain(NewtonRadiationIncrement(Q,S,IX,3))+Roots(2,B)*BatchIncrement(J,3,B))
            NewtonRadiationIncrement(Q,S,IX,4) = RetainedValue( &
              Retain(NewtonRadiationIncrement(Q,S,IX,4))+Roots(3,B)*BatchIncrement(J,4,B))
          END DO
          END DO
        END DO
#if defined(THORNADO_OMP_OL)
        !$OMP TARGET EXIT DATA MAP(release:Indices,Active,D0,I0,V0,Fluid0,W,U,Roots,BatchX,BatchIncrement, &
        !$OMP   BatchIterations,BatchQueries,BatchRadiationError,BatchConservationError,BatchEmAb,&
        !$OMP   LocalD,LocalT,LocalYe,LocalYm, &
        !$OMP   LocalE,MuE,MuM,LocalError,LocalElectron,LocalMuon,LocalOpacityMask,&
        !$OMP   DBScatter,DBPair)
#elif defined(THORNADO_OACC)
        !$ACC EXIT DATA DELETE(Indices,Active,D0,I0,V0,Fluid0,W,U,Roots,BatchX,BatchIncrement, &
        !$ACC   BatchIterations,BatchQueries,BatchRadiationError,BatchConservationError,BatchEmAb,&
        !$ACC   LocalD,LocalT,LocalYe,LocalYm, &
        !$ACC   LocalE,MuE,MuM,LocalError,LocalElectron,LocalMuon,LocalOpacityMask,&
        !$ACC   DBScatter,DBPair)
#endif
        DEALLOCATE(Indices,D0,I0,V0,Fluid0,W,U,Roots,BatchX,BatchIncrement,BatchIterations,BatchQueries,BatchFailure)
        DEALLOCATE(BatchRadiationError,BatchConservationError,BatchEmAb)
        IF( ALLOCATED(BatchCheck) )THEN
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET EXIT DATA MAP(release:BatchCheck)
#elif defined(THORNADO_OACC)
          !$ACC EXIT DATA DELETE(BatchCheck)
#endif
          DEALLOCATE(BatchCheck)
        END IF
        DEALLOCATE(LocalD,LocalT,LocalYe,LocalYm,LocalE,MuE,MuM,LocalError)
        DEALLOCATE(LocalElectron,LocalMuon,LocalOpacityMask,DBScatter,DBPair)
      END DO
      WHERE(Codes == MaskCode) Codes = 0
      DEALLOCATE(Active)
    END DO

    ! Every converged point's final kernel query was its accepted state:
    ! subsequent probes/trials exclude it. Reuse those kernels for the
    ! existing AP-transport and process-Q interfaces, without a new opacity
    ! query or a primitive-to-conserved subtraction for the actual increment.
    CALL ComputeRates_Packed(D,T,Dnu,I1,I2,I3,nX_P0=nX_G,Opacity_Mask_Input=Opacity_Mask)
    IF( TraceUnit /= 0 ) CLOSE(TraceUnit)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP(release:StageMask,QueryMask,QueryPack,QueryUnpack,WorkT,WorkYe,WorkYm,WorkE, &
    !$OMP   WorkSqrtG,OldE,OldYe,OldYm,OldV,KernelScale)
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA DELETE(StageMask,QueryMask,QueryPack,QueryUnpack,WorkT,WorkYe,WorkYm,WorkE, &
    !$ACC   WorkSqrtG,OldE,OldYe,OldYm,OldV,KernelScale)
#endif

  CONTAINS

    SUBROUTINE Fail(Message,Point,Code)
      CHARACTER(*), INTENT(in) :: Message
      INTEGER, INTENT(in) :: Point,Code
      WRITE(ERROR_UNIT,'(A,A,3I10)') 'Balanced Newton: ',TRIM(Message),Rank,Point,Code
      FLUSH(ERROR_UNIT)
      CALL MPI_ABORT(MPI_COMM_WORLD,-1,Ierr)
    END SUBROUTINE Fail

    SUBROUTINE QueryCoefficients &
      ( Mode,Points,Requested,Transfer,Process,Species,QA,QB,KT,Mu,Susceptibility,Chi,Sigma,Gain,Loss,Error )
      INTEGER, INTENT(in) :: Mode,Points(:),Process(:),Species(:),QA(:),QB(:)
      LOGICAL, INTENT(in) :: Requested(:)
      TYPE(RetainedReal), INTENT(in) :: Transfer(:,:)
      REAL(DP), INTENT(out) :: KT(:),Mu(:,:),Susceptibility(:),Chi(:,:),Sigma(:,:),Gain(:,:),Loss(:,:)
      INTEGER, INTENT(out) :: Error(:)
      INTEGER :: B,IX,Q,S,J,P,L,R,LO,HI,C,Count,Ierr
      REAL(DP) :: G,B0,V2,Fwd,Back,Diag,ProbeY,PairSusceptibility
      TYPE(RetainedReal) :: DV(3),Energy

      CALL TimersStart(Timer_Collisions_MatterFromRadiation)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(IX,DV,Energy)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(IX,DV,Energy)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(IX,DV,Energy)
#endif
      DO B = 1, SIZE(Points)
        IX = Points(B)
        LocalD(B) = D(IX)
        LocalT(B) = WorkT(IX)
        LocalYe(B) = OldYe(IX)
        LocalYm(B) = OldYm(IX)
        LocalE(B) = OldE(IX)
        LocalError(B) = 0
        Error(B) = 0
        LocalElectron(B) = Opacity_Mask(iOp_EmAb_Electron,IX) .OR. Opacity_Mask(iOp_ECTable,IX)
        LocalMuon(B) = Opacity_Mask(iOp_EmAb_Muon,IX)
        LocalOpacityMask(:,B) = Opacity_Mask(:,IX)
        IF( .NOT. Requested(B) ) CYCLE
        DV = -Transfer(NF+2:K,B)
        Energy = OldE(IX)-Transfer(NF+1,B)-OldV(1,IX)*DV(1)-OldV(2,IX)*DV(2)-OldV(3,IX)*DV(3) &
          -0.5_DP*(DV(1)*DV(1)+DV(2)*DV(2)+DV(3)*DV(3))
        LocalE(B) = RetainedValue(Energy)
        LocalYe(B) = RetainedValue(OldYe(IX)-Transfer(1,B))
        IF( NF >= 2 ) LocalYm(B) = RetainedValue(OldYm(IX)-Transfer(2,B))
        IF( .NOT. IEEE_IS_FINITE(LocalE(B)) .OR. LocalYe(B) <= Zero .OR. LocalYe(B) >= One &
            .OR. ABS(LocalYm(B)) >= One )THEN
          LocalError(B) = 1
          LocalYe(B) = OldYe(IX)
          LocalYm(B) = OldYm(IX)
          LocalE(B) = OldE(IX)
        END IF
      END DO
      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
        (LocalD,LocalE,LocalYe,LocalYm,LocalT,Error_Option=Error,DisplayErrors_Option=.FALSE.)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO B = 1, SIZE(Points)
        Error(B) = MAX(Error(B),LocalError(B))
      END DO
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM(Error)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE SELF(Error)
#endif
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(IX)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(IX)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(IX)
#endif
      DO B = 1, SIZE(Points)
        IF( Error(B) == 0 ) CYCLE
        IX = Points(B)
        LocalT(B) = T(IX)
        LocalYe(B) = OldYe(IX)
        LocalYm(B) = OldYm(IX)
      END DO
      CALL ComputeElectronNeutrinoChemicalPotential_TABLE &
        (LocalD,LocalT,LocalYe,LocalYm,MuE,ChargedCurrentMask_Option=LocalElectron)
      CALL ComputeMuonNeutrinoChemicalPotential_TABLE &
        (LocalD,LocalT,LocalYe,LocalYm,MuM,ChargedCurrentMask_Option=LocalMuon)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ProbeY,PairSusceptibility,Ierr,S)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(ProbeY,PairSusceptibility,Ierr,S)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(ProbeY,PairSusceptibility,Ierr,S)
#endif
      DO B = 1, SIZE(Points)
        IF( .NOT. Requested(B) .OR. Error(B) /= 0 ) CYCLE
        KT(B) = BoltzmannConstant*LocalT(B)
        Susceptibility(B) = Zero
        IF( LocalMuon(B) )THEN
          CALL ComputeMuonFractionFromDegeneracy_TABLE(LocalD(B),LocalT(B),LocalYe(B),Zero, &
            ProbeY,PairSusceptibility,Ierr)
          IF( Ierr == 0 ) Susceptibility(B) = PairSusceptibility
        END IF
        DO S = 1, nSpecies
          Mu(S,B) = Zero
          IF( S <= 2 ) Mu(S,B) = LeptonNumber(S)*MuE(B)
          IF( S == 3 .OR. S == 4 ) Mu(S,B) = LeptonNumber(S)*MuM(B)
        END DO
      END DO
      CALL TimersStop(Timer_Collisions_MatterFromRadiation)
      IF( Mode == 2 ) RETURN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO IX = 1, nX_G
        QueryMask(IX) = .FALSE.
      END DO
      Count = 0
      DO B = 1, SIZE(Points)
        IF( .NOT. Requested(B) .OR. Error(B) /= 0 ) CYCLE
        Count = Count+1
        QueryUnpack(Count) = Points(B)
        QueryPack(Points(B)) = Count
      END DO
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(QueryUnpack,QueryPack)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(IX)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(QueryUnpack,QueryPack)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(IX)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(IX)
#endif
      DO B = 1, SIZE(Points)
        IF( .NOT. Requested(B) .OR. Error(B) /= 0 ) CYCLE
        IX = Points(B)
        QueryMask(IX) = .TRUE.
        WorkT(IX) = LocalT(B)
        WorkYe(IX) = LocalYe(B)
        WorkYm(IX) = LocalYm(B)
        WorkE(IX) = LocalE(B)
      END DO
      IF( Count == 0 ) RETURN
      CALL TimersStart(Timer_Collisions_ComputeOpacity)
      CALL ComputeOpacities_Packed(D,WorkT,WorkYe,WorkYm,WorkSqrtG,QueryMask,Count,QueryPack,QueryUnpack,nX_G, &
        Opacity_Mask_Input=Opacity_Mask)
      CALL TimersStop(Timer_Collisions_ComputeOpacity)
      CALL ComputeDetailedBalanceFactors(1,nE_G,1,SIZE(Points),E_N,LocalT,.TRUE.,.TRUE.,DBScatter,DBPair)
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(IX,S,Q,R,J,P,L,C,Fwd,Back,Diag)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(IX,S,Q,R,J,P,L,C,Fwd,Back,Diag)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO PRIVATE(IX,S,Q,R,J,P,L,C,Fwd,Back,Diag)
#endif
      DO B = 1, SIZE(Points)
        IF( .NOT. Requested(B) .OR. Error(B) /= 0 ) CYCLE
        IX = Points(B)
        DO S = 1, nSpecies
        DO R = 1, nE_G
          J = (S-1)*nE_G+R
          C = nChirals-(LeptonNumber(S)+1)/nChirals
          Chi(J,B) = Chi_EmAb(R,S,IX)
          Sigma(J,B) = Sigma_Iso(R,C,IX)
          ! Diagonal self-scattering cancels identically in C0, but damps I.
          DO P = 1, 6
            IF( .NOT. Enabled(P) .OR. (P >= 3 .AND. P <= 5) ) CYCLE
            CALL ComputeCollisionKernel_OrderV_Point(P,S,R,R,IX,KernelScale,Fwd,Back)
            Sigma(J,B) = Sigma(J,B)+Back*W2_N(R)
          END DO
        END DO
        END DO
        DO L = 1, SIZE(Process)
          P = Process(L)
          S = Species(L)
          R = QA(L)
          Q = QB(L)
          CALL ComputeCollisionKernel_OrderV_Point(P,S,R,Q,IX,KernelScale,Fwd,Back)
          IF( P <= 2 .OR. P == 6 )THEN
            Fwd = Back*DBScatter(Q,R,B)
          ELSE IF( P <= 5 )THEN
            Fwd = Back*DBPair(Q,R,B)
          END IF
          Gain(L,B) = Fwd*W2_N(Q)
          Loss(L,B) = Back*W2_N(Q)
        END DO
      END DO
    END SUBROUTINE QueryCoefficients

    FUNCTION IntegerString(I) RESULT(S)
      INTEGER, INTENT(in) :: I
      CHARACTER(24) :: S
      WRITE(S,'(I0)') I
    END FUNCTION IntegerString

  END SUBROUTINE SolveNeutrinoMatterCoupling_Newton_OrderV

  SUBROUTINE ComputeCollisionKernel_OrderV_Point(P,S,R,Q,IX,KernelScale,Fwd,Back)
    USE OpacityModule_TABLE, ONLY: C1,C2,C1_NuPair,C2_NuPair
#if defined(THORNADO_OMP_OL)
      !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
      !$ACC ROUTINE SEQ
#endif
      INTEGER, INTENT(in) :: P,S,R,Q,IX
      REAL(DP), INTENT(in) :: KernelScale(5)
      REAL(DP), INTENT(out) :: Fwd,Back
      INTEGER :: LO,HI
      LO = MIN(R,Q)
      HI = MAX(R,Q)
      Back = Zero
      Fwd = Zero
      SELECT CASE(P)
        CASE(1)
          Back = Phi_NNS_0(LO,HI,IX)
          IF( MOD(S,2) == 0 ) Back = Phi_NbNS_0(LO,HI,IX)
        CASE(2)
          Back = C1(S)*H_I_0(LO,HI,IX)+C2(S)*H_II_0(LO,HI,IX)
        CASE(3)
          IF( Q <= R )THEN
            Back = C1(S)*J_I_0(LO,HI,IX)+C2(S)*J_II_0(LO,HI,IX)
          ELSE
            Back = C1(S)*J_II_0(LO,HI,IX)+C2(S)*J_I_0(LO,HI,IX)
          END IF
        CASE(4)
          IF( Q <= R )THEN
            Back = C1_NuPair(S)*Nu_J_I_0(LO,HI,IX)+C2_NuPair(S)*Nu_J_II_0(LO,HI,IX)
          ELSE
            Back = C1_NuPair(S)*Nu_J_II_0(LO,HI,IX)+C2_NuPair(S)*Nu_J_I_0(LO,HI,IX)
          END IF
        CASE(5)
          Back = S_Sigma(LO,HI,IX)
        CASE(6:9)
          CALL GetGeneralLeptonCoefficient(P-5,S,R,Q,IX,Fwd,Back)
      END SELECT
      IF( P <= 5 )THEN
        Back = Back*KernelScale(P)
        Fwd = Back
      END IF
  END SUBROUTINE ComputeCollisionKernel_OrderV_Point
#endif


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






  SUBROUTINE ComputeCollisionMoments_OrderV_Point &
    ( iN_E, iS, iN_X, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, RN, RG1, RG2, RG3, LinearRN_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER, INTENT(in) :: iN_E, iS, iN_X
    REAL(DP), INTENT(in) :: Dnu(:,:,:), Inu_u_1(:,:,:), Inu_u_2(:,:,:), Inu_u_3(:,:,:)
    REAL(DP), INTENT(in) :: Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP), INTENT(out) :: RN(nC), RG1(nC), RG2(nC), RG3(nC)
    REAL(DP), INTENT(out), OPTIONAL :: LinearRN_Option(nC)
    REAL(DP) :: LinearRN(nC)
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: I_d_1, I_d_2, I_d_3
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    INTEGER :: iG

    LinearRN = Zero
      I_d_1 = Gm_dd_11(iN_X) * Inu_u_1(iN_E,iS,iN_X)
      I_d_2 = Gm_dd_22(iN_X) * Inu_u_2(iN_E,iS,iN_X)
      I_d_3 = Gm_dd_33(iN_X) * Inu_u_3(iN_E,iS,iN_X)

      IF( Include_LinCorr ) CALL ComputeEddingtonTensorComponents_dd &
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
          LinearRN(iC_NNS) = -( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
          RN(iC_NNS) = RN(iC_NNS) + LinearRN(iC_NNS)
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
          LinearRN(iC_NES) = -( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
          RN(iC_NES) = RN(iC_NES) + LinearRN(iC_NES)
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
          LinearRN(iC_Pair) = -( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
          RN(iC_Pair) = RN(iC_Pair) + LinearRN(iC_Pair)
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
          LinearRN(iC_NuPair) = -( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
          RN(iC_NuPair) = RN(iC_NuPair) + LinearRN(iC_NuPair)
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
          LinearRN(iC_Brem) = -( L_u_1 * I_d_1 + L_u_2 * I_d_2 + L_u_3 * I_d_3 )
          RN(iC_Brem) = RN(iC_Brem) + LinearRN(iC_Brem)
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


    IF( PRESENT(LinearRN_Option) ) LinearRN_Option = LinearRN

  END SUBROUTINE ComputeCollisionMoments_OrderV_Point


  SUBROUTINE ComputeMatterRHS_Collisions_OrderV &
    ( Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
      V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33, MASK, Q_Matter, &
      NNS_RateDecomposition_Option, QH_ByProcessSpecies_Option, EmAb_RateDecomposition_Option, &
      EmAb_Spectrum_Option, RadiationCollisionRHS_Option, UseAcceptedEmAbRate_Option )

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
    LOGICAL, INTENT(in), OPTIONAL :: UseAcceptedEmAbRate_Option

    REAL(DP), ALLOCATABLE :: Q_SP(:,:,:,:), QLm_SP(:,:,:)
    REAL(DP), ALLOCATABLE :: EmAb_Spectrum(:,:,:,:)
    REAL(DP), ALLOCATABLE :: RadiationCollisionRHS(:,:,:,:)
    REAL(DP), ALLOCATABLE :: QH_ByProcessSpecies(:,:,:)

    INTEGER  :: ierr, myid
    LOGICAL  :: ComputeRadiationCollisionRHS, InvalidShape, UseAcceptedEmAbRate
    INTEGER  :: iG, iN_E, iN_X, iP, iQ, iS
    REAL(DP) :: k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33
    REAL(DP) :: I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3
    REAL(DP) :: RN(nC)
    REAL(DP) :: RG1(nC)
    REAL(DP) :: RG2(nC)
    REAL(DP) :: RG3(nC)
    REAL(DP) :: QL, QLm, QE, QM1, QM2, QM3
    REAL(DP) :: Q0(nC), Q0Correction(nC), Q0Term, LinearRN(nC)
    REAL(DP) :: RN_Total, RG1_Total, RG2_Total, RG3_Total
    REAL(DP) :: V2

    ComputeRadiationCollisionRHS = PRESENT( RadiationCollisionRHS_Option )
    ! Only the immediate post-solve caller may use the retained accepted
    ! source. Arbitrary-state diagnostics initialize their own collision
    ! storage and must evaluate the source of the state they supply.
    UseAcceptedEmAbRate = .FALSE.
    IF( PRESENT(UseAcceptedEmAbRate_Option) ) &
      UseAcceptedEmAbRate = UseAcceptedEmAbRate_Option

    IF( UseAcceptedEmAbRate )THEN
      IF( .NOT. ALLOCATED(AcceptedNewtonEmAbRate) )THEN
        WRITE(ERROR_UNIT,'(A)') 'ERROR: Accepted EmAb source requested without an accepted solve.'
        CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
      END IF
      IF( ANY(SHAPE(AcceptedNewtonEmAbRate) /= [nE_G,nSpecies,nX_G]) )THEN
        WRITE(ERROR_UNIT,'(A)') 'ERROR: Accepted EmAb source has inconsistent dimensions.'
        CALL MPI_ABORT( MPI_COMM_WORLD, -1, ierr )
      END IF
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE TO(AcceptedNewtonEmAbRate)
#elif defined(THORNADO_OACC)
      !$ACC UPDATE DEVICE(AcceptedNewtonEmAbRate)
#endif
    END IF

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
    DO iQ   = 1, SIZE(Q_SP,1)
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
    !$OMP FIRSTPRIVATE( ComputeRadiationCollisionRHS, UseAcceptedEmAbRate ) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$OMP          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$OMP          RG3_Total, QL, QLm, QE, QM1, QM2, QM3, Q0, Q0Correction, Q0Term, LinearRN, iN_E, iP, iG )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC FIRSTPRIVATE( ComputeRadiationCollisionRHS, UseAcceptedEmAbRate ) &
    !$ACC PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$ACC          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$ACC          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$ACC          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$ACC          RG3_Total, QL, QLm, QE, QM1, QM2, QM3, Q0, Q0Correction, Q0Term, LinearRN, iN_E, iP, iG ) &
    !$ACC PRESENT( Q_SP, QLm_SP, RadiationCollisionRHS )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP FIRSTPRIVATE( ComputeRadiationCollisionRHS, UseAcceptedEmAbRate ) &
    !$OMP PRIVATE( k_dd_11, k_dd_12, k_dd_13, k_dd_22, k_dd_23, k_dd_33, &
    !$OMP          I_d_1, I_d_2, I_d_3, V_d_1, V_d_2, V_d_3, &
    !$OMP          L_u_1, L_u_2, L_u_3, L_d_1, L_d_2, L_d_3, &
    !$OMP          RN, RG1, RG2, RG3, RN_Total, RG1_Total, RG2_Total, &
    !$OMP          RG3_Total, QL, QLm, QE, QM1, QM2, QM3, Q0, Q0Correction, Q0Term, LinearRN, iN_E, iP, iG )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies

      IF ( MASK(iN_X) ) THEN

        V_d_1 = Gm_dd_11(iN_X) * V_u_1(iN_X)
        V_d_2 = Gm_dd_22(iN_X) * V_u_2(iN_X)
        V_d_3 = Gm_dd_33(iN_X) * V_u_3(iN_X)
        Q0=Zero
        Q0Correction=Zero

        DO iN_E = 1, nE_G

          CALL ComputeCollisionMoments_OrderV_Point &
                 ( iN_E, iS, iN_X, Dnu, Inu_u_1, Inu_u_2, Inu_u_3, &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33, RN, RG1, RG2, RG3, LinearRN )

          IF( UseAcceptedEmAbRate ) RN(iC_EmAb) = AcceptedNewtonEmAbRate(iN_E,iS,iN_X)

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
    !$ACC PRESENT( Q_Matter, Q_SP )
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


  PURE SUBROUTINE ComputeChargedCurrentLeptonSource_Point &
    (Weights,Chi,Equilibrium,Distribution,Active,Source,AbsoluteSource)
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: Weights(:),Chi(:,:),Equilibrium(:,:),Distribution(:,:)
    LOGICAL, INTENT(in) :: Active
    REAL(DP), INTENT(out) :: Source,AbsoluteSource
    REAL(DP) :: Correction,Rate,SignedWeight
    INTEGER :: iE,iS

    Source=Zero
    Correction=Zero
    AbsoluteSource=Zero
    IF( .NOT. Active ) RETURN
    DO iE=1,SIZE(Weights)
    DO iS=1,2
      SignedWeight=REAL(2*iS-3,DP)*Weights(iE)
      Rate=(SignedWeight*Chi(iE,iS))*(Equilibrium(iE,iS)-Distribution(iE,iS))
      CALL AddMuonSourceCompensated(Rate,Source,Correction)
      AbsoluteSource=AbsoluteSource+ABS(Rate)
    END DO
    END DO
    Source=Source+Correction

  END SUBROUTINE ComputeChargedCurrentLeptonSource_Point


  PURE SUBROUTINE ComputeMuonLeptonSource_Point &
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

  END SUBROUTINE ComputeMuonLeptonSource_Point


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


  SUBROUTINE ComputeMuonLeptonSource_OrderV_Point &
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
    REAL(DP) :: Correction, Value
    INTEGER :: iProcess

    Source         = Zero
    AbsoluteSource = Zero
    IF( nSpecies < iNuM_Bar ) RETURN
    IF( .NOT. MatterSource_Mask(iYm,iN_X) ) RETURN

    TransferActive(1) = TransferRequested(1) .AND. Opacity_Mask(iOp_LFE,iN_X)
    TransferActive(2) = TransferRequested(2) .AND. Opacity_Mask(iOp_LFC,iN_X)
    TransferActive(3) = TransferRequested(3) .AND. Opacity_Mask(iOp_IMD,iN_X)

    CALL ComputeMuonLeptonSource_Point &
           ( W2_S, Chi_EmAb(:,iNuM:iNuM_Bar,iN_X), &
             Dnu_0(:,iNuM:iNuM_Bar,iN_X), Dnu(:,iNuM:iNuM_Bar,iN_X), &
             R_GL(:,iNuM:iNuM_Bar,iN_X,iGL_LFE:iGL_IMD), &
             Opacity_Mask(iOp_EmAb_Muon,iN_X), TransferActive, Source, AbsoluteSource )


  END SUBROUTINE ComputeMuonLeptonSource_OrderV_Point


  SUBROUTINE ComputeMuonLeptonSource_OrderV &
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
        CALL ComputeMuonLeptonSource_OrderV_Point &
               ( iN_X, Dnu, TransferRequested, Source(iN_X), AbsoluteSource(iN_X) )
    END DO

  END SUBROUTINE ComputeMuonLeptonSource_OrderV

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
            CALL ComputeMuonLeptonSource_OrderV_Point &
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

      IF( TraceNewtonSolve )THEN
        WRITE(FileName,'(A,".failure.rank",I6.6,".dat")') &
          TRIM(NewtonTraceFile), myid
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

          WRITE(*,*)                      '[SolveNeutrinoMatterCoupling] Non-coverged inner loop'
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

          WRITE(*,*)                      '[SolveNeutrinoMatterCoupling] Error'
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

    WRITE(*,*)                      '[SolveNeutrinoMatterCoupling] Convergence Status'
    WRITE(*,'(1x,4a8,5a23)')        'iQ', 'iN_X', 'k_outer', 'k_inner', 'Q_old', 'Q', 'U', 'F(iQ)', 'G(iQ)'
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iYe', iN_X, k_outer, k_inner, Ye0_P, Ye_P, U_Ye   (iN_X), Fm_outer(iYe,iN_X), G_Ye   (iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iYm', iN_X, k_outer, k_inner, Ym0_P, Ym_P, U_Ym   (iN_X), Fm_outer(iYm,iN_X), G_Ym   (iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iEf', iN_X, k_outer, k_inner,  E0_P,  E_P, U_Ef   (iN_X), Fm_outer(iEf,iN_X), G_Ef   (iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iV1', iN_X, k_outer, k_inner, V10_P, V1_P, U_V_d_1(iN_X), Fm_outer(iV1,iN_X), G_V_d_1(iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iV2', iN_X, k_outer, k_inner, V20_P, V2_P, U_V_d_2(iN_X), Fm_outer(iV2,iN_X), G_V_d_2(iN_X)
    WRITE(*,'(1x,a8,3i8,5es23.15)') 'iV3', iN_X, k_outer, k_inner, V30_P, V3_P, U_V_d_3(iN_X), Fm_outer(iV3,iN_X), G_V_d_3(iN_X)

  END SUBROUTINE PrintStatus_FP


END MODULE TwoMoment_NeutrinoMatterSolverModule
