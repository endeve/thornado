MODULE StandaloneMatterDiagnosticsModule

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Three, Half, FourPi, SqrtTiny
  USE UnitsModule, ONLY: &
    BoltzmannConstant, Centimeter, Erg, Gram, MeV, PlanckConstant, &
    Second, SpeedOfLight, UnitsActive
  USE ProgramHeaderModule, ONLY: &
    nDOFE, nDOFX, nDOFZ
  USE MeshModule, ONLY: &
    MeshE, MeshX, NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    WeightsX1, WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_X1_Dn, L_X1_Up
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep1, iGE_Ep2, iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_SqrtGm, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, iCF_Nm, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, &
    nAF, iAF_T, iAF_Ye, iAF_Ym
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, LimitEquilibriumDistributions_DG, &
    SetChargedCurrentMasks
  USE RadiationFieldsModule, ONLY: &
    nSpecies, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    iNuE_Bar, iNuM, iNuM_Bar, LeptonNumber, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    nGR, nDR, &
    iGR_D, iGR_I1, iGR_J, iGR_H1, &
    iDR_PL_Theta_1, iDR_PL_Theta_2, iDR_PL_dEnergy
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeMatterRHS_Collisions_OrderV, &
    ComputeMeanInverseMeanFreePaths_OrderV, &
    ComputeTransportOpacityFromConserved_OrderV, &
    ComputeIncrement_TwoMoment_Explicit_AP, &
    GetNeutrinoMatterSolverConvergenceDiagnostics, &
    nQ, nC, nConservationHandoffCheck, &
    iCHC_InputConversion_L, iCHC_PrimitivePair_L, &
    iCHC_FinalConversion_L, iCHC_MappedPair_L, iCHC_ReturnedPair_L, &
    iCHC_InputConversion_E, iCHC_PrimitivePair_E, &
    iCHC_FinalConversion_E, iCHC_MappedPair_E, iCHC_ReturnedPair_E
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    nLeptonFamilies, nLeptonFlavorCheck, &
    iLFC_InputConversion, iLFC_EntryProjection, &
    iLFC_CollisionSource, iLFC_RawMap, &
    iLFC_MapProjection, iLFC_FixedPointResidual, &
    iLFC_EquationResidual
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit, &
    UseSourceAwareLLF
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    iQ_Le, iQ_Lm, iQ_E, iQ_H, iQ_M1, iQ_M2, iQ_M3, &
    iC_EmAb, iC_Iso, &
    nEmAbRateDecomposition, &
    iEmAbRateDecomposition_NIn, iEmAbRateDecomposition_NOut, iEmAbRateDecomposition_NNet, &
    iEmAbRateDecomposition_EIn, iEmAbRateDecomposition_EOut, iEmAbRateDecomposition_QH, &
    nEmAbSpectrum, iEmAbSpectrum_Dnu, iEmAbSpectrum_DnuEq, &
    iEmAbSpectrum_Chi, iEmAbSpectrum_EIn, iEmAbSpectrum_EOut, &
    iEmAbSpectrum_QH, &
    nNNSRateDecomposition, &
    iNNSRateDecomposition_NIn, iNNSRateDecomposition_NOut, iNNSRateDecomposition_NNet, &
    iNNSRateDecomposition_EIn, iNNSRateDecomposition_EOut, iNNSRateDecomposition_QH
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment, ComputeConserved_TwoMoment, &
    FaceVelocity_X1, Flux_X1
  USE TwoMoment_TimeSteppingModule, ONLY: &
    iIMEX_Assembled, iIMEX_PreImplicit, iIMEX_PostImplicit, &
    iIMEX_PostImplicitLimiter, &
    iIMEX_CompletedStep, &
    GetImplicitStageCoefficient_IMEX_RK, &
    GetExplicitStageCoefficient_IMEX_RK, &
    GetExplicitStageWeight_IMEX_RK, GetImplicitStageWeight_IMEX_RK, &
    ExplicitRHSUsed_IMEX_RK, &
    HasFinalAssembly_IMEX_RK

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteStandaloneMatterDiagnostics
  PUBLIC :: InitializeMatterExchangeAttribution
  PUBLIC :: AccumulateMatterExchangeAttribution
  PUBLIC :: AccumulateActualMatterExchangeAttribution
  PUBLIC :: FinalizeMatterExchangeAttribution
  PUBLIC :: InitializeCumulativeConservationDiagnostics
  PUBLIC :: AccumulateCumulativeConservationDiagnostics
  PUBLIC :: FinalizeCumulativeConservationDiagnostics

  CHARACTER(6), PARAMETER :: ProcessName(nC) = [ CHARACTER(6) :: &
    'EmAb  ', 'Iso   ', 'NNS   ', 'NES   ', 'Pair  ', 'NuPair', 'Brem  ' ]

  INTEGER, PARAMETER :: iFace_EnergyLeft       = 1
  INTEGER, PARAMETER :: iFace_EnergyRight      = 2
  INTEGER, PARAMETER :: iFace_EnergyJump       = 3
  INTEGER, PARAMETER :: iFace_EnergyJumpAbs    = 4
  INTEGER, PARAMETER :: iFace_EnergyJumpRelL1  = 5
  INTEGER, PARAMETER :: iFace_FluxCentered     = 6
  INTEGER, PARAMETER :: iFace_FluxLLFPenalty   = 7
  INTEGER, PARAMETER :: iFace_FluxTotal        = 8
  INTEGER, PARAMETER :: iFace_EqEnergyLeft     = 9
  INTEGER, PARAMETER :: iFace_EqEnergyRight    = 10
  INTEGER, PARAMETER :: iFace_EqEnergyJump     = 11
  INTEGER, PARAMETER :: iFace_EqEnergyJumpAbs  = 12
  INTEGER, PARAMETER :: iFace_EqEnergyJumpRel  = 13
  INTEGER, PARAMETER :: iFace_NonEqEnergyLeft  = 14
  INTEGER, PARAMETER :: iFace_NonEqEnergyRight = 15
  INTEGER, PARAMETER :: iFace_NonEqEnergyJump  = 16
  INTEGER, PARAMETER :: iFace_NonEqEnergyJumpAbs = 17
  INTEGER, PARAMETER :: iFace_NonEqEnergyJumpRel = 18
  INTEGER, PARAMETER :: iFace_FluxLLFEq        = 19
  INTEGER, PARAMETER :: iFace_FluxLLFNonEq     = 20
  INTEGER, PARAMETER :: iFace_FluxLLFClosure   = 21
  INTEGER, PARAMETER :: nFaceDiagnostics       = 21

  INTEGER, PARAMETER :: iAttr_QL = 1
  INTEGER, PARAMETER :: iAttr_QH = 2
  INTEGER, PARAMETER :: nAttrQuantities = 2

  INTEGER, PARAMETER :: iAttr_Actual          = 1
  INTEGER, PARAMETER :: iAttr_NoStream        = 2
  INTEGER, PARAMETER :: iAttr_APMinusNoStream = 3
  INTEGER, PARAMETER :: iAttr_StandardMinusAP = 4
  INTEGER, PARAMETER :: iAttr_StandardTotal   = 5
  INTEGER, PARAMETER :: iAttr_APTotal         = 6
  INTEGER, PARAMETER :: iAttr_EmAbFull        = 7
  INTEGER, PARAMETER :: iAttr_EmAbMean        = 8
  INTEGER, PARAMETER :: iAttr_EmAbHigh        = 9
  INTEGER, PARAMETER :: nAttrComponents       = 9

  CHARACTER(20), PARAMETER :: AttributionComponentName(nAttrComponents) &
    = [ CHARACTER(20) :: &
        'actual', 'no_stream_shadow', 'ap_minus_no_stream', &
        'standard_minus_ap', &
        'standard_total', 'ap_total', 'emab_full', 'emab_mean', &
        'emab_high' ]

  LOGICAL, SAVE :: MatterExchangeAttributionInitialized = .FALSE.
  LOGICAL, SAVE :: MatterExchangeAttributionShadows = .TRUE.
  INTEGER, SAVE :: MatterExchangeAttributionCycle = 0
  REAL(DP), SAVE :: MatterExchangeAttributionTime = Zero
  CHARACTER(256), SAVE :: MatterExchangeAttributionFile = ''
  REAL(DP), ALLOCATABLE, SAVE :: &
    AttributionCumulative(:,:,:,:), AttributionAbsolute(:,:,:,:)

  LOGICAL, SAVE :: CumulativeConservationInitialized = .FALSE.
  CHARACTER(256), SAVE :: CumulativeConservationFile = ''
  REAL(DP), SAVE :: InitialMatterEnergy = Zero
  REAL(DP), SAVE :: InitialRadiationEnergy = Zero
  REAL(DP), SAVE :: InitialMatterLepton = Zero
  REAL(DP), SAVE :: InitialRadiationLepton = Zero
  REAL(DP), SAVE :: PreviousTotalEnergy = Zero
  REAL(DP), SAVE :: PreviousTotalLepton = Zero
  REAL(DP), SAVE :: PreviousMatterEnergy = Zero
  REAL(DP), SAVE :: PreviousRadiationEnergy = Zero
  REAL(DP), SAVE :: PreviousMatterLepton = Zero
  REAL(DP), SAVE :: PreviousRadiationLepton = Zero
  REAL(DP), SAVE :: CumulativeBoundaryEnergy = Zero
  REAL(DP), SAVE :: CumulativeBoundaryLepton = Zero

  INTEGER, PARAMETER :: iCheck_Energy = 1
  INTEGER, PARAMETER :: iCheck_Momentum1 = 2
  INTEGER, PARAMETER :: iCheck_Momentum2 = 3
  INTEGER, PARAMETER :: iCheck_Momentum3 = 4
  INTEGER, PARAMETER :: iCheck_LeptonElectron = 5
  INTEGER, PARAMETER :: iCheck_LeptonMuon = 6
  INTEGER, PARAMETER :: iCheck_LeptonTau = 7
  INTEGER, PARAMETER :: nCheckQuantities = 7
  CHARACTER(16), PARAMETER :: CheckQuantityName(nCheckQuantities) = [ &
    CHARACTER(16) :: 'energy', 'momentum1', 'momentum2', 'momentum3', &
                     'lepton_e', 'lepton_mu', 'lepton_tau' ]
  INTEGER, PARAMETER :: iCheck_MatterActual = 1
  INTEGER, PARAMETER :: iCheck_RadiationActual = 2
  INTEGER, PARAMETER :: iCheck_Boundary = 3
  INTEGER, PARAMETER :: iCheck_ImplicitMatter = 4
  INTEGER, PARAMETER :: iCheck_ImplicitRadiation = 5
  INTEGER, PARAMETER :: iCheck_ExplicitMatter = 6
  INTEGER, PARAMETER :: iCheck_ExplicitRadiation = 7
  INTEGER, PARAMETER :: iCheck_MatterRemainder = 8
  INTEGER, PARAMETER :: iCheck_RadiationRemainder = 9
  INTEGER, PARAMETER :: iCheck_CollisionResidual = 10
  INTEGER, PARAMETER :: iCheck_TransportBoundaryResidual = 11
  INTEGER, PARAMETER :: iCheck_ConservationResidual = 12
  INTEGER, PARAMETER :: iCheck_AlgebraicClosure = 13
  INTEGER, PARAMETER :: iCheck_PreLimiterMatterLocal = 14
  INTEGER, PARAMETER :: iCheck_PreLimiterRadiationLocal = 15
  INTEGER, PARAMETER :: iCheck_PostLimiterMatterLocal = 16
  INTEGER, PARAMETER :: iCheck_PostLimiterRadiationLocal = 17
  INTEGER, PARAMETER :: iCheck_DirectImplicitPairResidual = 18
  INTEGER, PARAMETER :: nCheckComponents = 18
  CHARACTER(32), PARAMETER :: CheckComponentName(nCheckComponents) = [ &
    CHARACTER(32) :: &
    'matter_actual', 'radiation_actual', 'boundary', &
    'implicit_matter', 'implicit_radiation', &
    'explicit_matter', 'explicit_radiation', &
    'matter_remainder', 'radiation_remainder', &
    'collision_pair_residual', 'transport_boundary_residual', &
    'conservation_residual', 'algebraic_closure', &
    'pre_stage_limiter_matter_local', &
    'pre_stage_limiter_rad_local', &
    'post_stage_limiter_matter_local', &
    'post_stage_limiter_rad_local', &
    'direct_implicit_pair_residual' ]
  CHARACTER(256), SAVE :: IMEXConservationCheckFile = ''
  CHARACTER(256), SAVE :: IMEXConservationStageDefectFile = ''
  CHARACTER(256), SAVE :: IMEXLeptonFlavorCheckFile = ''
  REAL(DP), SAVE :: &
    IMEXConservationCumulative(nCheckQuantities,nCheckComponents) = Zero

CONTAINS


  SUBROUTINE InitializeCumulativeConservationDiagnostics &
    ( FileName, Time, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP), INTENT(in) :: Time
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: Unit

    CALL ComputeConservedMatterRadiationTotals &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR, &
             InitialMatterEnergy, InitialRadiationEnergy, &
             InitialMatterLepton, InitialRadiationLepton )

    PreviousTotalEnergy = InitialMatterEnergy + InitialRadiationEnergy
    PreviousTotalLepton = InitialMatterLepton + InitialRadiationLepton
    PreviousMatterEnergy = InitialMatterEnergy
    PreviousRadiationEnergy = InitialRadiationEnergy
    PreviousMatterLepton = InitialMatterLepton
    PreviousRadiationLepton = InitialRadiationLepton
    CumulativeBoundaryEnergy = Zero
    CumulativeBoundaryLepton = Zero
    CumulativeConservationFile = TRIM( FileName )
    IMEXConservationCheckFile = TRIM( FileName ) // '.imex'
    IMEXConservationStageDefectFile &
      = TRIM( FileName ) // '.imex_implicit_pair_defect'
    IMEXLeptonFlavorCheckFile &
      = TRIM( FileName ) // '.imex_lepton_flavor_check'
    IMEXConservationCumulative = Zero

    OPEN( NEWUNIT = Unit, FILE = TRIM(CumulativeConservationFile), &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# cycle time_s dt_s matter_energy_erg radiation_energy_erg' // &
      ' total_energy_erg cumulative_boundary_energy_erg' // &
      ' cumulative_energy_residual_erg cumulative_energy_relative' // &
      ' cumulative_energy_exchange_relative' // &
      ' step_energy_change_erg step_boundary_energy_erg' // &
      ' step_energy_residual_erg step_energy_relative' // &
      ' step_energy_exchange_relative' // &
      ' matter_electron_lepton radiation_electron_lepton' // &
      ' total_electron_lepton cumulative_boundary_electron_lepton' // &
      ' cumulative_electron_lepton_residual' // &
      ' cumulative_electron_lepton_relative' // &
      ' cumulative_electron_lepton_exchange_relative' // &
      ' step_electron_lepton_change step_boundary_electron_lepton' // &
      ' step_electron_lepton_residual step_electron_lepton_relative' // &
      ' step_electron_lepton_exchange_relative'
    WRITE(Unit,'(A)') &
      '# residual = state change - signed off-grid increment'
    WRITE(Unit,'(A)') &
      '# cumulative relative uses the initial total; step relative uses' // &
      ' the adjacent-step total'
    WRITE(Unit,'(A)') &
      '# exchange relative uses the sum of absolute matter, radiation,' // &
      ' and boundary changes'
    WRITE(Unit,'(I8,1X,26(ES24.16E3,1X))') &
      0, Time / Second, Zero, &
      InitialMatterEnergy / Erg, InitialRadiationEnergy / Erg, &
      PreviousTotalEnergy / Erg, Zero, Zero, Zero, Zero, &
      Zero, Zero, Zero, Zero, Zero, &
      InitialMatterLepton, InitialRadiationLepton, PreviousTotalLepton, &
      Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero
    CLOSE(Unit)

    OPEN( NEWUNIT = Unit, FILE = TRIM(IMEXLeptonFlavorCheckFile), &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# cycle time_s dt_s stage diagonal_coefficient implicit_weight' // &
      ' family input_conversion entry_projection collision_source' // &
      ' raw_map map_projection' // &
      ' fixed_point_residual equation_residual' // &
      ' reconstructed_radiation_change' // &
      ' direct_radiation_change identity_closure direct_pair_residual'
    WRITE(Unit,'(A)') &
      '# changes are particle numbers over the diagonal implicit stage;' // &
      ' multiply by implicit_weight/diagonal_coefficient for the' // &
      ' full-step IMEX contribution.'
    WRITE(Unit,'(A)') &
      '# reconstructed = input_conversion + entry_projection + raw_map' // &
      ' + map_projection - fixed_point_residual; collision_source is an independent' // &
      ' final-state source integral.'
    CLOSE(Unit)

    OPEN( NEWUNIT = Unit, FILE = TRIM(IMEXConservationCheckFile), &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# cycle time_s dt_s quantity component step_value cumulative_value'
    WRITE(Unit,'(A)') &
      '# energy is in erg; momentum is in code/cgs momentum units;' // &
      ' lepton quantities are particle numbers.'
    WRITE(Unit,'(A)') &
      '# radiation energy and momentum use the conserved-moment order-v' // &
      ' combinations N+v^i G_i and G_i+v_i N used by the boundary tally.'
    WRITE(Unit,'(A)') &
      '# curvilinear momentum requires geometric source terms for a global' // &
      ' balance; collision_pair_residual remains a valid local-exchange test.'
    WRITE(Unit,'(A)') &
      '# remainders contain final assembly, limiter, and frame cross terms' // &
      ' not represented by the weighted stage RHS contributions.'
    WRITE(Unit,'(A)') &
      '# pre/post-stage limiter entries are local stage responses and are' // &
      ' not additional terms in the exact full-step decomposition.'
    CLOSE(Unit)

    OPEN( NEWUNIT = Unit, FILE = TRIM(IMEXConservationStageDefectFile), &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# cycle time_s dt_s stage diagonal_coefficient implicit_weight' // &
      ' quantity direct_pair_residual weighted_pair_residual'
    WRITE(Unit,'(A)') &
      '# pair defect is the change in matter plus radiation across the' // &
      ' accepted diagonal implicit stage; weighted values multiply it' // &
      ' by w_IM/a_IM(stage,stage).'
    CLOSE(Unit)

    CumulativeConservationInitialized = .TRUE.

  END SUBROUTINE InitializeCumulativeConservationDiagnostics


  SUBROUTINE AccumulateCumulativeConservationDiagnostics &
    ( Cycle, Time, TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
      GE, GX, CF, CR, OffGridFluidIncrement, OffGridRadiationIncrement, &
      OffGridRadiationLeptonFlavorIncrement, &
      StageFluidState, StageRadiationState, &
      StageFluidImplicitRate, StageFluidExplicitRate, &
      StageRadiationImplicitRate, StageRadiationExplicitRate, &
      StageLeptonFlavorCheckRate )

    INTEGER, INTENT(in) :: Cycle
    REAL(DP), INTENT(in) :: Time, TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: OffGridFluidIncrement(1:nCF)
    REAL(DP), INTENT(in) :: OffGridRadiationIncrement(1:2*nCR)
    REAL(DP), INTENT(in) :: OffGridRadiationLeptonFlavorIncrement(3)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageFluidImplicitRate(:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageFluidExplicitRate(:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationImplicitRate(:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationExplicitRate(:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageLeptonFlavorCheckRate(:,:,:,:,:,:,:)

    INTEGER :: Unit
    REAL(DP) :: MatterEnergy, RadiationEnergy, TotalEnergy
    REAL(DP) :: MatterLepton, RadiationLepton, TotalLepton
    REAL(DP) :: BoundaryEnergy, BoundaryLepton, hc3
    REAL(DP) :: StepEnergyChange, StepEnergyResidual
    REAL(DP) :: StepLeptonChange, StepLeptonResidual
    REAL(DP) :: CumulativeEnergyResidual, CumulativeLeptonResidual
    REAL(DP) :: StepEnergyRelative, StepLeptonRelative
    REAL(DP) :: CumulativeEnergyRelative, CumulativeLeptonRelative
    REAL(DP) :: StepEnergyExchangeRelative, StepLeptonExchangeRelative
    REAL(DP) :: &
      CumulativeEnergyExchangeRelative, CumulativeLeptonExchangeRelative
    REAL(DP) :: PreviousMatterEnergyLocal, PreviousRadiationEnergyLocal
    REAL(DP) :: PreviousMatterLeptonLocal, PreviousRadiationLeptonLocal
    REAL(DP) :: MatterComponents(nCheckQuantities)
    REAL(DP) :: RadiationComponents(nCheckQuantities)
    REAL(DP) :: MatterComponentsInitial(nCheckQuantities)
    REAL(DP) :: RadiationComponentsInitial(nCheckQuantities)
    REAL(DP) :: MatterChange(nCheckQuantities)
    REAL(DP) :: RadiationChange(nCheckQuantities)
    REAL(DP) :: BoundaryComponents(nCheckQuantities)

    IF( .NOT. CumulativeConservationInitialized ) &
      ERROR STOP 'Cumulative conservation diagnostics are not initialized.'

    PreviousMatterEnergyLocal    = PreviousMatterEnergy
    PreviousRadiationEnergyLocal = PreviousRadiationEnergy
    PreviousMatterLeptonLocal    = PreviousMatterLepton
    PreviousRadiationLeptonLocal = PreviousRadiationLepton

    CALL ComputeConservedMatterRadiationTotals &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR, &
             MatterEnergy, RadiationEnergy, MatterLepton, RadiationLepton )

    CALL ComputeConservedMatterRadiationComponents &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR, &
             MatterComponents, RadiationComponents )

    hc3 = One
    IF( UnitsActive ) hc3 = ( PlanckConstant * SpeedOfLight )**3

    BoundaryEnergy = OffGridFluidIncrement(iCF_E) &
      + FourPi * OffGridRadiationIncrement(nCR+iCR_N) / hc3
    BoundaryLepton = OffGridFluidIncrement(iCF_Ne) &
      + FourPi * OffGridRadiationLeptonFlavorIncrement(1) / hc3

    BoundaryComponents = Zero
    BoundaryComponents(iCheck_Energy) = BoundaryEnergy
    BoundaryComponents(iCheck_Momentum1) &
      = OffGridFluidIncrement(iCF_S1) &
          + FourPi * OffGridRadiationIncrement(nCR+iCR_G1) / hc3
    BoundaryComponents(iCheck_Momentum2) &
      = OffGridFluidIncrement(iCF_S2) &
          + FourPi * OffGridRadiationIncrement(nCR+iCR_G2) / hc3
    BoundaryComponents(iCheck_Momentum3) &
      = OffGridFluidIncrement(iCF_S3) &
          + FourPi * OffGridRadiationIncrement(nCR+iCR_G3) / hc3
    BoundaryComponents(iCheck_LeptonElectron) = BoundaryLepton
    BoundaryComponents(iCheck_LeptonMuon) &
      = OffGridFluidIncrement(iCF_Nm) &
          + FourPi * OffGridRadiationLeptonFlavorIncrement(2) / hc3
    BoundaryComponents(iCheck_LeptonTau) &
      = FourPi * OffGridRadiationLeptonFlavorIncrement(3) / hc3

    CALL ComputeConservedMatterRadiationComponents &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             StageFluidState(:,:,:,:,:,1,iIMEX_Assembled), &
             StageRadiationState(:,:,:,:,:,:,:,1,iIMEX_Assembled), &
             MatterComponentsInitial, RadiationComponentsInitial )
    MatterChange = MatterComponents - MatterComponentsInitial
    RadiationChange = RadiationComponents - RadiationComponentsInitial

    TotalEnergy = MatterEnergy + RadiationEnergy
    TotalLepton = MatterLepton + RadiationLepton
    StepEnergyChange = TotalEnergy - PreviousTotalEnergy
    StepLeptonChange = TotalLepton - PreviousTotalLepton
    StepEnergyResidual = StepEnergyChange - BoundaryEnergy
    StepLeptonResidual = StepLeptonChange - BoundaryLepton

    CumulativeBoundaryEnergy &
      = CumulativeBoundaryEnergy + BoundaryEnergy
    CumulativeBoundaryLepton &
      = CumulativeBoundaryLepton + BoundaryLepton
    CumulativeEnergyResidual = TotalEnergy &
      - ( InitialMatterEnergy + InitialRadiationEnergy ) &
      - CumulativeBoundaryEnergy
    CumulativeLeptonResidual = TotalLepton &
      - ( InitialMatterLepton + InitialRadiationLepton ) &
      - CumulativeBoundaryLepton

    StepEnergyRelative = ABS( StepEnergyResidual ) &
      / MAX( ABS(PreviousTotalEnergy), ABS(TotalEnergy), TINY(One) )
    StepLeptonRelative = ABS( StepLeptonResidual ) &
      / MAX( ABS(PreviousTotalLepton), ABS(TotalLepton), TINY(One) )
    CumulativeEnergyRelative = ABS( CumulativeEnergyResidual ) &
      / MAX( ABS(InitialMatterEnergy + InitialRadiationEnergy), TINY(One) )
    CumulativeLeptonRelative = ABS( CumulativeLeptonResidual ) &
      / MAX( ABS(InitialMatterLepton + InitialRadiationLepton), TINY(One) )
    StepEnergyExchangeRelative = ABS( StepEnergyResidual ) &
      / MAX( ABS(MatterEnergy - PreviousMatterEnergy) &
             + ABS(RadiationEnergy - PreviousRadiationEnergy) &
             + ABS(BoundaryEnergy), TINY(One) )
    StepLeptonExchangeRelative = ABS( StepLeptonResidual ) &
      / MAX( ABS(MatterLepton - PreviousMatterLepton) &
             + ABS(RadiationLepton - PreviousRadiationLepton) &
             + ABS(BoundaryLepton), TINY(One) )
    CumulativeEnergyExchangeRelative = ABS( CumulativeEnergyResidual ) &
      / MAX( ABS(MatterEnergy - InitialMatterEnergy) &
             + ABS(RadiationEnergy - InitialRadiationEnergy) &
             + ABS(CumulativeBoundaryEnergy), TINY(One) )
    CumulativeLeptonExchangeRelative = ABS( CumulativeLeptonResidual ) &
      / MAX( ABS(MatterLepton - InitialMatterLepton) &
             + ABS(RadiationLepton - InitialRadiationLepton) &
             + ABS(CumulativeBoundaryLepton), TINY(One) )

    OPEN( NEWUNIT = Unit, FILE = TRIM(CumulativeConservationFile), &
          STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE', &
          FORM = 'FORMATTED' )
    WRITE(Unit,'(I8,1X,26(ES24.16E3,1X))') &
      Cycle, Time / Second, TimeStep / Second, &
      MatterEnergy / Erg, RadiationEnergy / Erg, TotalEnergy / Erg, &
      CumulativeBoundaryEnergy / Erg, CumulativeEnergyResidual / Erg, &
      CumulativeEnergyRelative, CumulativeEnergyExchangeRelative, &
      StepEnergyChange / Erg, &
      BoundaryEnergy / Erg, StepEnergyResidual / Erg, StepEnergyRelative, &
      StepEnergyExchangeRelative, &
      MatterLepton, RadiationLepton, TotalLepton, &
      CumulativeBoundaryLepton, CumulativeLeptonResidual, &
      CumulativeLeptonRelative, CumulativeLeptonExchangeRelative, &
      StepLeptonChange, BoundaryLepton, StepLeptonResidual, &
      StepLeptonRelative, StepLeptonExchangeRelative
    CLOSE(Unit)

    CALL AccumulateIMEXConservationCheck &
           ( Cycle, Time, TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             GE, GX, &
             MatterChange, RadiationChange, BoundaryComponents, &
             StageFluidState, StageRadiationState, &
             StageFluidImplicitRate, StageFluidExplicitRate, &
             StageRadiationImplicitRate, StageRadiationExplicitRate, &
             StageLeptonFlavorCheckRate )

    PreviousTotalEnergy = TotalEnergy
    PreviousTotalLepton = TotalLepton
    PreviousMatterEnergy = MatterEnergy
    PreviousRadiationEnergy = RadiationEnergy
    PreviousMatterLepton = MatterLepton
    PreviousRadiationLepton = RadiationLepton

  END SUBROUTINE AccumulateCumulativeConservationDiagnostics


  SUBROUTINE AccumulateIMEXConservationCheck &
    ( Cycle, Time, TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      MatterChange, RadiationChange, BoundaryChange, &
      StageFluidState, StageRadiationState, &
      StageFluidImplicitRate, StageFluidExplicitRate, &
      StageRadiationImplicitRate, StageRadiationExplicitRate, &
      StageLeptonFlavorCheckRate )

    INTEGER, INTENT(in) :: Cycle
    REAL(DP), INTENT(in) :: Time, TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      MatterChange(nCheckQuantities), &
      RadiationChange(nCheckQuantities), &
      BoundaryChange(nCheckQuantities)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageFluidImplicitRate(:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageFluidExplicitRate(:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationImplicitRate(:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationExplicitRate(:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageLeptonFlavorCheckRate(:,:,:,:,:,:,:)

    INTEGER :: iComponent, iStage, iQuantity, Unit
    REAL(DP) :: aII, wEX, wIM
    REAL(DP) :: dMatter(nCheckQuantities)
    REAL(DP) :: dRadiation(nCheckQuantities)
    REAL(DP) :: DirectPairDefect(nCheckQuantities)
    REAL(DP) :: WeightedPairDefect(nCheckQuantities)
    REAL(DP) :: MatterA(nCheckQuantities)
    REAL(DP) :: RadiationA(nCheckQuantities)
    REAL(DP) :: MatterB(nCheckQuantities)
    REAL(DP) :: RadiationB(nCheckQuantities)
    REAL(DP) :: Step(nCheckQuantities,nCheckComponents)
    REAL(DP) :: LeptonCheck(nLeptonFamilies,nLeptonFlavorCheck)
    REAL(DP), ALLOCATABLE :: CF_Work(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR_Work(:,:,:,:,:,:,:)

    ALLOCATE &
      ( CF_Work(1:nDOFX, &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCF) )
    ALLOCATE &
      ( CR_Work(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

    Step = Zero
    Step(:,iCheck_MatterActual) = MatterChange
    Step(:,iCheck_RadiationActual) = RadiationChange
    Step(:,iCheck_Boundary) = BoundaryChange

    DO iStage = 1, SIZE( StageFluidState, 6 )

      wIM = GetImplicitStageWeight_IMEX_RK( iStage )
      aII = GetImplicitStageCoefficient_IMEX_RK( iStage, iStage )

      DirectPairDefect = Zero
      WeightedPairDefect = Zero

      IF( aII /= Zero )THEN

        CALL ComputeConservedMatterRadiationComponents &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                 StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
                 StageRadiationState &
                   (:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit), &
                 MatterA, RadiationA )
        CALL ComputeConservedMatterRadiationComponents &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                 StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
                 StageRadiationState &
                   (:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit), &
                 MatterB, RadiationB )

        DirectPairDefect = MatterB + RadiationB - MatterA - RadiationA
        WeightedPairDefect = ( wIM / aII ) * DirectPairDefect

        CALL ComputeIntegratedLeptonFlavorCheck &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
                 StageLeptonFlavorCheckRate(:,:,:,:,:,:,iStage), &
                 LeptonCheck )
        CALL WriteLeptonFlavorStageCheck &
               ( Cycle, Time, TimeStep, iStage, aII, wIM, &
                 LeptonCheck, MatterA, RadiationA, MatterB, RadiationB )

      END IF

      Step(:,iCheck_DirectImplicitPairResidual) &
        = Step(:,iCheck_DirectImplicitPairResidual) + WeightedPairDefect

      OPEN( NEWUNIT = Unit, FILE = TRIM(IMEXConservationStageDefectFile), &
            STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE', &
            FORM = 'FORMATTED' )
      DO iQuantity = 1, nCheckQuantities
        WRITE(Unit,'(I8,1X,2(ES24.16E3,1X),I4,1X,2(ES24.16E3,1X),' // &
                   'A,1X,2(ES24.16E3,1X))') &
          Cycle, Time / Second, TimeStep / Second, iStage, aII, wIM, &
          TRIM(CheckQuantityName(iQuantity)), &
          ConvertCheckQuantity( iQuantity, DirectPairDefect(iQuantity) ), &
          ConvertCheckQuantity( iQuantity, WeightedPairDefect(iQuantity) )
      END DO
      CLOSE(Unit)

      CALL ComputeWeightedStageConservationChange &
             ( TimeStep, wIM, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               StageFluidImplicitRate(:,:,:,:,:,iStage), &
               StageRadiationImplicitRate(:,:,:,:,:,:,:,iStage), &
               CF_Work, CR_Work, &
               dMatter, dRadiation )
      Step(:,iCheck_ImplicitMatter) &
        = Step(:,iCheck_ImplicitMatter) + dMatter
      Step(:,iCheck_ImplicitRadiation) &
        = Step(:,iCheck_ImplicitRadiation) + dRadiation

      wEX = GetExplicitStageWeight_IMEX_RK( iStage )
      CALL ComputeWeightedStageConservationChange &
             ( TimeStep, wEX, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               StageRadiationState &
                 (:,:,:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               StageFluidExplicitRate(:,:,:,:,:,iStage), &
               StageRadiationExplicitRate(:,:,:,:,:,:,:,iStage), &
               CF_Work, CR_Work, &
               dMatter, dRadiation )
      Step(:,iCheck_ExplicitMatter) &
        = Step(:,iCheck_ExplicitMatter) + dMatter
      Step(:,iCheck_ExplicitRadiation) &
        = Step(:,iCheck_ExplicitRadiation) + dRadiation

      CALL ComputeConservedMatterRadiationComponents &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_Assembled), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_Assembled), &
               MatterA, RadiationA )
      CALL ComputeConservedMatterRadiationComponents &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               MatterB, RadiationB )
      Step(:,iCheck_PreLimiterMatterLocal) &
        = Step(:,iCheck_PreLimiterMatterLocal) + MatterB - MatterA
      Step(:,iCheck_PreLimiterRadiationLocal) &
        = Step(:,iCheck_PreLimiterRadiationLocal) + RadiationB - RadiationA

      CALL ComputeConservedMatterRadiationComponents &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               MatterA, RadiationA )
      CALL ComputeConservedMatterRadiationComponents &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               StageFluidState &
                 (:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               StageRadiationState &
                 (:,:,:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               MatterB, RadiationB )
      Step(:,iCheck_PostLimiterMatterLocal) &
        = Step(:,iCheck_PostLimiterMatterLocal) + MatterB - MatterA
      Step(:,iCheck_PostLimiterRadiationLocal) &
        = Step(:,iCheck_PostLimiterRadiationLocal) + RadiationB - RadiationA

    END DO

    Step(:,iCheck_MatterRemainder) &
      = Step(:,iCheck_MatterActual) &
          - Step(:,iCheck_ImplicitMatter) &
          - Step(:,iCheck_ExplicitMatter)
    Step(:,iCheck_RadiationRemainder) &
      = Step(:,iCheck_RadiationActual) &
          - Step(:,iCheck_ImplicitRadiation) &
          - Step(:,iCheck_ExplicitRadiation)
    Step(:,iCheck_CollisionResidual) &
      = Step(:,iCheck_ImplicitMatter) &
          + Step(:,iCheck_ImplicitRadiation)
    Step(:,iCheck_TransportBoundaryResidual) &
      = Step(:,iCheck_ExplicitMatter) &
          + Step(:,iCheck_ExplicitRadiation) &
          - Step(:,iCheck_Boundary)
    Step(:,iCheck_ConservationResidual) &
      = Step(:,iCheck_MatterActual) &
          + Step(:,iCheck_RadiationActual) &
          - Step(:,iCheck_Boundary)
    Step(:,iCheck_AlgebraicClosure) &
      = Step(:,iCheck_ConservationResidual) &
          - Step(:,iCheck_CollisionResidual) &
          - Step(:,iCheck_TransportBoundaryResidual) &
          - Step(:,iCheck_MatterRemainder) &
          - Step(:,iCheck_RadiationRemainder)

    IMEXConservationCumulative = IMEXConservationCumulative + Step

    OPEN( NEWUNIT = Unit, FILE = TRIM(IMEXConservationCheckFile), &
          STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE', &
          FORM = 'FORMATTED' )
    DO iQuantity = 1, nCheckQuantities
      DO iComponent = 1, nCheckComponents
        WRITE(Unit,'(I8,1X,2(ES24.16E3,1X),2(A,1X),2(ES24.16E3,1X))') &
          Cycle, Time / Second, TimeStep / Second, &
          TRIM(CheckQuantityName(iQuantity)), &
          TRIM(CheckComponentName(iComponent)), &
          ConvertCheckQuantity( iQuantity, Step(iQuantity,iComponent) ), &
          ConvertCheckQuantity &
            ( iQuantity, IMEXConservationCumulative(iQuantity,iComponent) )
      END DO
    END DO
    CLOSE(Unit)

    DEALLOCATE( CF_Work, CR_Work )

  END SUBROUTINE AccumulateIMEXConservationCheck


  SUBROUTINE ComputeIntegratedLeptonFlavorCheck &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CheckRate, IntegratedRate )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CheckRate(1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
                iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                iZ_B1(4):iZ_E1(4))
    REAL(DP), INTENT(out) :: &
      IntegratedRate(nLeptonFamilies,nLeptonFlavorCheck)

    INTEGER :: iCheck, iFamily, iNodeX, iX1, iX2, iX3
    REAL(DP) :: SpaceWeight

    IntegratedRate = Zero

    DO iX3 = iZ_B0(4), iZ_E0(4)
    DO iX2 = iZ_B0(3), iZ_E0(3)
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      SpaceWeight = MeshX(1) % Width(iX1) &
                  * MeshX(2) % Width(iX2) &
                  * MeshX(3) % Width(iX3) &
                  * WeightsX_q(iNodeX) &
                  * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

      DO iCheck = 1, nLeptonFlavorCheck
      DO iFamily = 1, nLeptonFamilies
        IntegratedRate(iFamily,iCheck) &
          = IntegratedRate(iFamily,iCheck) &
              + SpaceWeight &
                  * CheckRate(iFamily,iCheck,iNodeX,iX1,iX2,iX3)
      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIntegratedLeptonFlavorCheck


  SUBROUTINE WriteLeptonFlavorStageCheck &
    ( Cycle, Time, TimeStep, iStage, aII, wIM, CheckRate, &
      MatterA, RadiationA, MatterB, RadiationB )

    INTEGER, INTENT(in) :: Cycle, iStage
    REAL(DP), INTENT(in) :: Time, TimeStep, aII, wIM
    REAL(DP), INTENT(in) :: &
      CheckRate(nLeptonFamilies,nLeptonFlavorCheck)
    REAL(DP), INTENT(in) :: &
      MatterA(nCheckQuantities), RadiationA(nCheckQuantities), &
      MatterB(nCheckQuantities), RadiationB(nCheckQuantities)

    CHARACTER(8), PARAMETER :: FamilyName(nLeptonFamilies) &
      = [ CHARACTER(8) :: 'electron', 'muon', 'tau' ]
    INTEGER :: iFamily, iQuantity, Unit
    REAL(DP) :: InputConversion, EntryProjection, CollisionSource
    REAL(DP) :: RawMap, MapProjection
    REAL(DP) :: FixedPointResidual, EquationResidual
    REAL(DP) :: Reconstructed, DirectRadiation
    REAL(DP) :: IdentityClosure, DirectPairResidual

    OPEN( NEWUNIT = Unit, FILE = TRIM(IMEXLeptonFlavorCheckFile), &
          STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE', &
          FORM = 'FORMATTED' )

    DO iFamily = 1, nLeptonFamilies

      iQuantity = iCheck_LeptonElectron + iFamily - 1
      InputConversion &
        = TimeStep * aII * CheckRate(iFamily,iLFC_InputConversion)
      EntryProjection &
        = TimeStep * aII * CheckRate(iFamily,iLFC_EntryProjection)
      CollisionSource &
        = TimeStep * aII * CheckRate(iFamily,iLFC_CollisionSource)
      RawMap = TimeStep * aII * CheckRate(iFamily,iLFC_RawMap)
      MapProjection &
        = TimeStep * aII * CheckRate(iFamily,iLFC_MapProjection)
      FixedPointResidual &
        = TimeStep * aII * CheckRate(iFamily,iLFC_FixedPointResidual)
      EquationResidual &
        = TimeStep * aII * CheckRate(iFamily,iLFC_EquationResidual)
      Reconstructed &
        = InputConversion + EntryProjection + RawMap + MapProjection &
            - FixedPointResidual
      DirectRadiation = RadiationB(iQuantity) - RadiationA(iQuantity)
      IdentityClosure = DirectRadiation - Reconstructed
      DirectPairResidual = DirectRadiation &
        + MatterB(iQuantity) - MatterA(iQuantity)

      WRITE(Unit,'(I8,1X,2(ES24.16E3,1X),I4,1X,2(ES24.16E3,1X),' // &
                 'A,1X,11(ES24.16E3,1X))') &
        Cycle, Time / Second, TimeStep / Second, iStage, aII, wIM, &
        TRIM(FamilyName(iFamily)), InputConversion, EntryProjection, &
        CollisionSource, &
        RawMap, MapProjection, FixedPointResidual, EquationResidual, &
        Reconstructed, &
        DirectRadiation, IdentityClosure, DirectPairResidual

    END DO

    CLOSE(Unit)

  END SUBROUTINE WriteLeptonFlavorStageCheck


  SUBROUTINE ComputeWeightedStageConservationChange &
    ( TimeStep, Weight, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      CF_Base, CR_Base, CF_Rate, CR_Rate, CF_Work, CR_Work, &
      MatterChange, RadiationChange )

    REAL(DP), INTENT(in) :: TimeStep, Weight
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_Base(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
              iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR_Base(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
              iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      CF_Rate(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
              iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR_Rate(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
              iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      CF_Work(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
              iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      CR_Work(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
              iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      MatterChange(nCheckQuantities), RadiationChange(nCheckQuantities)

    REAL(DP) :: MatterBase(nCheckQuantities)
    REAL(DP) :: RadiationBase(nCheckQuantities)
    REAL(DP) :: MatterWork(nCheckQuantities)
    REAL(DP) :: RadiationWork(nCheckQuantities)

    MatterChange = Zero
    RadiationChange = Zero
    IF( Weight == Zero ) RETURN

    CF_Work = CF_Base + TimeStep * Weight * CF_Rate
    CR_Work = CR_Base + TimeStep * Weight * CR_Rate

    CALL ComputeConservedMatterRadiationComponents &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF_Base, CR_Base, &
             MatterBase, RadiationBase )
    CALL ComputeConservedMatterRadiationComponents &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF_Work, CR_Work, &
             MatterWork, RadiationWork )

    MatterChange = MatterWork - MatterBase
    RadiationChange = RadiationWork - RadiationBase

  END SUBROUTINE ComputeWeightedStageConservationChange


  SUBROUTINE FinalizeCumulativeConservationDiagnostics

    CumulativeConservationInitialized = .FALSE.
    CumulativeConservationFile = ''
    InitialMatterEnergy = Zero
    InitialRadiationEnergy = Zero
    InitialMatterLepton = Zero
    InitialRadiationLepton = Zero
    PreviousTotalEnergy = Zero
    PreviousTotalLepton = Zero
    PreviousMatterEnergy = Zero
    PreviousRadiationEnergy = Zero
    PreviousMatterLepton = Zero
    PreviousRadiationLepton = Zero
    CumulativeBoundaryEnergy = Zero
    CumulativeBoundaryLepton = Zero
    IMEXConservationCheckFile = ''
    IMEXConservationStageDefectFile = ''
    IMEXConservationCumulative = Zero

  END SUBROUTINE FinalizeCumulativeConservationDiagnostics


  SUBROUTINE ComputeConservedMatterRadiationTotals &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR, &
      MatterEnergy, RadiationEnergy, MatterLepton, RadiationLepton )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      MatterEnergy, RadiationEnergy, MatterLepton, RadiationLepton

    REAL(DP) :: Matter(nCheckQuantities)
    REAL(DP) :: Radiation(nCheckQuantities)

    CALL ComputeConservedMatterRadiationComponents &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR, &
             Matter, Radiation )

    MatterEnergy = Matter(iCheck_Energy)
    RadiationEnergy = Radiation(iCheck_Energy)
    MatterLepton = Matter(iCheck_LeptonElectron)
    RadiationLepton = Radiation(iCheck_LeptonElectron)

  END SUBROUTINE ComputeConservedMatterRadiationTotals


  SUBROUTINE ComputeConservedMatterRadiationComponents &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF, CR, Matter, Radiation )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      Matter(nCheckQuantities), Radiation(nCheckQuantities)

    INTEGER :: iE, iFamily, iNodeE, iNodeX, iNodeZ, iS, iX1, iX2, iX3
    REAL(DP) :: hc3, PhaseWeight2, PhaseWeight3, SpaceWeight
    REAL(DP) :: V_u_1, V_u_2, V_u_3, V_d_1, V_d_2, V_d_3
    REAL(DP) :: Enu, Fnu_d_1, Fnu_d_2, Fnu_d_3

    hc3 = One
    IF( UnitsActive ) hc3 = ( PlanckConstant * SpeedOfLight )**3
    Matter = Zero
    Radiation = Zero

    DO iX3 = iZ_B0(4), iZ_E0(4)
    DO iX2 = iZ_B0(3), iZ_E0(3)
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      SpaceWeight = MeshX(1) % Width(iX1) &
                  * MeshX(2) % Width(iX2) &
                  * MeshX(3) % Width(iX3) &
                  * WeightsX_q(iNodeX) &
                  * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

      Matter(iCheck_Energy) = Matter(iCheck_Energy) &
        + SpaceWeight * CF(iNodeX,iX1,iX2,iX3,iCF_E)
      Matter(iCheck_Momentum1) = Matter(iCheck_Momentum1) &
        + SpaceWeight * CF(iNodeX,iX1,iX2,iX3,iCF_S1)
      Matter(iCheck_Momentum2) = Matter(iCheck_Momentum2) &
        + SpaceWeight * CF(iNodeX,iX1,iX2,iX3,iCF_S2)
      Matter(iCheck_Momentum3) = Matter(iCheck_Momentum3) &
        + SpaceWeight * CF(iNodeX,iX1,iX2,iX3,iCF_S3)
      Matter(iCheck_LeptonElectron) = Matter(iCheck_LeptonElectron) &
        + SpaceWeight * CF(iNodeX,iX1,iX2,iX3,iCF_Ne)
      Matter(iCheck_LeptonMuon) = Matter(iCheck_LeptonMuon) &
        + SpaceWeight * CF(iNodeX,iX1,iX2,iX3,iCF_Nm)

      V_u_1 = CF(iNodeX,iX1,iX2,iX3,iCF_S1) &
        / ( GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
            * CF(iNodeX,iX1,iX2,iX3,iCF_D) )
      V_u_2 = CF(iNodeX,iX1,iX2,iX3,iCF_S2) &
        / ( GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
            * CF(iNodeX,iX1,iX2,iX3,iCF_D) )
      V_u_3 = CF(iNodeX,iX1,iX2,iX3,iCF_S3) &
        / ( GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) &
            * CF(iNodeX,iX1,iX2,iX3,iCF_D) )
      V_d_1 = GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) * V_u_1
      V_d_2 = GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) * V_u_2
      V_d_3 = GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) * V_u_3

      DO iS = 1, nSpecies
        iFamily = ( iS + 1 ) / 2
        DO iE = iZ_B0(1), iZ_E0(1)
        DO iNodeE = 1, nDOFE

          iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
          PhaseWeight2 = FourPi / hc3 * MeshE % Width(iE) &
            * WeightsE(iNodeE) * GE(iNodeE,iE,iGE_Ep2)
          PhaseWeight3 = FourPi / hc3 * MeshE % Width(iE) &
            * WeightsE(iNodeE) * GE(iNodeE,iE,iGE_Ep3)

          Radiation(iCheck_LeptonElectron+iFamily-1) &
            = Radiation(iCheck_LeptonElectron+iFamily-1) &
                + SpaceWeight * PhaseWeight2 * LeptonNumber(iS) &
                  * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

          ! --- Evaluate the order-v Eulerian invariants directly from the
          !     evolved conserved moments.  These are the same linear
          !     combinations used to construct the energy and momentum
          !     boundary flux tallies.  In particular, do not apply the
          !     nonlinear conserved-to-primitive inversion to an assembled
          !     pre-limiter state or to a synthetic stage-weighted state.

          Enu &
            = CR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) &
                + V_u_1 * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS) &
                + V_u_2 * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS) &
                + V_u_3 * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS)
          Fnu_d_1 &
            = CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS) &
                + V_d_1 * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)
          Fnu_d_2 &
            = CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS) &
                + V_d_2 * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)
          Fnu_d_3 &
            = CR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS) &
                + V_d_3 * CR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

          Radiation(iCheck_Energy) = Radiation(iCheck_Energy) &
            + SpaceWeight * PhaseWeight3 * Enu
          Radiation(iCheck_Momentum1) = Radiation(iCheck_Momentum1) &
            + SpaceWeight * PhaseWeight3 * Fnu_d_1
          Radiation(iCheck_Momentum2) = Radiation(iCheck_Momentum2) &
            + SpaceWeight * PhaseWeight3 * Fnu_d_2
          Radiation(iCheck_Momentum3) = Radiation(iCheck_Momentum3) &
            + SpaceWeight * PhaseWeight3 * Fnu_d_3

        END DO
        END DO
      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeConservedMatterRadiationComponents


  REAL(DP) FUNCTION ConvertCheckQuantity( iQuantity, Value ) RESULT( Converted )

    INTEGER, INTENT(in) :: iQuantity
    REAL(DP), INTENT(in) :: Value

    Converted = Value
    IF( iQuantity == iCheck_Energy ) Converted = Value / Erg

  END FUNCTION ConvertCheckQuantity


  SUBROUTINE InitializeMatterExchangeAttribution &
    ( FileName, iZ_B0, iZ_E0, ComputeShadows_Option )

    CHARACTER(*), INTENT(in) :: FileName
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    LOGICAL, INTENT(in), OPTIONAL :: ComputeShadows_Option

    INTEGER :: Unit

    IF( MatterExchangeAttributionInitialized ) &
      CALL FinalizeMatterExchangeAttribution

    IF( iZ_E0(3) /= iZ_B0(3) .OR. iZ_E0(4) /= iZ_B0(4) ) &
      ERROR STOP 'Matter-exchange attribution supports only 1D runs.'

    ALLOCATE &
      ( AttributionCumulative &
          (nAttrQuantities,nAttrComponents,nDOFX, &
           iZ_B0(2):iZ_E0(2)) )
    ALLOCATE( AttributionAbsolute, MOLD = AttributionCumulative )

    AttributionCumulative = Zero
    AttributionAbsolute   = Zero
    MatterExchangeAttributionCycle = 0
    MatterExchangeAttributionTime  = Zero
    MatterExchangeAttributionFile  = TRIM( FileName )
    MatterExchangeAttributionShadows = .TRUE.
    IF( PRESENT( ComputeShadows_Option ) ) &
      MatterExchangeAttributionShadows = ComputeShadows_Option

    OPEN( NEWUNIT = Unit, &
          FILE = TRIM(MatterExchangeAttributionFile) // '.nodes', &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# cycle time_s dt_s r_cm volume_weight_cm3 element node quantity' // &
      ' component step_increment step_rate cumulative_signed' // &
      ' cumulative_absolute'
    IF( MatterExchangeAttributionShadows )THEN
      WRITE(Unit,'(A)') '# attribution_mode full_shadows'
    ELSE
      WRITE(Unit,'(A)') '# attribution_mode actual_only'
    END IF
    CLOSE(Unit)

    OPEN( NEWUNIT = Unit, &
          FILE = TRIM(MatterExchangeAttributionFile) // '.elements', &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# cycle time_s dt_s r_cm volume_cm3 element quantity component' // &
      ' step_increment_mean step_abs_node_mean step_rate_mean' // &
      ' cumulative_signed_mean cumulative_absolute_node_mean'
    IF( MatterExchangeAttributionShadows )THEN
      WRITE(Unit,'(A)') '# attribution_mode full_shadows'
    ELSE
      WRITE(Unit,'(A)') '# attribution_mode actual_only'
    END IF
    CLOSE(Unit)

    MatterExchangeAttributionInitialized = .TRUE.

  END SUBROUTINE InitializeMatterExchangeAttribution


  SUBROUTINE FinalizeMatterExchangeAttribution

    IF( ALLOCATED( AttributionCumulative ) ) &
      DEALLOCATE( AttributionCumulative )
    IF( ALLOCATED( AttributionAbsolute ) ) &
      DEALLOCATE( AttributionAbsolute )

    MatterExchangeAttributionInitialized = .FALSE.
    MatterExchangeAttributionShadows = .TRUE.
    MatterExchangeAttributionCycle = 0
    MatterExchangeAttributionTime  = Zero
    MatterExchangeAttributionFile  = ''

  END SUBROUTINE FinalizeMatterExchangeAttribution


  SUBROUTINE AccumulateMatterExchangeAttribution &
    ( TimeAfterStep, TimeStep, DnuMax, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      CF_Start, CF_Actual, CF_NoStream, CF_AP, CF_Standard, &
      StageFluidState, StageRadiationState )

    REAL(DP), INTENT(in) :: TimeAfterStep, TimeStep, DnuMax
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_Start(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CF_Actual(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CF_NoStream(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CF_AP(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CF_Standard(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)

    INTEGER :: iX1, iNodeX
    REAL(DP) :: E_Start, E_Actual, E_NoStream, E_AP, E_Standard
    REAL(DP), ALLOCATABLE :: StepIncrement(:,:,:,:)
    REAL(DP), ALLOCATABLE :: ModalIncrement(:,:,:,:)

    IF( .NOT. MatterExchangeAttributionInitialized ) &
      ERROR STOP 'Matter-exchange attribution is not initialized.'
    IF( .NOT. MatterExchangeAttributionShadows ) &
      ERROR STOP 'Full matter-exchange attribution called in actual-only mode.'
    IF( TimeStep <= Zero ) &
      ERROR STOP 'Matter-exchange attribution requires a positive timestep.'

    ALLOCATE( StepIncrement, MOLD = AttributionCumulative )
    ALLOCATE &
      ( ModalIncrement &
          (nAttrQuantities,3,nDOFX,iZ_B0(2):iZ_E0(2)) )
    StepIncrement  = Zero
    ModalIncrement = Zero

    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      StepIncrement(iAttr_QL,iAttr_Actual,iNodeX,iX1) &
        = CF_Actual(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_Start(iNodeX,iX1,1,1,iCF_Ne)
      StepIncrement(iAttr_QL,iAttr_NoStream,iNodeX,iX1) &
        = CF_NoStream(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_Start(iNodeX,iX1,1,1,iCF_Ne)
      StepIncrement(iAttr_QL,iAttr_APMinusNoStream,iNodeX,iX1) &
        = CF_AP(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_NoStream(iNodeX,iX1,1,1,iCF_Ne)
      StepIncrement(iAttr_QL,iAttr_StandardMinusAP,iNodeX,iX1) &
        = CF_Standard(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_AP(iNodeX,iX1,1,1,iCF_Ne)
      StepIncrement(iAttr_QL,iAttr_StandardTotal,iNodeX,iX1) &
        = CF_Standard(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_Start(iNodeX,iX1,1,1,iCF_Ne)
      StepIncrement(iAttr_QL,iAttr_APTotal,iNodeX,iX1) &
        = CF_AP(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_Start(iNodeX,iX1,1,1,iCF_Ne)

      E_Start = InternalEnergyDensity &
        ( CF_Start(iNodeX,iX1,1,1,:), GX(iNodeX,iX1,1,1,:) )
      E_Actual = InternalEnergyDensity &
        ( CF_Actual(iNodeX,iX1,1,1,:), GX(iNodeX,iX1,1,1,:) )
      E_NoStream = InternalEnergyDensity &
        ( CF_NoStream(iNodeX,iX1,1,1,:), &
          GX(iNodeX,iX1,1,1,:) )
      E_AP = InternalEnergyDensity &
        ( CF_AP(iNodeX,iX1,1,1,:), GX(iNodeX,iX1,1,1,:) )
      E_Standard = InternalEnergyDensity &
        ( CF_Standard(iNodeX,iX1,1,1,:), GX(iNodeX,iX1,1,1,:) )

      StepIncrement(iAttr_QH,iAttr_Actual,iNodeX,iX1) &
        = E_Actual - E_Start
      StepIncrement(iAttr_QH,iAttr_NoStream,iNodeX,iX1) &
        = E_NoStream - E_Start
      StepIncrement(iAttr_QH,iAttr_APMinusNoStream,iNodeX,iX1) &
        = E_AP - E_NoStream
      StepIncrement(iAttr_QH,iAttr_StandardMinusAP,iNodeX,iX1) &
        = E_Standard - E_AP
      StepIncrement(iAttr_QH,iAttr_StandardTotal,iNodeX,iX1) &
        = E_Standard - E_Start
      StepIncrement(iAttr_QH,iAttr_APTotal,iNodeX,iX1) &
        = E_AP - E_Start

    END DO
    END DO

    CALL ComputeStageWeightedEmAbModalIncrement &
           ( TimeStep, DnuMax, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             StageFluidState, StageRadiationState, ModalIncrement )

    StepIncrement(:,iAttr_EmAbFull,:,:) = ModalIncrement(:,1,:,:)
    StepIncrement(:,iAttr_EmAbMean,:,:) = ModalIncrement(:,2,:,:)
    StepIncrement(:,iAttr_EmAbHigh,:,:) = ModalIncrement(:,3,:,:)

    AttributionCumulative = AttributionCumulative + StepIncrement
    AttributionAbsolute   = AttributionAbsolute + ABS( StepIncrement )
    MatterExchangeAttributionCycle &
      = MatterExchangeAttributionCycle + 1
    MatterExchangeAttributionTime = TimeAfterStep

    CALL WriteMatterExchangeAttributionStep &
           ( TimeAfterStep, TimeStep, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, StepIncrement )

    DEALLOCATE( StepIncrement, ModalIncrement )

  END SUBROUTINE AccumulateMatterExchangeAttribution


  SUBROUTINE AccumulateActualMatterExchangeAttribution &
    ( TimeAfterStep, TimeStep, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF_Start, CF_Actual )

    REAL(DP), INTENT(in) :: TimeAfterStep, TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_Start(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CF_Actual(1:nDOFX, &
        iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
        iZ_B1(4):iZ_E1(4),1:nCF)

    INTEGER :: iX1, iNodeX
    REAL(DP), ALLOCATABLE :: StepIncrement(:,:,:,:)

    IF( .NOT. MatterExchangeAttributionInitialized ) &
      ERROR STOP 'Matter-exchange attribution is not initialized.'
    IF( MatterExchangeAttributionShadows ) &
      ERROR STOP 'Actual-only attribution called with shadows enabled.'
    IF( TimeStep <= Zero ) &
      ERROR STOP 'Matter-exchange attribution requires a positive timestep.'

    ALLOCATE( StepIncrement, MOLD = AttributionCumulative )
    StepIncrement = Zero

    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      StepIncrement(iAttr_QL,iAttr_Actual,iNodeX,iX1) &
        = CF_Actual(iNodeX,iX1,1,1,iCF_Ne) &
          - CF_Start(iNodeX,iX1,1,1,iCF_Ne)
      StepIncrement(iAttr_QH,iAttr_Actual,iNodeX,iX1) &
        = InternalEnergyDensity &
            ( CF_Actual(iNodeX,iX1,1,1,:), GX(iNodeX,iX1,1,1,:) ) &
          - InternalEnergyDensity &
            ( CF_Start(iNodeX,iX1,1,1,:), GX(iNodeX,iX1,1,1,:) )

    END DO
    END DO

    AttributionCumulative(:,iAttr_Actual,:,:) &
      = AttributionCumulative(:,iAttr_Actual,:,:) &
        + StepIncrement(:,iAttr_Actual,:,:)
    AttributionAbsolute(:,iAttr_Actual,:,:) &
      = AttributionAbsolute(:,iAttr_Actual,:,:) &
        + ABS( StepIncrement(:,iAttr_Actual,:,:) )
    MatterExchangeAttributionCycle &
      = MatterExchangeAttributionCycle + 1
    MatterExchangeAttributionTime = TimeAfterStep

    CALL WriteMatterExchangeAttributionStep &
           ( TimeAfterStep, TimeStep, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, StepIncrement )

    DEALLOCATE( StepIncrement )

  END SUBROUTINE AccumulateActualMatterExchangeAttribution


  FUNCTION InternalEnergyDensity( CF, GF ) RESULT( Eint )

    REAL(DP), INTENT(in) :: CF(1:nCF), GF(1:nGF)
    REAL(DP) :: Eint

    Eint = CF(iCF_E) &
      - Half * ( CF(iCF_S1)**2 / GF(iGF_Gm_dd_11) &
                 + CF(iCF_S2)**2 / GF(iGF_Gm_dd_22) &
                 + CF(iCF_S3)**2 / GF(iGF_Gm_dd_33) ) / CF(iCF_D)

  END FUNCTION InternalEnergyDensity


  SUBROUTINE ComputeStageWeightedEmAbModalIncrement &
    ( TimeStep, DnuMax, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      StageFluidState, StageRadiationState, ModalIncrement )

    REAL(DP), INTENT(in) :: TimeStep, DnuMax
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(out) :: &
      ModalIncrement(nAttrQuantities,3,nDOFX,iZ_B0(2):iZ_E0(2))

    INTEGER :: iE, iN_E, iN_X, iNodeE, iNodeX, iNodeZ
    INTEGER :: iS, iStage, iX1, nE_G, nX_G
    REAL(DP) :: ElementVolume, ModeMean, StageWeight, VolumeWeight
    REAL(DP), ALLOCATABLE :: PF_Stage(:,:,:,:,:), AF_Stage(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR_Full(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR_Mean(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: E_N(:), D_N(:), T_N(:), Ye_N(:), Ym_N(:)
    REAL(DP), ALLOCATABLE :: DnuEq(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Full(:,:,:), Q_Mean(:,:,:)
    LOGICAL, ALLOCATABLE :: ElectronChargedCurrentMask(:)
    LOGICAL, ALLOCATABLE :: MuonChargedCurrentMask(:)

    nE_G = nDOFE * ( iZ_E0(1) - iZ_B0(1) + 1 )
    nX_G = nDOFX * ( iZ_E0(2) - iZ_B0(2) + 1 )

    ALLOCATE &
      ( PF_Stage &
          (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4),1:nPF) )
    ALLOCATE &
      ( AF_Stage &
          (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4),1:nAF) )
    ALLOCATE &
      ( PR_Full &
          (1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies) )
    ALLOCATE( PR_Mean, MOLD = PR_Full )
    ALLOCATE( E_N(nE_G), D_N(nX_G), T_N(nX_G), Ye_N(nX_G), Ym_N(nX_G) )
    ALLOCATE( DnuEq(nE_G,nSpecies,nX_G) )
    ALLOCATE( ElectronChargedCurrentMask(nX_G) )
    ALLOCATE( MuonChargedCurrentMask(nX_G) )
    ALLOCATE( Q_Full(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( Q_Mean, MOLD = Q_Full )

    ModalIncrement = Zero

    iN_E = 0
    DO iE = iZ_B0(1), iZ_E0(1)
    DO iNodeE = 1, nDOFE
      iN_E = iN_E + 1
      E_N(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )
    END DO
    END DO

    DO iStage = 1, SIZE( StageFluidState, 6 )

      StageWeight = GetImplicitStageWeight_IMEX_RK( iStage )
      IF( StageWeight == Zero ) CYCLE

      CALL BuildPrimitiveState &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               StageRadiationState &
                 (:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               PF_Stage, AF_Stage, PR_Full )

      iN_X = 0
      DO iX1 = iZ_B0(2), iZ_E0(2)
      DO iNodeX = 1, nDOFX
        iN_X = iN_X + 1
        D_N(iN_X) = PF_Stage(iNodeX,iX1,1,1,iPF_D)
        T_N(iN_X) = AF_Stage(iNodeX,iX1,1,1,iAF_T)
        Ye_N(iN_X) = AF_Stage(iNodeX,iX1,1,1,iAF_Ye)
        Ym_N(iN_X) = AF_Stage(iNodeX,iX1,1,1,iAF_Ym)
      END DO
      END DO

      CALL SetChargedCurrentMasks &
             ( nX_G, nDOFX, D_N, T_N, Ye_N, Ym_N, &
               ElectronChargedCurrentMask, MuonChargedCurrentMask )
      CALL ComputeEquilibriumDistributions &
             ( 1, nE_G, 1, nSpecies, 1, nX_G, &
               E_N, D_N, T_N, Ye_N, Ym_N, DnuEq, &
               ElectronChargedCurrentMask_Option &
                 = ElectronChargedCurrentMask, &
               MuonChargedCurrentMask_Option = MuonChargedCurrentMask )
      CALL LimitEquilibriumDistributions_DG &
             ( 1, nE_G, 1, nSpecies, 1, nX_G, &
               E_N, DnuEq, DnuMax )

      PR_Mean = PR_Full
      DO iS = 1, nSpecies
      DO iE = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE
        iN_E = ( iE - iZ_B0(1) ) * nDOFE + iNodeE
        DO iX1 = iZ_B0(2), iZ_E0(2)

          ElementVolume = Zero
          ModeMean      = Zero
          DO iNodeX = 1, nDOFX
            iN_X = ( iX1 - iZ_B0(2) ) * nDOFX + iNodeX
            iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
            VolumeWeight = WeightsX1(iNodeX) &
              * GX(iNodeX,iX1,1,1,iGF_h_1) &
              * GX(iNodeX,iX1,1,1,iGF_h_2)**2
            ElementVolume = ElementVolume + VolumeWeight
            ModeMean = ModeMean + VolumeWeight &
              * ( PR_Full(iNodeZ,iE,iX1,1,1,iPR_D,iS) &
                  - DnuEq(iN_E,iS,iN_X) )
          END DO
          ModeMean = ModeMean / ElementVolume

          DO iNodeX = 1, nDOFX
            iN_X = ( iX1 - iZ_B0(2) ) * nDOFX + iNodeX
            iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
            PR_Mean(iNodeZ,iE,iX1,1,1,iPR_D,iS) &
              = DnuEq(iN_E,iS,iN_X) + ModeMean
          END DO

        END DO
      END DO
      END DO
      END DO

      CALL ComputeMatterRHS_Collisions_OrderV &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               PF_Stage, PR_Full, Q_Full )
      CALL ComputeMatterRHS_Collisions_OrderV &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               PF_Stage, PR_Mean, Q_Mean )

      iN_X = 0
      DO iX1 = iZ_B0(2), iZ_E0(2)
      DO iNodeX = 1, nDOFX
        iN_X = iN_X + 1
        ModalIncrement(iAttr_QL,1,iNodeX,iX1) &
          = ModalIncrement(iAttr_QL,1,iNodeX,iX1) &
            + TimeStep * StageWeight &
              * Q_Full(iQ_Le,1+nSpecies+iC_EmAb,iN_X)
        ModalIncrement(iAttr_QL,2,iNodeX,iX1) &
          = ModalIncrement(iAttr_QL,2,iNodeX,iX1) &
            + TimeStep * StageWeight &
              * Q_Mean(iQ_Le,1+nSpecies+iC_EmAb,iN_X)
        ModalIncrement(iAttr_QH,1,iNodeX,iX1) &
          = ModalIncrement(iAttr_QH,1,iNodeX,iX1) &
            + TimeStep * StageWeight &
              * Q_Full(iQ_H,1+nSpecies+iC_EmAb,iN_X)
        ModalIncrement(iAttr_QH,2,iNodeX,iX1) &
          = ModalIncrement(iAttr_QH,2,iNodeX,iX1) &
            + TimeStep * StageWeight &
              * Q_Mean(iQ_H,1+nSpecies+iC_EmAb,iN_X)
      END DO
      END DO

    END DO

    ModalIncrement(:,3,:,:) &
      = ModalIncrement(:,1,:,:) - ModalIncrement(:,2,:,:)

    DEALLOCATE &
      ( PF_Stage, AF_Stage, PR_Full, PR_Mean, E_N, D_N, T_N, Ye_N, Ym_N, &
        DnuEq, Q_Full, Q_Mean, &
        ElectronChargedCurrentMask, MuonChargedCurrentMask )

  END SUBROUTINE ComputeStageWeightedEmAbModalIncrement


  SUBROUTINE WriteMatterExchangeAttributionStep &
    ( TimeAfterStep, TimeStep, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, StepIncrement )

    REAL(DP), INTENT(in) :: TimeAfterStep, TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      StepIncrement &
        (nAttrQuantities,nAttrComponents,nDOFX,iZ_B0(2):iZ_E0(2))

    INTEGER :: iC, iC_End, iQ, iX1, iNodeX, NodeUnit, ElementUnit
    REAL(DP) :: ElementVolume, Radius, UnitIncrement, UnitRate
    REAL(DP) :: VolumeWeight
    REAL(DP) :: StepMean, StepAbsMean, CumulativeMean, AbsoluteMean
    CHARACTER(2) :: QuantityName

    OPEN( NEWUNIT = NodeUnit, &
          FILE = TRIM(MatterExchangeAttributionFile) // '.nodes', &
          STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE', &
          FORM = 'FORMATTED' )
    OPEN( NEWUNIT = ElementUnit, &
          FILE = TRIM(MatterExchangeAttributionFile) // '.elements', &
          STATUS = 'OLD', POSITION = 'APPEND', ACTION = 'WRITE', &
          FORM = 'FORMATTED' )

    IF( MatterExchangeAttributionShadows )THEN
      iC_End = nAttrComponents
    ELSE
      iC_End = iAttr_Actual
    END IF

    DO iQ = 1, nAttrQuantities
      IF( iQ == iAttr_QL )THEN
        QuantityName = 'QL'
        UnitIncrement = One / Centimeter**3
        UnitRate      = One / Centimeter**3 / Second
      ELSE
        QuantityName = 'QH'
        UnitIncrement = Erg / Centimeter**3
        UnitRate      = Erg / Centimeter**3 / Second
      END IF

      DO iC = 1, iC_End
      DO iX1 = iZ_B0(2), iZ_E0(2)

        ElementVolume = Zero
        StepMean      = Zero
        StepAbsMean   = Zero
        CumulativeMean = Zero
        AbsoluteMean   = Zero

        DO iNodeX = 1, nDOFX
          Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )
          VolumeWeight = FourPi * WeightsX1(iNodeX) &
            * MeshX(1) % Width(iX1) &
            * GX(iNodeX,iX1,1,1,iGF_h_1) &
            * GX(iNodeX,iX1,1,1,iGF_h_2)**2
          ElementVolume = ElementVolume + VolumeWeight
          StepMean = StepMean &
            + VolumeWeight * StepIncrement(iQ,iC,iNodeX,iX1)
          StepAbsMean = StepAbsMean &
            + VolumeWeight * ABS( StepIncrement(iQ,iC,iNodeX,iX1) )
          CumulativeMean = CumulativeMean &
            + VolumeWeight * AttributionCumulative(iQ,iC,iNodeX,iX1)
          AbsoluteMean = AbsoluteMean &
            + VolumeWeight * AttributionAbsolute(iQ,iC,iNodeX,iX1)

          WRITE(NodeUnit, &
            '(I0,1X,4(ES24.16E3,1X),2(I0,1X),A,1X,A,1X,' // &
            '4(ES24.16E3,1X))') &
            MatterExchangeAttributionCycle, TimeAfterStep / Second, &
            TimeStep / Second, Radius / Centimeter, &
            VolumeWeight / Centimeter**3, &
            iX1 - iZ_B0(2) + 1, iNodeX, QuantityName, &
            TRIM(AttributionComponentName(iC)), &
            StepIncrement(iQ,iC,iNodeX,iX1) / UnitIncrement, &
            ( StepIncrement(iQ,iC,iNodeX,iX1) / TimeStep ) / UnitRate, &
            AttributionCumulative(iQ,iC,iNodeX,iX1) / UnitIncrement, &
            AttributionAbsolute(iQ,iC,iNodeX,iX1) / UnitIncrement
        END DO

        StepMean       = StepMean       / ElementVolume
        StepAbsMean    = StepAbsMean    / ElementVolume
        CumulativeMean = CumulativeMean / ElementVolume
        AbsoluteMean   = AbsoluteMean   / ElementVolume

        WRITE(ElementUnit, &
          '(I0,1X,4(ES24.16E3,1X),I0,1X,A,1X,A,1X,' // &
          '5(ES24.16E3,1X))') &
          MatterExchangeAttributionCycle, TimeAfterStep / Second, &
          TimeStep / Second, MeshX(1) % Center(iX1) / Centimeter, &
          ElementVolume / Centimeter**3, iX1 - iZ_B0(2) + 1, &
          QuantityName, TRIM(AttributionComponentName(iC)), &
          StepMean / UnitIncrement, StepAbsMean / UnitIncrement, &
          ( StepMean / TimeStep ) / UnitRate, &
          CumulativeMean / UnitIncrement, AbsoluteMean / UnitIncrement

      END DO
      END DO
    END DO

    CLOSE(NodeUnit)
    CLOSE(ElementUnit)

  END SUBROUTINE WriteMatterExchangeAttributionStep


  SUBROUTINE WriteStandaloneMatterDiagnostics &
    ( FileName, CoordinateSystem, Time, TimeStep, DnuMax, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
      GE, GX, CF_BeforeStep, CR_BeforeStep, CF_AfterStep, &
      PF, AF, CR, PR, GR, DR, &
      MatterCollisionExchange_Rate, MatterUpdate_Rate, &
      MatterCollisionRHS_FinalStage, RadiationState_BeforeFinalLimiter, &
      IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

    CHARACTER(*), INTENT(in) :: FileName
    CHARACTER(*), INTENT(in) :: CoordinateSystem
    REAL(DP),     INTENT(in) :: Time
    REAL(DP),     INTENT(in) :: TimeStep, DnuMax
    INTEGER,      INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_BeforeStep(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR_BeforeStep(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      CF_AfterStep(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(in) :: &
      AF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nAF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      PR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      GR(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
         1:nGR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      DR(iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nDR)
    REAL(DP), INTENT(in) :: &
      MatterCollisionExchange_Rate &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4))
    REAL(DP), INTENT(in) :: &
      MatterUpdate_Rate &
        (1:nQ,1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4))
    REAL(DP), INTENT(in) :: &
      MatterCollisionRHS_FinalStage &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4))
    REAL(DP), INTENT(in) :: &
      RadiationState_BeforeFinalLimiter &
        (1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: IMEX_Stage_FluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: IMEX_Stage_RadiationState(:,:,:,:,:,:,:,:,:)

    INTEGER :: iC, iNodeX, iN_X, iS, iX1, nE_G, nX_G, Unit
    INTEGER :: InnerIterationsMin_Element, InnerIterationsMax_Element
    INTEGER :: OuterIterationsMax_Element, InnerUpdateCount_Element
    REAL(DP) :: ElementVolume, RadialWeight, Radius, VolumeWeight
    REAL(DP) :: E_BeforeStep, E_AfterStep
    REAL(DP) :: E_Before_Element, E_After_Element
    REAL(DP) :: Ne_Before_Element, Ne_After_Element
    REAL(DP) :: Nm_Before_Element, Nm_After_Element
    REAL(DP) :: FractionalEChange_Element
    REAL(DP) :: QH_Applied_Element
    REAL(DP) :: FluxFactor, T_Flux, T_Spectral
    REAL(DP) :: D_Element, T_Fluid_Element, Ye_Element, Ym_Element
    REAL(DP) :: FirstInnerRelativeResidualMax_Element
    REAL(DP) :: D_Nu_Element(nSpecies), I_Nu_Element(nSpecies)
    REAL(DP) :: J_Nu_Element(nSpecies), H_Nu_Element(nSpecies)
    REAL(DP) :: MeanInvMFP_Element(nC,nSpecies)
    REAL(DP) :: EmAb_RateDecomposition_Element(nEmAbRateDecomposition,nSpecies)
    REAL(DP) :: NNS_RateDecomposition_Element(nNNSRateDecomposition,nSpecies)
    REAL(DP) :: Q_Matter_Element(nQ,1+nSpecies+nC)
    REAL(DP) :: Q_Matter_PreLimiter_Element(nQ,1+nSpecies+nC)
    REAL(DP) :: MatterCollisionExchange_Rate_Element(nQ,1+nSpecies+nC)
    REAL(DP) :: MatterCollisionRHS_FinalStage_Element(nQ,1+nSpecies+nC)
    REAL(DP) :: MatterUpdate_Rate_Element(nQ)
    REAL(DP) :: MatterEquationResidual_Rate_Final_Element(nQ)
    REAL(DP) :: ConservativeProjectionMatterHeating_Rate_Final_Element
    REAL(DP) :: QH_ByProcessSpecies_Element(nC,nSpecies)
    REAL(DP) :: Q_Transport_Element(nQ,1+nSpecies)
    REAL(DP), ALLOCATABLE :: &
      MeanInvMFP(:,:,:), EmAb_RateDecomposition(:,:,:), NNS_RateDecomposition(:,:,:), &
      Q_Matter(:,:,:), Q_Matter_PreLimiter(:,:,:), &
      QH_ByProcessSpecies(:,:,:), PR_BeforeFinalLimiter(:,:,:,:,:,:,:), &
      EmAb_Spectrum(:,:,:,:), EmAb_Spectrum_PreLimiter(:,:,:,:), &
      RadiationCollisionRHS(:,:,:,:), RadiationTransportRHS(:,:,:,:,:,:,:), &
      Q_Transport(:,:,:), CF_Streaming(:,:,:,:,:), CR_Streaming(:,:,:,:,:,:,:)
    INTEGER, ALLOCATABLE :: &
      InnerIterationsFinal(:), OuterIterationsFinal(:), &
      ConservativeProjectionStatusFinal(:)
    REAL(DP), ALLOCATABLE :: FirstInnerRelativeResidual(:)
    REAL(DP), ALLOCATABLE :: MatterEquationResidual_Rate_Final(:,:)
    REAL(DP), ALLOCATABLE :: InitialRadiationLimiter_Rate_Final(:,:)
    REAL(DP), ALLOCATABLE :: RadiationResidualRelativeFinal(:)
    REAL(DP), ALLOCATABLE :: &
      ConservativeProjectionRadiationResidualRelativeFinal(:), &
      ConservativeProjectionMatterHeatingRateFinal(:)
    REAL(DP), ALLOCATABLE :: ConservationHandoffCheckRateFinal(:,:)

    IF( iZ_E0(3) /= iZ_B0(3) .OR. iZ_E0(4) /= iZ_B0(4) ) THEN
      ERROR STOP 'WriteStandaloneMatterDiagnostics supports only 1D runs.'
    END IF

    IF( TRIM( CoordinateSystem ) /= 'SPHERICAL' ) THEN
      ERROR STOP 'WriteStandaloneMatterDiagnostics requires SPHERICAL coordinates.'
    END IF

    nE_G = nDOFE * ( iZ_E0(1) - iZ_B0(1) + 1 )
    nX_G = nDOFX * PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

    CALL GetNeutrinoMatterSolverConvergenceDiagnostics &
           ( InnerIterationsFinal, OuterIterationsFinal, &
             FirstInnerRelativeResidual, &
             MatterEquationResidual_Rate_Final, &
             InitialRadiationLimiter_Rate_Final, &
             RadiationResidualRelativeFinal, &
             ConservativeProjectionRadiationResidualRelativeFinal, &
             ConservativeProjectionMatterHeatingRateFinal, &
             ConservativeProjectionStatusFinal, &
             ConservationHandoffCheckRateFinal )

    ALLOCATE( MeanInvMFP(nC,nSpecies,nX_G) )
    ALLOCATE( EmAb_RateDecomposition(nEmAbRateDecomposition,nSpecies,nX_G) )
    ALLOCATE( NNS_RateDecomposition(nNNSRateDecomposition,nSpecies,nX_G) )
    ALLOCATE( Q_Matter(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( Q_Matter_PreLimiter(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( QH_ByProcessSpecies(nC,nSpecies,nX_G) )
    ALLOCATE( EmAb_Spectrum(nEmAbSpectrum,nE_G,nSpecies,nX_G) )
    ALLOCATE &
      ( EmAb_Spectrum_PreLimiter(nEmAbSpectrum,nE_G,nSpecies,nX_G) )
    ALLOCATE( RadiationCollisionRHS(4,nE_G,nSpecies,nX_G) )
    ALLOCATE( RadiationTransportRHS, MOLD = CR )
    ALLOCATE( Q_Transport(nQ,1+nSpecies,nX_G) )
    ALLOCATE( CF_Streaming, MOLD = CF_AfterStep )
    ALLOCATE( CR_Streaming, MOLD = CR )
    ALLOCATE( PR_BeforeFinalLimiter, MOLD = PR )

    PR_BeforeFinalLimiter = Zero

    CALL ComputePrimitiveRadiation &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, &
             RadiationState_BeforeFinalLimiter, PR_BeforeFinalLimiter )

    CALL ComputeMatterRHS_Collisions_OrderV &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF, &
             PR_BeforeFinalLimiter, Q_Matter_PreLimiter, &
             EmAb_Spectrum_Option = EmAb_Spectrum_PreLimiter )

    CALL ComputeMatterRHS_Collisions_OrderV &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF, PR, Q_Matter, &
             NNS_RateDecomposition_Option = NNS_RateDecomposition, &
             QH_ByProcessSpecies_Option = QH_ByProcessSpecies, &
             EmAb_RateDecomposition_Option = EmAb_RateDecomposition, &
             EmAb_Spectrum_Option = EmAb_Spectrum, &
             RadiationCollisionRHS_Option = RadiationCollisionRHS )

    CF_Streaming = CF_AfterStep
    CR_Streaming = CR
    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             CF_Streaming, CR_Streaming, RadiationTransportRHS )

    CALL ComputeTransportEquivalentMatterRHS &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF, &
             RadiationTransportRHS, Q_Transport )

    CALL ComputeMeanInverseMeanFreePaths_OrderV &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF, PR, MeanInvMFP )

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )

    WRITE(Unit,'(A,1X,ES24.16E3)') '# time_s', Time / Second
    WRITE(Unit,'(A,1X,ES24.16E3)') '# applied_dt_s', TimeStep / Second
    WRITE(Unit,'(A,1X,I0,1X,A,1X,I0)') '# nE', iZ_E0(1) - iZ_B0(1) + 1, &
      'nDOFE', nDOFE
    CALL WriteHeader( Unit )

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)

      ElementVolume        = Zero
      RadialWeight         = Zero
      D_Element            = Zero
      T_Fluid_Element      = Zero
      Ye_Element           = Zero
      Ym_Element           = Zero
      D_Nu_Element         = Zero
      I_Nu_Element         = Zero
      J_Nu_Element         = Zero
      H_Nu_Element         = Zero
      MeanInvMFP_Element   = Zero
      EmAb_RateDecomposition_Element  = Zero
      NNS_RateDecomposition_Element   = Zero
      Q_Matter_Element     = Zero
      Q_Matter_PreLimiter_Element = Zero
      MatterCollisionExchange_Rate_Element = Zero
      MatterCollisionRHS_FinalStage_Element = Zero
      MatterUpdate_Rate_Element = Zero
      MatterEquationResidual_Rate_Final_Element = Zero
      ConservativeProjectionMatterHeating_Rate_Final_Element = Zero
      QH_ByProcessSpecies_Element = Zero
      Q_Transport_Element = Zero
      QH_Applied_Element   = Zero
      E_Before_Element     = Zero
      E_After_Element      = Zero
      Ne_Before_Element    = Zero
      Ne_After_Element     = Zero
      Nm_Before_Element    = Zero
      Nm_After_Element     = Zero
      InnerIterationsMin_Element = HUGE(0)
      InnerIterationsMax_Element = 0
      OuterIterationsMax_Element = 0
      InnerUpdateCount_Element = 0
      FirstInnerRelativeResidualMax_Element = Zero

      DO iNodeX = 1, nDOFX

        iN_X = iN_X + 1
        Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )
        VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                       * GX(iNodeX,iX1,1,1,iGF_h_1) &
                       * GX(iNodeX,iX1,1,1,iGF_h_2)**2

        ElementVolume    = ElementVolume    + VolumeWeight
        RadialWeight = RadialWeight + WeightsX1(iNodeX)

        D_Element       = D_Element       + VolumeWeight * PF(iNodeX,iX1,1,1,iPF_D)
        T_Fluid_Element = T_Fluid_Element + VolumeWeight * AF(iNodeX,iX1,1,1,iAF_T)
        Ye_Element      = Ye_Element      + VolumeWeight * AF(iNodeX,iX1,1,1,iAF_Ye)
        Ym_Element      = Ym_Element      + VolumeWeight * AF(iNodeX,iX1,1,1,iAF_Ym)

        E_BeforeStep = CF_BeforeStep(iNodeX,iX1,1,1,iCF_E) &
          - Half * ( &
              CF_BeforeStep(iNodeX,iX1,1,1,iCF_S1)**2 &
                / GX(iNodeX,iX1,1,1,iGF_Gm_dd_11) &
            + CF_BeforeStep(iNodeX,iX1,1,1,iCF_S2)**2 &
                / GX(iNodeX,iX1,1,1,iGF_Gm_dd_22) &
            + CF_BeforeStep(iNodeX,iX1,1,1,iCF_S3)**2 &
                / GX(iNodeX,iX1,1,1,iGF_Gm_dd_33) ) &
            / CF_BeforeStep(iNodeX,iX1,1,1,iCF_D)

        E_AfterStep = CF_AfterStep(iNodeX,iX1,1,1,iCF_E) &
          - Half * ( &
              CF_AfterStep(iNodeX,iX1,1,1,iCF_S1)**2 &
                / GX(iNodeX,iX1,1,1,iGF_Gm_dd_11) &
            + CF_AfterStep(iNodeX,iX1,1,1,iCF_S2)**2 &
                / GX(iNodeX,iX1,1,1,iGF_Gm_dd_22) &
            + CF_AfterStep(iNodeX,iX1,1,1,iCF_S3)**2 &
                / GX(iNodeX,iX1,1,1,iGF_Gm_dd_33) ) &
            / CF_AfterStep(iNodeX,iX1,1,1,iCF_D)

        QH_Applied_Element = QH_Applied_Element &
          + VolumeWeight * ( E_AfterStep - E_BeforeStep ) / TimeStep
        E_Before_Element = E_Before_Element + VolumeWeight * E_BeforeStep
        E_After_Element = E_After_Element + VolumeWeight * E_AfterStep
        Ne_Before_Element = Ne_Before_Element &
          + VolumeWeight * CF_BeforeStep(iNodeX,iX1,1,1,iCF_Ne)
        Ne_After_Element = Ne_After_Element &
          + VolumeWeight * CF_AfterStep(iNodeX,iX1,1,1,iCF_Ne)
        Nm_Before_Element = Nm_Before_Element &
          + VolumeWeight * CF_BeforeStep(iNodeX,iX1,1,1,iCF_Nm)
        Nm_After_Element = Nm_After_Element &
          + VolumeWeight * CF_AfterStep(iNodeX,iX1,1,1,iCF_Nm)

        IF( InnerIterationsFinal(iN_X) > 0 )THEN
          InnerIterationsMin_Element &
            = MIN( InnerIterationsMin_Element, InnerIterationsFinal(iN_X) )
          InnerIterationsMax_Element &
            = MAX( InnerIterationsMax_Element, InnerIterationsFinal(iN_X) )
          IF( InnerIterationsFinal(iN_X) > 1 ) &
            InnerUpdateCount_Element = InnerUpdateCount_Element + 1
        END IF
        OuterIterationsMax_Element &
          = MAX( OuterIterationsMax_Element, OuterIterationsFinal(iN_X) )
        FirstInnerRelativeResidualMax_Element &
          = MAX( FirstInnerRelativeResidualMax_Element, &
                 FirstInnerRelativeResidual(iN_X) )

        DO iS = 1, nSpecies
          D_Nu_Element(iS) = D_Nu_Element(iS) &
                             + VolumeWeight * GR(iNodeX,iX1,1,1,iGR_D ,iS)
          I_Nu_Element(iS) = I_Nu_Element(iS) &
                             + VolumeWeight * GR(iNodeX,iX1,1,1,iGR_I1,iS)
          J_Nu_Element(iS) = J_Nu_Element(iS) &
                             + VolumeWeight * GR(iNodeX,iX1,1,1,iGR_J ,iS)
          H_Nu_Element(iS) = H_Nu_Element(iS) &
                             + VolumeWeight * GR(iNodeX,iX1,1,1,iGR_H1,iS)
        END DO

        Q_Matter_Element = Q_Matter_Element + VolumeWeight * Q_Matter(:,:,iN_X)
        Q_Matter_PreLimiter_Element = Q_Matter_PreLimiter_Element &
          + VolumeWeight * Q_Matter_PreLimiter(:,:,iN_X)
        MatterCollisionExchange_Rate_Element = MatterCollisionExchange_Rate_Element &
          + VolumeWeight * MatterCollisionExchange_Rate(:,:,iNodeX,iX1,1,1)
        MatterCollisionRHS_FinalStage_Element &
          = MatterCollisionRHS_FinalStage_Element &
            + VolumeWeight * MatterCollisionRHS_FinalStage(:,:,iNodeX,iX1,1,1)
        MatterUpdate_Rate_Element = MatterUpdate_Rate_Element &
          + VolumeWeight * MatterUpdate_Rate(:,iNodeX,iX1,1,1)
        MatterEquationResidual_Rate_Final_Element &
          = MatterEquationResidual_Rate_Final_Element &
            + VolumeWeight * MatterEquationResidual_Rate_Final(:,iN_X)
        ConservativeProjectionMatterHeating_Rate_Final_Element &
          = ConservativeProjectionMatterHeating_Rate_Final_Element &
            + VolumeWeight &
              * ConservativeProjectionMatterHeatingRateFinal(iN_X)
        EmAb_RateDecomposition_Element = EmAb_RateDecomposition_Element &
                              + VolumeWeight * EmAb_RateDecomposition(:,:,iN_X)
        NNS_RateDecomposition_Element = NNS_RateDecomposition_Element &
                             + VolumeWeight * NNS_RateDecomposition(:,:,iN_X)
        QH_ByProcessSpecies_Element = QH_ByProcessSpecies_Element &
                                      + VolumeWeight * QH_ByProcessSpecies(:,:,iN_X)
        Q_Transport_Element = Q_Transport_Element &
                              + VolumeWeight * Q_Transport(:,:,iN_X)

        DO iS = 1, nSpecies
        DO iC = 1, nC
          MeanInvMFP_Element(iC,iS) = MeanInvMFP_Element(iC,iS) &
                                       + WeightsX1(iNodeX) * MeanInvMFP(iC,iS,iN_X)
        END DO
        END DO

      END DO

      Radius                    = MeshX(1) % Center(iX1)
      D_Element                 = D_Element       / ElementVolume
      T_Fluid_Element           = T_Fluid_Element / ElementVolume
      Ye_Element                = Ye_Element      / ElementVolume
      Ym_Element                = Ym_Element      / ElementVolume
      D_Nu_Element              = D_Nu_Element    / ElementVolume
      I_Nu_Element              = I_Nu_Element    / ElementVolume
      J_Nu_Element              = J_Nu_Element    / ElementVolume
      H_Nu_Element              = H_Nu_Element    / ElementVolume
      Q_Matter_Element          = Q_Matter_Element / ElementVolume
      Q_Matter_PreLimiter_Element &
        = Q_Matter_PreLimiter_Element / ElementVolume
      MatterCollisionExchange_Rate_Element = MatterCollisionExchange_Rate_Element / ElementVolume
      MatterCollisionRHS_FinalStage_Element &
        = MatterCollisionRHS_FinalStage_Element / ElementVolume
      MatterUpdate_Rate_Element &
        = MatterUpdate_Rate_Element / ElementVolume
      MatterEquationResidual_Rate_Final_Element &
        = MatterEquationResidual_Rate_Final_Element / ElementVolume
      ConservativeProjectionMatterHeating_Rate_Final_Element &
        = ConservativeProjectionMatterHeating_Rate_Final_Element / ElementVolume
      EmAb_RateDecomposition_Element       = EmAb_RateDecomposition_Element / ElementVolume
      NNS_RateDecomposition_Element        = NNS_RateDecomposition_Element / ElementVolume
      QH_ByProcessSpecies_Element = QH_ByProcessSpecies_Element / ElementVolume
      Q_Transport_Element          = Q_Transport_Element / ElementVolume
      QH_Applied_Element          = QH_Applied_Element / ElementVolume
      E_Before_Element            = E_Before_Element / ElementVolume
      E_After_Element             = E_After_Element / ElementVolume
      Ne_Before_Element           = Ne_Before_Element / ElementVolume
      Ne_After_Element            = Ne_After_Element / ElementVolume
      Nm_Before_Element           = Nm_Before_Element / ElementVolume
      Nm_After_Element            = Nm_After_Element / ElementVolume
      FractionalEChange_Element   = Zero
      IF( InnerIterationsMin_Element == HUGE(0) ) &
        InnerIterationsMin_Element = 0
      IF( E_After_Element /= Zero )THEN
        FractionalEChange_Element &
          = TimeStep * QH_Applied_Element / E_After_Element
      END IF
      ! Sources and moments are spatial-volume averages.  The inverse MFPs
      ! instead average over dr, since they are subsequently integrated as
      ! inverse-MFP times dr to form an optical depth.
      MeanInvMFP_Element        = MeanInvMFP_Element / RadialWeight

      T_Fluid_Element = BoltzmannConstant * T_Fluid_Element / MeV

      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Radius / Centimeter
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MeshX(1) % Width(iX1) / Centimeter
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ElementVolume / Centimeter**3
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') D_Element / ( Gram / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') T_Fluid_Element
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') Ye_Element
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') Ym_Element
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Ne_Before_Element / ( One / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Ne_After_Element / ( One / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Nm_Before_Element / ( One / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Nm_After_Element / ( One / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        E_Before_Element / ( Erg / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        E_After_Element / ( Erg / Centimeter**3 )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') FractionalEChange_Element
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') DR(iX1,1,1,iDR_PL_Theta_1)
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') DR(iX1,1,1,iDR_PL_Theta_2)
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        DR(iX1,1,1,iDR_PL_dEnergy) / Erg
      WRITE(Unit,'(I0,1X)',ADVANCE='NO') OuterIterationsMax_Element
      WRITE(Unit,'(I0,1X)',ADVANCE='NO') InnerIterationsMin_Element
      WRITE(Unit,'(I0,1X)',ADVANCE='NO') InnerIterationsMax_Element
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        FirstInnerRelativeResidualMax_Element
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        REAL( InnerUpdateCount_Element, DP ) / REAL( nDOFX, DP )

      DO iS = 1, nSpecies

        T_Spectral = Zero
        IF( D_Nu_Element(iS) > Zero ) THEN
          T_Spectral = J_Nu_Element(iS) &
                       / ( Three * D_Nu_Element(iS) * MeV )
        END IF

        T_Flux = Zero
        IF( I_Nu_Element(iS) /= Zero ) THEN
          T_Flux = H_Nu_Element(iS) &
                 / ( Three * I_Nu_Element(iS) * MeV )
        END IF

        FluxFactor = Zero
        IF( D_Nu_Element(iS) /= Zero ) THEN
          FluxFactor = I_Nu_Element(iS) / D_Nu_Element(iS)
        END IF

        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') T_Spectral
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') T_Flux
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') FluxFactor

      END DO

      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ( Ne_After_Element - Ne_Before_Element ) / TimeStep &
          / ( One / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ( Nm_After_Element - Nm_Before_Element ) / TimeStep &
          / ( One / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterCollisionExchange_Rate_Element(iQ_Le,1) &
          / ( One / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterCollisionExchange_Rate_Element(iQ_Lm,1) &
          / ( One / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterUpdate_Rate_Element(iQ_Le) &
          / ( One / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterUpdate_Rate_Element(iQ_Lm) &
          / ( One / Centimeter**3 / Second )

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MatterCollisionExchange_Rate_Element(iQ_Le,1+iS) &
            / ( One / Centimeter**3 / Second )
      END DO

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MatterCollisionExchange_Rate_Element(iQ_Lm,1+iS) &
            / ( One / Centimeter**3 / Second )
      END DO

      DO iC = 1, nC
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MatterCollisionExchange_Rate_Element(iQ_Le,1+nSpecies+iC) &
            / ( One / Centimeter**3 / Second )
      END DO

      DO iC = 1, nC
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MatterCollisionExchange_Rate_Element(iQ_Lm,1+nSpecies+iC) &
            / ( One / Centimeter**3 / Second )
      END DO

      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Q_Matter_Element(iQ_H,1) / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Q_Transport_Element(iQ_H,1) / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ( Q_Transport_Element(iQ_H,1) - Q_Matter_Element(iQ_H,1) ) &
          / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Q_Matter_PreLimiter_Element(iQ_H,1) &
          / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterCollisionRHS_FinalStage_Element(iQ_H,1) &
          / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        QH_Applied_Element / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterCollisionExchange_Rate_Element(iQ_H,1) / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterUpdate_Rate_Element(iQ_H) / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        MatterEquationResidual_Rate_Final_Element(iQ_H) &
          / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ConservativeProjectionMatterHeating_Rate_Final_Element &
          / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ( MatterCollisionRHS_FinalStage_Element(iQ_H,1) &
            - MatterUpdate_Rate_Element(iQ_H) ) &
          / ( Erg / Centimeter**3 / Second )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        ( MatterEquationResidual_Rate_Final_Element(iQ_H) &
            - MatterCollisionRHS_FinalStage_Element(iQ_H,1) &
            + MatterUpdate_Rate_Element(iQ_H) ) &
          / ( Erg / Centimeter**3 / Second )

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MatterCollisionExchange_Rate_Element(iQ_H,1+iS) &
            / ( Erg / Centimeter**3 / Second )
      END DO

      DO iC = 1, nC
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MatterCollisionExchange_Rate_Element(iQ_H,1+nSpecies+iC) &
            / ( Erg / Centimeter**3 / Second )
      END DO

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Q_Matter_Element(iQ_H,1+iS) / ( Erg / Centimeter**3 / Second )
      END DO

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Q_Transport_Element(iQ_H,1+iS) / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          ( Q_Transport_Element(iQ_H,1+iS) &
              - Q_Matter_Element(iQ_H,1+iS) ) &
            / ( Erg / Centimeter**3 / Second )
      END DO

      DO iC = 1, nC
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Q_Matter_Element(iQ_H,1+nSpecies+iC) / ( Erg / Centimeter**3 / Second )
      END DO

      DO iS = 1, nSpecies
      DO iC = 1, nC
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          QH_ByProcessSpecies_Element(iC,iS) / ( Erg / Centimeter**3 / Second )
      END DO
      END DO

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          EmAb_RateDecomposition_Element(iEmAbRateDecomposition_NIn,iS) / ( One / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          EmAb_RateDecomposition_Element(iEmAbRateDecomposition_NOut,iS) / ( One / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          EmAb_RateDecomposition_Element(iEmAbRateDecomposition_NNet,iS) / ( One / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          EmAb_RateDecomposition_Element(iEmAbRateDecomposition_EIn,iS) / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          EmAb_RateDecomposition_Element(iEmAbRateDecomposition_EOut,iS) / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          EmAb_RateDecomposition_Element(iEmAbRateDecomposition_QH,iS) / ( Erg / Centimeter**3 / Second )
      END DO

      DO iS = 1, nSpecies
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          NNS_RateDecomposition_Element(iNNSRateDecomposition_NIn,iS) / ( One / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          NNS_RateDecomposition_Element(iNNSRateDecomposition_NOut,iS) / ( One / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          NNS_RateDecomposition_Element(iNNSRateDecomposition_NNet,iS) / ( One / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          NNS_RateDecomposition_Element(iNNSRateDecomposition_EIn,iS) / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          NNS_RateDecomposition_Element(iNNSRateDecomposition_EOut,iS) / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          NNS_RateDecomposition_Element(iNNSRateDecomposition_QH,iS) / ( Erg / Centimeter**3 / Second )
      END DO

      DO iS = 1, nSpecies
      DO iC = 1, nC
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          MeanInvMFP_Element(iC,iS) * Centimeter
      END DO
      END DO

      WRITE(Unit,*)

    END DO

    CLOSE(Unit)

    CALL WriteEmAbSpectrum &
           ( TRIM( FileName ) // '.emab_spectrum', &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
             EmAb_Spectrum_PreLimiter, EmAb_Spectrum )

    CALL WriteTransportCollisionCheck &
           ( TRIM( FileName ) // '.transport_collision', &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, &
             Q_Matter, Q_Transport, RadiationCollisionRHS, &
             RadiationTransportRHS )

    CALL WriteIMEXStageCheck &
           ( TRIM( FileName ) // '.imex_stages', &
             TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

    CALL WriteIMEXStreamingCheck &
           ( TRIM( FileName ) // '.imex_streaming', &
             TRIM( FileName ) // '.imex_faces', &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             MeanInvMFP, DnuMax, &
             IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

    CALL WriteDiscreteStepBudget &
           ( TRIM( FileName ) // '.discrete_budget', TimeStep, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             CF_BeforeStep, CR_BeforeStep, CF_AfterStep, CR, &
             IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

    CALL WriteSolverConvergenceCheck &
           ( TRIM( FileName ) // '.solver_convergence', &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, InnerIterationsFinal, &
             OuterIterationsFinal, FirstInnerRelativeResidual, &
             MatterEquationResidual_Rate_Final, &
             InitialRadiationLimiter_Rate_Final, &
             RadiationResidualRelativeFinal, &
             ConservativeProjectionRadiationResidualRelativeFinal, &
             ConservativeProjectionMatterHeatingRateFinal, &
             ConservativeProjectionStatusFinal, &
             MatterCollisionRHS_FinalStage, MatterUpdate_Rate )

    CALL WriteSolverConservationCheck &
           ( TRIM( FileName ) // '.solver_conservation_check', TimeStep, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             IMEX_Stage_FluidState, IMEX_Stage_RadiationState, &
             MatterEquationResidual_Rate_Final, &
             InitialRadiationLimiter_Rate_Final, &
             ConservationHandoffCheckRateFinal )

    CALL WriteSolverConservationHandoffCheck &
           ( TRIM( FileName ) // '.solver_conservation_handoff', &
             iZ_B0, iZ_E0, ConservationHandoffCheckRateFinal, &
             InitialRadiationLimiter_Rate_Final, &
             MatterEquationResidual_Rate_Final )

    CALL WriteFrozenTransientOperator &
           ( TRIM( FileName ) // '.frozen_transient', &
             TRIM( FileName ) // '.frozen_transient_stages', &
             Time - TimeStep, TimeStep, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             CF_BeforeStep, CR_BeforeStep, CR, &
             IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

    DEALLOCATE &
      ( MeanInvMFP, EmAb_RateDecomposition, NNS_RateDecomposition, Q_Matter, Q_Matter_PreLimiter, &
        QH_ByProcessSpecies, PR_BeforeFinalLimiter, &
        EmAb_Spectrum, EmAb_Spectrum_PreLimiter, &
        RadiationCollisionRHS, RadiationTransportRHS, Q_Transport, &
        CF_Streaming, CR_Streaming, &
        InnerIterationsFinal, OuterIterationsFinal, &
        FirstInnerRelativeResidual, MatterEquationResidual_Rate_Final, &
        InitialRadiationLimiter_Rate_Final, &
        RadiationResidualRelativeFinal, &
        ConservativeProjectionRadiationResidualRelativeFinal, &
        ConservativeProjectionMatterHeatingRateFinal, &
        ConservativeProjectionStatusFinal, ConservationHandoffCheckRateFinal )

  END SUBROUTINE WriteStandaloneMatterDiagnostics


  SUBROUTINE WriteFrozenTransientOperator &
    ( FileName, StageFileName, StateTime, TimeStep, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      CF_BeforeStep, CR_BeforeStep, CR_AfterStep, &
      StageFluidState, StageRadiationState )

    CHARACTER(*), INTENT(in) :: FileName, StageFileName
    REAL(DP), INTENT(in) :: StateTime, TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_BeforeStep(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR_BeforeStep(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      CR_AfterStep(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)

    INTEGER, PARAMETER :: iRaw_Total    = 1
    INTEGER, PARAMETER :: iRaw_NoLLF    = 2
    INTEGER, PARAMETER :: iRaw_Centered = 3
    INTEGER, PARAMETER :: iRaw_Surface  = 4
    INTEGER, PARAMETER :: iRaw_Volume   = 5
    INTEGER, PARAMETER :: nRaw          = 5
    REAL(DP), PARAMETER :: LLFFactor(nRaw) &
      = [ One, Zero, Zero, Zero, Zero ]
    REAL(DP), PARAMETER :: SurfaceFactor(nRaw) &
      = [ One, One, One, One, Zero ]
    REAL(DP), PARAMETER :: VolumeFactor(nRaw) &
      = [ One, One, One, Zero, One ]
    LOGICAL, PARAMETER :: IncludeObserver(nRaw) &
      = [ .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. ]

    INTEGER :: iCR, iE, iN_E, iN_X, iNodeE, iNodeX, iNodeZ
    INTEGER :: iRaw, iS, iStage, iState, iStateX1, iStateX2, iStateX3
    INTEGER :: iX1, nE_G, nX_G, StageUnit, Unit
    REAL(DP) :: ActualRate, Energy, EulerPrediction, Radius
    REAL(DP) :: SemidiscreteRate
    REAL(DP), ALLOCATABLE :: AF_Initial(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF_Work(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR_Work(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: EmAbSpectrum(:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF_Initial(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR_Initial(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Dummy(:,:,:)
    REAL(DP), ALLOCATABLE :: RadiationCollisionRHS(:,:,:,:)
    REAL(DP), ALLOCATABLE :: StreamingRHS(:,:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: TransportOpacity(:,:,:)

    nE_G = nDOFE * ( iZ_E0(1) - iZ_B0(1) + 1 )
    nX_G = nDOFX * PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

    ALLOCATE &
      ( PF_Initial(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nPF) )
    ALLOCATE &
      ( AF_Initial(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nAF) )
    ALLOCATE &
      ( PR_Initial(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies) )
    ALLOCATE( CF_Work, MOLD = CF_BeforeStep )
    ALLOCATE( CR_Work, MOLD = CR_BeforeStep )
    ALLOCATE( Q_Dummy(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( EmAbSpectrum(nEmAbSpectrum,nE_G,nSpecies,nX_G) )
    ALLOCATE( RadiationCollisionRHS(nCR,nE_G,nSpecies,nX_G) )
    ALLOCATE( TransportOpacity(nE_G,nSpecies,nX_G) )
    ALLOCATE &
      ( StreamingRHS(1:nDOFZ,iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies,1:nRaw) )

    CALL BuildPrimitiveState &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
             CF_BeforeStep, CR_BeforeStep, PF_Initial, AF_Initial, PR_Initial )

    CALL ComputeMatterRHS_Collisions_OrderV &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF_Initial, PR_Initial, &
             Q_Dummy, EmAb_Spectrum_Option = EmAbSpectrum, &
             RadiationCollisionRHS_Option = RadiationCollisionRHS )

    DO iRaw = 1, nRaw
      CF_Work = CF_BeforeStep
      CR_Work = CR_BeforeStep
      IF( iRaw == iRaw_Total )THEN
        CALL ComputeIncrement_TwoMoment_Explicit_AP &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF_Work, CR_Work, &
                 StreamingRHS(:,:,:,:,:,:,:,iRaw), &
                 LLFPenaltyFactor_Option = LLFFactor(iRaw), &
                 SurfaceTermFactor_Option = SurfaceFactor(iRaw), &
                 VolumeTermFactor_Option = VolumeFactor(iRaw), &
                 IncludeObserverCorrections_Option = IncludeObserver(iRaw) )
      ELSE
        CALL ComputeIncrement_TwoMoment_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF_Work, CR_Work, &
                 StreamingRHS(:,:,:,:,:,:,:,iRaw), &
                 LLFPenaltyFactor_Option = LLFFactor(iRaw), &
                 SurfaceTermFactor_Option = SurfaceFactor(iRaw), &
                 VolumeTermFactor_Option = VolumeFactor(iRaw), &
                 IncludeObserverCorrections_Option = IncludeObserver(iRaw) )
      END IF
    END DO

    CF_Work = CF_BeforeStep
    CR_Work = CR_BeforeStep
    CALL ComputeTransportOpacityFromConserved_OrderV &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, CF_Work, CR_Work, &
             TransportOpacity_Option = TransportOpacity )

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A,1X,ES24.16E3)') '# state_time_code', StateTime
    WRITE(Unit,'(A,1X,ES24.16E3)') '# state_time_s', StateTime / Second
    WRITE(Unit,'(A,1X,ES24.16E3)') '# timestep_code', TimeStep
    WRITE(Unit,'(A,1X,ES24.16E3)') '# timestep_s', TimeStep / Second
    WRITE(Unit,'(A,1X,ES24.16E3)') '# centimeter_code', Centimeter
    WRITE(Unit,'(A,1X,ES24.16E3)') '# second_code', Second
    WRITE(Unit,'(A,1X,I0)') '# n_spatial_nodes', nDOFX
    WRITE(Unit,'(A,1X,I0)') '# n_energy_nodes', nDOFE
    WRITE(Unit,'(A,1X,I0)') '# source_aware_llf', &
      MERGE( 1, 0, UseSourceAwareLLF )
    WRITE(Unit,'(A)') &
      '# element spatial_node energy_element energy_node species' // &
      ' flat_spatial flat_energy r_code dr_code E_code dE_code' // &
      ' wx wE h1 h2 h3 sqrtgm gm11 gm22 gm33' // &
      ' rho_code T_code Ye v1 v2 v3' // &
      ' N_before G1_before G2_before G3_before' // &
      ' N_after G1_after G2_after G3_after' // &
      ' D_before I1_before I2_before I3_before' // &
      ' D0 chi_emab kappa_transport kappa_minus_chi_emab' // &
      ' coll_N coll_G1 coll_G2 coll_G3' // &
      ' stream_total_N stream_total_G1 stream_total_G2 stream_total_G3' // &
      ' stream_noLLF_N stream_noLLF_G1 stream_noLLF_G2 stream_noLLF_G3' // &
      ' stream_centered_N stream_centered_G1 stream_centered_G2 stream_centered_G3' // &
      ' stream_surface_N stream_surface_G1 stream_surface_G2 stream_surface_G3' // &
      ' stream_volume_N stream_volume_G1 stream_volume_G2 stream_volume_G3' // &
      ' actual_rate_N actual_rate_G1 actual_rate_G2 actual_rate_G3' // &
      ' semidiscrete_N semidiscrete_G1 semidiscrete_G2 semidiscrete_G3' // &
      ' forward_euler_N forward_euler_G1 forward_euler_G2 forward_euler_G3'

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
      iN_X = iN_X + 1
      Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )

      DO iS = 1, nSpecies
      DO iE = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE
        iN_E = iNodeE + ( iE - iZ_B0(1) ) * nDOFE
        iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
        Energy = NodeCoordinate( MeshE, iE, iNodeE )

        WRITE(Unit,'(7(I0,1X))',ADVANCE='NO') &
          iX1-iZ_B0(2)+1, iNodeX, iE-iZ_B0(1)+1, iNodeE, iS, iN_X, iN_E
        WRITE(Unit,'(*(ES24.16E3,1X))',ADVANCE='NO') &
          Radius, MeshX(1)%Width(iX1), Energy, MeshE%Width(iE), &
          WeightsX1(iNodeX), WeightsE(iNodeE), &
          GX(iNodeX,iX1,1,1,iGF_h_1), GX(iNodeX,iX1,1,1,iGF_h_2), &
          GX(iNodeX,iX1,1,1,iGF_h_3), GX(iNodeX,iX1,1,1,iGF_SqrtGm), &
          GX(iNodeX,iX1,1,1,iGF_Gm_dd_11), &
          GX(iNodeX,iX1,1,1,iGF_Gm_dd_22), &
          GX(iNodeX,iX1,1,1,iGF_Gm_dd_33), &
          PF_Initial(iNodeX,iX1,1,1,iPF_D), &
          AF_Initial(iNodeX,iX1,1,1,iAF_T), &
          AF_Initial(iNodeX,iX1,1,1,iAF_Ye), &
          PF_Initial(iNodeX,iX1,1,1,iPF_V1), &
          PF_Initial(iNodeX,iX1,1,1,iPF_V2), &
          PF_Initial(iNodeX,iX1,1,1,iPF_V3)
        DO iCR = 1, nCR
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            CR_BeforeStep(iNodeZ,iE,iX1,1,1,iCR,iS)
        END DO
        DO iCR = 1, nCR
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            CR_AfterStep(iNodeZ,iE,iX1,1,1,iCR,iS)
        END DO
        DO iCR = 1, nPR
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            PR_Initial(iNodeZ,iE,iX1,1,1,iCR,iS)
        END DO
        WRITE(Unit,'(*(ES24.16E3,1X))',ADVANCE='NO') &
          EmAbSpectrum(iEmAbSpectrum_DnuEq,iN_E,iS,iN_X), &
          EmAbSpectrum(iEmAbSpectrum_Chi,iN_E,iS,iN_X), &
          TransportOpacity(iN_E,iS,iN_X), &
          TransportOpacity(iN_E,iS,iN_X) &
            - EmAbSpectrum(iEmAbSpectrum_Chi,iN_E,iS,iN_X)
        DO iCR = 1, nCR
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            RadiationCollisionRHS(iCR,iN_E,iS,iN_X)
        END DO
        DO iRaw = 1, nRaw
        DO iCR = 1, nCR
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            StreamingRHS(iNodeZ,iE,iX1,1,1,iCR,iS,iRaw)
        END DO
        END DO
        DO iCR = 1, nCR
          ActualRate = &
            ( CR_AfterStep(iNodeZ,iE,iX1,1,1,iCR,iS) &
              - CR_BeforeStep(iNodeZ,iE,iX1,1,1,iCR,iS) ) / TimeStep
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') ActualRate
        END DO
        DO iCR = 1, nCR
          SemidiscreteRate = &
            StreamingRHS(iNodeZ,iE,iX1,1,1,iCR,iS,iRaw_Total) &
              + RadiationCollisionRHS(iCR,iN_E,iS,iN_X)
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') SemidiscreteRate
        END DO
        DO iCR = 1, nCR
          SemidiscreteRate = &
            StreamingRHS(iNodeZ,iE,iX1,1,1,iCR,iS,iRaw_Total) &
              + RadiationCollisionRHS(iCR,iN_E,iS,iN_X)
          EulerPrediction = CR_BeforeStep(iNodeZ,iE,iX1,1,1,iCR,iS) &
                              + TimeStep * SemidiscreteRate
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') EulerPrediction
        END DO
        WRITE(Unit,*)
      END DO
      END DO
      END DO
    END DO
    END DO
    CLOSE(Unit)

    OPEN( NEWUNIT = StageUnit, FILE = TRIM( StageFileName ), &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(StageUnit,'(A,1X,ES24.16E3)') '# timestep_code', TimeStep
    DO iStage = 1, SIZE( StageFluidState, 6 )
      WRITE(StageUnit,'(A,I0,A,3(ES24.16E3,1X))') &
        '# stage ', iStage, ' aii_wEX_wIM ', &
        GetImplicitStageCoefficient_IMEX_RK(iStage), &
        GetExplicitStageWeight_IMEX_RK(iStage), &
        GetImplicitStageWeight_IMEX_RK(iStage)
      WRITE(StageUnit,'(A,I0,A)',ADVANCE='NO') '# stage ', iStage, ' aEX '
      DO iState = 1, SIZE( StageFluidState, 6 )
        WRITE(StageUnit,'(ES24.16E3,1X)',ADVANCE='NO') &
          GetExplicitStageCoefficient_IMEX_RK(iStage,iState)
      END DO
      WRITE(StageUnit,*)
      WRITE(StageUnit,'(A,I0,A)',ADVANCE='NO') '# stage ', iStage, ' aIM '
      DO iState = 1, SIZE( StageFluidState, 6 )
        WRITE(StageUnit,'(ES24.16E3,1X)',ADVANCE='NO') &
          GetImplicitStageCoefficient_IMEX_RK(iStage,iState)
      END DO
      WRITE(StageUnit,*)
    END DO
    WRITE(StageUnit,'(A)') &
      '# stage state element spatial_node energy_element energy_node species' // &
      ' N G1 G2 G3 fluid_D fluid_S1 fluid_S2 fluid_S3 fluid_E fluid_Ne'
    DO iStage = 1, SIZE( StageFluidState, 6 )
    DO iState = 1, SIZE( StageFluidState, 7 )
    DO iX1 = iZ_B0(2), iZ_E0(2)
      iStateX1 = iX1 - iZ_B1(2) + 1
      iStateX2 = 1 - iZ_B1(3) + 1
      iStateX3 = 1 - iZ_B1(4) + 1
      DO iNodeX = 1, nDOFX
      DO iS = 1, nSpecies
      DO iE = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE
        iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE
        WRITE(StageUnit,'(7(I0,1X))',ADVANCE='NO') &
          iStage, iState, iX1-iZ_B0(2)+1, iNodeX, &
          iE-iZ_B0(1)+1, iNodeE, iS
        DO iCR = 1, nCR
          WRITE(StageUnit,'(ES24.16E3,1X)',ADVANCE='NO') &
            StageRadiationState &
              (iNodeZ,iE,iStateX1,iStateX2,iStateX3,iCR,iS,iStage,iState)
        END DO
        WRITE(StageUnit,'(*(ES24.16E3,1X))') &
          StageFluidState(iNodeX,iStateX1,iStateX2,iStateX3,iCF_D, &
                          iStage,iState), &
          StageFluidState(iNodeX,iStateX1,iStateX2,iStateX3,iCF_S1, &
                          iStage,iState), &
          StageFluidState(iNodeX,iStateX1,iStateX2,iStateX3,iCF_S2, &
                          iStage,iState), &
          StageFluidState(iNodeX,iStateX1,iStateX2,iStateX3,iCF_S3, &
                          iStage,iState), &
          StageFluidState(iNodeX,iStateX1,iStateX2,iStateX3,iCF_E, &
                          iStage,iState), &
          StageFluidState(iNodeX,iStateX1,iStateX2,iStateX3,iCF_Ne, &
                          iStage,iState)
      END DO
      END DO
      END DO
      END DO
    END DO
    END DO
    END DO
    CLOSE(StageUnit)

    DEALLOCATE &
      ( PF_Initial, AF_Initial, PR_Initial, CF_Work, CR_Work, Q_Dummy, &
        EmAbSpectrum, RadiationCollisionRHS, StreamingRHS, TransportOpacity )

  END SUBROUTINE WriteFrozenTransientOperator


  SUBROUTINE WriteSolverConvergenceCheck &
    ( FileName, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, InnerIterationsFinal, &
      OuterIterationsFinal, FirstInnerRelativeResidual, &
      MatterEquationResidual_Rate_Final, &
      InitialRadiationLimiter_Rate_Final, &
      RadiationResidualRelativeFinal, &
      ConservativeProjectionRadiationResidualRelativeFinal, &
      ConservativeProjectionMatterHeatingRateFinal, &
      ConservativeProjectionStatusFinal, &
      MatterCollisionRHS_FinalStage, MatterUpdate_Rate )

    CHARACTER(*), INTENT(in) :: FileName
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    INTEGER, INTENT(in) :: InnerIterationsFinal(:)
    INTEGER, INTENT(in) :: OuterIterationsFinal(:)
    REAL(DP), INTENT(in) :: FirstInnerRelativeResidual(:)
    REAL(DP), INTENT(in) :: MatterEquationResidual_Rate_Final(:,:)
    REAL(DP), INTENT(in) :: InitialRadiationLimiter_Rate_Final(:,:)
    REAL(DP), INTENT(in) :: RadiationResidualRelativeFinal(:)
    REAL(DP), INTENT(in) :: &
      ConservativeProjectionRadiationResidualRelativeFinal(:)
    REAL(DP), INTENT(in) :: &
      ConservativeProjectionMatterHeatingRateFinal(:)
    INTEGER,  INTENT(in) :: ConservativeProjectionStatusFinal(:)
    REAL(DP), INTENT(in) :: &
      MatterCollisionRHS_FinalStage &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4))
    REAL(DP), INTENT(in) :: &
      MatterUpdate_Rate &
        (1:nQ,1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4))

    INTEGER :: iN_X, iNodeX, iX1, Unit
    REAL(DP) :: CachedSource, ClosureError, MeasuredDefect
    REAL(DP) :: OuterDefect, Radius, SolverIncrement

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# element spatial_node r_cm outer_iterations final_inner_iterations' // &
      ' first_inner_relative_residual inner_update_applied' // &
      ' accepted_radiation_relative_residual' // &
      ' conservative_projection_radiation_relative_residual' // &
      ' conservative_projection_status' // &
      ' QH_conservative_projection_final_stage_erg_cm3_s' // &
      ' QL_initial_radiation_limiter_cm3_s' // &
      ' QE_initial_radiation_limiter_erg_cm3_s' // &
      ' QH_outer_equation_defect_erg_cm3_s' // &
      ' QH_cached_final_stage_erg_cm3_s' // &
      ' QH_solver_increment_erg_cm3_s' // &
      ' QH_cached_source_minus_solver_increment_erg_cm3_s' // &
      ' QH_outer_defect_closure_error_erg_cm3_s'

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      iN_X = iN_X + 1
      Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )

      OuterDefect = MatterEquationResidual_Rate_Final(iQ_H,iN_X) &
                      / ( Erg / Centimeter**3 / Second )
      CachedSource &
        = MatterCollisionRHS_FinalStage(iQ_H,1,iNodeX,iX1,1,1) &
            / ( Erg / Centimeter**3 / Second )
      SolverIncrement = MatterUpdate_Rate(iQ_H,iNodeX,iX1,1,1) &
                          / ( Erg / Centimeter**3 / Second )
      MeasuredDefect = CachedSource - SolverIncrement
      ClosureError = OuterDefect - MeasuredDefect

      WRITE(Unit,'(2(I0,1X),ES24.16E3,1X,2(I0,1X),ES24.16E3,1X,I0,1X,' // &
                 '2(ES24.16E3,1X),I0,1X,8(ES24.16E3,1X))') &
        iX1 - iZ_B0(2) + 1, iNodeX, Radius / Centimeter, &
        OuterIterationsFinal(iN_X), InnerIterationsFinal(iN_X), &
        FirstInnerRelativeResidual(iN_X), &
        MERGE( 1, 0, InnerIterationsFinal(iN_X) > 1 ), &
        RadiationResidualRelativeFinal(iN_X), &
        ConservativeProjectionRadiationResidualRelativeFinal(iN_X), &
        ConservativeProjectionStatusFinal(iN_X), &
        ConservativeProjectionMatterHeatingRateFinal(iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        InitialRadiationLimiter_Rate_Final(iQ_Le,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        InitialRadiationLimiter_Rate_Final(iQ_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        OuterDefect, CachedSource, SolverIncrement, MeasuredDefect, ClosureError

    END DO
    END DO

    CLOSE(Unit)

  END SUBROUTINE WriteSolverConvergenceCheck


  SUBROUTINE WriteSolverConservationCheck &
    ( FileName, TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      StageFluidState, StageRadiationState, &
      MatterEquationResidual_Rate_Final, &
      InitialRadiationLimiter_Rate_Final, ConservationHandoffCheckRateFinal )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP),     INTENT(in) :: TimeStep
    INTEGER,      INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: MatterEquationResidual_Rate_Final(:,:)
    REAL(DP), INTENT(in) :: InitialRadiationLimiter_Rate_Final(:,:)
    REAL(DP), INTENT(in) :: ConservationHandoffCheckRateFinal(:,:)

    INTEGER  :: iN_X, iNodeX, iStage, iStageCandidate, iX1, Unit
    REAL(DP) :: aII, DirectEnergy, DirectLepton
    REAL(DP) :: LimiterEnergy, LimiterLepton
    REAL(DP) :: MatterEnergyA, MatterEnergyB
    REAL(DP) :: MatterLeptonA, MatterLeptonB
    REAL(DP) :: OuterEnergy, OuterLepton
    REAL(DP) :: RadiationEnergyA, RadiationEnergyB
    REAL(DP) :: RadiationLeptonA, RadiationLeptonB
    REAL(DP) :: AngularCoordinateMeasure, StageTimeStep, VolumeWeight
    REAL(DP) :: Handoff(nConservationHandoffCheck)

    iStage = 0
    DO iStageCandidate = 1, SIZE( StageFluidState, 6 )
      IF( GetImplicitStageCoefficient_IMEX_RK &
            ( iStageCandidate, iStageCandidate ) /= Zero ) &
        iStage = iStageCandidate
    END DO
    IF( iStage == 0 ) RETURN

    aII = GetImplicitStageCoefficient_IMEX_RK( iStage, iStage )
    StageTimeStep = TimeStep * aII
    AngularCoordinateMeasure &
      = MeshX(2) % Width(iZ_B0(3)) * MeshX(3) % Width(iZ_B0(4))

    CALL ComputeConservedMatterRadiationTotals &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
             StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit), &
             MatterEnergyA, RadiationEnergyA, MatterLeptonA, RadiationLeptonA )
    CALL ComputeConservedMatterRadiationTotals &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
             StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
             StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit), &
             MatterEnergyB, RadiationEnergyB, MatterLeptonB, RadiationLeptonB )

    DirectEnergy = MatterEnergyB + RadiationEnergyB &
                     - MatterEnergyA - RadiationEnergyA
    DirectLepton = MatterLeptonB + RadiationLeptonB &
                     - MatterLeptonA - RadiationLeptonA

    LimiterEnergy = Zero
    LimiterLepton = Zero
    OuterEnergy = Zero
    OuterLepton = Zero
    Handoff = Zero
    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
      iN_X = iN_X + 1
      VolumeWeight = AngularCoordinateMeasure &
        * MeshX(1) % Width(iX1) * WeightsX_q(iNodeX) &
        * GX(iNodeX,iX1,1,1,iGF_SqrtGm)
      LimiterEnergy = LimiterEnergy + VolumeWeight &
        * InitialRadiationLimiter_Rate_Final(iQ_E,iN_X) * StageTimeStep
      LimiterLepton = LimiterLepton + VolumeWeight &
        * InitialRadiationLimiter_Rate_Final(iQ_Le,iN_X) * StageTimeStep
      OuterEnergy = OuterEnergy + VolumeWeight &
        * MatterEquationResidual_Rate_Final(iQ_E,iN_X) * StageTimeStep
      OuterLepton = OuterLepton + VolumeWeight &
        * MatterEquationResidual_Rate_Final(iQ_Le,iN_X) * StageTimeStep
      Handoff = Handoff + VolumeWeight &
        * ConservationHandoffCheckRateFinal(:,iN_X) * StageTimeStep
    END DO
    END DO

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# quantity stage_pair_change input_primitive_minus_conserved' // &
      ' solver_primitive_pair_change final_conserved_minus_primitive' // &
      ' mapped_pair_change returned_increment_pair_change' // &
      ' initial_radiation_limiter_change outer_equation_old_minus_new' // &
      ' solver_closure conversion_closure increment_formation_closure' // &
      ' retained_stage_minus_returned_closure'
    WRITE(Unit,'(A,1X,12(ES24.16E3,1X))') 'energy_erg', &
      DirectEnergy / Erg, &
      Handoff(iCHC_InputConversion_E) / Erg, &
      Handoff(iCHC_PrimitivePair_E) / Erg, &
      Handoff(iCHC_FinalConversion_E) / Erg, &
      Handoff(iCHC_MappedPair_E) / Erg, &
      Handoff(iCHC_ReturnedPair_E) / Erg, &
      LimiterEnergy / Erg, OuterEnergy / Erg, &
      ( Handoff(iCHC_PrimitivePair_E) + OuterEnergy ) / Erg, &
      ( Handoff(iCHC_MappedPair_E) - Handoff(iCHC_PrimitivePair_E) &
          - Handoff(iCHC_InputConversion_E) &
          - Handoff(iCHC_FinalConversion_E) ) / Erg, &
      ( Handoff(iCHC_ReturnedPair_E) - Handoff(iCHC_MappedPair_E) ) / Erg, &
      ( DirectEnergy - Handoff(iCHC_ReturnedPair_E) ) / Erg
    WRITE(Unit,'(A,1X,12(ES24.16E3,1X))') 'lepton', &
      DirectLepton, &
      Handoff(iCHC_InputConversion_L), &
      Handoff(iCHC_PrimitivePair_L), &
      Handoff(iCHC_FinalConversion_L), &
      Handoff(iCHC_MappedPair_L), &
      Handoff(iCHC_ReturnedPair_L), &
      LimiterLepton, OuterLepton, &
      Handoff(iCHC_PrimitivePair_L) + OuterLepton, &
      Handoff(iCHC_MappedPair_L) - Handoff(iCHC_PrimitivePair_L) &
        - Handoff(iCHC_InputConversion_L) - Handoff(iCHC_FinalConversion_L), &
      Handoff(iCHC_ReturnedPair_L) - Handoff(iCHC_MappedPair_L), &
      DirectLepton - Handoff(iCHC_ReturnedPair_L)
    CLOSE(Unit)

  END SUBROUTINE WriteSolverConservationCheck


  SUBROUTINE WriteSolverConservationHandoffCheck &
    ( FileName, iZ_B0, iZ_E0, ConservationHandoffCheckRateFinal, &
      InitialRadiationLimiterRateFinal, MatterEquationResidualRateFinal )

    CHARACTER(*), INTENT(in) :: FileName
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in) :: ConservationHandoffCheckRateFinal(:,:)
    REAL(DP), INTENT(in) :: InitialRadiationLimiterRateFinal(:,:)
    REAL(DP), INTENT(in) :: MatterEquationResidualRateFinal(:,:)

    INTEGER :: iN_X, iNodeX, iX1, Unit
    REAL(DP) :: Radius

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# element spatial_node r_cm' // &
      ' QL_input_primitive_minus_conserved_cm3_s' // &
      ' QL_solver_primitive_pair_cm3_s' // &
      ' QL_final_conserved_minus_primitive_cm3_s' // &
      ' QL_mapped_pair_cm3_s QL_returned_pair_cm3_s' // &
      ' QL_initial_limiter_cm3_s QL_outer_old_minus_new_cm3_s' // &
      ' QL_solver_closure_cm3_s QL_conversion_closure_cm3_s' // &
      ' QL_increment_formation_closure_cm3_s' // &
      ' QE_input_primitive_minus_conserved_erg_cm3_s' // &
      ' QE_solver_primitive_pair_erg_cm3_s' // &
      ' QE_final_conserved_minus_primitive_erg_cm3_s' // &
      ' QE_mapped_pair_erg_cm3_s QE_returned_pair_erg_cm3_s' // &
      ' QE_initial_limiter_erg_cm3_s QE_outer_old_minus_new_erg_cm3_s' // &
      ' QE_solver_closure_erg_cm3_s QE_conversion_closure_erg_cm3_s' // &
      ' QE_increment_formation_closure_erg_cm3_s'

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      iN_X = iN_X + 1
      Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )
      WRITE(Unit,'(2(I0,1X),21(ES24.16E3,1X))') &
        iX1 - iZ_B0(2) + 1, iNodeX, Radius / Centimeter, &
        ConservationHandoffCheckRateFinal(iCHC_InputConversion_L,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_PrimitivePair_L,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_FinalConversion_L,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_MappedPair_L,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_ReturnedPair_L,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        InitialRadiationLimiterRateFinal(iQ_Le,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        MatterEquationResidualRateFinal(iQ_Le,iN_X) &
          / ( One / Centimeter**3 / Second ), &
        ( ConservationHandoffCheckRateFinal(iCHC_PrimitivePair_L,iN_X) &
            + MatterEquationResidualRateFinal(iQ_Le,iN_X) ) &
          / ( One / Centimeter**3 / Second ), &
        ( ConservationHandoffCheckRateFinal(iCHC_MappedPair_L,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_PrimitivePair_L,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_InputConversion_L,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_FinalConversion_L,iN_X) ) &
          / ( One / Centimeter**3 / Second ), &
        ( ConservationHandoffCheckRateFinal(iCHC_ReturnedPair_L,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_MappedPair_L,iN_X) ) &
          / ( One / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_InputConversion_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_PrimitivePair_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_FinalConversion_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_MappedPair_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        ConservationHandoffCheckRateFinal(iCHC_ReturnedPair_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        InitialRadiationLimiterRateFinal(iQ_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        MatterEquationResidualRateFinal(iQ_E,iN_X) &
          / ( Erg / Centimeter**3 / Second ), &
        ( ConservationHandoffCheckRateFinal(iCHC_PrimitivePair_E,iN_X) &
            + MatterEquationResidualRateFinal(iQ_E,iN_X) ) &
          / ( Erg / Centimeter**3 / Second ), &
        ( ConservationHandoffCheckRateFinal(iCHC_MappedPair_E,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_PrimitivePair_E,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_InputConversion_E,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_FinalConversion_E,iN_X) ) &
          / ( Erg / Centimeter**3 / Second ), &
        ( ConservationHandoffCheckRateFinal(iCHC_ReturnedPair_E,iN_X) &
            - ConservationHandoffCheckRateFinal(iCHC_MappedPair_E,iN_X) ) &
          / ( Erg / Centimeter**3 / Second )

    END DO
    END DO

    CLOSE(Unit)

  END SUBROUTINE WriteSolverConservationHandoffCheck


  SUBROUTINE WriteIMEXStageCheck &
    ( FileName, TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      StageFluidState, StageRadiationState )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP), INTENT(in) :: TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)

    INTEGER :: iC, iE, iN_X, iNodeE, iNodeX, iNodeZ
    INTEGER :: iStateX1, iStateX2, iStateX3
    INTEGER :: iS, iStage, iState, iX1, nX_G, Unit
    REAL(DP) :: E, ElementVolume, Eint, ImplicitTimeStep
    REAL(DP) :: Radius, VolumeWeight, W3
    REAL(DP) :: D_Element, Eint_Element, T_Element, Ye_Element, Ym_Element
    REAL(DP) :: RadiationEnergy_Element(nSpecies)
    REAL(DP) :: Q_Element(nQ,1+nSpecies+nC)
    REAL(DP), ALLOCATABLE :: AF_Stage(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF_Stage(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR_Stage(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Stage(:,:,:)

    nX_G = nDOFX * PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

    ALLOCATE &
      ( PF_Stage(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nPF) )
    ALLOCATE &
      ( AF_Stage(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nAF) )
    ALLOCATE &
      ( PR_Stage(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nPR,1:nSpecies) )
    ALLOCATE( Q_Stage(1:nQ,1:1+nSpecies+nC,1:nX_G) )

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') '# state 1: assembled stage state before limiter'
    WRITE(Unit,'(A)') '# state 2: state entering the implicit collision solve'
    WRITE(Unit,'(A)') '# state 3: state after implicit solve, before limiter'
    WRITE(Unit,'(A)') '# state 4: state after the post-implicit limiter'
    WRITE(Unit,'(A)') '# state 5: completed step after any final assembly and limiter'
    WRITE(Unit,'(A)',ADVANCE='NO') &
      '# stage state element r_cm volume_cm3 D_g_cm3 T_MeV Ye Ym' // &
      ' implicit_dt_s eint_erg_cm3 QLe_total_cm3_s QLmu_total_cm3_s' // &
      ' QE_total_erg_cm3_s QH_total_erg_cm3_s'
    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' Erad_s', iS, '_erg_cm3'
    END DO
    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' QH_s', iS, '_erg_cm3_s'
    END DO
    DO iC = 1, nC
      WRITE(Unit,'(A)',ADVANCE='NO') &
        ' QH_' // TRIM( ProcessName(iC) ) // '_erg_cm3_s'
    END DO
    WRITE(Unit,*)

    DO iStage = 1, SIZE( StageFluidState, 6 )
    ImplicitTimeStep = TimeStep &
      * GetImplicitStageCoefficient_IMEX_RK( iStage )
    DO iState = iIMEX_Assembled, iIMEX_CompletedStep

      IF(      iState == iIMEX_CompletedStep &
          .AND. iStage /= SIZE( StageFluidState, 6 ) ) CYCLE
      IF(      iState == iIMEX_CompletedStep &
          .AND. .NOT. HasFinalAssembly_IMEX_RK() ) CYCLE

      PF_Stage = Zero
      AF_Stage = Zero
      PR_Stage = Zero
      Q_Stage  = Zero

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iZ_B0(2:4), iZ_E0(2:4), iZ_B1(2:4), iZ_E1(2:4), GX, &
               StageFluidState(:,:,:,:,:,iStage,iState), &
               PF_Stage, AF_Stage )

      CALL ComputePrimitiveRadiation &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF_Stage, &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iState), PR_Stage )

      CALL ComputeMatterRHS_Collisions_OrderV &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF_Stage, PR_Stage, &
               Q_Stage )

      iN_X = 0
      DO iX1 = iZ_B0(2), iZ_E0(2)

        iStateX1 = iX1 - iZ_B1(2) + 1
        iStateX2 = 1   - iZ_B1(3) + 1
        iStateX3 = 1   - iZ_B1(4) + 1

        ElementVolume           = Zero
        D_Element               = Zero
        T_Element               = Zero
        Ye_Element              = Zero
        Ym_Element              = Zero
        Eint_Element            = Zero
        RadiationEnergy_Element = Zero
        Q_Element               = Zero

        DO iNodeX = 1, nDOFX

          iN_X = iN_X + 1
          VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                         * GX(iNodeX,iX1,1,1,iGF_h_1) &
                         * GX(iNodeX,iX1,1,1,iGF_h_2)**2

          ElementVolume = ElementVolume + VolumeWeight
          D_Element = D_Element &
            + VolumeWeight * PF_Stage(iNodeX,iX1,1,1,iPF_D)
          T_Element = T_Element &
            + VolumeWeight * AF_Stage(iNodeX,iX1,1,1,iAF_T)
          Ye_Element = Ye_Element &
            + VolumeWeight * AF_Stage(iNodeX,iX1,1,1,iAF_Ye)
          Ym_Element = Ym_Element &
            + VolumeWeight * AF_Stage(iNodeX,iX1,1,1,iAF_Ym)

          Eint = StageFluidState &
            (iNodeX,iStateX1,iStateX2,iStateX3,iCF_E,iStage,iState) &
            - Half * ( &
                StageFluidState &
                  (iNodeX,iStateX1,iStateX2,iStateX3,iCF_S1,iStage,iState)**2 &
                  / GX(iNodeX,iX1,1,1,iGF_Gm_dd_11) &
              + StageFluidState &
                  (iNodeX,iStateX1,iStateX2,iStateX3,iCF_S2,iStage,iState)**2 &
                  / GX(iNodeX,iX1,1,1,iGF_Gm_dd_22) &
              + StageFluidState &
                  (iNodeX,iStateX1,iStateX2,iStateX3,iCF_S3,iStage,iState)**2 &
                  / GX(iNodeX,iX1,1,1,iGF_Gm_dd_33) ) &
              / StageFluidState &
                  (iNodeX,iStateX1,iStateX2,iStateX3,iCF_D,iStage,iState)
          Eint_Element = Eint_Element + VolumeWeight * Eint
          Q_Element = Q_Element + VolumeWeight * Q_Stage(:,:,iN_X)

          DO iS = 1, nSpecies
          DO iE = iZ_B0(1), iZ_E0(1)
          DO iNodeE = 1, nDOFE

            iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
            E = NodeCoordinate( MeshE, iE, iNodeE )
            W3 = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E**3 &
                 / ( PlanckConstant * SpeedOfLight )**3
            RadiationEnergy_Element(iS) = RadiationEnergy_Element(iS) &
              + VolumeWeight * W3 &
                * PR_Stage(iNodeZ,iE,iX1,1,1,iPR_D,iS)

          END DO
          END DO
          END DO

        END DO

        Radius                  = MeshX(1) % Center(iX1)
        D_Element               = D_Element / ElementVolume
        T_Element               = T_Element / ElementVolume
        Ye_Element              = Ye_Element / ElementVolume
        Ym_Element              = Ym_Element / ElementVolume
        Eint_Element            = Eint_Element / ElementVolume
        RadiationEnergy_Element = RadiationEnergy_Element / ElementVolume
        Q_Element               = Q_Element / ElementVolume

        WRITE(Unit,'(3(I0,1X),12(ES24.16E3,1X))',ADVANCE='NO') &
          iStage, iState, iX1 - iZ_B0(2) + 1, &
          Radius / Centimeter, ElementVolume / Centimeter**3, &
          D_Element / ( Gram / Centimeter**3 ), &
          BoltzmannConstant * T_Element / MeV, Ye_Element, Ym_Element, &
          ImplicitTimeStep / Second, &
          Eint_Element / ( Erg / Centimeter**3 ), &
          Q_Element(iQ_Le,1) / ( One / Centimeter**3 / Second ), &
          Q_Element(iQ_Lm,1) / ( One / Centimeter**3 / Second ), &
          Q_Element(iQ_E,1) / ( Erg / Centimeter**3 / Second ), &
          Q_Element(iQ_H,1) / ( Erg / Centimeter**3 / Second )
        DO iS = 1, nSpecies
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            RadiationEnergy_Element(iS) / ( Erg / Centimeter**3 )
        END DO
        DO iS = 1, nSpecies
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            Q_Element(iQ_H,1+iS) / ( Erg / Centimeter**3 / Second )
        END DO
        DO iC = 1, nC
          WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
            Q_Element(iQ_H,1+nSpecies+iC) &
              / ( Erg / Centimeter**3 / Second )
        END DO
        WRITE(Unit,*)

      END DO

    END DO
    END DO

    CLOSE(Unit)
    DEALLOCATE( PF_Stage, AF_Stage, PR_Stage, Q_Stage )

  END SUBROUTINE WriteIMEXStageCheck


  SUBROUTINE WriteIMEXStreamingCheck &
    ( FileName, FaceFileName, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      MeanInvMFP, DnuMax, &
      StageFluidState, StageRadiationState )

    CHARACTER(*), INTENT(in) :: FileName, FaceFileName
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: MeanInvMFP(:,:,:)
    REAL(DP), INTENT(in) :: DnuMax
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)

    INTEGER, PARAMETER :: iStreaming_Total   = 1
    INTEGER, PARAMETER :: iStreaming_NoLLF   = 2
    INTEGER, PARAMETER :: iStreaming_Spatial = 3
    INTEGER, PARAMETER :: iStreaming_Surface = 4
    INTEGER, PARAMETER :: iStreaming_Volume  = 5
    INTEGER, PARAMETER :: nStreamingRaw      = 5
    INTEGER, PARAMETER :: iComponent_Total    = 1
    INTEGER, PARAMETER :: iComponent_SpatialCentered = 2
    INTEGER, PARAMETER :: iComponent_Surface  = 3
    INTEGER, PARAMETER :: iComponent_Volume   = 4
    INTEGER, PARAMETER :: iComponent_LLF      = 5
    INTEGER, PARAMETER :: iComponent_Observer = 6
    INTEGER, PARAMETER :: iComponent_Closure  = 7
    INTEGER, PARAMETER :: iComponent_TotalClosure = 8
    INTEGER, PARAMETER :: nStreamingComponents = 8

    REAL(DP), PARAMETER :: LLFFactor(nStreamingRaw) &
      = [ One, Zero, Zero, Zero, Zero ]
    REAL(DP), PARAMETER :: SurfaceFactor(nStreamingRaw) &
      = [ One, One, One, One, Zero ]
    REAL(DP), PARAMETER :: VolumeFactor(nStreamingRaw) &
      = [ One, One, One, Zero, One ]
    LOGICAL, PARAMETER :: IncludeObserver(nStreamingRaw) &
      = [ .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. ]

    INTEGER :: iC, iE, iElement, iG, iN_E, iN_X, iNodeE, iNodeX
    INTEGER :: iRaw, iS, iStage, iX1
    INTEGER :: nElements, nE_G, nFaces, nX_G, FaceUnit, Unit
    REAL(DP) :: ElementVolume, ProperWidthWeight, Radius, VolumeWeight
    REAL(DP) :: Q_Element(nQ,1+nSpecies,nStreamingRaw)
    REAL(DP) :: Q_Component(nQ,1+nSpecies,nStreamingComponents)
    REAL(DP), ALLOCATABLE :: FaceDiagnostics_PreLimiter(:,:,:)
    REAL(DP), ALLOCATABLE :: FaceDiagnostics_PostLimiter(:,:,:)
    REAL(DP), ALLOCATABLE :: MeanInvMFP_Element(:,:,:)
    REAL(DP), ALLOCATABLE :: ProperWidth_Element(:)
    REAL(DP), ALLOCATABLE :: E_N(:), D_N(:), T_N(:), Ye_N(:), Ym_N(:)
    REAL(DP), ALLOCATABLE :: DnuEq(:,:,:)
    REAL(DP), ALLOCATABLE :: AF_Stage(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF_Work(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF_Stage(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR_Work(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: RadiationRHS(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Streaming(:,:,:,:)
    LOGICAL, ALLOCATABLE :: ElectronChargedCurrentMask(:)
    LOGICAL, ALLOCATABLE :: MuonChargedCurrentMask(:)

    nElements = iZ_E0(2) - iZ_B0(2) + 1
    nE_G      = nDOFE * ( iZ_E0(1) - iZ_B0(1) + 1 )
    nFaces    = nElements + 1
    nX_G = nDOFX * PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

    ALLOCATE &
      ( PF_Stage(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nPF) )
    ALLOCATE &
      ( AF_Stage(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nAF) )
    ALLOCATE &
      ( CF_Work(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nCF) )
    ALLOCATE &
      ( CR_Work(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
          1:nCR,1:nSpecies) )
    ALLOCATE( RadiationRHS, MOLD = CR_Work )
    ALLOCATE( Q_Streaming(1:nQ,1:1+nSpecies,1:nX_G,1:nStreamingRaw) )
    ALLOCATE &
      ( FaceDiagnostics_PreLimiter &
          (1:nFaceDiagnostics,1:nSpecies,1:nFaces) )
    ALLOCATE &
      ( FaceDiagnostics_PostLimiter &
          (1:nFaceDiagnostics,1:nSpecies,1:nFaces) )
    ALLOCATE( MeanInvMFP_Element(1:nC,1:nSpecies,1:nElements) )
    ALLOCATE( ProperWidth_Element(1:nElements) )
    ALLOCATE( E_N(1:nE_G), D_N(1:nX_G), T_N(1:nX_G), &
              Ye_N(1:nX_G), Ym_N(1:nX_G) )
    ALLOCATE( DnuEq(1:nE_G,1:nSpecies,1:nX_G) )
    ALLOCATE( ElectronChargedCurrentMask(1:nX_G) )
    ALLOCATE( MuonChargedCurrentMask(1:nX_G) )

    DO iE = iZ_B0(1), iZ_E0(1)
    DO iNodeE = 1, nDOFE

      iN_E = iNodeE + ( iE - iZ_B0(1) ) * nDOFE
      E_N(iN_E) = GE(iNodeE,iE,iGE_Ep1)

    END DO
    END DO

    MeanInvMFP_Element = Zero
    ProperWidth_Element = Zero
    iN_X = 0
    DO iElement = 1, nElements

      iX1 = iZ_B0(2) + iElement - 1
      ElementVolume = Zero
      ProperWidthWeight = Zero

      DO iNodeX = 1, nDOFX

        iN_X = iN_X + 1
        VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                       * GX(iNodeX,iX1,1,1,iGF_h_1) &
                       * GX(iNodeX,iX1,1,1,iGF_h_2)**2
        ElementVolume = ElementVolume + VolumeWeight
        ProperWidthWeight = ProperWidthWeight &
          + VolumeWeight * GX(iNodeX,iX1,1,1,iGF_h_1)

        DO iS = 1, nSpecies
        DO iC = 1, nC
          MeanInvMFP_Element(iC,iS,iElement) &
            = MeanInvMFP_Element(iC,iS,iElement) &
                + VolumeWeight * MeanInvMFP(iC,iS,iN_X)
        END DO
        END DO

      END DO

      MeanInvMFP_Element(:,:,iElement) &
        = MeanInvMFP_Element(:,:,iElement) / ElementVolume
      ProperWidth_Element(iElement) &
        = MeshX(1) % Width(iX1) * ProperWidthWeight / ElementVolume

    END DO

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# Explicit radiation RHS recomputed at each post-implicit-limiter stage state.'
    WRITE(Unit,'(A)') &
      '# total = spatial_centered + LLF_penalty + observer; ' // &
      'spatial_centered = surface_centered + volume.'
    WRITE(Unit,'(A)',ADVANCE='NO') &
      '# stage explicit_rhs_used element r_cm volume_cm3 species'
    CALL WriteStreamingHeader( Unit, 'QE'  )
    CALL WriteStreamingHeader( Unit, 'QH'  )
    CALL WriteStreamingHeader( Unit, 'QM1' )
    WRITE(Unit,*)

    OPEN( NEWUNIT = FaceUnit, FILE = TRIM( FaceFileName ), &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    CALL WriteIMEXFaceHeader( FaceUnit )

    DO iStage = 1, SIZE( StageFluidState, 6 )

      PF_Stage = Zero
      AF_Stage = Zero

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iZ_B0(2:4), iZ_E0(2:4), iZ_B1(2:4), iZ_E1(2:4), GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               PF_Stage, AF_Stage )

      iN_X = 0
      DO iX1 = iZ_B0(2), iZ_E0(2)
      DO iNodeX = 1, nDOFX

        iN_X = iN_X + 1
        D_N(iN_X) = PF_Stage(iNodeX,iX1,1,1,iPF_D)
        T_N(iN_X) = AF_Stage(iNodeX,iX1,1,1,iAF_T)
        Ye_N(iN_X) = AF_Stage(iNodeX,iX1,1,1,iAF_Ye)
        Ym_N(iN_X) = AF_Stage(iNodeX,iX1,1,1,iAF_Ym)

      END DO
      END DO

      CALL SetChargedCurrentMasks &
             ( nX_G, nDOFX, D_N, T_N, Ye_N, Ym_N, &
               ElectronChargedCurrentMask, MuonChargedCurrentMask )
      CALL ComputeEquilibriumDistributions &
             ( 1, nE_G, 1, nSpecies, 1, nX_G, &
               E_N, D_N, T_N, Ye_N, Ym_N, DnuEq, &
               ElectronChargedCurrentMask_Option &
                 = ElectronChargedCurrentMask, &
               MuonChargedCurrentMask_Option = MuonChargedCurrentMask )

      CALL LimitEquilibriumDistributions_DG &
             ( 1, nE_G, 1, nSpecies, 1, nX_G, &
               E_N, DnuEq, DnuMax )

      DO iRaw = 1, nStreamingRaw

        CF_Work &
          = StageFluidState &
              (:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter)
        CR_Work &
          = StageRadiationState &
              (:,:,:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter)

        CALL ComputeIncrement_TwoMoment_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                 CF_Work, CR_Work, RadiationRHS, &
                 LLFPenaltyFactor_Option = LLFFactor(iRaw), &
                 SurfaceTermFactor_Option = SurfaceFactor(iRaw), &
                 VolumeTermFactor_Option = VolumeFactor(iRaw), &
                 IncludeObserverCorrections_Option = IncludeObserver(iRaw) )

        CALL ComputeTransportEquivalentMatterRHS &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF_Stage, RadiationRHS, &
                 Q_Streaming(:,:,:,iRaw) )

        IF( iRaw == iStreaming_Total ) &
          CALL ComputeRadialFaceDiagnostics &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                   PF_Stage, CF_Work, CR_Work, DnuEq, &
                   FaceDiagnostics_PostLimiter )

      END DO

      CF_Work &
        = StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit)
      CR_Work &
        = StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit)

      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               CF_Work, CR_Work, RadiationRHS )

      CALL ComputeRadialFaceDiagnostics &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
               PF_Stage, CF_Work, CR_Work, DnuEq, &
               FaceDiagnostics_PreLimiter )

      CALL WriteIMEXFaceRows &
             ( FaceUnit, iStage, .FALSE., iZ_B0, &
               FaceDiagnostics_PreLimiter, MeanInvMFP_Element, &
               ProperWidth_Element )

      CALL WriteIMEXFaceRows &
             ( FaceUnit, iStage, .TRUE., iZ_B0, &
               FaceDiagnostics_PostLimiter, MeanInvMFP_Element, &
               ProperWidth_Element )

      iN_X = 0
      DO iX1 = iZ_B0(2), iZ_E0(2)

        ElementVolume = Zero
        Q_Element = Zero

        DO iNodeX = 1, nDOFX

          iN_X = iN_X + 1
          VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                         * GX(iNodeX,iX1,1,1,iGF_h_1) &
                         * GX(iNodeX,iX1,1,1,iGF_h_2)**2

          ElementVolume = ElementVolume + VolumeWeight
          Q_Element = Q_Element + VolumeWeight * Q_Streaming(:,:,iN_X,:)

        END DO

        Q_Element = Q_Element / ElementVolume
        Q_Component = Zero
        Q_Component(:,:,iComponent_Total) &
          = Q_Element(:,:,iStreaming_Total)
        Q_Component(:,:,iComponent_SpatialCentered) &
          = Q_Element(:,:,iStreaming_Spatial)
        Q_Component(:,:,iComponent_Surface) &
          = Q_Element(:,:,iStreaming_Surface)
        Q_Component(:,:,iComponent_Volume) &
          = Q_Element(:,:,iStreaming_Volume)
        Q_Component(:,:,iComponent_LLF) &
          = Q_Element(:,:,iStreaming_Total) &
              - Q_Element(:,:,iStreaming_NoLLF)
        Q_Component(:,:,iComponent_Observer) &
          = Q_Element(:,:,iStreaming_NoLLF) &
              - Q_Element(:,:,iStreaming_Spatial)
        Q_Component(:,:,iComponent_Closure) &
          = Q_Element(:,:,iStreaming_Spatial) &
              - Q_Element(:,:,iStreaming_Surface) &
              - Q_Element(:,:,iStreaming_Volume)
        Q_Component(:,:,iComponent_TotalClosure) &
          = Q_Component(:,:,iComponent_Total) &
              - Q_Component(:,:,iComponent_SpatialCentered) &
              - Q_Component(:,:,iComponent_LLF) &
              - Q_Component(:,:,iComponent_Observer)

        Radius = MeshX(1) % Center(iX1)

        DO iS = 0, nSpecies

          IF( iS == 0 )THEN
            iG = 1
          ELSE
            iG = 1 + iS
          END IF

          WRITE(Unit,'(3(I0,1X),2(ES24.16E3,1X),I0,1X)',ADVANCE='NO') &
            iStage, MERGE( 1, 0, ExplicitRHSUsed_IMEX_RK(iStage) ), &
            iX1 - iZ_B0(2) + 1, Radius / Centimeter, &
            ElementVolume / Centimeter**3, iS
          CALL WriteStreamingComponents &
                 ( Unit, Q_Component(iQ_E,iG,:), &
                   Erg / Centimeter**3 / Second )
          CALL WriteStreamingComponents &
                 ( Unit, Q_Component(iQ_H,iG,:), &
                   Erg / Centimeter**3 / Second )
          CALL WriteStreamingComponents &
                 ( Unit, Q_Component(iQ_M1,iG,:), &
                   Erg / Centimeter**3 / Second )
          WRITE(Unit,*)

        END DO

      END DO

    END DO

    CLOSE(Unit)
    CLOSE(FaceUnit)

    DEALLOCATE &
      ( PF_Stage, AF_Stage, CF_Work, CR_Work, RadiationRHS, Q_Streaming, &
        FaceDiagnostics_PreLimiter, FaceDiagnostics_PostLimiter, &
        MeanInvMFP_Element, ProperWidth_Element, &
        E_N, D_N, T_N, Ye_N, Ym_N, DnuEq, &
        ElectronChargedCurrentMask, MuonChargedCurrentMask )

  END SUBROUTINE WriteIMEXStreamingCheck


  SUBROUTINE ComputeRadialFaceDiagnostics &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF, CF, CR, DnuEq, &
      FaceDiagnostics )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(in) :: &
      CF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: DnuEq(:,:,:)
    REAL(DP), INTENT(out) :: FaceDiagnostics(:,:,:)

    INTEGER :: iCR, iE, iFace, iN_E, iN_X_Left, iN_X_Right
    INTEGER :: iNodeE, iNodeX, iNodeZ, iS, iCF
    INTEGER :: iXFace, iXLeft, iXRight
    REAL(DP) :: CF_Left(nCF), CF_Right(nCF)
    REAL(DP) :: ConservedFace_Left(nCR), ConservedFace_Right(nCR)
    REAL(DP) :: ConservedTrace_Left(nCR), ConservedTrace_Right(nCR)
    REAL(DP) :: EqConservedFace_Left(nCR), EqConservedFace_Right(nCR)
    REAL(DP) :: EqConservedNode_Left(nCR), EqConservedNode_Right(nCR)
    REAL(DP) :: EqConservedTrace_Left(nCR), EqConservedTrace_Right(nCR)
    REAL(DP) :: D_Left, D_Right, Denominator
    REAL(DP) :: EqD_Left, EqD_Right, EqI_Left(3), EqI_Right(3)
    REAL(DP) :: E, EnergyState_Left, EnergyState_Right
    REAL(DP) :: EqEnergyState_Left, EqEnergyState_Right
    REAL(DP) :: NonEqEnergyState_Left, NonEqEnergyState_Right
    REAL(DP) :: EqJump, NonEqJump
    REAL(DP) :: Flux_Left(nCR), Flux_Right(nCR)
    REAL(DP) :: FluxEnergy_Left, FluxEnergy_Right
    REAL(DP) :: Gm_dd(3), h_Left(3), h_Right(3)
    REAL(DP) :: Gm_dd_Node_Left(3), Gm_dd_Node_Right(3)
    REAL(DP) :: I_Left(3), I_Right(3)
    REAL(DP) :: V_Face(3), V_Left(3), V_Right(3)
    REAL(DP) :: V_Node_Left(3), V_Node_Right(3)
    REAL(DP) :: W3
    LOGICAL :: EquilibriumSplitAvailable

    FaceDiagnostics = Zero

    DO iFace = 1, SIZE( FaceDiagnostics, 3 )

      iXFace  = iZ_B0(2) + iFace - 1
      iXLeft  = iXFace - 1
      iXRight = iXFace
      EquilibriumSplitAvailable &
        = iFace > 1 .AND. iFace < SIZE( FaceDiagnostics, 3 )

      h_Left  = Zero
      h_Right = Zero
      CF_Left  = Zero
      CF_Right = Zero

      DO iNodeX = 1, nDOFX

        h_Left(1) = h_Left(1) &
          + LX_X1_Up(1,iNodeX) * GX(iNodeX,iXLeft,1,1,iGF_h_1)
        h_Left(2) = h_Left(2) &
          + LX_X1_Up(1,iNodeX) * GX(iNodeX,iXLeft,1,1,iGF_h_2)
        h_Left(3) = h_Left(3) &
          + LX_X1_Up(1,iNodeX) * GX(iNodeX,iXLeft,1,1,iGF_h_3)
        h_Right(1) = h_Right(1) &
          + LX_X1_Dn(1,iNodeX) * GX(iNodeX,iXRight,1,1,iGF_h_1)
        h_Right(2) = h_Right(2) &
          + LX_X1_Dn(1,iNodeX) * GX(iNodeX,iXRight,1,1,iGF_h_2)
        h_Right(3) = h_Right(3) &
          + LX_X1_Dn(1,iNodeX) * GX(iNodeX,iXRight,1,1,iGF_h_3)

        DO iCF = 1, nCF
          CF_Left(iCF) = CF_Left(iCF) &
            + LX_X1_Up(1,iNodeX) * CF(iNodeX,iXLeft,1,1,iCF)
          CF_Right(iCF) = CF_Right(iCF) &
            + LX_X1_Dn(1,iNodeX) * CF(iNodeX,iXRight,1,1,iCF)
        END DO

      END DO

      Gm_dd = MAX( ( Half * ( h_Left + h_Right ) )**2, SqrtTiny )

      V_Left(1)  = CF_Left(iCF_S1) &
        / ( Gm_dd(1) * CF_Left(iCF_D) )
      V_Left(2)  = CF_Left(iCF_S2) &
        / ( Gm_dd(2) * CF_Left(iCF_D) )
      V_Left(3)  = CF_Left(iCF_S3) &
        / ( Gm_dd(3) * CF_Left(iCF_D) )
      V_Right(1) = CF_Right(iCF_S1) &
        / ( Gm_dd(1) * CF_Right(iCF_D) )
      V_Right(2) = CF_Right(iCF_S2) &
        / ( Gm_dd(2) * CF_Right(iCF_D) )
      V_Right(3) = CF_Right(iCF_S3) &
        / ( Gm_dd(3) * CF_Right(iCF_D) )

      CALL FaceVelocity_X1 &
             ( V_Left(1), V_Left(2), V_Left(3), &
               V_Right(1), V_Right(2), V_Right(3), &
               V_Face(1), V_Face(2), V_Face(3) )

      DO iS = 1, nSpecies

        Denominator = Zero

        DO iE = iZ_B0(1), iZ_E0(1)
        DO iNodeE = 1, nDOFE

          ConservedTrace_Left  = Zero
          ConservedTrace_Right = Zero
          EqConservedTrace_Left  = Zero
          EqConservedTrace_Right = Zero

          iN_E = iNodeE + ( iE - iZ_B0(1) ) * nDOFE

          DO iNodeX = 1, nDOFX

            iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

            DO iCR = 1, nCR
              ConservedTrace_Left(iCR) = ConservedTrace_Left(iCR) &
                + L_X1_Up(iNodeE,iNodeZ) &
                    * CR(iNodeZ,iE,iXLeft,1,1,iCR,iS)
              ConservedTrace_Right(iCR) = ConservedTrace_Right(iCR) &
                + L_X1_Dn(iNodeE,iNodeZ) &
                    * CR(iNodeZ,iE,iXRight,1,1,iCR,iS)
            END DO

            IF( EquilibriumSplitAvailable )THEN

              iN_X_Left = iNodeX &
                + ( iXLeft - iZ_B0(2) ) * nDOFX
              iN_X_Right = iNodeX &
                + ( iXRight - iZ_B0(2) ) * nDOFX

              V_Node_Left(1) = PF(iNodeX,iXLeft,1,1,iPF_V1)
              V_Node_Left(2) = PF(iNodeX,iXLeft,1,1,iPF_V2)
              V_Node_Left(3) = PF(iNodeX,iXLeft,1,1,iPF_V3)
              V_Node_Right(1) = PF(iNodeX,iXRight,1,1,iPF_V1)
              V_Node_Right(2) = PF(iNodeX,iXRight,1,1,iPF_V2)
              V_Node_Right(3) = PF(iNodeX,iXRight,1,1,iPF_V3)

              Gm_dd_Node_Left(1) = GX(iNodeX,iXLeft,1,1,iGF_Gm_dd_11)
              Gm_dd_Node_Left(2) = GX(iNodeX,iXLeft,1,1,iGF_Gm_dd_22)
              Gm_dd_Node_Left(3) = GX(iNodeX,iXLeft,1,1,iGF_Gm_dd_33)
              Gm_dd_Node_Right(1) = GX(iNodeX,iXRight,1,1,iGF_Gm_dd_11)
              Gm_dd_Node_Right(2) = GX(iNodeX,iXRight,1,1,iGF_Gm_dd_22)
              Gm_dd_Node_Right(3) = GX(iNodeX,iXRight,1,1,iGF_Gm_dd_33)

              CALL ComputeConserved_TwoMoment &
                     ( DnuEq(iN_E,iS,iN_X_Left), Zero, Zero, Zero, &
                       EqConservedNode_Left(iCR_N), &
                       EqConservedNode_Left(iCR_G1), &
                       EqConservedNode_Left(iCR_G2), &
                       EqConservedNode_Left(iCR_G3), &
                       V_Node_Left(1), V_Node_Left(2), V_Node_Left(3), &
                       Gm_dd_Node_Left(1), Gm_dd_Node_Left(2), &
                       Gm_dd_Node_Left(3) )

              CALL ComputeConserved_TwoMoment &
                     ( DnuEq(iN_E,iS,iN_X_Right), Zero, Zero, Zero, &
                       EqConservedNode_Right(iCR_N), &
                       EqConservedNode_Right(iCR_G1), &
                       EqConservedNode_Right(iCR_G2), &
                       EqConservedNode_Right(iCR_G3), &
                       V_Node_Right(1), V_Node_Right(2), V_Node_Right(3), &
                       Gm_dd_Node_Right(1), Gm_dd_Node_Right(2), &
                       Gm_dd_Node_Right(3) )

              DO iCR = 1, nCR
                EqConservedTrace_Left(iCR) &
                  = EqConservedTrace_Left(iCR) &
                      + L_X1_Up(iNodeE,iNodeZ) &
                          * EqConservedNode_Left(iCR)
                EqConservedTrace_Right(iCR) &
                  = EqConservedTrace_Right(iCR) &
                      + L_X1_Dn(iNodeE,iNodeZ) &
                          * EqConservedNode_Right(iCR)
              END DO

            END IF

          END DO

          CALL ComputePrimitive_TwoMoment &
                 ( ConservedTrace_Left(iCR_N), &
                   ConservedTrace_Left(iCR_G1), &
                   ConservedTrace_Left(iCR_G2), &
                   ConservedTrace_Left(iCR_G3), &
                   D_Left, I_Left(1), I_Left(2), I_Left(3), &
                   V_Face(1), V_Face(2), V_Face(3), &
                   Gm_dd(1), Gm_dd(2), Gm_dd(3) )

          CALL ComputePrimitive_TwoMoment &
                 ( ConservedTrace_Right(iCR_N), &
                   ConservedTrace_Right(iCR_G1), &
                   ConservedTrace_Right(iCR_G2), &
                   ConservedTrace_Right(iCR_G3), &
                   D_Right, I_Right(1), I_Right(2), I_Right(3), &
                   V_Face(1), V_Face(2), V_Face(3), &
                   Gm_dd(1), Gm_dd(2), Gm_dd(3) )

          CALL ComputeConserved_TwoMoment &
                 ( D_Left, I_Left(1), I_Left(2), I_Left(3), &
                   ConservedFace_Left(iCR_N), &
                   ConservedFace_Left(iCR_G1), &
                   ConservedFace_Left(iCR_G2), &
                   ConservedFace_Left(iCR_G3), &
                   V_Face(1), Zero, Zero, &
                   Gm_dd(1), Gm_dd(2), Gm_dd(3) )

          CALL ComputeConserved_TwoMoment &
                 ( D_Right, I_Right(1), I_Right(2), I_Right(3), &
                   ConservedFace_Right(iCR_N), &
                   ConservedFace_Right(iCR_G1), &
                   ConservedFace_Right(iCR_G2), &
                   ConservedFace_Right(iCR_G3), &
                   V_Face(1), Zero, Zero, &
                   Gm_dd(1), Gm_dd(2), Gm_dd(3) )

          IF( EquilibriumSplitAvailable )THEN

            CALL ComputePrimitive_TwoMoment &
                   ( EqConservedTrace_Left(iCR_N), &
                     EqConservedTrace_Left(iCR_G1), &
                     EqConservedTrace_Left(iCR_G2), &
                     EqConservedTrace_Left(iCR_G3), &
                     EqD_Left, EqI_Left(1), EqI_Left(2), EqI_Left(3), &
                     V_Face(1), V_Face(2), V_Face(3), &
                     Gm_dd(1), Gm_dd(2), Gm_dd(3) )

            CALL ComputePrimitive_TwoMoment &
                   ( EqConservedTrace_Right(iCR_N), &
                     EqConservedTrace_Right(iCR_G1), &
                     EqConservedTrace_Right(iCR_G2), &
                     EqConservedTrace_Right(iCR_G3), &
                     EqD_Right, EqI_Right(1), EqI_Right(2), EqI_Right(3), &
                     V_Face(1), V_Face(2), V_Face(3), &
                     Gm_dd(1), Gm_dd(2), Gm_dd(3) )

            CALL ComputeConserved_TwoMoment &
                   ( EqD_Left, EqI_Left(1), EqI_Left(2), EqI_Left(3), &
                     EqConservedFace_Left(iCR_N), &
                     EqConservedFace_Left(iCR_G1), &
                     EqConservedFace_Left(iCR_G2), &
                     EqConservedFace_Left(iCR_G3), &
                     V_Face(1), Zero, Zero, &
                     Gm_dd(1), Gm_dd(2), Gm_dd(3) )

            CALL ComputeConserved_TwoMoment &
                   ( EqD_Right, EqI_Right(1), EqI_Right(2), EqI_Right(3), &
                     EqConservedFace_Right(iCR_N), &
                     EqConservedFace_Right(iCR_G1), &
                     EqConservedFace_Right(iCR_G2), &
                     EqConservedFace_Right(iCR_G3), &
                     V_Face(1), Zero, Zero, &
                     Gm_dd(1), Gm_dd(2), Gm_dd(3) )

          END IF

          Flux_Left = Flux_X1 &
            ( D_Left, I_Left(1), I_Left(2), I_Left(3), &
              V_Face(1), V_Face(2), V_Face(3), &
              Gm_dd(1), Gm_dd(2), Gm_dd(3) )
          Flux_Right = Flux_X1 &
            ( D_Right, I_Right(1), I_Right(2), I_Right(3), &
              V_Face(1), V_Face(2), V_Face(3), &
              Gm_dd(1), Gm_dd(2), Gm_dd(3) )

          EnergyState_Left &
            = ConservedFace_Left(iCR_N) &
                + V_Face(1) * ConservedFace_Left(iCR_G1) &
                + V_Face(2) * ConservedFace_Left(iCR_G2) &
                + V_Face(3) * ConservedFace_Left(iCR_G3)
          EnergyState_Right &
            = ConservedFace_Right(iCR_N) &
                + V_Face(1) * ConservedFace_Right(iCR_G1) &
                + V_Face(2) * ConservedFace_Right(iCR_G2) &
                + V_Face(3) * ConservedFace_Right(iCR_G3)

          IF( EquilibriumSplitAvailable )THEN

            EqEnergyState_Left &
              = EqConservedFace_Left(iCR_N) &
                  + V_Face(1) * EqConservedFace_Left(iCR_G1) &
                  + V_Face(2) * EqConservedFace_Left(iCR_G2) &
                  + V_Face(3) * EqConservedFace_Left(iCR_G3)
            EqEnergyState_Right &
              = EqConservedFace_Right(iCR_N) &
                  + V_Face(1) * EqConservedFace_Right(iCR_G1) &
                  + V_Face(2) * EqConservedFace_Right(iCR_G2) &
                  + V_Face(3) * EqConservedFace_Right(iCR_G3)
            NonEqEnergyState_Left = EnergyState_Left - EqEnergyState_Left
            NonEqEnergyState_Right = EnergyState_Right - EqEnergyState_Right
            EqJump = EqEnergyState_Right - EqEnergyState_Left
            NonEqJump = NonEqEnergyState_Right - NonEqEnergyState_Left

          END IF

          FluxEnergy_Left &
            = Flux_Left(iCR_N) &
                + V_Face(1) * Flux_Left(iCR_G1) &
                + V_Face(2) * Flux_Left(iCR_G2) &
                + V_Face(3) * Flux_Left(iCR_G3)
          FluxEnergy_Right &
            = Flux_Right(iCR_N) &
                + V_Face(1) * Flux_Right(iCR_G1) &
                + V_Face(2) * Flux_Right(iCR_G2) &
                + V_Face(3) * Flux_Right(iCR_G3)

          E = GE(iNodeE,iE,iGE_Ep1)
          W3 = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E**3 &
                 / ( PlanckConstant * SpeedOfLight )**3

          FaceDiagnostics(iFace_EnergyLeft,iS,iFace) &
            = FaceDiagnostics(iFace_EnergyLeft,iS,iFace) &
                + W3 * EnergyState_Left
          FaceDiagnostics(iFace_EnergyRight,iS,iFace) &
            = FaceDiagnostics(iFace_EnergyRight,iS,iFace) &
                + W3 * EnergyState_Right
          FaceDiagnostics(iFace_EnergyJump,iS,iFace) &
            = FaceDiagnostics(iFace_EnergyJump,iS,iFace) &
                + W3 * ( EnergyState_Right - EnergyState_Left )
          FaceDiagnostics(iFace_EnergyJumpAbs,iS,iFace) &
            = FaceDiagnostics(iFace_EnergyJumpAbs,iS,iFace) &
                + W3 * ABS( EnergyState_Right - EnergyState_Left )
          Denominator = Denominator &
            + Half * W3 &
                * ( ABS( EnergyState_Left ) + ABS( EnergyState_Right ) )
          FaceDiagnostics(iFace_FluxCentered,iS,iFace) &
            = FaceDiagnostics(iFace_FluxCentered,iS,iFace) &
                + Half * W3 * ( FluxEnergy_Left + FluxEnergy_Right )
          FaceDiagnostics(iFace_FluxLLFPenalty,iS,iFace) &
            = FaceDiagnostics(iFace_FluxLLFPenalty,iS,iFace) &
                - Half * W3 * ( EnergyState_Right - EnergyState_Left )

          IF( EquilibriumSplitAvailable )THEN

            FaceDiagnostics(iFace_EqEnergyLeft,iS,iFace) &
              = FaceDiagnostics(iFace_EqEnergyLeft,iS,iFace) &
                  + W3 * EqEnergyState_Left
            FaceDiagnostics(iFace_EqEnergyRight,iS,iFace) &
              = FaceDiagnostics(iFace_EqEnergyRight,iS,iFace) &
                  + W3 * EqEnergyState_Right
            FaceDiagnostics(iFace_EqEnergyJump,iS,iFace) &
              = FaceDiagnostics(iFace_EqEnergyJump,iS,iFace) &
                  + W3 * EqJump
            FaceDiagnostics(iFace_EqEnergyJumpAbs,iS,iFace) &
              = FaceDiagnostics(iFace_EqEnergyJumpAbs,iS,iFace) &
                  + W3 * ABS( EqJump )
            FaceDiagnostics(iFace_NonEqEnergyLeft,iS,iFace) &
              = FaceDiagnostics(iFace_NonEqEnergyLeft,iS,iFace) &
                  + W3 * NonEqEnergyState_Left
            FaceDiagnostics(iFace_NonEqEnergyRight,iS,iFace) &
              = FaceDiagnostics(iFace_NonEqEnergyRight,iS,iFace) &
                  + W3 * NonEqEnergyState_Right
            FaceDiagnostics(iFace_NonEqEnergyJump,iS,iFace) &
              = FaceDiagnostics(iFace_NonEqEnergyJump,iS,iFace) &
                  + W3 * NonEqJump
            FaceDiagnostics(iFace_NonEqEnergyJumpAbs,iS,iFace) &
              = FaceDiagnostics(iFace_NonEqEnergyJumpAbs,iS,iFace) &
                  + W3 * ABS( NonEqJump )
            FaceDiagnostics(iFace_FluxLLFEq,iS,iFace) &
              = FaceDiagnostics(iFace_FluxLLFEq,iS,iFace) &
                  - Half * W3 * EqJump
            FaceDiagnostics(iFace_FluxLLFNonEq,iS,iFace) &
              = FaceDiagnostics(iFace_FluxLLFNonEq,iS,iFace) &
                  - Half * W3 * NonEqJump

          END IF

        END DO
        END DO

        IF( Denominator > Zero ) &
          FaceDiagnostics(iFace_EnergyJumpRelL1,iS,iFace) &
            = FaceDiagnostics(iFace_EnergyJumpAbs,iS,iFace) / Denominator

        IF( EquilibriumSplitAvailable .AND. Denominator > Zero )THEN
          FaceDiagnostics(iFace_EqEnergyJumpRel,iS,iFace) &
            = FaceDiagnostics(iFace_EqEnergyJumpAbs,iS,iFace) / Denominator
          FaceDiagnostics(iFace_NonEqEnergyJumpRel,iS,iFace) &
            = FaceDiagnostics(iFace_NonEqEnergyJumpAbs,iS,iFace) / Denominator
        END IF

        FaceDiagnostics(iFace_FluxTotal,iS,iFace) &
          = FaceDiagnostics(iFace_FluxCentered,iS,iFace) &
              + FaceDiagnostics(iFace_FluxLLFPenalty,iS,iFace)
        FaceDiagnostics(iFace_FluxLLFClosure,iS,iFace) &
          = FaceDiagnostics(iFace_FluxLLFPenalty,iS,iFace) &
              - FaceDiagnostics(iFace_FluxLLFEq,iS,iFace) &
              - FaceDiagnostics(iFace_FluxLLFNonEq,iS,iFace)

      END DO

    END DO

  END SUBROUTINE ComputeRadialFaceDiagnostics


  SUBROUTINE WriteIMEXFaceHeader( Unit )

    INTEGER, INTENT(in) :: Unit

    WRITE(Unit,'(A)') &
      '# limiter_applied=0: post-implicit state before limiter; ' // &
      'limiter_applied=1: state entering the explicit streaming operator.'
    WRITE(Unit,'(A)') &
      '# Optical depths use final-state, element-averaged inverse mean free paths.'
    WRITE(Unit,'(A)') &
      '# LLF penalty is the penalty used by the configured streaming operator.'
    WRITE(Unit,'(A)') &
      '# The equilibrium split uses the nonlinear solver''s spectrally ' // &
      'limited D_0 with the same DnuMax and face construction.'
    WRITE(Unit,'(A)') &
      '# stage limiter_applied explicit_rhs_used face r_cm species ' // &
      'equilibrium_split_available ' // &
      'energy_left_erg_cm3 energy_right_erg_cm3 ' // &
      'energy_jump_signed_erg_cm3 energy_jump_abs_spectral_erg_cm3 ' // &
      'energy_jump_relative_L1 centered_flux_erg_cm2_s ' // &
      'LLF_penalty_flux_erg_cm2_s total_flux_erg_cm2_s ' // &
      'equilibrium_energy_left_erg_cm3 equilibrium_energy_right_erg_cm3 ' // &
      'equilibrium_energy_jump_signed_erg_cm3 ' // &
      'equilibrium_energy_jump_abs_spectral_erg_cm3 ' // &
      'equilibrium_jump_relative_to_actual_L1 ' // &
      'nonequilibrium_energy_left_erg_cm3 ' // &
      'nonequilibrium_energy_right_erg_cm3 ' // &
      'nonequilibrium_energy_jump_signed_erg_cm3 ' // &
      'nonequilibrium_energy_jump_abs_spectral_erg_cm3 ' // &
      'nonequilibrium_jump_relative_to_actual_L1 ' // &
      'equilibrium_LLF_penalty_flux_erg_cm2_s ' // &
      'nonequilibrium_LLF_penalty_flux_erg_cm2_s ' // &
      'LLF_decomposition_closure_erg_cm2_s ' // &
      'tau_iso_left tau_iso_right tau_noniso_left tau_noniso_right ' // &
      'tau_total_left tau_total_right'

  END SUBROUTINE WriteIMEXFaceHeader


  SUBROUTINE WriteIMEXFaceRows &
    ( Unit, iStage, LimiterApplied, iZ_B0, FaceDiagnostics, &
      MeanInvMFP_Element, ProperWidth_Element )

    INTEGER, INTENT(in) :: Unit, iStage, iZ_B0(4)
    LOGICAL, INTENT(in) :: LimiterApplied
    REAL(DP), INTENT(in) :: FaceDiagnostics(:,:,:)
    REAL(DP), INTENT(in) :: MeanInvMFP_Element(:,:,:)
    REAL(DP), INTENT(in) :: ProperWidth_Element(:)

    INTEGER :: iFace, iS, iXLeft
    REAL(DP) :: Radius
    REAL(DP) :: TauIso_Left, TauIso_Right
    REAL(DP) :: TauNonIso_Left, TauNonIso_Right
    REAL(DP) :: TauTotal_Left, TauTotal_Right

    DO iFace = 1, SIZE( FaceDiagnostics, 3 )

      IF( iFace == 1 )THEN
        Radius = MeshX(1) % Center(iZ_B0(2)) &
                   - Half * MeshX(1) % Width(iZ_B0(2))
      ELSE
        iXLeft = iZ_B0(2) + iFace - 2
        Radius = MeshX(1) % Center(iXLeft) &
                   + Half * MeshX(1) % Width(iXLeft)
      END IF

      DO iS = 1, nSpecies

        TauIso_Left = Zero
        TauIso_Right = Zero
        TauTotal_Left = Zero
        TauTotal_Right = Zero

        IF( iFace > 1 )THEN
          TauIso_Left = ProperWidth_Element(iFace-1) &
            * MeanInvMFP_Element(iC_Iso,iS,iFace-1)
          TauTotal_Left = ProperWidth_Element(iFace-1) &
            * SUM( MeanInvMFP_Element(:,iS,iFace-1) )
        END IF

        IF( iFace <= SIZE( ProperWidth_Element ) )THEN
          TauIso_Right = ProperWidth_Element(iFace) &
            * MeanInvMFP_Element(iC_Iso,iS,iFace)
          TauTotal_Right = ProperWidth_Element(iFace) &
            * SUM( MeanInvMFP_Element(:,iS,iFace) )
        END IF
        TauNonIso_Left = TauTotal_Left - TauIso_Left
        TauNonIso_Right = TauTotal_Right - TauIso_Right

        WRITE(Unit,'(4(I0,1X),ES24.16E3,1X,2(I0,1X))',ADVANCE='NO') &
          iStage, MERGE( 1, 0, LimiterApplied ), &
          MERGE( 1, 0, ExplicitRHSUsed_IMEX_RK(iStage) ), &
          iFace - 1, Radius / Centimeter, iS, &
          MERGE( 1, 0, iFace > 1 .AND. &
                        iFace < SIZE( FaceDiagnostics, 3 ) )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EnergyLeft,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EnergyRight,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EnergyJump,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EnergyJumpAbs,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EnergyJumpRelL1,iS,iFace)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_FluxCentered,iS,iFace) &
            / ( Erg / Centimeter**2 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_FluxLLFPenalty,iS,iFace) &
            / ( Erg / Centimeter**2 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_FluxTotal,iS,iFace) &
            / ( Erg / Centimeter**2 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EqEnergyLeft,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EqEnergyRight,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EqEnergyJump,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EqEnergyJumpAbs,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_EqEnergyJumpRel,iS,iFace)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_NonEqEnergyLeft,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_NonEqEnergyRight,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_NonEqEnergyJump,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_NonEqEnergyJumpAbs,iS,iFace) &
            / ( Erg / Centimeter**3 )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_NonEqEnergyJumpRel,iS,iFace)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_FluxLLFEq,iS,iFace) &
            / ( Erg / Centimeter**2 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_FluxLLFNonEq,iS,iFace) &
            / ( Erg / Centimeter**2 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          FaceDiagnostics(iFace_FluxLLFClosure,iS,iFace) &
            / ( Erg / Centimeter**2 / Second )
        WRITE(Unit,'(6(ES24.16E3,1X))') &
          TauIso_Left, TauIso_Right, TauNonIso_Left, TauNonIso_Right, &
          TauTotal_Left, TauTotal_Right

      END DO

    END DO

  END SUBROUTINE WriteIMEXFaceRows


  SUBROUTINE WriteStreamingHeader( Unit, Quantity )

    INTEGER,      INTENT(in) :: Unit
    CHARACTER(*), INTENT(in) :: Quantity

    WRITE(Unit,'(8(A))',ADVANCE='NO') &
      ' ' // TRIM( Quantity ) // '_total_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_spatial_centered_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_surface_centered_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_volume_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_LLF_penalty_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_observer_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_spatial_closure_erg_cm3_s', &
      ' ' // TRIM( Quantity ) // '_total_closure_erg_cm3_s'

  END SUBROUTINE WriteStreamingHeader


  SUBROUTINE WriteStreamingComponents( Unit, Components, UnitValue )

    INTEGER,  INTENT(in) :: Unit
    REAL(DP), INTENT(in) :: Components(:), UnitValue

    INTEGER :: iComponent

    DO iComponent = 1, SIZE( Components )
      WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
        Components(iComponent) / UnitValue
    END DO

  END SUBROUTINE WriteStreamingComponents


  SUBROUTINE WriteEmAbSpectrum &
    ( FileName, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
      Spectrum_PreLimiter, Spectrum )

    CHARACTER(*), INTENT(in) :: FileName
    INTEGER,      INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: Spectrum_PreLimiter(:,:,:,:), Spectrum(:,:,:,:)

    INTEGER :: iE, iN_E, iN_X, iNodeE, iNodeX, iS, iX1, Unit
    REAL(DP) :: Energy, Radius, VolumeWeight

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )

    WRITE(Unit,'(A)') &
      '# iX iNodeX r_cm volume_weight_cm3 iE iNodeE E_MeV species' // &
      ' Dnu_pre Dnu_post Dnu_eq_pre Dnu_eq_post' // &
      ' Chi_pre_s_inv Chi_post_s_inv' // &
      ' QH_pre_erg_cm3_s QH_post_erg_cm3_s' // &
      ' EIn_post_erg_cm3_s EOut_post_erg_cm3_s'

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      iN_X = iN_X + 1
      Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )
      VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                     * GX(iNodeX,iX1,1,1,iGF_h_1) &
                     * GX(iNodeX,iX1,1,1,iGF_h_2)**2

      DO iS = 1, nSpecies
      DO iE = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE

        iN_E = iNodeE + ( iE - iZ_B0(1) ) * nDOFE
        Energy = NodeCoordinate( MeshE, iE, iNodeE )

        WRITE(Unit,'(I0,1X)',ADVANCE='NO') iX1
        WRITE(Unit,'(I0,1X)',ADVANCE='NO') iNodeX
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') Radius / Centimeter
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          VolumeWeight / Centimeter**3
        WRITE(Unit,'(I0,1X)',ADVANCE='NO') iE
        WRITE(Unit,'(I0,1X)',ADVANCE='NO') iNodeE
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') Energy / MeV
        WRITE(Unit,'(I0,1X)',ADVANCE='NO') iS
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum_PreLimiter(iEmAbSpectrum_Dnu,iN_E,iS,iN_X)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum(iEmAbSpectrum_Dnu,iN_E,iS,iN_X)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum_PreLimiter(iEmAbSpectrum_DnuEq,iN_E,iS,iN_X)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum(iEmAbSpectrum_DnuEq,iN_E,iS,iN_X)
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum_PreLimiter(iEmAbSpectrum_Chi,iN_E,iS,iN_X) &
            / ( One / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum(iEmAbSpectrum_Chi,iN_E,iS,iN_X) / ( One / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum_PreLimiter(iEmAbSpectrum_QH,iN_E,iS,iN_X) &
            / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum(iEmAbSpectrum_QH,iN_E,iS,iN_X) &
            / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          Spectrum(iEmAbSpectrum_EIn,iN_E,iS,iN_X) &
            / ( Erg / Centimeter**3 / Second )
        WRITE(Unit,'(ES24.16E3,1X)') &
          Spectrum(iEmAbSpectrum_EOut,iN_E,iS,iN_X) &
            / ( Erg / Centimeter**3 / Second )

      END DO
      END DO
      END DO

    END DO
    END DO

    CLOSE(Unit)

  END SUBROUTINE WriteEmAbSpectrum


  SUBROUTINE WriteDiscreteStepBudget &
    ( FileName, TimeStep, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
      CF_BeforeStep, CR_BeforeStep, CF_AfterStep, CR_AfterStep, &
      StageFluidState, StageRadiationState )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP), INTENT(in) :: TimeStep
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_BeforeStep(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR_BeforeStep(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: &
      CF_AfterStep(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR_AfterStep(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in) :: StageFluidState(:,:,:,:,:,:,:)
    REAL(DP), INTENT(in) :: StageRadiationState(:,:,:,:,:,:,:,:,:)

    INTEGER, PARAMETER :: iRaw_Total   = 1
    INTEGER, PARAMETER :: iRaw_NoLLF   = 2
    INTEGER, PARAMETER :: iRaw_Centered = 3
    INTEGER, PARAMETER :: iRaw_Surface = 4
    INTEGER, PARAMETER :: iRaw_Volume  = 5
    INTEGER, PARAMETER :: nRaw = 5
    REAL(DP), PARAMETER :: LLFFactor(nRaw) &
      = [ One, Zero, Zero, Zero, Zero ]
    REAL(DP), PARAMETER :: SurfaceFactor(nRaw) &
      = [ One, One, One, One, Zero ]
    REAL(DP), PARAMETER :: VolumeFactor(nRaw) &
      = [ One, One, One, Zero, One ]
    LOGICAL, PARAMETER :: IncludeObserver(nRaw) &
      = [ .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. ]

    INTEGER :: iC, iRaw, iStage, jStage, nStages, nX_G, Unit
    REAL(DP) :: aii, wEX, wIM
    REAL(DP), ALLOCATABLE :: AF_Pre(:,:,:,:,:), AF_Post(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF_Pre(:,:,:,:,:), PF_Post(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR_Pre(:,:,:,:,:,:,:), PR_Post(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF_Work(:,:,:,:,:), CR_Work(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: RadiationRate(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Collision_Pre(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Collision_Post(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Collision_PostLimiter(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Streaming(:,:,:,:), Q_Component(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Matter_Actual(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Radiation_Actual(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Conservation(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Limiter_Matter(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Limiter_Radiation(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Weighted_Streaming(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Weighted_Collision_Matter(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Weighted_Collision_Radiation(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Full_Matter(:,:,:), Q_Full_Radiation(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Remainder_Matter(:,:,:)
    REAL(DP), ALLOCATABLE :: Q_Remainder_Radiation(:,:,:)
    REAL(DP), ALLOCATABLE :: QH_Exact(:), QH_Exact_Full(:)
    REAL(DP), ALLOCATABLE :: QH_Remainder_Exact(:)

    nStages = SIZE( StageFluidState, 6 )
    nX_G = nDOFX * PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

    ALLOCATE &
      ( PF_Pre(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nPF) )
    ALLOCATE( PF_Post, MOLD = PF_Pre )
    ALLOCATE &
      ( AF_Pre(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4),1:nAF) )
    ALLOCATE( AF_Post, MOLD = AF_Pre )
    ALLOCATE &
      ( PR_Pre(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies) )
    ALLOCATE( PR_Post, MOLD = PR_Pre )
    ALLOCATE( CF_Work, MOLD = CF_BeforeStep )
    ALLOCATE( CR_Work, MOLD = CR_BeforeStep )
    ALLOCATE( RadiationRate, MOLD = CR_BeforeStep )
    ALLOCATE( Q_Collision_Pre(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( Q_Collision_Post(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( Q_Collision_PostLimiter(nQ,1+nSpecies+nC,nX_G) )
    ALLOCATE( Q_Streaming(nQ,1+nSpecies,nX_G,nRaw) )
    ALLOCATE( Q_Component(nQ,1+nSpecies,nX_G) )
    ALLOCATE( Q_Matter_Actual(nQ,1,nX_G) )
    ALLOCATE( Q_Radiation_Actual(nQ,1+nSpecies,nX_G) )
    ALLOCATE( Q_Conservation(nQ,1,nX_G) )
    ALLOCATE( Q_Limiter_Matter(nQ,1,nX_G) )
    ALLOCATE( Q_Limiter_Radiation(nQ,1+nSpecies,nX_G) )
    ALLOCATE( Q_Weighted_Streaming(nQ,1+nSpecies,nX_G) )
    ALLOCATE( Q_Weighted_Collision_Matter(nQ,1,nX_G) )
    ALLOCATE( Q_Weighted_Collision_Radiation(nQ,1+nSpecies,nX_G) )
    ALLOCATE( Q_Full_Matter(nQ,1,nX_G), Q_Full_Radiation(nQ,1+nSpecies,nX_G) )
    ALLOCATE( Q_Remainder_Matter(nQ,1,nX_G) )
    ALLOCATE( Q_Remainder_Radiation(nQ,1+nSpecies,nX_G) )
    ALLOCATE( QH_Exact(nX_G), QH_Exact_Full(nX_G) )
    ALLOCATE( QH_Remainder_Exact(nX_G) )

    Q_Weighted_Streaming          = Zero
    Q_Weighted_Collision_Matter   = Zero
    Q_Weighted_Collision_Radiation = Zero

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') '# All rates are element-volume averages.'
    WRITE(Unit,'(A)') &
      '# Radiation rows use the sign of the radiation equation; collision-source rows use matter-gain sign.'
    WRITE(Unit,'(A)') &
      '# Hence COLLISION_CONSERVATION = COLLISION_MATTER_ACTUAL + COLLISION_RADIATION_ACTUAL.'
    WRITE(Unit,'(A,1X,ES24.16E3)') '# timestep_s', TimeStep / Second
    DO iStage = 1, nStages
      WRITE(Unit,'(A,I0,A,3(ES24.16E3,1X))') '# stage ', iStage, &
        ' aii_wEX_wIM ', GetImplicitStageCoefficient_IMEX_RK(iStage), &
        GetExplicitStageWeight_IMEX_RK(iStage), &
        GetImplicitStageWeight_IMEX_RK(iStage)
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') '# stage ', iStage, ' aEX_row '
      DO jStage = 1, nStages
        WRITE(Unit,'(ES24.16E3,1X)',ADVANCE='NO') &
          GetExplicitStageCoefficient_IMEX_RK(iStage,jStage)
      END DO
      WRITE(Unit,*)
    END DO
    WRITE(Unit,'(A)') &
      '# stage component element r_cm volume_cm3 species wEX wIM' // &
      ' QLe_cm-3_s-1 QLmu_cm-3_s-1 QE_erg_cm-3_s-1 QH_erg_cm-3_s-1' // &
      ' QM1_erg_cm-3_s-1 QH_exact_eint_erg_cm-3_s-1'

    DO iStage = 1, nStages

      aii = GetImplicitStageCoefficient_IMEX_RK( iStage )
      wEX = GetExplicitStageWeight_IMEX_RK( iStage )
      wIM = GetImplicitStageWeight_IMEX_RK( iStage )

      CALL BuildPrimitiveState &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               PF_Pre, AF_Pre, PR_Pre )
      CALL ComputeMatterRHS_Collisions_OrderV &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF_Pre, PR_Pre, &
               Q_Collision_Pre )

      CALL BuildPrimitiveState &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               PF_Post, AF_Post, PR_Post )

      CALL ComputeMatterRHS_Collisions_OrderV &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF_Post, PR_Post, &
               Q_Collision_Post )

      CALL WriteBudgetComponent &
             ( Unit, iStage, 'COLLISION_SOURCE_PRE', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Collision_Pre(:,1:1+nSpecies,:) )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'COLLISION_SOURCE_POST', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Collision_Post(:,1:1+nSpecies,:) )

      IF( aii /= Zero )THEN
        CALL ComputeFluidStateRate &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
                 StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
                 StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
                 PF_Post, TimeStep * aii, Q_Matter_Actual, QH_Exact )
        RadiationRate = &
          ( StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit) &
            - StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit) ) &
          / ( TimeStep * aii )
        CALL ComputeTransportEquivalentMatterRHS &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF_Post, RadiationRate, &
                 Q_Radiation_Actual )
      ELSE
        Q_Matter_Actual    = Zero
        Q_Radiation_Actual = Zero
        QH_Exact           = Zero
      END IF

      Q_Conservation(:,1,:) = Q_Matter_Actual(:,1,:) &
                               + Q_Radiation_Actual(:,1,:)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'COLLISION_MATTER_ACTUAL', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Matter_Actual, QH_Exact )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'COLLISION_RADIATION_ACTUAL', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Radiation_Actual )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'COLLISION_CONSERVATION', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Conservation )

      Q_Weighted_Collision_Matter(:,1,:) &
        = Q_Weighted_Collision_Matter(:,1,:) + wIM * Q_Matter_Actual(:,1,:)
      Q_Weighted_Collision_Radiation &
        = Q_Weighted_Collision_Radiation + wIM * Q_Radiation_Actual

      CALL BuildPrimitiveState &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               PF_Post, AF_Post, PR_Post )

      CALL ComputeMatterRHS_Collisions_OrderV &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, PF_Post, PR_Post, &
               Q_Collision_PostLimiter )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'COLLISION_SOURCE_POST_LIMITER', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               Q_Collision_PostLimiter(:,1:1+nSpecies,:) )
      DO iC = 1, nC
        CALL WriteBudgetComponent &
               ( Unit, iStage, &
                 'COLLISION_PROCESS_' // TRIM( ProcessName(iC) ), wEX, wIM, &
                 iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
                 Q_Collision_PostLimiter(:,1+nSpecies+iC:1+nSpecies+iC,:) )
      END DO

      DO iRaw = 1, nRaw
        CF_Work = StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter)
        CR_Work = StageRadiationState &
                    (:,:,:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter)
        CALL ComputeIncrement_TwoMoment_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, &
                 CF_Work, CR_Work, RadiationRate, &
                 LLFPenaltyFactor_Option = LLFFactor(iRaw), &
                 SurfaceTermFactor_Option = SurfaceFactor(iRaw), &
                 VolumeTermFactor_Option = VolumeFactor(iRaw), &
                 IncludeObserverCorrections_Option = IncludeObserver(iRaw) )
        CALL ComputeTransportEquivalentMatterRHS &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF_Post, RadiationRate, &
                 Q_Streaming(:,:,:,iRaw) )
      END DO

      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_TOTAL', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Streaming(:,:,:,iRaw_Total) )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_CENTERED', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Streaming(:,:,:,iRaw_Centered) )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_SURFACE', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Streaming(:,:,:,iRaw_Surface) )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_VOLUME', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Streaming(:,:,:,iRaw_Volume) )
      Q_Component = Q_Streaming(:,:,:,iRaw_Total) &
                      - Q_Streaming(:,:,:,iRaw_NoLLF)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_LLF', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Component )
      Q_Component = Q_Streaming(:,:,:,iRaw_NoLLF) &
                      - Q_Streaming(:,:,:,iRaw_Centered)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_OBSERVER', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Component )
      Q_Component = Q_Streaming(:,:,:,iRaw_Centered) &
                      - Q_Streaming(:,:,:,iRaw_Surface) &
                      - Q_Streaming(:,:,:,iRaw_Volume)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_DG_CLOSURE', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Component )
      Q_Component = Q_Streaming(:,:,:,iRaw_Total) &
                      - Q_Streaming(:,:,:,iRaw_Centered) &
                      - ( Q_Streaming(:,:,:,iRaw_Total) &
                          - Q_Streaming(:,:,:,iRaw_NoLLF) ) &
                      - ( Q_Streaming(:,:,:,iRaw_NoLLF) &
                          - Q_Streaming(:,:,:,iRaw_Centered) )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'STREAMING_TOTAL_CLOSURE', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Component )

      Q_Component = Q_Streaming(:,:,:,iRaw_Total) &
                      - Q_Collision_PostLimiter(:,1:1+nSpecies,:)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'SEMIDISCRETE_BALANCE_POST', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Component )

      Q_Weighted_Streaming = Q_Weighted_Streaming &
                               + wEX * Q_Streaming(:,:,:,iRaw_Total)

      CALL ComputeFluidStateRate &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_Assembled), &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PreImplicit), &
               PF_Pre, TimeStep, Q_Limiter_Matter, QH_Exact )
      RadiationRate = &
        ( StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PreImplicit) &
          - StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_Assembled) ) &
        / TimeStep
      CALL ComputeTransportEquivalentMatterRHS &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF_Pre, RadiationRate, &
               Q_Limiter_Radiation )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'PRE_STAGE_LIMITER_MATTER', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Limiter_Matter, QH_Exact )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'PRE_STAGE_LIMITER_RADIATION', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Limiter_Radiation )
      Q_Conservation(:,1,:) = Q_Limiter_Matter(:,1,:) &
                               + Q_Limiter_Radiation(:,1,:)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'PRE_STAGE_LIMITER_CONSERVATION', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Conservation )

      CALL ComputeFluidStateRate &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicit), &
               StageFluidState(:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter), &
               PF_Post, TimeStep, Q_Limiter_Matter, QH_Exact )
      RadiationRate = &
        ( StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicitLimiter) &
          - StageRadiationState(:,:,:,:,:,:,:,iStage,iIMEX_PostImplicit) ) &
        / TimeStep
      CALL ComputeTransportEquivalentMatterRHS &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF_Post, RadiationRate, &
               Q_Limiter_Radiation )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'POST_STAGE_LIMITER_MATTER', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Limiter_Matter, QH_Exact )
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'POST_STAGE_LIMITER_RADIATION', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Limiter_Radiation )
      Q_Conservation(:,1,:) = Q_Limiter_Matter(:,1,:) &
                               + Q_Limiter_Radiation(:,1,:)
      CALL WriteBudgetComponent &
             ( Unit, iStage, 'POST_STAGE_LIMITER_CONSERVATION', wEX, wIM, &
               iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Conservation )

    END DO

    CALL BuildPrimitiveState &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF_AfterStep, CR_AfterStep, &
             PF_Post, AF_Post, PR_Post )
    CALL ComputeFluidStateRate &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, &
             CF_BeforeStep, CF_AfterStep, PF_Post, TimeStep, &
             Q_Full_Matter, QH_Exact_Full )
    RadiationRate = ( CR_AfterStep - CR_BeforeStep ) / TimeStep
    CALL ComputeTransportEquivalentMatterRHS &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF_Post, RadiationRate, &
             Q_Full_Radiation )

    Q_Remainder_Matter = Q_Full_Matter - Q_Weighted_Collision_Matter
    QH_Remainder_Exact = QH_Exact_Full &
                           - Q_Weighted_Collision_Matter(iQ_H,1,:)
    Q_Remainder_Radiation = Q_Full_Radiation &
                              - Q_Weighted_Streaming &
                              - Q_Weighted_Collision_Radiation

    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_STEP_MATTER_ACTUAL', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Full_Matter, QH_Exact_Full )
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_STEP_RADIATION_ACTUAL', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Full_Radiation )
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_WEIGHTED_STREAMING', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Weighted_Streaming )
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_WEIGHTED_COLLISION_MATTER', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Weighted_Collision_Matter )
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_WEIGHTED_COLLISION_RADIATION', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Weighted_Collision_Radiation )
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_REMAINDER_MATTER', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Remainder_Matter, &
             QH_Remainder_Exact )
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_REMAINDER_RADIATION', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Remainder_Radiation )

    Q_Conservation(:,1,:) = Q_Full_Matter(:,1,:) &
      + Q_Full_Radiation(:,1,:) - Q_Weighted_Streaming(:,1,:) &
      - Q_Remainder_Matter(:,1,:) - Q_Remainder_Radiation(:,1,:)
    CALL WriteBudgetComponent &
           ( Unit, 0, 'FULL_ALGEBRAIC_CLOSURE', Zero, Zero, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q_Conservation )

    CLOSE(Unit)

    DEALLOCATE &
      ( PF_Pre, PF_Post, AF_Pre, AF_Post, PR_Pre, PR_Post, &
        CF_Work, CR_Work, RadiationRate, &
        Q_Collision_Pre, Q_Collision_Post, Q_Collision_PostLimiter, &
        Q_Streaming, Q_Component, &
        Q_Matter_Actual, Q_Radiation_Actual, Q_Conservation, &
        Q_Limiter_Matter, Q_Limiter_Radiation, &
        Q_Weighted_Streaming, Q_Weighted_Collision_Matter, &
        Q_Weighted_Collision_Radiation, Q_Full_Matter, Q_Full_Radiation, &
        Q_Remainder_Matter, Q_Remainder_Radiation, &
        QH_Exact, QH_Exact_Full, QH_Remainder_Exact )

  END SUBROUTINE WriteDiscreteStepBudget


  SUBROUTINE BuildPrimitiveState &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF, CR, PF, AF, PR )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(out) :: &
      AF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nAF)
    REAL(DP), INTENT(out) :: &
      PR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)

    PF = Zero
    AF = Zero
    PR = Zero
    CALL ComputeFromConserved_Euler_NonRelativistic &
           ( iZ_B0(2:4), iZ_E0(2:4), iZ_B1(2:4), iZ_E1(2:4), &
             GX, CF, PF, AF )
    CALL ComputePrimitiveRadiation &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CR, PR )

  END SUBROUTINE BuildPrimitiveState


  SUBROUTINE ComputeFluidStateRate &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, CF_Before, CF_After, PF_After, &
      TimeInterval, Q_Rate, QH_Exact )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      CF_Before(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      CF_After(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      PF_After(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(in) :: TimeInterval
    REAL(DP), INTENT(out) :: Q_Rate(:,:,:), QH_Exact(:)

    INTEGER :: iN_X, iNodeX, iX1, iX2, iX3
    REAL(DP) :: EintBefore, EintAfter

    Q_Rate = Zero
    QH_Exact = Zero
    iN_X = 0
    DO iX3 = iZ_B0(4), iZ_E0(4)
    DO iX2 = iZ_B0(3), iZ_E0(3)
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
      iN_X = iN_X + 1
      Q_Rate(iQ_Le,1,iN_X) &
        = ( CF_After(iNodeX,iX1,iX2,iX3,iCF_Ne) &
            - CF_Before(iNodeX,iX1,iX2,iX3,iCF_Ne) ) / TimeInterval
      Q_Rate(iQ_Lm,1,iN_X) &
        = ( CF_After(iNodeX,iX1,iX2,iX3,iCF_Nm) &
            - CF_Before(iNodeX,iX1,iX2,iX3,iCF_Nm) ) / TimeInterval
      Q_Rate(iQ_E,1,iN_X) &
        = ( CF_After(iNodeX,iX1,iX2,iX3,iCF_E) &
            - CF_Before(iNodeX,iX1,iX2,iX3,iCF_E) ) / TimeInterval
      Q_Rate(iQ_M1,1,iN_X) &
        = ( CF_After(iNodeX,iX1,iX2,iX3,iCF_S1) &
            - CF_Before(iNodeX,iX1,iX2,iX3,iCF_S1) ) / TimeInterval
      Q_Rate(iQ_M2,1,iN_X) &
        = ( CF_After(iNodeX,iX1,iX2,iX3,iCF_S2) &
            - CF_Before(iNodeX,iX1,iX2,iX3,iCF_S2) ) / TimeInterval
      Q_Rate(iQ_M3,1,iN_X) &
        = ( CF_After(iNodeX,iX1,iX2,iX3,iCF_S3) &
            - CF_Before(iNodeX,iX1,iX2,iX3,iCF_S3) ) / TimeInterval
      Q_Rate(iQ_H,1,iN_X) = Q_Rate(iQ_E,1,iN_X) &
        - PF_After(iNodeX,iX1,iX2,iX3,iPF_V1) * Q_Rate(iQ_M1,1,iN_X) &
        - PF_After(iNodeX,iX1,iX2,iX3,iPF_V2) * Q_Rate(iQ_M2,1,iN_X) &
        - PF_After(iNodeX,iX1,iX2,iX3,iPF_V3) * Q_Rate(iQ_M3,1,iN_X)
      EintBefore = CF_Before(iNodeX,iX1,iX2,iX3,iCF_E) &
        - Half * ( CF_Before(iNodeX,iX1,iX2,iX3,iCF_S1)**2 &
                    / GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                 + CF_Before(iNodeX,iX1,iX2,iX3,iCF_S2)**2 &
                    / GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                 + CF_Before(iNodeX,iX1,iX2,iX3,iCF_S3)**2 &
                    / GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) ) &
          / CF_Before(iNodeX,iX1,iX2,iX3,iCF_D)
      EintAfter = CF_After(iNodeX,iX1,iX2,iX3,iCF_E) &
        - Half * ( CF_After(iNodeX,iX1,iX2,iX3,iCF_S1)**2 &
                    / GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                 + CF_After(iNodeX,iX1,iX2,iX3,iCF_S2)**2 &
                    / GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                 + CF_After(iNodeX,iX1,iX2,iX3,iCF_S3)**2 &
                    / GX(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) ) &
          / CF_After(iNodeX,iX1,iX2,iX3,iCF_D)
      QH_Exact(iN_X) = ( EintAfter - EintBefore ) / TimeInterval
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFluidStateRate


  SUBROUTINE WriteBudgetComponent &
    ( Unit, iStage, Component, wEX, wIM, &
      iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, Q, &
      QH_Exact_Option )

    INTEGER, INTENT(in) :: Unit, iStage
    CHARACTER(*), INTENT(in) :: Component
    REAL(DP), INTENT(in) :: wEX, wIM
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: Q(:,:,:)
    REAL(DP), INTENT(in), OPTIONAL :: QH_Exact_Option(:)

    INTEGER :: iElement, iG, iN_X, iNodeX, iX1
    REAL(DP) :: ElementVolume, Q_Element(nQ), QH_Exact_Element
    REAL(DP) :: Radius, VolumeWeight

    iN_X = 0
    DO iElement = 1, iZ_E0(2) - iZ_B0(2) + 1
      iX1 = iZ_B0(2) + iElement - 1
      ElementVolume = Zero
      DO iNodeX = 1, nDOFX
        iN_X = iN_X + 1
        VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                       * GX(iNodeX,iX1,1,1,iGF_h_1) &
                       * GX(iNodeX,iX1,1,1,iGF_h_2)**2
        ElementVolume = ElementVolume + VolumeWeight
      END DO
      Radius = MeshX(1) % Center(iX1)

      DO iG = 1, SIZE(Q,2)
        Q_Element = Zero
        QH_Exact_Element = Zero
        DO iNodeX = 1, nDOFX
          iN_X = iNodeX + ( iElement - 1 ) * nDOFX
          VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                         * GX(iNodeX,iX1,1,1,iGF_h_1) &
                         * GX(iNodeX,iX1,1,1,iGF_h_2)**2
          Q_Element = Q_Element + VolumeWeight * Q(:,iG,iN_X)
          IF( PRESENT( QH_Exact_Option ) .AND. iG == 1 ) &
            QH_Exact_Element = QH_Exact_Element &
              + VolumeWeight * QH_Exact_Option(iN_X)
        END DO
        Q_Element = Q_Element / ElementVolume
        IF( PRESENT( QH_Exact_Option ) .AND. iG == 1 )THEN
          QH_Exact_Element = QH_Exact_Element / ElementVolume
        ELSE
          QH_Exact_Element = Q_Element(iQ_H)
        END IF
        WRITE(Unit,'(I0,1X,A,1X,I0,1X,2(ES24.16E3,1X),I0,1X,' // &
                   '2(ES24.16E3,1X),6(ES24.16E3,1X))') &
          iStage, TRIM(Component), iElement, Radius / Centimeter, &
          ElementVolume / Centimeter**3, iG - 1, wEX, wIM, &
          Q_Element(iQ_Le) / ( One / Centimeter**3 / Second ), &
          Q_Element(iQ_Lm) / ( One / Centimeter**3 / Second ), &
          Q_Element(iQ_E) / ( Erg / Centimeter**3 / Second ), &
          Q_Element(iQ_H) / ( Erg / Centimeter**3 / Second ), &
          Q_Element(iQ_M1) / ( Erg / Centimeter**3 / Second ), &
          QH_Exact_Element / ( Erg / Centimeter**3 / Second )
      END DO
    END DO

  END SUBROUTINE WriteBudgetComponent


  SUBROUTINE ComputeTransportEquivalentMatterRHS &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, PF, RadiationTransportRHS, Q_Transport )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(in) :: &
      RadiationTransportRHS(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: Q_Transport(:,:,:)

    INTEGER :: iE, iG, iN_X, iNodeE, iNodeX, iNodeZ, iS
    INTEGER :: iX1, iX2, iX3
    REAL(DP) :: E, W2, W3

    Q_Transport = Zero

    DO iX3 = iZ_B0(4), iZ_E0(4)
    DO iX2 = iZ_B0(3), iZ_E0(3)
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX

      iN_X = iNodeX &
             + ( iX1 - iZ_B0(2) ) * nDOFX &
             + ( iX2 - iZ_B0(3) ) * nDOFX &
                 * ( iZ_E0(2) - iZ_B0(2) + 1 ) &
             + ( iX3 - iZ_B0(4) ) * nDOFX &
                 * ( iZ_E0(2) - iZ_B0(2) + 1 ) &
                 * ( iZ_E0(3) - iZ_B0(3) + 1 )

      DO iS = 1, nSpecies
      DO iE = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE

        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
        E = NodeCoordinate( MeshE, iE, iNodeE )
        W2 = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E**2 &
             / ( PlanckConstant * SpeedOfLight )**3
        W3 = W2 * E

        IF( iS <= iNuE_Bar ) &
          Q_Transport(iQ_Le,1+iS,iN_X) &
            = Q_Transport(iQ_Le,1+iS,iN_X) &
              + LeptonNumber(iS) * W2 * RadiationTransportRHS &
                   (iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

        IF( iS == iNuM .OR. iS == iNuM_Bar ) &
          Q_Transport(iQ_Lm,1+iS,iN_X) &
            = Q_Transport(iQ_Lm,1+iS,iN_X) &
              + LeptonNumber(iS) * W2 * RadiationTransportRHS &
                   (iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)

        Q_Transport(iQ_E,1+iS,iN_X) = Q_Transport(iQ_E,1+iS,iN_X) &
          + W3 * RadiationTransportRHS &
                   (iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)
        Q_Transport(iQ_M1,1+iS,iN_X) = Q_Transport(iQ_M1,1+iS,iN_X) &
          + W3 * RadiationTransportRHS &
                   (iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS)
        Q_Transport(iQ_M2,1+iS,iN_X) = Q_Transport(iQ_M2,1+iS,iN_X) &
          + W3 * RadiationTransportRHS &
                   (iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS)
        Q_Transport(iQ_M3,1+iS,iN_X) = Q_Transport(iQ_M3,1+iS,iN_X) &
          + W3 * RadiationTransportRHS &
                   (iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS)

      END DO
      END DO

      Q_Transport(iQ_H,1+iS,iN_X) = Q_Transport(iQ_E,1+iS,iN_X) &
        - PF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            * Q_Transport(iQ_M1,1+iS,iN_X) &
        - PF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            * Q_Transport(iQ_M2,1+iS,iN_X) &
        - PF(iNodeX,iX1,iX2,iX3,iPF_V3) &
            * Q_Transport(iQ_M3,1+iS,iN_X)

      END DO

      DO iG = 1, nSpecies
        Q_Transport(:,1,iN_X) = Q_Transport(:,1,iN_X) &
                                 + Q_Transport(:,1+iG,iN_X)
      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeTransportEquivalentMatterRHS


  SUBROUTINE WriteTransportCollisionCheck &
    ( FileName, iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, &
      Q_Matter, Q_Transport, RadiationCollisionRHS, RadiationTransportRHS )

    CHARACTER(*), INTENT(in) :: FileName
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(in) :: Q_Matter(:,:,:), Q_Transport(:,:,:)
    REAL(DP), INTENT(in) :: RadiationCollisionRHS(:,:,:,:)
    REAL(DP), INTENT(in) :: &
      RadiationTransportRHS(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER :: iE, iG, iN_E, iN_X, iNodeE, iNodeX, iNodeZ, iS, iX1
    INTEGER :: Unit, UnitSpectrum
    REAL(DP) :: CollisionQE, CollisionQH, CollisionQM1
    REAL(DP) :: E, Radius, TransportQE, TransportQH, TransportQM1
    REAL(DP) :: VolumeWeight
    REAL(DP) :: V1, V2, V3, W3

    OPEN( NEWUNIT = Unit, FILE = TRIM( FileName ), STATUS = 'REPLACE', &
          ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(Unit,'(A)') &
      '# element spatial_node r_cm volume_weight_cm3 species' // &
      ' QE_collision_erg_cm3_s QE_transport_erg_cm3_s QE_difference_erg_cm3_s' // &
      ' QM1_collision_erg_cm3_s QM1_transport_erg_cm3_s' // &
      ' QM1_difference_erg_cm3_s' // &
      ' QH_collision_erg_cm3_s QH_transport_erg_cm3_s QH_difference_erg_cm3_s'

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
      iN_X = iN_X + 1
      Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )
      VolumeWeight = FourPi * WeightsX1(iNodeX) * MeshX(1) % Width(iX1) &
                     * GX(iNodeX,iX1,1,1,iGF_h_1) &
                     * GX(iNodeX,iX1,1,1,iGF_h_2)**2

      DO iS = 0, nSpecies
        IF( iS == 0 )THEN
          iG = 1
        ELSE
          iG = 1 + iS
        END IF

        WRITE(Unit, &
          '(I0,1X,I0,1X,2(ES24.16E3,1X),I0,1X)',ADVANCE='NO') &
          iX1 - iZ_B0(2) + 1, iNodeX, Radius / Centimeter, &
          VolumeWeight / Centimeter**3, iS
        WRITE(Unit,'(9(ES24.16E3,1X))') &
          Q_Matter(iQ_E,iG,iN_X) / ( Erg / Centimeter**3 / Second ), &
          Q_Transport(iQ_E,iG,iN_X) / ( Erg / Centimeter**3 / Second ), &
          ( Q_Transport(iQ_E,iG,iN_X) - Q_Matter(iQ_E,iG,iN_X) ) &
            / ( Erg / Centimeter**3 / Second ), &
          Q_Matter(iQ_M1,iG,iN_X) / ( Erg / Centimeter**3 / Second ), &
          Q_Transport(iQ_M1,iG,iN_X) / ( Erg / Centimeter**3 / Second ), &
          ( Q_Transport(iQ_M1,iG,iN_X) - Q_Matter(iQ_M1,iG,iN_X) ) &
            / ( Erg / Centimeter**3 / Second ), &
          Q_Matter(iQ_H,iG,iN_X) / ( Erg / Centimeter**3 / Second ), &
          Q_Transport(iQ_H,iG,iN_X) / ( Erg / Centimeter**3 / Second ), &
          ( Q_Transport(iQ_H,iG,iN_X) - Q_Matter(iQ_H,iG,iN_X) ) &
            / ( Erg / Centimeter**3 / Second )
      END DO
    END DO
    END DO
    CLOSE(Unit)

    OPEN( NEWUNIT = UnitSpectrum, FILE = TRIM( FileName ) // '_spectrum', &
          STATUS = 'REPLACE', ACTION = 'WRITE', FORM = 'FORMATTED' )
    WRITE(UnitSpectrum,'(A)') &
      '# element spatial_node r_cm species E_MeV' // &
      ' QE_collision_erg_cm3_s QE_transport_erg_cm3_s QE_difference_erg_cm3_s' // &
      ' QM1_collision_erg_cm3_s QM1_transport_erg_cm3_s' // &
      ' QM1_difference_erg_cm3_s' // &
      ' QH_collision_erg_cm3_s QH_transport_erg_cm3_s QH_difference_erg_cm3_s'

    iN_X = 0
    DO iX1 = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1, nDOFX
      iN_X = iN_X + 1
      Radius = NodeCoordinate( MeshX(1), iX1, iNodeX )
      V1 = PF(iNodeX,iX1,1,1,iPF_V1)
      V2 = PF(iNodeX,iX1,1,1,iPF_V2)
      V3 = PF(iNodeX,iX1,1,1,iPF_V3)

      DO iS = 1, nSpecies
      DO iE = iZ_B0(1), iZ_E0(1)
      DO iNodeE = 1, nDOFE
        iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE
        iN_E = iNodeE + ( iE - iZ_B0(1) ) * nDOFE
        E = NodeCoordinate( MeshE, iE, iNodeE )
        W3 = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E**3 &
             / ( PlanckConstant * SpeedOfLight )**3

        CollisionQE = -W3 * RadiationCollisionRHS(iCR_N,iN_E,iS,iN_X)
        CollisionQM1 = -W3 * RadiationCollisionRHS(iCR_G1,iN_E,iS,iN_X)
        CollisionQH = CollisionQE &
          + W3 * ( V1 * RadiationCollisionRHS(iCR_G1,iN_E,iS,iN_X) &
                 + V2 * RadiationCollisionRHS(iCR_G2,iN_E,iS,iN_X) &
                 + V3 * RadiationCollisionRHS(iCR_G3,iN_E,iS,iN_X) )

        TransportQE = W3 * RadiationTransportRHS &
                               (iNodeZ,iE,iX1,1,1,iCR_N,iS)
        TransportQM1 = W3 * RadiationTransportRHS &
                                (iNodeZ,iE,iX1,1,1,iCR_G1,iS)
        TransportQH = TransportQE &
          - W3 * ( V1 * RadiationTransportRHS &
                            (iNodeZ,iE,iX1,1,1,iCR_G1,iS) &
                 + V2 * RadiationTransportRHS &
                            (iNodeZ,iE,iX1,1,1,iCR_G2,iS) &
                 + V3 * RadiationTransportRHS &
                            (iNodeZ,iE,iX1,1,1,iCR_G3,iS) )

        WRITE(UnitSpectrum, &
          '(I0,1X,I0,1X,ES24.16E3,1X,I0,1X,ES24.16E3,1X)',ADVANCE='NO') &
          iX1 - iZ_B0(2) + 1, iNodeX, Radius / Centimeter, iS, E / MeV
        WRITE(UnitSpectrum,'(9(ES24.16E3,1X))') &
          CollisionQE / ( Erg / Centimeter**3 / Second ), &
          TransportQE / ( Erg / Centimeter**3 / Second ), &
          ( TransportQE - CollisionQE ) / ( Erg / Centimeter**3 / Second ), &
          CollisionQM1 / ( Erg / Centimeter**3 / Second ), &
          TransportQM1 / ( Erg / Centimeter**3 / Second ), &
          ( TransportQM1 - CollisionQM1 ) / ( Erg / Centimeter**3 / Second ), &
          CollisionQH / ( Erg / Centimeter**3 / Second ), &
          TransportQH / ( Erg / Centimeter**3 / Second ), &
          ( TransportQH - CollisionQH ) / ( Erg / Centimeter**3 / Second )
      END DO
      END DO
      END DO
    END DO
    END DO

    CLOSE(UnitSpectrum)

  END SUBROUTINE WriteTransportCollisionCheck


  SUBROUTINE WriteHeader( Unit )

    INTEGER, INTENT(in) :: Unit

    INTEGER :: iC, iS

    WRITE(Unit,'(A)',ADVANCE='NO') &
      '# r_cm dr_cm dV_cm3 D_g_cm3 T_fluid_MeV Ye Ym' // &
      ' ne_before_cm3 ne_after_cm3 nm_before_cm3 nm_after_cm3' // &
      ' eint_before_erg_cm3 eint_after_erg_cm3' // &
      ' fractional_eint_change_last_step' // &
      ' PL_theta1 PL_theta2 PL_dEnergy_erg' // &
      ' solver_outer_iterations_max' // &
      ' solver_final_inner_iterations_min' // &
      ' solver_final_inner_iterations_max' // &
      ' solver_first_inner_relative_residual_max' // &
      ' solver_final_inner_update_fraction'

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0)',ADVANCE='NO') ' T_spec_s', iS
      WRITE(Unit,'(A,I0)',ADVANCE='NO') ' T_flux_s', iS
      WRITE(Unit,'(A,I0)',ADVANCE='NO') ' flux_factor_s', iS
    END DO

    WRITE(Unit,'(A)',ADVANCE='NO') ' QLe_applied_last_step_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QLmu_applied_last_step_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QLe_IMEX_source_total_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QLmu_IMEX_source_total_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QLe_IMEX_increment_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QLmu_IMEX_increment_cm3_s'

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') &
        ' QLe_IMEX_source_species_s', iS, '_cm3_s'
    END DO

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') &
        ' QLmu_IMEX_source_species_s', iS, '_cm3_s'
    END DO

    DO iC = 1, nC
      WRITE(Unit,'(A)',ADVANCE='NO') &
        ' QLe_IMEX_source_process_' // TRIM( ProcessName(iC) ) // '_cm3_s'
    END DO

    DO iC = 1, nC
      WRITE(Unit,'(A)',ADVANCE='NO') &
        ' QLmu_IMEX_source_process_' // TRIM( ProcessName(iC) ) // '_cm3_s'
    END DO

    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_total_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_transport_equivalent_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_transport_minus_collision_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_fresh_pre_limiter_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_cached_final_stage_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_applied_last_step_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_IMEX_source_total_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_IMEX_increment_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') ' QH_outer_equation_defect_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') &
      ' QH_conservative_projection_final_stage_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') &
      ' QH_cached_source_minus_solver_increment_erg_cm3_s'
    WRITE(Unit,'(A)',ADVANCE='NO') &
      ' QH_outer_defect_closure_error_erg_cm3_s'

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') &
        ' QH_IMEX_source_species_s', iS, '_erg_cm3_s'
    END DO

    DO iC = 1, nC
      WRITE(Unit,'(A)',ADVANCE='NO') &
        ' QH_IMEX_source_process_' // TRIM( ProcessName(iC) ) // '_erg_cm3_s'
    END DO

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0)',ADVANCE='NO') ' QH_species_s', iS
    END DO

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0)',ADVANCE='NO') ' QH_transport_species_s', iS
      WRITE(Unit,'(A,I0)',ADVANCE='NO') &
        ' QH_transport_minus_collision_s', iS
    END DO

    DO iC = 1, nC
      WRITE(Unit,'(A)',ADVANCE='NO') ' QH_process_' // TRIM( ProcessName(iC) )
    END DO

    DO iS = 1, nSpecies
    DO iC = 1, nC
      WRITE(Unit,'(A,I0)',ADVANCE='NO') ' QH_' // TRIM( ProcessName(iC) ) // '_s', iS
    END DO
    END DO

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' EmAb_Nin_s', iS, '_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' EmAb_Nout_s', iS, '_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' EmAb_Nnet_s', iS, '_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' EmAb_Ein_s', iS, '_erg_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' EmAb_Eout_s', iS, '_erg_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' EmAb_QH_s', iS, '_erg_cm3_s'
    END DO

    DO iS = 1, nSpecies
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' NNS_Nin_s', iS, '_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' NNS_Nout_s', iS, '_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' NNS_Nnet_s', iS, '_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' NNS_Ein_s', iS, '_erg_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' NNS_Eout_s', iS, '_erg_cm3_s'
      WRITE(Unit,'(A,I0,A)',ADVANCE='NO') ' NNS_QH_s', iS, '_erg_cm3_s'
    END DO

    DO iS = 1, nSpecies
    DO iC = 1, nC
      WRITE(Unit,'(A)',ADVANCE='NO') ' inv_mfp_' // TRIM( ProcessName(iC) ) // '_s'
      WRITE(Unit,'(I0,A)',ADVANCE='NO') iS, '_cm_inv'
    END DO
    END DO

    WRITE(Unit,*)

  END SUBROUTINE WriteHeader


  SUBROUTINE ComputePrimitiveRadiation &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CR, PR )

    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      PF(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPF)
    REAL(DP), INTENT(in) :: &
      CR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      PR(1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)

    INTEGER :: iNodeX, iNodeZ, iS, iZ1, iZ2, iZ3, iZ4

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( ( iNodeZ - 1 ) / nDOFE, nDOFX ) + 1

        CALL ComputePrimitive_TwoMoment &
               ( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputePrimitiveRadiation


END MODULE StandaloneMatterDiagnosticsModule
