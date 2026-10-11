PROGRAM ApplicationDriver_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, Three, &
    Pi, TwoPi, SqrtTiny
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    Second, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFE, nDOFZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    WeightsX1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE, nGE
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF, uDF, &
    nCF, iCF_S1, iCF_S3, iCF_E
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR, uAR, uGR, uDR, &
    nCR, nDR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, &
    ReadFieldsHDF
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeFromConserved_TwoMoment, &
    ComputeTimeStep_TwoMoment
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos, ONLY: &
    ComputeIncrement_TwoMoment_Implicit, &
    ComputeIncrement_TwoMoment_Explicit_AP, &
    SetNeutrinoMatterSolverConvergenceDiagnostics, &
    GetNeutrinoMatterSolverConvergenceDiagnostics, &
    nLeptonFamilies, nLeptonFlavorCheck
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    InitializeEquilibriumDeviationLLF, &
    SetEquilibriumDeviationReferenceState, &
    ClearEquilibriumDeviationReferenceState, &
    InitializeSourceAwareLLF, &
    SetSourceAwareTransportOpacity_Flattened, &
    ClearSourceAwareTransportOpacity, &
    ComputeIncrement_TwoMoment_Explicit, &
    SetStreamingCheckSurfaceTermFactors, &
    ClearStreamingCheckSurfaceTermFactors
  USE TwoMoment_NeutrinoMatterSolverModule, ONLY: &
    InitializeNeutrinoMatterSolverParameters, &
    nQ, nC
  USE TwoMoment_TimeSteppingModule, ONLY: &
    Update_IMEX_RK, &
    GetNumberOfStages_IMEX_RK, &
    nIMEXStageStates
  USE InitializationModule_Neutrinos, ONLY: &
    InitializeFields, &
    ComputeError, &
    OpaqueDiffusionDomainLength, &
    OpaqueDiffusionOpacity, &
    OpaqueDiffusionDecayRate, &
    OpaqueDiffusionAmplitude, &
    OpaqueDiffusionBackground, &
    VariableOpacityDiffusionMeanOpacity, &
    SetVariableOpacityDiffusionMeanOpticalDepth, &
    VariableOpacityDiffusionOpacity, &
    VariableOpacityDiffusionDensitySource, &
    ScatteringDiffusionDensitySource, &
    AbsorptionDiffusionAbsorptionOpacity, &
    AbsorptionDiffusionEquilibriumDensity, &
    CoupledDiffusionAbsorptionOpacity, &
    CoupledDiffusionEquilibriumSlope, &
    CoupledDiffusionDecayRate, &
    SetMovingHomogeneousRelaxationParameters
  USE StandaloneMatterDiagnosticsModule, ONLY: &
    WriteStandaloneMatterDiagnostics, &
    InitializeMatterExchangeAttribution, &
    AccumulateMatterExchangeAttribution, &
    AccumulateActualMatterExchangeAttribution, &
    FinalizeMatterExchangeAttribution, &
    InitializeCumulativeConservationDiagnostics, &
    AccumulateCumulativeConservationDiagnostics, &
    FinalizeCumulativeConservationDiagnostics
  USE TwoMoment_TallyModule, ONLY: &
    ComputeTally
  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(512) :: ProfileName
  CHARACTER(64) :: EosTableName
  CHARACTER(64) :: OpacityTableName_EmAb
  CHARACTER(64) :: OpacityTableName_Iso
  CHARACTER(64) :: OpacityTableName_NNS
  CHARACTER(64) :: OpacityTableName_NES
  CHARACTER(64) :: OpacityTableName_Pair
  CHARACTER(64) :: OpacityTableName_Brem
  LOGICAL       :: EvolveEuler
  LOGICAL       :: UseSlopeLimiter_Euler
  LOGICAL       :: UseSlopeLimiter_TwoMoment
  LOGICAL       :: UsePositivityLimiter_Euler
  LOGICAL       :: UsePositivityLimiter_TwoMoment
  LOGICAL       :: UseEnergyLimiter_TwoMoment
  LOGICAL       :: PrescribedTimeStep
  LOGICAL       :: Include_NNS
  LOGICAL       :: Include_NES
  LOGICAL       :: Include_Pair
  LOGICAL       :: Include_NuPair
  LOGICAL       :: Include_Brem
  LOGICAL       :: Include_LinCorr
  LOGICAL       :: FreezeOpacities
  LOGICAL       :: UseConservativeProjection
  LOGICAL       :: UseRadiationOnlySolve
  LOGICAL       :: UseEquilibriumDeviationLLF
  LOGICAL       :: UseSourceAwareLLF
  LOGICAL       :: UseAPWellBalancedFluxSource
  LOGICAL       :: UseAbsorptionDiffusionScaledLLF
  LOGICAL       :: MovingAllProcesses
  LOGICAL       :: ProfileAllProcesses
  LOGICAL       :: ProfileMomentumCoupling
  LOGICAL       :: ProfileRadiationLimiters
  LOGICAL       :: WriteSolverHistory
  INTEGER, PARAMETER :: iLimiterPoint_Node = 1
  INTEGER, PARAMETER :: iLimiterPoint_E_Dn = 2
  INTEGER, PARAMETER :: iLimiterPoint_E_Up = 3
  INTEGER, PARAMETER :: iLimiterPoint_X_Dn = 4
  INTEGER, PARAMETER :: iLimiterPoint_X_Up = 5
  INTEGER       :: RestartFileNumber
  INTEGER       :: nSpecies
  INTEGER       :: nNodes
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: nEquidistantX
  INTEGER       :: iCycle, iCycleD, iCycleW, maxCycles
  INTEGER       :: nIMEXStages
  INTEGER       :: M_outer, MaxIter_outer
  INTEGER       :: M_inner, MaxIter_inner
  INTEGER       :: ScatteringDiffusionNX
  INTEGER       :: EnvironmentStatus, EnvironmentLength, ReadStatus
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: dEquidistantX
  REAL(DP)      :: t, dt, dt_CFL, dt_0, dt_MAX, dt_RATE, t_end
  REAL(DP)      :: MatterDiagnosticsTimeStep
  REAL(DP)      :: Rtol_outer, Rtol_inner
  REAL(DP)      :: wMatterRHS(6)
  REAL(DP)      :: DnuMax
  REAL(DP)      :: OffGridFluidIncrement(nCF)
  REAL(DP)      :: OffGridRadiationIncrement(2*nCR)
  REAL(DP)      :: OffGridRadiationLeptonFlavorIncrement(3)
  REAL(DP)      :: AbsorptionDiffusionMeanOpticalDepth
  REAL(DP)      :: AbsorptionDiffusionLLFPenaltyFactor
  REAL(DP)      :: AbsorptionDiffusionMeanCellWidth
  REAL(DP)      :: AbsorptionDiffusionMeanCellOpticalDepth
  CHARACTER(128) :: AbsorptionDiffusionOpticalDepthString
  CHARACTER(128) :: AbsorptionDiffusionLLFScalingString
  CHARACTER(128) :: ScatteringDiffusionNXString
  CHARACTER(128) :: SpatialElementsString
  CHARACTER(128) :: CoupledDiffusionSchemeString
  CHARACTER(128) :: NodesString
  CHARACTER(128) :: SourceAwareLLFString
  CHARACTER(128) :: ConservativeProjectionString
  CHARACTER(128) :: RadiationOnlySolveString
  CHARACTER(128) :: APWellBalancedFluxSourceString
  CHARACTER(128) :: MatterExchangeShadowsString
  CHARACTER(128) :: MaxCyclesString
  CHARACTER(128) :: EndTimeMillisecondsString
  CHARACTER(128) :: TimeStepMillisecondsString
  CHARACTER(128) :: MovingVelocityString
  CHARACTER(128) :: MovingVelocityDirectionString
  CHARACTER(128) :: MovingRadiationScaleString
  CHARACTER(128) :: MovingSpeciesString
  CHARACTER(128) :: MovingAllProcessesString
  CHARACTER(128) :: ProfileSpeciesString
  CHARACTER(128) :: ProfileAllProcessesString
  CHARACTER(128) :: ProfileMomentumCouplingString
  CHARACTER(128) :: ProfileRadiationLimitersString
  CHARACTER(128) :: WriteSolverHistoryString
  CHARACTER(128) :: WriteMatterDiagnosticsString
  CHARACTER(128) :: WriteMatterExchangeAttributionString
  CHARACTER(128) :: WriteCumulativeConservationString
  CHARACTER(512) :: ProfileNameString
  CHARACTER(128) :: AndersonDepthString
  REAL(DP)      :: MovingVelocityFraction
  REAL(DP)      :: MovingVelocityDirection(3)
  REAL(DP)      :: MovingRadiationScale
  REAL(DP), ALLOCATABLE :: CF_BeforeMatterDiagnosticsStep(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CR_BeforeMatterDiagnosticsStep(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: MatterCollisionExchange_Rate_LastStep(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: MatterCollisionRHS_FinalStage(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: MatterUpdate_Rate_LastStep(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: RadiationState_BeforeFinalLimiter(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_FluidState(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_RadiationState(:,:,:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_FluidImplicitRate(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_FluidExplicitRate(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_RadiationImplicitRate(:,:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_RadiationExplicitRate(:,:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: IMEX_Stage_LeptonFlavorCheckRate(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CF_AttributionStart(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CF_AttributionNoStream(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CF_AttributionAP(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CF_AttributionStandard(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CR_AttributionNoStream(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CR_AttributionAP(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: CR_AttributionStandard(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: EquilibriumTest_CF_Initial(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: EquilibriumTest_CR_Initial(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: EquilibriumTest_CF_Reference(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: EquilibriumTest_CR_Reference(:,:,:,:,:,:,:)
  LOGICAL       :: Relaxation_restart_from_file
  LOGICAL       :: WriteMatterDiagnostics
  LOGICAL       :: WriteMatterExchangeAttribution
  LOGICAL       :: WriteCumulativeConservation
  LOGICAL       :: ComputeMatterExchangeShadows
  LOGICAL       :: CaptureMatterDiagnosticsStep
  LOGICAL       :: CaptureSolverConvergenceStep
  LOGICAL       :: RunEquilibriumTest
  LOGICAL       :: RequireStationaryEquilibrium
  LOGICAL       :: RunManufacturedOpaqueDiffusion
  LOGICAL       :: RunManufacturedVariableOpacityDiffusion
  LOGICAL       :: RunScatteringDiffusionLimit
  LOGICAL       :: RunAbsorptionDiffusionBalance
  LOGICAL       :: RunCoupledAbsorptionDiffusion
  REAL(DP), PARAMETER :: ScatteringDiffusionEfoldings = 6.0_DP
  REAL(DP), PARAMETER :: CoupledDiffusionEfoldings = 2.0_DP
  CHARACTER(256) :: MatterDiagnosticsFileName
  CHARACTER(256) :: MatterDiagnosticsFileString
  CHARACTER(128) :: ProgramNameString

  !ProgramName = 'Relaxation'
  Relaxation_restart_from_file = .FALSE. 
  ProgramName = 'DeleptonizationWave1D'
  !ProgramName = 'HomogeneousEquilibrium'
  !ProgramName = 'VaryingLocalEquilibrium'
  !ProgramName = 'AbsorptionDiffusionBalance'
  !ProgramName = 'ScatteringDiffusionLimit'
  !ProgramName = 'CoupledAbsorptionDiffusion'

  CoordinateSystem = 'CARTESIAN'

  !EosTableName          = 'wl-EOS-SFHo-15-25-50.h5'
  !EosTableName          = 'BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
  EosTableName          = 'BaryonsPlusPhotonsPlusLeptonsEOS_1.h5'
  OpacityTableName_EmAb = 'wl-Op-SFHo-15-25-50-E40-EmAb.h5'
  OpacityTableName_Iso  = 'wl-Op-SFHo-15-25-50-E40-Iso.h5'
  OpacityTableName_NNS  = 'wl-Op-SFHo-15-25-50-E40-NNS.h5'
  OpacityTableName_NES  = 'wl-Op-SFHo-15-25-50-E40-NES.h5'
  OpacityTableName_Pair = 'wl-Op-SFHo-15-25-50-E40-Pair.h5'
  OpacityTableName_Brem = 'wl-Op-SFHo-15-25-50-E40-Brem.h5'

  PrescribedTimeStep = .FALSE.

  RestartFileNumber = - 1

  M_outer         = 2
  MaxIter_outer   = 1000
  Rtol_outer      = 1.0d-8
  M_inner         = 2
  MaxIter_inner   = 1000
  Rtol_inner      = 1.0d-8
  Include_NNS     = .TRUE.
  Include_NES     = .TRUE.
  Include_Pair    = .TRUE.
  Include_NuPair  = .TRUE.
  Include_Brem    = .TRUE.
  Include_LinCorr = .FALSE.
  wMatterRHS      = [ One, One, One, One, One, One ]
  DnuMax          = HUGE( One )
  FreezeOpacities = .FALSE.
  UseConservativeProjection = .FALSE.
  UseRadiationOnlySolve = .FALSE.
  UseEquilibriumDeviationLLF = .FALSE.
  UseSourceAwareLLF = .FALSE.
  UseAPWellBalancedFluxSource = .FALSE.
  UseAbsorptionDiffusionScaledLLF = .FALSE.
  AbsorptionDiffusionLLFPenaltyFactor = One
  WriteMatterDiagnostics = .FALSE.
  WriteMatterExchangeAttribution = .FALSE.
  WriteCumulativeConservation = .FALSE.
  WriteSolverHistory = .FALSE.
  ComputeMatterExchangeShadows = .FALSE.
  MatterDiagnosticsFileName = 'matter_diagnostics.dat'
  ProfileName = ''
  RunEquilibriumTest = .FALSE.
  RequireStationaryEquilibrium = .FALSE.
  RunManufacturedOpaqueDiffusion = .FALSE.
  RunManufacturedVariableOpacityDiffusion = .FALSE.
  RunScatteringDiffusionLimit = .FALSE.
  RunAbsorptionDiffusionBalance = .FALSE.
  RunCoupledAbsorptionDiffusion = .FALSE.

  nEquidistantX = 1
  dEquidistantX = 1.0_DP * Kilometer

  ProgramNameString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_PROGRAM_NAME', ProgramNameString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 ) &
    ProgramName = TRIM( ADJUSTL( ProgramNameString ) )

  SELECT CASE( TRIM( ProgramName ) )

    CASE( 'HomogeneousEquilibrium' )

      RunEquilibriumTest = .TRUE.
      RequireStationaryEquilibrium = .TRUE.

      CoordinateSystem = 'CARTESIAN'

      nSpecies = 2
      nNodes   = 2

      nX  = [ 4, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ 4.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 1, 1, 1 ]

      nE    = 16
      eL    = 0.0_DP * MeV
      eR    = 3.0d2 * MeV
      bcE   = 10
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'BackwardEuler'

      PrescribedTimeStep = .TRUE.
      dt_0               = 1.0d-4 * Millisecond
      dt_MAX             = dt_0
      dt_RATE            = One
      t_end              = 4.0d0 * dt_0

      iCycleD  = 1
      iCycleW  = 100
      maxCycles = 4

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .TRUE.

      Include_NNS     = .FALSE.
      Include_NES     = .FALSE.
      Include_Pair    = .FALSE.
      Include_NuPair  = .FALSE.
      Include_Brem    = .FALSE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .FALSE.

      wMatterRHS = [ One, One, One, One, One, One ]
      DnuMax     = One - 1.0d-6

      M_outer       = 2
      MaxIter_outer = 100
      Rtol_outer    = 1.0d-12
      M_inner       = 2
      MaxIter_inner = 1000
      Rtol_inner    = 1.0d-12
      WriteMatterDiagnostics = .FALSE.

    CASE( 'VaryingLocalEquilibrium' )

      RunEquilibriumTest = .TRUE.
      RequireStationaryEquilibrium = .FALSE.

      CoordinateSystem = 'SPHERICAL'

      nSpecies = 2
      nNodes   = 2

      nX    = [ 80, 1, 1 ]
      xL    = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR    = [ 2.0d1 * Kilometer, Pi, TwoPi ]
      bcX   = [ 31, 1, 1 ]
      ZoomX = [ One, One, One ]

      nEquidistantX = 0

      nE    = 16
      eL    = 0.0_DP * MeV
      eR    = 3.0d2 * MeV
      bcE   = 11
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'BackwardEuler'

      PrescribedTimeStep = .TRUE.
      dt_0               = 1.0d-4 * Millisecond
      dt_MAX             = dt_0
      dt_RATE            = One
      t_end              = 4.0d0 * dt_0

      iCycleD   = 1
      iCycleW   = 1
      maxCycles = 4

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .TRUE.

      Include_NNS     = .TRUE.
      Include_NES     = .TRUE.
      Include_Pair    = .TRUE.
      Include_NuPair  = .TRUE.
      Include_Brem    = .TRUE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .FALSE.

      wMatterRHS = [ One, One, One, One, One, One ]
      DnuMax     = One - 1.0d-6

      M_outer       = 2
      MaxIter_outer = 100
      Rtol_outer    = 1.0d-12
      M_inner       = 2
      MaxIter_inner = 1000
      Rtol_inner    = 1.0d-12

      ProfileName = 'input_thornado_VX_100ms.dat'
      WriteMatterDiagnostics = .FALSE.

    CASE( 'ManufacturedOpaqueDiffusion' )

      RunManufacturedOpaqueDiffusion = .TRUE.

      CoordinateSystem = 'CARTESIAN'

      nSpecies = 2
      nNodes   = 2

      nX  = [ 16, 1, 1 ]
      xL  = [ Zero, Zero, Zero ]
      xR  = [ OpaqueDiffusionDomainLength, One, One ]
      bcX = [ 1, 1, 1 ]

      nE    = 4
      eL    = 1.0_DP * MeV
      eR    = 5.0d1 * MeV
      bcE   = 0
      ZoomE = One

      TimeSteppingScheme = 'IMEX_PDARS'

      PrescribedTimeStep = .FALSE.
      dt_0               = Zero
      dt_MAX             = HUGE(One)
      dt_RATE            = One
      t_end              = 0.4_DP / ABS( OpaqueDiffusionDecayRate() )

      iCycleD   = 100
      iCycleW   = 1000000
      maxCycles = 100000

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .FALSE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

      Include_NNS     = .FALSE.
      Include_NES     = .FALSE.
      Include_Pair    = .FALSE.
      Include_NuPair  = .FALSE.
      Include_Brem    = .FALSE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .TRUE.

      wMatterRHS = [ Zero, Zero, Zero, Zero, Zero, Zero ]
      DnuMax = One
      WriteMatterDiagnostics = .FALSE.

    CASE( 'ManufacturedVariableDiffusion' )

      RunManufacturedVariableOpacityDiffusion = .TRUE.

      CoordinateSystem = 'CARTESIAN'

      nSpecies = 2
      nNodes   = 2

      nX  = [ 16, 1, 1 ]
      xL  = [ Zero, Zero, Zero ]
      xR  = [ OpaqueDiffusionDomainLength, One, One ]
      bcX = [ 1, 1, 1 ]

      nE    = 4
      eL    = 1.0_DP * MeV
      eR    = 5.0d1 * MeV
      bcE   = 0
      ZoomE = One

      TimeSteppingScheme = 'IMEX_PDARS'

      PrescribedTimeStep = .FALSE.
      dt_0               = Zero
      dt_MAX             = HUGE(One)
      dt_RATE            = One
      t_end              = 10.0_DP * OpaqueDiffusionDomainLength

      iCycleD   = 100
      iCycleW   = 1000000
      maxCycles = 100000

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .FALSE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

      Include_NNS     = .FALSE.
      Include_NES     = .FALSE.
      Include_Pair    = .FALSE.
      Include_NuPair  = .FALSE.
      Include_Brem    = .FALSE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .TRUE.

      wMatterRHS = [ Zero, Zero, Zero, Zero, Zero, Zero ]
      DnuMax = One
      WriteMatterDiagnostics = .FALSE.

    CASE( 'ScatteringDiffusionLimit' )

      RunScatteringDiffusionLimit = .TRUE.

      CoordinateSystem = 'CARTESIAN'

      nSpecies = 2
      nNodes   = 2

      nX  = [ 16, 1, 1 ]
      xL  = [ Zero, Zero, Zero ]
      xR  = [ OpaqueDiffusionDomainLength, One, One ]
      bcX = [ 1, 1, 1 ]

      nE    = 4
      eL    = 1.0_DP * MeV
      eR    = 5.0d1 * MeV
      bcE   = 0
      ZoomE = One

      TimeSteppingScheme = 'IMEX_PDARS'

      PrescribedTimeStep = .FALSE.
      dt_0               = Zero
      dt_MAX             = HUGE(One)
      dt_RATE            = One
      t_end              = Zero

      iCycleD   = 10000
      iCycleW   = 100000000
      maxCycles = 2000000

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .FALSE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

      Include_NNS     = .FALSE.
      Include_NES     = .FALSE.
      Include_Pair    = .FALSE.
      Include_NuPair  = .FALSE.
      Include_Brem    = .FALSE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .TRUE.

      wMatterRHS = [ Zero, Zero, Zero, Zero, Zero, Zero ]
      DnuMax = One
      WriteMatterDiagnostics = .FALSE.

    CASE( 'AbsorptionDiffusionBalance' )

      RunAbsorptionDiffusionBalance = .TRUE.

      CoordinateSystem = 'CARTESIAN'

      nSpecies = 2
      nNodes   = 2

      nX  = [ 16, 1, 1 ]
      xL  = [ Zero, Zero, Zero ]
      xR  = [ OpaqueDiffusionDomainLength, One, One ]
      bcX = [ 1, 1, 1 ]

      nE    = 4
      eL    = 1.0_DP * MeV
      eR    = 5.0d1 * MeV
      bcE   = 0
      ZoomE = One

      TimeSteppingScheme = 'IMEX_PDARS'

      PrescribedTimeStep = .FALSE.
      dt_0               = Zero
      dt_MAX             = HUGE(One)
      dt_RATE            = One
      t_end              = OpaqueDiffusionDomainLength

      iCycleD   = 1
      iCycleW   = 1000000
      maxCycles = 1

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .FALSE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

      Include_NNS     = .FALSE.
      Include_NES     = .FALSE.
      Include_Pair    = .FALSE.
      Include_NuPair  = .FALSE.
      Include_Brem    = .FALSE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .TRUE.

      wMatterRHS = [ Zero, Zero, Zero, Zero, Zero, Zero ]
      DnuMax = One
      WriteMatterDiagnostics = .FALSE.

    CASE( 'CoupledAbsorptionDiffusion' )

      RunCoupledAbsorptionDiffusion = .TRUE.

      CoordinateSystem = 'CARTESIAN'

      nSpecies = 2
      nNodes   = 2

      nX  = [ 16, 1, 1 ]
      xL  = [ Zero, Zero, Zero ]
      xR  = [ OpaqueDiffusionDomainLength, One, One ]
      bcX = [ 1, 1, 1 ]

      nE    = 4
      eL    = 1.0_DP * MeV
      eR    = 5.0d1 * MeV
      bcE   = 0
      ZoomE = One

      TimeSteppingScheme = 'IMEX_PDARS'

      PrescribedTimeStep = .FALSE.
      dt_0               = Zero
      dt_MAX             = HUGE(One)
      dt_RATE            = One
      t_end              = Zero

      iCycleD   = 10000
      iCycleW   = 100000000
      maxCycles = 2000000

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .FALSE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

      Include_NNS     = .FALSE.
      Include_NES     = .FALSE.
      Include_Pair    = .FALSE.
      Include_NuPair  = .FALSE.
      Include_Brem    = .FALSE.
      Include_LinCorr = .FALSE.
      FreezeOpacities = .TRUE.

      wMatterRHS = [ Zero, Zero, Zero, Zero, Zero, Zero ]
      DnuMax = One
      WriteMatterDiagnostics = .FALSE.

    CASE( 'Relaxation' )

      nSpecies = 6
      nNodes   = 2
      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 11
      ZoomE = 1.266038160710160_DP

      IF( Relaxation_restart_from_file ) THEN
        CoordinateSystem = 'SPHERICAL'

        nX    = [ 1, 1, 1 ]
        xL    = [ 0.0_DP            , 0.0_DP, 0.0_DP ]
        xR    = [ 8.0_DP * Kilometer, Pi    , TwoPi  ]
        bcX   = [ 0, 0, 0 ]
        ZoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]

        nEquidistantX = 0
        dEquidistantX = 1.0_DP * Kilometer

        nE    = 16
        eL    = 0.0d0 * MeV
        eR    = 3.0d2 * MeV
        bcE   = 11
        ZoomE = 1.266038160710160_DP

        TimeSteppingScheme = 'IMEX_PDARS'
        !TimeSteppingScheme = 'BackwardEuler'

        t_end = 1.0d-2 * Millisecond

        PrescribedTimeStep = .TRUE. ! If .FALSE., explicit CFL will be used.
        dt_0               = 8.895d-06 * Second !8.895d-3 * Millisecond
        dt_MAX             = 8.895d-06 * Second !8.895d-3 * Millisecond
        dt_RATE            = 1.01_DP

        dt_0               = 1.111880E-06 * Second !8.895d-3 * Millisecond
        dt_MAX             = 1.111880E-06 * Second !8.895d-3 * Millisecond

        iCycleD = 1
        iCycleW = 100
        maxCycles = 1000000

        EvolveEuler                    = .FALSE.
        UseSlopeLimiter_Euler          = .FALSE.
        UseSlopeLimiter_TwoMoment      = .FALSE.
        UsePositivityLimiter_Euler     = .FALSE.
        UsePositivityLimiter_TwoMoment = .TRUE.
        UseEnergyLimiter_TwoMoment     = .TRUE.

        Include_NNS     = .TRUE.
        Include_NES     = .TRUE.
        Include_Pair    = .TRUE.
        Include_NuPair  = .TRUE.
        Include_Brem    = .TRUE.
        Include_LinCorr = .FALSE.
        FreezeOpacities = .TRUE. ! --- Keep opacities fixed during iterations?
        M_outer         = 2
        MaxIter_outer   = 1000
        Rtol_outer      = 1.0d-8
        M_inner         = 2
        MaxIter_inner   = 1000
        Rtol_inner      = 1.0d-8
        Include_NNS     = .TRUE.
        Include_NES     = .TRUE.
        Include_Pair    = .TRUE.
        Include_NuPair  = .TRUE.
        Include_Brem    = .TRUE.
        Include_LinCorr = .FALSE.
        wMatterRHS      = [ One, One, One, One, One, One ]
        DnuMax          = One-1.0d-6
        FreezeOpacities = .FALSE.

      ELSE

        nX  = [ 1, 1, 1 ]
        xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
        xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
        bcX = [ 0, 0, 0 ]

        TimeSteppingScheme = 'BackwardEuler'

        t_end = 1.0d2 * Millisecond

        PrescribedTimeStep = .TRUE.
        dt_0               = 1.0d-4 * Millisecond
        dt_MAX             = 1.0d+1 * Millisecond
        dt_RATE            = 1.04_DP
        iCycleD            = 1
        iCycleW            = 1
        maxCycles          = 100000

        EvolveEuler                    = .FALSE.
        UseSlopeLimiter_Euler          = .FALSE.
        UseSlopeLimiter_TwoMoment      = .FALSE.
        UsePositivityLimiter_Euler     = .FALSE.
        UsePositivityLimiter_TwoMoment = .FALSE.
        UseEnergyLimiter_TwoMoment     = .FALSE.
      
      ENDIF

    CASE( 'DeleptonizationWave1D' )

      CoordinateSystem = 'SPHERICAL'

      nSpecies = 6
      nNodes   = 2

      nX    = [ 256, 1, 1 ]
      xL    = [ 0.0_DP           , 0.0_DP, 0.0_DP ]
      xR    = [ 5.12d2 * Kilometer, Pi    , TwoPi  ]
      bcX   = [ 31, 1, 1 ]
      ZoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]

      nEquidistantX = 0
      dEquidistantX = 1.0_DP * Kilometer

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 11
      !ZoomE = 1.473936813452146d0
      ZoomE = 1.266038160710160_DP
      !ZoomE = 1.119237083677839

      TimeSteppingScheme = 'IMEX_PDARS'
      !TimeSteppingScheme = 'BackwardEuler'
      !TimeSteppingScheme = 'IMEX_ARS_111'

      t_end = 4.0d-4 * Millisecond

      PrescribedTimeStep = .true. ! If .FALSE., explicit CFL will be used.
      dt_0               = 1d-3 * Millisecond
      dt_MAX             = 1d-3 * Millisecond
      dt_RATE            = 2.0d0 !1.01_DP

      iCycleD = 1
      iCycleW = 100
      maxCycles = 1000000

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .TRUE.

      !ProfileName = 'input_thornado_VX_100ms.dat'
      ProfileName = 'input_thornado_20M_100ms.dat'

      wMatterRHS = [ One, One, One, Zero, Zero, Zero ] ! --- Keep Velocity Fixed
      !wMatterRHS = [ Zero, Zero, Zero, Zero, Zero, Zero ] ! --- Keep Fluid Fixed

      Include_NNS     = .TRUE.
      Include_NES     = .TRUE.
      Include_Pair    = .TRUE.
      Include_NuPair  = .TRUE.
      Include_Brem    = .TRUE.
      Include_LinCorr = .false.
      FreezeOpacities = .FALSE. ! --- Keep opacities fixed during iterations?

      DnuMax          = One - 1.0d-6

      M_inner = 2
      MaxIter_inner   = 1000

      Rtol_outer = 1.0d-8
      Rtol_inner = 1.0d-8
  
  
      UseConservativeProjection = .false.
      UseEquilibriumDeviationLLF = .false.

      WriteMatterDiagnostics = .FALSE.
      WriteMatterExchangeAttribution = .FALSE.
      WriteCumulativeConservation = .FALSE.
      !MatterDiagnosticsFileName = 'md_PDARS_All_Iso_lim_tol_1e-8_1e-8_no_cons2_wb_all.dat'
      MatterDiagnosticsFileName = 'md_nX256_nE16_nNodes2_diagnostic.dat'

    CASE( 'EquilibriumAdvection' )

      nSpecies = 2
      nNodes   = 2

      nX  = [ 8, 1, 1 ]
      xL  = [ - 5.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ + 5.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 1, 1, 1 ]

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 10
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end = 1.0d1 * Millisecond

      iCycleD = 1
      iCycleW = 100
      maxCycles = 100000

      EvolveEuler                    = .TRUE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .TRUE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

    CASE( 'MovingHomogeneousRelaxation' )

      nSpecies = 2
      MovingSpeciesString = ''
      CALL GET_ENVIRONMENT_VARIABLE &
             ( 'THORNADO_MOVING_NSPECIES', MovingSpeciesString, &
               LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
      IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
        READ(MovingSpeciesString,*,IOSTAT=ReadStatus) nSpecies
        IF( ReadStatus /= 0 .OR. .NOT. ANY( nSpecies == [ 2, 6 ] ) ) &
          ERROR STOP 'THORNADO_MOVING_NSPECIES must be 2 or 6.'
      END IF

      nNodes   = 2

      nX  = [ 4, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ 4.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 1, 1, 1 ]

      nE    = 16
      eL    = 0.0_DP * MeV
      eR    = 3.0d2 * MeV
      bcE   = 10
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'BackwardEuler'

      PrescribedTimeStep = .TRUE.
      dt_0               = 1.0d-4 * Millisecond
      dt_MAX             = dt_0
      dt_RATE            = One
      t_end              = 4.0_DP * dt_0

      iCycleD  = 1
      iCycleW  = 1
      maxCycles = 4

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment     = .FALSE.

      MovingAllProcesses = .FALSE.
      MovingAllProcessesString = ''
      CALL GET_ENVIRONMENT_VARIABLE &
             ( 'THORNADO_MOVING_ALL_PROCESSES', &
               MovingAllProcessesString, &
               LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
      IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
        SELECT CASE( TRIM( ADJUSTL( MovingAllProcessesString ) ) )
          CASE( 'true', 'TRUE', 'on', 'ON', '1' )
            MovingAllProcesses = .TRUE.
          CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
            MovingAllProcesses = .FALSE.
          CASE DEFAULT
            ERROR STOP &
              'THORNADO_MOVING_ALL_PROCESSES must be true or false.'
        END SELECT
      END IF

      Include_NNS     = MovingAllProcesses
      Include_NES     = MovingAllProcesses
      Include_Pair    = MovingAllProcesses
      Include_NuPair  = MovingAllProcesses
      Include_Brem    = MovingAllProcesses
      Include_LinCorr = .FALSE.

      wMatterRHS = [ One, One, One, One, One, One ]
      DnuMax = One - 1.0d-6

      WriteMatterDiagnostics = .FALSE.
      WriteMatterExchangeAttribution = .FALSE.
      WriteCumulativeConservation = .TRUE.
      MatterDiagnosticsFileName = 'moving_homogeneous_relaxation.dat'
      WriteSolverHistory = .TRUE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  ProfileNameString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_PROFILE_NAME', ProfileNameString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    IF( EnvironmentLength > LEN( ProfileName ) ) &
      ERROR STOP 'THORNADO_PROFILE_NAME is too long.'
    ProfileName = TRIM( ProfileNameString )
  END IF

  IF( TRIM( ProgramName ) == 'DeleptonizationWave1D' )THEN

    ProfileSpeciesString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_PROFILE_NSPECIES', ProfileSpeciesString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
      READ(ProfileSpeciesString,*,IOSTAT=ReadStatus) nSpecies
      IF( ReadStatus /= 0 .OR. .NOT. ANY( nSpecies == [ 2, 6 ] ) ) &
        ERROR STOP 'THORNADO_PROFILE_NSPECIES must be 2 or 6.'
    END IF

    ProfileAllProcesses = .FALSE.
    ProfileAllProcessesString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_PROFILE_ALL_PROCESSES', &
             ProfileAllProcessesString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
      SELECT CASE( TRIM( ADJUSTL( ProfileAllProcessesString ) ) )
        CASE( 'true', 'TRUE', 'on', 'ON', '1' )
          ProfileAllProcesses = .TRUE.
        CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
          ProfileAllProcesses = .FALSE.
        CASE DEFAULT
          ERROR STOP &
            'THORNADO_PROFILE_ALL_PROCESSES must be true or false.'
      END SELECT
    END IF
    IF( ProfileAllProcesses )THEN
      Include_NNS    = .TRUE.
      Include_NES    = .TRUE.
      Include_Pair   = .TRUE.
      Include_NuPair = .TRUE.
      Include_Brem   = .TRUE.
    END IF

    ProfileMomentumCoupling = .FALSE.
    ProfileMomentumCouplingString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_PROFILE_MOMENTUM_COUPLING', &
             ProfileMomentumCouplingString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
      SELECT CASE( TRIM( ADJUSTL( ProfileMomentumCouplingString ) ) )
        CASE( 'true', 'TRUE', 'on', 'ON', '1' )
          ProfileMomentumCoupling = .TRUE.
        CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
          ProfileMomentumCoupling = .FALSE.
        CASE DEFAULT
          ERROR STOP &
            'THORNADO_PROFILE_MOMENTUM_COUPLING must be true or false.'
      END SELECT
    END IF
    IF( ProfileMomentumCoupling ) wMatterRHS = One

    ProfileRadiationLimiters = .FALSE.
    ProfileRadiationLimitersString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_PROFILE_RADIATION_LIMITERS', &
             ProfileRadiationLimitersString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
      SELECT CASE( TRIM( ADJUSTL( ProfileRadiationLimitersString ) ) )
        CASE( 'true', 'TRUE', 'on', 'ON', '1' )
          ProfileRadiationLimiters = .TRUE.
        CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
          ProfileRadiationLimiters = .FALSE.
        CASE DEFAULT
          ERROR STOP &
            'THORNADO_PROFILE_RADIATION_LIMITERS must be true or false.'
      END SELECT
    END IF
    IF( ProfileRadiationLimiters )THEN
      UsePositivityLimiter_TwoMoment = .TRUE.
      UseEnergyLimiter_TwoMoment = .TRUE.
    END IF

  END IF

  WriteSolverHistoryString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_WRITE_SOLVER_HISTORY', WriteSolverHistoryString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    SELECT CASE( TRIM( ADJUSTL( WriteSolverHistoryString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        WriteSolverHistory = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        WriteSolverHistory = .FALSE.
      CASE DEFAULT
        ERROR STOP 'THORNADO_WRITE_SOLVER_HISTORY must be true or false.'
    END SELECT
  END IF

  WriteMatterDiagnosticsString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_WRITE_MATTER_DIAGNOSTICS', &
           WriteMatterDiagnosticsString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    SELECT CASE( TRIM( ADJUSTL( WriteMatterDiagnosticsString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        WriteMatterDiagnostics = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        WriteMatterDiagnostics = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_WRITE_MATTER_DIAGNOSTICS must be true or false.'
    END SELECT
  END IF

  WriteMatterExchangeAttributionString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_WRITE_MATTER_EXCHANGE_ATTRIBUTION', &
           WriteMatterExchangeAttributionString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    SELECT CASE( TRIM( ADJUSTL( WriteMatterExchangeAttributionString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        WriteMatterExchangeAttribution = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        WriteMatterExchangeAttribution = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_WRITE_MATTER_EXCHANGE_ATTRIBUTION must be true or false.'
    END SELECT
  END IF

  WriteCumulativeConservationString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_WRITE_CUMULATIVE_CONSERVATION', &
           WriteCumulativeConservationString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    SELECT CASE( TRIM( ADJUSTL( WriteCumulativeConservationString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        WriteCumulativeConservation = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        WriteCumulativeConservation = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_WRITE_CUMULATIVE_CONSERVATION must be true or false.'
    END SELECT
  END IF

  IF( TRIM( ProgramName ) == 'MovingHomogeneousRelaxation' )THEN

    MovingVelocityFraction = 0.1_DP
    MovingVelocityDirection = [ One, Zero, Zero ]
    MovingRadiationScale = 0.9_DP

    MovingVelocityString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_MOVING_VELOCITY_FRACTION', MovingVelocityString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
      READ(MovingVelocityString,*,IOSTAT=ReadStatus) MovingVelocityFraction
      IF( ReadStatus /= 0 ) &
        ERROR STOP 'Invalid THORNADO_MOVING_VELOCITY_FRACTION.'
    END IF

    MovingVelocityDirectionString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_MOVING_VELOCITY_DIRECTION', &
             MovingVelocityDirectionString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      SELECT CASE( TRIM( ADJUSTL( MovingVelocityDirectionString ) ) )
        CASE( 'x1', 'X1' )
          MovingVelocityDirection = [ One, Zero, Zero ]
        CASE( 'x2', 'X2' )
          MovingVelocityDirection = [ Zero, One, Zero ]
        CASE( 'x3', 'X3' )
          MovingVelocityDirection = [ Zero, Zero, One ]
        CASE( 'diagonal', 'DIAGONAL' )
          MovingVelocityDirection = [ One, One, One ]
        CASE DEFAULT
          ERROR STOP &
            'THORNADO_MOVING_VELOCITY_DIRECTION must be x1, x2, x3, or diagonal.'
      END SELECT

    END IF

    MovingRadiationScaleString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_MOVING_RADIATION_SCALE', &
             MovingRadiationScaleString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
      READ(MovingRadiationScaleString,*,IOSTAT=ReadStatus) MovingRadiationScale
      IF( ReadStatus /= 0 ) &
        ERROR STOP 'Invalid THORNADO_MOVING_RADIATION_SCALE.'
    END IF

    CALL SetMovingHomogeneousRelaxationParameters &
           ( MovingVelocityFraction, MovingRadiationScale, &
             MovingVelocityDirection )

  END IF

  AndersonDepthString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_M_INNER', AndersonDepthString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    READ(AndersonDepthString,*,IOSTAT=ReadStatus) M_inner
    IF( ReadStatus /= 0 .OR. M_inner < 1 ) &
      ERROR STOP 'THORNADO_M_INNER must be a positive integer.'
  END IF

  EndTimeMillisecondsString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_T_END_MS', EndTimeMillisecondsString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    READ(EndTimeMillisecondsString,*,IOSTAT=ReadStatus) t_end
    IF( ReadStatus /= 0 .OR. t_end <= Zero ) &
      ERROR STOP 'THORNADO_T_END_MS must be positive.'
    t_end = t_end * Millisecond
  END IF

  TimeStepMillisecondsString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_DT_MS', TimeStepMillisecondsString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )
  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    READ(TimeStepMillisecondsString,*,IOSTAT=ReadStatus) dt_0
    IF( ReadStatus /= 0 .OR. dt_0 <= Zero ) &
      ERROR STOP 'THORNADO_DT_MS must be positive.'
    PrescribedTimeStep = .TRUE.
    dt_0 = dt_0 * Millisecond
    dt_MAX = dt_0
    dt_RATE = One
  END IF

  RadiationOnlySolveString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_RADIATION_ONLY_SOLVE', RadiationOnlySolveString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    SELECT CASE( TRIM( ADJUSTL( RadiationOnlySolveString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        UseRadiationOnlySolve = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        UseRadiationOnlySolve = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_RADIATION_ONLY_SOLVE must be true or false.'
    END SELECT

  END IF

  MatterExchangeShadowsString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_MATTER_EXCHANGE_SHADOWS', &
           MatterExchangeShadowsString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    SELECT CASE( TRIM( ADJUSTL( MatterExchangeShadowsString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        ComputeMatterExchangeShadows = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        ComputeMatterExchangeShadows = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_MATTER_EXCHANGE_SHADOWS must be true or false.'
    END SELECT

  END IF

  NodesString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_NNODES', NodesString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    READ(NodesString,*,IOSTAT=ReadStatus) nNodes
    IF( ReadStatus /= 0 .OR. nNodes < 1 ) &
      ERROR STOP 'THORNADO_NNODES must be a positive integer.'

  END IF

  SpatialElementsString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_NX1', SpatialElementsString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    READ(SpatialElementsString,*,IOSTAT=ReadStatus) nX(1)
    IF( ReadStatus /= 0 .OR. nX(1) < 1 ) &
      ERROR STOP 'THORNADO_NX1 must be a positive integer.'

  END IF

  SourceAwareLLFString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_WELL_BALANCED_TRANSPORT', SourceAwareLLFString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus /= 0 .OR. EnvironmentLength == 0 ) &
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_SOURCE_AWARE_LLF', SourceAwareLLFString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    SELECT CASE( TRIM( ADJUSTL( SourceAwareLLFString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        UseSourceAwareLLF = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        UseSourceAwareLLF = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_SOURCE_AWARE_LLF must be true or false.'
    END SELECT

  END IF

  ConservativeProjectionString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_CONSERVATIVE_PROJECTION', &
           ConservativeProjectionString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    SELECT CASE( TRIM( ADJUSTL( ConservativeProjectionString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        UseConservativeProjection = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        UseConservativeProjection = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_CONSERVATIVE_PROJECTION must be true or false.'
    END SELECT

  END IF

  MatterDiagnosticsFileString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_MATTER_DIAGNOSTICS_FILE', &
           MatterDiagnosticsFileString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    IF( EnvironmentLength > LEN( MatterDiagnosticsFileName ) ) &
      ERROR STOP 'THORNADO_MATTER_DIAGNOSTICS_FILE is too long.'

    MatterDiagnosticsFileName = TRIM( MatterDiagnosticsFileString )

  END IF

  MaxCyclesString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_MAX_CYCLES', MaxCyclesString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN
    READ(MaxCyclesString,*,IOSTAT=ReadStatus) maxCycles
    IF( ReadStatus /= 0 .OR. maxCycles < 1 ) &
      ERROR STOP 'THORNADO_MAX_CYCLES must be a positive integer.'
  END IF

  APWellBalancedFluxSourceString = ''
  CALL GET_ENVIRONMENT_VARIABLE &
         ( 'THORNADO_AP_WELL_BALANCED_FLUX_SOURCE', &
           APWellBalancedFluxSourceString, &
           LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

  IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

    SELECT CASE( TRIM( ADJUSTL( APWellBalancedFluxSourceString ) ) )
      CASE( 'true', 'TRUE', 'on', 'ON', '1' )
        UseAPWellBalancedFluxSource = .TRUE.
      CASE( 'false', 'FALSE', 'off', 'OFF', '0' )
        UseAPWellBalancedFluxSource = .FALSE.
      CASE DEFAULT
        ERROR STOP &
          'THORNADO_AP_WELL_BALANCED_FLUX_SOURCE must be true or false.'
    END SELECT

  END IF

  IF( UseAPWellBalancedFluxSource )THEN

    IF( .NOT. RunCoupledAbsorptionDiffusion ) &
      ERROR STOP &
        'The AP/well-balanced flux-source test requires ' // &
        'CoupledAbsorptionDiffusion.'

    UseSourceAwareLLF = .TRUE.
    UseEquilibriumDeviationLLF = .TRUE.

  END IF

  IF( RunAbsorptionDiffusionBalance )THEN

    AbsorptionDiffusionOpticalDepthString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_ABSORPTION_DIFFUSION_OPTICAL_DEPTH', &
             AbsorptionDiffusionOpticalDepthString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      READ(AbsorptionDiffusionOpticalDepthString,*,IOSTAT=ReadStatus) &
        AbsorptionDiffusionMeanOpticalDepth
      IF( ReadStatus /= 0 ) &
        ERROR STOP 'Invalid THORNADO_ABSORPTION_DIFFUSION_OPTICAL_DEPTH.'

      CALL SetVariableOpacityDiffusionMeanOpticalDepth &
             ( AbsorptionDiffusionMeanOpticalDepth )

    END IF

    AbsorptionDiffusionLLFScalingString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_ABSORPTION_DIFFUSION_LLF_SCALING', &
             AbsorptionDiffusionLLFScalingString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      SELECT CASE( TRIM( ADJUSTL( AbsorptionDiffusionLLFScalingString ) ) )
        CASE( 'standard' )
          UseAbsorptionDiffusionScaledLLF = .FALSE.
        CASE( 'scaled' )
          UseAbsorptionDiffusionScaledLLF = .TRUE.
        CASE DEFAULT
          ERROR STOP &
            'THORNADO_ABSORPTION_DIFFUSION_LLF_SCALING must be standard or scaled.'
      END SELECT

    END IF

  END IF

  IF( RunScatteringDiffusionLimit )THEN

    AbsorptionDiffusionOpticalDepthString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_SCATTERING_DIFFUSION_OPTICAL_DEPTH', &
             AbsorptionDiffusionOpticalDepthString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      READ(AbsorptionDiffusionOpticalDepthString,*,IOSTAT=ReadStatus) &
        AbsorptionDiffusionMeanOpticalDepth
      IF( ReadStatus /= 0 ) &
        ERROR STOP 'Invalid THORNADO_SCATTERING_DIFFUSION_OPTICAL_DEPTH.'

      CALL SetVariableOpacityDiffusionMeanOpticalDepth &
             ( AbsorptionDiffusionMeanOpticalDepth )

    END IF

    ScatteringDiffusionNXString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_SCATTERING_DIFFUSION_NX', &
             ScatteringDiffusionNXString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      READ(ScatteringDiffusionNXString,*,IOSTAT=ReadStatus) &
        ScatteringDiffusionNX
      IF( ReadStatus /= 0 .OR. ScatteringDiffusionNX < 2 ) &
        ERROR STOP 'THORNADO_SCATTERING_DIFFUSION_NX must be at least 2.'

      nX(1) = ScatteringDiffusionNX

    END IF

    t_end = ScatteringDiffusionEfoldings * Three &
      * VariableOpacityDiffusionMeanOpacity &
      * OpaqueDiffusionDomainLength**2 / TwoPi**2

    WRITE(*,* )
    WRITE(*,'(A,ES14.6E3)') &
      'Scattering-diffusion mean optical depth = ', &
      VariableOpacityDiffusionMeanOpacity * OpaqueDiffusionDomainLength
    WRITE(*,'(A,ES14.6E3)') &
      '  requested nominal diffusion e-foldings= ', &
      ScatteringDiffusionEfoldings
    WRITE(*,'(A,ES14.6E3)') &
      '  final time [light-crossing times]      = ', &
      t_end / OpaqueDiffusionDomainLength
    WRITE(*,'(A,I8)') &
      '  spatial elements                       = ', nX(1)
    WRITE(*,* )

  END IF

  IF( RunCoupledAbsorptionDiffusion )THEN

    AbsorptionDiffusionOpticalDepthString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_COUPLED_DIFFUSION_OPTICAL_DEPTH', &
             AbsorptionDiffusionOpticalDepthString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      READ(AbsorptionDiffusionOpticalDepthString,*,IOSTAT=ReadStatus) &
        AbsorptionDiffusionMeanOpticalDepth
      IF( ReadStatus /= 0 .OR. AbsorptionDiffusionMeanOpticalDepth <= Zero ) &
        ERROR STOP 'Invalid THORNADO_COUPLED_DIFFUSION_OPTICAL_DEPTH.'

      CALL SetVariableOpacityDiffusionMeanOpticalDepth &
             ( AbsorptionDiffusionMeanOpticalDepth )

    END IF

    ScatteringDiffusionNXString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_COUPLED_DIFFUSION_NX', &
             ScatteringDiffusionNXString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      READ(ScatteringDiffusionNXString,*,IOSTAT=ReadStatus) &
        ScatteringDiffusionNX
      IF( ReadStatus /= 0 .OR. ScatteringDiffusionNX < 2 ) &
        ERROR STOP 'THORNADO_COUPLED_DIFFUSION_NX must be at least 2.'

      nX(1) = ScatteringDiffusionNX

    END IF

    CoupledDiffusionSchemeString = ''
    CALL GET_ENVIRONMENT_VARIABLE &
           ( 'THORNADO_COUPLED_DIFFUSION_SCHEME', &
             CoupledDiffusionSchemeString, &
             LENGTH = EnvironmentLength, STATUS = EnvironmentStatus )

    IF( EnvironmentStatus == 0 .AND. EnvironmentLength > 0 )THEN

      SELECT CASE( TRIM( ADJUSTL( CoupledDiffusionSchemeString ) ) )
        CASE( 'IMEX_ARS_111', 'IMEX_PDARS' )
          TimeSteppingScheme &
            = TRIM( ADJUSTL( CoupledDiffusionSchemeString ) )
        CASE DEFAULT
          ERROR STOP &
            'THORNADO_COUPLED_DIFFUSION_SCHEME must be IMEX_ARS_111 or IMEX_PDARS.'
      END SELECT

    END IF

    t_end = CoupledDiffusionEfoldings &
      / ABS( CoupledDiffusionDecayRate() )

    WRITE(*,* )
    WRITE(*,'(A,ES14.6E3)') &
      'Coupled absorption-diffusion optical depth = ', &
      VariableOpacityDiffusionMeanOpacity * OpaqueDiffusionDomainLength
    WRITE(*,'(A,A)') '  time integrator                        = ', &
      TRIM( TimeSteppingScheme )
    WRITE(*,'(A,ES14.6E3)') &
      '  requested slow-mode e-foldings         = ', &
      CoupledDiffusionEfoldings
    WRITE(*,'(A,ES14.6E3)') &
      '  final time [light-crossing times]       = ', &
      t_end / OpaqueDiffusionDomainLength
    WRITE(*,'(A,I8)') &
      '  spatial elements                        = ', nX(1)
    WRITE(*,* )

  END IF

  CALL InitializeDriver

  IF( RunAbsorptionDiffusionBalance )THEN

    IF( UseSourceAwareLLF .AND. UseAbsorptionDiffusionScaledLLF ) &
      ERROR STOP &
        'Select either source-aware or globally scaled LLF, not both.'

    AbsorptionDiffusionMeanCellWidth &
      = SUM( MeshX(1) % Width(iX_B0(1):iX_E0(1)) ) &
          / REAL( iX_E0(1) - iX_B0(1) + 1, DP )
    AbsorptionDiffusionMeanCellOpticalDepth &
      = VariableOpacityDiffusionMeanOpacity &
          * AbsorptionDiffusionMeanCellWidth

    AbsorptionDiffusionLLFPenaltyFactor = One
    IF( UseAbsorptionDiffusionScaledLLF ) &
      AbsorptionDiffusionLLFPenaltyFactor &
        = MIN( One, One / AbsorptionDiffusionMeanCellOpticalDepth )

    WRITE(*,*)
    WRITE(*,'(A,L1)') &
      'Absorption-diffusion scaled LLF enabled = ', &
      UseAbsorptionDiffusionScaledLLF
    WRITE(*,'(A,ES14.6E3)') &
      '  mean cell optical depth              = ', &
      AbsorptionDiffusionMeanCellOpticalDepth
    WRITE(*,'(A,ES14.6E3)') &
      '  LLF dissipative factor               = ', &
      AbsorptionDiffusionLLFPenaltyFactor
    WRITE(*,*)

  END IF

  IF(      RunEquilibriumTest &
      .OR. RunManufacturedOpaqueDiffusion &
      .OR. RunManufacturedVariableOpacityDiffusion &
      .OR. RunScatteringDiffusionLimit &
      .OR. RunAbsorptionDiffusionBalance &
      .OR. RunCoupledAbsorptionDiffusion )THEN
    CALL InitializeFields( ProfileName, DnuMax_Option = DnuMax )
  ELSE IF (Relaxation_restart_from_file) THEN
    CALL InitializeFields( ProfileName, 'Relaxation_input.dat' )
  ELSE
    CALL InitializeFields( ProfileName )
  ENDIF


  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

  IF( RunEquilibriumTest )THEN
    ALLOCATE( EquilibriumTest_CF_Initial, MOLD = uCF )
    ALLOCATE( EquilibriumTest_CR_Initial, MOLD = uCR )
    EquilibriumTest_CF_Initial = uCF
    EquilibriumTest_CR_Initial = uCR
  END IF

  IF( RestartFileNumber .LT. 0 )THEN

    t = Zero

    ! --- Apply Slope Limiter to Initial Data ---

    CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    CALL ApplySlopeLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

    ! --- Apply Positivity Limiter to Initial Data ---

    CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

    IF( TRIM( ProgramName ) == 'VaryingLocalEquilibrium' ) &
      CALL ReportVaryingEquilibriumLimiterInput

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

    CALL WriteFieldsHDF &
           ( Time = 0.0_DP, &
             WriteGF_Option = .TRUE., &
             WriteFF_Option = .TRUE., &
             WriteRF_Option = .TRUE. )

    CALL ComputeTally &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR, &
             SetInitialValues_Option = .TRUE. )

    IF( RunEquilibriumTest )THEN

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( uCF, uCR )
#endif

      CALL BeginEquilibriumTest( dt_0, RequireStationaryEquilibrium )

    END IF

    IF( RunAbsorptionDiffusionBalance ) &
      CALL ReportAbsorptionDiffusionOperatorBalance

  ELSE

    CALL ReadFieldsHDF &
           ( RestartFileNumber, t, &
             ReadGF_Option = .TRUE., &
             ReadFF_Option = .TRUE., &
             ReadRF_Option = .TRUE. )

  END IF

  IF( WriteMatterDiagnostics )THEN
    ALLOCATE( CF_BeforeMatterDiagnosticsStep, MOLD = uCF )
    ALLOCATE( CR_BeforeMatterDiagnosticsStep, MOLD = uCR )
    ALLOCATE &
      ( MatterCollisionExchange_Rate_LastStep &
          (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
           iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), iX_B1(3):iX_E1(3)) )
    ALLOCATE &
      ( MatterUpdate_Rate_LastStep &
          (1:nQ,1:nDOFX, &
           iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), iX_B1(3):iX_E1(3)) )
    ALLOCATE &
      ( MatterCollisionRHS_FinalStage &
          (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
           iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), iX_B1(3):iX_E1(3)) )
    ALLOCATE( RadiationState_BeforeFinalLimiter, MOLD = uCR )
    CF_BeforeMatterDiagnosticsStep = uCF
    CR_BeforeMatterDiagnosticsStep = uCR
    MatterCollisionExchange_Rate_LastStep     = Zero
    MatterCollisionRHS_FinalStage   = Zero
    MatterUpdate_Rate_LastStep  = Zero
    RadiationState_BeforeFinalLimiter = Zero
    MatterDiagnosticsTimeStep = Zero
  END IF

  IF(      WriteMatterDiagnostics &
      .OR. WriteCumulativeConservation &
      .OR. ( WriteMatterExchangeAttribution &
             .AND. ComputeMatterExchangeShadows ) )THEN
    nIMEXStages = GetNumberOfStages_IMEX_RK()
    ALLOCATE &
      ( IMEX_Stage_FluidState &
          (1:nDOFX, &
           iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), iX_B1(3):iX_E1(3), &
           1:SIZE(uCF,5),1:nIMEXStages,1:nIMEXStageStates) )
    ALLOCATE &
      ( IMEX_Stage_RadiationState &
          (1:SIZE(uCR,1), &
           iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
           1:SIZE(uCR,6),1:SIZE(uCR,7), &
           1:nIMEXStages,1:nIMEXStageStates) )
    IMEX_Stage_FluidState           = Zero
    IMEX_Stage_RadiationState       = Zero
  END IF

  IF(      WriteMatterDiagnostics &
      .OR. WriteCumulativeConservation &
      .OR. ( WriteMatterExchangeAttribution &
             .AND. ComputeMatterExchangeShadows ) )THEN
    ALLOCATE &
      ( IMEX_Stage_FluidImplicitRate &
          (1:nDOFX, &
           iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), iX_B1(3):iX_E1(3), &
           1:SIZE(uCF,5),1:nIMEXStages) )
    ALLOCATE( IMEX_Stage_FluidExplicitRate, &
              MOLD = IMEX_Stage_FluidImplicitRate )
    ALLOCATE &
      ( IMEX_Stage_RadiationImplicitRate &
          (1:SIZE(uCR,1), &
           iZ_B1(1):iZ_E1(1), iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), iZ_B1(4):iZ_E1(4), &
           1:SIZE(uCR,6),1:SIZE(uCR,7),1:nIMEXStages) )
    ALLOCATE( IMEX_Stage_RadiationExplicitRate, &
              MOLD = IMEX_Stage_RadiationImplicitRate )
    ALLOCATE &
      ( IMEX_Stage_LeptonFlavorCheckRate &
          (1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
           iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3), &
           1:nIMEXStages) )
    IMEX_Stage_FluidImplicitRate     = Zero
    IMEX_Stage_FluidExplicitRate     = Zero
    IMEX_Stage_RadiationImplicitRate = Zero
    IMEX_Stage_RadiationExplicitRate = Zero
    IMEX_Stage_LeptonFlavorCheckRate = Zero
  END IF

  IF( WriteMatterExchangeAttribution )THEN
    ALLOCATE( CF_AttributionStart, MOLD = uCF )
    IF( ComputeMatterExchangeShadows )THEN
      ALLOCATE( CF_AttributionNoStream, MOLD = uCF )
      ALLOCATE( CF_AttributionAP, MOLD = uCF )
      ALLOCATE( CF_AttributionStandard, MOLD = uCF )
      ALLOCATE( CR_AttributionNoStream, MOLD = uCR )
      ALLOCATE( CR_AttributionAP, MOLD = uCR )
      ALLOCATE( CR_AttributionStandard, MOLD = uCR )
    END IF
    CALL InitializeMatterExchangeAttribution &
           ( TRIM(MatterDiagnosticsFileName) // '.exchange_attribution', &
             iZ_B0, iZ_E0, &
             ComputeShadows_Option = ComputeMatterExchangeShadows )
  END IF

  IF( WriteCumulativeConservation ) &
    CALL InitializeCumulativeConservationDiagnostics &
           ( TRIM(MatterDiagnosticsFileName) // '.cumulative_conservation', &
             t, iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  IF( WriteSolverHistory ) &
    CALL InitializeSolverDiagnostics &
           ( TRIM(MatterDiagnosticsFileName) // '.solver_history' )

  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end .AND. iCycle < maxCycles )

    iCycle = iCycle + 1

    IF( PrescribedTimeStep )THEN

      dt = MIN( ( dt_RATE )**( iCycle - 1 ) * dt_0, dt_MAX )

    ELSE

      CALL ComputeTimeStep_TwoMoment &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, &
               0.3_DP/(Two*DBLE(nNodes-1)+One), dt_CFL )

      dt = dt_CFL

    END IF

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( dt <= 100.0_DP * EPSILON(One) &
                  * MAX( ABS(t_end), TINY(One) ) )THEN

      t = t_end
      EXIT

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = '    ,  t / Millisecond, &
          '', 'dt = '   , dt / Millisecond

    END IF

    CaptureMatterDiagnosticsStep = WriteMatterDiagnostics .AND. &
      (      t + dt >= t_end &
        - 100.0_DP * EPSILON(One) * MAX( ABS(t_end), TINY(One) ) &
        .OR. iCycle == maxCycles )

    CaptureSolverConvergenceStep = CaptureMatterDiagnosticsStep .OR. &
      WriteSolverHistory

    CALL SetNeutrinoMatterSolverConvergenceDiagnostics &
           ( CaptureSolverConvergenceStep )

    IF( CaptureMatterDiagnosticsStep )THEN
      CF_BeforeMatterDiagnosticsStep = uCF
      CR_BeforeMatterDiagnosticsStep = uCR
      MatterDiagnosticsTimeStep = dt
      MatterCollisionExchange_Rate_LastStep       = Zero
      MatterUpdate_Rate_LastStep    = Zero
      MatterCollisionRHS_FinalStage     = Zero
      RadiationState_BeforeFinalLimiter = Zero
    END IF

    IF(      CaptureMatterDiagnosticsStep &
        .OR. ( WriteMatterExchangeAttribution &
               .AND. ComputeMatterExchangeShadows ) )THEN
      IMEX_Stage_FluidState     = Zero
      IMEX_Stage_RadiationState = Zero
    END IF

    IF( WriteCumulativeConservation )THEN
      IMEX_Stage_FluidState             = Zero
      IMEX_Stage_RadiationState         = Zero
      IMEX_Stage_FluidImplicitRate      = Zero
      IMEX_Stage_FluidExplicitRate      = Zero
      IMEX_Stage_RadiationImplicitRate  = Zero
      IMEX_Stage_RadiationExplicitRate  = Zero
      IMEX_Stage_LeptonFlavorCheckRate  = Zero
    END IF

    IF( WriteMatterExchangeAttribution )THEN

      CF_AttributionStart = uCF

    END IF

    IF( WriteMatterExchangeAttribution .AND. ComputeMatterExchangeShadows )THEN

      CF_AttributionNoStream  = uCF
      CF_AttributionAP        = uCF
      CF_AttributionStandard  = uCF
      CR_AttributionNoStream  = uCR
      CR_AttributionAP        = uCR
      CR_AttributionStandard  = uCR

      CALL InitializeSourceAwareLLF( Use_Option = .FALSE. )
      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, &
               CF_AttributionNoStream, CR_AttributionNoStream, &
               ComputeIncrement_TwoMoment_Implicit, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_TwoMoment_Explicit_Zero, &
               UpdateTallies_Option = .FALSE. )

      IF( UseSourceAwareLLF )THEN

        CALL Update_IMEX_RK &
               ( dt, uGE, uGF, &
                 CF_AttributionStandard, CR_AttributionStandard, &
                 ComputeIncrement_TwoMoment_Implicit, &
                 ComputeIncrement_TwoMoment_Explicit_Option &
                   = ComputeIncrement_TwoMoment_Explicit_AP, &
                 UpdateTallies_Option = .FALSE. )

      ELSE

        CALL InitializeSourceAwareLLF( Use_Option = .TRUE. )
        CALL Update_IMEX_RK &
               ( dt, uGE, uGF, CF_AttributionAP, CR_AttributionAP, &
                 ComputeIncrement_TwoMoment_Implicit, &
                 ComputeIncrement_TwoMoment_Explicit_Option &
                   = ComputeIncrement_TwoMoment_Explicit_AP, &
                 UpdateTallies_Option = .FALSE. )

      END IF

      CALL InitializeSourceAwareLLF( Use_Option = UseSourceAwareLLF )

    END IF

    OffGridFluidIncrement = Zero
    OffGridRadiationIncrement = Zero
    OffGridRadiationLeptonFlavorIncrement = Zero

    IF( RunCoupledAbsorptionDiffusion )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, &
               ComputeIncrement_CoupledAbsorptionDiffusion, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_CoupledAbsorptionDiffusion_Explicit, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE IF( RunAbsorptionDiffusionBalance )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, &
               ComputeIncrement_AbsorptionDiffusionBalance, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_AbsorptionDiffusionBalance_Explicit, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE IF(      RunManufacturedVariableOpacityDiffusion &
             .OR. RunScatteringDiffusionLimit )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, &
               ComputeIncrement_ManufacturedVariableOpacityDiffusion, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_ManufacturedVariableOpacityDiffusion_Explicit, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE IF( RunManufacturedOpaqueDiffusion )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, &
               ComputeIncrement_ManufacturedOpaqueDiffusion, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_ManufacturedOpaqueDiffusion_Explicit, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE IF( CaptureMatterDiagnosticsStep )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_TwoMoment_Explicit_AP, &
               MatterCollisionExchange_Rate_Option = MatterCollisionExchange_Rate_LastStep, &
               MatterUpdate_Rate_Option = MatterUpdate_Rate_LastStep, &
               MatterCollisionRHS_FinalStage_Option &
                 = MatterCollisionRHS_FinalStage, &
               RadiationState_BeforeFinalLimiter_Option &
                 = RadiationState_BeforeFinalLimiter, &
               IMEX_Stage_FluidState_Option &
                 = IMEX_Stage_FluidState, &
               IMEX_Stage_RadiationState_Option &
                 = IMEX_Stage_RadiationState, &
               IMEX_Stage_FluidImplicitRate_Option &
                 = IMEX_Stage_FluidImplicitRate, &
               IMEX_Stage_FluidExplicitRate_Option &
                 = IMEX_Stage_FluidExplicitRate, &
               IMEX_Stage_RadiationImplicitRate_Option &
                 = IMEX_Stage_RadiationImplicitRate, &
               IMEX_Stage_RadiationExplicitRate_Option &
                 = IMEX_Stage_RadiationExplicitRate, &
               IMEX_Stage_LeptonFlavorCheck_Rate_Option &
                 = IMEX_Stage_LeptonFlavorCheckRate, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE IF(      WriteMatterExchangeAttribution &
             .AND. ComputeMatterExchangeShadows )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_TwoMoment_Explicit_AP, &
               IMEX_Stage_FluidState_Option &
                 = IMEX_Stage_FluidState, &
               IMEX_Stage_RadiationState_Option &
                 = IMEX_Stage_RadiationState, &
               IMEX_Stage_FluidImplicitRate_Option &
                 = IMEX_Stage_FluidImplicitRate, &
               IMEX_Stage_FluidExplicitRate_Option &
                 = IMEX_Stage_FluidExplicitRate, &
               IMEX_Stage_RadiationImplicitRate_Option &
                 = IMEX_Stage_RadiationImplicitRate, &
               IMEX_Stage_RadiationExplicitRate_Option &
                 = IMEX_Stage_RadiationExplicitRate, &
               IMEX_Stage_LeptonFlavorCheck_Rate_Option &
                 = IMEX_Stage_LeptonFlavorCheckRate, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE IF( WriteCumulativeConservation )THEN

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_TwoMoment_Explicit_AP, &
               IMEX_Stage_FluidState_Option &
                 = IMEX_Stage_FluidState, &
               IMEX_Stage_RadiationState_Option &
                 = IMEX_Stage_RadiationState, &
               IMEX_Stage_FluidImplicitRate_Option &
                 = IMEX_Stage_FluidImplicitRate, &
               IMEX_Stage_FluidExplicitRate_Option &
                 = IMEX_Stage_FluidExplicitRate, &
               IMEX_Stage_RadiationImplicitRate_Option &
                 = IMEX_Stage_RadiationImplicitRate, &
               IMEX_Stage_RadiationExplicitRate_Option &
                 = IMEX_Stage_RadiationExplicitRate, &
               IMEX_Stage_LeptonFlavorCheck_Rate_Option &
                 = IMEX_Stage_LeptonFlavorCheckRate, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    ELSE

      CALL Update_IMEX_RK &
             ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit, &
               ComputeIncrement_TwoMoment_Explicit_Option &
                 = ComputeIncrement_TwoMoment_Explicit_AP, &
               OffGridFluidIncrement_Option = OffGridFluidIncrement, &
               OffGridRadiationIncrement_Option = OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement_Option &
                 = OffGridRadiationLeptonFlavorIncrement )

    END IF

    IF( WriteCumulativeConservation ) &
      CALL AccumulateCumulativeConservationDiagnostics &
             ( iCycle, t + dt, dt, iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
               uGE, uGF, uCF, uCR, &
               OffGridFluidIncrement, OffGridRadiationIncrement, &
               OffGridRadiationLeptonFlavorIncrement, &
               IMEX_Stage_FluidState, IMEX_Stage_RadiationState, &
               IMEX_Stage_FluidImplicitRate, &
               IMEX_Stage_FluidExplicitRate, &
               IMEX_Stage_RadiationImplicitRate, &
               IMEX_Stage_RadiationExplicitRate, &
               IMEX_Stage_LeptonFlavorCheckRate )

    IF( WriteMatterExchangeAttribution )THEN

      IF( ComputeMatterExchangeShadows )THEN

        IF( UseSourceAwareLLF )THEN
          CF_AttributionAP = uCF
          CR_AttributionAP = uCR
        ELSE
          CF_AttributionStandard = uCF
          CR_AttributionStandard = uCR
        END IF

        CALL AccumulateMatterExchangeAttribution &
               ( t + dt, dt, DnuMax, &
                 iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, &
                 CF_AttributionStart, uCF, &
                 CF_AttributionNoStream, &
                 CF_AttributionAP, CF_AttributionStandard, &
                 IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

      ELSE

        CALL AccumulateActualMatterExchangeAttribution &
               ( t + dt, dt, &
                 iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, &
                 CF_AttributionStart, uCF )

      END IF

    END IF

    IF( WriteSolverHistory ) &
      CALL WriteSolverDiagnostics &
             ( TRIM(MatterDiagnosticsFileName) // '.solver_history', &
               iCycle, t + dt, dt )

    CALL SetNeutrinoMatterSolverConvergenceDiagnostics( .FALSE. )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL ComputeFromConserved_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

      CALL ComputeTally &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

    END IF

  END DO

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

  IF( RunEquilibriumTest )THEN

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( uCF, uCR )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( uCF, uCR )
#endif

    CALL EndEquilibriumTest( RequireStationaryEquilibrium )

  END IF

  IF( WriteMatterDiagnostics )THEN
    IF( MatterDiagnosticsTimeStep <= Zero )THEN
      ERROR STOP 'No completed timestep is available for matter diagnostics.'
    END IF

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( uDR )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( uDR )
#endif

    CALL WriteStandaloneMatterDiagnostics &
           ( MatterDiagnosticsFileName, CoordinateSystem, t, &
             MatterDiagnosticsTimeStep, DnuMax, &
             iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, &
             CF_BeforeMatterDiagnosticsStep, CR_BeforeMatterDiagnosticsStep, &
             uCF, &
             uPF, uAF, uCR, uPR, uGR, uDR, &
             MatterCollisionExchange_Rate_LastStep, MatterUpdate_Rate_LastStep, &
             MatterCollisionRHS_FinalStage, &
             RadiationState_BeforeFinalLimiter, &
             IMEX_Stage_FluidState, IMEX_Stage_RadiationState )

    DEALLOCATE &
      ( CF_BeforeMatterDiagnosticsStep, CR_BeforeMatterDiagnosticsStep, &
        MatterCollisionExchange_Rate_LastStep, MatterUpdate_Rate_LastStep, &
        MatterCollisionRHS_FinalStage, &
        RadiationState_BeforeFinalLimiter )
  END IF

  IF( WriteMatterExchangeAttribution )THEN
    CALL FinalizeMatterExchangeAttribution
    DEALLOCATE( CF_AttributionStart )
    IF( ComputeMatterExchangeShadows )THEN
      DEALLOCATE &
        ( CF_AttributionNoStream, CF_AttributionAP, &
          CF_AttributionStandard, CR_AttributionNoStream, &
          CR_AttributionAP, CR_AttributionStandard )
    END IF
  END IF

  IF( WriteCumulativeConservation ) &
    CALL FinalizeCumulativeConservationDiagnostics

  IF(      WriteMatterDiagnostics &
      .OR. WriteCumulativeConservation &
      .OR. ( WriteMatterExchangeAttribution &
             .AND. ComputeMatterExchangeShadows ) ) &
    DEALLOCATE &
      ( IMEX_Stage_FluidState, IMEX_Stage_RadiationState, &
        IMEX_Stage_FluidImplicitRate, IMEX_Stage_FluidExplicitRate, &
        IMEX_Stage_RadiationImplicitRate, &
        IMEX_Stage_RadiationExplicitRate, &
        IMEX_Stage_LeptonFlavorCheckRate )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL ComputeTally &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, t, uGE, uGF, uCF, uCR )

  CALL ComputeError( t )

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeSolverDiagnostics( FileName )

    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: Unit

    OPEN( NEWUNIT = Unit, FILE = TRIM(FileName), STATUS = 'REPLACE', &
          ACTION = 'WRITE' )
    WRITE(Unit,'(A)') &
      '# cycle time_ms dt_ms point inner_iterations outer_iterations ' // &
      'first_relative_residual accepted_relative_residual ' // &
      'max_matter_residual_rate max_initial_projection_rate ' // &
      'max_conservation_handoff_rate'
    CLOSE(Unit)

  END SUBROUTINE InitializeSolverDiagnostics


  SUBROUTINE WriteSolverDiagnostics &
    ( FileName, Cycle, Time, TimeStep )

    CHARACTER(*), INTENT(in) :: FileName
    INTEGER,      INTENT(in) :: Cycle
    REAL(DP),     INTENT(in) :: Time, TimeStep

    INTEGER :: i, Unit
    INTEGER,  ALLOCATABLE :: InnerIterations(:), OuterIterations(:)
    INTEGER,  ALLOCATABLE :: ProjectionStatus(:)
    REAL(DP), ALLOCATABLE :: FirstResidual(:), AcceptedResidual(:)
    REAL(DP), ALLOCATABLE :: MatterResidualRate(:,:)
    REAL(DP), ALLOCATABLE :: InitialProjectionRate(:,:)
    REAL(DP), ALLOCATABLE :: ProjectionResidual(:)
    REAL(DP), ALLOCATABLE :: ProjectionHeatingRate(:)
    REAL(DP), ALLOCATABLE :: ConservationHandoffCheckRate(:,:)

    CALL GetNeutrinoMatterSolverConvergenceDiagnostics &
           ( InnerIterations, OuterIterations, FirstResidual, &
             MatterResidualRate, InitialProjectionRate, AcceptedResidual, &
             ProjectionResidual, ProjectionHeatingRate, ProjectionStatus, &
             ConservationHandoffCheckRate )

    OPEN( NEWUNIT = Unit, FILE = TRIM(FileName), STATUS = 'OLD', &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    DO i = 1, SIZE( InnerIterations )
      WRITE(Unit,'(I8,1X,2(ES24.16E3,1X),3(I8,1X),5(ES24.16E3,1X))') &
        Cycle, Time / Millisecond, TimeStep / Millisecond, &
        i, InnerIterations(i), OuterIterations(i), &
        FirstResidual(i), AcceptedResidual(i), &
        MAXVAL( ABS( MatterResidualRate(:,i) ) ), &
        MAXVAL( ABS( InitialProjectionRate(:,i) ) ), &
        MAXVAL( ABS( ConservationHandoffCheckRate(:,i) ) )
    END DO

    CLOSE(Unit)

    DEALLOCATE &
      ( InnerIterations, OuterIterations, ProjectionStatus, &
        FirstResidual, AcceptedResidual, MatterResidualRate, &
        InitialProjectionRate, ProjectionResidual, &
        ProjectionHeatingRate, ConservationHandoffCheckRate )

  END SUBROUTINE WriteSolverDiagnostics


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit_Zero &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: LLFPenaltyFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: SurfaceTermFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: VolumeTermFactor_Option
    LOGICAL, INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    dU_R = Zero

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit_Zero


  SUBROUTINE SetAnalyticSourceAwareTransportOpacity &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      UseVariableOpacity, &
      ConstantOpacity_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    LOGICAL, INTENT(in) :: UseVariableOpacity
    REAL(DP), INTENT(in), OPTIONAL :: ConstantOpacity_Option

    INTEGER :: iN_E, iN_X, iS
    INTEGER :: iNodeX, iNodeX1, iX1
    INTEGER :: iZ_B_Opacity(4), iZ_E_Opacity(4), nX_Opacity(3)
    INTEGER :: nE_Opacity, nE_G_Opacity, nX_G_Opacity
    REAL(DP) :: X, Opacity
    REAL(DP), ALLOCATABLE :: TransportOpacity(:,:,:)

    IF( .NOT. UseSourceAwareLLF ) RETURN

    IF( .NOT. UseVariableOpacity .AND. &
        .NOT. PRESENT( ConstantOpacity_Option ) ) &
      ERROR STOP 'A constant source-aware transport opacity is required.'

    iZ_B_Opacity(1)   = iZ_B0_Input(1)
    iZ_E_Opacity(1)   = iZ_E0_Input(1)
    iZ_B_Opacity(2:4) = iZ_B1_Input(2:4)
    iZ_E_Opacity(2:4) = iZ_E1_Input(2:4)

    nE_Opacity   = iZ_E_Opacity(1) - iZ_B_Opacity(1) + 1
    nX_Opacity   = iZ_E_Opacity(2:4) - iZ_B_Opacity(2:4) + 1
    nE_G_Opacity = nDOFE * nE_Opacity
    nX_G_Opacity = nDOFX * PRODUCT( nX_Opacity )

    ALLOCATE( TransportOpacity(nE_G_Opacity,nSpecies,nX_G_Opacity) )

    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G_Opacity

      iX1     = MOD( (iN_X-1) / nDOFX, nX_Opacity(1) ) &
                  + iZ_B_Opacity(2)
      iNodeX  = MOD(  iN_X-1,          nDOFX ) + 1

      IF( UseVariableOpacity )THEN
        iNodeX1 = NodeNumberTableX(1,iNodeX)
        X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        Opacity = VariableOpacityDiffusionOpacity( X )
      ELSE
        Opacity = ConstantOpacity_Option
      END IF

      DO iN_E = 1, nE_G_Opacity
        TransportOpacity(iN_E,iS,iN_X) = Opacity
      END DO

    END DO
    END DO

    CALL SetSourceAwareTransportOpacity_Flattened &
           ( iZ_B_Opacity, iZ_E_Opacity, nX_Opacity, TransportOpacity )

    DEALLOCATE( TransportOpacity )

  END SUBROUTINE SetAnalyticSourceAwareTransportOpacity


  SUBROUTINE ComputeIncrement_CoupledAbsorptionDiffusion_Explicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: LLFPenaltyFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: SurfaceTermFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: VolumeTermFactor_Option
    LOGICAL, INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    REAL(DP), ALLOCATABLE :: EquilibriumReference(:,:,:,:,:,:,:)

    IF( UseAPWellBalancedFluxSource )THEN

      ALLOCATE( EquilibriumReference, MOLD = U_R )
      EquilibriumReference = U_R

      ! This is the same stage-local equilibrium used by the implicit
      ! radiation--reservoir exchange below.
      EquilibriumReference(:,:,:,:,:,iCR_N,1) &
        = CoupledDiffusionEquilibriumSlope * U_R(:,:,:,:,:,iCR_N,2)
      EquilibriumReference(:,:,:,:,:,iCR_G1:iCR_G3,1) = Zero

      CALL SetEquilibriumDeviationReferenceState &
             ( iZ_B1_Input, iZ_E1_Input, EquilibriumReference )

    END IF

    CALL SetAnalyticSourceAwareTransportOpacity &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, .FALSE., &
             ConstantOpacity_Option = VariableOpacityDiffusionMeanOpacity )

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             GE, GX, U_F, U_R, dU_R, &
             LLFPenaltyFactor_Option = LLFPenaltyFactor_Option, &
             SurfaceTermFactor_Option = SurfaceTermFactor_Option, &
             VolumeTermFactor_Option = VolumeTermFactor_Option, &
             IncludeObserverCorrections_Option &
               = IncludeObserverCorrections_Option )

    IF( UseSourceAwareLLF ) CALL ClearSourceAwareTransportOpacity

    IF( UseAPWellBalancedFluxSource )THEN

      CALL ClearEquilibriumDeviationReferenceState
      DEALLOCATE( EquilibriumReference )

    END IF

    ! Species two is the non-streaming finite-capacity matter reservoir.
    dU_R(:,:,:,:,:,:,2) = Zero

  END SUBROUTINE ComputeIncrement_CoupledAbsorptionDiffusion_Explicit


  SUBROUTINE ComputeIncrement_CoupledAbsorptionDiffusion &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      dt_Input, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, uDR_Option, &
      MatterCoupling_Mask_Option, &
      MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
      IMEX_ImplicitWeight_Option, MatterCollisionRHS_Option, &
      LeptonFlavorCheck_Rate_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: dt_Input
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: AccumulateMatterCoupling_Mask
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nDR)
    LOGICAL, INTENT(inout), OPTIONAL :: &
      MatterCoupling_Mask_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterCollisionExchange_Rate_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterUpdate_Rate_Option &
        (1:nQ,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(in), OPTIONAL :: IMEX_ImplicitWeight_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      LeptonFlavorCheck_Rate_Option &
        (1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterCollisionRHS_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))

    INTEGER :: iCR, iE, iX1, iX2, iX3, iNodeZ
    REAL(DP) :: Radiation, Matter, ExchangeRate
    REAL(DP) :: AbsorptionOpacity, TransportOpacity

    dU_F = Zero
    dU_R = Zero

    AbsorptionOpacity = CoupledDiffusionAbsorptionOpacity()
    TransportOpacity = VariableOpacityDiffusionMeanOpacity

    DO iX3 = iZ_B1_Input(4), iZ_E1_Input(4)
    DO iX2 = iZ_B1_Input(3), iZ_E1_Input(3)
    DO iX1 = iZ_B1_Input(2), iZ_E1_Input(2)
    DO iE  = iZ_B1_Input(1), iZ_E1_Input(1)
    DO iNodeZ = 1, nDOFZ

      Radiation = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,1)
      Matter    = U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,2)

      ExchangeRate = AbsorptionOpacity &
        * ( CoupledDiffusionEquilibriumSlope * Matter - Radiation ) &
        / ( One + dt_Input * AbsorptionOpacity &
                    * ( One + CoupledDiffusionEquilibriumSlope ) )

      dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,1) =   ExchangeRate
      dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,2) = - ExchangeRate

      DO iCR = iCR_G1, iCR_G3
        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,1) &
          = - TransportOpacity * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,1) &
            / ( One + dt_Input * TransportOpacity )
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    IF( PRESENT( uDR_Option ) ) uDR_Option = Zero
    IF( PRESENT( MatterCoupling_Mask_Option ) ) &
      MatterCoupling_Mask_Option = .FALSE.
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) ) &
      MatterCollisionExchange_Rate_Option = Zero
    IF( PRESENT( MatterUpdate_Rate_Option ) ) &
      MatterUpdate_Rate_Option = Zero
    IF( PRESENT( MatterCollisionRHS_Option ) ) &
      MatterCollisionRHS_Option = Zero
    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) ) &
      LeptonFlavorCheck_Rate_Option = Zero

  END SUBROUTINE ComputeIncrement_CoupledAbsorptionDiffusion


  SUBROUTINE ComputeIncrement_AbsorptionDiffusionBalance_Explicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: LLFPenaltyFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: SurfaceTermFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: VolumeTermFactor_Option
    LOGICAL, INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    REAL(DP) :: LLFPenaltyFactor

    LLFPenaltyFactor = AbsorptionDiffusionLLFPenaltyFactor
    IF( PRESENT( LLFPenaltyFactor_Option ) ) &
      LLFPenaltyFactor &
        = LLFPenaltyFactor * LLFPenaltyFactor_Option

    CALL SetAnalyticSourceAwareTransportOpacity &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, .TRUE. )

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             GE, GX, U_F, U_R, dU_R, &
             LLFPenaltyFactor_Option = LLFPenaltyFactor, &
             SurfaceTermFactor_Option = SurfaceTermFactor_Option, &
             VolumeTermFactor_Option = VolumeTermFactor_Option, &
             IncludeObserverCorrections_Option &
               = IncludeObserverCorrections_Option )

    IF( UseSourceAwareLLF ) CALL ClearSourceAwareTransportOpacity

  END SUBROUTINE ComputeIncrement_AbsorptionDiffusionBalance_Explicit


  SUBROUTINE ComputeIncrement_AbsorptionDiffusionBalance &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      dt_Input, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, uDR_Option, &
      MatterCoupling_Mask_Option, &
      MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
      IMEX_ImplicitWeight_Option, MatterCollisionRHS_Option, &
      LeptonFlavorCheck_Rate_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: dt_Input
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: AccumulateMatterCoupling_Mask
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nDR)
    LOGICAL, INTENT(inout), OPTIONAL :: &
      MatterCoupling_Mask_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterCollisionExchange_Rate_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterUpdate_Rate_Option &
        (1:nQ,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(in), OPTIONAL :: IMEX_ImplicitWeight_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      LeptonFlavorCheck_Rate_Option &
        (1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterCollisionRHS_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))

    INTEGER :: iCR, iE, iS, iX1, iX2, iX3
    INTEGER :: iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, Dnu_0
    REAL(DP) :: AbsorptionOpacity, TransportOpacity

    dU_F = Zero
    dU_R = Zero

    DO iS  = 1, nSpecies
    DO iX3 = iZ_B1_Input(4), iZ_E1_Input(4)
    DO iX2 = iZ_B1_Input(3), iZ_E1_Input(3)
    DO iX1 = iZ_B1_Input(2), iZ_E1_Input(2)
    DO iE  = iZ_B1_Input(1), iZ_E1_Input(1)
    DO iNodeZ = 1, nDOFZ

      iNodeX = ( iNodeZ - 1 ) / nDOFE + 1
      iNodeX1 = NodeNumberTableX(1,iNodeX)
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

      AbsorptionOpacity = AbsorptionDiffusionAbsorptionOpacity( X )
      TransportOpacity = VariableOpacityDiffusionOpacity( X )

      IF( iS == 1 )THEN
        Dnu_0 = AbsorptionDiffusionEquilibriumDensity( X )
      ELSE
        Dnu_0 = Half * OpaqueDiffusionBackground
      END IF

      dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) &
        = AbsorptionOpacity &
            * ( Dnu_0 - U_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS) ) &
          / ( One + dt_Input * AbsorptionOpacity )

      DO iCR = iCR_G1, iCR_G3
        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
          = - TransportOpacity &
              * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
            / ( One + dt_Input * TransportOpacity )
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    IF( PRESENT( uDR_Option ) ) uDR_Option = Zero
    IF( PRESENT( MatterCoupling_Mask_Option ) ) &
      MatterCoupling_Mask_Option = .FALSE.
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) ) &
      MatterCollisionExchange_Rate_Option = Zero
    IF( PRESENT( MatterUpdate_Rate_Option ) ) &
      MatterUpdate_Rate_Option = Zero
    IF( PRESENT( MatterCollisionRHS_Option ) ) &
      MatterCollisionRHS_Option = Zero
    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) ) &
      LeptonFlavorCheck_Rate_Option = Zero

  END SUBROUTINE ComputeIncrement_AbsorptionDiffusionBalance


  SUBROUTINE ReportAbsorptionDiffusionOperatorBalance

    REAL(DP), ALLOCATABLE :: dU_F(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_R_EX(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_R_EX_Central(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_R_EX_Volume(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_R_EX_Face_Dn(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_R_EX_Face_Up(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU_R_IM(:,:,:,:,:,:,:)
    INTEGER :: iNodeX, iNodeX1, iNodeZ, iX1
    INTEGER :: OutputUnit
    REAL(DP) :: X, Weight, CellDensityResidual
    REAL(DP) :: DensityResidual, FluxResidual
    REAL(DP) :: DensityScale, FluxScale
    REAL(DP) :: DensityL1, DensityL2, DensityLInf
    REAL(DP) :: FluxL1, FluxL2, FluxLInf
    REAL(DP) :: DensityCheckerboard
    REAL(DP) :: StreamingDensity, StreamingDensityVolume
    REAL(DP) :: StreamingDensityCentral, StreamingDensityPenalty
    REAL(DP) :: StreamingDensityFaceDn, StreamingDensityFaceUp
    REAL(DP) :: StreamingFlux, StreamingFluxVolume
    REAL(DP) :: StreamingFluxCentral, StreamingFluxPenalty
    REAL(DP) :: StreamingFluxFaceDn, StreamingFluxFaceUp
    REAL(DP) :: DensityDecomposition, FluxDecomposition

    ALLOCATE( dU_F, MOLD = uCF )
    ALLOCATE( dU_R_EX, MOLD = uCR )
    ALLOCATE( dU_R_EX_Central, MOLD = uCR )
    ALLOCATE( dU_R_EX_Volume, MOLD = uCR )
    ALLOCATE( dU_R_EX_Face_Dn, MOLD = uCR )
    ALLOCATE( dU_R_EX_Face_Up, MOLD = uCR )
    ALLOCATE( dU_R_IM, MOLD = uCR )

    CALL ComputeIncrement_AbsorptionDiffusionBalance_Explicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, uCF, uCR, dU_R_EX, &
             IncludeObserverCorrections_Option = .FALSE. )

    CALL ComputeIncrement_AbsorptionDiffusionBalance_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, uCF, uCR, dU_R_EX_Central, &
             LLFPenaltyFactor_Option = Zero, &
             IncludeObserverCorrections_Option = .FALSE. )

    CALL ComputeIncrement_AbsorptionDiffusionBalance_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, uCF, uCR, dU_R_EX_Volume, &
             SurfaceTermFactor_Option = Zero, &
             VolumeTermFactor_Option = One, &
             IncludeObserverCorrections_Option = .FALSE. )

    CALL SetStreamingCheckSurfaceTermFactors( One, Zero )
    CALL ComputeIncrement_AbsorptionDiffusionBalance_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, uCF, uCR, dU_R_EX_Face_Dn, &
             SurfaceTermFactor_Option = Zero, &
             VolumeTermFactor_Option = Zero, &
             IncludeObserverCorrections_Option = .FALSE. )
    CALL ClearStreamingCheckSurfaceTermFactors

    CALL SetStreamingCheckSurfaceTermFactors( Zero, One )
    CALL ComputeIncrement_AbsorptionDiffusionBalance_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, uCF, uCR, dU_R_EX_Face_Up, &
             SurfaceTermFactor_Option = Zero, &
             VolumeTermFactor_Option = Zero, &
             IncludeObserverCorrections_Option = .FALSE. )
    CALL ClearStreamingCheckSurfaceTermFactors

    CALL ComputeIncrement_AbsorptionDiffusionBalance &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, Zero, &
             uGE, uGF, uCF, dU_F, uCR, dU_R_IM, &
             AccumulateMatterCoupling_Mask = .FALSE. )

    DensityScale = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude &
      * ( TwoPi / OpaqueDiffusionDomainLength )**2 &
      / ( Three * VariableOpacityDiffusionMeanOpacity )
    FluxScale = OpaqueDiffusionBackground * OpaqueDiffusionAmplitude &
      * ( TwoPi / OpaqueDiffusionDomainLength ) / Three

    DensityL1 = Zero
    DensityL2 = Zero
    DensityLInf = Zero
    FluxL1 = Zero
    FluxL2 = Zero
    FluxLInf = Zero
    DensityCheckerboard = Zero

    OPEN( NEWUNIT = OutputUnit, &
          FILE = 'absorption_diffusion_operator_balance.dat', &
          STATUS = 'REPLACE', ACTION = 'WRITE' )
    WRITE(OutputUnit,'(A)') &
      '# x_km element node streaming_D streaming_central_D ' // &
      'streaming_penalty_D streaming_volume_D ' // &
      'streaming_face_dn_D streaming_face_up_D absorption_D residual_D ' // &
      'decomposition_D streaming_I1 streaming_central_I1 ' // &
      'streaming_penalty_I1 streaming_volume_I1 ' // &
      'streaming_face_dn_I1 streaming_face_up_I1 collision_I1 ' // &
      'residual_I1 decomposition_I1'

    DO iX1 = iX_B0(1), iX_E0(1)

      CellDensityResidual = Zero

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeZ = ( iNodeX - 1 ) * nDOFE + 1
        X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        Weight = WeightsX1(iNodeX1) * MeshX(1) % Width(iX1)

        StreamingDensity &
          = dU_R_EX(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
        StreamingDensityCentral &
          = dU_R_EX_Central &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
        StreamingDensityPenalty &
          = StreamingDensity - StreamingDensityCentral
        StreamingDensityVolume &
          = dU_R_EX_Volume &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
        StreamingDensityFaceDn &
          = dU_R_EX_Face_Dn &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
        StreamingDensityFaceUp &
          = dU_R_EX_Face_Up &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
        StreamingFlux &
          = dU_R_EX(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)
        StreamingFluxCentral &
          = dU_R_EX_Central &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)
        StreamingFluxPenalty = StreamingFlux - StreamingFluxCentral
        StreamingFluxVolume &
          = dU_R_EX_Volume &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)
        StreamingFluxFaceDn &
          = dU_R_EX_Face_Dn &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)
        StreamingFluxFaceUp &
          = dU_R_EX_Face_Up &
              (iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)

        DensityDecomposition = StreamingDensity &
          - ( StreamingDensityVolume + StreamingDensityFaceDn &
                + StreamingDensityFaceUp )
        FluxDecomposition = StreamingFlux &
          - ( StreamingFluxVolume + StreamingFluxFaceDn &
                + StreamingFluxFaceUp )

        DensityResidual &
          = StreamingDensity &
              + dU_R_IM(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1)
        FluxResidual &
          = StreamingFlux &
              + dU_R_IM(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1)

        DensityL1 = DensityL1 + Weight * ABS( DensityResidual )
        DensityL2 = DensityL2 + Weight * DensityResidual**2
        DensityLInf = MAX( DensityLInf, ABS( DensityResidual ) )
        FluxL1 = FluxL1 + Weight * ABS( FluxResidual )
        FluxL2 = FluxL2 + Weight * FluxResidual**2
        FluxLInf = MAX( FluxLInf, ABS( FluxResidual ) )
        CellDensityResidual = CellDensityResidual + Weight * DensityResidual

        WRITE(OutputUnit,'(ES24.16E3,1X,2(I8,1X),18(ES24.16E3,1X))') &
          X / Kilometer, &
          iX1 - iX_B0(1) + 1, iNodeX1, &
          StreamingDensity, StreamingDensityCentral, &
          StreamingDensityPenalty, StreamingDensityVolume, &
          StreamingDensityFaceDn, StreamingDensityFaceUp, &
          dU_R_IM(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_N,1), &
          DensityResidual, DensityDecomposition, &
          StreamingFlux, StreamingFluxCentral, &
          StreamingFluxPenalty, StreamingFluxVolume, &
          StreamingFluxFaceDn, StreamingFluxFaceUp, &
          dU_R_IM(iNodeZ,iE_B0,iX1,iX_B0(2),iX_B0(3),iCR_G1,1), &
          FluxResidual, FluxDecomposition

      END DO

      DensityCheckerboard = DensityCheckerboard &
        + (-One)**(iX1-iX_B0(1)) * CellDensityResidual

    END DO

    CLOSE( OutputUnit )

    DensityL1 = DensityL1 &
      / ( DensityScale * OpaqueDiffusionDomainLength )
    DensityL2 = SQRT( DensityL2 / OpaqueDiffusionDomainLength ) &
      / DensityScale
    DensityLInf = DensityLInf / DensityScale
    FluxL1 = FluxL1 / ( FluxScale * OpaqueDiffusionDomainLength )
    FluxL2 = SQRT( FluxL2 / OpaqueDiffusionDomainLength ) / FluxScale
    FluxLInf = FluxLInf / FluxScale
    DensityCheckerboard = ABS( DensityCheckerboard ) &
      / ( DensityScale * OpaqueDiffusionDomainLength )

    WRITE(*,* )
    WRITE(*,'(A)') 'Absorption-diffusion semi-discrete balance:'
    WRITE(*,'(A,3(1X,ES14.6E3))') '  density rel. L1/L2/Linf     = ', &
      DensityL1, DensityL2, DensityLInf
    WRITE(*,'(A,3(1X,ES14.6E3))') '  flux rel. L1/L2/Linf        = ', &
      FluxL1, FluxL2, FluxLInf
    WRITE(*,'(A,ES14.6E3)') '  density checkerboard ampl.  = ', &
      DensityCheckerboard
    WRITE(*,'(A)') &
      '  profile: absorption_diffusion_operator_balance.dat'
    WRITE(*,* )

    DEALLOCATE( dU_F, dU_R_EX, dU_R_EX_Central, dU_R_EX_Volume, &
                dU_R_EX_Face_Dn, dU_R_EX_Face_Up, dU_R_IM )

  END SUBROUTINE ReportAbsorptionDiffusionOperatorBalance


  SUBROUTINE ComputeIncrement_ManufacturedOpaqueDiffusion_Explicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: LLFPenaltyFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: SurfaceTermFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: VolumeTermFactor_Option
    LOGICAL, INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    CALL SetAnalyticSourceAwareTransportOpacity &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, .FALSE., &
             ConstantOpacity_Option = OpaqueDiffusionOpacity )

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             GE, GX, U_F, U_R, dU_R, &
             LLFPenaltyFactor_Option = LLFPenaltyFactor_Option, &
             SurfaceTermFactor_Option = SurfaceTermFactor_Option, &
             VolumeTermFactor_Option = VolumeTermFactor_Option, &
             IncludeObserverCorrections_Option &
               = IncludeObserverCorrections_Option )

    IF( UseSourceAwareLLF ) CALL ClearSourceAwareTransportOpacity

  END SUBROUTINE ComputeIncrement_ManufacturedOpaqueDiffusion_Explicit


  SUBROUTINE ComputeIncrement_ManufacturedOpaqueDiffusion &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      dt_Input, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, uDR_Option, &
      MatterCoupling_Mask_Option, &
      MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
      IMEX_ImplicitWeight_Option, MatterCollisionRHS_Option, &
      LeptonFlavorCheck_Rate_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: dt_Input
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: AccumulateMatterCoupling_Mask
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nDR)
    LOGICAL, INTENT(inout), OPTIONAL :: &
      MatterCoupling_Mask_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterCollisionExchange_Rate_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterUpdate_Rate_Option &
        (1:nQ,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(in), OPTIONAL :: IMEX_ImplicitWeight_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      LeptonFlavorCheck_Rate_Option &
        (1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterCollisionRHS_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))

    INTEGER :: iCR, iE, iS, iX1, iX2, iX3, iNodeZ
    REAL(DP) :: ImplicitRate

    dU_F = Zero
    dU_R = Zero

    ! Return the collision RHS evaluated at the backward-Euler stage state:
    ! G^{n+1} = G^n / ( 1 + dt*kappa ).
    ImplicitRate = OpaqueDiffusionOpacity &
      / ( One + dt_Input * OpaqueDiffusionOpacity )

    DO iS  = 1, nSpecies
    DO iCR = iCR_G1, iCR_G3
    DO iX3 = iZ_B1_Input(4), iZ_E1_Input(4)
    DO iX2 = iZ_B1_Input(3), iZ_E1_Input(3)
    DO iX1 = iZ_B1_Input(2), iZ_E1_Input(2)
    DO iE  = iZ_B1_Input(1), iZ_E1_Input(1)
    DO iNodeZ = 1, nDOFZ

      dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
        = - ImplicitRate * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    IF( PRESENT( uDR_Option ) ) uDR_Option = Zero
    IF( PRESENT( MatterCoupling_Mask_Option ) ) &
      MatterCoupling_Mask_Option = .FALSE.
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) ) &
      MatterCollisionExchange_Rate_Option = Zero
    IF( PRESENT( MatterUpdate_Rate_Option ) ) &
      MatterUpdate_Rate_Option = Zero
    IF( PRESENT( MatterCollisionRHS_Option ) ) &
      MatterCollisionRHS_Option = Zero
    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) ) &
      LeptonFlavorCheck_Rate_Option = Zero

  END SUBROUTINE ComputeIncrement_ManufacturedOpaqueDiffusion


  SUBROUTINE ComputeIncrement_ManufacturedVariableOpacityDiffusion_Explicit &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      GE, GX, U_F, U_R, dU_R, &
      LLFPenaltyFactor_Option, SurfaceTermFactor_Option, &
      VolumeTermFactor_Option, IncludeObserverCorrections_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in), OPTIONAL :: LLFPenaltyFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: SurfaceTermFactor_Option
    REAL(DP), INTENT(in), OPTIONAL :: VolumeTermFactor_Option
    LOGICAL, INTENT(in), OPTIONAL :: IncludeObserverCorrections_Option

    INTEGER :: iE, iX1, iX2, iX3, iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, DensitySource

    CALL SetAnalyticSourceAwareTransportOpacity &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, .TRUE. )

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
             GE, GX, U_F, U_R, dU_R, &
             LLFPenaltyFactor_Option = LLFPenaltyFactor_Option, &
             SurfaceTermFactor_Option = SurfaceTermFactor_Option, &
             VolumeTermFactor_Option = VolumeTermFactor_Option, &
             IncludeObserverCorrections_Option &
               = IncludeObserverCorrections_Option )

    IF( UseSourceAwareLLF ) CALL ClearSourceAwareTransportOpacity

    DO iX3 = iZ_B0_Input(4), iZ_E0_Input(4)
    DO iX2 = iZ_B0_Input(3), iZ_E0_Input(3)
    DO iX1 = iZ_B0_Input(2), iZ_E0_Input(2)
    DO iE  = iZ_B0_Input(1), iZ_E0_Input(1)
    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1,iNodeX)
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      IF( RunScatteringDiffusionLimit )THEN
        DensitySource = ScatteringDiffusionDensitySource( X )
      ELSE
        DensitySource = VariableOpacityDiffusionDensitySource( X )
      END IF

      DO iNodeZ = (iNodeX-1)*nDOFE+1, iNodeX*nDOFE
        dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,1) &
          = dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR_N,1) &
              + DensitySource
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE &
    ComputeIncrement_ManufacturedVariableOpacityDiffusion_Explicit


  SUBROUTINE ComputeIncrement_ManufacturedVariableOpacityDiffusion &
    ( iZ_B0_Input, iZ_E0_Input, iZ_B1_Input, iZ_E1_Input, &
      dt_Input, GE, GX, U_F, dU_F, U_R, dU_R, &
      AccumulateMatterCoupling_Mask, uDR_Option, &
      MatterCoupling_Mask_Option, &
      MatterCollisionExchange_Rate_Option, MatterUpdate_Rate_Option, &
      IMEX_ImplicitWeight_Option, MatterCollisionRHS_Option, &
      LeptonFlavorCheck_Rate_Option )

    INTEGER, INTENT(in) :: &
      iZ_B0_Input(4), iZ_E0_Input(4), iZ_B1_Input(4), iZ_E1_Input(4)
    REAL(DP), INTENT(in) :: dt_Input
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE,iZ_B1_Input(1):iZ_E1_Input(1),1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(out) :: &
      dU_F(1:nDOFX, &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1_Input(1):iZ_E1_Input(1), &
          iZ_B1_Input(2):iZ_E1_Input(2), &
          iZ_B1_Input(3):iZ_E1_Input(3), &
          iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOFZ, &
           iZ_B1_Input(1):iZ_E1_Input(1), &
           iZ_B1_Input(2):iZ_E1_Input(2), &
           iZ_B1_Input(3):iZ_E1_Input(3), &
           iZ_B1_Input(4):iZ_E1_Input(4),1:nCR,1:nSpecies)
    LOGICAL, INTENT(in) :: AccumulateMatterCoupling_Mask
    REAL(DP), INTENT(inout), OPTIONAL :: &
      uDR_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4),1:nDR)
    LOGICAL, INTENT(inout), OPTIONAL :: &
      MatterCoupling_Mask_Option &
        (iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterCollisionExchange_Rate_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(inout), OPTIONAL :: &
      MatterUpdate_Rate_Option &
        (1:nQ,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(in), OPTIONAL :: IMEX_ImplicitWeight_Option
    REAL(DP), INTENT(out), OPTIONAL :: &
      LeptonFlavorCheck_Rate_Option &
        (1:nLeptonFamilies,1:nLeptonFlavorCheck,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2),iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))
    REAL(DP), INTENT(out), OPTIONAL :: &
      MatterCollisionRHS_Option &
        (1:nQ,1:1+nSpecies+nC,1:nDOFX, &
         iZ_B1_Input(2):iZ_E1_Input(2), &
         iZ_B1_Input(3):iZ_E1_Input(3), &
         iZ_B1_Input(4):iZ_E1_Input(4))

    INTEGER :: iCR, iE, iS, iX1, iX2, iX3
    INTEGER :: iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X, Opacity, ImplicitRate

    dU_F = Zero
    dU_R = Zero

    DO iS  = 1, nSpecies
    DO iCR = iCR_G1, iCR_G3
    DO iX3 = iZ_B1_Input(4), iZ_E1_Input(4)
    DO iX2 = iZ_B1_Input(3), iZ_E1_Input(3)
    DO iX1 = iZ_B1_Input(2), iZ_E1_Input(2)
    DO iE  = iZ_B1_Input(1), iZ_E1_Input(1)
    DO iNodeZ = 1, nDOFZ

      iNodeX = ( iNodeZ - 1 ) / nDOFE + 1
      iNodeX1 = NodeNumberTableX(1,iNodeX)
      X = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      Opacity = VariableOpacityDiffusionOpacity( X )
      ImplicitRate = Opacity / ( One + dt_Input * Opacity )

      dU_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS) &
        = - ImplicitRate * U_R(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    IF( PRESENT( uDR_Option ) ) uDR_Option = Zero
    IF( PRESENT( MatterCoupling_Mask_Option ) ) &
      MatterCoupling_Mask_Option = .FALSE.
    IF( PRESENT( MatterCollisionExchange_Rate_Option ) ) &
      MatterCollisionExchange_Rate_Option = Zero
    IF( PRESENT( MatterUpdate_Rate_Option ) ) &
      MatterUpdate_Rate_Option = Zero
    IF( PRESENT( MatterCollisionRHS_Option ) ) &
      MatterCollisionRHS_Option = Zero
    IF( PRESENT( LeptonFlavorCheck_Rate_Option ) ) &
      LeptonFlavorCheck_Rate_Option = Zero

  END SUBROUTINE ComputeIncrement_ManufacturedVariableOpacityDiffusion


  SUBROUTINE ReportVaryingEquilibriumLimiterInput

    USE ProgramHeaderModule, ONLY: &
      nDOFZ, nDOFE
    USE ReferenceElementModule, ONLY: &
      nDOF_E, nDOF_X1, &
      NodeNumberTable, NodeNumberTable_E, NodeNumberTable_X1, &
      Weights_Q
    USE ReferenceElementModule_Lagrange, ONLY: &
      L_E_Dn, L_E_Up, L_X1_Dn, L_X1_Up
    USE GeometryFieldsModuleE, ONLY: &
      iGE_Ep2
    USE GeometryFieldsModule, ONLY: &
      iGF_SqrtGm
    USE MeshModule, ONLY: &
      MeshE, MeshX, NodeCoordinate
    USE RadiationFieldsModule, ONLY: &
      iCR_G1, iCR_G2, iCR_G3

    INTEGER :: iE, iX1, iX2, iX3, iS, iNodeZ, iNodeE, iNodeX
    INTEGER :: iP
    INTEGER :: MinNLocation(7), MinGammaLocation(7)
    REAL(DP) :: N_P, MinN, MinGamma, MaxAbsG
    REAL(DP) :: SumN, SumW, Weight, N_K, NodalMin, NodalMax
    REAL(DP) :: R_MinN, E_MinN, R_MinGamma, E_MinGamma
    CHARACTER(16) :: PointName_MinN, PointName_MinGamma

    MinN     = HUGE( One )
    MinGamma = HUGE( One )
    MaxAbsG  = Zero

    DO iS  = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0

      DO iNodeZ = 1, nDOFZ

        N_P = uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N,iS)
        CALL RecordLimiterPoint &
               ( N_P, [ iLimiterPoint_Node, iNodeZ, iE, iX1, iX2, iX3, iS ], &
                 MinN, MinNLocation, MinGamma, MinGammaLocation )

        MaxAbsG = MAX( MaxAbsG, &
          ABS( uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS) ), &
          ABS( uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS) ), &
          ABS( uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS) ) )

      END DO

      DO iP = 1, nDOF_E

        N_P = DOT_PRODUCT &
          ( L_E_Dn(iP,1:nDOFZ), &
            uCR(1:nDOFZ,iE,iX1,iX2,iX3,iCR_N,iS) )
        CALL RecordLimiterPoint &
               ( N_P, [ iLimiterPoint_E_Dn, iP, iE, iX1, iX2, iX3, iS ], &
                 MinN, MinNLocation, MinGamma, MinGammaLocation )

        N_P = DOT_PRODUCT &
          ( L_E_Up(iP,1:nDOFZ), &
            uCR(1:nDOFZ,iE,iX1,iX2,iX3,iCR_N,iS) )
        CALL RecordLimiterPoint &
               ( N_P, [ iLimiterPoint_E_Up, iP, iE, iX1, iX2, iX3, iS ], &
                 MinN, MinNLocation, MinGamma, MinGammaLocation )

      END DO

      DO iP = 1, nDOF_X1

        N_P = DOT_PRODUCT &
          ( L_X1_Dn(iP,1:nDOFZ), &
            uCR(1:nDOFZ,iE,iX1,iX2,iX3,iCR_N,iS) )
        CALL RecordLimiterPoint &
               ( N_P, [ iLimiterPoint_X_Dn, iP, iE, iX1, iX2, iX3, iS ], &
                 MinN, MinNLocation, MinGamma, MinGammaLocation )

        N_P = DOT_PRODUCT &
          ( L_X1_Up(iP,1:nDOFZ), &
            uCR(1:nDOFZ,iE,iX1,iX2,iX3,iCR_N,iS) )
        CALL RecordLimiterPoint &
               ( N_P, [ iLimiterPoint_X_Up, iP, iE, iX1, iX2, iX3, iS ], &
                 MinN, MinNLocation, MinGamma, MinGammaLocation )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL LimiterPointCoordinates &
           ( MinNLocation(1), MinNLocation(2), &
             MinNLocation(3), MinNLocation(4), &
             R_MinN, E_MinN, PointName_MinN )
    CALL LimiterPointCoordinates &
           ( MinGammaLocation(1), MinGammaLocation(2), &
             MinGammaLocation(3), MinGammaLocation(4), &
             R_MinGamma, E_MinGamma, PointName_MinGamma )

    SumN = Zero
    SumW = Zero
    DO iNodeZ = 1, nDOFZ

      iNodeE = NodeNumberTable(1,iNodeZ)
      iNodeX = NodeNumberTable(2,iNodeZ)
      Weight = Weights_Q(iNodeZ) &
               * uGE(iNodeE,MinNLocation(3),iGE_Ep2) &
               * uGF(iNodeX,MinNLocation(4),MinNLocation(5), &
                     MinNLocation(6),iGF_SqrtGm)
      SumN = SumN + Weight &
             * uCR(iNodeZ,MinNLocation(3),MinNLocation(4), &
                   MinNLocation(5),MinNLocation(6),iCR_N,MinNLocation(7))
      SumW = SumW + Weight

    END DO
    N_K = SumN / SumW
    NodalMin = MINVAL &
      ( uCR(:,MinNLocation(3),MinNLocation(4),MinNLocation(5), &
              MinNLocation(6),iCR_N,MinNLocation(7)) )
    NodalMax = MAXVAL &
      ( uCR(:,MinNLocation(3),MinNLocation(4),MinNLocation(5), &
              MinNLocation(6),iCR_N,MinNLocation(7)) )

    WRITE(*,*)
    WRITE(*,'(A)') 'Varying-equilibrium limiter-input diagnostic:'
    WRITE(*,'(A,ES14.6E3)') '  density lower bound           = ', SqrtTiny
    WRITE(*,'(A,ES14.6E3)') '  realizability lower bound     = ', SqrtTiny
    WRITE(*,'(A,L1)')       '  density limiter triggered     = ', MinN < SqrtTiny
    WRITE(*,'(A,L1)')       '  nonzero flux in Gamma check    = ', &
      MaxAbsG > Zero
    WRITE(*,'(A,ES14.6E3)') '  minimum point value N         = ', MinN
    WRITE(*,'(A,A)')        '    point set                    = ', &
      TRIM( PointName_MinN )
    WRITE(*,'(A,4(I8,1X))') '    species, E/X element        = ', &
      MinNLocation(7), MinNLocation(3), MinNLocation(4), MinNLocation(2)
    WRITE(*,'(A,2(ES14.6E3,1X))') '    E [MeV], r [km]             = ', &
      E_MinN / MeV, R_MinN / Kilometer
    WRITE(*,'(A,2(ES14.6E3,1X))') '    nodal N min/max             = ', &
      NodalMin, NodalMax
    WRITE(*,'(A,ES14.6E3)') '    weighted cell-average N     = ', N_K
    WRITE(*,'(A,ES14.6E3)') '  pre-density-step min Gamma    = ', MinGamma
    WRITE(*,'(A,A)')        '    point set                    = ', &
      TRIM( PointName_MinGamma )
    WRITE(*,'(A,4(I8,1X))') '    species, E/X element        = ', &
      MinGammaLocation(7), MinGammaLocation(3), &
      MinGammaLocation(4), MinGammaLocation(2)
    WRITE(*,'(A,2(ES14.6E3,1X))') '    E [MeV], r [km]             = ', &
      E_MinGamma / MeV, R_MinGamma / Kilometer
    WRITE(*,'(A,ES14.6E3)') '  maximum nodal |G_i|           = ', MaxAbsG
    WRITE(*,*)

  END SUBROUTINE ReportVaryingEquilibriumLimiterInput


  SUBROUTINE RecordLimiterPoint &
    ( N_Value, Location, MinN, MinNLocation, MinGamma, MinGammaLocation )

    REAL(DP), INTENT(in)    :: N_Value
    INTEGER,  INTENT(in)    :: Location(7)
    REAL(DP), INTENT(inout) :: MinN, MinGamma
    INTEGER,  INTENT(inout) :: MinNLocation(7), MinGammaLocation(7)

    REAL(DP) :: Gamma_Value

    ! The varying-equilibrium initialization has G_i = 0 exactly, so the
    ! limiter's Gamma = N - sqrt(G_i G^i) reduces to N at every point in
    ! the expanded positivity set.  The density step is applied first and
    ! enforces the same SqrtTiny lower bound used by the Gamma step.
    Gamma_Value = N_Value

    IF( N_Value < MinN )THEN
      MinN         = N_Value
      MinNLocation = Location
    END IF

    IF( Gamma_Value < MinGamma )THEN
      MinGamma         = Gamma_Value
      MinGammaLocation = Location
    END IF

  END SUBROUTINE RecordLimiterPoint


  SUBROUTINE LimiterPointCoordinates &
    ( PointType, PointIndex, iE_Point, iX1_Point, R, E, PointName )

    USE ReferenceElementModule, ONLY: &
      NodeNumberTable, NodeNumberTable_E, NodeNumberTable_X1
    USE MeshModule, ONLY: &
      MeshE, MeshX, NodeCoordinate

    INTEGER,       INTENT(in)  :: PointType, PointIndex
    INTEGER,       INTENT(in)  :: iE_Point, iX1_Point
    REAL(DP),      INTENT(out) :: R, E
    CHARACTER(16), INTENT(out) :: PointName

    INTEGER :: iNodeE, iNodeX

    SELECT CASE( PointType )
    CASE( iLimiterPoint_Node )
      iNodeE = NodeNumberTable(1,PointIndex)
      iNodeX = NodeNumberTable(2,PointIndex)
      E = NodeCoordinate( MeshE,    iE_Point,  iNodeE )
      R = NodeCoordinate( MeshX(1), iX1_Point, iNodeX )
      PointName = 'node'
    CASE( iLimiterPoint_E_Dn )
      iNodeX = NodeNumberTable_E(1,PointIndex)
      E = NodeCoordinate &
            ( MeshE % Center(iE_Point), MeshE % Width(iE_Point), - Half )
      R = NodeCoordinate( MeshX(1), iX1_Point, iNodeX )
      PointName = 'lower E face'
    CASE( iLimiterPoint_E_Up )
      iNodeX = NodeNumberTable_E(1,PointIndex)
      E = NodeCoordinate &
            ( MeshE % Center(iE_Point), MeshE % Width(iE_Point), + Half )
      R = NodeCoordinate( MeshX(1), iX1_Point, iNodeX )
      PointName = 'upper E face'
    CASE( iLimiterPoint_X_Dn )
      iNodeE = NodeNumberTable_X1(1,PointIndex)
      E = NodeCoordinate( MeshE, iE_Point, iNodeE )
      R = NodeCoordinate &
            ( MeshX(1) % Center(iX1_Point), &
              MeshX(1) % Width(iX1_Point), - Half )
      PointName = 'lower X1 face'
    CASE( iLimiterPoint_X_Up )
      iNodeE = NodeNumberTable_X1(1,PointIndex)
      E = NodeCoordinate( MeshE, iE_Point, iNodeE )
      R = NodeCoordinate &
            ( MeshX(1) % Center(iX1_Point), &
              MeshX(1) % Width(iX1_Point), + Half )
      PointName = 'upper X1 face'
    END SELECT

  END SUBROUTINE LimiterPointCoordinates


  SUBROUTINE BeginEquilibriumTest &
    ( dt_Test, RequireZeroStreaming )

    REAL(DP), INTENT(in) :: dt_Test
    LOGICAL, INTENT(in) :: RequireZeroStreaming

    REAL(DP), PARAMETER :: TestTolerance = 1.0d-10
    REAL(DP), ALLOCATABLE :: CF_Work(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR_Work(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dCF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dCR(:,:,:,:,:,:,:)
    REAL(DP) :: LimiterFluidChange, LimiterRadiationChange
    REAL(DP) :: StreamingResponse_PreLimiter, StreamingResponse_PostLimiter
    REAL(DP) :: CollisionFluidUpdate_PreLimiter
    REAL(DP) :: CollisionRadiationUpdate_PreLimiter
    REAL(DP) :: CollisionFluidResponse_PostLimiter
    REAL(DP) :: CollisionRadiationResponse_PostLimiter

    LimiterFluidChange = RelativeFluidDifference &
      ( uCF, EquilibriumTest_CF_Initial )
    LimiterRadiationChange = RelativeRadiationDifference &
      ( uCR, EquilibriumTest_CR_Initial )

    ALLOCATE( EquilibriumTest_CF_Reference, MOLD = uCF )
    ALLOCATE( EquilibriumTest_CR_Reference, MOLD = uCR )
    EquilibriumTest_CF_Reference = uCF
    EquilibriumTest_CR_Reference = uCR

    ALLOCATE( CF_Work, MOLD = uCF )
    ALLOCATE( CR_Work, MOLD = uCR )
    ALLOCATE( dCF, MOLD = uCF )
    ALLOCATE( dCR, MOLD = uCR )

    CF_Work = EquilibriumTest_CF_Initial
    CR_Work = EquilibriumTest_CR_Initial
    dCR     = Zero

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, CF_Work, CR_Work, dCR )

    StreamingResponse_PreLimiter = dt_Test &
      * RadiationRateNorm( dCR, EquilibriumTest_CR_Initial )

    CF_Work = EquilibriumTest_CF_Initial
    CR_Work = EquilibriumTest_CR_Initial
    dCF     = Zero
    dCR     = Zero

    CALL ComputeIncrement_TwoMoment_Implicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt_Test, &
             uGE, uGF, CF_Work, dCF, CR_Work, dCR, &
             AccumulateMatterCoupling_Mask = .FALSE. )

    CollisionFluidUpdate_PreLimiter = dt_Test &
      * FluidRateNorm( dCF, EquilibriumTest_CF_Initial )
    CollisionRadiationUpdate_PreLimiter = dt_Test &
      * RadiationRateNorm( dCR, EquilibriumTest_CR_Initial )

    CF_Work = uCF
    CR_Work = uCR
    dCR     = Zero

    CALL ComputeIncrement_TwoMoment_Explicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
             uGE, uGF, CF_Work, CR_Work, dCR )

    StreamingResponse_PostLimiter = dt_Test &
      * RadiationRateNorm( dCR, uCR )

    CF_Work = uCF
    CR_Work = uCR
    dCF     = Zero
    dCR     = Zero

    CALL ComputeIncrement_TwoMoment_Implicit &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt_Test, &
             uGE, uGF, CF_Work, dCF, CR_Work, dCR, &
             AccumulateMatterCoupling_Mask = .FALSE. )

    CollisionFluidResponse_PostLimiter = dt_Test &
      * FluidRateNorm( dCF, uCF )
    CollisionRadiationResponse_PostLimiter = dt_Test &
      * RadiationRateNorm( dCR, uCR )

    WRITE(*,*)
    WRITE(*,'(A)') 'Local-equilibrium operator test:'
    WRITE(*,'(A,ES14.6E3)') '  limiter fluid rel. change    = ', &
      LimiterFluidChange
    WRITE(*,'(A,ES14.6E3)') '  limiter radiation rel. change= ', &
      LimiterRadiationChange
    WRITE(*,'(A,ES14.6E3)') '  pre-limiter streaming resp.  = ', &
      StreamingResponse_PreLimiter
    WRITE(*,'(A,ES14.6E3)') '  post-limiter streaming resp. = ', &
      StreamingResponse_PostLimiter
    IF( .NOT. RequireZeroStreaming ) &
      WRITE(*,'(A)') &
        '    (nonzero streaming is expected for the varying profile)'
    WRITE(*,'(A,ES14.6E3)') '  pre-limiter coll. fluid upd. = ', &
      CollisionFluidUpdate_PreLimiter
    WRITE(*,'(A,ES14.6E3)') '  pre-limiter coll. rad. upd.  = ', &
      CollisionRadiationUpdate_PreLimiter
    WRITE(*,'(A,ES14.6E3)') '  post-limiter coll. fluid resp= ', &
      CollisionFluidResponse_PostLimiter
    WRITE(*,'(A,ES14.6E3)') '  post-limiter coll. rad. resp.= ', &
      CollisionRadiationResponse_PostLimiter

    IF( MAX( CollisionFluidUpdate_PreLimiter, &
             CollisionRadiationUpdate_PreLimiter ) > TestTolerance )THEN
      ERROR STOP 'Local-equilibrium collision test failed.'
    END IF

    IF( RequireZeroStreaming )THEN
      IF( MAX( LimiterFluidChange, LimiterRadiationChange, &
               StreamingResponse_PreLimiter, StreamingResponse_PostLimiter, &
               CollisionFluidResponse_PostLimiter, &
               CollisionRadiationResponse_PostLimiter ) &
            > TestTolerance )THEN
        ERROR STOP 'Homogeneous-equilibrium operator test failed.'
      END IF
    END IF

    IF( RequireZeroStreaming )THEN
      WRITE(*,'(A)') '  result                        = PASS'
    ELSE
      WRITE(*,'(A)') &
        '  pre-limiter collision result = PASS'
    END IF
    WRITE(*,*)

    DEALLOCATE( EquilibriumTest_CF_Initial, EquilibriumTest_CR_Initial )
    DEALLOCATE( CF_Work, CR_Work, dCF, dCR )

  END SUBROUTINE BeginEquilibriumTest


  SUBROUTINE EndEquilibriumTest( RequireStationaryState )

    LOGICAL, INTENT(in) :: RequireStationaryState

    REAL(DP), PARAMETER :: TestTolerance = 1.0d-10
    REAL(DP) :: FluidRelativeChange, RadiationRelativeChange

    FluidRelativeChange = RelativeFluidDifference &
      ( uCF, EquilibriumTest_CF_Reference )
    RadiationRelativeChange = RelativeRadiationDifference &
      ( uCR, EquilibriumTest_CR_Reference )

    WRITE(*,*)
    WRITE(*,'(A)') 'Local-equilibrium evolution test:'
    WRITE(*,'(A,ES14.6E3)') '  fluid-state relative change  = ', &
      FluidRelativeChange
    WRITE(*,'(A,ES14.6E3)') '  radiation relative change    = ', &
      RadiationRelativeChange

    IF( RequireStationaryState )THEN
      IF( MAX( FluidRelativeChange, RadiationRelativeChange ) &
            > TestTolerance )THEN
        ERROR STOP 'Homogeneous-equilibrium evolution test failed.'
      END IF
      WRITE(*,'(A)') '  result                        = PASS'
    ELSE
      WRITE(*,'(A)') &
        '  result                        = STREAMING RESPONSE RECORDED'
    END IF
    WRITE(*,*)

    DEALLOCATE &
      ( EquilibriumTest_CF_Reference, EquilibriumTest_CR_Reference )

  END SUBROUTINE EndEquilibriumTest


  REAL(DP) FUNCTION RelativeFluidDifference( A, B ) RESULT( Difference )

    REAL(DP), INTENT(in) :: &
      A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      B(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF
    REAL(DP) :: Scale

    Difference = Zero

    DO iCF = 1, SIZE(A,5)
      Scale = MAX &
        ( MAXVAL &
            ( ABS( B(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCF) ) ), &
          TINY(One) )
      Difference = MAX &
        ( Difference, &
          MAXVAL &
            ( ABS( A(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCF) &
                 - B(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCF) ) ) / Scale )
    END DO

  END FUNCTION RelativeFluidDifference


  REAL(DP) FUNCTION RelativeRadiationDifference( A, B ) RESULT( Difference )

    REAL(DP), INTENT(in) :: &
      A(1:,iE_B1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:), &
      B(1:,iE_B1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

    INTEGER :: iCR, iS
    REAL(DP) :: Scale

    Difference = Zero

    DO iS = 1, SIZE(A,7)
    DO iCR = 1, SIZE(A,6)
      Scale = MAX &
        ( MAXVAL &
            ( ABS( B(:,iE_B0:iE_E0, &
                       iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCR,iS) ) ), &
          MAXVAL &
            ( ABS( B(:,iE_B0:iE_E0, &
                       iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCR_N,iS) ) ) &
              * EPSILON(One), &
          TINY(One) )
      Difference = MAX &
        ( Difference, &
          MAXVAL &
            ( ABS( A(:,iE_B0:iE_E0, &
                       iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCR,iS) &
                 - B(:,iE_B0:iE_E0, &
                       iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                       iX_B0(3):iX_E0(3),iCR,iS) ) ) / Scale )
    END DO
    END DO

  END FUNCTION RelativeRadiationDifference


  REAL(DP) FUNCTION FluidRateNorm( Rate, State ) RESULT( Norm )

    REAL(DP), INTENT(in) :: &
      Rate (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      State(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iCF
    REAL(DP) :: Scale

    Norm = Zero

    DO iCF = 1, SIZE(Rate,5)
      Scale = MAX &
        ( MAXVAL &
            ( ABS( State(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),iCF) ) ), &
          TINY(One) )
      IF( iCF >= iCF_S1 .AND. iCF <= iCF_S3 ) &
        Scale = MAX &
          ( MAXVAL &
              ( ABS( State(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                             iX_B0(3):iX_E0(3),iCF_E) ) ), &
            TINY(One) )
      Norm = MAX &
        ( Norm, &
          MAXVAL &
            ( ABS( Rate(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),iCF) ) ) / Scale )
    END DO

  END FUNCTION FluidRateNorm


  REAL(DP) FUNCTION RadiationRateNorm( Rate, State ) RESULT( Norm )

    REAL(DP), INTENT(in) :: &
      Rate &
        (1:,iE_B1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:), &
      State &
        (1:,iE_B1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)

    INTEGER :: iCR, iS
    REAL(DP) :: Scale

    Norm = Zero

    DO iS = 1, SIZE(Rate,7)
      Scale = MAX &
        ( MAXVAL &
            ( ABS( State(:,iE_B0:iE_E0, &
                         iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),iCR_N,iS) ) ), &
          TINY(One) )
      DO iCR = 1, SIZE(Rate,6)
        Norm = MAX &
          ( Norm, &
            MAXVAL &
              ( ABS( Rate(:,iE_B0:iE_E0, &
                           iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),iCR,iS) ) ) / Scale )
      END DO
    END DO

  END FUNCTION RadiationRateNorm


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule, ONLY: &
      InitializeTimers
    USE TimersModule_Euler, ONLY: &
      InitializeTimers_Euler
    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
    USE MeshModule, ONLY: &
      MeshX, &
      CreateMesh_Custom
    USE ReferenceElementModuleX, ONLY: &
      InitializeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      InitializeReferenceElementX_Lagrange
    USE GeometryComputationModule, ONLY: &
      ComputeGeometryX
    USE ReferenceElementModuleE, ONLY: &
      InitializeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      InitializeReferenceElementE_Lagrange
    USE GeometryComputationModuleE, ONLY: &
      ComputeGeometryE
    USE ReferenceElementModuleZ, ONLY: &
      InitializeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange
    USE EquationOfStateModule_TABLE, ONLY: &
      InitializeEquationOfState_TABLE, &
      Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y
    USE OpacityModule_TABLE, ONLY: &
      InitializeOpacities_TABLE
    USE TwoMoment_ClosureModule, ONLY: &
      InitializeClosure_TwoMoment
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      InitializeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      InitializePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      InitializeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule, ONLY: &
      InitializeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule, ONLY: &
      InitializePositivityLimiter_TwoMoment
    USE TwoMoment_TallyModule, ONLY: &
      InitializeTally
    USE TwoMoment_TimeSteppingModule, ONLY: &
      Initialize_IMEX_RK

    CALL InitializeTimers

    CALL InitializeTimers_Euler

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = [ 1, 0, 0 ], &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             zoomX_Option &
               = zoomX, &
             nE_Option &
               = nE, &
             swE_Option &
               = 1, &
             bcE_Option &
               = bcE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             zoomE_Option &
               = zoomE, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = TRIM( CoordinateSystem ), &
             ActivateUnits_Option &
               = .TRUE., &
             nSpecies_Option &
               = nSpecies, &
             BasicInitialization_Option &
               = .TRUE. )

    IF( nEquidistantX > 1 )THEN

      CALL CreateMesh_Custom &
             ( MeshX(1), nX(1), nNodes, 1, xL(1), xR(1), &
               nEquidistantX, dEquidistantX, Verbose_Option = .TRUE. )

    END IF

    ! --- Position Space Reference Element and Geometry ---

    CALL InitializeReferenceElementX

    CALL InitializeReferenceElementX_Lagrange

    CALL ComputeGeometryX &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    ! --- Energy Space Reference Element and Geometry ---

    CALL InitializeReferenceElementE

    CALL InitializeReferenceElementE_Lagrange

    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    ! --- Phase Space Reference Element ---

    CALL InitializeReferenceElementZ

    CALL InitializeReferenceElement

    CALL InitializeReferenceElement_Lagrange

    ! --- Initialize Equation of State ---

    CALL InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = EosTableName, &
             UseChemicalPotentialShift_Option &
               = .FALSE., &
             Verbose_Option &
               = .TRUE., &
             MuonEos_MinD_Option = 1.0d12 )

    ! --- Initialize Opacities ---

    CALL InitializeOpacities_TABLE &
           ( OpacityTableName_EmAb_Option &
               = TRIM( OpacityTableName_EmAb ), &
             OpacityTableName_Iso_Option  &
               = TRIM( OpacityTableName_Iso ), &
             OpacityTableName_NNS_Option &
               = TRIM( OpacityTableName_NNS ), &
             OpacityTableName_NES_Option &
               = TRIM( OpacityTableName_NES ), &
             OpacityTableName_Pair_Option &
               = TRIM( OpacityTableName_Pair ), &
             OpacityTableName_Brem_Option &
               = TRIM( OpacityTableName_Brem ), &
             EquationOfStateTableName_Option &
               = TRIM( EosTableName ), &
             Verbose_Option = .TRUE., &
             EmAb_Muon_MinD_Option = 1.0d12 )

    ! --- Initialize Moment Closure ---

    CALL InitializeClosure_TwoMoment

    ! --- Initialize Slope Limiter (Euler) ---

    CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
           ( BetaTVD_Option &
               = 1.75_DP, &
             SlopeTolerance_Option &
               = 1.0d-6, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_Euler, &
             UseTroubledCellIndicator_Option &
               = .FALSE., &
             LimiterThresholdParameter_Option &
               = Zero )

    ! --- Initialize Positivity Limiter (Euler) ---

    CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option &
               = UsePositivityLimiter_Euler, &
             Verbose_Option &
               = .TRUE., &
             Min_1_Option &
               = ( One + 1.0d3 * EPSILON( One ) ) * Min_D, &
             Min_2_Option &
               = ( One + 1.0d3 * EPSILON( One ) ) * Min_T, &
             Min_3_Option &
               = ( One + 1.0d3 * EPSILON( One ) ) * Min_Y, &
             Max_1_Option &
               = ( One - 1.0d3 * EPSILON( One ) ) * Max_D, &
             Max_2_Option &
               = ( One - 1.0d3 * EPSILON( One ) ) * Max_T, &
             Max_3_Option &
               = ( One - 1.0d3 * EPSILON( One ) ) * Max_Y )

    ! --- Initialize Troubled Cell Indicator (Two-Moment) ---

    CALL InitializeTroubledCellIndicator_TwoMoment &
           ( UseTroubledCellIndicator_Option &
               = .FALSE., &
             C_TCI_Option &
               = Zero, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Slope Limiter (Two-Moment) ---

    CALL InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = 1.75_DP, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_TwoMoment, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Positivity Limiter (Two-Moment) ---

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option &
               = SqrtTiny, &
             Min_2_Option &
               = SqrtTiny, &
             UsePositivityLimiter_Option &
               = UsePositivityLimiter_TwoMoment, &
             UseEnergyLimiter_Option &
               = UseEnergyLimiter_TwoMoment, &
             Verbose_Option &
               = .TRUE. )

    ! --- Set Neutrino-Matter Solver Parameters ---

    CALL InitializeNeutrinoMatterSolverParameters &
           ( M_outer_Option &
               = M_outer, &
             M_inner_Option &
               = M_inner, &
             MaxIter_outer_Option &
               = MaxIter_outer, &
             MaxIter_inner_Option &
               = MaxIter_inner, &
             Rtol_inner_Option &
               = Rtol_inner, &
             Rtol_outer_Option &
               = Rtol_outer, &
             Include_NNS_Option &
               = Include_NNS, &
             NNS_ApplyWeakMagnetism_Option &
               = .TRUE., &
             NNS_ApplyManyBodyCorrection_Option &
               = .TRUE., &
             Include_NES_Option &
               = Include_NES, &
             Include_Pair_Option &
               = Include_Pair, &
             Include_NuPair_Option &
               = Include_NuPair, &
             Include_Brem_Option &
               = Include_Brem, &
             Include_LinCorr_Option &
               = Include_LinCorr, &
             wMatrRHS_Option &
               = wMatterRHS, &
             DnuMax_Option &
               = DnuMax, &
             FreezeOpacities_Option &
               = FreezeOpacities, &
             UseConservativeProjection_Option &
               = UseConservativeProjection, &
             UseRadiationOnlySolve_Option &
               = UseRadiationOnlySolve, &
             Verbose_Option &
               = .TRUE. )

    CALL InitializeEquilibriumDeviationLLF &
           ( Use_Option = UseEquilibriumDeviationLLF, &
             DnuMax_Option = DnuMax, &
             MaxFluxFactor_Option = 0.1_DP, &
             Verbose_Option = .TRUE. )

    CALL InitializeSourceAwareLLF &
           ( Use_Option = UseSourceAwareLLF, &
             Verbose_Option = .TRUE. )

    ! --- Initialize Tally ---

    CALL InitializeTally

    ! --- Initialize Time Stepper ---

    CALL Initialize_IMEX_RK &
           ( TRIM( TimeSteppingScheme ), EvolveEuler_Option = EvolveEuler )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE TwoMoment_TimeSteppingModule, ONLY: &
      Finalize_IMEX_RK
    USE TwoMoment_TallyModule, ONLY: &
      FinalizeTally
    USE EquationOfStateModule_TABLE, ONLY: &
      FinalizeEquationOfState_TABLE
    USE OpacityModule_TABLE, ONLY: &
      FinalizeOpacities_TABLE
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      FinalizeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule, ONLY: &
      FinalizeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule, ONLY: &
      FinalizePositivityLimiter_TwoMoment
    USE ReferenceElementModuleX, ONLY: &
      FinalizeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      FinalizeReferenceElementX_Lagrange
    USE ReferenceElementModuleE, ONLY: &
      FinalizeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      FinalizeReferenceElementE_Lagrange
    USE ReferenceElementModuleZ, ONLY: &
      FinalizeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      FinalizeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      FinalizeReferenceElement_Lagrange
    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram
    USE TimersModule_Euler, ONLY: &
      FinalizeTimers_Euler
    USE TwoMoment_TimersModule, ONLY: &
      FinalizeTimers

    CALL Finalize_IMEX_RK

    CALL FinalizeTally

    CALL FinalizeEquationOfState_TABLE

    CALL FinalizeOpacities_TABLE

    CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizeTroubledCellIndicator_TwoMoment

    CALL FinalizeSlopeLimiter_TwoMoment

    CALL FinalizePositivityLimiter_TwoMoment

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

    CALL FinalizeTimers_Euler

    CALL FinalizeTimers

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver_Neutrinos
