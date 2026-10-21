#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_OPACITY
#endif

MODULE NeutrinoOpacitiesComputationModule

  USE KindModule, ONLY: &
    DP, Zero, One, SqrtTiny, Pi, TwoPi, FourPi
  USE PhysicalConstantsModule, ONLY: &
    AvogadroConstantMKS, &
    SpeedOfLightCGS, &
    hbarMeVs
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    PlanckConstant, &
    SpeedOfLight, &
    Centimeter, &
    Kilometer, &
    Gram, &
    Kelvin, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFE, nDOFX, &
    nE, nNodesE
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModule, ONLY: &
    Weights_Q
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InterpMat_E
  USE ReferenceElementModule_Lagrange, ONLY: &
    InterpMat
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeMuonChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeProtonMassFraction_TABLE, &
    ComputeNeutronMassFraction_TABLE, &
    ComputeHeavyMassFraction_TABLE, &
    ComputeHeavyMassNumber_TABLE, &
    ComputeElectronNeutrinoChemicalPotential_TABLE, &
    ComputeMuonNeutrinoChemicalPotential_TABLE, &
    ComputeNeutronEffectiveMass_TABLE, &
    ComputeProtonEffectiveMass_TABLE, &
    ComputeNeutronSelfEnergy_TABLE, &
    ComputeProtonSelfEnergy_TABLE
  USE OpacityModule_TABLE, ONLY: &
#ifdef MICROPHYSICS_WEAKLIB
    OS_EmAb, OS_Iso, OS_NNS, OS_NES, OS_Pair, OS_Brem, &
    EmAb_T, Iso_T, NNS_T, NES_T, Pair_T, Brem_T, &
    NNS_AT, NES_AT, Pair_AT, Brem_AT, &
    NNS_Xi_Nu_N, NNS_Xi_Nu_P, &
    NNS_Xi_NuBar_N, NNS_Xi_NuBar_P, &
    use_EC_table,             &
    OS_EmAb_EC_spec, EmAb_EC_spec_T, &
    OS_EmAb_EC_rate, EmAb_EC_rate_T, & 
    Ds_T, Ts_T, Ys_T, &
    Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T, &
    EC_nE, EC_dE, EC_iE_max, EC_iNodeE_max, &
    EC_kfmin, EC_kfmax, &
    EC_a, EC_b, EC_ak, EC_bk, &
#endif
    LogEs_T, LogDs_T, LogTs_T, Ys_T, &
    LogEtas_T, MuBs_T, &
    C1, C2, C1_NuPair, C2_NuPair, &
    QueryOpacity, &
    QueryOpacity_EmAb, &
    QueryOpacity_EmAb_Nucleon, &
    QueryOpacity_EmAb_Muon, &
    QueryOpacity_ECTable, &
    QueryOpacity_Iso, &
    QueryOpacity_NES, &
    QueryOpacity_NNS, &
    QueryOpacity_Pair, &
    QueryOpacity_Brem, &
    QueryOpacity_NuPair, &
    QueryOpacity_NMS, &
    QueryOpacity_LFE, &
    QueryOpacity_LFC, &
    QueryOpacity_IMD, &
    EmAb_Muon_Method, nE_Muon_2D, &
    MuonCC_QuadratureNodes, MuonCC_QuadratureWeights
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar, iNuM, iNuM_Bar, iNuT, iNuT_Bar, LeptonNumber, &
    iNu, iNu_Bar, nChirals

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules ---

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_3D_Custom,           &
    LogInterpolateSingleVariable_3D_Custom_Point,     &
    LogInterpolateSingleVariable_4D_Custom,           &
    LogInterpolateSingleVariable_4D_Custom_Point,     &
    LogInterpolateSingleVariable_1D3D_Custom,         &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned, &
    SumLogInterpolateSingleVariable_2D2D_Custom_Aligned

  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log, &
    GetIndexAndDelta_Lin, &
    GetIndexAndDelta_Log, &
    LinearInterp1D_1DArray_Point, &
    LinearInterp3D_3DArray_Point

  USE wlOpacityFieldsModule, ONLY: &
    iNeutron_NNS, iProton_NNS

  USE wlNeutralCurrentCorrectionsModule, ONLY: &
    ComputeNCManyBodyCorrection

  USE wlCalculateAbEmOpacityModule, ONLY: &
    ElasticAbsorptionOpacityNue, &
    ElasticAbsorptionOpacityNum, &
    CalculateHorowitzWeakMagRecoil
  USE wlEosConstantsModule, ONLY: &
    mn, mp, me, mmu, kmev
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    CC2D_CallData, CC2D_RowData, &
    CC2D_PrepareCall, CC2D_RowSetup, CC2D_NodeKernel, &
    CC2D_ApplyPrefactor
  ! ----------------------------------------------

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeNeutrinoOpacities_EC
  PUBLIC :: WriteMuonChargedCurrentOpacityDiagnostic
  PUBLIC :: ComputeNeutrinoOpacities_EC_Vector
  PUBLIC :: ComputeNeutrinoOpacities_ES
  PUBLIC :: ComputeNeutrinoOpacities_ES_Vector
  PUBLIC :: ComputeNeutrinoOpacities_NNS
  PUBLIC :: ComputeNeutrinoOpacityRates_NNS
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_NNS
  PUBLIC :: ComputeNeutrinoOpacities_NES
  PUBLIC :: ComputeNeutrinoOpacityRates_NES
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_NES
  PUBLIC :: ComputeNeutrinoOpacities_Pair
  PUBLIC :: ComputeNeutrinoOpacityRates_Pair
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_Pair
  PUBLIC :: ComputeNeutrinoOpacities_NuPair
  PUBLIC :: ComputeNeutrinoOpacityRates_NuPair
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_NuPair
  PUBLIC :: ComputeNeutrinoOpacities_Brem
  PUBLIC :: ComputeNeutrinoOpacityRates_Brem
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_Brem
  PUBLIC :: ComputeEquilibriumDistributions_Point
  PUBLIC :: ComputeEquilibriumDistributions
  PUBLIC :: ComputeEquilibriumDistributions_DG
  PUBLIC :: LimitEquilibriumDistributions_DG
  PUBLIC :: ComputeEquilibriumDistributionAndDerivatives
  PUBLIC :: FermiDirac
  PUBLIC :: dFermiDiracdT
  PUBLIC :: dFermiDiracdY
  PUBLIC :: SetOpacityMask
  PUBLIC :: SetChargedCurrentMasks
  PUBLIC :: ComputeDetailedBalanceFactors
  PUBLIC :: GetCollisionKernelScale

  PUBLIC :: iM0, iM1
  PUBLIC :: nOp
  PUBLIC :: iOp_EmAb
  PUBLIC :: iOp_EmAb_Electron
  PUBLIC :: iOp_EmAb_Muon
  PUBLIC :: iOp_ECTable
  PUBLIC :: iOp_Iso
  PUBLIC :: iOp_NNS
  PUBLIC :: iOp_NES
  PUBLIC :: iOp_Pair
  PUBLIC :: iOp_NuPair
  PUBLIC :: iOp_Brem
  PUBLIC :: iOp_NMS
  PUBLIC :: iOp_LFE
  PUBLIC :: iOp_LFC
  PUBLIC :: iOp_IMD

  REAL(DP), PARAMETER :: UnitD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitT    = Kelvin
  REAL(DP), PARAMETER :: UnitY    = One
  REAL(DP), PARAMETER :: UnitYe   = One
  REAL(DP), PARAMETER :: UnitYm   = One
  REAL(DP), PARAMETER :: UnitE    = MeV
  REAL(DP), PARAMETER :: UnitMn   = MeV
  REAL(DP), PARAMETER :: UnitMp   = MeV
  REAL(DP), PARAMETER :: UnitMe   = MeV
  REAL(DP), PARAMETER :: UnitMmu  = MeV
  REAL(DP), PARAMETER :: UnitEta  = One
  REAL(DP), PARAMETER :: UnitEC   = One / Centimeter
  REAL(DP), PARAMETER :: UnitEC2  = One / Kilometer
  REAL(DP), PARAMETER :: UnitES   = One / ( Centimeter * MeV**2 )
  REAL(DP), PARAMETER :: UnitNNS  = One / ( Centimeter * MeV**3 ) !Check this! TODO
  REAL(DP), PARAMETER :: UnitNES  = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitPair = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitBrem = One / ( Centimeter * MeV**3 )

  REAL(DP), PARAMETER :: C_A      = -1.26d0/2.0d0 ! C_A from HR98
  REAL(DP), PARAMETER :: G_F      = 1.166d-11 ! G_F from HR98 in [MeV**-2]
  REAL(DP), PARAMETER :: nb_14    = AvogadroConstantMKS * 1d14
  REAL(DP), PARAMETER :: Brem_const = C_A**2 * G_F**2 * nb_14 &
                                    * hbarMeVs**2 * SpeedOfLightCGS**2 &
                                    / TwoPi**3 * UnitBrem

  REAL(dp), PARAMETER :: Alpha_Brem(3) = [ 1.0d0, 1.0d0, 28.d0/3.d0 ]
  ! Weak magnetism and pseudoscalar terms; no form-factor correction.
  INTEGER,  PARAMETER :: MuonCCCorrections(3) = [ 1, 1, 0 ]
  REAL(dp), PARAMETER :: coeff = 2.0d0 &
                               * (Pi*SpeedOfLightCGS)**2 &
                               * hbarMeVs**3


  INTEGER,  PARAMETER :: iM0         = 1
  INTEGER,  PARAMETER :: iM1         = 2

  INTEGER,  PARAMETER :: nOp         = 14
  INTEGER,  PARAMETER :: iOp_EmAb    = 1
  INTEGER,  PARAMETER :: iOp_ECTable = 2
  INTEGER,  PARAMETER :: iOp_Iso     = 3
  INTEGER,  PARAMETER :: iOp_NNS     = 4
  INTEGER,  PARAMETER :: iOp_NES     = 5
  INTEGER,  PARAMETER :: iOp_Pair    = 6
  INTEGER,  PARAMETER :: iOp_NuPair  = 7
  INTEGER,  PARAMETER :: iOp_Brem    = 8
  INTEGER,  PARAMETER :: iOp_EmAb_Electron = 9
  INTEGER,  PARAMETER :: iOp_EmAb_Muon     = 10
  INTEGER,  PARAMETER :: iOp_NMS           = 11
  INTEGER,  PARAMETER :: iOp_LFE           = 12
  INTEGER,  PARAMETER :: iOp_LFC           = 13
  INTEGER,  PARAMETER :: iOp_IMD           = 14

  INTERFACE ComputeEquilibriumDistributions_DG
    MODULE PROCEDURE ComputeEquilibriumDistributions_DG_E
    MODULE PROCEDURE ComputeEquilibriumDistributions_DG_Z
  END INTERFACE ComputeEquilibriumDistributions_DG

  INTERFACE FermiDirac
    MODULE PROCEDURE FermiDirac_Scalar
    MODULE PROCEDURE FermiDirac_Vector
  END INTERFACE

  INTERFACE dFermiDiracdT
    MODULE PROCEDURE dFermiDiracdT_Scalar
    MODULE PROCEDURE dFermiDiracdT_Vector
  END INTERFACE

  INTERFACE dFermiDiracdY
    MODULE PROCEDURE dFermiDiracdY_Scalar
    MODULE PROCEDURE dFermiDiracdY_Vector
  END INTERFACE

CONTAINS


  PURE FUNCTION GetCollisionKernelScale( Process ) RESULT(Scale)
    INTEGER, INTENT(in) :: Process
    REAL(DP) :: Scale
    Scale = One
    SELECT CASE(Process)
      CASE(1)
        Scale = UnitNNS
      CASE(2)
        Scale = UnitNES
      CASE(3,4)
        Scale = UnitPair
      CASE(5)
        Scale = 3.0_DP*Brem_const
    END SELECT
  END FUNCTION GetCollisionKernelScale


  SUBROUTINE ComputeDetailedBalanceFactors &
    ( iE_B, iE_E, iX_B, iX_E, E, T, ComputeScattering, ComputePairs, &
      ScatteringFactor, PairFactor )

    INTEGER,  INTENT(in)  :: iE_B, iE_E, iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E), T(iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: ComputeScattering, ComputePairs
    REAL(DP), INTENT(out) :: &
      ScatteringFactor(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: &
      PairFactor(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    INTEGER  :: iE1, iE2, iX
    ! These factors depend only on the current opacity temperature and the
    ! fixed energy grid.  Build them once per opacity state and reuse them in
    ! every species, angular moment, and collision process.

    IF( ComputePairs )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E

      PairFactor(iE1,iE2,iX) &
        = DetailedBalance_Pair( E(iE1), E(iE2), T(iX) )

    END DO
    END DO
    END DO

    END IF

    IF( ComputeScattering )THEN

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E

      ScatteringFactor(iE1,iE2,iX) &
        = DetailedBalance_InelasticScattering &
            ( MIN( E(iE1), E(iE2) ), MAX( E(iE1), E(iE2) ), T(iX) )

    END DO
    END DO
    END DO

    END IF

  END SUBROUTINE ComputeDetailedBalanceFactors


  FUNCTION DetailedBalance_InelasticScattering( E_Low, E_High, T ) &
    RESULT( DetailedBalance )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! The kernel tables store the de-excitation direction.  Microscopic
    ! reversibility fixes the inverse direction through the matter bath and
    ! therefore depends on the energy transfer and temperature, not on a
    ! numerically limited representation of the equilibrium distribution.

    REAL(DP), INTENT(in) :: E_Low, E_High, T
    REAL(DP) :: DetailedBalance
    REAL(DP) :: Exponent

    Exponent = - ( E_High - E_Low ) / ( BoltzmannConstant * T )

    IF ( Exponent <= LOG( TINY( One ) ) ) THEN
      DetailedBalance = Zero
    ELSE
      DetailedBalance = EXP( MIN( Exponent, Zero ) )
    END IF

  END FUNCTION DetailedBalance_InelasticScattering


  FUNCTION DetailedBalance_Pair( E_1, E_2, T ) RESULT( DetailedBalance )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! For a neutrino-antineutrino pair the equilibrium chemical potentials
    ! cancel.  Production/annihilation detailed balance is consequently the
    ! Boltzmann factor for the total pair energy.

    REAL(DP), INTENT(in) :: E_1, E_2, T
    REAL(DP) :: DetailedBalance
    REAL(DP) :: Exponent

    Exponent = - ( E_1 + E_2 ) / ( BoltzmannConstant * T )

    IF ( Exponent <= LOG( TINY( One ) ) ) THEN
      DetailedBalance = Zero
    ELSE
      DetailedBalance = EXP( MIN( Exponent, Zero ) )
    END IF

  END FUNCTION DetailedBalance_Pair


  FUNCTION CachedDetailedBalanceScattering &
    ( iE1, iE2, iX, E, T, Cache_Option ) RESULT( DetailedBalance )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iE1, iE2, iX
    REAL(DP), INTENT(in) :: E(:), T(:)
    REAL(DP), INTENT(in), OPTIONAL :: Cache_Option(:,:,:)
    REAL(DP) :: DetailedBalance

    IF( PRESENT( Cache_Option ) )THEN
      DetailedBalance = Cache_Option(iE1,iE2,iX)
    ELSE
      DetailedBalance = DetailedBalance_InelasticScattering &
        ( MIN( E(iE1), E(iE2) ), MAX( E(iE1), E(iE2) ), T(iX) )
    END IF

  END FUNCTION CachedDetailedBalanceScattering


  FUNCTION CachedDetailedBalancePair &
    ( iE1, iE2, iX, E, T, Cache_Option ) RESULT( DetailedBalance )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iE1, iE2, iX
    REAL(DP), INTENT(in) :: E(:), T(:)
    REAL(DP), INTENT(in), OPTIONAL :: Cache_Option(:,:,:)
    REAL(DP) :: DetailedBalance

    IF( PRESENT( Cache_Option ) )THEN
      DetailedBalance = Cache_Option(iE1,iE2,iX)
    ELSE
      DetailedBalance = DetailedBalance_Pair( E(iE1), E(iE2), T(iX) )
    END IF

  END FUNCTION CachedDetailedBalancePair


  FUNCTION FermiDiracComplement_Scalar( E, Mu, kT ) RESULT( Complement )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! Evaluate 1-f without subtracting a value that may round to one.

    REAL(DP), INTENT(in) :: E, Mu, kT
    REAL(DP) :: Complement

    Complement = FermiDirac_Scalar( - E, - Mu, kT )

  END FUNCTION FermiDiracComplement_Scalar


  SUBROUTINE ComputeEquilibriumDistributions_Point &
    ( E, D, T, Ye, Ym, f0, iSpecies, &
      ElectronChargedCurrentMask_Option, MuonChargedCurrentMask_Option )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distributions (Single E,D,T,Y) ---

    REAL(DP), INTENT(in)  :: E
    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: f0
    INTEGER,  INTENT(in)  :: iSpecies
    LOGICAL,  INTENT(in), OPTIONAL :: ElectronChargedCurrentMask_Option
    LOGICAL,  INTENT(in), OPTIONAL :: MuonChargedCurrentMask_Option

    REAL(DP) :: M, Mnu, kT
    LOGICAL  :: ElectronChargedCurrentActive, MuonChargedCurrentActive

    ! --- Compute Chemical Potentials ---

    kT = BoltzmannConstant * T

    ElectronChargedCurrentActive = .TRUE.
    MuonChargedCurrentActive     = .TRUE.
    IF( PRESENT( ElectronChargedCurrentMask_Option ) ) &
      ElectronChargedCurrentActive = ElectronChargedCurrentMask_Option
    IF( PRESENT( MuonChargedCurrentMask_Option ) ) &
      MuonChargedCurrentActive = MuonChargedCurrentMask_Option

    IF (iSpecies == iNuE .OR. iSpecies == iNuE_Bar) THEN
      CALL ComputeElectronNeutrinoChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mnu, &
             ChargedCurrentActive_Option = ElectronChargedCurrentActive )
      M = Mnu * LeptonNumber(iSpecies)
    ELSE IF (iSpecies == iNuM .OR. iSpecies == iNuM_Bar) THEN
      CALL ComputeMuonNeutrinoChemicalPotential_TABLE &
             ( D, T, Ye, Ym, Mnu, &
               ChargedCurrentActive_Option = MuonChargedCurrentActive )
      M = Mnu * LeptonNumber(iSpecies)
    ELSE
      M = Zero
    END IF

    f0 = FermiDirac( E, M, kT )

  END SUBROUTINE ComputeEquilibriumDistributions_Point


  SUBROUTINE ComputeEquilibriumDistributions &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, f0, &
      OpacityMask_Option, ElectronChargedCurrentMask_Option, &
      MuonChargedCurrentMask_Option )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(out) :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL, INTENT(in), OPTIONAL :: OpacityMask_Option(:,iX_B:)
    LOGICAL, INTENT(in), OPTIONAL :: &
      ElectronChargedCurrentMask_Option(iX_B:iX_E)
    LOGICAL, INTENT(in), OPTIONAL :: &
      MuonChargedCurrentMask_Option(iX_B:iX_E)

    REAL(DP) :: Mnue(iX_B:iX_E)
    REAL(DP) :: Mnum(iX_B:iX_E)
    REAL(DP) :: M, kT
    LOGICAL  :: ElectronChargedCurrentMask(iX_B:iX_E)
    LOGICAL  :: MuonChargedCurrentMask(iX_B:iX_E)
    INTEGER  :: iE, iS, iX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: f0, Mnue, Mnum ) &
    !$OMP MAP( to: E, D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( f0, Mnue, Mnum ) &
    !$ACC COPYIN( E, D, T, Ye, Ym )
#endif

    ! --- Compute Chemical Potentials ---

    ElectronChargedCurrentMask = .TRUE.
    MuonChargedCurrentMask     = .TRUE.

    IF( PRESENT( OpacityMask_Option ) )THEN
      ElectronChargedCurrentMask &
        = OpacityMask_Option(iOp_EmAb_Electron,iX_B:iX_E) .OR. &
          OpacityMask_Option(iOp_ECTable,       iX_B:iX_E)
      MuonChargedCurrentMask &
        = OpacityMask_Option(iOp_EmAb_Muon,iX_B:iX_E)
    ELSE
      IF( PRESENT( ElectronChargedCurrentMask_Option ) ) &
        ElectronChargedCurrentMask = ElectronChargedCurrentMask_Option
      IF( PRESENT( MuonChargedCurrentMask_Option ) ) &
        MuonChargedCurrentMask = MuonChargedCurrentMask_Option
    END IF

    CALL ComputeElectronNeutrinoChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mnue, &
             ChargedCurrentMask_Option = ElectronChargedCurrentMask )
    CALL ComputeMuonNeutrinoChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mnum, &
             ChargedCurrentMask_Option = MuonChargedCurrentMask )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( M, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( M, kT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( M, kT )
#endif

    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E

      kT = BoltzmannConstant * T(iX)
      IF (iS <= iNuE_Bar) THEN
        M = Mnue(iX) * LeptonNumber(iS)
      ELSE IF(iS == iNuM .OR. iS == iNuM_Bar) THEN
        M = Mnum(iX) * LeptonNumber(iS)
      ELSE
        M = Zero
      END IF

      f0(iE,iS,iX) &
        = FermiDirac_Scalar( E(iE), M, kT )

    END DO
    END DO
    END DO
    !DO iS = iS_B, iS_E
      !DO iE = iE_B, iE_E
      !  WRITE(*,*)  'EQUI', iS, iE, E(iE) / UnitE, f0(iE,iS,1)
      !END DO
    !END DO
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Mnue, Mnum, E, D, T, Ye, Ym ) &
    !$OMP MAP( from: f0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Mnue, Mnum, E, D, T, Ye, Ym ) &
    !$ACC COPYOUT( f0 )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions


  SUBROUTINE LimitEquilibriumDistributions_DG &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, f0, f0_Max_Option )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)   , TARGET :: E(iE_B:iE_E)
    REAL(DP), INTENT(inout), TARGET :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in),  OPTIONAL :: f0_Max_Option

    REAL(DP), POINTER :: E_Q(:,:), f0_Q(:,:,:,:)

    REAL(DP), ALLOCATABLE :: f0_K(:,:,:)

    REAL(DP) :: N_K, V_K, f0_Min, f0_Max, f0_P, Min_K, Max_K, Theta
    INTEGER  :: iE, iS, iX, iNodeE, iP_E, nE, nS, nX

    IF ( PRESENT( f0_Max_Option ) ) THEN
      f0_Max = MIN( f0_Max_Option, One - EPSILON( One ) )
    ELSE
      f0_Max = One - EPSILON( One )
    ENDIF

    nE = ( iE_E - iE_B + 1 ) / nDOFE
    nS = ( iS_E - iS_B + 1 )
    nX = ( iX_E - iX_B + 1 )

    ALLOCATE( f0_K(nE,nS,nX) )

    ! --- Permute Data ---

    E_Q (1:nDOFE,1:nE)           => E (:)
    f0_Q(1:nDOFE,1:nE,1:nS,1:nX) => f0(:,:,:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: f0_K ) &
    !$OMP MAP( to: E, f0 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( f0_K ) &
    !$ACC COPYIN( E, f0 )
#endif

    ! --- Cell Average of Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( V_K, N_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( V_K, N_K ) &
    !$ACC PRESENT( WeightsE )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( V_K, N_K )
#endif
    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      V_K = Zero
      N_K = Zero

      DO iNodeE = 1, nDOFE
        V_K = V_K + WeightsE(iNodeE) * E_Q(iNodeE,iE)**2
        N_K = N_K + WeightsE(iNodeE) * E_Q(iNodeE,iE)**2 * f0_Q(iNodeE,iE,iS,iX)
      END DO

      f0_K(iE,iS,iX) = N_K / V_K

      IF( f0_K(iE,iS,iX) > f0_Max )THEN
        f0_K(iE,iS,iX) = f0_Max
        DO iNodeE = 1, nDOFE
          f0_Q(iNodeE,iE,iS,iX) = f0_Max
        END DO
      END IF

    END DO
    END DO
    END DO

    ! --- Limit Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta, f0_P )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( f0_Min, Min_K, Max_K, Theta, f0_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta, f0_P )
#endif
    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      IF( iE < nE ) THEN
        f0_Min = f0_K(iE+1,iS,iX)
      ELSE IF( nE > 1 )THEN
        ! --- Estimate Cell Average in Outer Ghost Element ---
        ! Ratio first avoids underflow of the square; a vanished tail is zero.
        f0_Min = Zero
        IF( f0_K(nE-1,iS,iX) > Zero ) &
          f0_Min = f0_K(nE,iS,iX) &
                     * ( f0_K(nE,iS,iX) / f0_K(nE-1,iS,iX) )
      ELSE
        f0_Min = Zero
      END IF

      Max_K = f0_Max
      Min_K = f0_Min
      DO iP_E = 1, nDOFE+2
        f0_P = Zero
        DO iNodeE = 1, nDOFE
          f0_P = f0_P + InterpMat_E(iP_E,iNodeE) * f0_Q(iNodeE,iE,iS,iX)
        END DO
        Max_K = MAX( Max_K, f0_P )
        Min_K = MIN( Min_K, f0_P )
      END DO

      IF( Min_K < f0_Min .OR. Max_K > f0_Max )THEN

        Theta &
          = MIN( One, &
                 ABS((f0_Min-f0_K(iE,iS,iX))/(Min_K-f0_K(iE,iS,iX)+SqrtTiny)), &
                 ABS((f0_Max-f0_K(iE,iS,iX))/(Max_K-f0_K(iE,iS,iX)+SqrtTiny)) )

        DO iNodeE = 1, nDOFE
          f0_Q(iNodeE,iE,iS,iX) &
            = ( One - Theta ) * f0_K(iE,iS,iX) &
              +       Theta   * f0_Q(iNodeE,iE,iS,iX)
        END DO

      END IF

    END DO
    END DO
    END DO
    !DO iS = iS_B, iS_E
    !  DO iE = iE_B, iE_E
    !    WRITE(*,*)  'LIMIT', iS, iE, E(iE) / UnitE, f0(iE,iS,1)
    !  END DO
    !END DO
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: f0_K, E ) &
    !$OMP MAP( from: f0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( f0_K, E ) &
    !$ACC COPYOUT( f0 )
#endif

  END SUBROUTINE LimitEquilibriumDistributions_DG


  SUBROUTINE ComputeEquilibriumDistributions_DG_E &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, f0, f0_Max_Option )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(out) :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: f0_Max_Option

    REAL(DP) :: f0_Max
    LOGICAL, ALLOCATABLE :: ElectronChargedCurrentMask(:)
    LOGICAL, ALLOCATABLE :: MuonChargedCurrentMask(:)

    ALLOCATE( ElectronChargedCurrentMask(iX_B:iX_E) )
    ALLOCATE( MuonChargedCurrentMask(iX_B:iX_E) )

    CALL SetChargedCurrentMasks &
           ( iX_E-iX_B+1, nDOFX, D, T, Ye, Ym, &
             ElectronChargedCurrentMask, MuonChargedCurrentMask )

    IF ( PRESENT( f0_Max_Option ) ) THEN
      f0_Max = MIN( f0_Max_Option, One - EPSILON( One ) )
    ELSE
      f0_Max = One - EPSILON( One )
    ENDIF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: f0 ) &
    !$OMP MAP( to: E, D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( f0 ) &
    !$ACC COPYIN( E, D, T, Ye, Ym )
#endif

    CALL ComputeEquilibriumDistributions &
             ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, f0, &
             ElectronChargedCurrentMask_Option = ElectronChargedCurrentMask, &
             MuonChargedCurrentMask_Option = MuonChargedCurrentMask )

    CALL LimitEquilibriumDistributions_DG &
           ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, f0, f0_Max )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E, D, T, Ye, Ym ) &
    !$OMP MAP( from: f0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( E, D, T, Ye, Ym ) &
    !$ACC COPYOUT( f0 )
#endif

    DEALLOCATE( ElectronChargedCurrentMask, MuonChargedCurrentMask )

  END SUBROUTINE ComputeEquilibriumDistributions_DG_E


  SUBROUTINE ComputeEquilibriumDistributions_DG_Z &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, SqrtGm, f0, f0_Max_Option )

    INTEGER,  INTENT(in)          :: iE_B, iE_E
    INTEGER,  INTENT(in)          :: iS_B, iS_E
    INTEGER,  INTENT(in)          :: iX_B, iX_E
    REAL(DP), INTENT(in) , TARGET :: E (iE_B:iE_E)
    REAL(DP), INTENT(in) , TARGET :: D                     (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: T                     (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: Ye                    (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: Ym                    (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: SqrtGm                (iX_B:iX_E)
    REAL(DP), INTENT(out), TARGET :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL:: f0_Max_Option

    REAL(DP), PARAMETER :: f0_Min = SqrtTiny
    REAL(DP), PARAMETER :: EXPP   = EXP( + One )
    REAL(DP), PARAMETER :: EXPM   = EXP( - One )

    INTEGER  :: nE, nX, nS
    INTEGER  :: iE, iX, iS, iS0, iP
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    REAL(DP) :: kT, M, N_K, V_K, f0_P
    REAL(DP) :: Min_K, Max_K, Theta, f0_Max
    REAL(DP), POINTER :: E_Q(:,:), D_Q(:,:), T_Q(:,:)
    REAL(DP), POINTER :: Mnue_Q(:,:), Mnum_Q(:,:), SqrtGm_Q(:,:)
    REAL(DP), ALLOCATABLE, TARGET :: Mnue(:)
    REAL(DP), ALLOCATABLE, TARGET :: Mnum(:)
    REAL(DP), ALLOCATABLE :: Tau_Q(:,:,:), f0_K(:,:,:), f0_Q(:,:,:,:)
    LOGICAL, ALLOCATABLE :: ElectronChargedCurrentMask(:)
    LOGICAL, ALLOCATABLE :: MuonChargedCurrentMask(:)

    IF ( PRESENT( f0_Max_Option ) ) THEN
      f0_Max = MIN( f0_Max_Option, One - EPSILON( One ) )
    ELSE
      f0_Max = One - EPSILON( One )
    ENDIF

    nE = (iE_E-iE_B+1)/nDOFE ! --- Number of energy  elements
    nX = (iX_E-iX_B+1)/nDOFX ! --- Number of spatial elements
    nS = (iS_E-iS_B+1)       ! --- Number of species

    ALLOCATE( Mnue(iX_B:iX_E) )
    ALLOCATE( Mnum(iX_B:iX_E) )
    ALLOCATE( Tau_Q(nDOFE*nDOFX,nE,nX) )
    ALLOCATE( f0_K(nE,nS,nX) )
    ALLOCATE( f0_Q(nDOFE*nDOFX,nE,nS,nX) )
    ALLOCATE( ElectronChargedCurrentMask(iX_B:iX_E) )
    ALLOCATE( MuonChargedCurrentMask(iX_B:iX_E) )

    ! --- Compute Chemical Potentials ---

    CALL SetChargedCurrentMasks &
           ( iX_E-iX_B+1, nDOFX, D, T, Ye, Ym, &
             ElectronChargedCurrentMask, MuonChargedCurrentMask )
    CALL ComputeElectronNeutrinoChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mnue, &
             ChargedCurrentMask_Option = ElectronChargedCurrentMask )
    CALL ComputeMuonNeutrinoChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mnum, &
             ChargedCurrentMask_Option = MuonChargedCurrentMask )

    E_Q   (1:nDOFE,1:nE) => E   (:)
    D_Q   (1:nDOFX,1:nX) => D   (:)
    T_Q   (1:nDOFX,1:nX) => T   (:)
    Mnue_Q(1:nDOFX,1:nX) => Mnue(:)
    Mnum_Q(1:nDOFX,1:nX) => Mnum(:)

    SqrtGm_Q(1:nDOFX,1:nX) => SqrtGm(:)

    ! --- Compute Fermi-Dirac Distribution in Elements ---

    DO iX = 1, nX
    DO iS = iS_B, iS_E
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        kT = BoltzmannConstant * T_Q(iNodeX,iX)

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        kT = BoltzmannConstant * T_Q(iNodeX,iX)
        IF(iS <= iNuE_Bar) THEN
          M = Mnue_Q(iNodeX,iX) * LeptonNumber(iS)
        ELSE IF (iS == iNuM .OR. iS == iNuM_Bar) THEN
          M = Mnum_Q(iNodeX,iX) * LeptonNumber(iS)
        ELSE
          M = Zero
        END IF

        f0_Q(iNodeZ,iE,iS,iX) &
          = FermiDirac_Scalar( E_Q(iNodeE,iE), M, kT )

      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Energy-Position Space Volume Jacobian ---

    DO iX = 1, nX
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        Tau_Q(iNodeZ,iE,iX) = E_Q(iNodeE,iE)**2 * SqrtGm_Q(iNodeX,iX)

      END DO
      END DO

    END DO
    END DO

    ! --- Cell Average of Equilibrium Distribution ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      V_K = Zero
      N_K = Zero

      DO iNodeZ = 1, nDOFE * nDOFX

        V_K = V_K + Weights_Q(iNodeZ) * Tau_Q(iNodeZ,iE,iX)
        N_K = N_K + Weights_Q(iNodeZ) * Tau_Q(iNodeZ,iE,iX) &
                      * f0_Q(iNodeZ,iE,iS,iX)

      END DO

      f0_K(iE,iS,iX) = MAX( MIN( N_K / V_K, f0_Max ), f0_Min )

    END DO
    END DO
    END DO

    ! --- Apply Bound-Enforcing Limiter to Equilibrium Distribution ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      Max_K = MIN( f0_Max, f0_K(iE,iS,iX) * EXPP )
      Min_K = MAX( f0_Min, f0_K(iE,iS,iX) * EXPM )
      DO iP = 1, SIZE( InterpMat, 1 )
        f0_P = Zero
        DO iNodeZ = 1, nDOFE * nDOFX
          f0_P = f0_P + InterpMat(iP,iNodeZ) * f0_Q(iNodeZ,iE,iS,iX)
        END DO
        Max_K = MAX( Max_K, f0_P )
        Min_K = MIN( Min_K, f0_P )
      END DO

      IF( Min_K < f0_Min .OR. Max_K > f0_Max )THEN

        Theta &
          = MIN( One, &
                 ABS((f0_Min-f0_K(iE,iS,iX))/(Min_K-f0_K(iE,iS,iX)+SqrtTiny)), &
                 ABS((f0_Max-f0_K(iE,iS,iX))/(Max_K-f0_K(iE,iS,iX)+SqrtTiny)) )

        Theta = ( One - 1.0d-2 ) * Theta

        DO iNodeZ = 1, nDOFE * nDOFX

          f0_Q(iNodeZ,iE,iS,iX) &
            = ( One - Theta ) * f0_K(       iE,iS,iX) &
              +       Theta   * f0_Q(iNodeZ,iE,iS,iX)

        END DO

      END IF

    END DO
    END DO
    END DO

    ! --- Permute for Output ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        f0((iE-1)*nDOFE+iNodeE,iS,(iX-1)*nDOFX+iNodeX) &
          = f0_Q(iNodeZ,iE,iS,iX)

      END DO
      END DO

    END DO
    END DO
    END DO

    DEALLOCATE &
      ( Mnue, Mnum, Tau_Q, f0_K, f0_Q, &
        ElectronChargedCurrentMask, MuonChargedCurrentMask )

  END SUBROUTINE ComputeEquilibriumDistributions_DG_Z


  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, f0, df0dY, df0dU )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Ye,Ym) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(out) :: f0   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: df0dY(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: df0dU(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP), ALLOCATABLE :: Me(:), dMedT(:), dMedYe(:), dMedYm(:), dMedD(:)
    REAL(DP), ALLOCATABLE :: Mm(:), dMmdT(:), dMmdYe(:), dMmdYm(:), dMmdD(:)
    REAL(DP), ALLOCATABLE :: Mp(:), dMpdT(:), dMpdYe(:), dMpdYm(:), dMpdD(:)
    REAL(DP), ALLOCATABLE :: Mn(:), dMndT(:), dMndYe(:), dMndYm(:), dMndD(:)
    REAL(DP), ALLOCATABLE :: U(:), dUdT(:), dUdYe(:), dUdYm(:), dUdD(:)
    REAL(DP) :: Mnu, dMnudT, dMnudY

    REAL(DP) :: kT, dUdL, df0dT_Y, df0dY_T, Vacancy, Susceptibility
    INTEGER  :: iE, iS, iX

    ALLOCATE( Me(iX_B:iX_E), dMedT(iX_B:iX_E), dMedYe(iX_B:iX_E), dMedYm(iX_B:iX_E), dMedD(iX_B:iX_E) )
    ALLOCATE( Mm(iX_B:iX_E), dMmdT(iX_B:iX_E), dMmdYe(iX_B:iX_E), dMmdYm(iX_B:iX_E), dMmdD(iX_B:iX_E) )
    ALLOCATE( Mp(iX_B:iX_E), dMpdT(iX_B:iX_E), dMpdYe(iX_B:iX_E), dMpdYm(iX_B:iX_E), dMpdD(iX_B:iX_E) )
    ALLOCATE( Mn(iX_B:iX_E), dMndT(iX_B:iX_E), dMndYe(iX_B:iX_E), dMndYm(iX_B:iX_E), dMndD(iX_B:iX_E) )
    ALLOCATE( U(iX_B:iX_E), dUdT(iX_B:iX_E), &
              dUdYe(iX_B:iX_E), dUdYm(iX_B:iX_E), dUdD(iX_B:iX_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Me,  dMedT,  dMedYe, dMedYm, dMedD, &
    !$OMP             Mm,  dMmdT,  dMmdYe, dMmdYm, dMmdD, &
    !$OMP             Mp,  dMpdT,  dMpdYe, dMpdYm, dMpdD, &
    !$OMP             Mn,  dMndT,  dMndYe, dMndYm, dMndD, &
    !$OMP             U, dUdT, dUdYe, dUdYm, dUdD, &
    !$OMP             f0,  df0dY, df0dU )          &
    !$OMP MAP( to: E, D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Me,  dMedT,  dMedYe, dMedYm, dMedD, &
    !$ACC         Mm,  dMmdT,  dMmdYe, dMmdYm, dMmdD, &
    !$ACC         Mp,  dMpdT,  dMpdYe, dMpdYm, dMpdD, &
    !$ACC         Mn,  dMndT,  dMndYe, dMndYm, dMndD, &
    !$ACC         U, dUdT, dUdYe, dUdYm, dUdD, &
    !$ACC         f0,  df0dY, df0dU )          &
    !$ACC COPYIN( E, D, T, Ye, Ym )
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Me, dMedD, dMedT, dMedYe, dMedYm )

    CALL ComputeMuonChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mm, dMmdD, dMmdT, dMmdYe, dMmdYm )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mp, dMpdD, dMpdT, dMpdYe, dMpdYm )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mn, dMndD, dMndT, dMndYe, dMndYm )

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Ye, Ym, U,  dUdD,  dUdT,  dUdYe,  dUdYm  )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( kT, dUdL, Mnu, dMnudT, dMnudY, &
    !$OMP          df0dT_Y, df0dY_T, Vacancy, Susceptibility )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( kT, dUdL, Mnu, dMnudT, dMnudY, &
    !$ACC          df0dT_Y, df0dY_T, Vacancy, Susceptibility )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( kT, dUdL, Mnu, dMnudT, dMnudY, &
    !$OMP          df0dT_Y, df0dY_T, Vacancy, Susceptibility )
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E

      kT = BoltzmannConstant * T(iX)

      IF ( iS == iNuE .OR. iS == iNuE_Bar ) THEN
        Mnu    = ( ( Me   (iX)  + Mp   (iX)  ) - Mn   (iX)  ) * LeptonNumber(iS)
        dMnudT = ( ( dMedT(iX)  + dMpdT(iX)  ) - dMndT(iX)  ) * LeptonNumber(iS)
        dMnudY = ( ( dMedYe(iX) + dMpdYe(iX) ) - dMndYe(iX) ) * LeptonNumber(iS)
        dUdL   = dUdYe(iX)
      ELSE IF ( iS == iNuM .OR. iS == iNuM_Bar ) THEN
        Mnu    = ( ( Mm   (iX)  + Mp   (iX)  ) - Mn   (iX)  ) * LeptonNumber(iS)
        dMnudT = ( ( dMmdT(iX)  + dMpdT(iX)  ) - dMndT(iX)  ) * LeptonNumber(iS)
        dMnudY = ( ( dMmdYm(iX) + dMpdYm(iX) ) - dMndYm(iX) ) * LeptonNumber(iS)
        dUdL   = dUdYm(iX)
      ELSE
        Mnu    = Zero
        dMnudT = Zero
        dMnudY = Zero
        dUdL   = Zero
      END IF

      ! One exponential supplies the occupation and both derivatives.
      CALL FermiDiracPair &
             ( E(iE), Mnu, kT, f0(iE,iS,iX), Vacancy, Susceptibility )
      df0dT_Y = Susceptibility * ( dMnudT + ( E(iE)-Mnu ) / T(iX) ) / kT
      df0dY_T = Susceptibility * dMnudY / kT

      df0dU(iE,iS,iX) = df0dT_Y / dUdT(iX)
      df0dY(iE,iS,iX) = df0dY_T - df0dU(iE,iS,iX) * dUdL

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Me,  dMedT,  dMedYe, dMedYm, dMedD, &
    !$OMP               Mm,  dMmdT,  dMmdYe, dMmdYm, dMmdD, &
    !$OMP               Mp,  dMpdT,  dMpdYe, dMpdYm, dMpdD, &
    !$OMP               Mn,  dMndT,  dMndYe, dMndYm, dMndD, &
    !$OMP               U, dUdT, dUdYe, dUdYm, dUdD, &
    !$OMP               E, D, T, Ye, Ym ) &
    !$OMP MAP( from: f0, df0dU, df0dY )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Me,  dMedT,  dMedYe, dMedYm, dMedD, &
    !$ACC         Mm,  dMmdT,  dMmdYe, dMmdYm, dMmdD, &
    !$ACC         Mp,  dMpdT,  dMpdYe, dMpdYm, dMpdD, &
    !$ACC         Mn,  dMndT,  dMndYe, dMndYm, dMndD, &
    !$ACC         U, dUdT, dUdYe, dUdYm, dUdD, &
    !$ACC         E, D, T, Ye, Ym ) &
    !$ACC COPYOUT( f0, df0dU, df0dY )
#endif

    DEALLOCATE( Me, dMedT, dMedYe, dMedYm, dMedD )
    DEALLOCATE( Mm, dMmdT, dMmdYe, dMmdYm, dMmdD )
    DEALLOCATE( Mp, dMpdT, dMpdYe, dMpdYm, dMpdD )
    DEALLOCATE( Mn, dMndT, dMndYe, dMndYm, dMndD )
    DEALLOCATE( U, dUdT, dUdYe, dUdYm, dUdD )

  END SUBROUTINE ComputeEquilibriumDistributionAndDerivatives


  SUBROUTINE ComputeNeutrinoOpacities_EC &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, f0, opEC, &
      Opacity_Mask, T_old, Ye_old, Ym_old, MuonMethod_Option )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL, TARGET :: T_old (iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL, TARGET :: Ye_old(iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL, TARGET :: Ym_old(iX_B:iX_E)
    INTEGER, INTENT(in), OPTIONAL :: MuonMethod_Option

    REAL(DP) :: LogE_P
    REAL(DP) :: LogD_P, LogT_P
    REAL(DP) :: E_P, D_P, T_P
    INTEGER  :: iE, iS, iX

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA             &
    !$OMP MAP( to:    E, D, T, Ye, Ym, Opacity_Mask ) &
    !$OMP MAP( alloc: opEC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA                  &
    !$ACC COPYIN( E, D, T, Ye, Ym, Opacity_Mask ) &
    !$ACC CREATE( opEC )
#endif

    IF ( iS_B <= iNuE .AND. iS_E >= iNuE_Bar ) THEN

      ! do EmAb on nucleons first
      CALL ComputeNeutrinoOpacities_EC_NuE_Nucleons &
             ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, &
               E, D, T, Ye, Ym, &
               Opacity_Mask, opEC )

      ! now add contribution from EC on heavy nuclei
      IF ( use_EC_table > 0 ) THEN

        CALL ComputeNeutrinoOpacities_EC_NuE_Nuclei &
              ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, &
                E, D, T, Ye, Ym, f0, Opacity_Mask, opEC, &
                T_old = T_old, Ye_old = Ye_old, Ym_old = Ym_old )

      END IF

    END IF

    IF ( iS_B <= iNuM .AND. iS_E >= iNuM_Bar ) THEN

      CALL ComputeNeutrinoOpacities_EC_NuM_Nucleons &
             ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, &
               E, D, T, Ye, Ym, Opacity_Mask, opEC, MuonMethod_Option )

    END IF

    IF ( iS_B <= iNuT .AND. iS_E >= iNuT_Bar ) THEN
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX = iX_B, iX_E
      DO iS = iNuT, iNuT_Bar
      DO iE = iE_B, iE_E
        opEC(iE,iS,iX) = Zero
      END DO
      END DO
      END DO
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA                &
    !$OMP MAP (release: E, D, T, Ye, Ym, Opacity_Mask ) &
    !$OMP MAP (from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA                 &
    !$ACC DELETE (E, D, T, Ye, Ym, Opacity_Mask ) &
    !$ACC COPYOUT ( opEC )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( opEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E
      opEC(iE,iS,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC


  SUBROUTINE ComputeNeutrinoOpacities_EC_NuE_Nucleons &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, &
      Opacity_Mask, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)    :: iE_B, iE_E
    INTEGER,  INTENT(in)    :: iS_B, iS_E
    INTEGER,  INTENT(in)    :: iX_B, iX_E
    REAL(DP), INTENT(in)    :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)    :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Ym(iX_B:iX_E)
    LOGICAL,  INTENT(in)    :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(inout) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: LogE_P
    REAL(DP) :: LogD_P, LogT_P
    REAL(DP) :: E_P, D_P, T_P, Ye_P
    INTEGER  :: iE, iS, iX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( LogE_P, LogD_P, LogT_P, Ye_P )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( LogE_P, LogD_P, LogT_P, Ye_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( LogE_P, LogD_P, LogT_P, Ye_P )
#endif
    DO iX = iX_B, iX_E
    DO iS = iNuE, iNuE_Bar
    DO iE = iE_B, iE_E

      IF ( Opacity_Mask(iOp_EmAb_Electron,iX) ) THEN

        LogE_P = LOG10( E(iE) / UnitE )

        LogD_P = LOG10( D(iX)  / UnitD )
        LogT_P = LOG10( T(iX)  / UnitT )
        Ye_P   =        Ye(iX) / UnitY

        CALL LogInterpolateSingleVariable_4D_Custom_Point &
                ( LogE_P , LogD_P , LogT_P , Ye_P, &
                  LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                  OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(iE,iS,iX) )

        opEC(iE,iS,iX) = opEC(iE,iS,iX) * UnitEC

      ELSE

        opEC(iE,iS,iX) = Zero

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacities_EC_NuE_Nucleons


  SUBROUTINE ComputeNeutrinoOpacities_EC_NuE_Nuclei &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, f0, &
      Opacity_Mask, opEC, T_old, Ye_old, Ym_old )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)    :: iE_B, iE_E
    INTEGER,  INTENT(in)    :: iS_B, iS_E
    INTEGER,  INTENT(in)    :: iX_B, iX_E
    REAL(DP), INTENT(in)    :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)    :: D(iX_B:iX_E)
    REAL(DP), INTENT(in), TARGET :: T(iX_B:iX_E)
    REAL(DP), INTENT(in), TARGET :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in), TARGET :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)    :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(inout) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in), TARGET, OPTIONAL :: T_old (iX_B:iX_E)
    REAL(DP), INTENT(in), TARGET, OPTIONAL :: Ye_old(iX_B:iX_E)
    REAL(DP), INTENT(in), TARGET, OPTIONAL :: Ym_old(iX_B:iX_E)

    REAL(DP), POINTER :: T0(:), Ye0(:), Ym0(:)

    REAL(DP) :: D_P, T_P, Ye_P, T0_P, Ye0_P, Ym0_P
    INTEGER  :: iE, iX

    REAL(dp) :: E_node
    REAL(dp) :: dE_node, dspec_fine, dEC_T
    REAL(dp) :: loctot
    REAL(dp) :: Xnuc, Mnu, Emissivity

    INTEGER  :: iE1, iNodeE1, iE2
    INTEGER  :: iD, iT, iY
    REAL(dp) :: dD, dT, dY
    REAL(dp) :: a, b, f_a, f_b

    REAL(dp), DIMENSION(:),   ALLOCATABLE :: Xh, Ah, EC_rate
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: spec_nodes, spec_fine
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: spec_elements, spec_elements_nodes

    IF( PRESENT ( T_old ) ) THEN
      T0 => T_old(iX_B:iX_E)
    ELSE
      T0 => T(iX_B:iX_E)
    END IF

    IF( PRESENT ( Ye_old ) ) THEN
      Ye0 => Ye_old(iX_B:iX_E)
    ELSE
      Ye0 => Ye(iX_B:iX_E)
    END IF

    IF( PRESENT ( Ym_old ) ) THEN
      Ym0 => Ym_old(iX_B:iX_E)
    ELSE
      Ym0 => Ym(iX_B:iX_E)
    END IF

    ALLOCATE( Xh     (iX_B:iX_E) )
    ALLOCATE( Ah     (iX_B:iX_E) )
    ALLOCATE( EC_rate(iX_B:iX_E) )

    ALLOCATE( spec_nodes         (iE_B:iE_E,iX_B:iX_E) )
    ALLOCATE( spec_fine          (EC_nE,    iX_B:iX_E) )
    ALLOCATE( spec_elements      (nE,       iX_B:iX_E) )
    ALLOCATE( spec_elements_nodes(nE,       iX_B:iX_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: f0 ) &
    !$OMP MAP( alloc: T0, Ye0, Ym0, &
    !$OMP             Xh, Ah, EC_rate, & 
    !$OMP             spec_nodes, spec_elements,      &
    !$OMP             spec_elements_nodes, spec_fine )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( f0 ) &    
    !$ACC CREATE( T0, Ye0, Ym0, &
    !$ACC         Xh, Ah, EC_rate, &
    !$ACC         spec_nodes, spec_elements, &
    !$ACC         spec_elements_nodes, spec_fine )
#endif

    ! --- Compute heavy mass fraction and number ---

    CALL ComputeHeavyMassFraction_TABLE( D, T , Ye , Ym , Xh )
    CALL ComputeHeavyMassNumber_TABLE  ( D, T , Ye , Ym , Ah )

    ASSOCIATE &
      ( CenterE => MeshE % Center, &
        WidthE  => MeshE % Width, &
        NodesE  => MeshE % Nodes )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: CenterE, WidthE, NodesE ) &
    !$OMP PRIVATE( D_P, T_P, Ye_P, T0_P, Ye0_P, Xnuc, loctot, & 
    !$OMP          dspec_fine, dEC_T, dE_node, E_node, Mnu, Emissivity, &
    !$OMP          iE, iE1, iE2, iNodeE1, a, b, f_a, f_b, & 
    !$OMP          iD, iT, iY, dD, dT, dY )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN ( CenterE, WidthE, NodesE ) &
    !$ACC PRIVATE( D_P, T_P, Ye_P, T0_P, Ye0_P, Xnuc, loctot, & 
    !$ACC          dspec_fine, dEC_T, dE_node, E_node, Mnu, Emissivity, &
    !$ACC          iE, iE1, iE2, iNodeE1, a, b, f_a, f_b, & 
    !$ACC          iD, iT, iY, dD, dT, dY ) & 
    !$ACC PRESENT( Xh, Ah, T0, Ye0, EC_rate, OS_EmAb_EC_rate, &
    !$ACC          OS_EmAb_EC_spec, EmAb_EC_spec_T, WeightsE,         &
    !$ACC          f0, EC_nE, EC_dE, EC_iE_max, EC_iNodeE_max,        &
    !$ACC          Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, T_P, Ye_P, T0_P, Ye0_P, Xnuc, loctot, & 
    !$OMP          dspec_fine, dEC_T, dE_node, E_node, Mnu, Emissivity, &
    !$OMP          iE, iE1, iE2, iNodeE1, a, b, f_a, f_b, & 
    !$OMP          iD, iT, iY, dD, dT, dY )
#endif
    DO iX = iX_B, iX_E

      IF ( Opacity_Mask(iOp_ECTable,iX) )THEN

        D_P  = D (iX) / UnitD
        T_P  = T (iX) / UnitT
        Ye_P = Ye(iX) / UnitY

        ! If T_old and Ye_old were present, use the old state to
        ! calculate the interpolation indices, otherwise the
        ! current state is used, ie T0 = T, Ye0 = Ye
        T0_P  = T0 (iX) / UnitT
        Ye0_P = Ye0(iX) / UnitY

        ! Determine indices for D, T, Ye using old state
        iD = Index1D_Log( D_P  , Ds_EC_T )
        iT = Index1D_Log( T0_P , Ts_EC_T )
        iY = Index1D_Lin( Ye0_P, Ys_EC_T )

        ! Use new state to calculate interpolation weights
        dD = LOG10( D_P / Ds_EC_T(iD) ) / LOG10( Ds_EC_T(iD+1) / Ds_EC_T(iD) )
        dT = LOG10( T_P / Ts_EC_T(iT) ) / LOG10( Ts_EC_T(iT+1) / Ts_EC_T(iT) )
        dY = ( Ye_P - Ys_EC_T(iY) ) / ( Ys_EC_T(iY+1) - Ys_EC_T(iY) )

        CALL LinearInterp3D_3DArray_Point &
             ( iD, iT, iY, dD, dT, dY, &
               OS_EmAb_EC_rate(iNuE), EmAb_EC_rate_T, &
               EC_rate(iX) )

        loctot = 0.0d0
        DO iE = 1, EC_nE

          CALL LinearInterp3D_3DArray_Point &
              ( iD, iT, iY, dD, dT, dY, &
                OS_EmAb_EC_spec(iNuE), EmAb_EC_spec_T(:,:,:,iE), &
                spec_fine(iE,iX) )

          loctot = loctot + spec_fine(iE,iX) * EC_dE

        END DO

        DO iE = iE_B, iE_E
          spec_nodes(iE,iX) = 0.0d0
        END DO

        DO iE = 1, nE
          spec_elements(iE,iX) = 0.0d0
          spec_elements_nodes(iE,iX) = 0.0d0
        END DO

        !renormalise interpolated fine energy grid spectrum to 1
        DO iE = 1, EC_nE
          spec_fine(iE,iX) = spec_fine(iE,iX) / loctot
        END DO
      
        !calculate spectrum in elements from interpolated fine spectrum 
        DO iE1 = 1, EC_iE_max
          loctot = 0.0d0
          DO iNodeE1 = 1, nNodesE
            iE = (iE1-1) * nNodesE + iNodeE1
            IF ( iE <= EC_iNodeE_max ) THEN
              E_node = NodeCoordinate( CenterE(iE1), WidthE(iE1), NodesE(iNodeE1) ) / UnitE

              iE2 = FLOOR(E_node/EC_dE) + 1

              dspec_fine = spec_fine(iE2+1,iX) - spec_fine(iE2,iX)
              dEC_T      = Es_EC_T  (iE2+1   ) - Es_EC_T  (iE2   ) ! AH: is this just EC_dE?
              dE_node    = E_node              - Es_EC_T  (iE2   )

              spec_nodes(iE,iX) = spec_fine(iE2,iX) + dspec_fine * dE_node / dEC_T

              loctot = loctot &
                     + spec_nodes(iE,iX) * ( WidthE(iE1) / UnitE ) * WeightsE(iNodeE1)

              spec_nodes(iE,iX) = spec_nodes(iE,iX) / E_node**2
            END IF
          END DO
          spec_elements_nodes(iE1,iX) = loctot
        END DO

        !integrate fine spectrum over elements for renormlisation
        DO iE = 1, EC_iE_max

          IF ( EC_kfmin(iE) >= EC_kfmax(iE) ) THEN

            a = EC_a(iE,1)
            b = EC_b(iE,2)

            f_a =  spec_fine(EC_ak(iE,1),iX) + (a-Es_EC_T(EC_ak(iE,1)))    &
                * (spec_fine(EC_ak(iE,2),iX) -  spec_fine(EC_ak(iE,1),iX)) &
                / (  Es_EC_T(EC_ak(iE,2))    -    Es_EC_T(EC_ak(iE,1)))
            f_b =  spec_fine(EC_bk(iE,1),iX) + (b-Es_EC_T(EC_bk(iE,1)))    &
                * (spec_fine(EC_bk(iE,2),iX) -  spec_fine(EC_bk(iE,1),iX)) &
                / (  Es_EC_T(EC_bk(iE,2))    -    Es_EC_T(EC_bk(iE,1)))
 
            spec_elements(iE,iX) = spec_elements(iE,iX) &
                                 + 0.5d0 * (b - a) * (f_a + f_b)
          ELSE
            
            !integrate from lower thornado element face to next
            !EC table energy bin face  
            a = EC_a(iE,1)
            b = EC_b(iE,1)

            f_a =  spec_fine(EC_ak(iE,1),iX) + (a-Es_EC_T(EC_ak(iE,1)))    &
                * (spec_fine(EC_ak(iE,2),iX) -  spec_fine(EC_ak(iE,1),iX)) &
                / (  Es_EC_T(EC_ak(iE,2))    -    Es_EC_T(EC_ak(iE,1)))

            f_b = spec_fine(EC_kfmin(iE)+1,iX)

            spec_elements(iE,iX) = spec_elements(iE,iX) &
                                 + 0.5d0 * (b - a) * (f_a + f_b)

            !integrate over the fully contained EC table energy bins
            DO iE1 = EC_kfmin(iE)+1, EC_kfmax(iE)-1
              spec_elements(iE,iX) = spec_elements(iE,iX) &
                                   + 0.5d0 * EC_dE * (spec_fine(iE1,iX) + spec_fine(iE1+1,iX))
            END DO

            !integrate from EC table energy bin face
            !to upper thornado element face  
            a = EC_a(iE,2)            
            b = EC_b(iE,2)

            f_a = spec_fine(EC_kfmax(iE),iX)            

            f_b =  spec_fine(EC_bk(iE,1),iX) + (b-Es_EC_T(EC_bk(iE,1)))    &
                * (spec_fine(EC_bk(iE,2),iX) -  spec_fine(EC_bk(iE,1),iX)) &
                / (  Es_EC_T(EC_bk(iE,2))    -    Es_EC_T(EC_bk(iE,1)))

            spec_elements(iE,iX) = spec_elements(iE,iX) &
                                 + 0.5d0 * (b - a) * (f_a + f_b)            

          END IF

        END DO

        Xnuc = Xh(iX) / (Ah(iX) / AvogadroConstantMKS) * D_P 

        !renormalise spectrum at the nodes so that it integrates to 1
        DO iE = iE_B, EC_iNodeE_max
          iE1     = MOD( (iE-1) / nNodesE, nE ) + 1

          spec_nodes(iE,iX) = spec_nodes(iE,iX) * spec_elements(iE1,iX) &
                            / spec_elements_nodes(iE1,iX)

        END DO

        DO iE = iE_B, iE_E

          !now add Chi from EC table on heavy nuclei to Chi from EmAb on nucleons
          Emissivity = Xnuc * coeff * EC_rate(iX) * spec_nodes(iE,iX)
          IF( Emissivity == Zero )THEN
            CYCLE
          ELSE IF( f0(iE,iNuE,iX) >= TINY(One) )THEN
            opEC(iE,iNuE,iX) = opEC(iE,iNuE,iX) &
                                + Emissivity / f0(iE,iNuE,iX) * UnitEC
          ELSE
            ! Recover j/f from its exponent, not a fabricated occupation floor.
            CALL ComputeElectronNeutrinoChemicalPotential_TABLE &
                   ( D(iX), T(iX), Ye(iX), Ym(iX), Mnu )
            opEC(iE,iNuE,iX) = opEC(iE,iNuE,iX) &
              + DivideByFermiDirac &
                  ( Emissivity, ( E(iE)-Mnu ) / ( BoltzmannConstant*T(iX) ) ) * UnitEC
          END IF

        END DO

      END IF

    END DO

    END ASSOCIATE

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: f0, T0, Ye0, Ym0, &
    !$OMP               Xh, Ah, EC_rate, &
    !$OMP               spec_nodes, spec_elements, &
    !$OMP               spec_elements_nodes, spec_fine )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( f0, T0, Ye0, Ym0, &
    !$ACC         Xh, Ah, EC_rate, &
    !$ACC         spec_nodes, spec_elements, &
    !$ACC         spec_elements_nodes, spec_fine )
#endif

    DEALLOCATE( Xh      )
    DEALLOCATE( Ah      )
    DEALLOCATE( EC_rate )

    DEALLOCATE( spec_nodes          )
    DEALLOCATE( spec_fine           )
    DEALLOCATE( spec_elements       )
    DEALLOCATE( spec_elements_nodes )

  END SUBROUTINE ComputeNeutrinoOpacities_EC_NuE_Nuclei


  SUBROUTINE ComputeNeutrinoOpacities_EC_NuM_Nucleons &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, &
      Opacity_Mask, opEC, MuonMethod_Option )

    ! --- Full-kinematics muon charged-current opacities ---

    INTEGER,  INTENT(in)    :: iE_B, iE_E
    INTEGER,  INTENT(in)    :: iS_B, iS_E
    INTEGER,  INTENT(in)    :: iX_B, iX_E
    REAL(DP), INTENT(in)    :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)    :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Ym(iX_B:iX_E)
    LOGICAL,  INTENT(in)    :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(inout) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    INTEGER, INTENT(in), OPTIONAL :: MuonMethod_Option
    INTEGER :: LocalMethod
    INTEGER  :: iE, iS, iX

    REAL(DP), ALLOCATABLE :: Mumu(:)
    REAL(DP), ALLOCATABLE :: Mun(:), Mn_eff(:), Un(:)
    REAL(DP), ALLOCATABLE :: Mup(:), Mp_eff(:), Up(:)

#ifndef EOSMODE_COMPOSE
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iS = iNuM, iNuM_Bar
    DO iE = iE_B, iE_E
      opEC(iE,iS,iX) = Zero
    END DO
    END DO
    END DO
    RETURN
#endif

    IF( .NOT. ANY( Opacity_Mask(iOp_EmAb_Muon,iX_B:iX_E) ) )THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX = iX_B, iX_E
      DO iS = iNuM, iNuM_Bar
      DO iE = iE_B, iE_E
        opEC(iE,iS,iX) = Zero
      END DO
      END DO
      END DO
      RETURN
    END IF

    ALLOCATE( Mumu  (iX_B:iX_E) )
    ALLOCATE( Mun   (iX_B:iX_E) )
    ALLOCATE( Mn_eff(iX_B:iX_E) )
    ALLOCATE( Un    (iX_B:iX_E) )
    ALLOCATE( Mup   (iX_B:iX_E) )
    ALLOCATE( Mp_eff(iX_B:iX_E) )
    ALLOCATE( Up    (iX_B:iX_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Mun, Mn_eff, Un, Mup, Mp_eff, Up, Mumu )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA  &
    !$ACC CREATE( Mun, Mn_eff, Un, Mup, Mp_eff, Up, Mumu )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( D, T, Ye, Ym, Opacity_Mask, &
    !$ACC          Mumu, Mun, Mn_eff, Un, Mup, Mp_eff, Up )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E
      IF( Opacity_Mask(iOp_EmAb_Muon,iX) )THEN
        CALL ComputeMuonChemicalPotential_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Mumu(iX) )
        CALL ComputeNeutronChemicalPotential_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Mun(iX) )
        CALL ComputeNeutronEffectiveMass_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Mn_eff(iX) )
        CALL ComputeNeutronSelfEnergy_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Un(iX) )
        CALL ComputeProtonChemicalPotential_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Mup(iX) )
        CALL ComputeProtonEffectiveMass_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Mp_eff(iX) )
        CALL ComputeProtonSelfEnergy_TABLE &
               ( D(iX), T(iX), Ye(iX), Ym(iX), Up(iX) )
      ELSE
        Mumu  (iX) = Zero
        Mun   (iX) = Zero
        Mn_eff(iX) = Zero
        Un    (iX) = Zero
        Mup   (iX) = Zero
        Mp_eff(iX) = Zero
        Up    (iX) = Zero
      END IF
    END DO

    LocalMethod = EmAb_Muon_Method
    IF( PRESENT(MuonMethod_Option) ) LocalMethod = MuonMethod_Option
    SELECT CASE( LocalMethod )
    CASE( 1 )
      CALL ComputeNeutrinoOpacities_EC_NuM_Nucleons_Elastic &
             ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, &
               E, D, T, Ye, Ym, Mumu, &
               Mun, Mn_eff, Un, Mup, Mp_eff, Up, &
               Opacity_Mask, opEC )
    CASE( 2 )
      CALL ComputeNeutrinoOpacities_EC_NuM_Nucleons_2D &
             ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, &
               E, T, Mumu, &
               Mun, Mn_eff, Un, Mup, Mp_eff, Up, &
               Opacity_Mask, opEC )
    END SELECT

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Mun, Mn_eff, Un, Mup, Mp_eff, Up, Mumu )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Mun, Mn_eff, Un, Mup, Mp_eff, Up, Mumu )
#endif

    DEALLOCATE( Mumu   )
    DEALLOCATE( Mun    )
    DEALLOCATE( Mn_eff )
    DEALLOCATE( Un     )
    DEALLOCATE( Mup    )
    DEALLOCATE( Mp_eff )
    DEALLOCATE( Up     )

  END SUBROUTINE ComputeNeutrinoOpacities_EC_NuM_Nucleons


  SUBROUTINE ComputeNeutrinoOpacities_EC_NuM_Nucleons_Elastic &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Ye, Ym, &
      Mumu, Mun, Mn_eff, Un, Mup, Mp_eff, Up, Opacity_Mask, opEC )

    ! --- Muon Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)    :: iE_B, iE_E
    INTEGER,  INTENT(in)    :: iS_B, iS_E
    INTEGER,  INTENT(in)    :: iX_B, iX_E
    REAL(DP), INTENT(in)    :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)    :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mumu(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mun(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mn_eff(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Un(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mup(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mp_eff(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Up(iX_B:iX_E)
    LOGICAL,  INTENT(in)    :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(inout) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: E_P, D_P, T_P
    INTEGER  :: iE, iS, iX

    REAL(DP), ALLOCATABLE :: Xn(:), Xp(:)

    REAL(DP) :: Mumu_P
    REAL(DP) :: Mun_P, Mn_eff_P, Un_P, Xn_P
    REAL(DP) :: Mup_P, Mp_eff_P, Up_P, Xp_P
    REAL(DP) :: WeakMagRecNuLep
    REAL(DP) :: WeakMagRecNuLepBar
    REAL(DP) :: Mnum, kT

    ALLOCATE( Xn(iX_B:iX_E) )
    ALLOCATE( Xp(iX_B:iX_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Xn, Xp )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA  &
    !$ACC CREATE( Xn, Xp )
#endif

    CALL ComputeNeutronMassFraction_TABLE( D, T, Ye, Ym, Xn )
    CALL ComputeProtonMassFraction_TABLE ( D, T, Ye, Ym, Xp )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Mun_P, Mn_eff_P, Un_P, Xn_P, &
    !$OMP          Mup_P, Mp_eff_P, Up_P, Xp_P, Mumu_P, &
    !$OMP          E_P, D_P, T_P, WeakMagRecNuLep, WeakMagRecNuLepBar, &
    !$OMP          Mnum, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Mun_P, Mn_eff_P, Un_P, Xn_P, &
    !$ACC          Mup_P, Mp_eff_P, Up_P, Xp_P, Mumu_P, &
    !$ACC          E_P, D_P, T_P, WeakMagRecNuLep, WeakMagRecNuLepBar, &
    !$ACC          Mnum, kT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( Mun_P, Mn_eff_P, Un_P, Xn_P, &
    !$OMP          Mup_P, Mp_eff_P, Up_P, Xp_P, Mumu_P, &
    !$OMP          E_P, D_P, T_P, WeakMagRecNuLep, WeakMagRecNuLepBar, &
    !$OMP          Mnum, kT )
#endif
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E

      IF ( Opacity_Mask(iOp_EmAb_Muon,iX) ) THEN

        E_P  = E(iE)  / UnitE
        D_P  = D(iX)  / UnitD
        T_P  = T(iX)  / UnitT
        Mun_P    = Mun   (iX) / UnitMn
        Mn_eff_P = Mn_eff(iX) / UnitMn
        Un_P     = Un    (iX) / UnitMn
        Xn_P     = Xn    (iX)
        Mup_P    = Mup   (iX) / UnitMp
        Mp_eff_P = Mp_eff(iX) / UnitMp
        Up_P     = Up    (iX) / UnitMp
        Xp_P     = Xp    (iX)
        Mumu_P   = Mumu  (iX) / UnitMmu

        CALL ElasticAbsorptionOpacityNum &
            ( D_P, T_P, E_P, &
              Mun_P, Mn_eff_P, Un_P, Xn_P, &
              Mup_P, Mp_eff_P, Up_P, Xp_P, &
              Mumu_P, opEC(iE,iNuM,iX), opEC(iE,iNuM_Bar,iX) )

        CALL CalculateHorowitzWeakMagRecoil(E_P, WeakMagRecNuLep, WeakMagRecNuLepBar)

        Mnum = ( Mumu(iX) + Mup(iX) ) - Mun(iX)
        kT   = BoltzmannConstant * T(iX)
        opEC(iE,iNuM,iX) &
          = DivideByFermiDirac &
              ( opEC(iE,iNuM,iX) * WeakMagRecNuLep, &
                ( Mnum-E(iE) ) / kT ) * UnitEC
        opEC(iE,iNuM_Bar,iX) &
          = DivideByFermiDirac &
              ( opEC(iE,iNuM_Bar,iX) * WeakMagRecNuLepBar, &
                ( -Mnum-E(iE) ) / kT ) * UnitEC

      ELSE

        opEC(iE,iNuM    ,iX) = Zero
        opEC(iE,iNuM_Bar,iX) = Zero

      END IF

    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Xn, Xp )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Xn, Xp )
#endif

    DEALLOCATE( Xn )
    DEALLOCATE( Xp )

  END SUBROUTINE ComputeNeutrinoOpacities_EC_NuM_Nucleons_Elastic


  SUBROUTINE ComputeNeutrinoOpacities_EC_NuM_Nucleons_2D &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, T, &
      Mumu, Mun, Mn_eff, Un, Mup, Mp_eff, Up, Opacity_Mask, opEC )

    ! --- Muon Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)    :: iE_B, iE_E
    INTEGER,  INTENT(in)    :: iS_B, iS_E
    INTEGER,  INTENT(in)    :: iX_B, iX_E
    REAL(DP), INTENT(in)    :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)    :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mumu(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mun(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mn_eff(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Un(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mup(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Mp_eff(iX_B:iX_E)
    REAL(DP), INTENT(in)    :: Up(iX_B:iX_E)
    LOGICAL,  INTENT(in)    :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(inout) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: E_P, T_P
    INTEGER  :: iE, iX

    REAL(DP) :: Mumu_P
    REAL(DP) :: Mun_P, Mn_eff_P, Un_P
    REAL(DP) :: Mup_P, Mp_eff_P, Up_P
    REAL(DP) :: ReactionOpacity, NodeContribution
    REAL(DP) :: Mnum, kT
    INTEGER  :: iReaction, iQ1, iQ2
    TYPE(CC2D_CallData) :: CallData
    TYPE(CC2D_RowData)  :: RowData

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
    !$OMP PRIVATE( Mun_P, Mn_eff_P, Un_P, Mup_P, Mp_eff_P, Up_P, Mumu_P, &
    !$OMP          E_P, T_P, ReactionOpacity, NodeContribution, &
    !$OMP          iReaction, iQ1, iQ2, CallData, RowData, &
    !$OMP          Mnum, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG COLLAPSE(2) &
    !$ACC PRIVATE( Mun_P, Mn_eff_P, Un_P, Mup_P, Mp_eff_P, Up_P, Mumu_P, &
    !$ACC          E_P, T_P, ReactionOpacity, NodeContribution, &
    !$ACC          iReaction, iQ1, iQ2, CallData, RowData, &
    !$ACC          Mnum, kT ) &
    !$ACC PRESENT( MuonCC_QuadratureNodes, MuonCC_QuadratureWeights )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( Mun_P, Mn_eff_P, Un_P, Mup_P, Mp_eff_P, Up_P, Mumu_P, &
    !$OMP          E_P, T_P, ReactionOpacity, NodeContribution, &
    !$OMP          iReaction, iQ1, iQ2, CallData, RowData, &
    !$OMP          Mnum, kT )
#endif
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E

      IF ( Opacity_Mask(iOp_EmAb_Muon,iX) ) THEN

        E_P  = E(iE)  / UnitE
        T_P  = T(iX)  / UnitT * kmev

        Mun_P  = Mun (iX) / UnitMn
        Mup_P  = Mup (iX) / UnitMp
        Mumu_P = Mumu(iX) / UnitMmu
        Mn_eff_P = Mn_eff(iX) / UnitMn
        Mp_eff_P = Mp_eff(iX) / UnitMp
        Un_P     = Un(iX) / UnitMn
        Up_P     = Up(iX) / UnitMp

        opEC(iE,iNuM    ,iX) = Zero
        opEC(iE,iNuM_Bar,iX) = Zero

        DO iReaction = 1, 3

          CALL CC2D_PrepareCall &
                 ( iReaction, E_P, T_P, Mumu_P, Mun_P, Mup_P, &
                   mmu, Mn_eff_P, Mp_eff_P, Un_P, Up_P, &
                   MuonCCCorrections, CallData )

          IF( .NOT. CallData % HasPhaseSpace ) CYCLE

          ReactionOpacity = Zero
#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO SIMD REDUCTION(+:ReactionOpacity) &
          !$OMP PRIVATE( RowData, iQ2, NodeContribution )
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR REDUCTION(+:ReactionOpacity) &
          !$ACC PRIVATE( RowData, iQ2, NodeContribution )
#endif
          DO iQ1 = 1, nE_Muon_2D
            CALL CC2D_RowSetup &
                   ( CallData, MuonCC_QuadratureNodes(iQ1), &
                     MuonCC_QuadratureWeights(iQ1), RowData )
            IF( RowData % HasWindow )THEN
              DO iQ2 = 1, nE_Muon_2D
                CALL CC2D_NodeKernel &
                       ( CallData, RowData, MuonCC_QuadratureNodes(iQ2), &
                         MuonCC_QuadratureWeights(iQ2), NodeContribution )
                ReactionOpacity = ReactionOpacity + NodeContribution
              END DO
            END IF
          END DO
          CALL CC2D_ApplyPrefactor( CallData, ReactionOpacity )

          ! Reactions 1--3 are nu_mu absorption, anti-nu_mu absorption,
          ! and inverse neutron decay.  The corresponding inverse
          ! emissivities follow from detailed balance in the relaxation
          ! source Chi_EmAb * ( f_eq - f ).  Consequently opEC stores the
          ! corrected absorption opacity chi_abs / ( 1 - f_eq ), consistent
          ! with the electron charged-current opacity tables.
          IF( iReaction == 1 )THEN
            opEC(iE,iNuM,iX) = ReactionOpacity
          ELSE
            opEC(iE,iNuM_Bar,iX) &
              = opEC(iE,iNuM_Bar,iX) + ReactionOpacity
          END IF

        END DO

        Mnum = ( Mumu(iX) + Mup(iX) ) - Mun(iX)
        kT   = BoltzmannConstant * T(iX)
        opEC(iE,iNuM,iX) &
          = DivideByFermiDirac &
              ( opEC(iE,iNuM,iX), ( Mnum-E(iE) ) / kT ) * UnitEC
        opEC(iE,iNuM_Bar,iX) &
          = DivideByFermiDirac &
              ( opEC(iE,iNuM_Bar,iX), ( -Mnum-E(iE) ) / kT ) * UnitEC

      ELSE

        opEC(iE,iNuM    ,iX) = Zero
        opEC(iE,iNuM_Bar,iX) = Zero

      END IF

    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacities_EC_NuM_Nucleons_2D


  SUBROUTINE WriteMuonChargedCurrentOpacityDiagnostic &
    ( iFile, E, D, T, Ye, Ym, f0 )

    INTEGER,  INTENT(in) :: iFile
    REAL(DP), INTENT(in) :: E(:)
    REAL(DP), INTENT(in) :: D, T, Ye, Ym
    REAL(DP), INTENT(in) :: f0(:,:)

#ifdef MICROPHYSICS_WEAKLIB
    REAL(DP) :: Mumu, Mun, Mn_eff, Un
    REAL(DP) :: Mup, Mp_eff, Up
    REAL(DP) :: D_P, T_P, E_P
    REAL(DP) :: Mumu_P, Mun_P, Mn_eff_P, Un_P
    REAL(DP) :: Mup_P, Mp_eff_P, Up_P
    REAL(DP) :: RawIntegral, FinalOpacity, CorrectedOpacity
    REAL(DP) :: Mnum, kT, Exponent
    REAL(DP) :: NodeContribution
    INTEGER  :: iE, iReaction, iQ1, iQ2, iS, nRows
    TYPE(CC2D_CallData) :: CallData
    TYPE(CC2D_RowData)  :: RowData

    IF( EmAb_Muon_Method /= 2 )THEN
      WRITE(iFile,'(/,A,I0,A)') &
        '# detailed muon charged-current opacity diagnostic omitted: method ', &
        EmAb_Muon_Method, ' does not use the 2D quadrature'
      RETURN
    END IF

    IF( .NOT. ALLOCATED(MuonCC_QuadratureNodes) &
        .OR. .NOT. ALLOCATED(MuonCC_QuadratureWeights) )THEN
      WRITE(iFile,'(/,A)') &
        '# detailed muon charged-current opacity diagnostic omitted: ' // &
        '2D quadrature is not allocated'
      RETURN
    END IF

    CALL ComputeMuonChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mumu )
    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mun )
    CALL ComputeNeutronEffectiveMass_TABLE &
           ( D, T, Ye, Ym, Mn_eff )
    CALL ComputeNeutronSelfEnergy_TABLE &
           ( D, T, Ye, Ym, Un )
    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mup )
    CALL ComputeProtonEffectiveMass_TABLE &
           ( D, T, Ye, Ym, Mp_eff )
    CALL ComputeProtonSelfEnergy_TABLE &
           ( D, T, Ye, Ym, Up )

    D_P      = D      / UnitD
    T_P      = T      / UnitT * kmev
    Mumu_P   = Mumu   / UnitMmu
    Mun_P    = Mun    / UnitMn
    Mn_eff_P = Mn_eff / UnitMn
    Un_P     = Un     / UnitMn
    Mup_P    = Mup    / UnitMp
    Mp_eff_P = Mp_eff / UnitMp
    Up_P     = Up     / UnitMp
    Mnum     = ( Mumu + Mup ) - Mun
    kT       = BoltzmannConstant * T

    WRITE(iFile,'(/,A)') '# muon charged-current opacity diagnostic'
    WRITE(iFile,'(A,5ES26.16E3)') &
      '# rho_gcm3 T_MeV Ye Ym mu_mu_MeV: ', D_P, T_P, Ye, Ym, Mumu_P
    WRITE(iFile,'(A,6ES26.16E3)') &
      '# mu_n mstar_n U_n mu_p mstar_p U_p [MeV]: ', &
      Mun_P, Mn_eff_P, Un_P, Mup_P, Mp_eff_P, Up_P
    WRITE(iFile,'(A)') &
      '# iE reaction E_MeV phase_space E2_min E2_max active_rows ' // &
      'raw_integral final_cm^-1 corrected_cm^-1'

    DO iE = 1, SIZE(E)
      E_P = E(iE) / UnitE
      DO iReaction = 1, 3

        CALL CC2D_PrepareCall &
               ( iReaction, E_P, T_P, Mumu_P, Mun_P, Mup_P, &
                 mmu, Mn_eff_P, Mp_eff_P, Un_P, Up_P, &
                 MuonCCCorrections, CallData )

        RawIntegral = Zero
        nRows        = 0
        IF( CallData % HasPhaseSpace )THEN
          DO iQ1 = 1, nE_Muon_2D
            CALL CC2D_RowSetup &
                   ( CallData, MuonCC_QuadratureNodes(iQ1), &
                     MuonCC_QuadratureWeights(iQ1), RowData )
            IF( .NOT. RowData % HasWindow ) CYCLE
            nRows = nRows + 1
            DO iQ2 = 1, nE_Muon_2D
              CALL CC2D_NodeKernel &
                     ( CallData, RowData, MuonCC_QuadratureNodes(iQ2), &
                       MuonCC_QuadratureWeights(iQ2), NodeContribution )
              RawIntegral = RawIntegral + NodeContribution
            END DO
          END DO
        END IF

        FinalOpacity = RawIntegral
        CALL CC2D_ApplyPrefactor( CallData, FinalOpacity )
        IF( iReaction == 1 )THEN
          iS = iNuM
          Exponent = ( Mnum-E(iE) ) / kT
        ELSE
          iS = iNuM_Bar
          Exponent = ( -Mnum-E(iE) ) / kT
        END IF
        CorrectedOpacity = DivideByFermiDirac( FinalOpacity, Exponent )

        WRITE(iFile,'(2I6,ES26.16E3,1X,L1,2ES26.16E3,I8,3ES26.16E3)') &
          iE, iReaction, E_P, CallData % HasPhaseSpace, &
          CallData % E2_min, CallData % E2_max, nRows, &
          RawIntegral, FinalOpacity, CorrectedOpacity

      END DO
    END DO
#else
    WRITE(iFile,'(/,A)') &
      '# muon charged-current opacity diagnostic unavailable without WeakLib'
#endif

  END SUBROUTINE WriteMuonChargedCurrentOpacityDiagnostic

  SUBROUTINE ComputeNeutrinoOpacities_EC_Vector &
    ( iP_B, iP_E, iS_B, iS_E, E, D, T, Ye, Ym, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E (iP_B:iP_E)
    REAL(DP), INTENT(in)  :: D (iP_B:iP_E)
    REAL(DP), INTENT(in)  :: T (iP_B:iP_E)
    REAL(DP), INTENT(in)  :: Ye(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: Ym(iP_B:iP_E)
    REAL(DP), INTENT(out) :: opEC(iP_B:iP_E,iS_B:iS_E)

    REAL(DP), ALLOCATABLE :: LogE_P(:)
    REAL(DP), ALLOCATABLE :: LogD_P(:), LogT_P(:), Y_P(:)
    INTEGER  :: iP, iS

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( LogE_P(iP_B:iP_E) )
    ALLOCATE( LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P, opEC ) &
    !$OMP MAP( to: E, D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P, opEC ) &
    !$ACC COPYIN( E, D, T, Ye, Ym )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP)      / UnitD )
      LogT_P(iP) = LOG10( T(iP)      / UnitT )
      Y_P   (iP) = (Ye(iP) + Ym(iP)) / UnitY
      LogE_P(iP) = LOG10( E(iP)      / UnitE )
    END DO

    DO iS = iS_B, iS_E

      CALL LogInterpolateSingleVariable_4D_Custom &
             ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(:,iS) )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS = iS_B, iS_E
    DO iP = iP_B, iP_E
      IF ( QueryOpacity_EmAb( D(iP) / UnitD ) ) THEN
        opEC(iP,iS) = opEC(iP,iS) * UnitEC
      ELSE
        opEC(iP,iS) = Zero
      END IF
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Ye, Ym ) &
    !$OMP MAP( from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Ye, Ym ) &
    !$ACC COPYOUT( opEC )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP MAP( from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC COPYOUT( opEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS = iS_B, iS_E
    DO iP = iP_B, iP_E
      opEC(iP,iS) = Zero
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Vector


  SUBROUTINE ComputeNeutrinoOpacities_ES &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Ye, Ym, iMoment, opES, &
      Opacity_Mask )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(iE_B:iE_E,iNu:iNu_Bar,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)

    REAL(DP) :: LogE_P, LogD_P, LogT_P, Yp_P
    INTEGER  :: iX, iE, iC

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA          &
    !$OMP MAP( to:    E, D, T, Ye, Ym, &
    !$OMP             Opacity_Mask ) &
    !$OMP MAP( alloc: opES )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA                 &
    !$ACC COPYIN  ( E, D, T, Ye, Ym, &
    !$ACC           Opacity_Mask )   &
    !$ACC CREATE  ( opES )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( LogE_P, LogD_P, LogT_P, Yp_P )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( LogEs_T, LogDs_T, LogTs_T, Ys_T, &
    !$ACC          OS_Iso, Iso_T, Opacity_Mask ) &
    !$ACC PRIVATE( LogE_P, LogD_P, LogT_P, Yp_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( LogE_P, LogD_P, LogT_P, Yp_P )
#endif
    DO iX = iX_B, iX_E
    DO iC = iNu,  iNu_Bar
    DO iE = iE_B, iE_E

      IF ( Opacity_Mask(iOp_Iso,iX) ) THEN

        LogE_P = LOG10( E(iE) / UnitE )

        LogD_P = LOG10( D (iX) / UnitD )
        LogT_P = LOG10( T (iX) / UnitT )
        Yp_P   = (Ye(iX) + Ym(iX)) / UnitY

        CALL LogInterpolateSingleVariable_4D_Custom_Point &
               ( LogE_P , LogD_P , LogT_P , Yp_P, &
                 LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                 OS_Iso(iC,iMoment), Iso_T(:,:,:,:,iMoment,iC), opES(iE,iC,iX) )

        opES(iE,iC,iX) = opES(iE,iC,iX) * UnitES

      ELSE

        opES(iE,iC,iX) = Zero

      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA             &
    !$OMP MAP (release: E, D, T, Ye, Ym, &
    !$OMP               Opacity_Mask ) &
    !$OMP MAP (from: opES )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA                    &
    !$ACC DELETE (      E, D, T, Ye, Ym, &
    !$ACC               Opacity_Mask ) &
    !$ACC COPYOUT ( opES )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( opES )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iC = iNu,  iNu_Bar
    DO iE = iE_B, iE_E
      opES(iE,iC,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES

  
  SUBROUTINE ComputeNeutrinoOpacities_ES_Vector &
    ( iP_B, iP_E, E, D, T, Ye, Ym, iMoment, opES )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    REAL(DP), INTENT(in)  :: E(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: D(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: T(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: Ye(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: Ym(iP_B:iP_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(iP_B:iP_E)

    REAL(DP), ALLOCATABLE :: LogE_P(:)
    REAL(DP), ALLOCATABLE :: LogD_P(:), LogT_P(:), Yp_P(:)
    INTEGER  :: iP

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( LogE_P(iP_B:iP_E) )
    ALLOCATE( LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Yp_P(iP_B:iP_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Yp_P, opES ) &
    !$OMP MAP( to: E, D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Yp_P, opES ) &
    !$ACC COPYIN( E, D, T, Ye, Ym )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Yp_P(iP)   = (Ye(iP) + Ym(iP)) / UnitY
      LogE_P(iP) = LOG10( E(iP) / UnitE )
    END DO

    CALL LogInterpolateSingleVariable_4D_Custom &
           ( LogE_P, LogD_P, LogT_P, Yp_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(1,iMoment), Iso_T(:,:,:,:,iMoment,1), opES )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      IF ( QueryOpacity_Iso( D(iP) / UnitD ) ) THEN
        opES(iP) = opES(iP) * UnitES
      ELSE
        opES(iP) = Zero
      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Yp_P, E, D, T, Ye, Ym ) &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Yp_P, E, D, T, Ye, Ym ) &
    !$ACC COPYOUT( opES )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( opES )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      opES(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Vector

  SUBROUTINE ComputeNeutrinoOpacities_NNS &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Ye, Ym, &
      iMoment, Phi_NNS_0, Phi_NbNS_0, &
      ApplyWeakMagnetism, ApplyManyBodyCorrection )

    ! --- Brem Opacities (Multiple D,T) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E), Ym(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_NNS_0 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Phi_NbNS_0(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: ApplyWeakMagnetism
    LOGICAL,  INTENT(in)  :: ApplyManyBodyCorrection

    INTEGER  :: iX, iE1, iE2
    REAL(DP), ALLOCATABLE :: Mup(:), Mun(:) !Proton and neutron chemical potentials
    REAL(DP), ALLOCATABLE :: LogT_P(:)
    REAL(DP), ALLOCATABLE :: Mun_P(:), Mup_P(:) !neutron and proton chemical potential
    REAL(DP), ALLOCATABLE :: Phi_n(:,:,:), Phi_p(:,:,:) !Scattering kernels on n and p

    REAL(DP), ALLOCATABLE :: S_tot(:)

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( Mup(iX_B:iX_E), Mun(iX_B:iX_E) )
    ALLOCATE( LogT_P(iX_B:iX_E) )
    ALLOCATE( Mun_P(iX_B:iX_E) )
    ALLOCATE( Mup_P(iX_B:iX_E) )
    ALLOCATE( Phi_n(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E) )
    ALLOCATE( Phi_p(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E) )
    ALLOCATE( S_tot(iX_B:iX_E) )

    S_tot = One

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Mun, Mup, LogT_P, Mun_P, Mup_P, &
    !$OMP             Phi_n, Phi_p, Phi_NNS_0, Phi_NbNS_0 ) &
    !$OMP MAP( to: D, T, Ye, Ym, S_tot )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Mun, Mup, LogT_P, Mun_P, Mup_P, &
    !$ACC         Phi_n, Phi_p, Phi_NNS_0, Phi_NbNS_0 ) &
    !$ACC COPYIN( D, T, Ye, Ym, S_tot )
#endif

    ! --- Compute neutron and proton chemical potentials ---

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mun )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Ye, Ym, Mup )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      Mup_P(iX) = Mup(iX) / UnitMp
      Mun_P(iX) = Mun(iX) / UnitMn

      LogT_P(iX)   = LOG10( T(iX) / UnitT )

    END DO

    IF ( ApplyManyBodyCorrection ) THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
      DO iX = iX_B, iX_E

        CALL ComputeNCManyBodyCorrection &
               ( ( D(iX) / UnitD ) * AvogadroConstantMKS * 1.0d-39, &
                 T(iX) / MeV, Ye(iX), S_tot(iX) )

      END DO

    END IF

    ! --- Interpolate Phi_n  ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, Mun_P, LogTs_T, MuBs_T, &
             OS_NNS(iNeutron_NNS,iMoment), NNS_AT(:,:,:,:,iMoment,iNeutron_NNS), Phi_n )

    ! --- Interpolate Phi_p ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, Mup_P, LogTs_T, MuBs_T, &
             OS_NNS(iProton_NNS,iMoment), NNS_AT(:,:,:,:,iMoment,iProton_NNS), Phi_p )

    IF( ApplyWeakMagnetism ) THEN
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX  = iX_B, iX_E
      DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        Phi_NNS_0 (iE1,iE2,iX) = ( Phi_n(iE1,iE2,iX) * NNS_Xi_Nu_N   (iE2)  &
                               + Phi_p(iE1,iE2,iX) * NNS_Xi_Nu_P   (iE2)  ) * S_tot(iX)

        Phi_NbNS_0(iE1,iE2,iX) = ( Phi_n(iE1,iE2,iX) * NNS_Xi_NuBar_N(iE2) &
                               + Phi_p(iE1,iE2,iX) * NNS_Xi_NuBar_P(iE2) ) * S_tot(iX)

      END DO
      END DO
      END DO

    ELSE
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX  = iX_B, iX_E
      DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        Phi_NNS_0 (iE1,iE2,iX) = ( Phi_n(iE1,iE2,iX) + Phi_p(iE1,iE2,iX) ) * S_tot(iX)

        Phi_NbNS_0(iE1,iE2,iX) = ( Phi_n(iE1,iE2,iX) + Phi_p(iE1,iE2,iX) ) * S_tot(iX)

      END DO
      END DO
      END DO

    ENDIF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Mup, Mun, LogT_P, Mun_P, Mup_P, &
    !$OMP      D, T, Ye, Ym, Phi_n, Phi_p, S_tot ) &
    !$OMP MAP( from: Phi_NNS_0, Phi_NbNS_0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Mup, Mun, LogT_P, Mun_P, Mup_P, &
    !$ACC         D, T, Ye, Ym, Phi_n, Phi_p, S_tot ) &
    !$ACC COPYOUT( Phi_NNS_0, Phi_NbNS_0 )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: Phi_NNS_0, Phi_NbNS_0 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( Phi_NNS_0, Phi_NbNS_0 )
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      Phi_NNS_0 (iE1,iE2,iX) = Zero
      Phi_NbNS_0(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NNS

  SUBROUTINE ComputeNeutrinoOpacityRates_NNS &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      J, Phi_NNS_0, Phi_NbNS_0, Eta, Chi, R, &
      Opacity_Mask, DetailedBalanceScattering_Option, EnergyTransfer_Option )

    ! --- Neutrino-Nucleon Scattering Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D       (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T       (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E       (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2      (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J       (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Phi_NNS_0 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Phi_NbNS_0(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta     (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi     (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: R       (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalanceScattering_Option &
        (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    ! Fully weighted fluid-energy source for this spectral row. Half of
    ! each ordered event is counted; sum rows/species for the total.
    REAL(DP), INTENT(out), OPTIONAL :: &
      EnergyTransfer_Option(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_Out, Phi_In
    REAL(DP) :: SUM1, SUM2, SUMR, SUMQ, Event
    INTEGER  :: iE1, iE2, iS, iX, iC

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_In, Phi_Out, iC ) &
    !$OMP MAP( to: Phi_NNS_0, Phi_NbNS_0, T, E, W2, J, D, &
    !$OMP          Opacity_Mask ) &
    !$OMP MAP( from: Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_In, Phi_Out, iC ) &
    !$ACC COPYIN( Phi_NNS_0, Phi_NbNS_0, T, E, W2, J, D, &
    !$ACC         Opacity_Mask ) &
    !$ACC COPYOUT( Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_In, Phi_Out, iC )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero
      SUMR = Zero
      SUMQ = Zero
      iC   = nChirals - (LeptonNumber(iS) + 1) / nChirals

      IF ( Opacity_Mask(iOp_NNS,iX) ) THEN

        IF ( iC == iNu ) THEN !neutrino-nucleon scattering

          DO iE1 = iE_B, iE_E

            IF ( iE1 <= iE2 ) THEN
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_Out = Phi_NNS_0(iE1,iE2,iX) * UnitNNS
              Phi_In  = Phi_Out * DetBal
            ELSE
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_In  = Phi_NNS_0(iE2,iE1,iX) * UnitNNS
              Phi_Out = Phi_In * DetBal
            END IF

            SUM1 = SUM1 + Phi_In  * W2(iE1) * J(iE1,iS,iX)
            SUM2 = SUM2 + Phi_Out * W2(iE1) * ( One - J(iE1,iS,iX) )
            Event =   Phi_In * J(iE1,iS,iX) &
                             * ( One - J(iE2,iS,iX) ) &
                         - Phi_Out * ( One - J(iE1,iS,iX) ) &
                             * J(iE2,iS,iX)
          SUMR = SUMR + W2(iE1) * Event
          IF( PRESENT(EnergyTransfer_Option) ) &
            SUMQ = SUMQ - 0.5_DP * ( E(iE2) - E(iE1) ) &
              * (W2(iE2)/(PlanckConstant*SpeedOfLight)**3) * W2(iE1) * Event

          END DO

        ELSE !antineutrino-nucleon scattering

          DO iE1 = iE_B, iE_E

            IF ( iE1 <= iE2 ) THEN
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_Out = Phi_NbNS_0(iE1,iE2,iX) * UnitNNS
              Phi_In  = Phi_Out * DetBal
            ELSE
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_In  = Phi_NbNS_0(iE2,iE1,iX) * UnitNNS
              Phi_Out = Phi_In * DetBal
            END IF

            SUM1 = SUM1 + Phi_In  * W2(iE1) * J(iE1,iS,iX)
            SUM2 = SUM2 + Phi_Out * W2(iE1) * ( One - J(iE1,iS,iX) )

            Event =   Phi_In * J(iE1,iS,iX) &
                             * ( One - J(iE2,iS,iX) ) &
                         - Phi_Out * ( One - J(iE1,iS,iX) ) &
                             * J(iE2,iS,iX)
          SUMR = SUMR + W2(iE1) * Event
          IF( PRESENT(EnergyTransfer_Option) ) &
            SUMQ = SUMQ - 0.5_DP * ( E(iE2) - E(iE1) ) &
              * (W2(iE2)/(PlanckConstant*SpeedOfLight)**3) * W2(iE1) * Event

          END DO

        END IF

      END IF

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2
      R  (iE2,iS,iX) = SUMR
      IF( PRESENT(EnergyTransfer_Option) ) &
        EnergyTransfer_Option(iE2,iS,iX) = SUMQ

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_NNS


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NNS &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      H_1, H_2, H_3, &
      Phi_NNS_1, Phi_NbNS_1, &
      L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3, Opacity_Mask, &
      DetailedBalanceScattering_Option )

    ! --- Neutrino-Nucleon Scattering Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D         (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T         (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E         (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2        (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1       (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2       (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3       (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Phi_NNS_1 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Phi_NbNS_1(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_1   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_2   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_3   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalanceScattering_Option &
        (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_In, Phi_1_Out
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iC, iE1, iE2, iS, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iC, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_In, Phi_1_Out ) &
    !$OMP MAP( to: Phi_NNS_1, Phi_NbNS_1, T, E, W2, H_1, H_2, H_3, &
    !$OMP          D, Opacity_Mask ) &
    !$OMP MAP( from: L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iC, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_In, Phi_1_Out ) &
    !$ACC COPYIN( Phi_NNS_1, Phi_NbNS_1, T, E, W2, H_1, H_2, H_3, &
    !$ACC         D, Opacity_Mask ) &
    !$ACC COPYOUT( L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iC, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_In, Phi_1_Out )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      iC = nChirals - ( LeptonNumber(iS) + 1 ) / nChirals

      IF ( Opacity_Mask(iOp_NNS,iX) ) THEN

        DO iE1 = iE_B, iE_E

          IF ( iC == iNu ) THEN
            IF ( iE1 <= iE2 ) THEN
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_1_Out = Phi_NNS_1(iE1,iE2,iX) * UnitNNS
              Phi_1_In  = Phi_1_Out * DetBal
            ELSE
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_1_In  = Phi_NNS_1(iE2,iE1,iX) * UnitNNS
              Phi_1_Out = Phi_1_In * DetBal
            END IF
          ELSE
            IF ( iE1 <= iE2 ) THEN
              DetBal = DetailedBalance_InelasticScattering &
                         ( E(iE1), E(iE2), T(iX) )
              Phi_1_Out = Phi_NbNS_1(iE1,iE2,iX) * UnitNNS
              Phi_1_In  = Phi_1_Out * DetBal
            ELSE
              DetBal = CachedDetailedBalanceScattering &
                ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
              Phi_1_In  = Phi_NbNS_1(iE2,iE1,iX) * UnitNNS
              Phi_1_Out = Phi_1_In * DetBal
            END IF
          END IF

          SUM1 = SUM1 + Phi_1_In  * W2(iE1) * H_1(iE1,iS,iX)
          SUM2 = SUM2 + Phi_1_Out * W2(iE1) * H_1(iE1,iS,iX)
          SUM3 = SUM3 + Phi_1_In  * W2(iE1) * H_2(iE1,iS,iX)
          SUM4 = SUM4 + Phi_1_Out * W2(iE1) * H_2(iE1,iS,iX)
          SUM5 = SUM5 + Phi_1_In  * W2(iE1) * H_3(iE1,iS,iX)
          SUM6 = SUM6 + Phi_1_Out * W2(iE1) * H_3(iE1,iS,iX)

        END DO

      END IF

      L_In_1 (iE2,iS,iX) = SUM1
      L_Out_1(iE2,iS,iX) = SUM2
      L_In_2 (iE2,iS,iX) = SUM3
      L_Out_2(iE2,iS,iX) = SUM4
      L_In_3 (iE2,iS,iX) = SUM5
      L_Out_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NNS


  SUBROUTINE ComputeNeutrinoOpacities_NES &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Ye, Ym, iMoment, H_I, H_II )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: H_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: H_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP), ALLOCATABLE :: LogT_P(:), LogEta_P(:)
    INTEGER  :: iX, iE1, iE2, iH_I, iH_II

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E) )

    iH_I  = ( iMoment - 1 ) * 2 + 1
    iH_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogT_P, LogEta_P, H_I, H_II ) &
    !$OMP MAP( to: D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogT_P, LogEta_P, H_I, H_II ) &
    !$ACC COPYIN( D, T, Ye, Ym )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Ye, Ym, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX)   = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate HI  ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH_I), NES_AT(:,:,:,:,iH_I,1), H_I )

    ! --- Interpolate HII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH_II), NES_AT(:,:,:,:,iH_II,1), H_II )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogT_P, LogEta_P, D, T, Ye, Ym ) &
    !$OMP MAP( from: H_I, H_II )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogT_P, LogEta_P, D, T, Ye, Ym ) &
    !$ACC COPYOUT( H_I, H_II )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: H_I, H_II )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( H_I, H_II )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      H_I (iE1,iE2,iX) = Zero
      H_II(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES


  SUBROUTINE ComputeNeutrinoOpacityRates_NES &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      J, H_I, H_II, Eta, Chi, R, Opacity_Mask, &
      DetailedBalanceScattering_Option, EnergyTransfer_Option )

    ! --- Neutrino-Electron Scattering Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D   (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T   (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E   (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2  (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: R   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalanceScattering_Option &
        (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    ! Fully weighted fluid-energy source for this spectral row. Half of
    ! each ordered event is counted; sum rows/species for the total.
    REAL(DP), INTENT(out), OPTIONAL :: &
      EnergyTransfer_Option(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_Out, Phi_In
    REAL(DP) :: SUM1, SUM2, SUMR, SUMQ, Event
    INTEGER  :: iE1, iE2, iS, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_In, Phi_Out ) &
    !$OMP MAP( to: H_I, H_II, T, E, W2, J, D, Opacity_Mask ) &
    !$OMP MAP( from: Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_In, Phi_Out ) &
    !$ACC COPYIN( H_I, H_II, T, E, W2, J, D, Opacity_Mask ) &
    !$ACC COPYOUT( Eta, Chi, R, EnergyTransfer_Option ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_In, Phi_Out )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero
      SUMR = Zero
      SUMQ = Zero

      IF ( Opacity_Mask(iOp_NES,iX) ) THEN

        DO iE1 = iE_B, iE_E

          IF ( iE1 <= iE2 ) THEN
            DetBal = CachedDetailedBalanceScattering &
              ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
            Phi_Out = (   C1(iS) * H_I (iE1,iE2,iX) &
                        + C2(iS) * H_II(iE1,iE2,iX) ) * UnitNES
            Phi_In  = Phi_Out * DetBal
          ELSE
            DetBal = CachedDetailedBalanceScattering &
              ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
            Phi_In  = (   C1(iS) * H_I (iE2,iE1,iX) &
                        + C2(iS) * H_II(iE2,iE1,iX) ) * UnitNES
            Phi_Out = Phi_In * DetBal
          END IF

          SUM1 = SUM1 + Phi_In  * W2(iE1) * J(iE1,iS,iX)
          SUM2 = SUM2 + Phi_Out * W2(iE1) * ( One - J(iE1,iS,iX) )

          Event =   Phi_In * J(iE1,iS,iX) &
                           * ( One - J(iE2,iS,iX) ) &
                       - Phi_Out * ( One - J(iE1,iS,iX) ) &
                           * J(iE2,iS,iX)
          SUMR = SUMR + W2(iE1) * Event
          IF( PRESENT(EnergyTransfer_Option) ) &
            SUMQ = SUMQ - 0.5_DP * ( E(iE2) - E(iE1) ) &
              * (W2(iE2)/(PlanckConstant*SpeedOfLight)**3) * W2(iE1) * Event

        END DO

      END IF

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2
      R  (iE2,iS,iX) = SUMR
      IF( PRESENT(EnergyTransfer_Option) ) &
        EnergyTransfer_Option(iE2,iS,iX) = SUMQ

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_NES


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NES &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      H_1, H_2, H_3, &
      H_I_1, H_II_1, L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3, &
      Opacity_Mask, DetailedBalanceScattering_Option )

    ! --- Neutrino-Electron Scattering Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E      (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_I_1  (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_II_1 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_1 (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_2 (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_3 (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalanceScattering_Option &
        (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Out, Phi_1_In
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_In, Phi_1_Out ) &
    !$OMP MAP( to: H_I_1, H_II_1, T, E, W2, H_1, H_2, H_3, D, &
    !$OMP          Opacity_Mask ) &
    !$OMP MAP( from: L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_In, Phi_1_Out ) &
    !$ACC COPYIN( H_I_1, H_II_1, T, E, W2, H_1, H_2, H_3, D, &
    !$ACC         Opacity_Mask ) &
    !$ACC COPYOUT( L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_In, Phi_1_Out )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      IF ( Opacity_Mask(iOp_NES,iX) ) THEN

        DO iE1 = iE_B, iE_E

          IF ( iE1 <= iE2 ) THEN
            DetBal = CachedDetailedBalanceScattering &
              ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
            Phi_1_Out = (   C1(iS) * H_I_1 (iE1,iE2,iX) &
                          + C2(iS) * H_II_1(iE1,iE2,iX) ) * UnitNES
            Phi_1_In  = Phi_1_Out * DetBal
          ELSE
            DetBal = CachedDetailedBalanceScattering &
              ( iE1, iE2, iX, E, T, DetailedBalanceScattering_Option )
            Phi_1_In  = (   C1(iS) * H_I_1 (iE2,iE1,iX) &
                          + C2(iS) * H_II_1(iE2,iE1,iX) ) * UnitNES
            Phi_1_Out = Phi_1_In * DetBal
          END IF

          SUM1 = SUM1 + Phi_1_In  * W2(iE1) * H_1(iE1,iS,iX)
          SUM2 = SUM2 + Phi_1_Out * W2(iE1) * H_1(iE1,iS,iX)
          SUM3 = SUM3 + Phi_1_In  * W2(iE1) * H_2(iE1,iS,iX)
          SUM4 = SUM4 + Phi_1_Out * W2(iE1) * H_2(iE1,iS,iX)
          SUM5 = SUM5 + Phi_1_In  * W2(iE1) * H_3(iE1,iS,iX)
          SUM6 = SUM6 + Phi_1_Out * W2(iE1) * H_3(iE1,iS,iX)

        END DO

      END IF

      L_In_1 (iE2,iS,iX) = SUM1
      L_Out_1(iE2,iS,iX) = SUM2
      L_In_2 (iE2,iS,iX) = SUM3
      L_Out_2(iE2,iS,iX) = SUM4
      L_In_3 (iE2,iS,iX) = SUM5
      L_Out_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NES


  SUBROUTINE ComputeNeutrinoOpacities_Pair &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Ye, Ym, iMoment, J_I, J_II )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: J_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: J_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP), ALLOCATABLE :: LogT_P(:), LogEta_P(:)
    INTEGER  :: iX, iE1, iE2, iJ_I, iJ_II

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E) )

    iJ_I  = ( iMoment - 1 ) * 2 + 1
    iJ_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogT_P, LogEta_P, J_I, J_II ) &
    !$OMP MAP( to: D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogT_P, LogEta_P, J_I, J_II ) &
    !$ACC COPYIN( D, T, Ye, Ym )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Ye, Ym, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      LogT_P  (iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) &
                            / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate JI  ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_I), Pair_AT(:,:,:,:,iJ_I,1), J_I )

    ! --- Interpolate JII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_II), Pair_AT(:,:,:,:,iJ_II,1), J_II )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogT_P, LogEta_P, D, T, Ye, Ym ) &
    !$OMP MAP( from: J_I, J_II )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogT_P, LogEta_P, D, T, Ye, Ym ) &
    !$ACC COPYOUT( J_I, J_II )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: J_I, J_II)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( J_I, J_II)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      J_I (iE1,iE2,iX) = Zero
      J_II(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair




  SUBROUTINE ComputeNeutrinoOpacityRates_Pair &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      J, J_I, J_II, Eta, Chi, R, Opacity_Mask, &
      DetailedBalancePair_Option, EnergyTransfer_Option )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D   (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T   (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E   (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2  (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: R   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalancePair_Option(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    ! Fully weighted fluid-energy source for this spectral row. Half of
    ! each ordered event is counted; sum rows/species for the total.
    REAL(DP), INTENT(out), OPTIONAL :: &
      EnergyTransfer_Option(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2, SUMR, SUMQ, Event
    INTEGER  :: iX, iE1, iE2, iS, iS_A

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Pro, &
    !$OMP          Phi_0_Ann ) &
    !$OMP MAP( to: J_I, J_II, T, E, W2, J, D, Opacity_Mask ) &
    !$OMP MAP( from: Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Pro, &
    !$ACC          Phi_0_Ann ) &
    !$ACC COPYIN( J_I, J_II, T, E, W2, J, D, Opacity_Mask ) &
    !$ACC COPYOUT( Eta, Chi, R, EnergyTransfer_Option ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Pro, &
    !$OMP          Phi_0_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUMR = Zero
      SUMQ = Zero

      IF ( Opacity_Mask(iOp_Pair,iX) ) THEN

        DO iE1 = iE_B, iE_E

          DetBal = CachedDetailedBalancePair &
            ( iE1, iE2, iX, E, T, DetailedBalancePair_Option )

          IF ( iE1 <= iE2 ) THEN
            Phi_0_Ann = (   C1(iS) * J_I (iE1,iE2,iX) &
                          + C2(iS) * J_II(iE1,iE2,iX) ) * UnitPair
          ELSE
            Phi_0_Ann = (   C1(iS) * J_II(iE2,iE1,iX) &
                          + C2(iS) * J_I (iE2,iE1,iX) ) * UnitPair
          END IF
          Phi_0_Pro = Phi_0_Ann * DetBal

          SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iS_A,iX) )
          SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iS_A,iX)

          Event =   Phi_0_Pro * ( One - J(iE1,iS_A,iX) ) &
                           * ( One - J(iE2,iS,iX) ) &
                       - Phi_0_Ann * J(iE1,iS_A,iX) &
                           * J(iE2,iS,iX)
          SUMR = SUMR + W2(iE1) * Event
          IF( PRESENT(EnergyTransfer_Option) ) &
            SUMQ = SUMQ - 0.5_DP * ( E(iE2) + E(iE1) ) &
              * (W2(iE2)/(PlanckConstant*SpeedOfLight)**3) * W2(iE1) * Event

        END DO

      END IF

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2
      R  (iE2,iS,iX) = SUMR
      IF( PRESENT(EnergyTransfer_Option) ) &
        EnergyTransfer_Option(iE2,iS,iX) = SUMQ

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_Pair


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      H_1, H_2, H_3, &
      J_I_1, J_II_1, L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3, &
      Opacity_Mask, DetailedBalancePair_Option )

    ! --- e^+ e^- Pair Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E      (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_I_1  (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_II_1 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalancePair_Option(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Ann, Phi_1_Pro
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iS_A, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann ) &
    !$OMP MAP( to: J_I_1, J_II_1, T, E, W2, H_1, H_2, H_3, D, &
    !$OMP          Opacity_Mask ) &
    !$OMP MAP( from: L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_Pro, Phi_1_Ann ) &
    !$ACC COPYIN( J_I_1, J_II_1, T, E, W2, H_1, H_2, H_3, D, &
    !$ACC         Opacity_Mask ) &
    !$ACC COPYOUT( L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      IF ( Opacity_Mask(iOp_Pair,iX) ) THEN

        DO iE1 = iE_B, iE_E

          DetBal = CachedDetailedBalancePair &
            ( iE1, iE2, iX, E, T, DetailedBalancePair_Option )

          IF ( iE1 <= iE2 ) THEN
            Phi_1_Ann = (   C1(iS) * J_I_1 (iE1,iE2,iX) &
                          + C2(iS) * J_II_1(iE1,iE2,iX) ) * UnitPair
          ELSE
            Phi_1_Ann = (   C1(iS) * J_II_1(iE2,iE1,iX) &
                          + C2(iS) * J_I_1 (iE2,iE1,iX) ) * UnitPair
          END IF
          Phi_1_Pro = Phi_1_Ann * DetBal

          SUM1 = SUM1 + Phi_1_Pro * W2(iE1) * H_1(iE1,iS_A,iX)
          SUM2 = SUM2 + Phi_1_Ann * W2(iE1) * H_1(iE1,iS_A,iX)
          SUM3 = SUM3 + Phi_1_Pro * W2(iE1) * H_2(iE1,iS_A,iX)
          SUM4 = SUM4 + Phi_1_Ann * W2(iE1) * H_2(iE1,iS_A,iX)
          SUM5 = SUM5 + Phi_1_Pro * W2(iE1) * H_3(iE1,iS_A,iX)
          SUM6 = SUM6 + Phi_1_Ann * W2(iE1) * H_3(iE1,iS_A,iX)

        END DO

      END IF

      L_Pro_1(iE2,iS,iX) = SUM1
      L_Ann_1(iE2,iS,iX) = SUM2
      L_Pro_2(iE2,iS,iX) = SUM3
      L_Ann_2(iE2,iS,iX) = SUM4
      L_Pro_3(iE2,iS,iX) = SUM5
      L_Ann_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Pair


  SUBROUTINE ComputeNeutrinoOpacities_NuPair &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Ye, Ym, iMoment, J_I, J_II )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: J_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: J_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP), ALLOCATABLE :: LogT_P(:), LogEta_P(:), SignEta_P(:)
    INTEGER  :: iX, iE1, iE2, iJ_I, iJ_II
    REAL(DP) :: J_I_tmp, J_II_tmp

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E), SignEta_P(iX_B:iX_E) )

    iJ_I  = ( iMoment - 1 ) * 2 + 1
    iJ_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogT_P, LogEta_P, SignEta_P, J_I, J_II ) &
    !$OMP MAP( to: D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogT_P, LogEta_P, SignEta_P, J_I, J_II ) &
    !$ACC COPYIN( D, T, Ye, Ym )
#endif

    ! --- Compute Electron Neutrino Chemical Potential ---

    CALL ComputeElectronNeutrinoChemicalPotential_TABLE &
           ( D, T, Ye, Ym, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      LogT_P   (iX) = LOG10( T(iX) / UnitT )
      !The J_I and J_II have the following symmetry w.r.t. the sign of Eta:
      !J_I(-Eta) = J_II(Eta), so we store a sign to select the correct J_I/J_II 
      !depending on the sign of Eta
      SignEta_P(iX) = MAX(SIGN(One,LogEta_P(iX)),Zero)
      LogEta_P (iX) = LOG10( ABS(LogEta_P(iX)) &
                            / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate JI  ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_I), Pair_AT(:,:,:,:,iJ_I,1), J_I )

    ! --- Interpolate JII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_II), Pair_AT(:,:,:,:,iJ_II,1), J_II )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE ( J_I_tmp, J_II_tmp )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE ( J_I_tmp, J_II_tmp )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE ( J_I_tmp, J_II_tmp )
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E

      J_I_tmp  = J_I (iE1,iE2,iX)
      J_II_tmp = J_II(iE1,iE2,iX)

      J_I (iE1,iE2,iX) = SignEta_P(iX)*J_I_tmp  + (One-SignEta_P(iX))*J_II_tmp
      J_II(iE1,iE2,iX) = SignEta_P(iX)*J_II_tmp + (One-SignEta_P(iX))*J_I_tmp
    
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogT_P, LogEta_P, SignEta_P, D, T, Ye, Ym ) &
    !$OMP MAP( from: J_I, J_II )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogT_P, LogEta_P, SignEta_P, D, T, Ye, Ym ) &
    !$ACC COPYOUT( J_I, J_II )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: J_I, J_II)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( J_I, J_II)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      J_I (iE1,iE2,iX) = Zero
      J_II(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NuPair


  SUBROUTINE ComputeNeutrinoOpacityRates_NuPair &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      J, J_I, J_II, Eta, Chi, R, Opacity_Mask, &
      DetailedBalancePair_Option, EnergyTransfer_Option )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D   (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T   (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E   (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2  (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: R   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalancePair_Option(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    ! Fully weighted fluid-energy source for this spectral row. Half of
    ! each ordered event is counted; sum rows/species for the total.
    REAL(DP), INTENT(out), OPTIONAL :: &
      EnergyTransfer_Option(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2, SUMR, SUMQ, Event
    INTEGER  :: iX, iE1, iE2, iS, iS_A

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Pro, &
    !$OMP          Phi_0_Ann ) &
    !$OMP MAP( to: J_I, J_II, T, E, W2, J, D, Opacity_Mask ) &
    !$OMP MAP( from: Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Pro, &
    !$ACC          Phi_0_Ann ) &
    !$ACC COPYIN( J_I, J_II, T, E, W2, J, D, Opacity_Mask ) &
    !$ACC COPYOUT( Eta, Chi, R, EnergyTransfer_Option ) &
    !$ACC PRESENT( C1_NuPair, C2_NuPair )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Pro, &
    !$OMP          Phi_0_Ann )
#endif
      DO iX  = iX_B, iX_E
      DO iS  = iS_B, iS_E
      DO iE2 = iE_B, iE_E

        ! Get index for corresponding anti-neutrino
        iS_A = iS + 2*MOD(iS,2) - 1

        SUM1 = Zero
        SUM2 = Zero
        SUMR = Zero
        SUMQ = Zero

        IF ( Opacity_Mask(iOp_NuPair,iX) ) THEN

          DO iE1 = iE_B, iE_E

            DetBal = CachedDetailedBalancePair &
              ( iE1, iE2, iX, E, T, DetailedBalancePair_Option )

            IF ( iE1 <= iE2 ) THEN
              Phi_0_Ann = (   C1_NuPair(iS) * J_I (iE1,iE2,iX) &
                            + C2_NuPair(iS) * J_II(iE1,iE2,iX) ) * UnitPair
            ELSE
              Phi_0_Ann = (   C1_NuPair(iS) * J_II(iE2,iE1,iX) &
                            + C2_NuPair(iS) * J_I (iE2,iE1,iX) ) * UnitPair
            END IF
            Phi_0_Pro = Phi_0_Ann * DetBal

            SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iS_A,iX) )
            SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iS_A,iX)

            Event =   Phi_0_Pro * ( One - J(iE1,iS_A,iX) ) &
                             * ( One - J(iE2,iS,iX) ) &
                         - Phi_0_Ann * J(iE1,iS_A,iX) &
                             * J(iE2,iS,iX)
          SUMR = SUMR + W2(iE1) * Event
          IF( PRESENT(EnergyTransfer_Option) ) &
            SUMQ = SUMQ - 0.5_DP * ( E(iE2) + E(iE1) ) &
              * (W2(iE2)/(PlanckConstant*SpeedOfLight)**3) * W2(iE1) * Event

          END DO

        END IF

        Eta(iE2,iS,iX) = SUM1
        Chi(iE2,iS,iX) = SUM1 + SUM2
        R  (iE2,iS,iX) = SUMR
      IF( PRESENT(EnergyTransfer_Option) ) &
        EnergyTransfer_Option(iE2,iS,iX) = SUMQ

      END DO
      END DO
      END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_NuPair


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NuPair &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      H_1, H_2, H_3, &
      J_I_1, J_II_1, L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3, &
      Opacity_Mask, DetailedBalancePair_Option )

    ! --- nue+nuebar Pair Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E      (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_I_1  (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_II_1 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalancePair_Option(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Ann, Phi_1_Pro
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iS_A, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann ) &
    !$OMP MAP( to: J_I_1, J_II_1, T, E, W2, H_1, H_2, H_3, D, &
    !$OMP          Opacity_Mask ) &
    !$OMP MAP( from: L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_Pro, Phi_1_Ann ) &
    !$ACC COPYIN( J_I_1, J_II_1, T, E, W2, H_1, H_2, H_3, D, &
    !$ACC         Opacity_Mask ) &
    !$ACC COPYOUT( L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 ) &
    !$ACC PRESENT( C1_NuPair, C2_NuPair )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      IF ( Opacity_Mask(iOp_NuPair,iX) ) THEN

        DO iE1 = iE_B, iE_E

          DetBal = CachedDetailedBalancePair &
            ( iE1, iE2, iX, E, T, DetailedBalancePair_Option )

          IF ( iE1 <= iE2 ) THEN
            Phi_1_Ann = (   C1_NuPair(iS) * J_I_1 (iE1,iE2,iX) &
                          + C2_NuPair(iS) * J_II_1(iE1,iE2,iX) ) * UnitPair
          ELSE
            Phi_1_Ann = (   C1_NuPair(iS) * J_II_1(iE2,iE1,iX) &
                          + C2_NuPair(iS) * J_I_1 (iE2,iE1,iX) ) * UnitPair
          END IF
          Phi_1_Pro = Phi_1_Ann * DetBal

          SUM1 = SUM1 + Phi_1_Pro * W2(iE1) * H_1(iE1,iS_A,iX)
          SUM2 = SUM2 + Phi_1_Ann * W2(iE1) * H_1(iE1,iS_A,iX)
          SUM3 = SUM3 + Phi_1_Pro * W2(iE1) * H_2(iE1,iS_A,iX)
          SUM4 = SUM4 + Phi_1_Ann * W2(iE1) * H_2(iE1,iS_A,iX)
          SUM5 = SUM5 + Phi_1_Pro * W2(iE1) * H_3(iE1,iS_A,iX)
          SUM6 = SUM6 + Phi_1_Ann * W2(iE1) * H_3(iE1,iS_A,iX)

        END DO

      END IF

      L_Pro_1(iE2,iS,iX) = SUM1
      L_Ann_1(iE2,iS,iX) = SUM2
      L_Pro_2(iE2,iS,iX) = SUM3
      L_Ann_2(iE2,iS,iX) = SUM4
      L_Pro_3(iE2,iS,iX) = SUM5
      L_Ann_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NuPair


  SUBROUTINE ComputeNeutrinoOpacities_Brem &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Ye, Ym, S_sigma )

    ! --- Brem Opacities (Multiple D,T) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ye(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Ym(iX_B:iX_E)
    REAL(DP), INTENT(out) :: S_Sigma(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    INTEGER  :: iX, iE1, iE2
    REAL(DP), ALLOCATABLE :: Xp(:), Xn(:) !Proton and neutron mass fractions
    REAL(DP), ALLOCATABLE :: LogT_P(:)
    REAL(DP), ALLOCATABLE :: LogDX_P(:,:)

#ifdef MICROPHYSICS_WEAKLIB

    ALLOCATE( Xp(iX_B:iX_E), Xn(iX_B:iX_E) )
    ALLOCATE( LogT_P(iX_B:iX_E) )
    ALLOCATE( LogDX_P(3,iX_B:iX_E) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Xp, Xn, LogT_P, LogDX_P, S_sigma ) &
    !$OMP MAP( to: D, T, Ye, Ym )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Xp, Xn, LogT_P, LogDX_P, S_sigma ) &
    !$ACC COPYIN( D, T, Ye, Ym )
#endif

    ! --- Compute proton and neutron fractions ---

    CALL ComputeProtonMassFraction_TABLE &
           ( D, T, Ye, Ym, Xp )

    CALL ComputeNeutronMassFraction_TABLE &
           ( D, T, Ye, Ym, Xn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      LogDX_P(1,iX) = LOG10( D(iX) * Xp(iX) / UnitD )
      LogDX_P(2,iX) = LOG10( D(iX) * Xn(iX) / UnitD )
      LogDX_P(3,iX) = LOG10( D(iX) * SQRT(ABS(Xp(iX)*Xn(iX))) / UnitD )

      LogT_P(iX)    = LOG10( T(iX) / UnitT )    

    END DO

    CALL SumLogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogDX_P, LogT_P, LogDs_T, LogTs_T, Alpha_Brem, &
             OS_Brem(1,1), Brem_AT(:,:,:,:,1,1), S_sigma )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Xp, Xn, LogT_P, LogDX_P, D, T, Ye, Ym ) &
    !$OMP MAP( from: S_sigma )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Xp, Xn, LogT_P, LogDX_P, D, T, Ye, Ym ) &
    !$ACC COPYOUT( S_sigma )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: S_sigma )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( S_sigma )
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      S_sigma(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Brem

  
  SUBROUTINE ComputeNeutrinoOpacityRates_Brem &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      J, S_Sigma, Eta, Chi, R, Opacity_Mask, &
      DetailedBalancePair_Option, EnergyTransfer_Option )

    ! --- Brem Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E      (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J      (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: S_Sigma(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: R      (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalancePair_Option(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    ! Fully weighted fluid-energy source for this spectral row. Half of
    ! each ordered event is counted; sum rows/species for the total.
    REAL(DP), INTENT(out), OPTIONAL :: &
      EnergyTransfer_Option(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2, SUMR, SUMQ, Event
    INTEGER  :: iX, iE1, iE2, iS, iS_A

    !For medium modification from Fischer A&A 593, A103 (2016)
    !nuclear matter density, 0.16 fm^-3
    REAL(DP), PARAMETER :: D_0 = 2.656862401594383d14
    REAL(DP)            :: medium_fac

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Ann, &
    !$OMP          Phi_0_Pro, medium_fac ) &
    !$OMP MAP( to: S_Sigma, T, E, W2, J, D, Opacity_Mask ) &
    !$OMP MAP( from: Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Ann, &
    !$ACC          Phi_0_Pro, medium_fac ) &
    !$ACC COPYIN( S_Sigma, T, E, W2, J, D, Opacity_Mask ) &
    !$ACC COPYOUT( Eta, Chi, R, EnergyTransfer_Option )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUMR, SUMQ, Event, DetBal, Phi_0_Ann, &
    !$OMP          Phi_0_Pro, medium_fac )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUMR = Zero
      SUMQ = Zero

      !Equation 11 from Fischer Fischer A&A 593, A103 (2016)
      !medium_fac = (1.0d0 / ( 1.0d0 + 1.0d0/3.0d0 * (D(iX) / UnitD / D_0)**(1.0d0/3.0d0) ) )**6
      !Set it to one to switch off
      medium_fac = 1.0d0

      IF ( Opacity_Mask(iOp_Brem,iX) ) THEN

        DO iE1 = iE_B, iE_E

          DetBal = CachedDetailedBalancePair &
            ( iE1, iE2, iX, E, T, DetailedBalancePair_Option )

          IF ( iE1 <= iE2 ) THEN
            Phi_0_Ann = S_Sigma(iE1,iE2,iX) * 3.0d0 * Brem_const * medium_fac
          ELSE
            Phi_0_Ann = S_Sigma(iE2,iE1,iX) * 3.0d0 * Brem_const * medium_fac
          END IF
          Phi_0_Pro = Phi_0_Ann * DetBal

          SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iS_A,iX) )
          SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iS_A,iX)

          Event =   Phi_0_Pro * ( One - J(iE1,iS_A,iX) ) &
                           * ( One - J(iE2,iS,iX) ) &
                       - Phi_0_Ann * J(iE1,iS_A,iX) &
                           * J(iE2,iS,iX)
          SUMR = SUMR + W2(iE1) * Event
          IF( PRESENT(EnergyTransfer_Option) ) &
            SUMQ = SUMQ - 0.5_DP * ( E(iE2) + E(iE1) ) &
              * (W2(iE2)/(PlanckConstant*SpeedOfLight)**3) * W2(iE1) * Event

        END DO

      END IF

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2
      R  (iE2,iS,iX) = SUMR
      IF( PRESENT(EnergyTransfer_Option) ) &
        EnergyTransfer_Option(iE2,iS,iX) = SUMQ

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_Brem


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Brem &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, D, T, E, W2, &
      H_1, H_2, H_3, &
      S_Sigma, L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3, &
      Opacity_Mask, DetailedBalancePair_Option )

    ! --- N N Brem Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T      (iX_B:iX_E)
    REAL(DP), INTENT(in)  :: E      (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: S_Sigma(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    LOGICAL,  INTENT(in)  :: Opacity_Mask(nOp,iX_B:iX_E)
    REAL(DP), INTENT(in), OPTIONAL :: &
      DetailedBalancePair_Option(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Ann, Phi_1_Pro
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iS_A, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann ) &
    !$OMP MAP( to: S_Sigma, T, E, W2, H_1, H_2, H_3, D, &
    !$OMP          Opacity_Mask ) &
    !$OMP MAP( from: L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_Pro, Phi_1_Ann ) &
    !$ACC COPYIN( S_Sigma, T, E, W2, H_1, H_2, H_3, D, &
    !$ACC         Opacity_Mask ) &
    !$ACC COPYOUT( L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      IF ( Opacity_Mask(iOp_Brem,iX) ) THEN

        DO iE1 = iE_B, iE_E

          DetBal = CachedDetailedBalancePair &
            ( iE1, iE2, iX, E, T, DetailedBalancePair_Option )

          IF ( iE1 <= iE2 ) THEN
            Phi_1_Ann = - S_Sigma(iE1,iE2,iX) * Brem_const
          ELSE
            Phi_1_Ann = - S_Sigma(iE2,iE1,iX) * Brem_const
          END IF
          Phi_1_Pro = Phi_1_Ann * DetBal

          SUM1 = SUM1 + Phi_1_Pro * W2(iE1) * H_1(iE1,iS_A,iX)
          SUM2 = SUM2 + Phi_1_Ann * W2(iE1) * H_1(iE1,iS_A,iX)
          SUM3 = SUM3 + Phi_1_Pro * W2(iE1) * H_2(iE1,iS_A,iX)
          SUM4 = SUM4 + Phi_1_Ann * W2(iE1) * H_2(iE1,iS_A,iX)
          SUM5 = SUM5 + Phi_1_Pro * W2(iE1) * H_3(iE1,iS_A,iX)
          SUM6 = SUM6 + Phi_1_Ann * W2(iE1) * H_3(iE1,iS_A,iX)

        END DO

      END IF

      L_Pro_1(iE2,iS,iX) = SUM1
      L_Ann_1(iE2,iS,iX) = SUM2
      L_Pro_2(iE2,iS,iX) = SUM3
      L_Ann_2(iE2,iS,iX) = SUM4
      L_Pro_3(iE2,iS,iX) = SUM5
      L_Ann_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Brem

  
  FUNCTION FermiDirac_Scalar( E, Mu, kT ) RESULT( FermiDirac )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT
    REAL(DP) :: FermiDirac

    REAL(DP) :: X, Z

    X = ( E - Mu ) / kT
    Z = EXP( - ABS(X) )
    FermiDirac = MERGE( Z, One, X >= Zero ) / ( One + Z )

    RETURN
  END FUNCTION FermiDirac_Scalar


  FUNCTION FermiDirac_Vector( E, Mu, kT ) RESULT( FermiDirac )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT
    REAL(DP) :: FermiDirac(SIZE(E))

    REAL(DP) :: X(SIZE(E)), Z(SIZE(E))

    X = ( E - Mu ) / kT
    Z = EXP( - ABS(X) )
    FermiDirac = MERGE( Z, One, X >= Zero ) / ( One + Z )

    RETURN
  END FUNCTION FermiDirac_Vector


  FUNCTION dFermiDiracdT_Scalar( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudT, T
    REAL(DP) :: dFermiDiracdT

    REAL(DP) :: Z, Susceptibility

    Z = EXP( - ABS( ( E - Mu ) / kT ) )
    Susceptibility = ( Z / ( One + Z ) ) / ( One + Z )

    dFermiDiracdT &
      = Susceptibility * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT_Scalar


  FUNCTION dFermiDiracdT_Vector( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudT, T
    REAL(DP) :: dFermiDiracdT(SIZE(E))

    REAL(DP) :: Z(SIZE(E)), Susceptibility(SIZE(E))

    Z = EXP( - ABS( ( E - Mu ) / kT ) )
    Susceptibility = ( Z / ( One + Z ) ) / ( One + Z )

    dFermiDiracdT &
      = Susceptibility * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT_Vector


  FUNCTION dFermiDiracdY_Scalar( E, Mu, kT, dMudY, T ) RESULT( dFermiDiracdY )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY

    REAL(DP) :: Z, Susceptibility

    Z = EXP( - ABS( ( E - Mu ) / kT ) )
    Susceptibility = ( Z / ( One + Z ) ) / ( One + Z )

    dFermiDiracdY &
      = Susceptibility * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Scalar


  FUNCTION dFermiDiracdY_Vector( E, Mu, kT, dMudY, T ) RESULT( dFermiDiracdY )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY(SIZE(E))

    REAL(DP) :: Z(SIZE(E)), Susceptibility(SIZE(E))

    Z = EXP( - ABS( ( E - Mu ) / kT ) )
    Susceptibility = ( Z / ( One + Z ) ) / ( One + Z )

    dFermiDiracdY &
      = Susceptibility * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Vector


  SUBROUTINE FermiDiracPair( E, Mu, kT, Occupation, Vacancy, Susceptibility )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: E, Mu, kT
    REAL(DP), INTENT(out) :: Occupation, Vacancy, Susceptibility
    REAL(DP) :: X, Z, Inverse

    ! No positive exponential and no artificial occupation floor. MERGE has
    ! safe arguments on both sides, also under vectorized/device evaluation.
    X = ( E - Mu ) / kT
    Z = EXP( - ABS(X) )
    Inverse = One / ( One + Z )
    Occupation     = MERGE( Z, One, X >= Zero ) * Inverse
    Vacancy        = MERGE( One, Z, X >= Zero ) * Inverse
    Susceptibility = ( Z * Inverse ) * Inverse

  END SUBROUTINE FermiDiracPair


  FUNCTION DivideByFermiDirac( Value, X ) RESULT( Quotient )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: Value, X
    REAL(DP) :: Quotient, Z

    ! Value / f(X). Preserve a representable quotient even if f underflows.
    ! Value is a nonnegative emissivity/absorption coefficient. If that input
    ! has already underflowed, its missing information cannot be recovered.
    IF( Value == Zero )THEN
      Quotient = Zero
    ELSE IF( X <= Zero )THEN
      Quotient = Value * ( One + EXP(X) )
    ELSE IF( X < - LOG(TINY(One)) )THEN
      Z = EXP( - X )
      Quotient = ( Value / Z ) * ( One + Z )
    ELSE
      Quotient = Value + EXP( LOG(Value) + X )
    END IF

  END FUNCTION DivideByFermiDirac

  SUBROUTINE SetChargedCurrentMasks &
    ( nX, nDOF_X, D, T, Ye, Ym, &
      ElectronChargedCurrentMask, MuonChargedCurrentMask )

    INTEGER,  INTENT(in)  :: nX, nDOF_X
    REAL(DP), INTENT(in)  :: D(nX), T(nX), Ye(nX), Ym(nX)
    LOGICAL,  INTENT(out) :: ElectronChargedCurrentMask(nX)
    LOGICAL,  INTENT(out) :: MuonChargedCurrentMask(nX)

    INTEGER :: iN_X, iNodeX, iX
    LOGICAL :: ElectronElementActive, MuonElementActive
    REAL(DP), ALLOCATABLE :: Ah(:)

    ALLOCATE( Ah(nX) )
    Ah = Zero

#ifdef MICROPHYSICS_WEAKLIB
    IF( use_EC_table > 0 ) &
      CALL ComputeHeavyMassNumber_TABLE( D, T, Ye, Ym, Ah )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET DATA &
    !$OMP MAP( to: D, T, Ye, Ym, Ah ) &
    !$OMP MAP( from: ElectronChargedCurrentMask, MuonChargedCurrentMask )
#elif defined(THORNADO_OACC)
    !$ACC DATA COPYIN( D, T, Ye, Ym, Ah ) &
    !$ACC COPYOUT( ElectronChargedCurrentMask, MuonChargedCurrentMask )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX
      ElectronChargedCurrentMask(iN_X) &
        = QueryOpacity_EmAb_Nucleon( D(iN_X) / UnitD ) .OR. &
          QueryOpacity_ECTable &
            ( D(iN_X) / UnitD, T(iN_X) / UnitT, &
              Ye(iN_X) / UnitYe, Ah(iN_X) )
      MuonChargedCurrentMask(iN_X) &
        = QueryOpacity_EmAb_Muon( D(iN_X) / UnitD )
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( iNodeX, iN_X, ElectronElementActive, MuonElementActive )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( iNodeX, iN_X, ElectronElementActive, MuonElementActive )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iNodeX, iN_X, ElectronElementActive, MuonElementActive )
#endif
    DO iX = 1, nX / nDOF_X

      ElectronElementActive = .TRUE.
      MuonElementActive     = .TRUE.
      DO iNodeX = 1, nDOF_X
        iN_X = iNodeX + ( iX - 1 ) * nDOF_X
        ElectronElementActive &
          = ElectronElementActive .AND. ElectronChargedCurrentMask(iN_X)
        MuonElementActive &
          = MuonElementActive .AND. MuonChargedCurrentMask(iN_X)
      END DO

      DO iNodeX = 1, nDOF_X
        iN_X = iNodeX + ( iX - 1 ) * nDOF_X
        ElectronChargedCurrentMask(iN_X) = ElectronElementActive
        MuonChargedCurrentMask(iN_X)     = MuonElementActive
      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
    !$ACC END DATA
#endif

    DEALLOCATE( Ah )

  END SUBROUTINE SetChargedCurrentMasks


  SUBROUTINE SetOpacityMask &
    ( nX, nDOF_X, D, T, Ye, Ym, Opacity_Mask )

    INTEGER,  INTENT(in)  :: nX, nDOF_X
    REAL(DP), INTENT(in)  :: D(nX), T(nX), Ye(nX), Ym(nX)
    LOGICAL,  INTENT(out) :: Opacity_Mask(nOp,nX)

    INTEGER :: iN_X, iNodeX, iX, iOp
    LOGICAL :: OP_ACTIVE

    REAL(DP), ALLOCATABLE :: Ah(:)

    ALLOCATE( Ah(nX) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET DATA               &
    !$OMP MAP( to: D, T, Ye, Ym )   &
    !$OMP MAP( from: Opacity_Mask ) &
    !$OMP MAP( alloc: Ah )
#elif defined(THORNADO_OACC)
    !$ACC DATA                      &
    !$ACC COPYIN( D, T, Ye, Ym )    &
    !$ACC COPYOUT( Opacity_Mask )   &
    !$ACC CREATE( Ah )
#endif

    CALL ComputeHeavyMassNumber_TABLE &
           ( D, T, Ye, Ym, Ah )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX

      Opacity_Mask(iOp_EmAb_Electron,iN_X) &
        = QueryOpacity_EmAb_Nucleon( D(iN_X) / UnitD )

      Opacity_Mask(iOp_EmAb_Muon,iN_X) &
        = QueryOpacity_EmAb_Muon( D(iN_X) / UnitD )

      Opacity_Mask(iOp_EmAb,iN_X) &
        = Opacity_Mask(iOp_EmAb_Electron,iN_X) .OR. &
          Opacity_Mask(iOp_EmAb_Muon,    iN_X)

      Opacity_Mask(iOp_ECTable,iN_X) &
        = QueryOpacity_ECTable &
            ( D(iN_X) / UnitD, &
              T(iN_X) / Kelvin, &
              Ye(iN_X) / UnitYe, &
              Ah(iN_X) )

      Opacity_Mask(iOp_Iso,iN_X) &
        = QueryOpacity_Iso( D(iN_X) / UnitD )

      Opacity_Mask(iOp_NNS,iN_X) &
        = QueryOpacity_NNS( D(iN_X) / UnitD )

      Opacity_Mask(iOp_NES,iN_X) &
        = QueryOpacity_NES( D(iN_X) / UnitD )

      Opacity_Mask(iOp_Pair,iN_X) &
        = QueryOpacity_Pair( D(iN_X) / UnitD )

      Opacity_Mask(iOp_NuPair,iN_X) &
        = QueryOpacity_NuPair( D(iN_X) / UnitD )

      Opacity_Mask(iOp_Brem,iN_X) &
        = QueryOpacity_Brem( D(iN_X) / UnitD )

      Opacity_Mask(iOp_NMS,iN_X) &
        = QueryOpacity_NMS( D(iN_X) / UnitD )
      Opacity_Mask(iOp_LFE,iN_X) &
        = QueryOpacity_LFE( D(iN_X) / UnitD )
      Opacity_Mask(iOp_LFC,iN_X) &
        = QueryOpacity_LFC( D(iN_X) / UnitD )
      Opacity_Mask(iOp_IMD,iN_X) &
        = QueryOpacity_IMD( D(iN_X) / UnitD )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE(iNodeX,iN_X,OP_ACTIVE)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE(iNodeX,iN_X,OP_ACTIVE)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE(iNodeX,iN_X,OP_ACTIVE)
#endif
    DO iX  = 1, nX / nDOF_X
    DO iOp = 1, nOp

      OP_ACTIVE = .TRUE.

      DO iNodeX = 1, nDOF_X

        iN_X = iNodeX + ( iX - 1 ) * nDOF_X

        OP_ACTIVE = OP_ACTIVE .AND. Opacity_Mask(iOp,iN_X)

      END DO

      DO iNodeX = 1, nDOF_X

        iN_X = iNodeX + ( iX - 1 ) * nDOF_X

        Opacity_Mask(iOp,iN_X) = OP_ACTIVE

      END DO

    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD PRIVATE(iNodeX,iN_X)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(iNodeX,iN_X)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO PRIVATE(iNodeX,iN_X)
#endif
    DO iX = 1, nX / nDOF_X
      DO iNodeX = 1, nDOF_X
        iN_X = iNodeX + ( iX - 1 ) * nDOF_X
        Opacity_Mask(iOp_EmAb,iN_X) = &
          Opacity_Mask(iOp_EmAb_Electron,iN_X) .OR. &
          Opacity_Mask(iOp_EmAb_Muon,    iN_X)
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(THORNADO_OACC)
    !$ACC END DATA
#endif

    DEALLOCATE( Ah )

  END SUBROUTINE SetOpacityMask


END MODULE NeutrinoOpacitiesComputationModule
