MODULE NeutrinoLeptonInteractionsModule

  USE KindModule, ONLY: DP, Zero, One, FourPi
  USE UnitsModule, ONLY: BoltzmannConstant, Centimeter, MeV, PlanckConstant, SpeedOfLight
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeMuonChemicalPotential_TABLE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, iNuM, iNuM_Bar, iNuT, iNuT_Bar
  USE TwoMoment_RetainedArithmeticModule, ONLY: &
    RetainedReal, Retain, RetainedValue, &
    OPERATOR(+), OPERATOR(-), OPERATOR(*), ASSIGNMENT(=)
#ifdef MICROPHYSICS_WEAKLIB
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    InitGeneralScatteringKernels, FinalizeGeneralScatteringKernels, &
    CalculateAllPhout
#endif

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: iGL_NMS  = 1
  INTEGER, PARAMETER, PUBLIC :: iGL_LFE  = 2
  INTEGER, PARAMETER, PUBLIC :: iGL_LFC  = 3
  INTEGER, PARAMETER, PUBLIC :: iGL_IMD  = 4 ! Both muon charges, resolved by neutrino species.
  INTEGER, PARAMETER, PUBLIC :: nGL      = 4

  INTEGER, PARAMETER :: iProcessMin = 13
  INTEGER, PARAMETER :: iProcessMax = 34
  INTEGER, PARAMETER :: nProcesses  = iProcessMax - iProcessMin + 1
  INTEGER, PARAMETER :: iChannelNu  = 1
  INTEGER, PARAMETER :: iChannelNub = 2

  REAL(DP), PARAMETER :: UnitKernel = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: MaxExponent = LOG( HUGE( One ) ) - One
  REAL(DP), PARAMETER :: PhaseNormalization = One/(PlanckConstant*SpeedOfLight)**3

  LOGICAL :: Initialized = .FALSE.
  LOGICAL :: GeneralKernelWorkspaceInitialized = .FALSE.
  INTEGER :: nE_Global = 0, nX_Global = 0, nSpecies_Global = 0
  INTEGER :: GeneralKernelQuadratureOrder = 0
  INTEGER,  ALLOCATABLE :: IdentityIndex(:)
  REAL(DP), ALLOCATABLE :: EnergyGrid(:)
  REAL(DP), ALLOCATABLE :: GeneralKernelEnergyGrid(:)
  REAL(DP), ALLOCATABLE :: KernelTemperature(:)
  REAL(DP), ALLOCATABLE :: KernelMuElectron(:), KernelMuMuon(:)

  ! Phi_NMS stores the de-excitation kernel for the ordered pair
  ! (min(iE1,iE2),max(iE1,iE2)).  Flavor-changing and inverse-muon-decay
  ! kernels retain both forward and reverse directions explicitly because
  ! their detailed-balance factors contain charged-lepton chemical
  ! potentials.
  REAL(DP), ALLOCATABLE, PUBLIC :: Phi_NMS(:,:,:,:)
  REAL(DP), ALLOCATABLE :: Phi_LFE_F(:,:,:,:), Phi_LFE_R(:,:,:,:)
  REAL(DP), ALLOCATABLE :: Phi_LFC_F(:,:,:,:), Phi_LFC_R(:,:,:,:)
  REAL(DP), ALLOCATABLE :: Phi_IMD_F(:,:,:,:), Phi_IMD_R(:,:,:,:)

  REAL(DP), ALLOCATABLE, PUBLIC :: Chi_GL(:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: R_GL(:,:,:,:)
  ! Fluid-gain electron/muon number and comoving energy in code units.
  ! Accumulated from shared events, not from separately reduced species.
  REAL(DP), ALLOCATABLE, PUBLIC :: Source_GL(:,:,:)

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET( nE_Global, nX_Global, nSpecies_Global, &
  !$OMP   IdentityIndex, EnergyGrid, KernelTemperature, &
  !$OMP   KernelMuElectron, KernelMuMuon, &
  !$OMP   Phi_NMS, Phi_LFE_F, Phi_LFE_R, Phi_LFC_F, Phi_LFC_R, &
  !$OMP   Phi_IMD_F, Phi_IMD_R, Chi_GL, R_GL, Source_GL )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE( nE_Global, nX_Global, nSpecies_Global, &
  !$ACC   IdentityIndex, EnergyGrid, KernelTemperature, &
  !$ACC   KernelMuElectron, KernelMuMuon, &
  !$ACC   Phi_NMS, Phi_LFE_F, Phi_LFE_R, Phi_LFC_F, Phi_LFC_R, &
  !$ACC   Phi_IMD_F, Phi_IMD_R, Chi_GL, R_GL, Source_GL )
#endif

  PUBLIC :: InitializeNeutrinoLeptonInteractions
  PUBLIC :: FinalizeNeutrinoLeptonInteractions
  PUBLIC :: FinalizeNeutrinoLeptonInteractionKernels
  PUBLIC :: ComputeNeutrinoLeptonKernels
  PUBLIC :: ComputeNeutrinoLeptonRates
  PUBLIC :: ComputeNeutrinoLeptonPairingResiduals
  PUBLIC :: GetGeneralLeptonCoefficient

CONTAINS

  PURE SUBROUTINE GetGeneralLeptonCoefficient( Process, Species, R, Q, X, Gain, Loss )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(in) :: Process, Species, R, Q, X
    REAL(DP), INTENT(out) :: Gain, Loss
    Gain = Zero
    Loss = Zero
    SELECT CASE(Process)
      CASE(iGL_NMS)
        Loss = Phi_NMS(MIN(R,Q),MAX(R,Q),Species,X)
        Gain = Loss
      CASE(iGL_LFE)
        Gain = Phi_LFE_F(Q,R,Species,X)
        Loss = Phi_LFE_R(Q,R,Species,X)
      CASE(iGL_LFC)
        Gain = Phi_LFC_F(Q,R,Species,X)
        Loss = Phi_LFC_R(Q,R,Species,X)
      CASE(iGL_IMD)
        Gain = Phi_IMD_R(R,Q,3-Species,X)
        Loss = Phi_IMD_F(R,Q,3-Species,X)
    END SELECT
  END SUBROUTINE GetGeneralLeptonCoefficient

  SUBROUTINE InitializeNeutrinoLeptonInteractions &
    ( E, nX, nTheta_Option, EnableGeneralInteractions_Option )

    REAL(DP), INTENT(in) :: E(:)
    INTEGER,  INTENT(in) :: nX
    INTEGER,  INTENT(in), OPTIONAL :: nTheta_Option
    LOGICAL,  INTENT(in), OPTIONAL :: EnableGeneralInteractions_Option

    INTEGER :: iX, nTheta
    LOGICAL :: EnableGeneralInteractions

    IF( Initialized ) RETURN

    nE_Global = SIZE(E)
    nX_Global = nX
    nSpecies_Global = nSpecies
    nTheta = 16
    IF( PRESENT(nTheta_Option) ) nTheta = nTheta_Option
    EnableGeneralInteractions = .TRUE.
    IF( PRESENT(EnableGeneralInteractions_Option) ) &
      EnableGeneralInteractions = EnableGeneralInteractions_Option

    ALLOCATE( IdentityIndex(nX_Global), EnergyGrid(nE_Global) )
    EnergyGrid = E
    DO iX = 1, nX_Global
      IdentityIndex(iX) = iX
    END DO
    ALLOCATE( Chi_GL(nE_Global,nSpecies_Global,nX_Global,nGL) )
    ALLOCATE( R_GL  (nE_Global,nSpecies_Global,nX_Global,nGL) )
    ALLOCATE( Source_GL(3,nGL,nX_Global) )

    Chi_GL    = Zero; R_GL      = Zero
    Source_GL = Zero

    IF( EnableGeneralInteractions )THEN
      ALLOCATE( KernelTemperature(nX_Global) )
      ALLOCATE( KernelMuElectron(nX_Global), KernelMuMuon(nX_Global) )
      ALLOCATE( Phi_NMS  (nE_Global,nE_Global,nSpecies_Global,nX_Global) )
      ALLOCATE( Phi_LFE_F(nE_Global,nE_Global,2,nX_Global) )
      ALLOCATE( Phi_LFE_R(nE_Global,nE_Global,2,nX_Global) )
      ALLOCATE( Phi_LFC_F(nE_Global,nE_Global,2,nX_Global) )
      ALLOCATE( Phi_LFC_R(nE_Global,nE_Global,2,nX_Global) )
      ALLOCATE( Phi_IMD_F(nE_Global,nE_Global,2,nX_Global) )
      ALLOCATE( Phi_IMD_R(nE_Global,nE_Global,2,nX_Global) )
      KernelTemperature = Zero
      KernelMuElectron  = Zero
      KernelMuMuon      = Zero
      Phi_NMS   = Zero
      Phi_LFE_F = Zero; Phi_LFE_R = Zero
      Phi_LFC_F = Zero; Phi_LFC_R = Zero
      Phi_IMD_F = Zero; Phi_IMD_R = Zero
#ifdef MICROPHYSICS_WEAKLIB
      CALL EnsureGeneralKernelWorkspace( EnergyGrid, nTheta )
#endif
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( to: IdentityIndex, EnergyGrid, &
    !$OMP   nE_Global, nX_Global, nSpecies_Global, Chi_GL, R_GL, Source_GL )
    IF( EnableGeneralInteractions )THEN
      !$OMP TARGET ENTER DATA &
      !$OMP MAP( alloc: KernelTemperature, KernelMuElectron, KernelMuMuon ) &
      !$OMP MAP( alloc: Phi_NMS, Phi_LFE_F, Phi_LFE_R, &
      !$OMP             Phi_LFC_F, Phi_LFC_R, Phi_IMD_F, Phi_IMD_R )
    END IF
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( IdentityIndex, EnergyGrid, &
    !$ACC   nE_Global, nX_Global, nSpecies_Global, Chi_GL, R_GL, Source_GL )
#endif

    Initialized = .TRUE.

  END SUBROUTINE InitializeNeutrinoLeptonInteractions


  SUBROUTINE FinalizeNeutrinoLeptonInteractions

    IF( .NOT. Initialized ) RETURN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( delete: IdentityIndex, EnergyGrid, &
    !$OMP   Chi_GL, R_GL, Source_GL )
    IF( ALLOCATED(KernelTemperature) )THEN
      !$OMP TARGET EXIT DATA MAP( delete: KernelTemperature, &
      !$OMP   KernelMuElectron, KernelMuMuon, Phi_NMS, &
      !$OMP   Phi_LFE_F, Phi_LFE_R, Phi_LFC_F, Phi_LFC_R, &
      !$OMP   Phi_IMD_F, Phi_IMD_R )
    END IF
#elif defined(THORNADO_OACC)
#endif

    DEALLOCATE( IdentityIndex, EnergyGrid, Chi_GL, R_GL, Source_GL )
    IF( ALLOCATED(KernelTemperature) )THEN
      DEALLOCATE( KernelTemperature, KernelMuElectron, KernelMuMuon )
      DEALLOCATE( Phi_NMS, Phi_LFE_F, Phi_LFE_R )
      DEALLOCATE( Phi_LFC_F, Phi_LFC_R, Phi_IMD_F, Phi_IMD_R )
    END IF
    Initialized = .FALSE.
    nE_Global = 0; nX_Global = 0; nSpecies_Global = 0

  END SUBROUTINE FinalizeNeutrinoLeptonInteractions


  SUBROUTINE EnsureGeneralKernelWorkspace( E, nTheta )

    REAL(DP), INTENT(in) :: E(:)
    INTEGER,  INTENT(in) :: nTheta

#ifdef MICROPHYSICS_WEAKLIB
    LOGICAL :: ConfigurationChanged

    ConfigurationChanged = GeneralKernelWorkspaceInitialized
    IF( GeneralKernelWorkspaceInitialized )THEN
      ConfigurationChanged &
        =      GeneralKernelQuadratureOrder /= nTheta &
          .OR. SIZE(GeneralKernelEnergyGrid) /= SIZE(E)
      IF( .NOT. ConfigurationChanged ) &
        ConfigurationChanged = ANY( GeneralKernelEnergyGrid /= E )
    END IF

    IF( ConfigurationChanged )THEN
      CALL FinalizeGeneralScatteringKernels
      DEALLOCATE( GeneralKernelEnergyGrid )
      GeneralKernelWorkspaceInitialized = .FALSE.
      GeneralKernelQuadratureOrder = 0
    END IF

    IF( .NOT. GeneralKernelWorkspaceInitialized )THEN
      ALLOCATE( GeneralKernelEnergyGrid(SIZE(E)) )
      GeneralKernelEnergyGrid = E
      GeneralKernelQuadratureOrder = nTheta
      CALL InitGeneralScatteringKernels &
             ( E / MeV, E / MeV, SIZE(E), SIZE(E), nTheta, &
               iProcessMin_Option = iProcessMin, &
               iProcessMax_Option = iProcessMax )
      GeneralKernelWorkspaceInitialized = .TRUE.
    END IF
#endif

  END SUBROUTINE EnsureGeneralKernelWorkspace


  SUBROUTINE FinalizeNeutrinoLeptonInteractionKernels

#ifdef MICROPHYSICS_WEAKLIB
    IF( GeneralKernelWorkspaceInitialized )THEN
      CALL FinalizeGeneralScatteringKernels
      DEALLOCATE( GeneralKernelEnergyGrid )
      GeneralKernelWorkspaceInitialized = .FALSE.
      GeneralKernelQuadratureOrder = 0
    END IF
#endif

  END SUBROUTINE FinalizeNeutrinoLeptonInteractionKernels


  SUBROUTINE ComputeNeutrinoLeptonKernels &
    ( D, T, Ye, Ym, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
      IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
      nActive, UnpackIndex )

    REAL(DP), INTENT(in) :: D(:), T(:), Ye(:), Ym(:)
    LOGICAL,  INTENT(in) :: ActiveNMS(:), ActiveLFE(:)
    LOGICAL,  INTENT(in) :: ActiveLFC(:), ActiveIMD(:)
    LOGICAL,  INTENT(in) :: IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD
    INTEGER,  INTENT(in) :: nActive
    INTEGER,  INTENT(in), OPTIONAL :: UnpackIndex(:)

    IF( PRESENT(UnpackIndex) )THEN
      CALL ComputeNeutrinoLeptonKernels_Mapped &
             ( D, T, Ye, Ym, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
               IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
               nActive, UnpackIndex )
    ELSE
      CALL ComputeNeutrinoLeptonKernels_Mapped &
             ( D, T, Ye, Ym, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
               IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
               nActive, IdentityIndex )
    END IF

  END SUBROUTINE ComputeNeutrinoLeptonKernels


  SUBROUTINE ComputeNeutrinoLeptonKernels_Mapped &
    ( D, T, Ye, Ym, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
      IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
      nActive, UnpackIndex )

    REAL(DP), INTENT(in) :: D(:), T(:), Ye(:), Ym(:)
    LOGICAL,  INTENT(in) :: ActiveNMS(:), ActiveLFE(:)
    LOGICAL,  INTENT(in) :: ActiveLFC(:), ActiveIMD(:)
    LOGICAL,  INTENT(in) :: IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD
    INTEGER,  INTENT(in) :: nActive
    INTEGER,  INTENT(in) :: UnpackIndex(:)

    INTEGER  :: iActive, iE1, iE3, iS, iX
    INTEGER  :: iP13, iP25, iP27, iP29, iP31, iP33, iP34
    REAL(DP) :: E1, E3, kT, MuE, MuM, Scale
    REAL(DP) :: Phout(nProcesses,1)

    iP13 = 13 - iProcessMin + 1
    iP25 = 25 - iProcessMin + 1
    iP27 = 27 - iProcessMin + 1
    iP29 = 29 - iProcessMin + 1
    iP31 = 31 - iProcessMin + 1
    iP33 = 33 - iProcessMin + 1
    iP34 = 34 - iProcessMin + 1
    Scale = FourPi * UnitKernel

    ! Kernel construction is deliberately fused across all requested
    ! reactions.  CalculateAllPhout evaluates the distinct charged-lepton
    ! phase-space integrals once and maps them to every process.
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE( iX, MuE, MuM ) &
    !$OMP FIRSTPRIVATE( nActive, Scale, IncludeNMS, IncludeLFE, &
    !$OMP               IncludeLFC, IncludeIMD )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( iX, MuE, MuM ) &
    !$ACC FIRSTPRIVATE( nActive, Scale, IncludeNMS, IncludeLFE, &
    !$ACC               IncludeLFC, IncludeIMD ) &
    !$ACC PRESENT( D, T, Ye, Ym, ActiveNMS, ActiveLFE, ActiveLFC, &
    !$ACC          ActiveIMD, UnpackIndex, KernelTemperature, &
    !$ACC          KernelMuElectron, KernelMuMuon )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iX, MuE, MuM ) &
    !$OMP FIRSTPRIVATE( nActive, Scale, IncludeNMS, IncludeLFE, &
    !$OMP               IncludeLFC, IncludeIMD )
#endif
    DO iActive = 1, nActive
      iX = UnpackIndex(iActive)

      IF( .NOT. ( (IncludeNMS .AND. ActiveNMS(iX)) .OR. &
                  (IncludeLFE .AND. ActiveLFE(iX)) .OR. &
                  (IncludeLFC .AND. ActiveLFC(iX)) .OR. &
                  (IncludeIMD .AND. ActiveIMD(iX)) ) )THEN
        Phi_NMS(:,:,:,iX)   = Zero
        Phi_LFE_F(:,:,:,iX) = Zero; Phi_LFE_R(:,:,:,iX) = Zero
        Phi_LFC_F(:,:,:,iX) = Zero; Phi_LFC_R(:,:,:,iX) = Zero
        Phi_IMD_F(:,:,:,iX) = Zero; Phi_IMD_R(:,:,:,iX) = Zero
        KernelTemperature(iX) = Zero
        KernelMuElectron (iX) = Zero
        KernelMuMuon     (iX) = Zero
        CYCLE
      END IF

      CALL ComputeElectronChemicalPotential_TABLE &
             ( D(iX), T(iX), Ye(iX), Ym(iX), MuE )
      CALL ComputeMuonChemicalPotential_TABLE &
             ( D(iX), T(iX), Ye(iX), Ym(iX), MuM )
      KernelMuElectron (iX) = MuE / MeV
      KernelMuMuon     (iX) = MuM / MeV
      KernelTemperature(iX) = BoltzmannConstant * T(iX) / MeV
    END DO

    ! The thermodynamic work above is point-local.  The costly angular and
    ! charged-lepton phase-space integrals are independent for every
    ! (point,E1,E3) tuple and must expose that full parallelism to accelerators.
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iX, iS, E1, E3, kT, MuE, MuM, Phout ) &
    !$OMP FIRSTPRIVATE( nActive, Scale, IncludeNMS, IncludeLFE, &
    !$OMP               IncludeLFC, IncludeIMD )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iX, iS, E1, E3, kT, MuE, MuM, Phout ) &
    !$ACC FIRSTPRIVATE( nActive, Scale, IncludeNMS, IncludeLFE, &
    !$ACC               IncludeLFC, IncludeIMD ) &
    !$ACC PRESENT( ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
    !$ACC          UnpackIndex, EnergyGrid, KernelTemperature, &
    !$ACC          KernelMuElectron, KernelMuMuon, Phi_NMS, &
    !$ACC          Phi_LFE_F, Phi_LFE_R, Phi_LFC_F, Phi_LFC_R, &
    !$ACC          Phi_IMD_F, Phi_IMD_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iX, iS, E1, E3, kT, MuE, MuM, Phout ) &
    !$OMP FIRSTPRIVATE( nActive, Scale, IncludeNMS, IncludeLFE, &
    !$OMP               IncludeLFC, IncludeIMD )
#endif
    DO iActive = 1, nActive
    DO iE3 = 1, nE_Global
    DO iE1 = 1, nE_Global
        iX = UnpackIndex(iActive)
        kT  = KernelTemperature(iX)
        MuE = KernelMuElectron(iX)
        MuM = KernelMuMuon(iX)
        IF( kT <= Zero ) CYCLE
        E1 = EnergyGrid(iE1) / MeV
        E3 = EnergyGrid(iE3) / MeV
#ifdef MICROPHYSICS_WEAKLIB
        CALL CalculateAllPhout &
               ( iE1, iE3, E1, E3, kT, MuE, MuM, Phout, 1 )
#else
        Phout = Zero
#endif

        IF( IncludeNMS .AND. ActiveNMS(iX) )THEN
          DO iS = 1, nSpecies_Global
            IF( iE1 <= iE3 ) &
              Phi_NMS(iE1,iE3,iS,iX) = Scale &
                * ( Phout(iP13 + 2 * (iS-1),1) &
                    + Phout(iP13 + 2 * (iS-1) + 1,1) )
          END DO
        ELSE IF( iE1 <= iE3 )THEN
          Phi_NMS(iE1,iE3,:,iX) = Zero
        END IF

        IF( IncludeLFE .AND. ActiveLFE(iX) )THEN
          Phi_LFE_F(iE1,iE3,iChannelNu, iX) = Scale * Phout(iP25,1)
          Phi_LFE_F(iE1,iE3,iChannelNub,iX) = Scale * Phout(iP27,1)
          Phi_LFE_R(iE1,iE3,iChannelNu, iX) = Scale &
            * ScaleByExponential &
                ( Phout(iP25,1), ( E3 - E1 - MuE + MuM ) / kT )
          Phi_LFE_R(iE1,iE3,iChannelNub,iX) = Scale &
            * ScaleByExponential &
                ( Phout(iP27,1), ( E3 - E1 + MuE - MuM ) / kT )
        ELSE
          Phi_LFE_F(iE1,iE3,:,iX) = Zero
          Phi_LFE_R(iE1,iE3,:,iX) = Zero
        END IF

        IF( IncludeLFC .AND. ActiveLFC(iX) )THEN
          Phi_LFC_F(iE1,iE3,iChannelNu, iX) = Scale * Phout(iP29,1)
          Phi_LFC_F(iE1,iE3,iChannelNub,iX) = Scale * Phout(iP31,1)
          Phi_LFC_R(iE1,iE3,iChannelNu, iX) = Scale &
            * ScaleByExponential &
                ( Phout(iP29,1), ( E3 - E1 + MuM - MuE ) / kT )
          Phi_LFC_R(iE1,iE3,iChannelNub,iX) = Scale &
            * ScaleByExponential &
                ( Phout(iP31,1), ( E3 - E1 - MuM + MuE ) / kT )
        ELSE
          Phi_LFC_F(iE1,iE3,:,iX) = Zero
          Phi_LFC_R(iE1,iE3,:,iX) = Zero
        END IF

        IF( IncludeIMD .AND. ActiveIMD(iX) )THEN
          Phi_IMD_F(iE1,iE3,iChannelNu, iX) = Scale * Phout(iP33,1)
          Phi_IMD_F(iE1,iE3,iChannelNub,iX) = Scale * Phout(iP34,1)
          Phi_IMD_R(iE1,iE3,iChannelNu, iX) = Scale &
            * ScaleByExponential &
                ( Phout(iP33,1), ( MuM - MuE - E1 - E3 ) / kT )
          Phi_IMD_R(iE1,iE3,iChannelNub,iX) = Scale &
            * ScaleByExponential &
                ( Phout(iP34,1), ( MuE - MuM - E1 - E3 ) / kT )
        ELSE
          Phi_IMD_F(iE1,iE3,:,iX) = Zero
          Phi_IMD_R(iE1,iE3,:,iX) = Zero
        END IF
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoLeptonKernels_Mapped


  SUBROUTINE ComputeNeutrinoLeptonRates &
    ( T, W2, Dnu, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
      IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
      nActive, UnpackIndex )

    REAL(DP), INTENT(in) :: T(:), W2(:), Dnu(:,:,:)
    LOGICAL,  INTENT(in) :: ActiveNMS(:), ActiveLFE(:)
    LOGICAL,  INTENT(in) :: ActiveLFC(:), ActiveIMD(:)
    LOGICAL,  INTENT(in) :: IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD
    INTEGER,  INTENT(in) :: nActive
    INTEGER,  INTENT(in), OPTIONAL :: UnpackIndex(:)

    IF( PRESENT(UnpackIndex) )THEN
      CALL ComputeNeutrinoLeptonRates_Mapped &
             ( T, W2, Dnu, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
               IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
               nActive, UnpackIndex )
    ELSE
      CALL ComputeNeutrinoLeptonRates_Mapped &
             ( T, W2, Dnu, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
               IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
               nActive, IdentityIndex )
    END IF

  END SUBROUTINE ComputeNeutrinoLeptonRates


  SUBROUTINE ComputeNeutrinoLeptonRates_Mapped &
    ( T, W2, Dnu, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
      IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD, &
      nActive, UnpackIndex )

    REAL(DP), INTENT(in) :: T(:), W2(:), Dnu(:,:,:)
    LOGICAL,  INTENT(in) :: ActiveNMS(:), ActiveLFE(:)
    LOGICAL,  INTENT(in) :: ActiveLFC(:), ActiveIMD(:)
    LOGICAL,  INTENT(in) :: IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD
    INTEGER,  INTENT(in) :: nActive
    INTEGER,  INTENT(in) :: UnpackIndex(:)

    INTEGER :: iActive, iX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE( iX ) FIRSTPRIVATE( nActive, IncludeNMS, IncludeLFE, &
    !$OMP   IncludeLFC, IncludeIMD )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG PRIVATE( iX ) &
    !$ACC FIRSTPRIVATE( nActive, IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD ) &
    !$ACC PRESENT( T, W2, Dnu, ActiveNMS, ActiveLFE, ActiveLFC, ActiveIMD, &
    !$ACC          UnpackIndex, &
    !$ACC          Chi_GL, R_GL, Source_GL, Phi_NMS, Phi_LFE_F, Phi_LFE_R, &
    !$ACC          Phi_LFC_F, Phi_LFC_R, Phi_IMD_F, Phi_IMD_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO PRIVATE( iX ) &
    !$OMP FIRSTPRIVATE( nActive, IncludeNMS, IncludeLFE, IncludeLFC, IncludeIMD )
#endif
    DO iActive = 1, nActive
      iX = UnpackIndex(iActive)
      Chi_GL(:,:,iX,:) = Zero
      R_GL  (:,:,iX,:) = Zero
      Source_GL(:,:,iX) = Zero

      IF( IncludeNMS .AND. ActiveNMS(iX) ) &
        CALL ComputeNMSPoint( iX, T(iX), W2, Dnu )
      IF( IncludeLFE .AND. ActiveLFE(iX) ) &
        CALL ComputeFlavorPoint &
               ( iX, W2, Dnu, Phi_LFE_F, Phi_LFE_R, iGL_LFE )
      IF( IncludeLFC .AND. ActiveLFC(iX) ) &
        CALL ComputeFlavorPoint &
               ( iX, W2, Dnu, Phi_LFC_F, Phi_LFC_R, iGL_LFC )
      IF( IncludeIMD .AND. ActiveIMD(iX) ) &
        CALL ComputeIMDPoint( iX, W2, Dnu )
    END DO

  END SUBROUTINE ComputeNeutrinoLeptonRates_Mapped


  SUBROUTINE ComputeNMSPoint( iX, T, W2, Dnu )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iX
    REAL(DP), INTENT(in) :: T, W2(:), Dnu(:,:,:)
    INTEGER :: iE1, iE2, iS
    REAL(DP) :: DetBal, PhiIn, PhiOut, SumChi, SumR, Event
    REAL(DP) :: EnergySource, EnergyCorrection

    EnergySource = Zero
    EnergyCorrection = Zero
    DO iS = 1, nSpecies_Global
    DO iE2 = 1, nE_Global
      SumChi = Zero; SumR = Zero
      DO iE1 = 1, nE_Global
        DetBal = SafeExponential &
          ( -ABS(EnergyGrid(iE2)-EnergyGrid(iE1)) &
              / (BoltzmannConstant*T) )
        IF( iE1 <= iE2 )THEN
          PhiOut = Phi_NMS(iE1,iE2,iS,iX)
          PhiIn  = PhiOut * DetBal
        ELSE
          PhiIn  = Phi_NMS(iE2,iE1,iS,iX)
          PhiOut = PhiIn * DetBal
        END IF
        SumChi = SumChi + W2(iE1) &
          * ( PhiIn*Dnu(iE1,iS,iX) &
              + PhiOut*(One-Dnu(iE1,iS,iX)) )
        Event = PhiIn*Dnu(iE1,iS,iX)*(One-Dnu(iE2,iS,iX)) &
                  - PhiOut*(One-Dnu(iE1,iS,iX))*Dnu(iE2,iS,iX)
        SumR = SumR + W2(iE1)*Event
        CALL AddInteractionSource &
               ( -0.5_DP*(EnergyGrid(iE2)-EnergyGrid(iE1)) &
                   * ((W2(iE2)*PhaseNormalization)*W2(iE1))*Event, &
                 EnergySource, EnergyCorrection )
      END DO
      Chi_GL(iE2,iS,iX,iGL_NMS) = SumChi
      R_GL  (iE2,iS,iX,iGL_NMS) = SumR
    END DO
    END DO
    Source_GL(3,iGL_NMS,iX) = EnergySource + EnergyCorrection

  END SUBROUTINE ComputeNMSPoint


  SUBROUTINE ComputeFlavorPoint &
    ( iX, W2, Dnu, PhiF, PhiR, iProcess )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iX, iProcess
    REAL(DP), INTENT(in) :: W2(:), Dnu(:,:,:)
    REAL(DP), INTENT(in) :: PhiF(:,:,:,:), PhiR(:,:,:,:)
    INTEGER :: iChannel, iEa, iEb, iSa, iSb
    REAL(DP) :: Event, Kf, Kr, Fa, Fb, NumberEvent, Sign
    REAL(DP) :: NumberSource, NumberCorrection, EnergySource, EnergyCorrection

    NumberSource = Zero; NumberCorrection = Zero
    EnergySource = Zero; EnergyCorrection = Zero
    DO iChannel = 1, 2
      Sign = REAL(3-2*iChannel,DP)
      IF( iChannel == iChannelNu )THEN
        iSa = iNuM;     iSb = iNuE
      ELSE
        iSa = iNuM_Bar; iSb = iNuE_Bar
      END IF
      IF( iSa > nSpecies_Global .OR. iSb > nSpecies_Global ) CYCLE

      DO iEa = 1, nE_Global
        Fa = Dnu(iEa,iSa,iX)
        DO iEb = 1, nE_Global
          Fb = Dnu(iEb,iSb,iX)
          Kf = PhiF(iEa,iEb,iChannel,iX)
          Kr = PhiR(iEa,iEb,iChannel,iX)
          Event = Kf*Fa*(One-Fb) - Kr*Fb*(One-Fa)
          NumberEvent = ((W2(iEa)*PhaseNormalization)*W2(iEb))*Event
          CALL AddInteractionSource &
                 ( Sign*NumberEvent, NumberSource, NumberCorrection )
          CALL AddInteractionSource &
                 ( (EnergyGrid(iEa)-EnergyGrid(iEb))*NumberEvent, &
                   EnergySource, EnergyCorrection )
          R_GL(iEa,iSa,iX,iProcess) = R_GL(iEa,iSa,iX,iProcess) &
            - W2(iEb)*Event
          R_GL(iEb,iSb,iX,iProcess) = R_GL(iEb,iSb,iX,iProcess) &
            + W2(iEa)*Event
          Chi_GL(iEa,iSa,iX,iProcess) = Chi_GL(iEa,iSa,iX,iProcess) &
            + W2(iEb)*( Kf*(One-Fb) + Kr*Fb )
          Chi_GL(iEb,iSb,iX,iProcess) = Chi_GL(iEb,iSb,iX,iProcess) &
            + W2(iEa)*( Kr*(One-Fa) + Kf*Fa )
        END DO
      END DO
    END DO
    Source_GL(2,iProcess,iX) = NumberSource + NumberCorrection
    Source_GL(1,iProcess,iX) = -Source_GL(2,iProcess,iX)
    Source_GL(3,iProcess,iX) = EnergySource + EnergyCorrection

  END SUBROUTINE ComputeFlavorPoint


  SUBROUTINE ComputeIMDPoint( iX, W2, Dnu )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iX
    REAL(DP), INTENT(in) :: W2(:), Dnu(:,:,:)
    INTEGER :: iChannel, iE1, iE3, iS1, iS3
    REAL(DP) :: Sign
    TYPE(RetainedReal) :: Event, Kf, Kr, F1, F3, NumberEvent
    TYPE(RetainedReal) :: NumberSource, EnergySource, Rate1, Chi1
    TYPE(RetainedReal) :: Rate3(nE_Global), Chi3(nE_Global)

    NumberSource = Zero
    EnergySource = Zero
    DO iChannel = 1, 2
      Sign = REAL(3-2*iChannel,DP)
      IF( iChannel == iChannelNu )THEN
        iS1 = iNuE_Bar; iS3 = iNuM
      ELSE
        iS1 = iNuE;     iS3 = iNuM_Bar
      END IF
      IF( iS1 > nSpecies_Global .OR. iS3 > nSpecies_Global ) CYCLE
      Rate3 = Zero
      Chi3 = Zero
      DO iE1 = 1, nE_Global
        F1 = Dnu(iE1,iS1,iX)
        Rate1 = Zero
        Chi1 = Zero
        DO iE3 = 1, nE_Global
          F3 = Dnu(iE3,iS3,iX)
          Kf = Phi_IMD_F(iE1,iE3,iChannel,iX)
          Kr = Phi_IMD_R(iE1,iE3,iChannel,iX)
          Event = Kf*F1*F3 - Kr*(One-F1)*(One-F3)
          NumberEvent = (Retain(W2(iE1))*PhaseNormalization)*W2(iE3)*Event
          NumberSource = NumberSource + Sign*NumberEvent
          EnergySource = EnergySource &
            + (Retain(EnergyGrid(iE1))+EnergyGrid(iE3))*NumberEvent
          ! One retained event supplies both species and their fluid source.
          Rate1 = Rate1 - W2(iE3)*Event
          Rate3(iE3) = Rate3(iE3) - W2(iE1)*Event
          Chi1 = Chi1 + W2(iE3)*( Kf*F3 + Kr*(One-F3) )
          Chi3(iE3) = Chi3(iE3) + W2(iE1)*( Kf*F1 + Kr*(One-F1) )
        END DO
        R_GL(iE1,iS1,iX,iGL_IMD) = RetainedValue(Rate1)
        Chi_GL(iE1,iS1,iX,iGL_IMD) = RetainedValue(Chi1)
      END DO
      DO iE3 = 1, nE_Global
        R_GL(iE3,iS3,iX,iGL_IMD) = RetainedValue(Rate3(iE3))
        Chi_GL(iE3,iS3,iX,iGL_IMD) = RetainedValue(Chi3(iE3))
      END DO
    END DO
    Source_GL(2,iGL_IMD,iX) = RetainedValue(NumberSource)
    Source_GL(1,iGL_IMD,iX) = -Source_GL(2,iGL_IMD,iX)
    Source_GL(3,iGL_IMD,iX) = RetainedValue(EnergySource)

  END SUBROUTINE ComputeIMDPoint


  PURE SUBROUTINE AddInteractionSource( Value, Total, Correction )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: Value
    REAL(DP), INTENT(inout) :: Total, Correction
    REAL(DP) :: Next

    Next = Total + Value
    IF( ABS(Total) >= ABS(Value) )THEN
      Correction = Correction + ((Total-Next)+Value)
    ELSE
      Correction = Correction + ((Value-Next)+Total)
    END IF
    Total = Next

  END SUBROUTINE AddInteractionSource


  SUBROUTINE ComputeNeutrinoLeptonPairingResiduals( W2, Residual, Scale )

    REAL(DP), INTENT(in)  :: W2(:)
    REAL(DP), INTENT(out) :: Residual(:,:), Scale(:,:)

    INTEGER :: iE, iProcess, iS, iX, nX_Output
    REAL(DP) :: A, B

    nX_Output = SIZE( Residual, 2 )

    IF( SIZE(W2)         /= nE_Global .OR. &
        SIZE(Residual,1) /= nGL .OR. &
        SIZE(Scale,1)    /= nGL .OR. &
        SIZE(Scale,2)    /= nX_Output .OR. &
        nX_Output > nX_Global ) &
      ERROR STOP 'Invalid general-lepton pairing diagnostic shape.'

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
    !$OMP PRIVATE( iE, iProcess, iS, A, B ) MAP( to: W2 ) &
    !$OMP MAP( from: Residual, Scale )
#elif defined(THORNADO_OACC)
    !$ACC DATA COPYIN( W2 ) CREATE( Residual, Scale )
    !$ACC PARALLEL LOOP GANG PRIVATE( iE, iProcess, iS, A, B ) &
    !$ACC PRESENT( W2, R_GL, Residual, Scale )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO PRIVATE( iE, iProcess, iS, A, B )
#endif
    DO iX = 1, nX_Output

      Residual(:,iX) = Zero
      Scale   (:,iX) = Zero

      ! NMS must conserve the number of every neutrino species.
      DO iS = 1, nSpecies_Global
        A = Zero; B = Zero
        DO iE = 1, nE_Global
          A = A + W2(iE) * R_GL(iE,iS,iX,iGL_NMS)
          B = B + ABS( W2(iE) * R_GL(iE,iS,iX,iGL_NMS) )
        END DO
        IF( ABS(A) > ABS(Residual(iGL_NMS,iX)) )THEN
          Residual(iGL_NMS,iX) = A
          Scale   (iGL_NMS,iX) = B
        END IF
      END DO

      IF( nSpecies_Global >= iNuM_Bar )THEN
        ! LFE and LFC each move one neutrino between electron and muon
        ! flavor.  Check neutrino and antineutrino channels separately so
        ! cancellation between charge-conjugate reactions cannot hide an
        ! error.
        DO iProcess = iGL_LFE, iGL_LFC
          A = Zero; B = Zero
          DO iE = 1, nE_Global
            A = A + W2(iE) &
                  * ( R_GL(iE,iNuE,iX,iProcess) &
                      + R_GL(iE,iNuM,iX,iProcess) )
            B = B + ABS(W2(iE)*R_GL(iE,iNuE,iX,iProcess)) &
                  + ABS(W2(iE)*R_GL(iE,iNuM,iX,iProcess))
          END DO
          Residual(iProcess,iX) = A
          Scale   (iProcess,iX) = B

          A = Zero; B = Zero
          DO iE = 1, nE_Global
            A = A + W2(iE) &
                  * ( R_GL(iE,iNuE_Bar,iX,iProcess) &
                      + R_GL(iE,iNuM_Bar,iX,iProcess) )
            B = B + ABS(W2(iE)*R_GL(iE,iNuE_Bar,iX,iProcess)) &
                  + ABS(W2(iE)*R_GL(iE,iNuM_Bar,iX,iProcess))
          END DO
          IF( ABS(A) > ABS(Residual(iProcess,iX)) )THEN
            Residual(iProcess,iX) = A
            Scale   (iProcess,iX) = B
          END IF
        END DO

        ! Each IMD charge channel preserves L_e+L_mu.  Again check the two
        ! channels separately.
        A = Zero; B = Zero
        DO iE = 1, nE_Global
          A = A + W2(iE) &
                * ( -R_GL(iE,iNuE_Bar,iX,iGL_IMD) &
                    +R_GL(iE,iNuM,    iX,iGL_IMD) )
          B = B + ABS(W2(iE)*R_GL(iE,iNuE_Bar,iX,iGL_IMD)) &
                + ABS(W2(iE)*R_GL(iE,iNuM,iX,iGL_IMD))
        END DO
        Residual(iGL_IMD,iX) = A
        Scale   (iGL_IMD,iX) = B

        A = Zero; B = Zero
        DO iE = 1, nE_Global
          A = A + W2(iE) &
                * ( R_GL(iE,iNuE,    iX,iGL_IMD) &
                    -R_GL(iE,iNuM_Bar,iX,iGL_IMD) )
          B = B + ABS(W2(iE)*R_GL(iE,iNuE,iX,iGL_IMD)) &
                + ABS(W2(iE)*R_GL(iE,iNuM_Bar,iX,iGL_IMD))
        END DO
        IF( ABS(A) > ABS(Residual(iGL_IMD,iX)) )THEN
          Residual(iGL_IMD,iX) = A
          Scale   (iGL_IMD,iX) = B
        END IF
      END IF

    END DO
#if defined(THORNADO_OACC)
    !$ACC UPDATE HOST( Residual, Scale )
    !$ACC END DATA
#endif

  END SUBROUTINE ComputeNeutrinoLeptonPairingResiduals


  PURE FUNCTION SafeExponential( X ) RESULT( Value )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: X
    REAL(DP) :: Value

    IF( X <= LOG(TINY(One)) )THEN
      Value = Zero
    ELSE
      Value = EXP( MIN(X,MaxExponent) )
    END IF

  END FUNCTION SafeExponential


  PURE FUNCTION ScaleByExponential( Scale, X ) RESULT( Value )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: Scale, X
    REAL(DP) :: LogValue, Value

    IF( Scale <= Zero )THEN
      Value = Zero
      RETURN
    END IF

    LogValue = LOG( Scale ) + X
    IF( LogValue <= LOG(TINY(One)) )THEN
      Value = Zero
    ELSE
      Value = EXP( MIN(LogValue,MaxExponent) )
    END IF

  END FUNCTION ScaleByExponential

END MODULE NeutrinoLeptonInteractionsModule
