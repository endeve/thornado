PROGRAM NeutrinoOpacities

  USE KindModule, ONLY: &
    DP, FourPi
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV, &
    BoltzmannConstant, &
    PlanckConstant, &
    AtomicMassUnit, &
    Second
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    WeightsE
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE, &
    C1, C2, C1_NuPair, C2_NuPair
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar, iNuM, iNuM_Bar, &
    nChirals, LeptonNumber
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, &
    ComputeEquilibriumDistributions_DG, &
    ComputeNeutrinoOpacities_EC, &
    ComputeNeutrinoOpacities_ES, &
    ComputeNeutrinoOpacities_NES, &
    ComputeNeutrinoOpacityRates_NES, &
    ComputeNeutrinoOpacityRates_LinearCorrections_NES, &
    ComputeNeutrinoOpacities_Pair, &
    ComputeNeutrinoOpacityRates_Pair, &
    ComputeNeutrinoOpacityRates_LinearCorrections_Pair, &
    ComputeNeutrinoOpacities_NuPair, &
    ComputeNeutrinoOpacityRates_NuPair, &
    ComputeNeutrinoOpacities_Brem, &
    ComputeNeutrinoOpacityRates_Brem
  USE DeviceModule, ONLY: &
    InitializeDevice, &
    FinalizeDevice
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop
  USE OpacityModule_TABLE, ONLY: &
#ifdef MICROPHYSICS_WEAKLIB
    EmAb_EC_spec_T, Es_T, Ds_T, Ts_T, Ys_T, &
    OS_EmAb_EC_rate, EmAb_EC_rate_T, &
    Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T, EC_nE, &
    use_EC_table
#endif

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
    GetIndexAndDelta_Lin, &
    GetIndexAndDelta_Log

  ! ----------------------------------------------

#endif

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: &
    nNodes   = 2, & !2, &
    nE       = 16, & !2**4, &
    nX1      = 2**6, &
    nPointsX = nNodes * nX1, &
    nPointsE = nNodes * nE, &
    nSpecies = 6
  REAL(DP), PARAMETER :: &
    Unit_D     = Gram / Centimeter**3, &
    Unit_T     = Kelvin, &
    Unit_Y     = 1.0_DP, &
    Unit_E     = MeV, &
    Unit_Chi   = 1.0_DP / Centimeter, &
    Unit_Sigma = 1.0_DP / Centimeter, &
    UnitNES    = 1.0_DP / ( Centimeter * MeV**3 ), &
    UnitPair   = 1.0_DP / ( Centimeter * MeV**3 ), &
    BaryonMass = AtomicMassUnit, &
    Unit_Qdot  = 1.0_DP / ( BaryonMass * Second ), &
    Unit_Edot  = 1.0_DP / ( BaryonMass * Second ), &
    eL         = 0.0e0_DP * Unit_E, &
    eR         = 3.0e2_DP * Unit_E, &
    !ZoomE      = 1.183081754893913_DP
    ZoomE      = 1.266038160710160d0

  INTEGER :: &
    mpierr, iE, iX, iS, iNodeE, iN_E, iE1, iE2, iC
  REAL(DP) :: &
    kT, DetBal, &
    Timer_ReadEos, &
    Timer_ReadOpacities, &
    Timer_ComputeEquilibrium, &
    Timer_ComputeEquilibrium_DG, &
    Timer_Compute_EC, &
    Timer_Compute_ES, &
    Timer_ComputeKrnl_NES, &
    Timer_ComputeRate_NES, &
    Timer_Compute_NES, &
    Timer_ComputeCorr_NES, &
    Timer_ComputeKrnl_Pair, &
    Timer_ComputeRate_Pair, &
    Timer_Compute_Pair, &
    Timer_ComputeCorr_Pair, &
    Timer_ComputeKrnl_NuPair, &
    Timer_ComputeRate_NuPair, &
    Timer_Compute_NuPair, &
    Timer_ComputeKrnl_Brem, &
    Timer_ComputeRate_Brem, &
    Timer_Compute_Brem, &
    Timer_Total
  REAL(DP), DIMENSION(nPointsX) :: &
    D, T, Y
  REAL(DP), DIMENSION(nE) :: &
    dE
  REAL(DP), DIMENSION(nE,nSpecies,nPointsX) :: Edot_EmAb_element, &
    Edot_Iso_element,  &
    Qdot_Pair_element, &
    Qdot_NuPair_element, &
    Qdot_Brem_element, &
    Edot_NES_element
  REAL(DP), DIMENSION(nPointsE) :: &
    E, W2
  REAL(DP), DIMENSION(nPointsE,nChirals,nPointsX) :: &
    Sigma_Iso    ! --- Iso-energertic Kernel
  REAL(DP), DIMENSION(nPointsE,nSpecies,nPointsX) :: &
    f0       , & ! --- Equilibrium Distribution
    f0_DG    , & ! --- Equilibrium Distribution (DG Approximation)
    J        , & ! --- Neutrino Number Density
    Zero_J   , & ! --- All zeros
    H_1      , & ! --- Neutrino Flux (x-direction)
    H_2      , & ! --- Neutrino Flux (y-direction)
    H_3      , & ! --- Neutrino Flux (z-direction)
    Eta_EmAb , & ! --- Electron Capture Emissivity
    Chi_EmAb , & ! --- Electron Capture Opacity
    Edot_EmAb, & ! --- FourPi / PlanckConstant**3 / D * |Chi_EmAb| * E**3
    Eta_Iso  , & ! --- Iso Emissivity
    Chi_Iso  , & ! --- Iso Opacity
    Edot_Iso , & ! --- FourPi / PlanckConstant**3 / D * Chi_Iso * E**3
    Eta_NES  , & ! --- NES Emissivity
    Chi_NES  , & ! --- NES Opacity
    Edot_NES , & ! --- FourPi / PlanckConstant**3 / D * Chi_NES * E**3
    Eta_Pair , & ! --- Pair Emissivity
    Chi_Pair , & ! --- Pair Opacity
    Qdot_Pair, & ! --- Pair Heating Rate
    Eta_NuPair , & ! --- electron neutrino Pair annihilation Emissivity
    Chi_NuPair , & ! --- electron neutrino Pair annihilation Opacity
    Qdot_NuPair, & ! --- electron neutrino Pair annihilation Heating Rate
    Eta_Brem , & ! --- Brem Emissivity
    Chi_Brem , & ! --- Brem Opacity
    Qdot_Brem, & ! --- Brem Heating Rate
    A_In_1   , & ! --- Linear Corrections to NES rates
    A_In_2   , & ! --- Linear Corrections to NES rates
    A_In_3   , & ! --- Linear Corrections to NES rates
    A_Out_1  , & ! --- Linear Corrections to NES rates
    A_Out_2  , & ! --- Linear Corrections to NES rates
    A_Out_3  , & ! --- Linear Corrections to NES rates
    A_Pro_1  , & ! --- Linear Corrections to Pair rates
    A_Pro_2  , & ! --- Linear Corrections to Pair rates
    A_Pro_3  , & ! --- Linear Corrections to Pair rates
    A_Ann_1  , & ! --- Linear Corrections to Pair rates
    A_Ann_2  , & ! --- Linear Corrections to Pair rates
    A_Ann_3      ! --- Linear Corrections to Pair rates
  REAL(DP), DIMENSION(nPointsE,nPointsE,nPointsX) :: &
    H_I_0, H_II_0, &  ! --- NES  Scattering Functions (0th moment)
    H_I_1, H_II_1, &  ! --- NES  Scattering Functions (1st moment)
    J_I_0, J_II_0, &  ! --- Pair Scattering Functions (0th moment)
    J_I_1, J_II_1, &  ! --- Pair Scattering Functions (1st moment)
    Nu_J_I_0, Nu_J_II_0, &  ! --- electron neutrino Pair annihilation Scattering Functions (0th moment)
    S_sigma       ! --- Brem Scattering Kernel

  REAL(dp) :: loctot

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'NeutrinoOpacities', &
           nX_Option &
             = [ nX1, 1, 1 ], &
           swX_Option &
             = [ 01, 00, 00 ], &
           bcX_Option &
             = [ 32, 00, 00 ], &
           nE_Option &
             = nE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           ZoomE_Option &
             = ZoomE, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'NeutrinoOpacities'
  WRITE(*,*)
  WRITE(*,'(A6,A,I8.8)') '', 'nPointsX = ', nPointsX
  WRITE(*,'(A6,A,I8.8)') '', 'nPointsE = ', nPointsE
  WRITE(*,*)

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  ! --- Thermodynamic State ---

!  D = 5.6d13 * Unit_D
!  T = 3.0d11 * Unit_T
!  Y = 2.9d-1 * Unit_Y
!  D = 1.0d10 * Unit_D
!  T = 10444070263.062922d0 * Unit_T !3.481d10 * Unit_T
!  Y = 0.46d0 * Unit_Y
!  D = 10.0d0**10.2000000001d0 * Unit_D !2.0d10 * Unit_D
!  T = 0.861733d0 / 8.61733d-11 * Unit_T !0.9d0 / 8.61733d-11 * Unit_T
!  Y = 0.46d0 * Unit_Y

   !D = 2.0d10 * Unit_D
   !T = 0.9d0 / 8.61733d-11 * Unit_T
   !Y = 0.455 * Unit_Y
   !D = 2.0d10 * Unit_D
   !T = 0.95d0 / 8.61733d-11 * Unit_T
   !Y = 0.465 * Unit_Y
   !D = 1.0d12 * Unit_D
   !T = 2.0d0 / 8.61733d-11 * Unit_T
   !Y = 0.33 * Unit_Y

   !D = 0.70428d10 * Unit_D
   D = 0.70428d13 * Unit_D
   T = 0.66006d0 / 8.61733d-11 * Unit_T
   Y = 0.43642d0

  ! --- Energy Grid ---

  DO iN_E = 1, nPointsE
    iE       = MOD( (iN_E-1) / nNodes, nE     ) + 1
    iNodeE   = MOD( (iN_E-1)         , nNodes ) + 1
    dE(iE)   = MeshE % Width(iE)
    E(iN_E)  = NodeCoordinate( MeshE, iE, iNodeE )
    W2(iN_E) = FourPi * WeightsE(iNodeE) * E(iN_E)**2 * dE(iE)
    WRITE(*,'(A6,A2,I3.3,A10,ES8.2E2)') &
      '', 'E(',iN_E,') [MeV] = ', E(iN_E) / Unit_E
  END DO

  ! --- Initialize Equation of State ---

  Timer_ReadEos = 0.0d0
  CALL TimersStart( Timer_ReadEos )
  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option &
             = 'wl-EOS-SFHo-15-25-50.h5', &
           Verbose_Option = .TRUE. )
  CALL TimersStop( Timer_ReadEos )

  ! --- Initialize Opacities ---
 
  Timer_ReadOpacities = 0.0d0
  CALL TimersStart( Timer_ReadOpacities )
  !CALL InitializeOpacities_TABLE &
  !       ( OpacityTableName_EmAb_Option &
  !           = 'wl-Op-LS220-25-50-100-E40-B85-EmAb.h5', &
  !         OpacityTableName_Iso_Option  &
  !           = 'wl-Op-LS220-25-50-100-E40-B85-Iso.h5',  &
  !         OpacityTableName_NES_Option &
  !           = 'wl-Op-LS220-25-50-100-E40-B85-NES.h5',  &
  !         OpacityTableName_Pair_Option &
  !           = 'wl-Op-LS220-25-50-100-E40-B85-Pair.h5', &
  !         OpacityTableName_Brem_Option &
  !           = 'wl-Op-LS220-25-50-100-E40-HR98-Brem.h5', &
  !         Verbose_Option = .TRUE. ) 
  CALL InitializeOpacities_TABLE &
         ( EquationOfStateTableName_Option &
             = 'wl-EOS-SFHo-15-25-50.h5', &
           OpacityTableName_EmAb_Option &
             = 'wl-Op-SFHo-15-25-50-E40-EmAb.h5', &
           OpacityTableName_Iso_Option  &
             = 'wl-Op-SFHo-15-25-50-E40-Iso.h5',  &
           OpacityTableName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-NES.h5',  &
           OpacityTableName_Pair_Option &
             = 'wl-Op-SFHo-15-25-50-E40-Pair.h5', &
           OpacityTableName_Brem_Option &
             = 'wl-Op-SFHo-15-25-50-E40-Brem.h5', &
           Verbose_Option = .TRUE. ) 
  CALL TimersStop( Timer_ReadOpacities )

  ! --- Initialize distributions to zero ---

  f0 = 0.0d0
  f0_DG = 0.0d0
  J = 0.0d0
  H_1 = 0.0d0
  H_2 = 0.0d0
  H_3 = 0.0d0

  ! --- Initialize scattering functions to zero ---

  H_I_0 = 0.0d0
  H_II_0 = 0.0d0
  H_I_1 = 0.0d0
  H_II_1 = 0.0d0

  J_I_0 = 0.0d0
  J_II_0 = 0.0d0
  J_I_1 = 0.0d0
  J_II_1 = 0.0d0

  S_sigma = 0.0d0

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: E, D, T, Y, W2, f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$OMP          H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1 ) &
  !$OMP MAP( alloc: Chi_EmAb, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, Chi_Iso, Sigma_Iso, &
  !$OMP             Eta_EmAb, Eta_NES, Eta_Pair, Eta_NuPair, Eta_Brem, Eta_Iso, &
  !$OMP             A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$OMP             A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( E, D, T, Y, W2, f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$ACC         H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1 ) &
  !$ACC CREATE( Chi_EmAb, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, Chi_Iso, Sigma_Iso, &
  !$ACC         Eta_EmAb, Eta_NES, Eta_Pair, Eta_NuPair, Eta_Brem, Eta_Iso, &
  !$ACC         A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$ACC         A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 )
#endif

  ! --- Compute Equilibrium Distributions ---

  Timer_ComputeEquilibrium = 0.0d0
  CALL TimersStart( Timer_ComputeEquilibrium )
  CALL ComputeEquilibriumDistributions &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Y, f0 )
  CALL TimersStop( Timer_ComputeEquilibrium )

  ! --- Compute Equilibrium Distributions (DG) ---

  Timer_ComputeEquilibrium_DG = 0.0d0
  CALL TimersStart( Timer_ComputeEquilibrium_DG )
  CALL ComputeEquilibriumDistributions_DG &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Y, f0_DG )
  CALL TimersStop( Timer_ComputeEquilibrium_DG )

  ! --- Compute Neutrino Number Density ---

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
  !$ACC PRESENT( J, H_1, H_2, H_3 )
#endif
  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    !J(iE,iS,iX) = 0.0d0
    J(iE,iS,iX) = f0_DG(iE,iS,iX)
    Zero_J(iE,iS,iX) = 0.0d0

    !H_1(iE,iS,iX) = 0.0d0
    !H_2(iE,iS,iX) = 0.0d0
    !H_3(iE,iS,iX) = 0.0d0
    H_1(iE,iS,iX) = 0.1d0
    H_2(iE,iS,iX) = 0.1d0
    H_3(iE,iS,iX) = 0.1d0

  END DO
  END DO
  END DO
  
  ! --- Compute Electron Capture Opacities ---

  Timer_Compute_EC = 0.0d0
  CALL TimersStart( Timer_Compute_EC )
  CALL ComputeNeutrinoOpacities_EC &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Y, f0_DG, Chi_EmAb )
  CALL TimersStop( Timer_Compute_EC )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
  !$ACC PRESENT( Eta_EmAb, Chi_EmAb, f0_DG ) 
#endif
  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    Eta_EmAb(iE,iS,iX) = Chi_EmAb(iE,iS,iX) * f0_DG(iE,iS,iX)

  END DO
  END DO
  END DO

  ! --- Compute Elastic Scattering Opacities ---

  Timer_Compute_ES = 0.0d0
  CALL TimersStart( Timer_Compute_ES )
  CALL ComputeNeutrinoOpacities_ES &
         ( 1, nPointsE, 1, nPointsX, E, D, T, Y, 1, Sigma_Iso )
  CALL TimersStop( Timer_Compute_ES )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
  !$OMP PRIVATE( iC )
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
  !$ACC PRESENT( Eta_Iso, Chi_Iso, Sigma_Iso, f0_DG, E  ) &
  !$ACC PRIVATE( iC )
#endif
  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    iC = nChirals - (LeptonNumber(iS) + 1) / nChirals

    Chi_Iso(iE,iS,iX) = FourPi * E(iE)**2 * Sigma_Iso(iE,iC,iX)
    Eta_Iso(iE,iS,iX) = Chi_Iso(iE,iS,iX) * f0_DG(iE,iS,iX)

  END DO
  END DO
  END DO

  ! --- Compute NES Opacities ---

  Timer_Compute_NES = 0.0d0
  CALL TimersStart( Timer_Compute_NES )

  Timer_ComputeKrnl_NES = 0.0d0
  CALL TimersStart( Timer_ComputeKrnl_NES )
  CALL ComputeNeutrinoOpacities_NES &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, H_I_0, H_II_0 )
  CALL TimersStop( Timer_ComputeKrnl_NES )

  Timer_ComputeRate_NES = 0.0d0
  CALL TimersStart( Timer_ComputeRate_NES )
  CALL ComputeNeutrinoOpacityRates_NES &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, H_I_0, H_II_0, Eta_NES, Chi_NES )
  CALL TimersStop( Timer_ComputeRate_NES )

  ! --- Compute NES Linear Corrections ---

  Timer_ComputeCorr_NES = 0.0d0
  CALL TimersStart( Timer_ComputeCorr_NES )
  CALL ComputeNeutrinoOpacities_NES &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 2, H_I_1, H_II_1 )
  CALL TimersStop( Timer_ComputeCorr_NES )

  Timer_ComputeCorr_NES = 0.0d0
  CALL TimersStart( Timer_ComputeCorr_NES )
  CALL ComputeNeutrinoOpacityRates_LinearCorrections_NES &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, H_1, H_2, H_3, f0_DG, &
           H_I_1, H_II_1, A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3 )
  CALL TimersStop( Timer_ComputeCorr_NES )

  CALL TimersStop( Timer_Compute_NES )

  ! --- Compute Pair Opacities ---

  Timer_Compute_Pair = 0.0d0
  CALL TimersStart( Timer_Compute_Pair )

  Timer_ComputeKrnl_Pair = 0.0d0
  CALL TimersStart( Timer_ComputeKrnl_Pair )
  CALL ComputeNeutrinoOpacities_Pair &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, J_I_0, J_II_0 )
  CALL TimersStop( Timer_ComputeKrnl_Pair )

  Timer_ComputeRate_Pair = 0.0d0
  CALL TimersStart( Timer_ComputeRate_Pair )
  CALL ComputeNeutrinoOpacityRates_Pair &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, J_I_0, J_II_0, Eta_Pair, Chi_Pair )
  CALL TimersStop( Timer_ComputeRate_Pair )

  ! --- Compute Pair Linear Corrections ---

  Timer_ComputeCorr_Pair = 0.0d0
  CALL TimersStart( Timer_ComputeCorr_Pair )
  CALL ComputeNeutrinoOpacities_Pair &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 2, J_I_1, J_II_1 )
  CALL TimersStop( Timer_ComputeCorr_Pair )

  Timer_ComputeCorr_Pair = 0.0d0
  CALL TimersStart( Timer_ComputeCorr_Pair )
  CALL ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, H_1, H_2, H_3, f0_DG, &
           J_I_1, J_II_1, A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 )
  CALL TimersStop( Timer_ComputeCorr_Pair )

  CALL TimersStop( Timer_Compute_Pair )

  ! --- Compute NuPair Opacities ---

  Timer_Compute_NuPair = 0.0d0
  CALL TimersStart( Timer_Compute_NuPair )

  Timer_ComputeKrnl_NuPair = 0.0d0
  CALL TimersStart( Timer_ComputeKrnl_NuPair )
  CALL ComputeNeutrinoOpacities_NuPair &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, Nu_J_I_0, Nu_J_II_0 )
  CALL TimersStop( Timer_ComputeKrnl_NuPair )

  Timer_ComputeRate_NuPair = 0.0d0
  CALL TimersStart( Timer_ComputeRate_NuPair )
  CALL ComputeNeutrinoOpacityRates_NuPair &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, Nu_J_I_0, Nu_J_II_0, Eta_NuPair, Chi_NuPair )
  CALL TimersStop( Timer_ComputeRate_NuPair )

  CALL TimersStop( Timer_Compute_NuPair )

  ! --- Compute Brem Opacities ---
  
  Timer_Compute_Brem = 0.0d0
  CALL TimersStart( Timer_Compute_Brem )

  Timer_ComputeKrnl_Brem = 0.0d0
  CALL TimersStart( Timer_ComputeKrnl_Brem )
  CALL ComputeNeutrinoOpacities_Brem &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, S_sigma )
  CALL TimersStop( Timer_ComputeKrnl_Brem )

  Timer_ComputeRate_Brem = 0.0d0
  CALL TimersStart( Timer_ComputeRate_Brem )
  CALL ComputeNeutrinoOpacityRates_Brem &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, S_sigma, Eta_Brem, Chi_Brem )
  CALL TimersStop( Timer_ComputeRate_Brem )

  CALL TimersStop( Timer_Compute_Brem )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from: f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$OMP            Chi_EmAb, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, Chi_Iso, &
  !$OMP            Eta_EmAb, Eta_NES, Eta_Pair, Eta_NuPair, Eta_Brem, Eta_Iso, &
  !$OMP            H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1, &
  !$OMP            A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$OMP            A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 ) &
  !$OMP MAP( release: E, D, T, Y, W2, Sigma_Iso )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC COPYOUT( f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$ACC          Chi_EmAb, Chi_NES, Chi_Pair, Chi_NuPair, Chi_Brem, Chi_Iso, &
  !$ACC          Eta_EmAb, Eta_NES, Eta_Pair, Eta_NuPair, Eta_Brem, Eta_Iso, &
  !$ACC          H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1, &
  !$ACC          A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$ACC          A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 ) &
  !$ACC DELETE( E, D, T, Y, W2, Sigma_Iso )
#endif

  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    Edot_EmAb(iE,iS,iX)   = FourPi / PlanckConstant**3 / D(iX) * abs(Chi_EmAb(iE,iS,iX)) * E(iE)**3
    !Edot_EmAb(iE,iS,iX) = abs(Chi_EmAb(iE,iS,iX))
    Edot_Iso(iE,iS,iX)    = FourPi / PlanckConstant**3 / D(iX) * abs(Chi_Iso(iE,iS,iX))  * E(iE)**3
    Qdot_Pair(iE,iS,iX)   = FourPi / PlanckConstant**3 / D(iX) * Eta_Pair(iE,iS,iX)      * E(iE)**3
    Qdot_NuPair(iE,iS,iX) = FourPi / PlanckConstant**3 / D(iX) * Eta_NuPair(iE,iS,iX)    * E(iE)**3
    Qdot_Brem(iE,iS,iX)   = FourPi / PlanckConstant**3 / D(iX) * Eta_Brem(iE,iS,iX)      * E(iE)**3
    Edot_NES(iE,iS,iX)    = FourPi / PlanckConstant**3 / D(iX) * abs(Chi_NES(iE,iS,iX))  * E(iE)**3

  END DO
  END DO
  END DO

  Edot_EmAb_element   = 0.0d0
  Edot_Iso_element    = 0.0d0
  Qdot_Pair_element   = 0.0d0
  Qdot_NuPair_element = 0.0d0
  Qdot_Brem_element   = 0.0d0
  Edot_NES_element    = 0.0d0

  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iN_E = 1, nPointsE

    iE       = MOD( (iN_E-1) / nNodes, nE     ) + 1
    iNodeE   = MOD( (iN_E-1)         , nNodes ) + 1

    Edot_EmAb_element(iE,iS,iX)   = Edot_EmAb_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Edot_EmAb(iN_E,iS,iX) &
                                  * (E(iN_E) / Unit_E)**2
    Edot_Iso_element(iE,iS,iX)    = Edot_Iso_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Edot_Iso(iN_E,iS,iX) &
                                  * (E(iN_E) / Unit_E)**2
    Qdot_Pair_element(iE,iS,iX)   = Qdot_Pair_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Qdot_Pair(iN_E,iS,iX) &
                                  * (E(iN_E) / Unit_E)**2
    Qdot_NuPair_element(iE,iS,iX) = Qdot_NuPair_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Qdot_NuPair(iN_E,iS,iX) &
                                  * (E(iN_E) / Unit_E)**2
    Qdot_Brem_element(iE,iS,iX)   = Qdot_Brem_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Qdot_Brem(iN_E,iS,iX) &
                                  * (E(iN_E) / Unit_E)**2
    Edot_NES_element(iE,iS,iX)    = Edot_NES_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Edot_NES(iN_E,iS,iX) &
                                  * (E(iN_E) / Unit_E)**2
  END DO
  END DO
  END DO

  CALL WriteVector &
         ( nPointsE, E / Unit_E, 'E.dat' )

  CALL WriteVector &
         ( nE, MeshE % Center / Unit_E, 'EC.dat' )
  CALL WriteVector & ! --- NuE
         ( nE, Edot_EmAb_element(:,iNuE, 1) / Unit_Edot, 'Edot_EmAb_NuE_element.dat' )
  CALL WriteVector & ! --- NuE
         ( nE, Edot_Iso_element (:,iNuE, 1) / Unit_Edot, 'Edot_Iso_NuE_element.dat'  )
  CALL WriteVector & ! --- NuE
         ( nE, Qdot_Pair_element(:,iNuE, 1) / Unit_Edot, 'Qdot_Pair_NuE_element.dat' )
  CALL WriteVector & ! --- NuE
         ( nE, Qdot_NuPair_element(:,iNuE, 1) / Unit_Edot, 'Qdot_NuPair_NuE_element.dat' )
  CALL WriteVector & ! --- NuM
         ( nE, Qdot_NuPair_element(:,iNuM, 1) / Unit_Edot, 'Qdot_NuPair_NuM_element.dat' )
  CALL WriteVector & ! --- NuE
         ( nE, Qdot_Brem_element(:,iNuE, 1) / Unit_Edot, 'Qdot_Brem_NuE_element.dat' )
  CALL WriteVector & ! --- NuE
         ( nE, Edot_NES_element (:,iNuE, 1) / Unit_Edot, 'Edot_NES_NuE_element.dat'  )

  CALL WriteVector & ! --- NuE
         ( nPointsE, f0   (:,iNuE    ,1), 'f0_NuE.dat'        )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, f0   (:,iNuE_Bar,1), 'f0_NuE_Bar.dat'    )
  CALL WriteVector & ! --- NuE
         ( nPointsE, f0_DG(:,iNuE    ,1), 'f0_DG_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, f0_DG(:,iNuE_Bar,1), 'f0_DG_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_EmAb(:,iNuE    ,1) / Unit_Chi, 'Chi_EmAb_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_EmAb(:,iNuE    ,1) / Unit_Chi, 'Eta_EmAb_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_EmAb(:,iNuE_Bar,1) / Unit_Chi, 'Chi_EmAb_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_EmAb(:,iNuE_Bar,1) / Unit_Chi, 'Eta_EmAb_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_Iso(:,iNuE    ,1) / Unit_Chi, 'Chi_Iso_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_Iso(:,iNuE    ,1) / Unit_Chi, 'Eta_Iso_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_Iso(:,iNuE_Bar,1) / Unit_Chi, 'Chi_Iso_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_Iso(:,iNuE_Bar,1) / Unit_Chi, 'Eta_Iso_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_NES(:,iNuE    ,1) / Unit_Chi, 'Chi_NES_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_NES(:,iNuE    ,1) / Unit_Chi, 'Eta_NES_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_NES(:,iNuE_Bar,1) / Unit_Chi, 'Chi_NES_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_NES(:,iNuE_Bar,1) / Unit_Chi, 'Eta_NES_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, A_In_1 (:,iNuE    ,1) / Unit_Chi, 'A_In_NuE.dat'        )
  CALL WriteVector & ! --- NuE
         ( nPointsE, A_Out_1(:,iNuE    ,1) / Unit_Chi, 'A_Out_NuE.dat'       )
  CALL WriteVector & ! --- NuE
         ( nPointsE, A_In_1 (:,iNuE_Bar,1) / Unit_Chi, 'A_In_NuE_Bar.dat'    )
  CALL WriteVector & ! --- NuE
         ( nPointsE, A_Out_1(:,iNuE_Bar,1) / Unit_Chi, 'A_Out_NuE_Bar.dat'   )

  CALL WriteVector & ! --- NuE
         ( nPointsE, A_Pro_1 (:,iNuE   ,1) / Unit_Chi, 'A_Pro_NuE.dat'       )
  CALL WriteVector & ! --- NuE
         ( nPointsE, A_Ann_1(:,iNuE    ,1) / Unit_Chi, 'A_Ann_NuE.dat'       )
  CALL WriteVector & ! --- NuE
         ( nPointsE, A_Pro_1(:,iNuE_Bar,1) / Unit_Chi, 'A_Pro_NuE_Bar.dat'   )
  CALL WriteVector & ! --- NuE
         ( nPointsE, A_Ann_1(:,iNuE_Bar,1) / Unit_Chi, 'A_Ann_NuE_Bar.dat'   )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_Pair(:,iNuE    ,1) / Unit_Chi, 'Chi_Pair_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_Pair(:,iNuE    ,1) / Unit_Chi, 'Eta_Pair_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_Pair(:,iNuE_Bar,1) / Unit_Chi, 'Chi_Pair_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_Pair(:,iNuE_Bar,1) / Unit_Chi, 'Eta_Pair_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_NuPair(:,iNuE    ,1) / Unit_Chi, 'Chi_NuPair_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_NuPair(:,iNuE    ,1) / Unit_Chi, 'Eta_NuPair_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_NuPair(:,iNuE_Bar,1) / Unit_Chi, 'Chi_NuPair_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_NuPair(:,iNuE_Bar,1) / Unit_Chi, 'Eta_NuPair_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuM
         ( nPointsE, Chi_NuPair(:,iNuM    ,1) / Unit_Chi, 'Chi_NuPair_NuM.dat'     )
  CALL WriteVector & ! --- NuM
         ( nPointsE, Eta_NuPair(:,iNuM    ,1) / Unit_Chi, 'Eta_NuPair_NuM.dat'     )
  CALL WriteVector & ! --- NuM_Bar
         ( nPointsE, Chi_NuPair(:,iNuM_Bar,1) / Unit_Chi, 'Chi_NuPair_NuM_Bar.dat' )
  CALL WriteVector & ! --- NuM_Bar
         ( nPointsE, Eta_NuPair(:,iNuM_Bar,1) / Unit_Chi, 'Eta_NuPair_NuM_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_Brem(:,iNuE    ,1) / Unit_Chi, 'Chi_Brem_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_Brem(:,iNuE    ,1) / Unit_Chi, 'Eta_Brem_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_Brem(:,iNuE_Bar,1) / Unit_Chi, 'Chi_Brem_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_Brem(:,iNuE_Bar,1) / Unit_Chi, 'Eta_Brem_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Edot_EmAb(:,iNuE    ,1) / Unit_Edot, 'Edot_EmAb_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Edot_EmAb(:,iNuE_Bar,1) / Unit_Edot, 'Edot_EmAb_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Qdot_Pair(:,iNuE    ,1) / Unit_Qdot, 'Qdot_Pair_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Qdot_Pair(:,iNuE_Bar,1) / Unit_Qdot, 'Qdot_Pair_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Qdot_NuPair(:,iNuE    ,1) / Unit_Qdot, 'Qdot_NuPair_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Qdot_NuPair(:,iNuE_Bar,1) / Unit_Qdot, 'Qdot_NuPair_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuM
         ( nPointsE, Qdot_NuPair(:,iNuM    ,1) / Unit_Qdot, 'Qdot_NuPair_NuM.dat'     )
  CALL WriteVector & ! --- NuM_Bar
         ( nPointsE, Qdot_NuPair(:,iNuM_Bar,1) / Unit_Qdot, 'Qdot_NuPair_NuM_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Qdot_Brem(:,iNuE    ,1) / Unit_Qdot, 'Qdot_Brem_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Qdot_Brem(:,iNuE_Bar,1) / Unit_Qdot, 'Qdot_Brem_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Edot_NES(:,iNuE    ,1)  / Unit_Qdot, 'Edot_NES_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Edot_NES(:,iNuE_Bar,1)  / Unit_Qdot, 'Edot_NES_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Edot_Iso(:,iNuE    ,1)  / Unit_Qdot, 'Edot_Iso_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Edot_Iso(:,iNuE_Bar,1)  / Unit_Qdot, 'Edot_Iso_NuE_Bar.dat' )

  CALL WriteMatrix &
         ( nPointsE, nPointsE, H_I_0(:,:,1), 'H1.dat'  )
  CALL WriteMatrix &
         ( nPointsE, nPointsE, H_I_1(:,:,1), 'H1_1.dat'  )

  CALL WriteMatrix &
         ( nPointsE, nPointsE, J_I_0(:,:,1), 'J1.dat' )
  CALL WriteMatrix &
         ( nPointsE, nPointsE, J_I_1(:,:,1), 'J1_1.dat' )

  CALL WriteMatrix &
         ( nPointsE, nPointsE, S_sigma(:,:,1), 'S_sigma.dat' )

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  Timer_Total &
    = Timer_Compute_EC + Timer_Compute_ES &
      + Timer_Compute_NES + Timer_Compute_Pair &
      + Timer_Compute_NuPair + Timer_Compute_Brem

  WRITE(*,*)
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ReadEos = ',       &
    Timer_ReadEos
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ReadOpacities = ', &
    Timer_ReadOpacities
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ComputeEquil = ', &
    Timer_ComputeEquilibrium
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ComputeEquil_DG = ', &
    Timer_ComputeEquilibrium_DG
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_EC = ',    &
    Timer_Compute_EC, Timer_Compute_EC / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_ES = ',    &
    Timer_Compute_ES, Timer_Compute_ES / Timer_Total

  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_NES = ',   &
    Timer_Compute_NES, Timer_Compute_NES / Timer_Total
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Kernel)',   &
    Timer_ComputeKrnl_NES
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Rate)',   &
    Timer_ComputeRate_NES
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Corrections)',   &
    Timer_ComputeCorr_NES

  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_Pair = ',  &
    Timer_Compute_Pair, Timer_Compute_Pair / Timer_Total
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Kernel)',   &
    Timer_ComputeKrnl_Pair
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Rate)',   &
    Timer_ComputeRate_Pair
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Corrections)',   &
    Timer_ComputeCorr_Pair

  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_NuPair = ',  &
    Timer_Compute_NuPair, Timer_Compute_NuPair / Timer_Total
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Kernel)',   &
    Timer_ComputeKrnl_NuPair
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Rate)',   &
    Timer_ComputeRate_NuPair

  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_Brem = ',  &
    Timer_Compute_Brem, Timer_Compute_Brem / Timer_Total
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Kernel)',   &
    Timer_ComputeKrnl_Brem
  WRITE(*,'(A4,A22,1ES10.2E2)') '', '(Rate)',   &
    Timer_ComputeRate_Brem
  WRITE(*,*)

  CALL FinalizeDevice

  CALL MPI_FINALIZE( mpierr )


CONTAINS


  PURE REAL(DP) FUNCTION TRAPEZ( n, x, y )

    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(in) :: x(n), y(n)

    INTEGER :: i

    TRAPEZ = 0.0_DP
    DO i = 1, n - 1
      TRAPEZ = TRAPEZ + 0.5_dp * ( x(i+1) - x(i) ) * ( y(i) + y(i+1) )
    END DO

    RETURN
  END FUNCTION TRAPEZ


END PROGRAM NeutrinoOpacities
