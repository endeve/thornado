PROGRAM NeutrinoOpacities

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
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
    C1, C2
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar
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
    ComputeNeutrinoOpacities_Brem, &
    ComputeNeutrinoOpacityRates_Brem
  USE DeviceModule, ONLY: &
    InitializeDevice, &
    FinalizeDevice
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: &
    nNodes   = 2, &
    nE       = 2**4, &
    nX1      = 2**11, &
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
    eL         = 0.0e0_DP * Unit_E, &
    eR         = 3.0e2_DP * Unit_E, &
    ZoomE      = 1.183081754893913_DP

  INTEGER :: &
    mpierr, iE, iX, iS, iNodeE, iN_E, iE1, iE2
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
    Timer_ComputeKrnl_Brem, &
    Timer_ComputeRate_Brem, &
    Timer_Compute_Brem, &
    Timer_Total
  REAL(DP), DIMENSION(nPointsX) :: &
    D, T, Y
  REAL(DP), DIMENSION(nE) :: &
    dE
  REAL(DP), DIMENSION(nPointsE) :: &
    E, W2
  REAL(DP), DIMENSION(nPointsE,nPointsX) :: &
    Sigma_Iso    ! --- Iso-energertic Kernel
  REAL(DP), DIMENSION(nPointsE,nSpecies,nPointsX) :: &
    f0       , & ! --- Equilibrium Distribution
    f0_DG    , & ! --- Equilibrium Distribution (DG Approximation)
    J        , & ! --- Neutrino Number Density
    H_1      , & ! --- Neutrino Flux (x-direction)
    H_2      , & ! --- Neutrino Flux (y-direction)
    H_3      , & ! --- Neutrino Flux (z-direction)
    Eta_EmAb , & ! --- Electron Capture Emissivity
    Chi_EmAb , & ! --- Electron Capture Opacity
    Eta_Iso  , & ! --- Iso Emissivity
    Chi_Iso  , & ! --- Iso Opacity
    Eta_NES  , & ! --- NES Emissivity
    Chi_NES  , & ! --- NES Opacity
    Eta_Pair , & ! --- Pair Emissivity
    Chi_Pair , & ! --- Pair Opacity
    Qdot_Pair, & ! --- Pair Heating Rate
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
    S_sigma       ! --- Brem Scattering Kernel

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

!  D = 1.3d14 * Unit_D
  D = 5.6d13 * Unit_D
  T = 3.0d11 * Unit_T
  Y = 2.9d-1 * Unit_Y
!  D = 1.050d10 * Unit_D
!  T = 3.481d10 * Unit_T
!  Y = 2.530d-1 * Unit_Y

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

  CALL TimersStart( Timer_ReadEos )
  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5', &
           Verbose_Option = .TRUE. )
  CALL TimersStop( Timer_ReadEos )

  ! --- Initialize Opacities ---

  CALL TimersStart( Timer_ReadOpacities )
  CALL InitializeOpacities_TABLE &
         ( OpacityTableName_EmAb_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &
           OpacityTableName_Iso_Option  &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5',  &
           OpacityTableName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5',  &
           OpacityTableName_Pair_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &
           OpacityTableName_Brem_Option &
             = 'wl-Op-SFHo-15-25-50-E40-HR98-Brem.h5', &
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
  !$OMP MAP( alloc: Chi_EmAb, Chi_NES, Chi_Pair, Chi_Brem, Chi_Iso, Sigma_Iso, &
  !$OMP             Eta_EmAb, Eta_NES, Eta_Pair, Eta_Brem, Eta_Iso, &
  !$OMP             A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$OMP             A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( E, D, T, Y, W2, f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$ACC         H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1 ) &
  !$ACC CREATE( Chi_EmAb, Chi_NES, Chi_Pair, Chi_Brem, Chi_Iso, Sigma_Iso, &
  !$ACC         Eta_EmAb, Eta_NES, Eta_Pair, Eta_Brem, Eta_Iso, &
  !$ACC         A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$ACC         A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 )
#endif

  ! --- Compute Equilibrium Distributions ---

  CALL TimersStart( Timer_ComputeEquilibrium )
  CALL ComputeEquilibriumDistributions &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Y, f0 )
  CALL TimersStop( Timer_ComputeEquilibrium )

  ! --- Compute Equilibrium Distributions (DG) ---

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

    J(iE,iS,iX) = 0.0d0
    !J(iE,iS,iX) = f0_DG(iE,iS,iX)

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

  CALL TimersStart( Timer_Compute_EC )
  CALL ComputeNeutrinoOpacities_EC &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Y, Chi_EmAb )
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

  CALL TimersStart( Timer_Compute_ES )
  CALL ComputeNeutrinoOpacities_ES &
         ( 1, nPointsE, 1, nPointsX, E, D, T, Y, 1, Sigma_Iso )
  CALL TimersStop( Timer_Compute_ES )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
  !$ACC PRESENT( Eta_Iso, Chi_Iso, Sigma_Iso, f0_DG, E  )
#endif
  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    Chi_Iso(iE,iS,iX) = FourPi * E(iE)**2 * Sigma_Iso(iE,iX)
    Eta_Iso(iE,iS,iX) = Chi_Iso(iE,iS,iX) * f0_DG(iE,iS,iX)

  END DO
  END DO
  END DO

  ! --- Compute NES Opacities ---

  CALL TimersStart( Timer_Compute_NES )

  CALL TimersStart( Timer_ComputeKrnl_NES )
  CALL ComputeNeutrinoOpacities_NES &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, H_I_0, H_II_0 )
  CALL TimersStop( Timer_ComputeKrnl_NES )

  CALL TimersStart( Timer_ComputeRate_NES )
  CALL ComputeNeutrinoOpacityRates_NES &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, H_I_0, H_II_0, Eta_NES, Chi_NES )
  CALL TimersStop( Timer_ComputeRate_NES )

  ! --- Compute NES Linear Corrections ---

  CALL TimersStart( Timer_ComputeCorr_NES )
  CALL ComputeNeutrinoOpacities_NES &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 2, H_I_1, H_II_1 )
  CALL TimersStop( Timer_ComputeCorr_NES )

  CALL TimersStart( Timer_ComputeCorr_NES )
  CALL ComputeNeutrinoOpacityRates_LinearCorrections_NES &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, H_1, H_2, H_3, f0_DG, &
           H_I_1, H_II_1, A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3 )
  CALL TimersStop( Timer_ComputeCorr_NES )

  CALL TimersStop( Timer_Compute_NES )

  ! --- Compute Pair Opacities ---

  CALL TimersStart( Timer_Compute_Pair )

  CALL TimersStart( Timer_ComputeKrnl_Pair )
  CALL ComputeNeutrinoOpacities_Pair &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, J_I_0, J_II_0 )
  CALL TimersStop( Timer_ComputeKrnl_Pair )

  CALL TimersStart( Timer_ComputeRate_Pair )
  CALL ComputeNeutrinoOpacityRates_Pair &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, J_I_0, J_II_0, Eta_Pair, Chi_Pair )
  CALL TimersStop( Timer_ComputeRate_Pair )

  ! --- Compute Pair Linear Corrections ---

  CALL TimersStart( Timer_ComputeCorr_Pair )
  CALL ComputeNeutrinoOpacities_Pair &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 2, J_I_1, J_II_1 )
  CALL TimersStop( Timer_ComputeCorr_Pair )

  CALL TimersStart( Timer_ComputeCorr_Pair )
  CALL ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, H_1, H_2, H_3, f0_DG, &
           J_I_1, J_II_1, A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 )
  CALL TimersStop( Timer_ComputeCorr_Pair )

  CALL TimersStop( Timer_Compute_Pair )

  ! --- Compute Brem Opacities ---

  CALL TimersStart( Timer_Compute_Brem )

  CALL TimersStart( Timer_ComputeKrnl_Brem )
  CALL ComputeNeutrinoOpacities_Brem &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, S_sigma )
  CALL TimersStop( Timer_ComputeKrnl_Brem )

  CALL TimersStart( Timer_ComputeRate_Brem )
  CALL ComputeNeutrinoOpacityRates_Brem &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, D, W2, &
           J, f0_DG, S_sigma, Eta_Brem, Chi_Brem )
  CALL TimersStop( Timer_ComputeRate_Brem )

  CALL TimersStop( Timer_Compute_Brem )

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from: f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$OMP            Chi_EmAb, Chi_NES, Chi_Pair, Chi_Brem, Chi_Iso, &
  !$OMP            Eta_EmAb, Eta_NES, Eta_Pair, Eta_Brem, Eta_Iso, &
  !$OMP            H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1, &
  !$OMP            A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$OMP            A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 ) &
  !$OMP MAP( release: E, D, T, Y, W2, Sigma_Iso )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC COPYOUT( f0, f0_DG, J, H_1, H_2, H_3, S_sigma, &
  !$ACC          Chi_EmAb, Chi_NES, Chi_Pair, Chi_Brem, Chi_Iso, &
  !$ACC          Eta_EmAb, Eta_NES, Eta_Pair, Eta_Brem, Eta_Iso, &
  !$ACC          H_I_0, H_II_0, H_I_1, H_II_1, J_I_0, J_II_0, J_I_1, J_II_1, &
  !$ACC          A_In_1, A_In_2, A_In_3, A_Out_1, A_Out_2, A_Out_3, &
  !$ACC          A_Pro_1, A_Pro_2, A_Pro_3, A_Ann_1, A_Ann_2, A_Ann_3 ) &
  !$ACC DELETE( E, D, T, Y, W2, Sigma_Iso )
#endif

  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    Qdot_Pair(iE,iS,iX) = FourPi / PlanckConstant**3 / D(iX) * Eta_Pair(iE,iS,iX) * E(iE)**3
    Qdot_Brem(iE,iS,iX) = FourPi / PlanckConstant**3 / D(iX) * Eta_Brem(iE,iS,iX) * E(iE)**3

  END DO
  END DO
  END DO

  CALL WriteVector &
         ( nPointsE, E / Unit_E, 'E.dat' )

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
         ( nPointsE, Chi_Brem(:,iNuE    ,1) / Unit_Chi, 'Chi_Brem_NuE.dat'     )
  CALL WriteVector & ! --- NuE
         ( nPointsE, Eta_Brem(:,iNuE    ,1) / Unit_Chi, 'Eta_Brem_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_Brem(:,iNuE_Bar,1) / Unit_Chi, 'Chi_Brem_NuE_Bar.dat' )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Eta_Brem(:,iNuE_Bar,1) / Unit_Chi, 'Eta_Brem_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Qdot_Pair(:,iNuE    ,1) / Unit_Qdot, 'Qdot_Pair_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Qdot_Pair(:,iNuE_Bar,1) / Unit_Qdot, 'Qdot_Pair_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Qdot_Brem(:,iNuE    ,1) / Unit_Qdot, 'Qdot_Brem_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Qdot_Brem(:,iNuE_Bar,1) / Unit_Qdot, 'Qdot_Brem_NuE_Bar.dat' )

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
      + Timer_Compute_Brem

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
