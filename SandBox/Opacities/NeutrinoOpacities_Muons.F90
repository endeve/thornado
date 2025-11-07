PROGRAM NeutrinoOpacities_Muons

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV, &
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
    WriteVector
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE, &
    WeightsE
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar, iNuM, iNuM_Bar
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, &
    ComputeEquilibriumDistributions_DG, &
    ComputeNeutrinoOpacities_EC
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop

  USE mpi

  IMPLICIT NONE

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
    Unit_Ye    = 1.0_DP, &
    Unit_Ym    = 1.0_DP, &
    Unit_E     = MeV, &
    Unit_Chi   = 1.0_DP / Centimeter, &
    BaryonMass = AtomicMassUnit, &
    Unit_Edot  = 1.0_DP / ( BaryonMass * Second ), &
    eL         = 0.0e0_DP * Unit_E, &
    eR         = 3.0e2_DP * Unit_E, &
    ZoomE      = 1.266038160710160d0

#ifdef EOSMODE_3D
    CHARACTER(LEN=128) :: FileName_EOS = "3DEOSTable.h5"
#elif defined(EOSMODE_4D)
    CHARACTER(LEN=128) :: FileName_EOS = "4DEOSTable.h5"
#elif defined(EOSMODE_COMPOSE)
    CHARACTER(LEN=128) :: FileName_EOS = "BaryonsPlusPhotonsPlusLeptonsEOS.h5"
#endif

    CHARACTER(LEN=128) :: FileName_EmAb = "wl-Op-SFHo-15-25-50-E40-EmAb.h5"

  INTEGER :: &
    iE, iX, iS, iNodeE, iN_E
  REAL(DP) :: &
    Timer_ReadEos, &
    Timer_ReadOpacities, &
    Timer_ComputeEquilibrium, &
    Timer_ComputeEquilibrium_DG, &
    Timer_Compute_EC, &
    Timer_Total
  REAL(DP), DIMENSION(nPointsX) :: &
    D, T, Ye, Ym
  REAL(DP), DIMENSION(nE) :: &
    dE
  REAL(DP), DIMENSION(nE,nSpecies,nPointsX) :: &
    Edot_EmAb_element
  REAL(DP), DIMENSION(nPointsE) :: &
    E, W2
  REAL(DP), DIMENSION(nPointsE,nSpecies,nPointsX) :: &
    f0       , & ! --- Equilibrium Distribution
    f0_DG    , & ! --- Equilibrium Distribution (DG Approximation)
    J        , & ! --- Neutrino Number Density
    Zero_J   , & ! --- All zeros
    Eta_EmAb , & ! --- Electron Capture Emissivity
    Chi_EmAb , & ! --- Electron Capture Opacity
    Edot_EmAb    ! --- FourPi / PlanckConstant**3 / D * |Chi_EmAb| * E**3

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'NeutrinoOpacities_Muons', &
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
  WRITE(*,'(A4,A)') '', 'NeutrinoOpacities_Muons'
  WRITE(*,*)
  WRITE(*,'(A6,A,I8.8)') '', 'nPointsX = ', nPointsX
  WRITE(*,'(A6,A,I8.8)') '', 'nPointsE = ', nPointsE
  WRITE(*,*)

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  ! --- Thermodynamic State ---

  ! Case (a) of Table II in Fischer 2020
  D  = 1.0E+13_DP * Gram / Centimeter**3
  T  = 10.0_DP * MeV
  Ye = 0.2_DP
  Ym = 1.0E-04_DP

  ! Case (b) of Table II in Fischer 2020
  !D  = 1.0E+14_DP * Gram / Centimeter**3
  !T  = 25.0_DP * MeV
  !Ye = 0.15_DP
  !Ym = 0.05_DP

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
             = FileName_EOS, &
           Verbose_Option = .TRUE. )
  CALL TimersStop( Timer_ReadEos )

  ! --- Initialize Opacities ---
 
  Timer_ReadOpacities = 0.0d0
  CALL TimersStart( Timer_ReadOpacities )
  CALL InitializeOpacities_TABLE &
         ( EquationOfStateTableName_Option &
             = FileName_EOS, &
           OpacityTableName_EmAb_Option &
             = FileName_EmAb, &
           Verbose_Option = .TRUE. ) 
  CALL TimersStop( Timer_ReadOpacities )

  ! --- Initialize distributions to zero ---

  f0 = 0.0d0
  f0_DG = 0.0d0
  J = 0.0d0

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: E, D, T, Ye, Ym, W2, f0, f0_DG, J ) &
  !$OMP MAP( alloc: Chi_EmAb, Eta_EmAb )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( E, D, T, Ye, Ym, W2, f0, f0_DG, J ) &
  !$ACC CREATE( Chi_EmAb, Eta_EmAb )
#endif

  ! --- Compute Equilibrium Distributions ---

  Timer_ComputeEquilibrium = 0.0d0
  CALL TimersStart( Timer_ComputeEquilibrium )
  CALL ComputeEquilibriumDistributions &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Ye, Ym, f0 )
  CALL TimersStop( Timer_ComputeEquilibrium )

  ! --- Compute Equilibrium Distributions (DG) ---

  Timer_ComputeEquilibrium_DG = 0.0d0
  CALL TimersStart( Timer_ComputeEquilibrium_DG )
  CALL ComputeEquilibriumDistributions_DG &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Ye, Ym, f0_DG )
  CALL TimersStop( Timer_ComputeEquilibrium_DG )

  ! --- Compute Neutrino Number Density ---

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
  !$ACC PRESENT( J )
#endif
  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    !J(iE,iS,iX) = 0.0d0
    J(iE,iS,iX) = f0_DG(iE,iS,iX)
    Zero_J(iE,iS,iX) = 0.0d0

  END DO
  END DO
  END DO
  
  ! --- Compute Electron/Muon Capture Opacities ---

  Timer_Compute_EC = 0.0d0
  CALL TimersStart( Timer_Compute_EC )
  CALL ComputeNeutrinoOpacities_EC &
         ( 1, nPointsE, 1, nSpecies, 1, nPointsX, E, D, T, Ye, Ym, f0_DG, Chi_EmAb )
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


#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from: f0, f0_DG, J, &
  !$OMP            Chi_EmAb, Eta_EmAb ) &
  !$OMP MAP( release: E, D, T, Ye, Ym, W2 )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC COPYOUT( f0, f0_DG, J, &
  !$ACC          Chi_EmAb, Eta_EmAb ) &
  !$ACC DELETE( E, D, T, Ye, Ym, W2 )
#endif

  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iE = 1, nPointsE

    Edot_EmAb(iE,iS,iX)   = FourPi / PlanckConstant**3 / D(iX) * abs(Chi_EmAb(iE,iS,iX)) * E(iE)**3
    !Edot_EmAb(iE,iS,iX) = abs(Chi_EmAb(iE,iS,iX))

  END DO
  END DO
  END DO

  Edot_EmAb_element   = 0.0d0

  DO iX = 1, nPointsX
  DO iS = 1, nSpecies
  DO iN_E = 1, nPointsE

    iE       = MOD( (iN_E-1) / nNodes, nE     ) + 1
    iNodeE   = MOD( (iN_E-1)         , nNodes ) + 1

    Edot_EmAb_element(iE,iS,iX)   = Edot_EmAb_element(iE,iS,iX) &
                                  + dE(iE) / Unit_E * WeightsE(iNodeE) * Edot_EmAb(iN_E,iS,iX) &
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

  CALL WriteVector & ! --- NuM
         ( nPointsE, Chi_EmAb(:,iNuM    ,1) / Unit_Chi, 'Chi_EmAb_NuM.dat'     )
  CALL WriteVector & ! --- NuM
         ( nPointsE, Eta_EmAb(:,iNuM    ,1) / Unit_Chi, 'Eta_EmAb_NuM.dat'     )
  CALL WriteVector & ! --- NuM_Bar
         ( nPointsE, Chi_EmAb(:,iNuM_Bar,1) / Unit_Chi, 'Chi_EmAb_NuM_Bar.dat' )
  CALL WriteVector & ! --- NuM_Bar
         ( nPointsE, Eta_EmAb(:,iNuM_Bar,1) / Unit_Chi, 'Eta_EmAb_NuM_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Edot_EmAb(:,iNuE    ,1) / Unit_Edot, 'Edot_EmAb_NuE.dat'     )
  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Edot_EmAb(:,iNuE_Bar,1) / Unit_Edot, 'Edot_EmAb_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuM
         ( nPointsE, Edot_EmAb(:,iNuM    ,1) / Unit_Edot, 'Edot_EmAb_NuM.dat'     )
  CALL WriteVector & ! --- NuM_Bar
         ( nPointsE, Edot_EmAb(:,iNuM_Bar,1) / Unit_Edot, 'Edot_EmAb_NuM_Bar.dat' )

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  Timer_Total &
    = Timer_Compute_EC

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

  CALL FinalizeReferenceElementE
  CALL FinalizeReferenceElementE_Lagrange
  CALL FinalizeProgram


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


END PROGRAM NeutrinoOpacities_Muons