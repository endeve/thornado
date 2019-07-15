#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_NeutrinoMatterSolverModule

  USE KindModule, ONLY: &
    DP, Zero, One, FourPi
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    BoltzmannConstant, &
    AtomicMassUnit, &
    Centimeter, &
    Gram, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
    nNodesZ, &
    nDOFE
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Im_ComputeOpacity, &
    Timer_Im_ComputeRate, &
    Timer_Im_ComputeLS, &
    Timer_Im_UpdateFP
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply, &
    LinearLeastSquares, &
    LinearLeastSquares_LWORK
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_Points, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacities_NES_Point, &
    ComputeNeutrinoOpacities_Pair_Points, &
    ComputeNeutrinoOpacities_Pair_Point, &
    ComputeNeutrinoOpacitiesRates_NES_Points, &
    ComputeNeutrinoOpacitiesRates_Pair_Points

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: SolveMatterEquations_EmAb
  !PUBLIC :: SolveMatterEquations_EmAb_NuE
  PUBLIC :: SolveMatterEquations_FP_Coupled

  ! --- Units Only for Displaying to Screen ---

  PUBLIC :: SolveMatterEquations_EmAb
  !PUBLIC :: SolveMatterEquations_EmAb_NuE
  PUBLIC :: SolveMatterEquations_FP_Coupled

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

  LOGICAL, PARAMETER :: SolveMatter = .TRUE.
  LOGICAL, PARAMETER :: UsePreconditionerEmAb = .FALSE.
  LOGICAL, PARAMETER :: UsePreconditionerPair = .FALSE.
  LOGICAL, PARAMETER :: UsePreconditionerPairLagAllButJ0 = .FALSE.

  INTEGER,  PARAMETER :: M_FP = 3
  REAL(DP), PARAMETER :: WFactor_FP = FourPi / PlanckConstant**3

  INTEGER  :: nE_G, nX_G, n_FP, nZ(4)
  INTEGER  :: iE_B, iE_E

  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:)
  REAL(DP), ALLOCATABLE :: W3_S(:)

  REAL(DP), ALLOCATABLE :: AMAT(:,:,:)
  REAL(DP), ALLOCATABLE :: BVEC(:,:)
  REAL(DP), ALLOCATABLE :: TAU (:,:)
  REAL(DP), ALLOCATABLE :: WORK(:,:)
  INTEGER,  ALLOCATABLE :: INFO(:)  
  INTEGER               :: LWORK

  PUBLIC :: E_N

CONTAINS


  SUBROUTINE InitializeNeutrinoMatterSolver( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    REAL(DP) :: TMP(1)

    iE_B = iZ_B(1)
    iE_E = iZ_E(1)
    nZ   = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )
    n_FP = 2 + 2*nE_G

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S(:) = WFactor_FP * W2_N(:)
    W3_S(:) = WFactor_FP * W3_N(:)

    ALLOCATE( AMAT(n_FP,M_FP,nX_G) )
    ALLOCATE( BVEC(n_FP,nX_G) )
    ALLOCATE( TAU (n_FP,nX_G) )
    ALLOCATE( INFO(nX_G) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$OMP MAP( alloc: AMAT, BVEC, TAU, INFO )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$ACC CREATE( AMAT, BVEC, TAU, INFO )
#endif

    IF ( M_FP > 3 ) THEN

      CALL LinearLeastSquares_LWORK( 'N', n_FP, M_FP-1, 1, AMAT, n_FP, BVEC, n_FP, TMP, LWORK )
      ALLOCATE( WORK(LWORK,nX_G) )

    ELSE

      ALLOCATE( WORK(1,1) )

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: WORK )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( WORK )
#endif

  END SUBROUTINE InitializeNeutrinoMatterSolver


  SUBROUTINE FinalizeNeutrinoMatterSolver

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, &
    !$OMP               AMAT, BVEC, TAU, WORK, INFO )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, &
    !$ACC         AMAT, BVEC, TAU, WORK, INFO )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S )
    DEALLOCATE( AMAT, BVEC, TAU, WORK, INFO )

  END SUBROUTINE FinalizeNeutrinoMatterSolver


  SUBROUTINE ComputePointsAndWeightsE( E, W2, W3 )

    REAL(DP), INTENT(out) :: E (1:nE_G)
    REAL(DP), INTENT(out) :: W2(1:nE_G)
    REAL(DP), INTENT(out) :: W3(1:nE_G)

    INTEGER  :: iN_E, iE, iNodeE

    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      E (iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

      W2(iN_E) = WeightsE(iNodeE) * E(iN_E)**2 * MeshE % Width(iE)
      W3(iN_E) = WeightsE(iNodeE) * E(iN_E)**3 * MeshE % Width(iE)

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveMatterEquations_EmAb &
    ( J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

    ! --- Electron Neutrinos (1) and Electron Antineutrinos (2) ---

    REAL(DP), INTENT(in)    :: J_1  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J_2  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: D    (1:nX_G)
    REAL(DP), INTENT(in)    :: T    (1:nX_G)
    REAL(DP), INTENT(in)    :: Y    (1:nX_G)
    REAL(DP), INTENT(in)    :: E    (1:nX_G)

    REAL(DP), PARAMETER :: Log1d100 = LOG( 1.0d100 )
    REAL(DP) :: Mnu_1(1:nX_G), dMnudT_1(1:nX_G), dMnudY_1(1:nX_G)
    REAL(DP) :: Mnu_2(1:nX_G), dMnudT_2(1:nX_G), dMnudY_2(1:nX_G)
    REAL(DP) :: FD1_Exp, FD2_Exp

    INTEGER  :: iN_X, iN_E

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#endif

    ! --- Neutrino Chemical Potentials and Derivatives ---

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu_1, dMnudT_1, dMnudY_1, iSpecies = iNuE )

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu_2, dMnudT_2, dMnudY_2, iSpecies = iNuE_Bar )

    ! --- Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( FD1_Exp, FD2_Exp )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2), &
    !$ACC PRIVATE( FD1_Exp, FD2_Exp ) &
    !$ACC PRESENT( E_N, T, Mnu_1, Mnu_2, J0_1, J0_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( FD1_Exp, FD2_Exp )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        FD1_Exp = ( E_N(iN_E) - Mnu_1(iN_X) ) / ( BoltzmannConstant * T(iN_X) )
        FD1_Exp = MIN( MAX( FD1_Exp, - Log1d100 ), + Log1d100 )
        J0_1(iN_E,iN_X) = One / ( EXP( FD1_Exp ) + One )

        FD2_Exp = ( E_N(iN_E) - Mnu_2(iN_X) ) / ( BoltzmannConstant * T(iN_X) )
        FD2_Exp = MIN( MAX( FD2_Exp, - Log1d100 ), + Log1d100 )
        J0_2(iN_E,iN_X) = One / ( EXP( FD2_Exp ) + One )

      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Mnu_1, dMnudT_1, dMnudY_1, Mnu_2, dMnudT_2, dMnudY_2 )
#endif


  END SUBROUTINE SolveMatterEquations_EmAb


  SUBROUTINE SolveMatterEquations_FP_Coupled &
    ( dt, iS_1, iS_2, J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      D, T, Y, E, nIterations )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP), INTENT(in)    :: dt
    INTEGER,  INTENT(in)    :: iS_1, iS_2
    REAL(DP), INTENT(inout) :: J_1        (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: J_2        (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_1      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_2      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_1       (1:nE_G,1:nX_G)
    REAL(DP), INTENT(out)   :: J0_2       (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_NES_1  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_NES_2  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_NES_1  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_NES_2  (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_Pair_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Chi_Pair_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_Pair_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Eta_Pair_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: D          (1:nX_G)
    REAL(DP), INTENT(inout) :: T          (1:nX_G)
    REAL(DP), INTENT(inout) :: Y          (1:nX_G)
    REAL(DP), INTENT(inout) :: E          (1:nX_G)
    INTEGER,  INTENT(out)   :: nIterations(1:nX_G)

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = iE
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(              1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(              1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E
    REAL(DP), DIMENSION(1:nE_G,       1:nX_G) :: Jold_1, Jnew_1
    REAL(DP), DIMENSION(1:nE_G,       1:nX_G) :: Jold_2, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2

    INTEGER,  DIMENSION(       1:nX_G) :: PackedToUnpackedTable, UnpackedToPackedTable
    REAL(DP), DIMENSION(       1:nX_G) :: D_P, T_P, Y_P, E_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: J0_1_P, Jnew_1_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: J0_2_P, Jnew_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Eta_NES_1_P, Eta_NES_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Chi_NES_1_P, Chi_NES_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Eta_Pair_1_P, Eta_Pair_2_P
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Chi_Pair_1_P, Chi_Pair_2_P

    REAL(DP), DIMENSION(1:n_FP,1:M_FP,1:nX_G) :: GVEC, FVEC
    REAL(DP), DIMENSION(1:n_FP,       1:nX_G) :: GVECm, FVECm
    REAL(DP), DIMENSION(       1:M_FP,1:nX_G) :: Alpha

    REAL(DP), DIMENSION(    1:nX_G) :: JNRM_1, JNRM_2
    LOGICAL,  DIMENSION(    1:nX_G) :: CONVERGED, ITERATE
    INTEGER,  DIMENSION(    1:nX_G) :: Error

    REAL(DP) :: AERR_Y, AERR_E, AERR_J1, AERR_J2
    REAL(DP) :: RERR_Y, RERR_E, RERR_J1, RERR_J2
    REAL(DP) :: SUM1, SUM2
    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: i, k, iFP, iM, Mk, iN_X, iN_E, iX_P, nX_P
    INTEGER  :: OS_2

    OS_2 = OS_1 + nE_G
    ITERATE(:) = .TRUE.
    CONVERGED(:) = .FALSE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: CONVERGED, ITERATE ) &
    !$OMP MAP( alloc: Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP             Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP             Jold_1, Jnew_1, JNRM_1, &
    !$OMP             Jold_2, Jnew_2, JNRM_2, &
    !$OMP             Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP             Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP             PackedToUnpackedTable, UnpackedToPackedTable, &
    !$OMP             D_P, T_P, Y_P, E_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$OMP             Eta_NES_1_P, Eta_NES_2_P, &
    !$OMP             Chi_NES_1_P, Chi_NES_2_P, &
    !$OMP             Eta_Pair_1_P, Eta_Pair_2_P, &
    !$OMP             Chi_Pair_1_P, Chi_Pair_2_P, &
    !$OMP             GVEC, FVEC, GVECm, FVECm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( CONVERGED, ITERATE ) &
    !$ACC CREATE( Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, JNRM_1, &
    !$ACC         Jold_2, Jnew_2, JNRM_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackedToUnpackedTable, UnpackedToPackedTable, &
    !$ACC         D_P, T_P, Y_P, E_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$ACC         Eta_NES_1_P, Eta_NES_2_P, &
    !$ACC         Chi_NES_1_P, Chi_NES_2_P, &
    !$ACC         Eta_Pair_1_P, Eta_Pair_2_P, &
    !$ACC         Chi_Pair_1_P, Chi_Pair_2_P, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( Y, Yold, E, Eold )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      Yold(iN_X) = Y(iN_X)
      Eold(iN_X) = E(iN_X)
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( Jold_1, Jold_2, J_1, J_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        Jold_1(iN_E,iN_X) = J_1(iN_E,iN_X)
        Jold_2(iN_E,iN_X) = J_2(iN_E,iN_X)
      END DO
    END DO

    !IF( UsePreconditionerEmAb )THEN
    !  CALL TimersStart( Timer_Im_EmAb_FP )
    !  CALL SolveMatterEquations_EmAb_FP &
    !         ( dt, iS_1, iS_2, J_1, J_2, &
    !           Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E)
    !  CALL TimersStop( Timer_Im_EmAb_FP )
    !END IF

    CALL TimersStart( Timer_Im_ComputeOpacity )

    ! --- Compute Opacity Kernels ---

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, J0_1, iS_1 )

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, J0_2, iS_2 )

    ! --- NES Kernels ---

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_1, 1, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1 )

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_2, 1, &
             Phi_0_In_NES_2, Phi_0_Ot_NES_2 )

    ! --- Pair Kernels ---

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_1, 1, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1 )

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX_G, E_N, D, T, Y, iS_2, 1, &
             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

    CALL TimersStop( Timer_Im_ComputeOpacity )

    !IF( UsePreconditionerPair .OR. UsePreconditionerPairLagAllButJ0 )THEN
    !  CALL TimersStart( Timer_Im_Presolve )
    !  CALL SolveMatterEquations_Presolve &
    !         ( dt, iS_1, iS_2, J_1, J_2, &
    !           Chi_1, Chi_2, J0_1, J0_2, &
    !           Phi_0_In_NES_1, Phi_0_Ot_NES_1, 
    !           Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !           Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !           Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !           D, T, Y, E, &
    !           UsePreconditionerPairLagAllButJ0 )
    !  CALL TimersStop( Timer_Im_Presolve )
    !END IF

    ! --- Initial RHS ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_1, nE_G, W2_S, 1, Zero, C_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, Jold_2, nE_G, W2_S, 1,  One, C_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_1, nE_G, W3_S, 1, Zero, C_E, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_2, nE_G, W3_S, 1,  One, C_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1, SUM2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( SUM1, SUM2 ) &
    !$ACC PRESENT( Unew_Y, Y, Yold, S_Y, C_Y, &
    !$ACC          Unew_E, E, Eold, S_E, C_E, D, &
    !$ACC          JNRM_1, JNRM_2, Jold_1, Jold_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1, SUM2 )
#endif
    DO iN_X = 1, nX_G

      S_Y(iN_X) = One / ( D(iN_X) * Yold(iN_X) / AtomicMassUnit )
      S_E(iN_X) = One / ( D(iN_X) * Eold(iN_X) )

      C_Y(iN_X) = C_Y(iN_X) * S_Y(iN_X)
      C_E(iN_X) = C_E(iN_X) * S_E(iN_X)

      Unew_Y(iN_X) = Y(iN_X) / Yold(iN_X) ! --- Initial Guess
      Unew_E(iN_X) = E(iN_X) / Eold(iN_X) ! --- Initial Guess

      SUM1 = Zero
      SUM2 = Zero
      DO iN_E = 1, nE_G
        SUM1 = SUM1 + Jold_1(iN_E,iN_X) * Jold_1(iN_E,iN_X)
        SUM2 = SUM2 + Jold_2(iN_E,iN_X) * Jold_2(iN_E,iN_X)
      END DO
      JNRM_1(iN_X) = SQRT( SUM1 )
      JNRM_2(iN_X) = SQRT( SUM2 )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( Jnew_1, Jnew_2, J_1, J_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        Jnew_1(iN_E,iN_X) = J_1(iN_E,iN_X) ! --- Initial Guess
        Jnew_2(iN_E,iN_X) = J_2(iN_E,iN_X) ! --- Initial Guess
      END DO
    END DO

    !IF( .NOT. (UsePreconditionerPair .OR. UsePreconditionerPairLagAllButJ0 ))THEN

      ! --- Compute Neutrino Rates ---

      CALL TimersStart( Timer_Im_ComputeRate )

      ! --- NES Emissivities and Opacities ---

      CALL ComputeNeutrinoOpacitiesRates_NES_Points &
             ( 1, nE_G, 1, nX_G, W2_N, Jnew_1, &
               Phi_0_In_NES_1, Phi_0_Ot_NES_1, Eta_NES_1, Chi_NES_1 )

      CALL ComputeNeutrinoOpacitiesRates_NES_Points &
             ( 1, nE_G, 1, nX_G, W2_N, Jnew_2, &
               Phi_0_In_NES_2, Phi_0_Ot_NES_2, Eta_NES_2, Chi_NES_2 )

      ! --- Pair Emissivities and Opacities ---

      CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
             ( 1, nE_G, 1, nX_G, W2_N, Jnew_2, &
               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Eta_Pair_1, Chi_Pair_1 )

      CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
             ( 1, nE_G, 1, nX_G, W2_N, Jnew_1, &
               Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, Eta_Pair_2, Chi_Pair_2 )

      CALL TimersStop( Timer_Im_ComputeRate )

      ! --- Update Neutrino Densities ---

      CALL UpdateJ_FP &
             ( ITERATE, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
               Chi_1, Chi_2, J0_1, J0_2, &
               Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
               Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    !END IF

    k = 0
    DO WHILE( ANY( ITERATE(:) ) .AND. k < MaxIter )

      k  = k + 1
      Mk = MIN( M_FP, k )
      iM = Mk
      !iM = 1 + MOD( k-1, M_FP )

      CALL CreatePackTable_FP &
             ( ITERATE, nX_P, PackedToUnpackedTable, UnpackedToPackedTable )

      IF ( k > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Im_ComputeOpacity )

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
        !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR &
        !$ACC PRIVATE( iN_X ) &
        !$ACC PRESENT( UnpackedToPackedTable, D, T, Y, D_P, T_P, Y_P )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO SIMD &
        !$OMP PRIVATE( iN_X )
#endif
        DO iX_P = 1, nX_P
          iN_X = UnpackedToPackedTable(iX_P)
          D_P(iX_P) = D(iN_X)
          T_P(iX_P) = T(iN_X)
          Y_P(iX_P) = Y(iN_X)
        END DO

        ! --- Equilibrium Distributions ---

        CALL ComputeEquilibriumDistributions_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), &
                 J0_1_P(:,1:nX_P), iS_1 )

        CALL ComputeEquilibriumDistributions_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), &
                 J0_2_P(:,1:nX_P), iS_2 )

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( iX_P )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
        !$ACC PRIVATE( iX_P ) &
        !$ACC PRESENT( ITERATE, PackedToUnpackedTable, J0_1, J0_2, J0_1_P, J0_2_P )
#elif defined(THORNADO_OMP)
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( iX_P )
#endif
        DO iN_X = 1, nX_G
          DO iN_E = 1, nE_G
            IF ( ITERATE(iN_X) ) THEN
              iX_P = PackedToUnpackedTable(iN_X)
              J0_1(iN_E,iN_X) = J0_1_P(iN_E,iX_P)
              J0_2(iN_E,iN_X) = J0_2_P(iN_E,iX_P)
            END IF
          END DO
        END DO

        ! --- NES Kernels ---

        CALL ComputeNeutrinoOpacities_NES_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_1, 1, &
                 Phi_0_In_NES_1(:,:,1:nX_P), Phi_0_Ot_NES_1(:,:,1:nX_P) )

        CALL ComputeNeutrinoOpacities_NES_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_2, 1, &
                 Phi_0_In_NES_2(:,:,1:nX_P), Phi_0_Ot_NES_2(:,:,1:nX_P) )

        ! --- Pair Kernels ---

        CALL ComputeNeutrinoOpacities_Pair_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_1, 1, &
                 Phi_0_In_Pair_1(:,:,1:nX_P), Phi_0_Ot_Pair_1(:,:,1:nX_P) )

        CALL ComputeNeutrinoOpacities_Pair_Points &
               ( 1, nE_G, 1, nX_P, E_N, &
                 D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), iS_2, 1, &
                 Phi_0_In_Pair_2(:,:,1:nX_P), Phi_0_Ot_Pair_2(:,:,1:nX_P) )

        CALL TimersStop( Timer_Im_ComputeOpacity )

      END IF

      ! --- NES Emissivities and Opacities ---

      CALL TimersStart( Timer_Im_ComputeRate )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( UnpackedToPackedTable, Jnew_1, Jnew_2, Jnew_1_P, Jnew_2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iN_X )
#endif
      DO iX_P = 1, nX_P
        DO iN_E = 1, nE_G
          iN_X = UnpackedToPackedTable(iX_P)
          Jnew_1_P(iN_E,iX_P) = Jnew_1(iN_E,iN_X)
          Jnew_2_P(iN_E,iX_P) = Jnew_2(iN_E,iN_X)
        END DO
      END DO

      CALL ComputeNeutrinoOpacitiesRates_NES_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_1_P(:,1:nX_P), &
               Phi_0_In_NES_1(:,:,1:nX_P), Phi_0_Ot_NES_1(:,:,1:nX_P), &
               Eta_NES_1_P(:,1:nX_P), Chi_NES_1_P(:,1:nX_P) )

      CALL ComputeNeutrinoOpacitiesRates_NES_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_2_P(:,1:nX_P), &
               Phi_0_In_NES_2(:,:,1:nX_P), Phi_0_Ot_NES_2(:,:,1:nX_P), &
               Eta_NES_2_P(:,1:nX_P), Chi_NES_2_P(:,1:nX_P) )

      ! --- Pair Emissivities and Opacities ---

      CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_2_P(:,1:nX_P), &
               Phi_0_In_Pair_1(:,:,1:nX_P), Phi_0_Ot_Pair_1(:,:,1:nX_P), &
               Eta_Pair_1_P(:,1:nX_P), Chi_Pair_1_P(:,1:nX_P) )

      CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
             ( 1, nE_G, 1, nX_P, W2_N, Jnew_1_P(:,1:nX_P), &
               Phi_0_In_Pair_2(:,:,1:nX_P), Phi_0_Ot_Pair_2(:,:,1:nX_P), &
               Eta_Pair_2_P(:,1:nX_P), Chi_Pair_2_P(:,1:nX_P) )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iX_P )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iX_P ) &
      !$ACC PRESENT( ITERATE, PackedToUnpackedTable, &
      !$ACC          Chi_NES_1, Chi_NES_2, Chi_NES_1_P, Chi_NES_2_P, &
      !$ACC          Eta_NES_1, Eta_NES_2, Eta_NES_1_P, Eta_NES_2_P, &
      !$ACC          Chi_Pair_1, Chi_Pair_2, Chi_Pair_1_P, Chi_Pair_2_P, &
      !$ACC          Eta_Pair_1, Eta_Pair_2, Eta_Pair_1_P, Eta_Pair_2_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iX_P )
#endif
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G
          IF ( ITERATE(iN_X) ) THEN
            iX_P = PackedToUnpackedTable(iN_X)
            Chi_NES_1(iN_E,iN_X) = Chi_NES_1_P(iN_E,iX_P)
            Chi_NES_2(iN_E,iN_X) = Chi_NES_2_P(iN_E,iX_P)
            Eta_NES_1(iN_E,iN_X) = Eta_NES_1_P(iN_E,iX_P)
            Eta_NES_2(iN_E,iN_X) = Eta_NES_2_P(iN_E,iX_P)
            Chi_Pair_1(iN_E,iN_X) = Chi_Pair_1_P(iN_E,iX_P)
            Chi_Pair_2(iN_E,iN_X) = Chi_Pair_2_P(iN_E,iX_P)
            Eta_Pair_1(iN_E,iN_X) = Eta_Pair_1_P(iN_E,iX_P)
            Eta_Pair_2(iN_E,iN_X) = Eta_Pair_2_P(iN_E,iX_P)
          END IF
        END DO
      END DO

      CALL TimersStop( Timer_Im_ComputeRate )

      ! --- Right-Hand Side Vectors and Residuals ---

      CALL ComputeRHS_FP &
             ( ITERATE, Mk, iY, iE, OS_1, OS_2, &
               C_Y, S_Y, Unew_Y, GVEC_Y, C_E, S_E, Unew_E, GVEC_E, &
               dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
               Chi_1, Chi_2, J0_1, J0_2, &
               Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
               Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
               FVECm, GVECm, FVEC, GVEC )

      ! --- Anderson Acceleration ---

      CALL TimersStart( Timer_Im_ComputeLS )

      IF ( Mk > 1 ) THEN

        CALL BuildLS_FP &
               ( ITERATE, Mk, FVECm, FVEC, AMAT, BVEC )

        CALL SolveLS_FP &
              ( ITERATE, Mk, AMAT, BVEC, Alpha )

        CALL UpdateRHS_AA_FP &
               ( ITERATE, Mk, Alpha, GVEC, GVECm )

      END IF

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors---

      CALL TimersStart( Timer_Im_UpdateFP )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( ITERATE, FVECm, GVECm, Y, E, Yold, Eold, &
      !$ACC          Unew_Y, Unew_E, Jnew_1, Jnew_2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          FVECm(iY,iN_X) = GVECm(iY,iN_X) - Unew_Y(iN_X)
          FVECm(iE,iN_X) = GVECm(iE,iN_X) - Unew_E(iN_X)

          DO iFP = OS_1+1, OS_1+nE_G
            FVECm(iFP,iN_X) = GVECm(iFP,iN_X) - Jnew_1(iFP-OS_1,iN_X)
          END DO

          DO iFP = OS_2+1, OS_2+nE_G
            FVECm(iFP,iN_X) = GVECm(iFP,iN_X) - Jnew_2(iFP-OS_2,iN_X)
          END DO

          Unew_Y(iN_X) = GVECm(iY,iN_X)
          Y(iN_X) = Unew_Y(iN_X) * Yold(iN_X)

          Unew_E(iN_X) = GVECm(iE,iN_X)
          E(iN_X) = Unew_E(iN_X) * Eold(iN_X)

          DO iN_E = 1, nE_G
            Jnew_1(iN_E,iN_X) = GVECm(iN_E+OS_1,iN_X)
            Jnew_2(iN_E,iN_X) = GVECm(iN_E+OS_2,iN_X)
          END DO

        END IF
      END DO

      ! --- Update Temperature ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( UnpackedToPackedTable, D, T, Y, E, D_P, T_P, Y_P, E_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( iN_X )
#endif
      DO iX_P = 1, nX_P
        iN_X = UnpackedToPackedTable(iX_P)
        D_P(iX_P) = D(iN_X)
        T_P(iX_P) = T(iN_X)
        Y_P(iX_P) = Y(iN_X)
        E_P(iX_P) = E(iN_X)
      END DO

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D_P(1:nX_P), E_P(1:nX_P), Y_P(1:nX_P), T_P(1:nX_P) )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX_P )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iX_P ) &
      !$ACC PRESENT( ITERATE, PackedToUnpackedTable, T, T_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( iX_P )
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN
          iX_P = PackedToUnpackedTable(iN_X)
          T(iN_X) = T_P(iX_P)
        END IF
      END DO

      ! --- Check Convergence ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( AERR_Y, AERR_E, AERR_J1, AERR_J2, &
      !$OMP          RERR_Y, RERR_E, RERR_J1, RERR_J2 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( AERR_Y, AERR_E, AERR_J1, AERR_J2, &
      !$ACC          RERR_Y, RERR_E, RERR_J1, RERR_J2 ) &
      !$ACC PRESENT( ITERATE, CONVERGED, nIterations, JNRM_1, JNRM_2, FVECm )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( AERR_Y, AERR_E, AERR_J1, AERR_J2, &
      !$OMP          RERR_Y, RERR_E, RERR_J1, RERR_J2 )
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          AERR_Y = ABS( FVECm(iY,iN_X) )
          AERR_E = ABS( FVECm(iE,iN_X) )

          AERR_J1 = Zero
          AERR_J2 = Zero
          DO iN_E = 1, nE_G
            AERR_J1 = AERR_J1 + FVECm(OS_1+iN_E,iN_X)**2
            AERR_J2 = AERR_J2 + FVECm(OS_2+iN_E,iN_X)**2
          END DO
          AERR_J1 = SQRT( AERR_J1 )
          AERR_J2 = SQRT( AERR_J2 )

          RERR_Y  = AERR_Y
          RERR_E  = AERR_E
          RERR_J1 = AERR_J1 / JNRM_1(iN_X)
          RERR_J2 = AERR_J2 / JNRM_2(iN_X)

          CONVERGED(iN_X) = RERR_Y  <= Rtol &
                      .AND. RERR_E  <= Rtol &
                      .AND. RERR_J1 <= Rtol &
                      .AND. RERR_J2 <= Rtol

          ITERATE(iN_X) = .NOT. CONVERGED(iN_X)

          IF ( CONVERGED(iN_X) ) nIterations(iN_X) = k

        END IF
      END DO

      ! --- Shift History Arrays ---

      IF ( Mk == M_FP ) THEN

        CALL ShiftRHS_FP &
               ( ITERATE, Mk, FVEC, GVEC )

      END IF

      CALL TimersStop( Timer_Im_UpdateFP )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( ITERATE, CONVERGED )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( ITERATE, CONVERGED )
#endif

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( Jnew_1, Jnew_2, J_1, J_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        J_1(iN_E,iN_X) = Jnew_1(iN_E,iN_X)
        J_2(iN_E,iN_X) = Jnew_2(iN_E,iN_X)
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: CONVERGED, ITERATE, &
    !$OMP               Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP               Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP               Jold_1, Jnew_1, JNRM_1, &
    !$OMP               Jold_2, Jnew_2, JNRM_2, &
    !$OMP               Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP               Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP               Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP               PackedToUnpackedTable, UnpackedToPackedTable, &
    !$OMP               D_P, T_P, Y_P, E_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$OMP               Eta_NES_1_P, Eta_NES_2_P, &
    !$OMP               Chi_NES_1_P, Chi_NES_2_P, &
    !$OMP               Eta_Pair_1_P, Eta_Pair_2_P, &
    !$OMP               Chi_Pair_1_P, Chi_Pair_2_P, &
    !$OMP               GVEC, FVEC, GVECm, FVECm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( CONVERGED, ITERATE, &
    !$ACC         Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, JNRM_1, &
    !$ACC         Jold_2, Jnew_2, JNRM_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackedToUnpackedTable, UnpackedToPackedTable, &
    !$ACC         D_P, T_P, Y_P, E_P, J0_1_P, Jnew_1_P, J0_2_P, Jnew_2_P, &
    !$ACC         Eta_NES_1_P, Eta_NES_2_P, &
    !$ACC         Chi_NES_1_P, Chi_NES_2_P, &
    !$ACC         Eta_Pair_1_P, Eta_Pair_2_P, &
    !$ACC         Chi_Pair_1_P, Chi_Pair_2_P, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha )
#endif

  END SUBROUTINE SolveMatterEquations_FP_Coupled


  SUBROUTINE CreatePackTable_FP &
    ( MASK, nX_P, PackedToUnpackedTable, UnpackedToPackedTable )

    LOGICAL,  INTENT(in)  :: MASK(1:nX_G)
    INTEGER,  INTENT(out) :: nX_P
    INTEGER,  INTENT(out) :: PackedToUnpackedTable(1:nX_G)
    INTEGER,  INTENT(out) :: UnpackedToPackedTable(1:nX_G)

    INTEGER  :: iN_X, iX_P, iX_U

    ! --- Build Lookup Tables ---

    nX_P = COUNT( MASK(:) )
    iX_P = 0
    iX_U = nX_P
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN
        iX_P = iX_P + 1
        PackedToUnpackedTable(iN_X) = iX_P
        UnpackedToPackedTable(iX_P) = iN_X
      ELSE
        iX_U = iX_U + 1
        PackedToUnpackedTable(iN_X) = iX_U
        UnpackedToPackedTable(iX_U) = iN_X
      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( PackedToUnpackedTable, UnpackedToPackedTable )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( PackedToUnpackedTable, UnpackedToPackedTable )
#endif

  END SUBROUTINE CreatePackTable_FP


  SUBROUTINE UpdateJ_FP &
    ( MASK, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    LOGICAL,  INTENT(in)    :: MASK      (1:nX_G)
    REAL(DP), INTENT(in)    :: dt
    REAL(DP), INTENT(in)    :: Jold_1    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Jold_2    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Jnew_1    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Jnew_2    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_1     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_2     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J0_1      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J0_2      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_NES_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_NES_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_NES_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_NES_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_Pair_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_Pair_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_Pair_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_Pair_2(1:nE_G,1:nX_G)

    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T ) &
    !$ACC PRESENT( MASK, &
    !$ACC          Jold_1, Jnew_1, J0_1, Chi_1, Chi_NES_1, Eta_NES_1, Chi_Pair_1, Eta_Pair_1, &
    !$ACC          Jold_2, Jnew_2, J0_2, Chi_2, Chi_NES_2, Eta_NES_2, Chi_Pair_2, Eta_Pair_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Eta = Chi_1(iN_E,iN_X) * J0_1(iN_E,iN_X)
          Eta_T = Eta + Eta_NES_1(iN_E,iN_X) + Eta_Pair_1(iN_E,iN_X)
          Chi_T = Chi_1(iN_E,iN_X) + Chi_NES_1(iN_E,iN_X) + Chi_Pair_1(iN_E,iN_X)
          Jnew_1(iN_E,iN_X) = ( Jold_1(iN_E,iN_X) + dt * Eta_T ) / ( One + dt * Chi_T )

          Eta = Chi_2(iN_E,iN_X) * J0_2(iN_E,iN_X)
          Eta_T = Eta + Eta_NES_2(iN_E,iN_X) + Eta_Pair_2(iN_E,iN_X)
          Chi_T = Chi_2(iN_E,iN_X) + Chi_NES_2(iN_E,iN_X) + Chi_Pair_2(iN_E,iN_X)
          Jnew_2(iN_E,iN_X) = ( Jold_2(iN_E,iN_X) + dt * Eta_T ) / ( One + dt * Chi_T )

        END IF
      END DO
    END DO

  END SUBROUTINE UpdateJ_FP


  SUBROUTINE ComputeRHS_FP &
    ( MASK, Mk, iY, iE, OS_1, OS_2, &
      C_Y, S_Y, Unew_Y, G_Y, C_E, S_E, Unew_E, G_E, &
      dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      Fm, Gm, F, G )

    LOGICAL,  INTENT(in)    :: MASK      (1:nX_G)
    INTEGER,  INTENT(in)    :: Mk, iY, iE, OS_1, OS_2
    REAL(DP), INTENT(in)    :: C_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: S_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: Unew_Y    (1:nX_G)
    REAL(DP), INTENT(out)   :: G_Y       (1:nX_G)
    REAL(DP), INTENT(in)    :: C_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: S_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: Unew_E    (1:nX_G)
    REAL(DP), INTENT(out)   :: G_E       (1:nX_G)
    REAL(DP), INTENT(in)    :: dt
    REAL(DP), INTENT(in)    :: Jold_1    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Jold_2    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Jnew_1    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Jnew_2    (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_1     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_2     (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J0_1      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: J0_2      (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_NES_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_NES_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_NES_1 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_NES_2 (1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_Pair_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Chi_Pair_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_Pair_1(1:nE_G,1:nX_G)
    REAL(DP), INTENT(in)    :: Eta_Pair_2(1:nE_G,1:nX_G)
    REAL(DP), INTENT(inout) :: Fm        (1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: Gm        (1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: F         (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: G         (1:n_FP,1:M_FP,1:nX_G)

    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_X, iFP

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jnew_1, nE_G, W2_S, 1, Zero, G_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, Jnew_2, nE_G, W2_S, 1,  One, G_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jnew_1, nE_G, W3_S, 1, Zero, G_E, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jnew_2, nE_G, W3_S, 1,  One, G_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, C_Y, S_Y, Unew_Y, G_Y, C_E, S_E, Unew_E, G_E, Fm, Gm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_Y(iN_X) = One + C_Y(iN_X) - G_Y(iN_X) * S_Y(iN_X)
        G_E(iN_X) = One + C_E(iN_X) - G_E(iN_X) * S_E(iN_X)

        Gm(iY,iN_X) = G_Y(iN_X)
        Gm(iE,iN_X) = G_E(iN_X)

        Fm(iY,iN_X) = G_Y(iN_X) - Unew_Y(iN_X)
        Fm(iE,iN_X) = G_E(iN_X) - Unew_E(iN_X)

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T ) &
    !$ACC PRESENT( MASK, Fm, Gm, &
    !$ACC          Jold_1, Jnew_1, J0_1, Chi_1, Chi_NES_1, Eta_NES_1, Chi_Pair_1, Eta_Pair_1, &
    !$ACC          Jold_2, Jnew_2, J0_2, Chi_2, Chi_NES_2, Eta_NES_2, Chi_Pair_2, Eta_Pair_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Eta = Chi_1(iN_E,iN_X) * J0_1(iN_E,iN_X)
          Eta_T = Eta + Eta_NES_1(iN_E,iN_X) + Eta_Pair_1(iN_E,iN_X)
          Chi_T = Chi_1(iN_E,iN_X) + Chi_NES_1(iN_E,iN_X) + Chi_Pair_1(iN_E,iN_X)
          Gm(OS_1+iN_E,iN_X) = ( Jold_1(iN_E,iN_X) + dt * Eta_T ) / ( One + dt * Chi_T )
          Fm(OS_1+iN_E,iN_X) = Gm(OS_1+iN_E,iN_X) - Jnew_1(iN_E,iN_X)

          Eta = Chi_2(iN_E,iN_X) * J0_2(iN_E,iN_X)
          Eta_T = Eta + Eta_NES_2(iN_E,iN_X) + Eta_Pair_2(iN_E,iN_X)
          Chi_T = Chi_2(iN_E,iN_X) + Chi_NES_2(iN_E,iN_X) + Chi_Pair_2(iN_E,iN_X)
          Gm(OS_2+iN_E,iN_X) = ( Jold_2(iN_E,iN_X) + dt * Eta_T ) / ( One + dt * Chi_T )
          Fm(OS_2+iN_E,iN_X) = Gm(OS_2+iN_E,iN_X) - Jnew_2(iN_E,iN_X)

        END IF
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, Fm, Gm, F, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iFP = 1, n_FP
        IF ( MASK(iN_X) ) THEN
          F(iFP,Mk,iN_X) = Fm(iFP,iN_X)
          G(iFP,Mk,iN_X) = Gm(iFP,iN_X)
        END IF
      END DO
    END DO

  END SUBROUTINE ComputeRHS_FP


  SUBROUTINE BuildLS_FP &
    ( MASK, Mk, Fm, F, AMAT, BVEC )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(in)    :: Fm  (1:n_FP,       1:nX_G)
    REAL(DP), INTENT(in)    :: F   (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: AMAT(1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: BVEC(1:n_FP,       1:nX_G)

    INTEGER  :: iN_X, iFP, i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, Fm, BVEC )
#elif defined(THORNADO_OMP)
    !$OMP DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iFP = 1, n_FP
        IF ( MASK(iN_X) ) THEN
          BVEC(iFP,iN_X) = - Fm(iFP,iN_X)
        END IF
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( MASK, F, Fm, AMAT )
#elif defined(THORNADO_OMP)
    !$OMP DO SIMD COLLAPSE(3)
#endif
    DO iN_X = 1, nX_G
      DO i = 1, Mk-1
        DO iFP = 1, n_FP
          IF ( MASK(iN_X) ) THEN
            AMAT(iFP,i,iN_X) = F(iFP,i,iN_X) - Fm(iFP,iN_X)
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE BuildLS_FP


  SUBROUTINE SolveLS_FP &
    ( MASK, Mk, AMAT, BVEC, Alpha )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(inout) :: AMAT (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: BVEC (1:n_FP,       1:nX_G)
    REAL(DP), INTENT(inout) :: Alpha(       1:M_FP,1:nX_G)

    REAL(DP) :: AA11, AA12, AA22, AB1, AB2, DET_AA, SUM1
    INTEGER  :: iN_X, iFP, i

    IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( AA11, AB1 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( AA11, AB1 ) &
      !$ACC PRESENT( MASK, AMAT, BVEC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( AA11, AB1 )
#endif
      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          AA11 = Zero
          AB1 = Zero

#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO SIMD &
          !$OMP REDUCTION( +: AA11, AB1 )
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR &
          !$ACC REDUCTION( +: AA11, AB1 )
#endif
          DO iFP = 1, n_FP

            AA11 = AA11 + AMAT(iFP,1,iN_X) * AMAT(iFP,1,iN_X)
            AB1  = AB1  + AMAT(iFP,1,iN_X) * BVEC(iFP,iN_X)

          END DO

          BVEC(1,iN_X) = AB1 / AA11

        END IF
      END DO

    ELSE IF ( Mk == 3 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA ) &
      !$ACC PRESENT( MASK, AMAT, BVEC )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#endif
      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          AA11 = Zero
          AA12 = Zero
          AA22 = Zero
          AB1  = Zero
          AB2  = Zero

#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO SIMD &
          !$OMP REDUCTION( +: AA11, AA12, AA22, AB1, AB2 )
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR &
          !$ACC REDUCTION( +: AA11, AA12, AA22, AB1, AB2 )
#endif
          DO iFP = 1, n_FP

            AA11 = AA11 + AMAT(iFP,1,iN_X) * AMAT(iFP,1,iN_X)
            AA12 = AA12 + AMAT(iFP,1,iN_X) * AMAT(iFP,2,iN_X)
            AA22 = AA22 + AMAT(iFP,2,iN_X) * AMAT(iFP,2,iN_X)

            AB1  = AB1  + AMAT(iFP,1,iN_X) * BVEC(iFP,iN_X)
            AB2  = AB2  + AMAT(iFP,2,iN_X) * BVEC(iFP,iN_X)

          END DO

          DET_AA = AA11*AA22 - AA12*AA12

          BVEC(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
          BVEC(2,iN_X) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

        END IF
      END DO

    ELSE

      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          CALL LinearLeastSquares &
            ( 'N', n_FP, Mk-1, 1, AMAT(1,1,iN_X), n_FP, &
              BVEC(1,iN_X), n_FP, TAU(1,iN_X), WORK(1,iN_X), LWORK, INFO(iN_X) )

        END IF
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( SUM1 ) &
    !$ACC PRESENT( MASK, Alpha, BVEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        SUM1 = Zero
        DO i = 1, Mk-1
          Alpha(i,iN_X) = BVEC(i,iN_X)
          SUM1 = SUM1 + BVEC(i,iN_X)
        END DO
        Alpha(Mk,iN_X) = One - SUM1

      END IF
    END DO

  END SUBROUTINE SolveLS_FP


  SUBROUTINE UpdateRHS_AA_FP &
    ( Mask, Mk, Alpha, G, Gm )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(in)    :: Alpha(       1:M_FP,1:nX_G)
    REAL(DP), INTENT(in)    :: G    (1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: Gm   (1:n_FP,       1:nX_G)

    REAL(DP) :: SUM1
    INTEGER  :: iN_X, iFP, i

    !CALL MatrixMatrixMultiplyBatched &
    !       ( 'N', 'N', n_FP, 1, Mk, Zero, G, n_FP, n_FP*M_FP, &
    !         Alpha, M_FP, M_FP, Gm, n_FP, n_FP, nX_G )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( SUM1 ) &
    !$ACC PRESENT( MASK, G, Gm, Alpha )
#elif defined(THORNADO_OMP)
    !$OMP DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( SUM1 )
#endif
    DO iN_X = 1, nX_G
      DO iFP = 1, n_FP
        IF ( MASK(iN_X) ) THEN

          SUM1 = Zero
          DO i = 1, Mk
            SUM1 = SUM1 + G(iFP,i,iN_X) * Alpha(i,iN_X)
          END DO
          Gm(iFP,iN_X) = SUM1

        END IF
      END DO
    END DO

  END SUBROUTINE UpdateRHS_AA_FP


  SUBROUTINE ShiftRHS_FP &
    ( MASK, Mk, F, G )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(inout) :: F(1:n_FP,1:M_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: G(1:n_FP,1:M_FP,1:nX_G)

    REAL(DP) :: FTMP(1:n_FP,1:M_FP)
    REAL(DP) :: GTMP(1:n_FP,1:M_FP)
    INTEGER  :: iN_X, iFP, i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( FTMP, GTMP )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( FTMP, GTMP ) &
    !$ACC PRESENT( MASK, F, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( FTMP, GTMP )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
        !$ACC LOOP VECTOR COLLAPSE(2)
#endif
        DO i = 1, Mk-1
          DO iFP = 1, n_FP
            FTMP(iFP,i) = F(iFP,i+1,iN_X)
            GTMP(iFP,i) = G(iFP,i+1,iN_X)
          END DO
        END DO

#if defined(THORNADO_OMP_OL)
        !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
        !$ACC LOOP VECTOR COLLAPSE(2)
#endif
        DO i = 1, Mk-1
          DO iFP = 1, n_FP
            F(iFP,i,iN_X) = FTMP(iFP,i)
            G(iFP,i,iN_X) = GTMP(iFP,i)
          END DO
        END DO

      END IF
    END DO

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE ComputeNeutrinoChemicalPotentials &
    ( D, T, Y, M, dMdT, dMdY, iSpecies )

    REAL(DP), INTENT(in)  :: D   (1:nX_G)
    REAL(DP), INTENT(in)  :: T   (1:nX_G)
    REAL(DP), INTENT(in)  :: Y   (1:nX_G)
    REAL(DP), INTENT(out) :: M   (1:nX_G)
    REAL(DP), INTENT(out) :: dMdT(1:nX_G)
    REAL(DP), INTENT(out) :: dMdY(1:nX_G)
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: dMdD(1:nX_G)
    REAL(DP) :: Me(1:nX_G), dMedT(1:nX_G), dMedY(1:nX_G)
    REAL(DP) :: Mp(1:nX_G), dMpdT(1:nX_G), dMpdY(1:nX_G)
    REAL(DP) :: Mn(1:nX_G), dMndT(1:nX_G), dMndY(1:nX_G)

    INTEGER :: i

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#endif

    ! --- Matter Chemical Potentials and Derivatives ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, M = Me, dMdD_Option = dMdD, dMdT_Option = dMedT, dMdY_Option = dMedY )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, M = Mp, dMdD_Option = dMdD, dMdT_Option = dMpdT, dMdY_Option = dMpdY )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, M = Mn, dMdD_Option = dMdD, dMdT_Option = dMndT, dMdY_Option = dMndY )

    ! --- Neutrino Chemical Potential and Derivatives ---

    IF ( iSpecies == iNuE )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( M, dMdT, dMdY, Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO i = 1, nX_G
        M   (i) = ( Me   (i) + Mp   (i) ) - Mn   (i)
        dMdT(i) = ( dMedT(i) + dMpdT(i) ) - dMndT(i)
        dMdY(i) = ( dMedY(i) + dMpdY(i) ) - dMndY(i)
      END DO

    ELSE IF ( iSpecies == iNuE_Bar )THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( M, dMdT, dMdY, Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO i = 1, nX_G
        M   (i) = Mn   (i) - ( Me   (i) + Mp   (i) )
        dMdT(i) = dMndT(i) - ( dMedT(i) + dMpdT(i) )
        dMdY(i) = dMndY(i) - ( dMedY(i) + dMpdY(i) )
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Me, dMedT, dMedY, Mp, dMpdT, dMpdY, Mn, dMndT, dMndY, dMdD )
#endif

  END SUBROUTINE ComputeNeutrinoChemicalPotentials


END MODULE TwoMoment_NeutrinoMatterSolverModule