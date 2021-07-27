#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_NeutrinoMatterSolverModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, One, Two, FourPi, Half
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    AtomicMassUnit, &
    Centimeter, &
    Gram, &
    MeV, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, &
    nDOFE
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Collisions_InnerLoop, &
    Timer_Collisions_ComputeOpacity, &
    Timer_Collisions_ComputeRates, &
    Timer_Collisions_InitializeRHS, &
    Timer_Collisions_NeutrinoRHS, &
    Timer_Collisions_MatterRHS, &
    Timer_Collisions_SolveLS, &
    Timer_Collisions_UpdateFP, &
    Timer_Collisions_CheckOuter, &
    Timer_Collisions_CheckInner
  USE ArrayUtilitiesModule, ONLY: &
    CreatePackIndex, &
    ArrayPack, &
    ArrayUnpack, &
    ArrayCopy
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    MatrixVectorMultiply, &
    LinearLeastSquares, &
    LinearLeastSquares_LWORK, &
    LinearSolveBatched
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    iNuE, &
    iNuE_Bar, &
    LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_DG_Points, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacities_Pair_Points, &
    ComputeNeutrinoOpacitiesRates_NES_Points, &
    ComputeNeutrinoOpacitiesRates_Pair_Points
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    EddingtonTensorComponents_dd

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: InitializeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling_FP_Nested_AA

  LOGICAL, PARAMETER :: Include_NES  = .FALSE.
  LOGICAL, PARAMETER :: Include_Pair = .FALSE.

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

  REAL(DP), PARAMETER :: WFactor_FP = FourPi / PlanckConstant**3

  INTEGER  :: nE_G, nX_G, nZ(4)
  INTEGER  :: n_FP_inner, n_FP_outer
  INTEGER  :: iE_B, iE_E

  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:)
  REAL(DP), ALLOCATABLE :: W3_S(:)

  INTEGER,  ALLOCATABLE :: INFO(:)

  INTEGER               :: LWORK_outer
  REAL(DP), ALLOCATABLE :: WORK_outer(:,:)
  REAL(DP), ALLOCATABLE :: TAU_outer (:,:)
  REAL(DP), ALLOCATABLE :: BVEC_outer(:,:)
  REAL(DP), ALLOCATABLE :: AMAT_outer(:,:,:)

  INTEGER               :: LWORK_inner
  REAL(DP), ALLOCATABLE :: WORK_inner(:,:)
  REAL(DP), ALLOCATABLE :: TAU_inner (:,:)
  REAL(DP), ALLOCATABLE :: BVEC_inner(:,:)
  REAL(DP), ALLOCATABLE :: AMAT_inner(:,:,:)

  ! --- Solver Parameters to be initialized

  INTEGER  :: M_FP, M_outer, M_inner
  INTEGER  :: MaxIter_outer, MaxIter_inner
  REAL(DP) :: Rtol_outer, Rtol_inner

  INTEGER :: iS_1 = iNuE
  INTEGER :: iS_2 = iNuE_Bar

  INTEGER, PARAMETER :: iY  = 1
  INTEGER, PARAMETER :: iEf = 2
  INTEGER, PARAMETER :: iV1 = 3
  INTEGER, PARAMETER :: iV2 = 4
  INTEGER, PARAMETER :: iV3 = 5

  ! --- Temporary arrays for scatter/gather (packing)

  REAL(DP), ALLOCATABLE, TARGET :: P1D(:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D(:,:,:,:)

  INTEGER, PARAMETER :: iP1D_D = 1
  INTEGER, PARAMETER :: iP1D_T = 2
  INTEGER, PARAMETER :: iP1D_Y = 3
  INTEGER, PARAMETER :: iP1D_E = 4
  INTEGER, PARAMETER :: nP1D   = 4

  INTEGER, PARAMETER :: iP2D_Chi_NES_1  = 1
  INTEGER, PARAMETER :: iP2D_Chi_NES_2  = 2
  INTEGER, PARAMETER :: iP2D_Eta_NES_1  = 3
  INTEGER, PARAMETER :: iP2D_Eta_NES_2  = 4
  INTEGER, PARAMETER :: iP2D_Chi_Pair_1 = 5
  INTEGER, PARAMETER :: iP2D_Chi_Pair_2 = 6
  INTEGER, PARAMETER :: iP2D_Eta_Pair_1 = 7
  INTEGER, PARAMETER :: iP2D_Eta_Pair_2 = 8
  INTEGER, PARAMETER :: iP2D_Chi_1      = 9
  INTEGER, PARAMETER :: iP2D_Chi_2      = 10
  INTEGER, PARAMETER :: iP2D_J0_1       = 11
  INTEGER, PARAMETER :: iP2D_J0_2       = 12
  INTEGER, PARAMETER :: iP2D_J_1        = 13
  INTEGER, PARAMETER :: iP2D_J_2        = 14
  INTEGER, PARAMETER :: iP2D_Sig_1      = 15
  INTEGER, PARAMETER :: iP2D_Sig_2      = 16
  INTEGER, PARAMETER :: nP2D            = 16

  INTEGER, PARAMETER :: iP3D_Phi_0_In_NES_1  = 1
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_NES_1  = 2
  INTEGER, PARAMETER :: iP3D_Phi_0_In_NES_2  = 3
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_NES_2  = 4
  INTEGER, PARAMETER :: iP3D_Phi_0_In_Pair_1 = 5
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_Pair_1 = 6
  INTEGER, PARAMETER :: iP3D_Phi_0_In_Pair_2 = 7
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_Pair_2 = 8
  INTEGER, PARAMETER :: iP3D_WORK1           = 9
  INTEGER, PARAMETER :: iP3D_WORK2           = 10
  INTEGER, PARAMETER :: nP3D                 = 10

CONTAINS


  SUBROUTINE InitializeNeutrinoMatterSolver( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    REAL(DP) :: TMP(1)

    iE_B = iZ_B(1)
    iE_E = iZ_E(1)
    nZ   = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    n_FP_outer = 5
    n_FP_inner = nE_G * nCR * nSpecies

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S(:) = WFactor_FP * W2_N(:)
    W3_S(:) = WFactor_FP * W3_N(:)

    ALLOCATE( INFO(nX_G) )

    ALLOCATE( P1D(          nX_G,nP1D) )
    ALLOCATE( P2D(nE_G     ,nX_G,nP2D) )
    ALLOCATE( P3D(nE_G,nE_G,nX_G,nP3D) )

    ALLOCATE( TAU_outer (n_FP_outer,        nX_G) )
    ALLOCATE( BVEC_outer(n_FP_outer,        nX_G) )
    ALLOCATE( AMAT_outer(n_FP_outer,M_outer,nX_G) )

    ALLOCATE( TAU_inner (n_FP_inner,        nX_G) )
    ALLOCATE( BVEC_inner(n_FP_inner,        nX_G) )
    ALLOCATE( AMAT_inner(n_FP_inner,M_inner,nX_G) )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$OMP MAP( alloc: INFO, P1D, P2D, P3D, &
    !$OMP             AMAT_outer, BVEC_outer, TAU_outer, &
    !$OMP             AMAT_inner, BVEC_inner, TAU_inner )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$ACC CREATE( INFO, P1D, P2D, P3D, &
    !$ACC         AMAT_outer, BVEC_outer, TAU_outer, &
    !$ACC         AMAT_inner, BVEC_inner, TAU_inner )
#endif

    IF( M_outer > 3 )THEN

      CALL LinearLeastSquares_LWORK &
             ( 'N', n_FP_outer, M_outer-1, 1, AMAT_outer, n_FP_outer, &
               BVEC_outer, n_FP_outer, TMP, LWORK_outer )

    ELSE

      LWORK_outer = 1

    END IF

    IF( M_inner > 3 )THEN

      CALL LinearLeastSquares_LWORK &
             ( 'N', n_FP_inner, M_inner-1, 1, AMAT_inner, n_FP_inner, &
               BVEC_inner, n_FP_inner, TMP, LWORK_inner )

    ELSE

      LWORK_inner = 1

    END IF

    ALLOCATE( WORK_outer(LWORK_outer,nX_G) )
    ALLOCATE( WORK_inner(LWORK_inner,nX_G) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: WORK_outer, WORK_inner )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( WORK_outer, WORK_inner )
#endif

  END SUBROUTINE InitializeNeutrinoMatterSolver


  SUBROUTINE InitializeNeutrinoMatterSolverParameters &
    ( M_outer_Option, M_inner_Option, MaxIter_outer_Option, &
      MaxIter_inner_Option, Rtol_inner_Option, Rtol_outer_Option )

    INTEGER , INTENT(in), OPTIONAL :: M_outer_Option
    INTEGER , INTENT(in), OPTIONAL :: M_inner_Option
    INTEGER , INTENT(in), OPTIONAL :: MaxIter_outer_Option
    INTEGER , INTENT(in), OPTIONAL :: MaxIter_inner_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rtol_inner_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rtol_outer_Option

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

    IF( PRESENT( Rtol_inner_Option ) )THEN
      Rtol_inner = Rtol_inner_Option
    ELSE
      Rtol_inner = 1.0d-08
    END IF

    IF( PRESENT( Rtol_outer_Option ) )THEN
      Rtol_outer = Rtol_outer_Option
    ELSE
      Rtol_outer = 1.0d-08
    END IF

    M_FP = MAX( M_outer, M_inner )

  END SUBROUTINE InitializeNeutrinoMatterSolverParameters


  SUBROUTINE FinalizeNeutrinoMatterSolver

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, &
    !$OMP               INFO, P1D, P2D, P3D, &
    !$OMP               AMAT_outer, BVEC_outer, TAU_outer, WORK_outer, &
    !$OMP               AMAT_inner, BVEC_inner, TAU_inner, WORK_inner )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, &
    !$ACC         INFO, P1D, P2D, P3D, &
    !$ACC         AMAT_outer, BVEC_outer, TAU_outer, WORK_outer, &
    !$ACC         AMAT_inner, BVEC_inner, TAU_inner, WORK_inner )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S )
    DEALLOCATE( INFO )
    DEALLOCATE( TAU_outer, BVEC_outer, AMAT_outer, WORK_outer )
    DEALLOCATE( TAU_inner, BVEC_inner, AMAT_inner, WORK_inner )
    DEALLOCATE( P1D, P2D, P3D )

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

      W2(iN_E) = WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**2
      W3(iN_E) = WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**3

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA &
    ( dt, J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, D, T, Y, E, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Inner, nIterations_Outer )

    REAL(DP),                                      INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(inout) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(inout) :: D, T, Y, E
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER,  DIMENSION(1:nX_G),                   INTENT(out)   :: nIterations_Inner, nIterations_Outer

    ! --- Local Variables ---

    INTEGER  :: iN_E, iN_X, iS
    INTEGER  :: k_outer, Mk_outer, nX_P_outer
    INTEGER  :: k_inner, Mk_inner, nX_P_inner

    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE_OUTER, ITERATE_INNER
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex_outer, UnpackIndex_outer
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex_inner, UnpackIndex_inner

    REAL(DP), DIMENSION(1:nX_G) :: Omega, NormVsq

    REAL(DP), DIMENSION(1:nSpecies,1:nX_G) :: Jnorm

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: H_d_1, H_d_2, H_d_3
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: J_old, H_d_1_old, H_d_2_old, H_d_3_old
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: J_new, H_d_1_new, H_d_2_new, H_d_3_new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: C_J, C_H_d_1, C_H_d_2, C_H_d_3

    REAL(DP), DIMENSION(1:nX_G) :: Ef, V_d_1, V_d_2, V_d_3
    REAL(DP), DIMENSION(1:nX_G) :: Y_old, Ef_old, V_d_1_old, V_d_2_old, V_d_3_old
    REAL(DP), DIMENSION(1:nX_G) :: S_Y, S_ef, S_V_d_1, S_V_d_2, S_V_d_3
    REAL(DP), DIMENSION(1:nX_G) :: C_Y, C_ef, C_V_d_1, C_V_d_2, C_V_d_3
    REAL(DP), DIMENSION(1:nX_G) :: U_Y, U_ef, U_V_d_1, U_V_d_2, U_V_d_3

    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: J0, Chi, Eta
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair

    REAL(DP), DIMENSION(1:n_FP_outer,1:M_outer,1:nX_G) :: GVEC_outer, FVEC_outer
    REAL(DP), DIMENSION(1:n_FP_outer,          1:nX_G) :: GVECm_outer, FVECm_outer
    REAL(DP), DIMENSION(             1:M_outer,1:nX_G) :: Alpha_outer

    REAL(DP), DIMENSION(1:n_FP_inner,1:M_inner,1:nX_G) :: GVEC_inner, FVEC_inner
    REAL(DP), DIMENSION(1:n_FP_inner,          1:nX_G) :: GVECm_inner, FVECm_inner
    REAL(DP), DIMENSION(             1:M_inner,1:nX_G) :: Alpha_inner

    ITERATE_OUTER(:) = .TRUE.
    ITERATE_INNER(:) = .TRUE.

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE_OUTER, ITERATE_INNER ) &
    !$OMP MAP( alloc: Omega, NormVsq, Jnorm, &
    !$OMP             H_d_1, H_d_2, H_d_3, &
    !$OMP             J_old, H_d_1_old, H_d_2_old, H_d_3_old, &
    !$OMP             J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
    !$OMP             C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
    !$OMP             Ef, V_d_1, V_d_2, V_d_3, &
    !$OMP             Y_old, Ef_old, V_d_1_old, V_d_2_old, V_d_3_old, &
    !$OMP             S_Y, S_ef, S_V_d_1, S_V_d_2, S_V_d_3, &
    !$OMP             C_Y, C_ef, C_V_d_1, C_V_d_2, C_V_d_3, &
    !$OMP             U_Y, U_ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$OMP             J0, Chi, Eta, &
    !$OMP             Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$OMP             Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$OMP             PackIndex_outer, UnpackIndex_outer, GVEC_outer, FVEC_outer, &
    !$OMP             GVECm_outer, FVECm_outer, Alpha_outer, &
    !$OMP             PackIndex_inner, UnpackIndex_inner, GVEC_inner, FVEC_inner, &
    !$OMP             GVECm_inner, FVECm_inner, Alpha_inner )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE_OUTER, ITERATE_INNER ) &
    !$ACC CREATE( Omega, NormVsq, Jnorm, &
    !$ACC         H_d_1, H_d_2, H_d_3, &
    !$ACC         J_old, H_d_1_old, H_d_2_old, H_d_3_old, &
    !$ACC         J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
    !$ACC         C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
    !$ACC         Ef, V_d_1, V_d_2, V_d_3, &
    !$ACC         Y_old, Ef_old, V_d_1_old, V_d_2_old, V_d_3_old, &
    !$ACC         S_Y, S_ef, S_V_d_1, S_V_d_2, S_V_d_3, &
    !$ACC         C_Y, C_ef, C_V_d_1, C_V_d_2, C_V_d_3, &
    !$ACC         U_Y, U_ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$ACC         J0, Chi, Eta, &
    !$ACC         Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$ACC         Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, GVEC_outer, FVEC_outer, &
    !$ACC         GVECm_outer, FVECm_outer, Alpha_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, GVEC_inner, FVEC_inner, &
    !$ACC         GVECm_inner, FVECm_inner, Alpha_inner )
#endif

    ! --- Initial RHS ---

    CALL TimersStart( Timer_Collisions_InitializeRHS )

    CALL InitializeRHS_FP &
           ( J, H_u_1, H_u_2, H_u_3, &
             H_d_1, H_d_2, H_d_3, &
             C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
             D, Y, E, V_u_1, V_u_2, V_u_3, &
             Ef, V_d_1, V_d_2, V_d_3, &
             C_Y, C_Ef, C_V_d_1, C_V_d_2, C_V_d_3, &
             S_Y, S_Ef, S_V_d_1, S_V_d_2, S_V_d_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    CALL TimersStop( Timer_Collisions_InitializeRHS )

    ! --- Store Initial Matter State ---

    CALL ArrayCopy &
           ( Y    , Ef    , V_d_1    , V_d_2    , V_d_3, &
             Y_old, Ef_old, V_d_1_old, V_d_2_old, V_d_3_old )

    ! --- Store Initial Neutrino State ---

    CALL ArrayCopy &
           ( J    , H_d_1    , H_d_2    , H_d_3, &
             J_old, H_d_1_old, H_d_2_old, H_d_3_old )

    ! --- Initial Guess for Neutrino State ---

    CALL ArrayCopy &
           ( J    , H_d_1    , H_d_2    , H_d_3, &
             J_new, H_d_1_new, H_d_2_new, H_d_3_new )

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Collisions_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( D, T, Y, J0, Chi, Sig, Phi_0_In_NES, Phi_0_Ot_NES, &
             Phi_0_In_Pair, Phi_0_Ot_Pair )

    CALL TimersStop( Timer_Collisions_ComputeOpacity )

    ! --- Start Outer Loop ---

    k_outer = 0
    DO WHILE( ANY( ITERATE_OUTER(:) ) .AND. k_outer < MaxIter_outer )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      CALL ComputeJNorm &
             ( ITERATE_OUTER, J_new, Jnorm )

      CALL CreatePackIndex &
             ( ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      IF ( k_outer > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Collisions_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( D, T, Y, J0, Chi, Sig, Phi_0_In_NES, Phi_0_Ot_NES, &
                 Phi_0_In_Pair, Phi_0_Ot_Pair, &
                 ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        CALL TimersStop( Timer_Collisions_ComputeOpacity )

      END IF

      ! --- Compute Omega (Richardson damping coeff.) based on velocity ---
      ! --- For first interation: must be consistent with initial guess ---

      DO iN_X = 1, nX_G

        Omega(iN_X) = One / ( One + SQRT( NormVsq(iN_X) ) )

      END DO

      ! --- Start Inner Loop ---

      CALL TimersStart( Timer_Collisions_InnerLoop )

      k_inner = 0
      DO WHILE( ANY( ITERATE_INNER(:) ) .AND. k_inner < MaxIter_inner )

        k_inner  = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

        CALL TimersStart( Timer_Collisions_ComputeRates )

        CALL ComputeRates_Packed &
               ( J_new, &
                 Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, Phi_0_Ot_Pair, &
                 Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
                 ITERATE_INNER, nX_P_inner, PackIndex_inner, &
                 UnpackIndex_inner, nX_P_outer )

        CALL TimersStop( Timer_Collisions_ComputeRates )

        ! --- Right-Hand Side Vectors and Residuals (inner) ---

        CALL TimersStart( Timer_Collisions_NeutrinoRHS )

        CALL ComputeNeutrinoRHS_FP &
               ( ITERATE_INNER, n_FP_inner, FVECm_inner, GVECm_inner, &
                 dt, Omega, V_u_1, V_u_2, V_u_3, &
                 C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
                 J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
                 J0, Chi, Sig, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        CALL TimersStop( Timer_Collisions_NeutrinoRHS )

        ! --- Anderson Acceleration (inner) ---

        CALL TimersStart( Timer_Collisions_SolveLS )

        CALL SolveLS_FP &
               ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                 FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner, &
                 AMAT_inner, BVEC_inner, Alpha_inner, TAU_inner, &
                 LWORK_inner, WORK_inner )

        CALL TimersStop( Timer_Collisions_SolveLS )

        ! --- Update Residuals and Solution Vectors (inner) ---

        CALL TimersStart( Timer_Collisions_UpdateFP )

        CALL UpdateNeutrinoRHS_FP &
               ( ITERATE_INNER, n_FP_inner, FVECm_inner, GVECm_inner, &
                 J_new, H_d_1_new, H_d_2_new, H_d_3_new )

        ! --- Check Convergence (inner) ---

        CALL TimersStart( Timer_Collisions_CheckInner )

        CALL CheckConvergence_Inner &
               ( ITERATE_INNER, n_FP_inner, k_inner, &
                 nIterations_Inner, FVECm_inner, Jnorm )

        CALL TimersStop( Timer_Collisions_CheckInner )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Collisions_UpdateFP )

      END DO ! --- Inner Loop ---

      CALL TimersStop( Timer_Collisions_InnerLoop )

      ! --- Right-Hand Side Vectors and Residuals (outer) ---

      CALL TimersStart( Timer_Collisions_MatterRHS )

      CALL ComputeMatterRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, FVECm_outer, GVECm_outer, &
               J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
                      C_Y    , S_Y    , U_Y    , &
                      C_Ef   , S_Ef   , U_Ef   , &
               V_d_1, C_V_d_1, S_V_d_1, U_V_d_1, &
               V_d_2, C_V_d_2, S_V_d_2, U_V_d_2, &
               V_d_3, C_V_d_3, S_V_d_3, U_V_d_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL TimersStop( Timer_Collisions_MatterRHS )

      ! --- Anderson Acceleration (outer) ---

      CALL TimersStart( Timer_Collisions_SolveLS )

      CALL SolveLS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
               FVECm_outer, GVECm_outer, FVEC_outer, GVEC_outer, &
               AMAT_outer, BVEC_outer, Alpha_outer, TAU_outer, &
               LWORK_outer, WORK_outer )

      CALL TimersStop( Timer_Collisions_SolveLS )

      ! --- Update Residuals and Solution Vectors (outer) ---

      CALL TimersStart( Timer_Collisions_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, &
               Y    , Y_old , U_Y , &
               Ef   , Ef_old, U_Ef, &
               V_d_1, V_d_1_old, U_V_d_1, &
               V_d_2, V_d_2_old, U_V_d_2, &
               V_d_3, V_d_3_old, U_V_d_3, &
               FVECm_outer, GVECm_outer )

      ! --- Update Three-Velocity (Index Up) ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        V_u_1(iN_X) = V_d_1(iN_X) / Gm_dd_11(iN_X)
        V_u_2(iN_X) = V_d_2(iN_X) / Gm_dd_22(iN_X)
        V_u_3(iN_X) = V_d_3(iN_X) / Gm_dd_33(iN_X)

      END DO

      ! --- Compute E from Ef ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G

        NormVsq(iN_X) =   V_u_1(iN_X) * V_d_1(iN_X) &
                        + V_u_2(iN_X) * V_d_2(iN_X) &
                        + V_u_3(iN_X) * V_d_3(iN_X)

        E(iN_X) = Ef(iN_X) - Half * NormVsq(iN_X)

      END DO

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      ! --- Check Convergence (outer) ---

      CALL TimersStart( Timer_Collisions_CheckOuter )

      CALL CheckConvergence_Outer &
             ( ITERATE_OUTER, ITERATE_INNER, n_FP_outer, k_outer, &
               FVECm_outer, nIterations_Outer )

      CALL TimersStop( Timer_Collisions_CheckOuter )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
               FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Collisions_UpdateFP )

    END DO ! --- Outer Loop ---

    nIterations_Inner &
      = FLOOR( DBLE( nIterations_Inner ) / DBLE( nIterations_Outer ) )

    CALL ArrayCopy &
           ( J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
             J    , H_d_1    , H_d_2    , H_d_3 )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      H_u_1(iN_E,iN_X,iS) = H_d_1(iN_E,iN_X,iS) / Gm_dd_11(iN_X)
      H_u_2(iN_E,iN_X,iS) = H_d_2(iN_E,iN_X,iS) / Gm_dd_22(iN_X)
      H_u_3(iN_E,iN_X,iS) = H_d_3(iN_E,iN_X,iS) / Gm_dd_33(iN_X)

    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE_OUTER, ITERATE_INNER, &
    !$OMP               Omega, NormVsq, Jnorm, &
    !$OMP               H_d_1, H_d_2, H_d_3, &
    !$OMP               J_old, H_d_1_old, H_d_2_old, H_d_3_old, &
    !$OMP               J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
    !$OMP               C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
    !$OMP               Ef, V_d_1, V_d_2, V_d_3, &
    !$OMP               Y_old, Ef_old, V_d_1_old, V_d_2_old, V_d_3_old, &
    !$OMP               S_Y, S_ef, S_V_d_1, S_V_d_2, S_V_d_3, &
    !$OMP               C_Y, C_ef, C_V_d_1, C_V_d_2, C_V_d_3, &
    !$OMP               U_Y, U_ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$OMP               J0, Chi, Eta, &
    !$OMP               Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$OMP               Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$OMP               PackIndex_outer, UnpackIndex_outer, GVEC_outer, FVEC_outer, &
    !$OMP               GVECm_outer, FVECm_outer, Alpha_outer, &
    !$OMP               PackIndex_inner, UnpackIndex_inner, GVEC_inner, FVEC_inner, &
    !$OMP               GVECm_inner, FVECm_inner, Alpha_inner )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE_OUTER, ITERATE_INNER, &
    !$ACC         Omega, NormVsq, Jnorm, &
    !$ACC         H_d_1, H_d_2, H_d_3, &
    !$ACC         J_old, H_d_1_old, H_d_2_old, H_d_3_old, &
    !$ACC         J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
    !$ACC         C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
    !$ACC         Ef, V_d_1, V_d_2, V_d_3, &
    !$ACC         Y_old, Ef_old, V_d_1_old, V_d_2_old, V_d_3_old, &
    !$ACC         S_Y, S_ef, S_V_d_1, S_V_d_2, S_V_d_3, &
    !$ACC         C_Y, C_ef, C_V_d_1, C_V_d_2, C_V_d_3, &
    !$ACC         U_Y, U_ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$ACC         J0, Chi, Eta, &
    !$ACC         Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$ACC         Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, GVEC_outer, FVEC_outer, &
    !$ACC         GVECm_outer, FVECm_outer, Alpha_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, GVEC_inner, FVEC_inner, &
    !$ACC         GVECm_inner, FVECm_inner, Alpha_inner )
#endif

  END SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA


  SUBROUTINE ComputeOpacities_Packed &
    ( D, T, Y, J0, Chi, Sig, Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, &
      Phi_0_Ot_Pair, MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), INTENT(in)   , TARGET   :: D(1:nX_G)
    REAL(DP), INTENT(in)   , TARGET   :: T(1:nX_G)
    REAL(DP), INTENT(in)   , TARGET   :: Y(1:nX_G)
    REAL(DP), INTENT(inout), TARGET   :: J0 (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: Chi(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: Sig(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: &
      Phi_0_In_NES (1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: &
      Phi_0_Ot_NES (1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: &
      Phi_0_In_Pair(1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: &
      Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    LOGICAL,  INTENT(in)   , OPTIONAL :: MASK(1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: nX_P
    INTEGER,  INTENT(in)   , OPTIONAL :: PackIndex  (1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: UnpackIndex(1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: nX_P0

    REAL(DP), DIMENSION(:),     POINTER :: D_P, T_P, Y_P
    ! --- to be changed to handle cases when iSpecies > 2 ---
    REAL(DP), DIMENSION(:,:),   POINTER :: J0_1_P, J0_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Chi_1_P, Chi_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Sig_1_P, Sig_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_1_P,  Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_2_P,  Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P
    ! --- to be changed to handle cases when iSpecies > 2 ---

    INTEGER :: nX, nX0, iX, iE1, iE2

    IF( PRESENT( nX_P ) )THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF( PRESENT( nX_P0 ) )THEN
      nX0 = nX_P0
    ELSE
      nX0 = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P => P1D(1:nX,iP1D_D)
      T_P => P1D(1:nX,iP1D_T)
      Y_P => P1D(1:nX,iP1D_Y)

      CALL ArrayPack( nX, UnpackIndex, D, T, Y, D_P, T_P, Y_P )

      J0_1_P  => P2D(:,1:nX,iP2D_J0_1)
      J0_2_P  => P2D(:,1:nX,iP2D_J0_2)
      Chi_1_P => P2D(:,1:nX,iP2D_Chi_1)
      Chi_2_P => P2D(:,1:nX,iP2D_Chi_2)
      Sig_1_P => P2D(:,1:nX,iP2D_Sig_1)
      Sig_2_P => P2D(:,1:nX,iP2D_Sig_2)

      Phi_0_In_NES_1_P  => P3D(:,:,1:nX,iP3D_Phi_0_In_NES_1)
      Phi_0_Ot_NES_1_P  => P3D(:,:,1:nX,iP3D_Phi_0_Ot_NES_1)
      Phi_0_In_NES_2_P  => P3D(:,:,1:nX,iP3D_Phi_0_In_NES_2)
      Phi_0_Ot_NES_2_P  => P3D(:,:,1:nX,iP3D_Phi_0_Ot_NES_2)
      Phi_0_In_Pair_1_P => P3D(:,:,1:nX,iP3D_Phi_0_In_Pair_1)
      Phi_0_Ot_Pair_1_P => P3D(:,:,1:nX,iP3D_Phi_0_Ot_Pair_1)
      Phi_0_In_Pair_2_P => P3D(:,:,1:nX,iP3D_Phi_0_In_Pair_2)
      Phi_0_Ot_Pair_2_P => P3D(:,:,1:nX,iP3D_Phi_0_Ot_Pair_2)

    ELSE

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      J0_1_P  => J0 (:,:,iS_1)
      J0_2_P  => J0 (:,:,iS_2)
      Chi_1_P => Chi(:,:,iS_1)
      Chi_2_P => Chi(:,:,iS_2)
      Sig_1_P => Sig(:,:,iS_1)
      Sig_2_P => Sig(:,:,iS_2)

      Phi_0_In_NES_1_P  => Phi_0_In_NES (:,:,:,iS_1)
      Phi_0_Ot_NES_1_P  => Phi_0_Ot_NES (:,:,:,iS_1)
      Phi_0_In_NES_2_P  => Phi_0_In_NES (:,:,:,iS_2)
      Phi_0_Ot_NES_2_P  => Phi_0_Ot_NES (:,:,:,iS_2)
      Phi_0_In_Pair_1_P => Phi_0_In_Pair(:,:,:,iS_1)
      Phi_0_Ot_Pair_1_P => Phi_0_Ot_Pair(:,:,:,iS_1)
      Phi_0_In_Pair_2_P => Phi_0_In_Pair(:,:,:,iS_2)
      Phi_0_Ot_Pair_2_P => Phi_0_Ot_Pair(:,:,:,iS_2)

    END IF

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_DG_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, J0_2_P, iS_1, iS_2 )

    ! --- EmAb ---

    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, Chi_1_P )

    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_2, Chi_2_P )

    ! --- Isoenergetic scattering ---

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, 1, Sig_1_P )

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_2, 1, Sig_2_P )

    IF( Include_NES )THEN

      ! --- NES Kernels ---

      CALL ComputeNeutrinoOpacities_NES_Points &
             ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, iS_2, 1, &
               Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, &
               Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
               P3D(:,:,:,iP3D_WORK1), P3D(:,:,:,iP3D_WORK2) )

      ! --- Enforce Detailed Balance (Again) Based on J0 ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX  = 1, nX
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        IF( iE1 <= iE2 )THEN

          Phi_0_In_NES_1_P(iE1,iE2,iX) &
            = Phi_0_Ot_NES_1_P(iE1,iE2,iX) &
              * ( J0_1_P(iE2,iX) * ( One - J0_1_P(iE1,iX) ) ) &
              / ( J0_1_P(iE1,iX) * ( One - J0_1_P(iE2,iX) ) )

          Phi_0_In_NES_2_P(iE1,iE2,iX) &
            = Phi_0_Ot_NES_2_P(iE1,iE2,iX) &
              * ( J0_2_P(iE2,iX) * ( One - J0_2_P(iE1,iX) ) ) &
              / ( J0_2_P(iE1,iX) * ( One - J0_2_P(iE2,iX) ) )

        ELSE

          Phi_0_Ot_NES_1_P(iE1,iE2,iX) &
            = Phi_0_In_NES_1_P(iE1,iE2,iX) &
              * ( J0_1_P(iE1,iX) * ( One - J0_1_P(iE2,iX) ) ) &
              / ( J0_1_P(iE2,iX) * ( One - J0_1_P(iE1,iX) ) )

          Phi_0_Ot_NES_2_P(iE1,iE2,iX) &
            = Phi_0_In_NES_2_P(iE1,iE2,iX) &
              * ( J0_2_P(iE1,iX) * ( One - J0_2_P(iE2,iX) ) ) &
              / ( J0_2_P(iE2,iX) * ( One - J0_2_P(iE1,iX) ) )

        END IF

      END DO
      END DO
      END DO

    ELSE

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX  = 1, nX
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        Phi_0_In_NES_1_P(iE1,iE2,iX) = Zero
        Phi_0_Ot_NES_1_P(iE1,iE2,iX) = Zero
        Phi_0_In_NES_2_P(iE1,iE2,iX) = Zero
        Phi_0_Ot_NES_2_P(iE1,iE2,iX) = Zero

      END DO
      END DO
      END DO

    END IF

    IF( Include_Pair )THEN

      ! --- Pair Kernels ---

      CALL ComputeNeutrinoOpacities_Pair_Points &
             ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, iS_2, 1, &
               Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
               Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
               P3D(:,:,:,iP3D_WORK1), P3D(:,:,:,iP3D_WORK2) )

      ! --- Enforce Detailed Balance (Again) Based on J0 ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX  = 1, nX
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        Phi_0_In_Pair_1_P(iE1,iE2,iX) &
          = Phi_0_Ot_Pair_1_P(iE1,iE2,iX) &
            * (         J0_1_P(iE2,iX)   *         J0_2_P(iE1,iX)   ) &
            / ( ( One - J0_1_P(iE2,iX) ) * ( One - J0_2_P(iE1,iX) ) )

        Phi_0_In_Pair_2_P(iE1,iE2,iX) &
          = Phi_0_Ot_Pair_2_P(iE1,iE2,iX) &
            * (         J0_2_P(iE2,iX)   *         J0_1_P(iE1,iX)   ) &
            / ( ( One - J0_2_P(iE2,iX) ) * ( One - J0_1_P(iE1,iX) ) )

      END DO
      END DO
      END DO

    ELSE

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iX  = 1, nX
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        Phi_0_In_Pair_1_P(iE1,iE2,iX) = Zero
        Phi_0_Ot_Pair_1_P(iE1,iE2,iX) = Zero
        Phi_0_In_Pair_2_P(iE1,iE2,iX) = Zero
        Phi_0_Ot_Pair_2_P(iE1,iE2,iX) = Zero

      END DO
      END DO
      END DO

    END IF

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               J0_1_P, J0_2_P, J0(:,:,iS_1), J0(:,:,iS_2) )
      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_1_P, Chi_2_P, Sig_1_P, Sig_2_P, &
               Chi(:,:,iS_1), Chi(:,:,iS_2), Sig(:,:,iS_1), Sig(:,:,iS_2) )

      IF ( nX < nX0 ) THEN

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 Phi_0_In_NES_1_P , Phi_0_Ot_NES_1_P, &
                 Phi_0_In_NES_2_P , Phi_0_Ot_NES_2_P, &
                 Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
                 Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
                 Phi_0_In_NES (:,:,:,iS_1), Phi_0_Ot_NES (:,:,:,iS_1), &
                 Phi_0_In_NES (:,:,:,iS_2), Phi_0_Ot_NES (:,:,:,iS_2), &
                 Phi_0_In_Pair(:,:,:,iS_1), Phi_0_Ot_Pair(:,:,:,iS_1), &
                 Phi_0_In_Pair(:,:,:,iS_2), Phi_0_Ot_Pair(:,:,:,iS_2) )

      END IF

    END IF

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( J, Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, Phi_0_Ot_Pair, &
      Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, MASK, nX_P, PackIndex, &
      UnpackIndex, nX_P0 )

    REAL(DP), INTENT(in)   , TARGET   :: J(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)   , TARGET   :: &
      Phi_0_In_NES (1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)   , TARGET   :: &
      Phi_0_Ot_NES (1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)   , TARGET   :: &
      Phi_0_In_Pair(1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)   , TARGET   :: &
      Phi_0_Ot_Pair(1:nE_G,1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: Chi_NES (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: Eta_NES (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: Chi_Pair(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout), TARGET   :: Eta_Pair(1:nE_G,1:nX_G,1:nSpecies)
    LOGICAL,  INTENT(in)   , OPTIONAL :: MASK(1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: nX_P
    INTEGER,  INTENT(in)   , OPTIONAL :: PackIndex  (1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: UnpackIndex(1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: nX_P0

    REAL(DP), DIMENSION(:,:),   POINTER :: J_1_P, J_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Chi_NES_1_P, Chi_NES_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Eta_NES_1_P, Eta_NES_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Chi_Pair_1_P, Chi_Pair_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Eta_Pair_1_P, Eta_Pair_2_P

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

      Chi_NES_1_P  => P2D(:,1:nX,iP2D_Chi_NES_1)
      Chi_NES_2_P  => P2D(:,1:nX,iP2D_Chi_NES_2)
      Eta_NES_1_P  => P2D(:,1:nX,iP2D_Eta_NES_1)
      Eta_NES_2_P  => P2D(:,1:nX,iP2D_Eta_NES_2)
      Chi_Pair_1_P => P2D(:,1:nX,iP2D_Chi_Pair_1)
      Chi_Pair_2_P => P2D(:,1:nX,iP2D_Chi_Pair_2)
      Eta_Pair_1_P => P2D(:,1:nX,iP2D_Eta_Pair_1)
      Eta_Pair_2_P => P2D(:,1:nX,iP2D_Eta_Pair_2)

      J_1_P => P2D(:,1:nX,iP2D_J_1)
      J_2_P => P2D(:,1:nX,iP2D_J_2)

      CALL ArrayPack &
             ( nX, UnpackIndex, J(:,:,iS_1), J(:,:,iS_2), J_1_P, J_2_P )

      Phi_0_In_NES_1_P  => P3D(:,:,1:nX,iP3D_Phi_0_In_NES_1)
      Phi_0_Ot_NES_1_P  => P3D(:,:,1:nX,iP3D_Phi_0_Ot_NES_1)
      Phi_0_In_NES_2_P  => P3D(:,:,1:nX,iP3D_Phi_0_In_NES_2)
      Phi_0_Ot_NES_2_P  => P3D(:,:,1:nX,iP3D_Phi_0_Ot_NES_2)
      Phi_0_In_Pair_1_P => P3D(:,:,1:nX,iP3D_Phi_0_In_Pair_1)
      Phi_0_Ot_Pair_1_P => P3D(:,:,1:nX,iP3D_Phi_0_Ot_Pair_1)
      Phi_0_In_Pair_2_P => P3D(:,:,1:nX,iP3D_Phi_0_In_Pair_2)
      Phi_0_Ot_Pair_2_P => P3D(:,:,1:nX,iP3D_Phi_0_Ot_Pair_2)

      IF ( nX < nX0 ) THEN

        CALL ArrayPack &
               ( nX, UnpackIndex, &
                 Phi_0_In_NES (:,:,:,iS_1), Phi_0_Ot_NES (:,:,:,iS_1), &
                 Phi_0_In_NES (:,:,:,iS_2), Phi_0_Ot_NES (:,:,:,iS_2), &
                 Phi_0_In_Pair(:,:,:,iS_1), Phi_0_Ot_Pair(:,:,:,iS_1), &
                 Phi_0_In_Pair(:,:,:,iS_2), Phi_0_Ot_Pair(:,:,:,iS_2), &
                 Phi_0_In_NES_1_P , Phi_0_Ot_NES_1_P , &
                 Phi_0_In_NES_2_P , Phi_0_Ot_NES_2_P , &
                 Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
                 Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P )

      END IF

    ELSE

      Chi_NES_1_P  => Chi_NES (:,:,iS_1)
      Chi_NES_2_P  => Chi_NES (:,:,iS_2)
      Eta_NES_1_P  => Eta_NES (:,:,iS_1)
      Eta_NES_2_P  => Eta_NES (:,:,iS_2)
      Chi_Pair_1_P => Chi_Pair(:,:,iS_1)
      Chi_Pair_2_P => Chi_Pair(:,:,iS_2)
      Eta_Pair_1_P => Eta_Pair(:,:,iS_1)
      Eta_Pair_2_P => Eta_Pair(:,:,iS_2)

      J_1_P => J(:,:,iS_1)
      J_2_P => J(:,:,iS_2)

      Phi_0_In_NES_1_P  => Phi_0_In_NES (:,:,:,iS_1)
      Phi_0_Ot_NES_1_P  => Phi_0_Ot_NES (:,:,:,iS_1)
      Phi_0_In_NES_2_P  => Phi_0_In_NES (:,:,:,iS_2)
      Phi_0_Ot_NES_2_P  => Phi_0_Ot_NES (:,:,:,iS_2)
      Phi_0_In_Pair_1_P => Phi_0_In_Pair(:,:,:,iS_1)
      Phi_0_Ot_Pair_1_P => Phi_0_Ot_Pair(:,:,:,iS_1)
      Phi_0_In_Pair_2_P => Phi_0_In_Pair(:,:,:,iS_2)
      Phi_0_Ot_Pair_2_P => Phi_0_Ot_Pair(:,:,:,iS_2)

    END IF

    ! --- NES Emissivities and Opacities ---

    CALL ComputeNeutrinoOpacitiesRates_NES_Points &
           ( 1, nE_G, 1, nX, W2_N, J_1_P, &
             Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, &
             Eta_NES_1_P, Chi_NES_1_P )

    CALL ComputeNeutrinoOpacitiesRates_NES_Points &
           ( 1, nE_G, 1, nX, W2_N, J_2_P, &
             Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
             Eta_NES_2_P, Chi_NES_2_P )

    ! --- Pair Emissivities and Opacities ---

    CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
           ( 1, nE_G, 1, nX, W2_N, J_2_P, &
             Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
             Eta_Pair_1_P, Chi_Pair_1_P )

    CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
           ( 1, nE_G, 1, nX, W2_N, J_1_P, &
             Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
             Eta_Pair_2_P, Chi_Pair_2_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NES_1_P , Chi_NES_2_P , &
               Eta_NES_1_P , Eta_NES_2_P , &
               Chi_Pair_1_P, Chi_Pair_2_P, &
               Eta_Pair_1_P, Eta_Pair_2_P, &
               Chi_NES (:,:,iS_1), Chi_NES (:,:,iS_2), &
               Eta_NES (:,:,iS_1), Eta_NES (:,:,iS_2), &
               Chi_Pair(:,:,iS_1), Chi_Pair(:,:,iS_2), &
               Eta_Pair(:,:,iS_1), Eta_Pair(:,:,iS_2) )

    END IF

  END SUBROUTINE ComputeRates_Packed


  SUBROUTINE UpdateTemperature_Packed &
    ( D, E, Y, T, MASK, nX_P, PackIndex, UnpackIndex )

    REAL(DP), INTENT(in)   , TARGET   :: D(1:nX_G)
    REAL(DP), INTENT(in)   , TARGET   :: E(1:nX_G)
    REAL(DP), INTENT(in)   , TARGET   :: Y(1:nX_G)
    REAL(DP), INTENT(inout), TARGET   :: T(1:nX_G)
    LOGICAL,  INTENT(in)   , OPTIONAL :: MASK(1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: nX_P
    INTEGER,  INTENT(in)   , OPTIONAL :: PackIndex  (1:nX_G)
    INTEGER,  INTENT(in)   , OPTIONAL :: UnpackIndex(1:nX_G)

    INTEGER           :: nX
    REAL(DP), POINTER :: D_P(:), E_P(:), Y_P(:), T_P(:)

    IF( PRESENT( nX_P ) )THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P => P1D(1:nX,iP1D_D)
      Y_P => P1D(1:nX,iP1D_Y)
      T_P => P1D(1:nX,iP1D_T)
      E_P => P1D(1:nX,iP1D_E)

      CALL ArrayPack( nX, UnpackIndex, D, Y, T, E, D_P, Y_P, T_P, E_P )

    ELSE

      D_P => D(:)
      Y_P => Y(:)
      T_P => T(:)
      E_P => E(:)

    END IF

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D_P, E_P, Y_P, T_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack( nX, MASK, PackIndex, T_P, T )

    END IF

  END SUBROUTINE UpdateTemperature_Packed


  SUBROUTINE InitializeRHS_FP &
    ( J, H_u_1, H_u_2, H_u_3, &
      H_d_1, H_d_2, H_d_3, &
      C_J, C_H_d_1, C_H_d_2, C_H_d_3, &
      D, Y, E, V_u_1, V_u_2, V_u_3, &
      Ef, V_d_1, V_d_2, V_d_3, &
      C_Y, C_Ef, C_V_d_1, C_V_d_2, C_V_d_3, &
      S_Y, S_Ef, S_V_d_1, S_V_d_2, S_V_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)  :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(out) :: H_d_1, H_d_2, H_d_3
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(out) :: C_J, C_H_d_1, C_H_2_d, C_H_3_d
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)  :: D, Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(out) :: Ef, V_d_1, V_d_2, V_d_3
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(out) :: C_Y, C_Ef, C_V_d_1, C_V_d_2, C_V_d_3
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(out) :: S_Y, S_Ef, S_V_d_1, S_V_d_2, S_V_d_3
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC CREATE( N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 )
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( V_u_1, V_u_2, V_u_3, NormVsq, E, &
    !$ACC          Y, Ef, V_d_1, V_d_2, V_d_3, &
    !$ACC          U_Y, U_Ef, U_V_d_1, U_V_d_2, U_V_d_3, &
    !$ACC          S_Y, S_Ef, S_V_d_1, S_V_d_2, S_V_d_3, &
    !$ACC          Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      V_d_1(iN_X) = Gm_dd_11(iN_X) * V_u_1(iN_X)
      V_d_2(iN_X) = Gm_dd_22(iN_X) * V_u_2(iN_X)
      V_d_3(iN_X) = Gm_dd_33(iN_X) * V_u_3(iN_X)

      ! --- Specific Fluid Energy ---

      NormVsq(iN_X) &
        =    V_u_1(iN_X) * V_d_1(iN_X) &
           + V_u_2(iN_X) * V_d_2(iN_X) &
           + V_u_3(iN_X) * V_d_3(iN_X)

      Ef(iN_X) = E(iN_X) + Half * NormVsq(iN_X)

      ! --- Scaling Factors ---

      S_Y    (iN_X) = One / ( D(iN_X) * Y (iN_X) / AtomicMassUnit )
      S_Ef   (iN_X) = One / ( D(iN_X) * Ef(iN_X) )
      S_V_d_1(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_2(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_3(iN_X) = One / ( D(iN_X) * SpeedOfLight )

      ! --- Initial Guess for Matter State ---

      U_Y    (iN_X) = One
      U_Ef   (iN_X) = One
      U_V_d_1(iN_X) = V_d_1(iN_X) / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2(iN_X) / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3(iN_X) / SpeedOfLight

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_Y    (iN_X) = Zero
      C_Ef   (iN_X) = Zero
      C_V_d_1(iN_X) = Zero
      C_V_d_2(iN_X) = Zero
      C_V_d_3(iN_X) = Zero

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      H_d_1(iN_E,iN_X,iS) = Gm_dd_11(iN_X) * H_u_1(iN_E,iN_X,iS)
      H_d_2(iN_E,iN_X,iS) = Gm_dd_22(iN_X) * H_u_2(iN_E,iN_X,iS)
      H_d_3(iN_E,iN_X,iS) = Gm_dd_33(iN_X) * H_u_3(iN_E,iN_X,iS)

      vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iN_X,iS) &
              + V_u_2(iN_X) * H_d_2(iN_E,iN_X,iS) &
              + V_u_3(iN_X) * H_d_3(iN_E,iN_X,iS)

      k_dd = EddingtonTensorComponents_dd &
               ( J    (iN_E,iN_X,iS), H_u_1(iN_E,iN_X,iS), &
                 H_u_2(iN_E,iN_X,iS), H_u_3(iN_E,iN_X,iS), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

      vDotK_d_1 &
        = ( V_u_1(iN_X) * k_dd(1,1) &
          + V_u_2(iN_X) * k_dd(2,1) &
          + V_u_3(iN_X) * k_dd(3,1) ) * J(iN_E,iN_X,iS)
      vDotK_d_2 &                                                    
        = ( V_u_1(iN_X) * k_dd(1,2) &
          + V_u_2(iN_X) * k_dd(2,2) &
          + V_u_3(iN_X) * k_dd(3,2) ) * J(iN_E,iN_X,iS)
      vDotK_d_3 &                                                    
        = ( V_u_1(iN_X) * k_dd(1,3) &
          + V_u_2(iN_X) * k_dd(2,3) &
          + V_u_3(iN_X) * k_dd(3,3) ) * J(iN_E,iN_X,iS)

      ! --- Eulerian Neutrino Number Density ---

      N_nu(iN_E,iN_X,iS) = J(iN_E,iN_X,iS) + vDotH

      ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

      E_nu(iN_E,iN_X,iS) = J(iN_E,iN_X,iS) + Two * vDotH

      ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

      F_nu_d_1(iN_E,iN_X,iS) &
        = H_d_1(iN_E,iN_X,iS) + V_d_1(iN_X) * J(iN_E,iN_X,iS) + vDotK_d_1

      F_nu_d_2(iN_E,iN_X,iS) &
        = H_d_2(iN_E,iN_X,iS) + V_d_2(iN_X) * J(iN_E,iN_X,iS) + vDotK_d_2

      F_nu_d_3(iN_E,iN_X,iS) &
        = H_d_3(iN_E,iN_X,iS) + V_d_3(iN_X) * J(iN_E,iN_X,iS) + vDotK_d_3

      ! --- Old States for Neutrino Number Density and Flux ---

      C_J    (iN_E,iN_X,iS) = J    (iN_E,iN_X,iS) + vDotH
      C_H_d_1(iN_E,iN_X,iS) = H_d_1(iN_E,iN_X,iS) + vDotK_d_1
      C_H_d_2(iN_E,iN_X,iS) = H_d_2(iN_E,iN_X,iS) + vDotK_d_2
      C_H_d_3(iN_E,iN_X,iS) = H_d_3(iN_E,iN_X,iS) + vDotK_d_3

    END DO
    END DO
    END DO

    DO iS = iNuE, iNuE_Bar

      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, LeptonNumber(iS), N_nu(:,:,iS), &
               nE_G, W2_S, 1, One, C_Y, 1 )

    END DO

    DO iS = 1, nSpecies

      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, E_nu    (:,:,iS), nE_G, W3_S, 1, &
               One, C_Ef   , 1 )
      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, F_nu_d_1(:,:,iS), nE_G, W3_S, 1, &
               One, C_V_d_1, 1 )
      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, F_nu_d_2(:,:,iS), nE_G, W3_S, 1, &
               One, C_V_d_2, 1 )
      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, F_nu_d_3(:,:,iS), nE_G, W3_S, 1, &
               One, C_V_d_3, 1 )

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      !SUM_Y  = Zero
      !DO iS = iNuE, iNuE_Bar
      !DO iN_E = 1, nE_G
      !  SUM_Y  = SUM_Y  + N_nu    (iN_E,iN_X,iS) * W2_S(iN_E) * LeptonNumber(iS)
      !END DO
      !END DO
      !C_Y    (iN_X) = SUM_Y

      !SUM_Ef = Zero
      !SUM_V1 = Zero
      !SUM_V2 = Zero
      !SUM_V3 = Zero
      !DO iS = 1, nSpecies
      !DO iN_E = 1, nE_G
      !  SUM_Ef = SUM_Ef + E_nu    (iN_E,iN_X,iS) * W3_S(iN_E)
      !  SUM_V1 = SUM_V1 + F_nu_d_1(iN_E,iN_X,iS) * W3_S(iN_E)
      !  SUM_V2 = SUM_V2 + F_nu_d_2(iN_E,iN_X,iS) * W3_S(iN_E)
      !  SUM_V3 = SUM_V3 + F_nu_d_3(iN_E,iN_X,iS) * W3_S(iN_E)
      !END DO
      !END DO

      !C_Ef   (iN_X) = SUM_Ef
      !C_V_d_1(iN_X) = SUM_V1
      !C_V_d_2(iN_X) = SUM_V2
      !C_V_d_3(iN_X) = SUM_V3

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_Y    (iN_X) &
        = U_Y    (iN_X) + C_Y    (iN_X) * S_Y    (iN_X)
      C_Ef   (iN_X) &
        = U_Ef   (iN_X) + C_Ef   (iN_X) * S_Ef   (iN_X)
      C_V_d_1(iN_X) &
        = U_V_d_1(iN_X) + C_V_d_1(iN_X) * S_V_d_1(iN_X)
      C_V_d_2(iN_X) &
        = U_V_d_2(iN_X) + C_V_d_2(iN_X) * S_V_d_2(iN_X)
      C_V_d_3(iN_X) &
        = U_V_d_3(iN_X) + C_V_d_3(iN_X) * S_V_d_3(iN_X)

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 )
#endif

  END SUBROUTINE InitializeRHS_FP


  SUBROUTINE ComputeJNorm( MASK, J, Jnorm )

    LOGICAL , INTENT(in)    :: MASK(1:nX_G)
    REAL(DP), INTENT(in)    :: J(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout) :: Jnorm(1:nSpecies,1:nX_G)

    INTEGER :: iN_X, iS

    ! --- Fix me: GPU pragmas are stale
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        Jnorm(iS,iN_X) = WNORM( J(:,iN_X,iS), W2_S )

      END IF

    END DO
    END DO

  END SUBROUTINE ComputeJNorm


  SUBROUTINE ComputeMatterRHS_FP &
    ( MASK, n_FP, Fm, Gm, J, H_d_1, H_d_2, H_d_3, &
             C_Y    , S_Y    , U_Y    , &
             C_Ef   , S_Ef   , U_Ef   , &
      V_d_1, C_V_d_1, S_V_d_1, U_V_d_1, &
      V_d_2, C_V_d_2, S_V_d_2, U_V_d_2, &
      V_d_3, C_V_d_3, S_V_d_3, U_V_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    REAL(DP), INTENT(inout) :: Fm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: Gm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(in)    :: J    (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: H_d_1(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: H_d_2(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: H_d_3(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: C_Y     (1:nX_G)
    REAL(DP), INTENT(in)    :: S_Y     (1:nX_G)
    REAL(DP), INTENT(in)    :: U_Y     (1:nX_G)
    REAL(DP), INTENT(in)    :: C_Ef    (1:nX_G)
    REAL(DP), INTENT(in)    :: S_Ef    (1:nX_G)
    REAL(DP), INTENT(in)    :: U_Ef    (1:nX_G)
    REAL(DP), INTENT(in)    :: V_d_1   (1:nX_G)
    REAL(DP), INTENT(in)    :: C_V_d_1 (1:nX_G)
    REAL(DP), INTENT(in)    :: S_V_d_1 (1:nX_G)
    REAL(DP), INTENT(in)    :: U_V_d_1 (1:nX_G)
    REAL(DP), INTENT(in)    :: V_d_2   (1:nX_G)
    REAL(DP), INTENT(in)    :: C_V_d_2 (1:nX_G)
    REAL(DP), INTENT(in)    :: S_V_d_2 (1:nX_G)
    REAL(DP), INTENT(in)    :: U_V_d_2 (1:nX_G)
    REAL(DP), INTENT(in)    :: V_d_3   (1:nX_G)
    REAL(DP), INTENT(in)    :: C_V_d_3 (1:nX_G)
    REAL(DP), INTENT(in)    :: S_V_d_3 (1:nX_G)
    REAL(DP), INTENT(in)    :: U_V_d_3 (1:nX_G)
    REAL(DP), INTENT(in)    :: Gm_dd_11(1:nX_G)
    REAL(DP), INTENT(in)    :: Gm_dd_22(1:nX_G)
    REAL(DP), INTENT(in)    :: Gm_dd_33(1:nX_G)

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: V_u_1, V_u_2, V_u_3, k_dd(3,3)
    REAL(DP) :: H_u_1, H_u_2, H_u_3, vDotH
    REAL(DP) :: vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: G_Y    (1:nX_G)
    REAL(DP) :: G_Ef   (1:nX_G)
    REAL(DP) :: G_V_d_1(1:nX_G)
    REAL(DP) :: G_V_d_2(1:nX_G)
    REAL(DP) :: G_V_d_3(1:nX_G)
    REAL(DP) :: N_nu    (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP) :: E_nu    (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP) :: F_nu_d_1(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP) :: F_nu_d_2(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP) :: F_nu_d_3(1:nE_G,1:nX_G,1:nSpecies)

    ! --- Fix me: GPU pragmas are stale
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( V_u_1, V_u_2, V_u_3, vDotH, H_u_1, H_u_2, H_u_3, &
    !$OMP          k_dd, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      V_u_1 = V_d_1(iN_X) / Gm_dd_11(iN_X)
      V_u_2 = V_d_2(iN_X) / Gm_dd_22(iN_X)
      V_u_3 = V_d_3(iN_X) / Gm_dd_33(iN_X)

      vDotH =   V_u_1 * H_d_1(iN_E,iN_X,iS) &
              + V_u_2 * H_d_2(iN_E,iN_X,iS) &
              + V_u_3 * H_d_3(iN_E,iN_X,iS)

      H_u_1 = H_d_1(iN_E,iN_X,iS) / Gm_dd_11(iN_X)
      H_u_2 = H_d_2(iN_E,iN_X,iS) / Gm_dd_22(iN_X)
      H_u_3 = H_d_3(iN_E,iN_X,iS) / Gm_dd_33(iN_X)

      k_dd = EddingtonTensorComponents_dd &
               ( J(iN_E,iN_X,iS), H_u_1, H_u_2, H_u_3, &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

      vDotK_d_1 &
        = ( V_u_1 * k_dd(1,1) + V_u_2 * k_dd(2,1) + V_u_3 * k_dd(3,1) )
      vDotK_d_2 &
        = ( V_u_1 * k_dd(1,2) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(3,2) )
      vDotK_d_3 &
        = ( V_u_1 * k_dd(1,3) + V_u_2 * k_dd(2,3) + V_u_3 * k_dd(3,3) )

      vDotK_d_1 = vDotK_d_1 * J(iN_E,iN_X,iS)
      vDotK_d_2 = vDotK_d_2 * J(iN_E,iN_X,iS)
      vDotK_d_3 = vDotK_d_3 * J(iN_E,iN_X,iS)

      ! --- Eulerian Neutrino Number Density ---

      N_nu(iN_E,iN_X,iS) = J(iN_E,iN_X,iS) + vDotH

      ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

      E_nu(iN_E,iN_X,iS) = J(iN_E,iN_X,iS) + Two * vDotH

      ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

      F_nu_d_1(iN_E,iN_X,iS) &
        = H_d_1(iN_E,iN_X,iS) + V_d_1(iN_X) * J(iN_E,iN_X,iS) + vDotK_d_1

      F_nu_d_2(iN_E,iN_X,iS) &
        = H_d_2(iN_E,iN_X,iS) + V_d_2(iN_X) * J(iN_E,iN_X,iS) + vDotK_d_2

      F_nu_d_3(iN_E,iN_X,iS) &
        = H_d_3(iN_E,iN_X,iS) + V_d_3(iN_X) * J(iN_E,iN_X,iS) + vDotK_d_3

    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)

#elif defined(THORNADO_OACC  )

#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      G_Y    (iN_X) = Zero
      G_Ef   (iN_X) = Zero
      G_V_d_1(iN_X) = Zero
      G_V_d_2(iN_X) = Zero
      G_V_d_3(iN_X) = Zero

    END DO

    DO iS = iNuE, iNuE_Bar

      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, LeptonNumber(iS), N_nu(:,:,iS), &
               nE_G, W2_S, 1, One, G_Y, 1 )

    END DO

    DO iS = 1, nSpecies

      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, E_nu    (:,:,iS), nE_G, W3_S, 1, &
               One, G_Ef   , 1 )
      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, F_nu_d_1(:,:,iS), nE_G, W3_S, 1, &
               One, G_V_d_1, 1 )
      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, F_nu_d_2(:,:,iS), nE_G, W3_S, 1, &
               One, G_V_d_2, 1 )
      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, F_nu_d_3(:,:,iS), nE_G, W3_S, 1, &
               One, G_V_d_3, 1 )

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_Y    (iN_X)  = C_Y    (iN_X) - G_Y    (iN_X) * S_Y    (iN_X)
        G_Ef   (iN_X)  = C_Ef   (iN_X) - G_Ef   (iN_X) * S_Ef   (iN_X)
        G_V_d_1(iN_X)  = C_V_d_1(iN_X) - G_V_d_1(iN_X) * S_V_d_1(iN_X)
        G_V_d_2(iN_X)  = C_V_d_2(iN_X) - G_V_d_2(iN_X) * S_V_d_2(iN_X)
        G_V_d_3(iN_X)  = C_V_d_3(iN_X) - G_V_d_3(iN_X) * S_V_d_3(iN_X)

        Gm(iY ,iN_X) = G_Y    (iN_X)
        Gm(iEf,iN_X) = G_Ef   (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

        Fm(iY ,iN_X) = G_Y    (iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = G_Ef   (iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = G_V_d_1(iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = G_V_d_2(iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = G_V_d_3(iN_X) - U_V_d_3(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_FP


  SUBROUTINE ComputeNeutrinoRHS_FP &
    ( MASK, n_FP, FVECm, GVECm, dt, Omega, V_u_1, V_u_2, V_u_3, &
      C_J  , C_H_d_1  , C_H_d_2  , C_H_d_3  , &
      J_new, H_d_1_new, H_d_2_new, H_d_3_new, &
      J0, Chi, Sig, Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    REAL(DP), INTENT(inout) :: FVECm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: GVECm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(in)    :: dt
    REAL(DP), INTENT(in)    :: Omega(1:nX_G)
    REAL(DP), INTENT(in)    :: V_u_1(1:nX_G)
    REAL(DP), INTENT(in)    :: V_u_2(1:nX_G)
    REAL(DP), INTENT(in)    :: V_u_3(1:nX_G)
    REAL(DP), INTENT(in)    :: C_J      (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: C_H_d_1  (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: C_H_d_2  (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: C_H_d_3  (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: J_new    (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: H_d_1_new(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: H_d_2_new(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: H_d_3_new(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: J0       (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Chi      (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Sig      (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Chi_NES  (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Eta_NES  (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Chi_Pair (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Eta_Pair (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(in)    :: Gm_dd_11(1:nX_G)
    REAL(DP), INTENT(in)    :: Gm_dd_22(1:nX_G)
    REAL(DP), INTENT(in)    :: Gm_dd_33(1:nX_G)

    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: H_u_1, H_u_2, H_u_3, vDotH, k_dd(3,3)
    REAL(DP) :: vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: Eta_T, Chi_T, Kappa

    ! --- Fix me: GPU Pragmas are Stale ---
#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T, Kappa )
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T, Kappa )
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( vDotH, H_u_1, H_u_2, H_u_3, k_dd, vDotK_d_1, &
    !$OMP          vDotK_d_2, vDotK_d_3, Eta_T, Chi_T, Kappa, iOS )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        vDotH =   V_u_1(iN_X) * H_d_1_new(iN_E,iN_X,iS) &
                + V_u_2(iN_X) * H_d_2_new(iN_E,iN_X,iS) &
                + V_u_3(iN_X) * H_d_3_new(iN_E,iN_X,iS)

        H_u_1 = H_d_1_new(iN_E,iN_X,iS) / Gm_dd_11(iN_X)
        H_u_2 = H_d_2_new(iN_E,iN_X,iS) / Gm_dd_22(iN_X)
        H_u_3 = H_d_3_new(iN_E,iN_X,iS) / Gm_dd_33(iN_X)

        k_dd = EddingtonTensorComponents_dd &
                 (  J_new(iN_E,iN_X,iS), H_u_1, H_u_2, H_u_3, &
                    Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vDotK_d_1 = (   V_u_1(iN_X) * k_dd(1,1) &
                      + V_u_2(iN_X) * k_dd(2,1) &
                      + V_u_3(iN_X) * k_dd(3,1) ) * J_new(iN_E,iN_X,iS)

        vDotK_d_2 = (   V_u_1(iN_X) * k_dd(1,2) &
                      + V_u_2(iN_X) * k_dd(2,2) &
                      + V_u_3(iN_X) * k_dd(3,2) ) * J_new(iN_E,iN_X,iS)

        vDotK_d_3 = (   V_u_1(iN_X) * k_dd(1,3) &
                      + V_u_2(iN_X) * k_dd(2,3) &
                      + V_u_3(iN_X) * k_dd(3,3) ) * J_new(iN_E,iN_X,iS)

        ! --- Emissivity ---

        Eta_T = Chi       (iN_E,iN_X,iS) * J0(iN_E,iN_X,iS) &
                + Eta_NES (iN_E,iN_X,iS) &
                + Eta_Pair(iN_E,iN_X,iS)

        ! --- Number Opacity ---

        Chi_T = Chi       (iN_E,iN_X,iS) &
                + Chi_NES (iN_E,iN_X,iS) &
                + Chi_Pair(iN_E,iN_X,iS)

        ! --- Number Flux Opacity ---

        Kappa = Chi_T + Sig(iN_E,iN_X,iS)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! --- Number Equation ---

        GVECm(iOS+iCR_N,iN_X) &
          = ( One - Omega(iN_X) ) * J_new(iN_E,iN_X,iS) &
            + Omega(iN_X) * ( C_J(iN_E,iN_X,iS) - vDotH + dt * Eta_T ) &
                            / ( One + dt * Chi_T )

        FVECm(iOS+iCR_N,iN_X) &
          = GVECm(iOS+iCR_N,iN_X) - J_new(iN_E,iN_X,iS)

        ! --- Number Flux 1 Equation ---

        GVECm(iOS+iCR_G1,iN_X) &
          = ( One - Omega(iN_X) ) * H_d_1_new(iN_E,iN_X,iS) &
            + Omega(iN_X) * ( C_H_d_1(iN_E,iN_X,iS) - vDotK_d_1 ) &
                            / ( One + dt * Kappa )

        FVECm(iOS+iCR_G1,iN_X) &
          = GVECm(iOS+iCR_G1,iN_X) - H_d_1_new(iN_E,iN_X,iS)

        ! --- Number Flux 2 Equation ---

        GVECm(iOS+iCR_G2,iN_X) &
          = ( One - Omega(iN_X) ) * H_d_2_new(iN_E,iN_X,iS) &
            + Omega(iN_X) * ( C_H_d_2(iN_E,iN_X,iS) - vDotK_d_2 ) &
                            / ( One + dt * Kappa )

        FVECm(iOS+iCR_G2,iN_X) &
          = GVECm(iOS+iCR_G2,iN_X) - H_d_2_new(iN_E,iN_X,iS)

        ! --- Number Flux 3 Equation ---

        GVECm(iOS+iCR_G3,iN_X) &
          = ( One - Omega(iN_X) ) * H_d_3_new(iN_E,iN_X,iS) &
            + Omega(iN_X) * ( C_H_d_3(iN_E,iN_X,iS) - vDotK_d_3 ) &
                            / ( One + dt * Kappa )

        FVECm(iOS+iCR_G3,iN_X) &
          = GVECm(iOS+iCR_G3,iN_X) - H_d_3_new(iN_E,iN_X,iS)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_FP


  SUBROUTINE UpdateMatterRHS_FP &
    ( MASK, n_FP, Y, Y_old, U_Y, Ef, Ef_old, U_Ef, &
      V_d_1, V_d_1_old, U_V_d_1, V_d_2, V_d_2_old, U_V_d_2, &
      V_d_3, V_d_3_old, U_V_d_3, FVECm, GVECm )

    LOGICAL,  INTENT(in)    :: MASK     (1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    REAL(DP), INTENT(inout) :: Y        (1:nX_G)
    REAL(DP), INTENT(in)    :: Y_old    (1:nX_G)
    REAL(DP), INTENT(inout) :: U_Y      (1:nX_G)
    REAL(DP), INTENT(inout) :: Ef       (1:nX_G)
    REAL(DP), INTENT(in)    :: Ef_old   (1:nX_G)
    REAL(DP), INTENT(inout) :: U_Ef     (1:nX_G)
    REAL(DP), INTENT(inout) :: V_d_1    (1:nX_G)
    REAL(DP), INTENT(in)    :: V_d_1_old(1:nX_G)
    REAL(DP), INTENT(inout) :: U_V_d_1  (1:nX_G)
    REAL(DP), INTENT(inout) :: V_d_2    (1:nX_G)
    REAL(DP), INTENT(in)    :: V_d_2_old(1:nX_G)
    REAL(DP), INTENT(inout) :: U_V_d_2  (1:nX_G)
    REAL(DP), INTENT(inout) :: V_d_3    (1:nX_G)
    REAL(DP), INTENT(in)    :: V_d_3_old(1:nX_G)
    REAL(DP), INTENT(inout) :: U_V_d_3  (1:nX_G)
    REAL(DP), INTENT(inout) :: FVECm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: GVECm(1:n_FP,1:nX_G)

    INTEGER :: iN_X

#if defined  ( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        FVECm(iY ,iN_X) = GVECm(iY ,iN_X) - U_Y    (iN_X)
        FVECm(iEf,iN_X) = GVECm(iEf,iN_X) - U_Ef   (iN_X)
        FVECm(iV1,iN_X) = GVECm(iV1,iN_X) - U_V_d_1(iN_X)
        FVECm(iV2,iN_X) = GVECm(iV2,iN_X) - U_V_d_2(iN_X)
        FVECm(iV3,iN_X) = GVECm(iV3,iN_X) - U_V_d_3(iN_X)

        U_Y    (iN_X) = GVECm(iY ,iN_X)
        U_Ef   (iN_X) = GVECm(iEf,iN_X)
        U_V_d_1(iN_X) = GVECm(iV1,iN_X)
        U_V_d_2(iN_X) = GVECm(iV2,iN_X)
        U_V_d_3(iN_X) = GVECm(iV3,iN_X)

        Y    (iN_X) = U_Y    (iN_X) * Y_old (iN_X)
        Ef   (iN_X) = U_Ef   (iN_X) * Ef_old(iN_X)
        V_d_1(iN_X) = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2(iN_X) = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3(iN_X) = U_V_d_3(iN_X) * SpeedOfLight

      END IF

    END DO

  END SUBROUTINE UpdateMatterRHS_FP


  SUBROUTINE UpdateNeutrinoRHS_FP &
    ( MASK, n_FP, FVECm, GVECm, J_new, H_d_1_new, H_d_2_new, H_d_3_new )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    REAL(DP), INTENT(inout) :: FVECm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: GVECm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(inout) :: J_new    (1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout) :: H_d_1_new(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout) :: H_d_2_new(1:nE_G,1:nX_G,1:nSpecies)
    REAL(DP), INTENT(inout) :: H_d_3_new(1:nE_G,1:nX_G,1:nSpecies)

    INTEGER  :: iN_E, iN_X, iS, iOS

    ! --- Fix me: GPU Pragmas are Stale ---
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iOS )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        FVECm    (iOS+iCR_N ,iN_X) &
          = GVECm(iOS+iCR_N ,iN_X) - J_new    (iN_E,iN_X,iS)
        FVECm    (iOS+iCR_G1,iN_X) &
          = GVECm(iOS+iCR_G1,iN_X) - H_d_1_new(iN_E,iN_X,iS)
        FVECm    (iOS+iCR_G2,iN_X) &
          = GVECm(iOS+iCR_G2,iN_X) - H_d_2_new(iN_E,iN_X,iS)
        FVECm    (iOS+iCR_G3,iN_X) &
          = GVECm(iOS+iCR_G3,iN_X) - H_d_3_new(iN_E,iN_X,iS)

        J_new    (iN_E,iN_X,iS) = GVECm(iOS+iCR_N ,iN_X)
        H_d_1_new(iN_E,iN_X,iS) = GVECm(iOS+iCR_G1,iN_X)
        H_d_2_new(iN_E,iN_X,iS) = GVECm(iOS+iCR_G2,iN_X)
        H_d_3_new(iN_E,iN_X,iS) = GVECm(iOS+iCR_G3,iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE UpdateNeutrinoRHS_FP


  SUBROUTINE SolveLS_FP &
    ( MASK, n_FP, M, Mk, Fm, Gm, F, G, A, B, Alpha, TAU, LWORK, WORK )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    INTEGER,  INTENT(in)    :: M
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(inout) :: Fm   (1:n_FP    ,1:nX_G)
    REAL(DP), INTENT(inout) :: Gm   (1:n_FP    ,1:nX_G)
    REAL(DP), INTENT(inout) :: B    (1:n_FP    ,1:nX_G)
    REAL(DP), INTENT(inout) :: F    (1:n_FP,1:M,1:nX_G)
    REAL(DP), INTENT(inout) :: G    (1:n_FP,1:M,1:nX_G)
    REAL(DP), INTENT(inout) :: A    (1:n_FP,1:M,1:nX_G)
    REAL(DP), INTENT(inout) :: Alpha(       1:M,1:nX_G)
    REAL(DP), INTENT(inout) :: TAU  (1:n_FP    ,1:nX_G)
    INTEGER,  INTENT(inout) :: LWORK
    REAL(DP), INTENT(inout) :: WORK(1:LWORK    ,1:nX_G)

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
      !$OMP DO COLLAPSE(2)
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
      !$OMP DO COLLAPSE(3)
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
        !$OMP PARALLEL DO SIMD &
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

            B(1,iN_X) = AB1 / AA11

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

            B(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
            B(2,iN_X) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

          END IF
        END DO

      ELSE IF ( Mk > 3 ) THEN

        DO iN_X = 1, nX_G
          IF( MASK(iN_X) )THEN

            CALL LinearLeastSquares &
                   ( 'N', n_FP, Mk-1, 1, A(:,:,iN_X), n_FP, B(:,iN_X), n_FP, &
                     TAU(1,iN_X), WORK(1,iN_X), LWORK, INFO(iN_X) )

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
      !$OMP PARALLEL DO SIMD &
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
      !$OMP DO SIMD COLLAPSE(2) &
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

    END IF

  END SUBROUTINE SolveLS_FP


  SUBROUTINE ShiftRHS_FP( MASK, n_FP, M, Mk, F, G )

    LOGICAL,  INTENT(in)    :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    INTEGER,  INTENT(in)    :: M
    INTEGER,  INTENT(in)    :: Mk
    REAL(DP), INTENT(inout) :: F(1:n_FP,1:M,1:nX_G)
    REAL(DP), INTENT(inout) :: G(1:n_FP,1:M,1:nX_G)

    INTEGER  :: iN_X, iFP, iM
    REAL(DP) :: FTMP(1:n_FP,1:M)
    REAL(DP) :: GTMP(1:n_FP,1:M)

    IF ( Mk == M ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( FTMP, GTMP )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( FTMP, GTMP )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( FTMP, GTMP )
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
            FTMP(iFP,iM) = F(iFP,iM+1,iN_X)
            GTMP(iFP,iM) = G(iFP,iM+1,iN_X)
          END DO
          END DO

#if   defined( THORNADO_OMP_OL )
          !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM  = 1, Mk-1
          DO iFP = 1, n_FP
            F(iFP,iM,iN_X) = FTMP(iFP,iM)
            G(iFP,iM,iN_X) = GTMP(iFP,iM)
          END DO
          END DO

        END IF
      END DO

    END IF

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE CheckConvergence_Inner &
    ( MASK, n_FP, k_inner, nIterations_Inner, Fm, Jnorm )

    LOGICAL,  INTENT(inout) :: MASK(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    INTEGER,  INTENT(in)    :: k_inner
    INTEGER,  INTENT(inout) :: nIterations_Inner(1:nX_G)
    REAL(DP), INTENT(in)    :: Fm(1:n_FP,1:nX_G)
    REAL(DP), INTENT(in)    :: Jnorm(1:nSpecies,1:nX_G)

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X, iS, iCR, iN_E, iOS
    REAL(DP) :: Fnorm(1:nSpecies,1:nX_G)
    REAL(DP) :: F(1:nE_G,1:nCR,1:nSpecies,1:nX_G)

    ! --- Permute Fm ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( iOS )
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iCR  = 1, nCR
    DO iN_E = 1, nE_G

      iOS = ( (iS-1) * nE_G + (iN_E-1) ) * nCR

      F(iN_E,iCR,iS,iN_X) = Fm(iOS+iCR,iN_X)

    END DO
    END DO
    END DO
    END DO

    ! --- Compute Norm of F (Use Maximum Across Moments) ---

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies

      IF( MASK(iN_X) )THEN
      
        Fnorm(iS,iN_X) = Zero      
        DO iCR  = 1, nCR

          Fnorm(iS,iN_X) &
            = MAX( Fnorm(iS,iN_X), WNORM( F(:,iCR,iS,iN_X), W2_S ) )

        END DO

      END IF

    END DO
    END DO

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )

#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED )
#endif
    DO iN_X = 1, nX_G

      IF( MASK(iN_X) )THEN

        CONVERGED = ALL( Fnorm(:,iN_X) <= Rtol_inner * Jnorm(:,iN_X) )

        IF( CONVERGED )THEN

          MASK(iN_X) = .FALSE.

          nIterations_Inner(iN_X) &
            = nIterations_Inner(iN_X) + k_inner

        END IF

      END IF

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( MASK )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( MASK )
#endif

  END SUBROUTINE CheckConvergence_Inner


  SUBROUTINE CheckConvergence_Outer &
    ( MASK_OUTER, MASK_INNER, n_FP, k_outer, Fm, nIterations_Outer )

    LOGICAL,  INTENT(inout) :: MASK_OUTER(1:nX_G)
    LOGICAL,  INTENT(inout) :: MASK_INNER(1:nX_G)
    INTEGER,  INTENT(in)    :: n_FP
    INTEGER,  INTENT(in)    :: k_outer
    REAL(DP), INTENT(in)    :: Fm(1:n_FP,1:nX_G)
    INTEGER,  INTENT(inout) :: nIterations_Outer(1:nX_G)

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X
    REAL(DP) :: Fnorm_Y, Fnorm_Ef, Fnorm_V

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( Fnorm_Y, Fnorm_Ef, Fnorm_V, CONVERGED )
#endif
    DO iN_X = 1, nX_G
      IF( MASK_OUTER(iN_X) )THEN

        Fnorm_Y  =      ABS( Fm(iY ,iN_X) )
        Fnorm_Ef =      ABS( Fm(iEf,iN_X) )
        Fnorm_V  = MAX( ABS( Fm(iV1,iN_X) ), &
                        ABS( Fm(iV2,iN_X) ), &
                        ABS( Fm(iV3,iN_X) ) )

        CONVERGED = Fnorm_Y  <= Rtol_outer .AND. &
                    Fnorm_Ef <= Rtol_outer .AND. &
                    Fnorm_V  <= Rtol_outer

        IF( CONVERGED )THEN
          MASK_OUTER(iN_X) = .FALSE.
          nIterations_Outer(iN_X) = k_outer
        END IF

        MASK_INNER(iN_X) = MASK_OUTER(iN_X)

      END IF
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE FROM( MASK_OUTER, MASK_INNER )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE HOST( MASK_OUTER, MASK_INNER )
#endif

  END SUBROUTINE CheckConvergence_Outer


  FUNCTION ENORM( X )
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


END MODULE TwoMoment_NeutrinoMatterSolverModule_OrderV
