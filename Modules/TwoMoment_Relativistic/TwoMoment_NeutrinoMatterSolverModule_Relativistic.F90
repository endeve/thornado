MODULE TwoMoment_NeutrinoMatterSolverModule_Relativistic



  USE KindModule, ONLY: &
    DP, Zero, One, Two, FourPi, Half
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    AtomicMassUnit, &
    Centimeter, &
    Gram, &
    MeV, &
    Erg, &
    Kelvin, &
    Dyne, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, &
    nDOFE
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
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputePressure_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_DG_Points, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacities_Pair_Points, &
    ComputeNeutrinoOpacitiesRates_NES_Points, &
    ComputeNeutrinoOpacitiesRates_Pair_Points
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputeEddingtonTensorComponents_dd

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: InitializeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling_FP_Nested_AA

  LOGICAL, PARAMETER :: Include_NES  = .TRUE.
  LOGICAL, PARAMETER :: Include_Pair = .TRUE.

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

  ! --- Solver scratch arrays ---

  LOGICAL,  DIMENSION(:), ALLOCATABLE :: ITERATE_outer, ITERATE_inner
  INTEGER,  DIMENSION(:), ALLOCATABLE :: PackIndex_outer, UnpackIndex_outer
  INTEGER,  DIMENSION(:), ALLOCATABLE :: PackIndex_inner, UnpackIndex_inner

  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: Jnorm
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_J
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_H_d_1, H_d_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_H_d_2, H_d_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_H_d_3, H_d_3

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3

  REAL(DP), DIMENSION(:), ALLOCATABLE :: Omega
  REAL(DP), DIMENSION(:), ALLOCATABLE :: S_Eps, C_Eps, U_Eps, G_Eps, Eps_old
  REAL(DP), DIMENSION(:), ALLOCATABLE :: cD_old, Y_old, C_Y, S_Y, G_Y, U_Y
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rho_old, C_rho, S_rho, G_rho, U_rho
  REAL(DP), DIMENSION(:), ALLOCATABLE :: C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3

  REAL(DP), DIMENSION(:,:,:)  , ALLOCATABLE, TARGET :: J0, Chi, Sig
  REAL(DP), DIMENSION(:,:,:)  , ALLOCATABLE, TARGET :: Chi_NES, Eta_NES
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: Phi_0_In_NES, Phi_0_Ot_NES
  REAL(DP), DIMENSION(:,:,:)  , ALLOCATABLE, TARGET :: Chi_Pair, Eta_Pair
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: Phi_0_In_Pair, Phi_0_Ot_Pair

  ! --- Least-squares scratch arrays ---

  INTEGER                                 :: LWORK_outer
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: AMAT_outer, GVEC_outer, FVEC_outer
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: BVEC_outer, GVECm_outer, FVECm_outer
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: WORK_outer, TAU_outer, Alpha_outer

  INTEGER                                 :: LWORK_inner
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: AMAT_inner, GVEC_inner, FVEC_inner
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: BVEC_inner, GVECm_inner, FVECm_inner
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: WORK_inner, TAU_inner, Alpha_inner

  INTEGER,  DIMENSION(:)    , ALLOCATABLE :: INFO

  ! --- Solver Parameters to be initialized

  INTEGER  :: M_FP, M_outer, M_inner
  INTEGER  :: MaxIter_outer, MaxIter_inner
  REAL(DP) :: Rtol_outer, Rtol_inner

  INTEGER :: iS_1 = iNuE
  INTEGER :: iS_2 = iNuE_Bar

  INTEGER, PARAMETER :: iY  = 1
  INTEGER, PARAMETER :: iEps= 2
  INTEGER, PARAMETER :: iV1 = 3
  INTEGER, PARAMETER :: iV2 = 4
  INTEGER, PARAMETER :: iV3 = 5
  INTEGER, PARAMETER :: iD  = 6

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

    n_FP_outer = 6 
    n_FP_inner = nE_G * nCR * nSpecies

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S(:) = WFactor_FP * W2_N(:)
    W3_S(:) = WFactor_FP * W3_N(:)

    ALLOCATE(     ITERATE_outer(nX_G) )
    ALLOCATE(   PackIndex_outer(nX_G) )
    ALLOCATE( UnpackIndex_outer(nX_G) )

    ALLOCATE(     ITERATE_inner(nX_G) )
    ALLOCATE(   PackIndex_inner(nX_G) )
    ALLOCATE( UnpackIndex_inner(nX_G) )

    ITERATE_outer(:) = .TRUE.
    ITERATE_inner(:) = .TRUE.

    ALLOCATE( Jnorm(nSpecies,nX_G) )

    ALLOCATE( C_J  (nE_G,nX_G,nSpecies) )

    ALLOCATE( C_H_d_1    (nE_G,nX_G,nSpecies) )
    ALLOCATE(   H_d_1    (nE_G,nX_G,nSpecies) )

    ALLOCATE( C_H_d_2    (nE_G,nX_G,nSpecies) )
    ALLOCATE(   H_d_2    (nE_G,nX_G,nSpecies) )

    ALLOCATE( C_H_d_3    (nE_G,nX_G,nSpecies) )
    ALLOCATE(   H_d_3    (nE_G,nX_G,nSpecies) )

    ALLOCATE( N_nu    (nE_G,nX_G,nSpecies) )
    ALLOCATE( E_nu    (nE_G,nX_G,nSpecies) )
    ALLOCATE( F_nu_d_1(nE_G,nX_G,nSpecies) )
    ALLOCATE( F_nu_d_2(nE_G,nX_G,nSpecies) )
    ALLOCATE( F_nu_d_3(nE_G,nX_G,nSpecies) )

    ALLOCATE( Omega  (nX_G) )

    ALLOCATE( Eps_old   (nX_G) )
    ALLOCATE( S_Eps     (nX_G) )
    ALLOCATE( C_Eps     (nX_G) )
    ALLOCATE( G_Eps     (nX_G) )
    ALLOCATE( U_Eps     (nX_G) )

    ALLOCATE(   Y_old(nX_G) )
    ALLOCATE( C_Y    (nX_G) )
    ALLOCATE( S_Y    (nX_G) )
    ALLOCATE( G_Y    (nX_G) )
    ALLOCATE( U_Y    (nX_G) )

    ALLOCATE(   rho_old(nX_G) )
    ALLOCATE(   cD_old (nX_G) )
    ALLOCATE( C_rho    (nX_G) )
    ALLOCATE( S_rho    (nX_G) )
    ALLOCATE( G_rho    (nX_G) )
    ALLOCATE( U_rho    (nX_G) )

    ALLOCATE( C_V_d_1    (nX_G) )
    ALLOCATE( S_V_d_1    (nX_G) )
    ALLOCATE( G_V_d_1    (nX_G) )
    ALLOCATE( U_V_d_1    (nX_G) )
    ALLOCATE(   V_d_1    (nX_G) )

    ALLOCATE( C_V_d_2    (nX_G) )
    ALLOCATE( S_V_d_2    (nX_G) )
    ALLOCATE( G_V_d_2    (nX_G) )
    ALLOCATE( U_V_d_2    (nX_G) )
    ALLOCATE(   V_d_2    (nX_G) )

    ALLOCATE( C_V_d_3    (nX_G) )
    ALLOCATE( S_V_d_3    (nX_G) )
    ALLOCATE( G_V_d_3    (nX_G) )
    ALLOCATE( U_V_d_3    (nX_G) )
    ALLOCATE(   V_d_3    (nX_G) )

    ALLOCATE(  J0(nE_G,nX_G,nSpecies) )
    ALLOCATE( Chi(nE_G,nX_G,nSpecies) )
    ALLOCATE( Sig(nE_G,nX_G,nSpecies) )

    ALLOCATE(      Chi_NES(     nE_G,nX_G,nSpecies) )
    ALLOCATE(      Eta_NES(     nE_G,nX_G,nSpecies) )
    ALLOCATE( Phi_0_In_NES(nE_G,nE_G,nX_G,nSpecies) )
    ALLOCATE( Phi_0_Ot_NES(nE_G,nE_G,nX_G,nSpecies) )

    ALLOCATE(      Chi_Pair(     nE_G,nX_G,nSpecies) )
    ALLOCATE(      Eta_Pair(     nE_G,nX_G,nSpecies) )
    ALLOCATE( Phi_0_In_Pair(nE_G,nE_G,nX_G,nSpecies) )
    ALLOCATE( Phi_0_Ot_Pair(nE_G,nE_G,nX_G,nSpecies) )

    ALLOCATE(  AMAT_outer(n_FP_outer,M_outer,nX_G) )
    ALLOCATE(  GVEC_outer(n_FP_outer,M_outer,nX_G) )
    ALLOCATE(  FVEC_outer(n_FP_outer,M_outer,nX_G) )
    ALLOCATE(  BVEC_outer(n_FP_outer,        nX_G) )
    ALLOCATE( GVECm_outer(n_FP_outer,        nX_G) )
    ALLOCATE( FVECm_outer(n_FP_outer,        nX_G) )
    ALLOCATE(   TAU_outer(n_FP_outer,        nX_G) )
    ALLOCATE( Alpha_outer(           M_outer,nX_G) )

    ALLOCATE(  AMAT_inner(n_FP_inner,M_inner,nX_G) )
    ALLOCATE(  GVEC_inner(n_FP_inner,M_inner,nX_G) )
    ALLOCATE(  FVEC_inner(n_FP_inner,M_inner,nX_G) )
    ALLOCATE(  BVEC_inner(n_FP_inner,        nX_G) )
    ALLOCATE( GVECm_inner(n_FP_inner,        nX_G) )
    ALLOCATE( FVECm_inner(n_FP_inner,        nX_G) )
    ALLOCATE(   TAU_inner(n_FP_inner,        nX_G) )
    ALLOCATE( Alpha_inner(           M_inner,nX_G) )

    ALLOCATE( INFO(nX_G) )

    ALLOCATE( P1D(          nX_G,nP1D) )
    ALLOCATE( P2D(nE_G     ,nX_G,nP2D) )
    ALLOCATE( P3D(nE_G,nE_G,nX_G,nP3D) )

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S, &
    !$OMP          ITERATE_outer, ITERATE_inner ) &
    !$OMP MAP( alloc: INFO, P1D, P2D, P3D, &
    !$OMP             Omega, Jnorm, &
    !$OMP             C_J, &
    !$OMP             C_H_d_1, H_d_1, &
    !$OMP             C_H_d_2, H_d_2, &
    !$OMP             C_H_d_3, H_d_3, &
    !$OMP             N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3, &
    !$OMP             Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$OMP             Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$OMP             C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$OMP             C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$OMP             C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$OMP             J0, Chi, Sig, &
    !$OMP             Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$OMP             Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$OMP             PackIndex_outer, UnpackIndex_outer, &
    !$OMP             AMAT_outer, GVEC_outer, FVEC_outer, &
    !$OMP             BVEC_outer, GVECm_outer, FVECm_outer, &
    !$OMP             TAU_outer, Alpha_outer, &
    !$OMP             PackIndex_inner, UnpackIndex_inner, &
    !$OMP             AMAT_inner, GVEC_inner, FVEC_inner, &
    !$OMP             BVEC_inner, GVECm_inner, FVECm_inner, &
    !$OMP             TAU_inner, Alpha_inner )
#elif defined(THORNADO_OACC  )
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S, &
    !$ACC         ITERATE_outer, ITERATE_inner ) &
    !$ACC CREATE( INFO, P1D, P2D, P3D, &
    !$ACC         Omega, Jnorm, &
    !$ACC         C_J, &
    !$ACC         C_H_d_1, H_d_1, &
    !$ACC         C_H_d_2, H_d_2, &
    !$ACC         C_H_d_3, H_d_3, &
    !$ACC         N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3, &
    !$ACC         Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$ACC         Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$ACC         C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$ACC         C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$ACC         C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$ACC         J0, Chi, Sig, &
    !$ACC         Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$ACC         Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, &
    !$ACC         AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC         BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC         TAU_outer, Alpha_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC         BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC         TAU_inner, Alpha_inner )
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
      M_outer = 2
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
    !$OMP               Omega, Jnorm, &
    !$OMP               C_J, &
    !$OMP               C_H_d_1, H_d_1, &
    !$OMP               C_H_d_2, H_d_2, &
    !$OMP               C_H_d_3, H_d_3, &
    !$OMP               N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3, &
    !$OMP               Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$OMP               Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$OMP               C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$OMP               C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$OMP               C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$OMP               J0, Chi, Sig, &
    !$OMP               Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$OMP               Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$OMP               ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
    !$OMP               AMAT_outer, GVEC_outer, FVEC_outer, &
    !$OMP               BVEC_outer, GVECm_outer, FVECm_outer, &
    !$OMP               WORK_outer, TAU_outer, Alpha_outer, &
    !$OMP               ITERATE_inner, PackIndex_inner, UnpackIndex_inner, &
    !$OMP               AMAT_inner, GVEC_inner, FVEC_inner, &
    !$OMP               BVEC_inner, GVECm_inner, FVECm_inner, &
    !$OMP               WORK_inner, TAU_inner, Alpha_inner )
#elif defined(THORNADO_OACC  )
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, &
    !$ACC         INFO, P1D, P2D, P3D, &
    !$ACC         Omega, Jnorm, &
    !$ACC         C_J, &
    !$ACC         C_H_d_1, H_d_1, &
    !$ACC         C_H_d_2, H_d_2, &
    !$ACC         C_H_d_3, H_d_3, &
    !$ACC         N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3, &
    !$ACC         Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$ACC         Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$ACC         C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$ACC         C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$ACC         C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$ACC         J0, Chi, Sig, &
    !$ACC         Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES, &
    !$ACC         Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair, &
    !$ACC         ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
    !$ACC         AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC         BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC         WORK_outer, TAU_outer, Alpha_outer, &
    !$ACC         ITERATE_inner, PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC         BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC         WORK_inner, TAU_inner, Alpha_inner )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S )
    DEALLOCATE( INFO, P1D, P2D, P3D )
    DEALLOCATE( Omega, Jnorm )
    DEALLOCATE( C_J )
    DEALLOCATE( C_H_d_1, H_d_1 )
    DEALLOCATE( C_H_d_2, H_d_2 )
    DEALLOCATE( C_H_d_3, H_d_3 )
    DEALLOCATE( N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 )
    DEALLOCATE( Eps_old, C_Eps, S_Eps, U_Eps, G_Eps )
    DEALLOCATE( Y_old, C_Y, S_Y, G_Y, U_Y )
    DEALLOCATE( cD_old, rho_old, C_rho, S_rho, G_rho, U_rho )
    DEALLOCATE( C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1 )
    DEALLOCATE( C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2 )
    DEALLOCATE( C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3 )
    DEALLOCATE( J0, Chi, Sig )
    DEALLOCATE( Chi_NES, Eta_NES, Phi_0_In_NES, Phi_0_Ot_NES )
    DEALLOCATE( Chi_Pair, Eta_Pair, Phi_0_In_Pair, Phi_0_Ot_Pair )
    DEALLOCATE( ITERATE_outer, PackIndex_outer, UnpackIndex_outer )
    DEALLOCATE( AMAT_outer, GVEC_outer, FVEC_outer )
    DEALLOCATE( BVEC_outer, GVECm_outer, FVECm_outer )
    DEALLOCATE( WORK_outer, TAU_outer, Alpha_outer )
    DEALLOCATE( ITERATE_inner, PackIndex_inner, UnpackIndex_inner )
    DEALLOCATE( AMAT_inner, GVEC_inner, FVEC_inner )
    DEALLOCATE( BVEC_inner, GVECm_inner, FVECm_inner )
    DEALLOCATE( WORK_inner, TAU_inner, Alpha_inner )

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
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3, &
      nIterations_Inner, nIterations_Outer )

    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: D, T, Y, E
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:),     INTENT(inout) :: alp, B_u_1, B_u_2, B_u_3
    INTEGER,  DIMENSION(:),     INTENT(out)   :: nIterations_Inner, nIterations_Outer

    ! --- Local Variables ---

    INTEGER  :: k_outer, Mk_outer, nX_P_outer
    INTEGER  :: k_inner, Mk_inner, nX_P_inner
    REAL(DP) :: P(SIZE(D))

    ! --- Initial RHS ---
    CALL ComputePressure_TABLE & 
           ( D, T, Y, P )
    
    CALL InitializeRHS_FP &
           ( J, H_u_1, H_u_2, H_u_3, D, Y, E, P, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )

    ! --- Compute Opacity Kernels ---


    CALL ComputeOpacities_Packed &
           ( D, T, Y )


    ! --- Start Outer Loop ---

    k_outer = 0
    DO WHILE( ANY( ITERATE_outer(:) ) .AND. k_outer < MaxIter_outer )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      CALL ComputeJNorm &
             ( ITERATE_outer, J )

      CALL CreatePackIndex &
             ( ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      IF ( k_outer > 1 ) THEN

        ! --- Recompute Opacity Kernels ---


        CALL ComputeOpacities_Packed &
               ( D, T, Y, &
                 ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )


      END IF

      ! --- Start Inner Loop ---


      k_inner = 0
      DO WHILE( ANY( ITERATE_inner(:) ) .AND. k_inner < MaxIter_inner )

        k_inner  = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_inner, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---


        CALL ComputeRates_Packed &
               ( J, ITERATE_inner, nX_P_inner, PackIndex_inner, &
                 UnpackIndex_inner, nX_P_outer )

        ! --- Right-Hand Side Vectors and Residuals (inner) ---


        CALL ComputeNeutrinoRHS_FP &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, dt, &
                 J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3  )


        ! --- Anderson Acceleration (inner) ---


        CALL SolveLS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner, &
                 AMAT_inner, BVEC_inner, Alpha_inner, TAU_inner, &
                 LWORK_inner, WORK_inner )


        ! --- Update Residuals and Solution Vectors (inner) ---


        CALL UpdateNeutrinoRHS_FP &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                 J, H_u_1, H_u_2, H_u_3, &
                 V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 alp, B_u_1, B_u_2, B_u_3 )


        ! --- Check Convergence (inner) ---


        CALL CheckConvergence_Inner &
               ( ITERATE_inner, n_FP_inner, k_inner, &
                 nIterations_Inner, FVECm_inner )
        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )


      END DO ! --- Inner Loop ---


      ! --- Right-Hand Side Vectors and Residuals (outer) ---


      CALL ComputeMatterRHS_FP &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               J, H_u_1, H_u_2, H_u_3, D, E, P, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )


      ! --- Anderson Acceleration (outer) ---


      CALL SolveLS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVECm_outer, GVECm_outer, FVEC_outer, GVEC_outer, &
               AMAT_outer, BVEC_outer, Alpha_outer, TAU_outer, &
               LWORK_outer, WORK_outer )


      ! --- Update Residuals and Solution Vectors (outer) ---


      CALL UpdateMatterRHS_FP &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               D, Y, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

 

      CALL ComputePressure_TABLE & 
           ( D, T, Y, P )

      ! --- Check Convergence (outer) ---


      CALL CheckConvergence_Outer &
             ( ITERATE_outer, ITERATE_inner, n_FP_outer, k_outer, &
               nIterations_Outer, FVECm_outer )


      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVEC_outer, GVEC_outer )


    END DO ! --- Outer Loop ---

    nIterations_Inner &
      = FLOOR( DBLE( nIterations_Inner ) / DBLE( nIterations_Outer ) )

  END SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA


  SUBROUTINE ComputeOpacities_Packed &
    ( D, T, Y, MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(:), INTENT(in), TARGET   :: D, T, Y
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: MASK
    INTEGER,                INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                INTENT(in), OPTIONAL :: nX_P0

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
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
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
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
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
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
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
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
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
    ( J, MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET   :: J
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: MASK
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:),     INTENT(in), OPTIONAL :: PackIndex, UnpackIndex
    INTEGER,                    INTENT(in), OPTIONAL :: nX_P0

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

    REAL(DP), DIMENSION(:), INTENT(in)   , TARGET   :: D, E, Y
    REAL(DP), DIMENSION(:), INTENT(inout), TARGET   :: T
    LOGICAL,  DIMENSION(:), INTENT(in)   , OPTIONAL :: MASK
    INTEGER,                INTENT(in)   , OPTIONAL :: nX_P
    INTEGER,  DIMENSION(:), INTENT(in)   , OPTIONAL :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:), POINTER :: D_P, E_P, Y_P, T_P

    INTEGER :: nX

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
    ( J, H_u_1, H_u_2, H_u_3, D, Y, E, P, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Y, E, P, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: alp, B_u_1, B_u_2, B_u_3

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotV, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, W, vvDotK
    REAL(DP) :: DT, B_d_1, B_d_2, B_d_3, cD, SJ(3), Ef
#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( vDotV )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( vDotV )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( vDotV )
#endif
    DO iN_X = 1, nX_G

      V_d_1(iN_X) = Gm_dd_11(iN_X) * V_u_1(iN_X)
      V_d_2(iN_X) = Gm_dd_22(iN_X) * V_u_2(iN_X)
      V_d_3(iN_X) = Gm_dd_33(iN_X) * V_u_3(iN_X)

      ! --- Specific Fluid Energy ---

      vDotV =   V_u_1(iN_X) * V_d_1(iN_X) &
              + V_u_2(iN_X) * V_d_2(iN_X) &
              + V_u_3(iN_X) * V_d_3(iN_X)


      W = 1.0_DP / SQRT( 1.0_DP - vDotV)

      
      ! --- Compute Omega (Richardson damping coeff.) based on velocity ---
      ! --- For first interation: must be consistent with initial guess ---

      Omega(iN_X) = One / ( W + SQRT( vDotV ) )
      
      ! --- Store Initial Matter State ---

      rho_old(iN_X)  = D(iN_X)
      Y_old (iN_X)   = Y(iN_X)
      Eps_old (iN_X) = E(iN_X)
      cD_old(iN_X)   = W * D(iN_X)

      ! --- Scaling Factors ---

      S_rho  (iN_X) = One / ( D (iN_X) )
      S_Y    (iN_X) = One / ( Y (iN_X) )
      S_Eps  (iN_X) = One / ( E(iN_X) )
      S_V_d_1(iN_X) = One / ( SpeedOfLight )
      S_V_d_2(iN_X) = One / ( SpeedOfLight )
      S_V_d_3(iN_X) = One / ( SpeedOfLight )
      ! --- Initial Guess for Matter State ---

      U_rho  (iN_X) = One
      U_Y    (iN_X) = One
      U_Eps  (iN_X) = One
      U_V_d_1(iN_X) = V_d_1(iN_X) / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2(iN_X) / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3(iN_X) / SpeedOfLight

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_rho  (iN_X) = Zero
      C_Y    (iN_X) = Zero
      C_Eps  (iN_X) = Zero
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

      vDotV =   V_u_1(iN_X) * V_d_1(iN_X) &
              + V_u_2(iN_X) * V_d_2(iN_X) &
              + V_u_3(iN_X) * V_d_3(iN_X)

      W = 1.0_DP / SQRT(1.0_DP - vDotV)


      B_d_1 = Gm_dd_11(iN_X) * B_u_1(iN_X)
      B_d_2 = Gm_dd_22(iN_X) * B_u_2(iN_X)
      B_d_3 = Gm_dd_33(iN_X) * B_u_3(iN_X)

      DT = 1.0_DP / ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - alp(iN_X) )

      H_d_1(iN_E,iN_X,iS) = DT * ( B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - alp(iN_X) ) & 
                          * Gm_dd_11(iN_X) * H_u_1(iN_E,iN_X,iS) &
                          - DT * ( B_d_1 * V_u_2(iN_X) *Gm_dd_22(iN_X) ) * H_u_2(iN_E,iN_X,iS)   &
                          - DT * ( B_d_1 * V_u_3(iN_X) * Gm_dd_33(iN_X) ) * H_u_3(iN_E,iN_X,iS) 

      H_d_2(iN_E,iN_X,iS) = DT * ( B_d_1 * V_u_1(iN_X) + B_d_3 * V_u_3(iN_X) - alp(iN_X) ) &
                          * Gm_dd_22(iN_X) * H_u_2(iN_E,iN_X,iS) &
                          - DT * ( B_d_2 * V_u_1(iN_X) * Gm_dd_11(iN_X) ) * H_u_1(iN_E,iN_X,iS) &
                          - DT * ( Gm_dd_33(iN_X) * H_u_3(iN_E,iN_X,iS) * B_d_2 * V_u_3(iN_X) ) 

      H_d_3(iN_E,iN_X,iS) = DT * ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) - alp(iN_X) ) &
                          * Gm_dd_33(iN_X) * H_u_3(iN_E,iN_X,iS) &
                          - DT * ( Gm_dd_11(iN_X) * H_u_1(iN_E,iN_X,iS) * B_d_3 * V_u_1(iN_X) ) &
                          - DT * ( Gm_dd_22(iN_X) * H_u_2(iN_E,iN_X,iS) * B_d_3 * V_u_2(iN_X) )

      vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iN_X,iS) &
              + V_u_2(iN_X) * H_d_2(iN_E,iN_X,iS) &
              + V_u_3(iN_X) * H_d_3(iN_E,iN_X,iS)

       CALL ComputeEddingtonTensorComponents_dd &
               ( J(iN_E,iN_X,iS), H_u_1(iN_E,iN_X,iS), &
                 H_u_2(iN_E,iN_X,iS), H_u_3(iN_E,iN_X,iS), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 alp(iN_X), B_u_1(iN_X), B_u_2(iN_X), B_u_3(iN_X), & 
                 V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), k_dd  )

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


       vvDotK = V_u_1(iN_X) * vDotK_d_1 &
              + V_u_2(iN_X) * vDotK_d_2 & 
              + V_u_3(iN_X) * vDotK_d_3 

      ! --- Eulerian Neutrino Number Density ---

      N_nu(iN_E,iN_X,iS) = W * J(iN_E,iN_X,iS) + vDotH

      ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

      E_nu(iN_E,iN_X,iS) = W**2 * J(iN_E,iN_X,iS) +  W * Two * vDotH + vvDotK 

      ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

      F_nu_d_1(iN_E,iN_X,iS) &
        = W * H_d_1(iN_E,iN_X,iS) + W**2 * V_d_1(iN_X) * J(iN_E,iN_X,iS) &
        + W * vDotH * V_d_1(iN_X)  + vDotK_d_1

      F_nu_d_2(iN_E,iN_X,iS) &
        = W * H_d_2(iN_E,iN_X,iS) + W**2 * V_d_2(iN_X) * J(iN_E,iN_X,iS) &
        + W * vDotH * V_d_2(iN_X)  + vDotK_d_2

      F_nu_d_3(iN_E,iN_X,iS) &
        = W * H_d_3(iN_E,iN_X,iS) + W**2 * V_d_3(iN_X) * J(iN_E,iN_X,iS) &
        + W * vDotH * V_d_3(iN_X)  + vDotK_d_3

      ! --- Old States for Neutrino Number Density and Flux ---

      C_J    (iN_E,iN_X,iS) = W * J    (iN_E,iN_X,iS) + vDotH
      C_H_d_1(iN_E,iN_X,iS) = W * H_d_1(iN_E,iN_X,iS) + vDotK_d_1
      C_H_d_2(iN_E,iN_X,iS) = W * H_d_2(iN_E,iN_X,iS) + vDotK_d_2
      C_H_d_3(iN_E,iN_X,iS) = W * H_d_3(iN_E,iN_X,iS) + vDotK_d_3

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
               One, C_Eps  , 1 )
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


      vDotV =   V_u_1(iN_X) * V_d_1(iN_X) &
              + V_u_2(iN_X) * V_d_2(iN_X) &
              + V_u_3(iN_X) * V_d_3(iN_X)

      W = 1.0_DP / SQRT(1.0_DP - vDotV)
 
      C_Y(iN_X) = ( AtomicMassUnit / cD_old(iN_X) ) * C_Y(iN_X)

      SJ(1) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_1(iN_X)
      SJ(2) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_2(iN_X)
      SJ(3) = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 * V_d_3(iN_X)




      Ef = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) * D(iN_X) * W**2 - P(iN_X)
    

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_rho      (iN_X) &
        = W * U_rho  (iN_X) 
      C_Y    (iN_X) &
        = U_Y    (iN_X) + C_Y    (iN_X) * S_Y    (iN_X)
      C_Eps   (iN_X) &
        = ( Ef + C_Eps  (iN_X) ) * S_Eps  (iN_X) / cD_old(iN_X)
      C_V_d_1(iN_X) &
        = ( SJ(1) + C_V_d_1(iN_X) ) * S_V_d_1(iN_X) / cD_old(iN_X)
      C_V_d_2(iN_X) &
        = ( SJ(2) + C_V_d_2(iN_X) ) * S_V_d_2(iN_X) / cD_old(iN_X)
      C_V_d_3(iN_X) &
        = ( SJ(3) + C_V_d_3(iN_X) ) * S_V_d_3(iN_X) / cD_old(iN_X)
    END DO

  END SUBROUTINE InitializeRHS_FP


  SUBROUTINE ComputeJNorm( MASK, J )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: J

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

        Jnorm(iS,iN_X) = WNORM( J(:,iN_X,iS), W2_S )

      END IF

    END DO
    END DO

  END SUBROUTINE ComputeJNorm


  SUBROUTINE ComputeMatterRHS_FP &
    ( MASK, Fm, Gm, &
      J, H_u_1, H_u_2, H_u_3, D, E, P, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: D, E, P, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: alp, B_u_1, B_u_2, B_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotV, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, vvDotK, W, h
    REAL(DP) :: DT, B_d_1, B_d_2, B_d_3

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

      IF ( MASK(iN_X) ) THEN

        vDotV =   V_u_1(iN_X) * V_d_1(iN_X) &
                + V_u_2(iN_X) * V_d_2(iN_X) &
                + V_u_3(iN_X) * V_d_3(iN_X)

        W = 1.0_DP / SQRT(1.0_DP - vDotV)


        B_d_1 = Gm_dd_11(iN_X) * B_u_1(iN_X)
        B_d_2 = Gm_dd_22(iN_X) * B_u_2(iN_X)
        B_d_3 = Gm_dd_33(iN_X) * B_u_3(iN_X)

        DT = 1.0_DP / ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - alp(iN_X) )

        H_d_1(iN_E,iN_X,iS) = DT * ( B_d_2 * V_u_2(iN_X) + B_d_3 * V_u_3(iN_X) - alp(iN_X) ) & 
                            * Gm_dd_11(iN_X) * H_u_1(iN_E,iN_X,iS) &
                            - DT * ( B_d_1 * V_u_2(iN_X) *Gm_dd_22(iN_X) ) * H_u_2(iN_E,iN_X,iS)   &
                            - DT * ( B_d_1 * V_u_3(iN_X) * Gm_dd_33(iN_X) ) * H_u_3(iN_E,iN_X,iS) 

        H_d_2(iN_E,iN_X,iS) = DT * ( B_d_1 * V_u_1(iN_X) + B_d_3 * V_u_3(iN_X) - alp(iN_X) ) &
                            * Gm_dd_22(iN_X) * H_u_2(iN_E,iN_X,iS) &
                            - DT * ( B_d_2 * V_u_1(iN_X) * Gm_dd_11(iN_X) ) * H_u_1(iN_E,iN_X,iS) &
                            - DT * ( Gm_dd_33(iN_X) * H_u_3(iN_E,iN_X,iS) * B_d_2 * V_u_3(iN_X) ) 

        H_d_3(iN_E,iN_X,iS) = DT * ( B_d_1 * V_u_1(iN_X) + B_d_2 * V_u_2(iN_X) - alp(iN_X) ) &
                            * Gm_dd_33(iN_X) * H_u_3(iN_E,iN_X,iS) &
                            - DT * ( Gm_dd_11(iN_X) * H_u_1(iN_E,iN_X,iS) * B_d_3 * V_u_1(iN_X) ) &
                            - DT * ( Gm_dd_22(iN_X) * H_u_2(iN_E,iN_X,iS) * B_d_3 * V_u_2(iN_X) )


   
       vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iN_X,iS) &
               + V_u_2(iN_X) * H_d_2(iN_E,iN_X,iS) &
               + V_u_3(iN_X) * H_d_3(iN_E,iN_X,iS)

       CALL ComputeEddingtonTensorComponents_dd &
               ( J(iN_E,iN_X,iS), H_u_1(iN_E,iN_X,iS), &
                 H_u_2(iN_E,iN_X,iS), H_u_3(iN_E,iN_X,iS), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 alp(iN_X), B_u_1(iN_X), B_u_2(iN_X), B_u_3(iN_X), & 
                 V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), k_dd  )

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


        vvDotK = V_u_1(iN_X) * vDotK_d_1 &
              + V_u_2(iN_X) * vDotK_d_2 & 
              + V_u_3(iN_X) * vDotK_d_3 


      ! --- Eulerian Neutrino Number Density ---

      N_nu(iN_E,iN_X,iS) = W * J(iN_E,iN_X,iS) + vDotH

      ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

      E_nu(iN_E,iN_X,iS) = W**2 * J(iN_E,iN_X,iS) +  W * Two * vDotH + vvDotK 

      ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

      F_nu_d_1(iN_E,iN_X,iS) &
        = W * H_d_1(iN_E,iN_X,iS) + W**2 * V_d_1(iN_X) * J(iN_E,iN_X,iS) &
        + W * vDotH * V_d_1(iN_X)  + vDotK_d_1

      F_nu_d_2(iN_E,iN_X,iS) &
        = W * H_d_2(iN_E,iN_X,iS) + W**2 * V_d_2(iN_X) * J(iN_E,iN_X,iS) &
        + W * vDotH * V_d_2(iN_X)  + vDotK_d_2

      F_nu_d_3(iN_E,iN_X,iS) &
        = W * H_d_3(iN_E,iN_X,iS) + W**2 * V_d_3(iN_X) * J(iN_E,iN_X,iS) &
        + W * vDotH * V_d_3(iN_X)  + vDotK_d_3
      END IF

    END DO
    END DO
    END DO

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC  )
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP   )
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_rho  (iN_X) = Zero
        G_Y    (iN_X) = Zero
        G_Eps  (iN_X) = Zero
        G_V_d_1(iN_X) = Zero
        G_V_d_2(iN_X) = Zero
        G_V_d_3(iN_X) = Zero

      END IF
    END DO

    DO iS = iNuE, iNuE_Bar

      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, LeptonNumber(iS), N_nu(:,:,iS), &
               nE_G, W2_S, 1, One, G_Y, 1 )

    END DO

    DO iS = 1, nSpecies

      CALL MatrixVectorMultiply &
             ( 'T', nE_G, nX_G, One, E_nu    (:,:,iS), nE_G, W3_S, 1, &
               One, G_Eps  , 1 )
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

        vDotV =   V_u_1(iN_X) * V_u_1(iN_X) * Gm_dd_11(iN_X) &
                + V_u_2(iN_X) * V_u_2(iN_X) * Gm_dd_22(iN_X)&
                + V_u_3(iN_X) * V_u_3(iN_X) * Gm_dd_33(iN_X)

        W = 1.0_DP / SQRT(1.0_DP - vDotV)
       
        h = ( 1.0_DP + E(iN_X) + P(iN_X) / D(iN_X) ) 

        G_rho  (iN_X)  = ( 1.0_DP / W ) * C_rho  (iN_X) 
        G_Y    (iN_X)  = C_Y    (iN_X) - ( AtomicMassUnit / cD_old(iN_X) ) * G_Y(iN_X) * S_Y(iN_X)
        G_Eps  (iN_X)  = ( 1.0_DP / W ) * ( C_Eps(iN_X) - G_Eps(iN_X) * S_Eps(iN_X) / cD_old(iN_X) ) &
                       + S_Eps(iN_X)    * ( P(iN_X) / ( W * cD_old(iN_X) ) - ( 1.0_DP + P(iN_X) / D(iN_X) ))
      
        G_V_d_1(iN_X)  = 1.0_DP / ( h * W ) * (  C_V_d_1(iN_X) - G_V_d_1(iN_X) * S_V_d_1(iN_X) / cD_old(iN_X) ) 
        G_V_d_2(iN_X)  = 1.0_DP / ( h * W ) * (  C_V_d_2(iN_X) - G_V_d_2(iN_X) * S_V_d_2(iN_X) / cD_old(iN_X) ) 
        G_V_d_3(iN_X)  = 1.0_DP / ( h * W ) * (  C_V_d_3(iN_X) - G_V_d_3(iN_X) * S_V_d_3(iN_X) / cD_old(iN_X) ) 
        



        Gm(iD ,iN_X) = G_rho  (iN_X)
        Gm(iY ,iN_X) = G_Y    (iN_X)
        Gm(iEps,iN_X)= G_Eps  (iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)

        Fm(iD ,iN_X) = G_rho  (iN_X) - U_rho  (iN_X)
        Fm(iY ,iN_X) = G_Y    (iN_X) - U_Y    (iN_X)
        Fm(iEps ,iN_X) = G_Eps  (iN_X) - U_Eps  (iN_X)
        Fm(iV1,iN_X) = G_V_d_1(iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = G_V_d_2(iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = G_V_d_3(iN_X) - U_V_d_3(iN_X)
      END IF
    END DO
  END SUBROUTINE ComputeMatterRHS_FP


  SUBROUTINE ComputeNeutrinoRHS_FP &
    ( MASK, Fm, Gm, dt, &
      J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: alp, B_u_1, B_u_2, B_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: k_dd(3,3), vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, W, vDotV 
    REAL(DP) :: Eta_T, Chi_T, Kappa

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          Eta_T, Chi_T, Kappa, iOS )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$ACC          Eta_T, Chi_T, Kappa, iOS )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
    !$OMP          Eta_T, Chi_T, Kappa, iOS )
#endif
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iN_X,iS) &
                + V_u_2(iN_X) * H_d_2(iN_E,iN_X,iS) &
                + V_u_3(iN_X) * H_d_3(iN_E,iN_X,iS)

       CALL ComputeEddingtonTensorComponents_dd &
               ( J(iN_E,iN_X,iS), H_u_1(iN_E,iN_X,iS), &
                 H_u_2(iN_E,iN_X,iS), H_u_3(iN_E,iN_X,iS), &
                 Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X), &
                 alp(iN_X), B_u_1(iN_X), B_u_2(iN_X), B_u_3(iN_X), & 
                 V_u_1(iN_X), V_u_2(iN_X), V_u_3(iN_X), k_dd  )

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

        vDotV =   V_u_1(iN_X) * V_u_1(iN_X) * Gm_dd_11(iN_X) &
                + V_u_2(iN_X) * V_u_2(iN_X) * Gm_dd_22(iN_X)&
                + V_u_3(iN_X) * V_u_3(iN_X) * Gm_dd_33(iN_X)

        W = 1.0_DP / SQRT(1.0_DP - vDotV) 
        ! --- Emissivity ---

        Eta_T =   Chi     (iN_E,iN_X,iS) * J0(iN_E,iN_X,iS) &
                + Eta_NES (iN_E,iN_X,iS) &
                + Eta_Pair(iN_E,iN_X,iS)

        ! --- Number Opacity ---

        Chi_T =   Chi     (iN_E,iN_X,iS) &
                + Chi_NES (iN_E,iN_X,iS) &
                + Chi_Pair(iN_E,iN_X,iS)

        ! --- Number Flux Opacity ---

        Kappa = Chi_T + Sig(iN_E,iN_X,iS)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! --- Number Equation ---

        Gm(iOS+iCR_N,iN_X) &
          = ( One - Omega(iN_X) ) *     J(iN_E,iN_X,iS) &
           + Omega(iN_X)   *  ( C_J(iN_E,iN_X,iS) - vDotH + dt * Eta_T ) &
              / ( W + dt * Chi_T )

        Fm(iOS+iCR_N,iN_X) &
          = Gm(iOS+iCR_N,iN_X) - J(iN_E,iN_X,iS)

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( One - Omega(iN_X) ) *     H_d_1(iN_E,iN_X,iS) &
            +       Omega(iN_X)   * ( C_H_d_1(iN_E,iN_X,iS) - vDotK_d_1 ) &
                                    / ( W + dt * Kappa )

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - H_d_1(iN_E,iN_X,iS)

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( One - Omega(iN_X) ) *     H_d_2(iN_E,iN_X,iS) &
            +       Omega(iN_X)   * ( C_H_d_2(iN_E,iN_X,iS) - vDotK_d_2 ) &
                                    / ( W + dt * Kappa )

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - H_d_2(iN_E,iN_X,iS)

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( One - Omega(iN_X) ) *     H_d_3(iN_E,iN_X,iS) &
            +       Omega(iN_X)   * ( C_H_d_3(iN_E,iN_X,iS) - vDotK_d_3 ) &
                                    / ( W + dt * Kappa )

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - H_d_3(iN_E,iN_X,iS)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_FP


  SUBROUTINE UpdateMatterRHS_FP &
    ( MASK, Fm, Gm, D, Y, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: D, Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_X
    REAL(DP) :: vDotV, W

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( vDotV )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( vDotV )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( vDotV )
#endif
    DO iN_X = 1, nX_G

      IF ( MASK(iN_X) ) THEN

        Fm(iD ,iN_X) = Gm(iD ,iN_X) - U_rho  (iN_X)
        Fm(iY ,iN_X) = Gm(iY ,iN_X) - U_Y    (iN_X)
        Fm(iEps,iN_X)= Gm(iEps ,iN_X) - U_Eps  (iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)
        
        U_rho  (iN_X) = Gm(iD ,iN_X)
        U_Y    (iN_X) = Gm(iY ,iN_X)
        U_Eps  (iN_X) = Gm(iEps ,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        D    (iN_X) = U_rho  (iN_X) * rho_old(iN_X)
        Y    (iN_X) = U_Y    (iN_X) * Y_old  (iN_X)
        E    (iN_X) = U_Eps  (iN_X) * Eps_old(iN_X)
        V_d_1(iN_X) = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2(iN_X) = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3(iN_X) = U_V_d_3(iN_X) * SpeedOfLight

        ! --- Update Three-Velocity (Index Up) ---

        V_u_1(iN_X) = V_d_1(iN_X) / Gm_dd_11(iN_X)
        V_u_2(iN_X) = V_d_2(iN_X) / Gm_dd_22(iN_X)
        V_u_3(iN_X) = V_d_3(iN_X) / Gm_dd_33(iN_X)

        ! --- Compute E from Ef ---

        vDotV =   V_u_1(iN_X) * V_d_1(iN_X) &
                + V_u_2(iN_X) * V_d_2(iN_X) &
                + V_u_3(iN_X) * V_d_3(iN_X)
        W = 1.0_DP / SQRT( 1.0_DP - vDotV)

        ! --- Compute Omega (Richardson damping coeff.) based on velocity ---

        Omega(iN_X) = One / ( W + SQRT( vDotV ) )
      END IF

    END DO

  END SUBROUTINE UpdateMatterRHS_FP


  SUBROUTINE UpdateNeutrinoRHS_FP &
    ( MASK, Fm, Gm, J, H_u_1, H_u_2, H_u_3, &
      V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B_u_1, B_u_2, B_u_3 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: alp, B_u_1, B_u_2, B_u_3

    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: B_d_1, B_d_2, B_d_3

 
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
    DO iS   = 1, nSpecies
    DO iN_X = 1, nX_G
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        Fm(iOS+iCR_N ,iN_X) = Gm(iOS+iCR_N ,iN_X) - J    (iN_E,iN_X,iS)
        Fm(iOS+iCR_G1,iN_X) = Gm(iOS+iCR_G1,iN_X) - H_d_1(iN_E,iN_X,iS)
        Fm(iOS+iCR_G2,iN_X) = Gm(iOS+iCR_G2,iN_X) - H_d_2(iN_E,iN_X,iS)
        Fm(iOS+iCR_G3,iN_X) = Gm(iOS+iCR_G3,iN_X) - H_d_3(iN_E,iN_X,iS)

        J    (iN_E,iN_X,iS) = Gm(iOS+iCR_N ,iN_X)
        H_d_1(iN_E,iN_X,iS) = Gm(iOS+iCR_G1,iN_X)
        H_d_2(iN_E,iN_X,iS) = Gm(iOS+iCR_G2,iN_X)
        H_d_3(iN_E,iN_X,iS) = Gm(iOS+iCR_G3,iN_X)

        B_d_1 = Gm_dd_11(iN_X) * B_u_1(iN_X)
        B_d_2 = Gm_dd_22(iN_X) * B_u_2(iN_X)
        B_d_3 = Gm_dd_33(iN_X) * B_u_3(iN_X)

        H_u_1(iN_E,iN_X,iS) = ( 1.0_DP - B_d_1 * V_u_1(iN_X) / alp(iN_X) ) * H_d_1(iN_E,iN_X,iS) / Gm_dd_11(iN_X)  &
                 - H_d_2(iN_E,iN_X,iS) * B_d_1 * V_u_2(iN_X) / ( alp(iN_X) *Gm_dd_11(iN_X) ) &
                 - H_d_3(iN_E,iN_X,iS) * B_d_1 * V_u_3(iN_X) / ( Gm_dd_11(iN_X) * alp(iN_X) )

        H_u_2(iN_E,iN_X,iS) =( 1.0_DP -   B_d_2 * V_u_2(iN_X) / alp(iN_X) ) * H_d_2(iN_E,iN_X,iS) / Gm_dd_22(iN_X)  &
                  - H_d_1(iN_E,iN_X,iS) * B_d_2 * V_u_1(iN_X) / ( alp(iN_X) *Gm_dd_22(iN_X) ) &
                  - H_d_3(iN_E,iN_X,iS) * B_d_2 * V_u_3(iN_X) / ( Gm_dd_22(iN_X) * alp(iN_X) )

        H_u_3(iN_E,iN_X,iS) =( 1.0_DP - B_d_3 * V_u_3(iN_X) / alp(iN_X) ) * H_d_3(iN_E,iN_X,iS) / Gm_dd_33(iN_X)  &
                  - H_d_1(iN_E,iN_X,iS) * B_d_3 * V_u_1(iN_X) / ( alp(iN_X) *Gm_dd_33(iN_X) ) &
                  - H_d_2(iN_E,iN_X,iS) * B_d_3 * V_u_2(iN_X) / ( Gm_dd_33(iN_X) * alp(iN_X) )




      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE UpdateNeutrinoRHS_FP


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

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F, G

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
    ( MASK, n_FP, k_inner, nIterations_Inner, Fm )

    LOGICAL,  DIMENSION(:)    , INTENT(inout) :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, k_inner
    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Inner
    REAL(DP), DIMENSION(:,:)  , INTENT(in)    :: Fm

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X, iS, iN_E, iOS
    REAL(DP) :: Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3

    ! --- Compute Norm of F (Use Maximum Across Moments) ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, iOS, Fnorm, Fnorm_N, Fnorm_G1, Fnorm_G2, Fnorm_G3 )
#endif
    DO iN_X = 1, nX_G

      IF( MASK(iN_X) )THEN

        CONVERGED = .TRUE.
 
        DO iS   = 1, nSpecies

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

          CONVERGED = CONVERGED .AND. ( Fnorm <= Rtol_inner * Jnorm(iS,iN_X) )

        END DO

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
    ( MASK_OUTER, MASK_INNER, n_FP, k_outer, nIterations_Outer, Fm )

    LOGICAL,  DIMENSION(:)    , INTENT(inout) :: MASK_OUTER, MASK_INNER
    INTEGER,                    INTENT(in)    :: n_FP, k_outer
    INTEGER,  DIMENSION(:)    , INTENT(inout) :: nIterations_Outer
    REAL(DP), DIMENSION(:,:)  , INTENT(in)    :: Fm

    LOGICAL  :: CONVERGED
    INTEGER  :: iN_X
    REAL(DP) :: Fnorm_Y, Fnorm_E, Fnorm_V, Fnorm_D

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_Ef, Fnorm_V )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_Y, Fnorm_Ef, Fnorm_V )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_Ef, Fnorm_V )
#endif
    DO iN_X = 1, nX_G
      IF( MASK_OUTER(iN_X) )THEN


        Fnorm_D  =      ABS( Fm(iD ,iN_X) )
        Fnorm_Y  =      ABS( Fm(iY ,iN_X) )
        Fnorm_E  =      ABS( Fm(iEps ,iN_X) )
        Fnorm_V  = MAX( ABS( Fm(iV1,iN_X) ), &
                        ABS( Fm(iV2,iN_X) ), &
                        ABS( Fm(iV3,iN_X) ) )

        CONVERGED = Fnorm_Y  <= Rtol_outer .AND. &
                    Fnorm_E  <= Rtol_outer .AND. &
                    Fnorm_D  <= Rtol_outer .AND. &
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


END MODULE TwoMoment_NeutrinoMatterSolverModule_Relativistic
