#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_NeutrinoMatterSolverModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, FourPi, SqrtTiny
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
    Timer_Collisions_OuterLoop, &
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
    ComputeEquilibriumDistributions_DG, &
    ComputeNeutrinoOpacities_EC, &
    ComputeNeutrinoOpacities_ES, &
    ComputeNeutrinoOpacities_NES, &
    ComputeNeutrinoOpacities_Pair, &
    ComputeNeutrinoOpacities_Brem, &
    ComputeNeutrinoOpacityRates_NES, &
    ComputeNeutrinoOpacityRates_Pair, &
    ComputeNeutrinoOpacityRates_Brem
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    EddingtonTensorComponents_dd

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: InitializeNeutrinoMatterSolverParameters
  PUBLIC :: SolveNeutrinoMatterCoupling_FP_Nested_AA

  LOGICAL, PARAMETER :: Include_NES  = .TRUE.
  LOGICAL, PARAMETER :: Include_Pair = .TRUE.
  LOGICAL, PARAMETER :: Include_Brem = .TRUE.

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

  INTEGER  :: nE_G, nX_G, nZ(4)
  INTEGER  :: n_FP_inner, n_FP_outer
  INTEGER  :: iE_B, iE_E

  REAL(DP), ALLOCATABLE :: E_N(:)  ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:) ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:) ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:) ! --- Integration Weights Scaled by (hc)^3
  REAL(DP), ALLOCATABLE :: W3_S(:) ! --- Integration Weights Scaled by (hc)^3
  REAL(DP), ALLOCATABLE :: FourPiEp2(:)

  ! --- Solver scratch arrays ---

  LOGICAL,  DIMENSION(:), ALLOCATABLE :: ITERATE_outer, ITERATE_inner
  INTEGER,  DIMENSION(:), ALLOCATABLE :: PackIndex_outer, UnpackIndex_outer
  INTEGER,  DIMENSION(:), ALLOCATABLE :: PackIndex_inner, UnpackIndex_inner

  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE :: Jnorm
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_J
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_H_d_1, H_d_1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_H_d_2, H_d_2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C_H_d_3, H_d_3

  REAL(DP), DIMENSION(:), ALLOCATABLE :: Omega
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Y_old, C_Y, S_Y, G_Y, U_Y
  REAL(DP), DIMENSION(:), ALLOCATABLE :: C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J0
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Sigma_Iso
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_EmAb, Eta_EmAb
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NES, Eta_NES
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Pair, Eta_Pair
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Brem, Eta_Brem

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_0, H_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_0, J_II_0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: S_Sigma

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
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wMatterRHS

  INTEGER, PARAMETER :: iY  = 1
  INTEGER, PARAMETER :: iEf = 2
  INTEGER, PARAMETER :: iV1 = 3
  INTEGER, PARAMETER :: iV2 = 4
  INTEGER, PARAMETER :: iV3 = 5

  ! --- Temporary arrays for scatter/gather (packing)

  REAL(DP), DIMENSION(:)    , ALLOCATABLE, TARGET :: D_T, T_T, Y_T, E_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J0_T
  REAL(DP), DIMENSION(:,:)  , ALLOCATABLE, TARGET :: Sigma_Iso_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_EmAb_T, Eta_EmAb_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_NES_T, Eta_NES_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Pair_T, Eta_Pair_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: Chi_Brem_T, Eta_Brem_T

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: H_I_0_T, H_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: J_I_0_T, J_II_0_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: S_Sigma_T

CONTAINS


  SUBROUTINE InitializeNeutrinoMatterSolver( iZ_B, iZ_E )

    INTEGER, INTENT(in) :: iZ_B(4), iZ_E(4)

    REAL(DP) :: TMP(1)

    INTEGER :: iE1, iE2, iN_X

    iE_B = iZ_B(1)
    iE_E = iZ_E(1)
    nZ   = iZ_E - iZ_B + 1

    nE_G = nZ(1) * nNodesZ(1)
    nX_G = PRODUCT( nZ(2:4) * nNodesZ(2:4) )

    n_FP_outer = 5
    n_FP_inner = nE_G * nCR * nSpecies

    ALLOCATE( wMatterRHS(n_FP_outer) )

    wMatterRHS(iY ) = One  ! --- One = On, Zero = Off
    wMatterRHS(iEf) = One
    wMatterRHS(iV1) = One
    wMatterRHS(iV2) = One
    wMatterRHS(iV3) = One

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )
    ALLOCATE( FourPiEp2(nE_G) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S      = W2_N / ( PlanckConstant * SpeedOfLight )**3
    W3_S      = W3_N / ( PlanckConstant * SpeedOfLight )**3
    FourPiEp2 = FourPi * E_N * E_N

    ALLOCATE(     ITERATE_outer(nX_G) )
    ALLOCATE(   PackIndex_outer(nX_G) )
    ALLOCATE( UnpackIndex_outer(nX_G) )

    ALLOCATE(     ITERATE_inner(nX_G) )
    ALLOCATE(   PackIndex_inner(nX_G) )
    ALLOCATE( UnpackIndex_inner(nX_G) )

    ITERATE_outer = .TRUE.
    ITERATE_inner = .TRUE.

    ALLOCATE( Jnorm(nSpecies,nX_G) )

    ALLOCATE( C_J(nE_G,nSpecies,nX_G) )

    ALLOCATE( C_H_d_1(nE_G,nSpecies,nX_G) )
    ALLOCATE(   H_d_1(nE_G,nSpecies,nX_G) )

    ALLOCATE( C_H_d_2(nE_G,nSpecies,nX_G) )
    ALLOCATE(   H_d_2(nE_G,nSpecies,nX_G) )

    ALLOCATE( C_H_d_3(nE_G,nSpecies,nX_G) )
    ALLOCATE(   H_d_3(nE_G,nSpecies,nX_G) )

    ALLOCATE( Omega(nX_G) )

    ALLOCATE(   Ef_old(nX_G) )
    ALLOCATE( C_Ef    (nX_G) )
    ALLOCATE( S_Ef    (nX_G) )
    ALLOCATE( G_Ef    (nX_G) )
    ALLOCATE( U_Ef    (nX_G) )
    ALLOCATE(   Ef    (nX_G) )

    ALLOCATE(   Y_old(nX_G) )
    ALLOCATE( C_Y    (nX_G) )
    ALLOCATE( S_Y    (nX_G) )
    ALLOCATE( G_Y    (nX_G) )
    ALLOCATE( U_Y    (nX_G) )

    ALLOCATE( C_V_d_1(nX_G) )
    ALLOCATE( S_V_d_1(nX_G) )
    ALLOCATE( G_V_d_1(nX_G) )
    ALLOCATE( U_V_d_1(nX_G) )
    ALLOCATE(   V_d_1(nX_G) )

    ALLOCATE( C_V_d_2(nX_G) )
    ALLOCATE( S_V_d_2(nX_G) )
    ALLOCATE( G_V_d_2(nX_G) )
    ALLOCATE( U_V_d_2(nX_G) )
    ALLOCATE(   V_d_2(nX_G) )

    ALLOCATE( C_V_d_3(nX_G) )
    ALLOCATE( S_V_d_3(nX_G) )
    ALLOCATE( G_V_d_3(nX_G) )
    ALLOCATE( U_V_d_3(nX_G) )
    ALLOCATE(   V_d_3(nX_G) )

    ALLOCATE(        J0(nE_G,nSpecies,nX_G) )
    ALLOCATE( Sigma_Iso(nE_G,         nX_G) )
    ALLOCATE(  Chi_EmAb(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Eta_EmAb(nE_G,nSpecies,nX_G) )
    ALLOCATE(   Chi_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(   Eta_NES(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Chi_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Eta_Pair(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Chi_Brem(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Eta_Brem(nE_G,nSpecies,nX_G) )

    ALLOCATE(   H_I_0(nE_G,nE_G,nX_G) )
    ALLOCATE(  H_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE(   J_I_0(nE_G,nE_G,nX_G) )
    ALLOCATE(  J_II_0(nE_G,nE_G,nX_G) )
    ALLOCATE( S_Sigma(nE_G,nE_G,nX_G) )

    ALLOCATE( D_T(nX_G) )
    ALLOCATE( T_T(nX_G) )
    ALLOCATE( Y_T(nX_G) )
    ALLOCATE( E_T(nX_G) )

    ALLOCATE( J_T(nE_G,nSpecies,nX_G) )

    ALLOCATE(        J0_T(nE_G,nSpecies,nX_G) )
    ALLOCATE( Sigma_Iso_T(nE_G,         nX_G) )
    ALLOCATE(  Chi_EmAb_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Eta_EmAb_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(   Chi_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(   Eta_NES_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Chi_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Eta_Pair_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Chi_Brem_T(nE_G,nSpecies,nX_G) )
    ALLOCATE(  Eta_Brem_T(nE_G,nSpecies,nX_G) )

    ALLOCATE(   H_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(  H_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(   J_I_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE(  J_II_0_T(nE_G,nE_G,nX_G) )
    ALLOCATE( S_Sigma_T(nE_G,nE_G,nX_G) )

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

#if   defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatterRHS, &
    !$OMP          ITERATE_outer, ITERATE_inner ) &
    !$OMP MAP( alloc: INFO, &
    !$OMP             Jnorm, &
    !$OMP             C_J, &
    !$OMP             C_H_d_1, H_d_1, &
    !$OMP             C_H_d_2, H_d_2, &
    !$OMP             C_H_d_3, H_d_3, &
    !$OMP             Omega, &
    !$OMP             Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$OMP             Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$OMP             C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$OMP             C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$OMP             C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$OMP             J0, Sigma_Iso, &
    !$OMP             Chi_EmAb, Eta_EmAb, &
    !$OMP             Chi_NES, Eta_NES, &
    !$OMP             Chi_Pair, Eta_Pair, &
    !$OMP             Chi_Brem, Eta_Brem, &
    !$OMP             H_I_0, H_II_0, J_I_0, J_II_0, S_Sigma, &
    !$OMP             D_T, T_T, Y_T, E_T, J_T, &
    !$OMP             J0_T, Sigma_Iso_T, &
    !$OMP             Chi_EmAb_T, Eta_EmAb_T, &
    !$OMP             Chi_NES_T, Eta_NES_T, &
    !$OMP             Chi_Pair_T, Eta_Pair_T, &
    !$OMP             Chi_Brem_T, Eta_Brem_T, &
    !$OMP             H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, S_Sigma_T, &
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
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatterRHS, &
    !$ACC         ITERATE_outer, ITERATE_inner ) &
    !$ACC CREATE( INFO, &
    !$ACC         Jnorm, &
    !$ACC         C_J, &
    !$ACC         C_H_d_1, H_d_1, &
    !$ACC         C_H_d_2, H_d_2, &
    !$ACC         C_H_d_3, H_d_3, &
    !$ACC         Omega, &
    !$ACC         Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$ACC         Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$ACC         C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$ACC         C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$ACC         C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$ACC         J0, Sigma_Iso, &
    !$ACC         Chi_EmAb, Eta_EmAb, &
    !$ACC         Chi_NES, Eta_NES, &
    !$ACC         Chi_Pair, Eta_Pair, &
    !$ACC         Chi_Brem, Eta_Brem, &
    !$ACC         H_I_0, H_II_0, J_I_0, J_II_0, S_Sigma, &
    !$ACC         D_T, T_T, Y_T, E_T, J_T, &
    !$ACC         J0_T, Sigma_Iso_T, &
    !$ACC         Chi_EmAb_T, Eta_EmAb_T, &
    !$ACC         Chi_NES_T, Eta_NES_T, &
    !$ACC         Chi_Pair_T, Eta_Pair_T, &
    !$ACC         Chi_Brem_T, Eta_Brem_T, &
    !$ACC         H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, S_Sigma_T, &
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

    IF ( .NOT. Include_NES ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X  = 1, nX_G
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        H_I_0   (iE1,iE2,iN_X) = Zero
        H_II_0  (iE1,iE2,iN_X) = Zero
        H_I_0_T (iE1,iE2,iN_X) = Zero
        H_II_0_T(iE1,iE2,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

    IF ( .NOT. Include_Pair ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X  = 1, nX_G
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        J_I_0   (iE1,iE2,iN_X) = Zero
        J_II_0  (iE1,iE2,iN_X) = Zero
        J_I_0_T (iE1,iE2,iN_X) = Zero
        J_II_0_T(iE1,iE2,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

    IF ( .NOT. Include_Brem ) THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iN_X  = 1, nX_G
      DO iE2 = 1, nE_G
      DO iE1 = 1, nE_G

        S_Sigma  (iE1,iE2,iN_X) = Zero
        S_Sigma_T(iE1,iE2,iN_X) = Zero

      END DO
      END DO
      END DO

    END IF

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
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatterRHS, &
    !$OMP               INFO, &
    !$OMP               Jnorm, &
    !$OMP               C_J, &
    !$OMP               C_H_d_1, H_d_1, &
    !$OMP               C_H_d_2, H_d_2, &
    !$OMP               C_H_d_3, H_d_3, &
    !$OMP               Omega, &
    !$OMP               Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$OMP               Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$OMP               C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$OMP               C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$OMP               C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$OMP               J0, Sigma_Iso, &
    !$OMP               Chi_EmAb, Eta_EmAb, &
    !$OMP               Chi_NES, Eta_NES, &
    !$OMP               Chi_Pair, Eta_Pair, &
    !$OMP               Chi_Brem, Eta_Brem, &
    !$OMP               H_I_0, H_II_0, J_I_0, J_II_0, S_Sigma, &
    !$OMP               D_T, T_T, Y_T, E_T, J_T, &
    !$OMP               J0_T, Sigma_Iso_T, &
    !$OMP               Chi_EmAb_T, Eta_EmAb_T, &
    !$OMP               Chi_NES_T, Eta_NES_T, &
    !$OMP               Chi_Pair_T, Eta_Pair_T, &
    !$OMP               Chi_Brem_T, Eta_Brem_T, &
    !$OMP               H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, S_Sigma_T, &
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
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2, wMatterRHS, &
    !$ACC         INFO, &
    !$ACC         Jnorm, &
    !$ACC         C_J, &
    !$ACC         C_H_d_1, H_d_1, &
    !$ACC         C_H_d_2, H_d_2, &
    !$ACC         C_H_d_3, H_d_3, &
    !$ACC         Omega, &
    !$ACC         Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef, &
    !$ACC         Y_old, C_Y, S_Y, G_Y, U_Y, &
    !$ACC         C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1, &
    !$ACC         C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2, &
    !$ACC         C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3, &
    !$ACC         J0, Sigma_Iso, &
    !$ACC         Chi_EmAb, Eta_EmAb, &
    !$ACC         Chi_NES, Eta_NES, &
    !$ACC         Chi_Pair, Eta_Pair, &
    !$ACC         Chi_Brem, Eta_Brem, &
    !$ACC         H_I_0, H_II_0, J_I_0, J_II_0, S_Sigma, &
    !$ACC         D_T, T_T, Y_T, E_T, J_T, &
    !$ACC         J0_T, Sigma_Iso_T, &
    !$ACC         Chi_EmAb_T, Eta_EmAb_T, &
    !$ACC         Chi_NES_T, Eta_NES_T, &
    !$ACC         Chi_Pair_T, Eta_Pair_T, &
    !$ACC         Chi_Brem_T, Eta_Brem_T, &
    !$ACC         H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, S_Sigma_T, &
    !$ACC         ITERATE_outer, PackIndex_outer, UnpackIndex_outer, &
    !$ACC         AMAT_outer, GVEC_outer, FVEC_outer, &
    !$ACC         BVEC_outer, GVECm_outer, FVECm_outer, &
    !$ACC         WORK_outer, TAU_outer, Alpha_outer, &
    !$ACC         ITERATE_inner, PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_inner, GVEC_inner, FVEC_inner, &
    !$ACC         BVEC_inner, GVECm_inner, FVECm_inner, &
    !$ACC         WORK_inner, TAU_inner, Alpha_inner )
#endif

    DEALLOCATE( wMatterRHS )
    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S, FourPiEp2 )
    DEALLOCATE( INFO )
    DEALLOCATE( Jnorm )
    DEALLOCATE( C_J )
    DEALLOCATE( C_H_d_1, H_d_1 )
    DEALLOCATE( C_H_d_2, H_d_2 )
    DEALLOCATE( C_H_d_3, H_d_3 )
    DEALLOCATE( Omega )
    DEALLOCATE( Ef_old, C_Ef, S_Ef, G_Ef, U_Ef, Ef )
    DEALLOCATE( Y_old, C_Y, S_Y, G_Y, U_Y )
    DEALLOCATE( C_V_d_1, S_V_d_1, G_V_d_1, U_V_d_1, V_d_1 )
    DEALLOCATE( C_V_d_2, S_V_d_2, G_V_d_2, U_V_d_2, V_d_2 )
    DEALLOCATE( C_V_d_3, S_V_d_3, G_V_d_3, U_V_d_3, V_d_3 )
    DEALLOCATE( J0, Sigma_Iso )
    DEALLOCATE( Chi_EmAb, Eta_EmAb )
    DEALLOCATE( Chi_NES, Eta_NES )
    DEALLOCATE( Chi_Pair, Eta_Pair )
    DEALLOCATE( Chi_Brem, Eta_Brem )
    DEALLOCATE( H_I_0, H_II_0, J_I_0, J_II_0, S_Sigma )
    DEALLOCATE( D_T, T_T, Y_T, E_T, J_T )
    DEALLOCATE( J0_T, Sigma_Iso_T )
    DEALLOCATE( Chi_EmAb_T, Eta_EmAb_T )
    DEALLOCATE( Chi_NES_T, Eta_NES_T )
    DEALLOCATE( Chi_Pair_T, Eta_Pair_T )
    DEALLOCATE( Chi_Brem_T, Eta_Brem_T )
    DEALLOCATE( H_I_0_T, H_II_0_T, J_I_0_T, J_II_0_T, S_Sigma_T )
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

      W2(iN_E) = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**2
      W3(iN_E) = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**3

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveNeutrinoMatterCoupling_FP_Nested_AA &
    ( dt, J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, D, T, Y, E, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, nIterations_Inner, nIterations_Outer )

    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:),     INTENT(inout) :: D, T, Y, E
    REAL(DP), DIMENSION(:),     INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER,  DIMENSION(:),     INTENT(out)   :: nIterations_Inner
    INTEGER,  DIMENSION(:),     INTENT(out)   :: nIterations_Outer

    ! --- Local Variables ---

    INTEGER  :: k_outer, Mk_outer, nX_P_outer
    INTEGER  :: k_inner, Mk_inner, nX_P_inner

    ! --- Initial RHS ---

    CALL TimersStart( Timer_Collisions_InitializeRHS )

    CALL InitializeRHS_FP &
           ( J, H_u_1, H_u_2, H_u_3, D, Y, E, V_u_1, V_u_2, V_u_3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    CALL TimersStop( Timer_Collisions_InitializeRHS )

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Collisions_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( D, T, Y )

    CALL TimersStop( Timer_Collisions_ComputeOpacity )

    ! --- Start Outer Loop ---

    CALL TimersStart( Timer_Collisions_OuterLoop )

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

        CALL TimersStart( Timer_Collisions_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( D, T, Y, &
                 ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        CALL TimersStop( Timer_Collisions_ComputeOpacity )

      END IF

      ! --- Start Inner Loop ---

      CALL TimersStart( Timer_Collisions_InnerLoop )

      k_inner = 0
      DO WHILE( ANY( ITERATE_inner(:) ) .AND. k_inner < MaxIter_inner )

        k_inner  = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_inner, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

        CALL TimersStart( Timer_Collisions_ComputeRates )

        CALL ComputeRates_Packed &
               ( J, ITERATE_inner, nX_P_inner, PackIndex_inner, &
                 UnpackIndex_inner, nX_P_outer )

        CALL TimersStop( Timer_Collisions_ComputeRates )

        ! --- Right-Hand Side Vectors and Residuals (inner) ---

        CALL TimersStart( Timer_Collisions_NeutrinoRHS )

        CALL ComputeNeutrinoRHS_FP &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, dt, &
                 J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        CALL TimersStop( Timer_Collisions_NeutrinoRHS )

        ! --- Anderson Acceleration (inner) ---

        CALL TimersStart( Timer_Collisions_SolveLS )

        CALL SolveLS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner, &
                 AMAT_inner, BVEC_inner, Alpha_inner, TAU_inner, &
                 LWORK_inner, WORK_inner )

        CALL TimersStop( Timer_Collisions_SolveLS )

        ! --- Update Residuals and Solution Vectors (inner) ---

        CALL TimersStart( Timer_Collisions_UpdateFP )

        CALL UpdateNeutrinoRHS_FP &
               ( ITERATE_inner, FVECm_inner, GVECm_inner, &
                 J, H_u_1, H_u_2, H_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        ! --- Check Convergence (inner) ---

        CALL TimersStart( Timer_Collisions_CheckInner )

        CALL CheckConvergence_Inner &
               ( ITERATE_inner, n_FP_inner, k_inner, &
                 nIterations_Inner, FVECm_inner )

        CALL TimersStop( Timer_Collisions_CheckInner )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_inner, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Collisions_UpdateFP )

      END DO ! --- Inner Loop ---

      CALL TimersStop( Timer_Collisions_InnerLoop )

      ! --- Right-Hand Side Vectors and Residuals (outer) ---

      CALL TimersStart( Timer_Collisions_MatterRHS )

      CALL ComputeMatterRHS_FP &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      CALL TimersStop( Timer_Collisions_MatterRHS )

      ! --- Anderson Acceleration (outer) ---

      CALL TimersStart( Timer_Collisions_SolveLS )

      CALL SolveLS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVECm_outer, GVECm_outer, FVEC_outer, GVEC_outer, &
               AMAT_outer, BVEC_outer, Alpha_outer, TAU_outer, &
               LWORK_outer, WORK_outer )

      CALL TimersStop( Timer_Collisions_SolveLS )

      ! --- Update Residuals and Solution Vectors (outer) ---

      CALL TimersStart( Timer_Collisions_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE_outer, FVECm_outer, GVECm_outer, &
               Y, E, V_u_1, V_u_2, V_u_3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      ! --- Check Convergence (outer) ---

      CALL TimersStart( Timer_Collisions_CheckOuter )

      CALL CheckConvergence_Outer &
             ( ITERATE_outer, ITERATE_inner, n_FP_outer, k_outer, &
               nIterations_Outer, FVECm_outer )

      CALL TimersStop( Timer_Collisions_CheckOuter )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_outer, n_FP_outer, M_outer, Mk_outer, &
               FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Collisions_UpdateFP )

    END DO ! --- Outer Loop ---

    CALL TimersStop( Timer_Collisions_OuterLoop )

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
    REAL(DP), DIMENSION(:,:,:), POINTER :: J0_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Sigma_Iso_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_EmAb_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_0_P, H_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_0_P, J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: S_Sigma_P

    INTEGER :: nX, nX0, iX, iE

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

      D_P => D_T(1:nX)
      T_P => T_T(1:nX)
      Y_P => Y_T(1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, D, T, Y, D_P, T_P, Y_P )

      J0_P        => J0_T       (:,:,1:nX)
      Sigma_Iso_P => Sigma_Iso_T(  :,1:nX)
      Chi_EmAb_P  => Chi_EmAb_T (:,:,1:nX)

      H_I_0_P     => H_I_0_T    (:,:,1:nX)
      H_II_0_P    => H_II_0_T   (:,:,1:nX)

      J_I_0_P     => J_I_0_T    (:,:,1:nX)
      J_II_0_P    => J_II_0_T   (:,:,1:nX)

      S_Sigma_P   => S_Sigma_T  (:,:,1:nX)

    ELSE

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      J0_P        => J0       (:,:,:)
      Sigma_Iso_P => Sigma_Iso(  :,:)
      Chi_EmAb_P  => Chi_EmAb (:,:,:)

      H_I_0_P     => H_I_0    (:,:,:)
      H_II_0_P    => H_II_0   (:,:,:)

      J_I_0_P     => J_I_0    (:,:,:)
      J_II_0_P    => J_II_0   (:,:,:)

      S_Sigma_P   => S_Sigma  (:,:,:)

    END IF

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_DG &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Y_P, J0_P )

    ! --- EmAb ---

    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE_G, 1, nSpecies, 1, nX, E_N, D_P, T_P, Y_P, Chi_EmAb_P )

    ! --- Isoenergetic scattering ---

    CALL ComputeNeutrinoOpacities_ES &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, 1, Sigma_Iso_P )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX = 1, nX
    DO iE = 1, nE_G

      Sigma_Iso_P(iE,iX) = FourPiEp2(iE) * Sigma_Iso_P(iE,iX)

    END DO
    END DO

    IF( Include_NES )THEN

      ! --- NES Scattering Functions ---

      CALL ComputeNeutrinoOpacities_NES &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, 1, H_I_0_P, H_II_0_P )


    END IF

    IF( Include_Pair )THEN

      ! --- Pair Kernels ---

      CALL ComputeNeutrinoOpacities_Pair &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, 1, J_I_0_P, J_II_0_P )


    END IF

    IF( Include_Brem )THEN

      ! --- Brem Kernels ---

      CALL ComputeNeutrinoOpacities_Brem &
             ( 1, nE_G, 1, nX, D_P, T_P, Y_P, S_Sigma_P )


    END IF

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               J0_P, J0 )
      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_EmAb_P, Chi_EmAb )
      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Sigma_Iso_P, Sigma_Iso )

      IF ( nX < nX0 ) THEN

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 H_I_0_P, H_II_0_P, H_I_0, H_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 J_I_0_P, J_II_0_P, J_I_0, J_II_0 )

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 S_Sigma_P, S_Sigma )

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

    REAL(DP), DIMENSION(:,:,:), POINTER :: J_P, J0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_NES_P, Eta_NES_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_Pair_P, Eta_Pair_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Chi_Brem_P, Eta_Brem_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: H_I_0_P, H_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: J_I_0_P, J_II_0_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: S_Sigma_P

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

      J_P        => J_T       (:,:,1:nX)
      J0_P       => J0_T      (:,:,1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, J, J0, J_P, J0_P )

      Chi_NES_P  => Chi_NES_T (:,:,1:nX)
      Eta_NES_P  => Eta_NES_T (:,:,1:nX)
      Chi_Pair_P => Chi_Pair_T(:,:,1:nX)
      Eta_Pair_P => Eta_Pair_T(:,:,1:nX)
      Chi_Brem_P => Chi_Brem_T(:,:,1:nX)
      Eta_Brem_P => Eta_Brem_T(:,:,1:nX)

      H_I_0_P    => H_I_0_T   (:,:,1:nX)
      H_II_0_P   => H_II_0_T  (:,:,1:nX)

      J_I_0_P    => J_I_0_T   (:,:,1:nX)
      J_II_0_P   => J_II_0_T  (:,:,1:nX)

      S_Sigma_P  => S_Sigma_T (:,:,1:nX)

      IF ( nX < nX0 ) THEN

        CALL ArrayPack &
               ( nX, UnpackIndex, &
                 H_I_0, H_II_0, H_I_0_P, H_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, &
                 J_I_0, J_II_0, J_I_0_P, J_II_0_P )

        CALL ArrayPack &
               ( nX, UnpackIndex, &
                 S_Sigma, S_Sigma_P )

      END IF

    ELSE

      J_P        => J       (:,:,:)
      J0_P       => J0      (:,:,:)

      Chi_NES_P  => Chi_NES (:,:,:)
      Eta_NES_P  => Eta_NES (:,:,:)
      Chi_Pair_P => Chi_Pair(:,:,:)
      Eta_Pair_P => Eta_Pair(:,:,:)
      Chi_Brem_P => Chi_Brem(:,:,:)
      Eta_Brem_P => Eta_Brem(:,:,:)

      H_I_0_P    => H_I_0   (:,:,:)
      H_II_0_P   => H_II_0  (:,:,:)

      J_I_0_P    => J_I_0   (:,:,:)
      J_II_0_P   => J_II_0  (:,:,:)

      S_Sigma_P  => S_Sigma (:,:,:)

    END IF

    ! --- NES Emissivities and Opacities ---

    CALL ComputeNeutrinoOpacityRates_NES &
           ( 1, nE_G, 1, nSpecies, 1, nX, W2_N, J_P, J0_P, H_I_0_P, H_II_0_P, &
             Eta_NES_P, Chi_NES_P )

    ! --- Pair Emissivities and Opacities ---

    CALL ComputeNeutrinoOpacityRates_Pair &
           ( 1, nE_G, 1, nSpecies, 1, nX, W2_N, J_P, J0_P, J_I_0_P, J_II_0_P, &
             Eta_Pair_P, Chi_Pair_P )

    ! --- Brem Emissivities and Opacities ---

    CALL ComputeNeutrinoOpacityRates_Brem &
           ( 1, nE_G, 1, nSpecies, 1, nX, W2_N, J_P, J0_P, S_Sigma_P, &
             Eta_Brem_P, Chi_Brem_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NES_P, Eta_NES_P, Chi_NES, Eta_NES )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_Pair_P, Eta_Pair_P, Chi_Pair, Eta_Pair )

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_Brem_P, Eta_Brem_P, Chi_Brem, Eta_Brem )

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

      D_P => D_T(1:nX)
      T_P => T_T(1:nX)
      Y_P => Y_T(1:nX)
      E_P => E_T(1:nX)

      CALL ArrayPack &
             ( nX, UnpackIndex, D, Y, T, E, D_P, Y_P, T_P, E_P )

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

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, T_P, T )

    END IF

  END SUBROUTINE UpdateTemperature_Packed


  SUBROUTINE InitializeRHS_FP &
    ( J, H_u_1, H_u_2, H_u_3, D, Y, E, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: D, Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotV, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3
    REAL(DP) :: SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( vDotV, SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( vDotV, SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( vDotV, SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3, &
    !$OMP          k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iN_X = 1, nX_G

      V_d_1(iN_X) = Gm_dd_11(iN_X) * V_u_1(iN_X)
      V_d_2(iN_X) = Gm_dd_22(iN_X) * V_u_2(iN_X)
      V_d_3(iN_X) = Gm_dd_33(iN_X) * V_u_3(iN_X)

      ! --- Specific Fluid Energy ---

      vDotV =   V_u_1(iN_X) * V_d_1(iN_X) &
              + V_u_2(iN_X) * V_d_2(iN_X) &
              + V_u_3(iN_X) * V_d_3(iN_X)

      Ef(iN_X) = E(iN_X) + Half * vDotV

      ! --- Compute Omega (Richardson damping coeff.) based on velocity ---
      ! --- For first interation: must be consistent with initial guess ---

      Omega(iN_X) = One / ( One + SQRT( vDotV ) )

      ! --- Store Initial Matter State ---

      Y_old (iN_X) = Y (iN_X)
      Ef_old(iN_X) = Ef(iN_X)

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

      SUM_Y  = Zero
      SUM_Ef = Zero
      SUM_V1 = Zero
      SUM_V2 = Zero
      SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
      !$OMP          N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 ) &
      !$OMP REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
      !$ACC LOOP VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
      !$ACC          N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 ) &
      !$ACC REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#endif
      DO iS   = 1, nSpecies
      DO iN_E = 1, nE_G

        H_d_1(iN_E,iS,iN_X) = Gm_dd_11(iN_X) * H_u_1(iN_E,iS,iN_X)
        H_d_2(iN_E,iS,iN_X) = Gm_dd_22(iN_X) * H_u_2(iN_E,iS,iN_X)
        H_d_3(iN_E,iS,iN_X) = Gm_dd_33(iN_X) * H_u_3(iN_E,iS,iN_X)

        vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iS,iN_X) &
                + V_u_2(iN_X) * H_d_2(iN_E,iS,iN_X) &
                + V_u_3(iN_X) * H_d_3(iN_E,iS,iN_X)

        k_dd = EddingtonTensorComponents_dd &
                 ( J    (iN_E,iS,iN_X), H_u_1(iN_E,iS,iN_X), &
                   H_u_2(iN_E,iS,iN_X), H_u_3(iN_E,iS,iN_X), &
                   Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vDotK_d_1 &
          = ( V_u_1(iN_X) * k_dd(1,1) &
            + V_u_2(iN_X) * k_dd(2,1) &
            + V_u_3(iN_X) * k_dd(3,1) ) * J(iN_E,iS,iN_X)
        vDotK_d_2 &
          = ( V_u_1(iN_X) * k_dd(1,2) &
            + V_u_2(iN_X) * k_dd(2,2) &
            + V_u_3(iN_X) * k_dd(3,2) ) * J(iN_E,iS,iN_X)
        vDotK_d_3 &
          = ( V_u_1(iN_X) * k_dd(1,3) &
            + V_u_2(iN_X) * k_dd(2,3) &
            + V_u_3(iN_X) * k_dd(3,3) ) * J(iN_E,iS,iN_X)

        ! --- Eulerian Neutrino Number Density ---

        N_nu = J(iN_E,iS,iN_X) + vDotH

        ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

        E_nu = J(iN_E,iS,iN_X) + Two * vDotH

        ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

        F_nu_d_1 &
          = H_d_1(iN_E,iS,iN_X) + V_d_1(iN_X) * J(iN_E,iS,iN_X) + vDotK_d_1

        F_nu_d_2 &
          = H_d_2(iN_E,iS,iN_X) + V_d_2(iN_X) * J(iN_E,iS,iN_X) + vDotK_d_2

        F_nu_d_3 &
          = H_d_3(iN_E,iS,iN_X) + V_d_3(iN_X) * J(iN_E,iS,iN_X) + vDotK_d_3

        ! --- Old States for Neutrino Number Density and Flux ---

        C_J    (iN_E,iS,iN_X) = J    (iN_E,iS,iN_X) + vDotH
        C_H_d_1(iN_E,iS,iN_X) = H_d_1(iN_E,iS,iN_X) + vDotK_d_1
        C_H_d_2(iN_E,iS,iN_X) = H_d_2(iN_E,iS,iN_X) + vDotK_d_2
        C_H_d_3(iN_E,iS,iN_X) = H_d_3(iN_E,iS,iN_X) + vDotK_d_3

        IF ( iS <= iNuE_Bar ) THEN
        SUM_Y  = SUM_Y  + N_nu     * W2_S(iN_E) * LeptonNumber(iS)
        END IF

        SUM_Ef = SUM_Ef + E_nu     * W3_S(iN_E)
        SUM_V1 = SUM_V1 + F_nu_d_1 * W3_S(iN_E)
        SUM_V2 = SUM_V2 + F_nu_d_2 * W3_S(iN_E)
        SUM_V3 = SUM_V3 + F_nu_d_3 * W3_S(iN_E)

      END DO
      END DO

      ! --- Include Old Matter State in Constant (C) Terms ---

      C_Y    (iN_X) = U_Y    (iN_X) + wMatterRHS(iY ) * SUM_Y  * S_Y    (iN_X)
      C_Ef   (iN_X) = U_Ef   (iN_X) + wMatterRHS(iEf) * SUM_Ef * S_Ef   (iN_X)
      C_V_d_1(iN_X) = U_V_d_1(iN_X) + wMatterRHS(iV1) * SUM_V1 * S_V_d_1(iN_X)
      C_V_d_2(iN_X) = U_V_d_2(iN_X) + wMatterRHS(iV2) * SUM_V2 * S_V_d_2(iN_X)
      C_V_d_3(iN_X) = U_V_d_3(iN_X) + wMatterRHS(iV3) * SUM_V3 * S_V_d_3(iN_X)

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

        Jnorm(iS,iN_X) = WNORM( J(:,iS,iN_X), W2_S )

      END IF

    END DO
    END DO

  END SUBROUTINE ComputeJNorm


  SUBROUTINE ComputeMatterRHS_FP &
    ( MASK, Fm, Gm, &
      J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS
    REAL(DP) :: k_dd(3,3), vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3
    REAL(DP) :: N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3
    REAL(DP) :: SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3, &
    !$OMP          N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3, &
    !$OMP          k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        SUM_Y  = Zero
        SUM_Ef = Zero
        SUM_V1 = Zero
        SUM_V2 = Zero
        SUM_V3 = Zero

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
        !$OMP          N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 ) &
        !$OMP REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR COLLAPSE(2) &
        !$ACC PRIVATE( k_dd, vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3, &
        !$ACC          N_nu, E_nu, F_nu_d_1, F_nu_d_2, F_nu_d_3 ) &
        !$ACC REDUCTION( + : SUM_Y, SUM_Ef, SUM_V1, SUM_V2, SUM_V3 )
#endif
        DO iS   = 1, nSpecies
        DO iN_E = 1, nE_G

          vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iS,iN_X) &
                  + V_u_2(iN_X) * H_d_2(iN_E,iS,iN_X) &
                  + V_u_3(iN_X) * H_d_3(iN_E,iS,iN_X)

          k_dd = EddingtonTensorComponents_dd &
                   ( J    (iN_E,iS,iN_X), H_u_1(iN_E,iS,iN_X), &
                     H_u_2(iN_E,iS,iN_X), H_u_3(iN_E,iS,iN_X), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

          vDotK_d_1 &
            = ( V_u_1(iN_X) * k_dd(1,1) &
              + V_u_2(iN_X) * k_dd(2,1) &
              + V_u_3(iN_X) * k_dd(3,1) ) * J(iN_E,iS,iN_X)
          vDotK_d_2 &
            = ( V_u_1(iN_X) * k_dd(1,2) &
              + V_u_2(iN_X) * k_dd(2,2) &
              + V_u_3(iN_X) * k_dd(3,2) ) * J(iN_E,iS,iN_X)
          vDotK_d_3 &
            = ( V_u_1(iN_X) * k_dd(1,3) &
              + V_u_2(iN_X) * k_dd(2,3) &
              + V_u_3(iN_X) * k_dd(3,3) ) * J(iN_E,iS,iN_X)

          ! --- Eulerian Neutrino Number Density ---

          N_nu = J(iN_E,iS,iN_X) + vDotH

          ! --- Eulerian Neutrino Energy Density (Scaled by Neutrino Energy) ---

          E_nu = J(iN_E,iS,iN_X) + Two * vDotH

          ! --- Eulerian Neutrino Momentum Density (Scaled by Neutrino Energy) ---

          F_nu_d_1 &
            = H_d_1(iN_E,iS,iN_X) + V_d_1(iN_X) * J(iN_E,iS,iN_X) + vDotK_d_1

          F_nu_d_2 &
            = H_d_2(iN_E,iS,iN_X) + V_d_2(iN_X) * J(iN_E,iS,iN_X) + vDotK_d_2

          F_nu_d_3 &
            = H_d_3(iN_E,iS,iN_X) + V_d_3(iN_X) * J(iN_E,iS,iN_X) + vDotK_d_3

          IF ( iS <= iNuE_Bar ) THEN
          SUM_Y  = SUM_Y  + N_nu     * W2_S(iN_E) * LeptonNumber(iS)
          END IF

          SUM_Ef = SUM_Ef + E_nu     * W3_S(iN_E)
          SUM_V1 = SUM_V1 + F_nu_d_1 * W3_S(iN_E)
          SUM_V2 = SUM_V2 + F_nu_d_2 * W3_S(iN_E)
          SUM_V3 = SUM_V3 + F_nu_d_3 * W3_S(iN_E)

        END DO
        END DO

        G_Y    (iN_X) = C_Y    (iN_X) - wMatterRHS(iY ) * SUM_Y  * S_Y    (iN_X)
        G_Ef   (iN_X) = C_Ef   (iN_X) - wMatterRHS(iEf) * SUM_Ef * S_Ef   (iN_X)
        G_V_d_1(iN_X) = C_V_d_1(iN_X) - wMatterRHS(iV1) * SUM_V1 * S_V_d_1(iN_X)
        G_V_d_2(iN_X) = C_V_d_2(iN_X) - wMatterRHS(iV2) * SUM_V2 * S_V_d_2(iN_X)
        G_V_d_3(iN_X) = C_V_d_3(iN_X) - wMatterRHS(iV3) * SUM_V3 * S_V_d_3(iN_X)

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
    ( MASK, Fm, Gm, dt, &
      J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP),                   INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS
    REAL(DP) :: k_dd(3,3), vDotH, vDotK_d_1, vDotK_d_2, vDotK_d_3
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
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        vDotH =   V_u_1(iN_X) * H_d_1(iN_E,iS,iN_X) &
                + V_u_2(iN_X) * H_d_2(iN_E,iS,iN_X) &
                + V_u_3(iN_X) * H_d_3(iN_E,iS,iN_X)

        k_dd = EddingtonTensorComponents_dd &
                 ( J    (iN_E,iS,iN_X), H_u_1(iN_E,iS,iN_X), &
                   H_u_2(iN_E,iS,iN_X), H_u_3(iN_E,iS,iN_X), &
                   Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vDotK_d_1 &
          = ( V_u_1(iN_X) * k_dd(1,1) &
            + V_u_2(iN_X) * k_dd(2,1) &
            + V_u_3(iN_X) * k_dd(3,1) ) * J(iN_E,iS,iN_X)
        vDotK_d_2 &                                                  
          = ( V_u_1(iN_X) * k_dd(1,2) &
            + V_u_2(iN_X) * k_dd(2,2) &
            + V_u_3(iN_X) * k_dd(3,2) ) * J(iN_E,iS,iN_X)
        vDotK_d_3 &                                                  
          = ( V_u_1(iN_X) * k_dd(1,3) &
            + V_u_2(iN_X) * k_dd(2,3) &
            + V_u_3(iN_X) * k_dd(3,3) ) * J(iN_E,iS,iN_X)

        ! --- Emissivity ---

        Eta_T =   Chi_EmAb(iN_E,iS,iN_X) * J0(iN_E,iS,iN_X) &
                + Eta_NES (iN_E,iS,iN_X) &
                + Eta_Pair(iN_E,iS,iN_X) &
                + Eta_Brem(iN_E,iS,iN_X)

        ! --- Number Opacity ---

        Chi_T =   Chi_EmAb(iN_E,iS,iN_X) &
                + Chi_NES (iN_E,iS,iN_X) &
                + Chi_Pair(iN_E,iS,iN_X) &
                + Chi_Brem(iN_E,iS,iN_X)

        ! --- Number Flux Opacity ---

        Kappa = Chi_T + Sigma_Iso(iN_E,iN_X)

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        ! --- Number Equation ---

        Gm(iOS+iCR_N,iN_X) &
          = ( One - Omega(iN_X) ) *     J(iN_E,iS,iN_X) &
            +       Omega(iN_X)   * ( C_J(iN_E,iS,iN_X) - vDotH + dt * Eta_T ) &
                                    / ( One + dt * Chi_T )

        Fm(iOS+iCR_N,iN_X) &
          = Gm(iOS+iCR_N,iN_X) - J(iN_E,iS,iN_X)

        ! --- Number Flux 1 Equation ---

        Gm(iOS+iCR_G1,iN_X) &
          = ( One - Omega(iN_X) ) *     H_d_1(iN_E,iS,iN_X) &
            +       Omega(iN_X)   * ( C_H_d_1(iN_E,iS,iN_X) - vDotK_d_1 ) &
                                    / ( One + dt * Kappa )

        Fm(iOS+iCR_G1,iN_X) &
          = Gm(iOS+iCR_G1,iN_X) - H_d_1(iN_E,iS,iN_X)

        ! --- Number Flux 2 Equation ---

        Gm(iOS+iCR_G2,iN_X) &
          = ( One - Omega(iN_X) ) *     H_d_2(iN_E,iS,iN_X) &
            +       Omega(iN_X)   * ( C_H_d_2(iN_E,iS,iN_X) - vDotK_d_2 ) &
                                    / ( One + dt * Kappa )

        Fm(iOS+iCR_G2,iN_X) &
          = Gm(iOS+iCR_G2,iN_X) - H_d_2(iN_E,iS,iN_X)

        ! --- Number Flux 3 Equation ---

        Gm(iOS+iCR_G3,iN_X) &
          = ( One - Omega(iN_X) ) *     H_d_3(iN_E,iS,iN_X) &
            +       Omega(iN_X)   * ( C_H_d_3(iN_E,iS,iN_X) - vDotK_d_3 ) &
                                    / ( One + dt * Kappa )

        Fm(iOS+iCR_G3,iN_X) &
          = Gm(iOS+iCR_G3,iN_X) - H_d_3(iN_E,iS,iN_X)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_FP


  SUBROUTINE UpdateMatterRHS_FP &
    ( MASK, Fm, Gm, Y, E, V_u_1, V_u_2, V_u_3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:)    , INTENT(inout) :: Y, E, V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_X
    REAL(DP) :: vDotV

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

        Fm(iY ,iN_X) = Gm(iY ,iN_X) - U_Y    (iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef   (iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)

        U_Y    (iN_X) = Gm(iY ,iN_X)
        U_Ef   (iN_X) = Gm(iEf,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        Y    (iN_X) = U_Y    (iN_X) * Y_old (iN_X)
        Ef   (iN_X) = U_Ef   (iN_X) * Ef_old(iN_X)
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

        E(iN_X) = Ef(iN_X) - Half * vDotV

        ! --- Compute Omega (Richardson damping coeff.) based on velocity ---

        Omega(iN_X) = One / ( One + SQRT( vDotV ) )

      END IF

    END DO

  END SUBROUTINE UpdateMatterRHS_FP


  SUBROUTINE UpdateNeutrinoRHS_FP &
    ( MASK, Fm, Gm, J, H_u_1, H_u_2, H_u_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(:)    , INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:)  , INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(:)    , INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    INTEGER  :: iN_E, iN_X, iS, iOS

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
    DO iN_X = 1, nX_G
    DO iS   = 1, nSpecies
    DO iN_E = 1, nE_G

      IF( MASK(iN_X) )THEN

        iOS = ( (iN_E-1) + (iS-1) * nE_G ) * nCR

        Fm(iOS+iCR_N ,iN_X) = Gm(iOS+iCR_N ,iN_X) - J    (iN_E,iS,iN_X)
        Fm(iOS+iCR_G1,iN_X) = Gm(iOS+iCR_G1,iN_X) - H_d_1(iN_E,iS,iN_X)
        Fm(iOS+iCR_G2,iN_X) = Gm(iOS+iCR_G2,iN_X) - H_d_2(iN_E,iS,iN_X)
        Fm(iOS+iCR_G3,iN_X) = Gm(iOS+iCR_G3,iN_X) - H_d_3(iN_E,iS,iN_X)

        J    (iN_E,iS,iN_X) = Gm(iOS+iCR_N ,iN_X)
        H_d_1(iN_E,iS,iN_X) = Gm(iOS+iCR_G1,iN_X)
        H_d_2(iN_E,iS,iN_X) = Gm(iOS+iCR_G2,iN_X)
        H_d_3(iN_E,iS,iN_X) = Gm(iOS+iCR_G3,iN_X)

        H_u_1(iN_E,iS,iN_X) = H_d_1(iN_E,iS,iN_X) / Gm_dd_11(iN_X)
        H_u_2(iN_E,iS,iN_X) = H_d_2(iN_E,iS,iN_X) / Gm_dd_22(iN_X)
        H_u_3(iN_E,iS,iN_X) = H_d_3(iN_E,iS,iN_X) / Gm_dd_33(iN_X)

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

            IF( ABS( AA11 ) < SqrtTiny )THEN

              B(1,iN_X) = Zero

            ELSE

              B(1,iN_X) = AB1 / AA11

            END IF

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

            IF( ABS( DET_AA ) < SqrtTiny )THEN

              B(1,iN_X) = Zero
              B(2,iN_X) = Zero

            ELSE

              B(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
              B(2,iN_X) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

            END IF

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

        DO iS = 1, nSpecies

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
    REAL(DP) :: Fnorm_Y, Fnorm_Ef, Fnorm_V

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
