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
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Im_EmAb_FP, &
    Timer_Im_NestInner, &
    Timer_Im_ComputeOpacity, &
    Timer_Im_ComputeRate, &
    Timer_Im_ComputeLS, &
    Timer_Im_UpdateFP
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
    nCR, &
    iNuE, &
    iNuE_Bar
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_Points, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeEquilibriumDistributionAndDerivatives_Points, &
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
  ! PUBLIC :: SolveMatterEquations_EmAb_FP
  PUBLIC :: SolveMatterEquations_FP_NestedAA

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

  REAL(DP), PARAMETER :: WFactor_FP = FourPi / PlanckConstant**3

  INTEGER  :: nE_G, nX_G, nZ(4)
  INTEGER  :: n_FP, n_FP_inner, n_FP_outer
  INTEGER  :: iE_B, iE_E
  INTEGER  :: iX_B(3), iX_E(3)

  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:)
  REAL(DP), ALLOCATABLE :: W3_S(:)

  REAL(DP), ALLOCATABLE :: AMAT(:,:,:)
  REAL(DP), ALLOCATABLE :: BVEC(:,:)
  REAL(DP), ALLOCATABLE :: TAU (:,:)
  REAL(DP), ALLOCATABLE :: WORK(:,:)
  INTEGER,  ALLOCATABLE :: IPIV(:,:)
  INTEGER,  ALLOCATABLE :: INFO(:)
  INTEGER               :: LWORK


  ! --- Solver Parameters to be initialized
  LOGICAL  :: UsePreconditionerEmAb
  INTEGER  :: M_FP, M_outer, M_inner
  INTEGER  :: MaxIter_outer, MaxIter_inner
  REAL(DP) :: Rtol, Utol


  INTEGER :: iS_1 = iNuE
  INTEGER :: iS_2 = iNuE_Bar
  INTEGER :: OS_JNuE, OS_H1NuE, OS_H2NuE, OS_H3NuE
  INTEGER :: OS_JNuE_Bar, OS_H1NuE_Bar, OS_H2NuE_Bar, OS_H3NuE_Bar
  INTEGER :: iY  = 1
  INTEGER :: iEf = 2
  INTEGER :: iV1 = 3
  INTEGER :: iV2 = 4
  INTEGER :: iV3 = 5

  ! --- Temporary arrays for scatter/gather (packing)

  REAL(DP), ALLOCATABLE, TARGET :: P1D(:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D(:,:,:,:)

  INTEGER, PARAMETER :: nP1D = 4
  INTEGER, PARAMETER :: iP1D_D = 1
  INTEGER, PARAMETER :: iP1D_T = 2
  INTEGER, PARAMETER :: iP1D_Y = 3
  INTEGER, PARAMETER :: iP1D_E = 4

  INTEGER, PARAMETER :: nP2D = 16
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
  INTEGER, PARAMETER :: iP2D_Sig_1        = 15
  INTEGER, PARAMETER :: iP2D_Sig_2        = 16

  INTEGER, PARAMETER :: nP3D = 10
  INTEGER, PARAMETER :: iP3D_Phi_0_In_NES_1  = 1
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_NES_1  = 2
  INTEGER, PARAMETER :: iP3D_Phi_0_In_NES_2  = 3
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_NES_2  = 4
  INTEGER, PARAMETER :: iP3D_Phi_0_In_Pair_1 = 5
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_Pair_1 = 6
  INTEGER, PARAMETER :: iP3D_Phi_0_In_Pair_2 = 7
  INTEGER, PARAMETER :: iP3D_Phi_0_Ot_Pair_2 = 8
  INTEGER, PARAMETER :: iP3D_WORK1 = 9
  INTEGER, PARAMETER :: iP3D_WORK2 = 10

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

    n_FP       = 5 + 2*(nE_G + 3*nE_G)
    n_FP_inner = 2*(nE_G + 3*nE_G)
    n_FP_outer = 5

    OS_JNuE = 0
    OS_H1NuE = OS_JNuE + nE_G
    OS_H2NuE = OS_H1NuE + nE_G
    OS_H3NuE = OS_H2NuE + nE_G
    OS_JNuE_Bar = OS_H3NuE + nE_G
    OS_H1NuE_Bar = OS_JNuE_Bar + nE_G
    OS_H2NuE_Bar = OS_H1NuE_Bar + nE_G
    OS_H3NuE_Bar = OS_H2NuE_Bar + nE_G

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
    ALLOCATE( IPIV(n_FP,nX_G) )
    ALLOCATE( INFO(nX_G) )

    ALLOCATE( P1D(nX_G,nP1D) )
    ALLOCATE( P2D(nE_G,nX_G,nP2D) )
    ALLOCATE( P3D(nE_G,nE_G,nX_G,nP3D) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$OMP MAP( alloc: AMAT, BVEC, TAU, IPIV, INFO, P1D, P2D, P3D )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$ACC CREATE( AMAT, BVEC, TAU, IPIV, INFO, P1D, P2D, P3D )
#endif

    IF ( M_FP> 3 ) THEN

      CALL LinearLeastSquares_LWORK &
             ( 'N', n_FP, M_FP-1, 1, &
               AMAT, n_FP, BVEC, n_FP, TMP, LWORK )
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


  SUBROUTINE InitializeNeutrinoMatterSolverParameters &
    ( Mout, Min, MaxItout, MaxItin, RelTol, AbsTol, EmAbPrecond)

    INTEGER,  OPTIONAL, INTENT(in) :: Mout, Min
    INTEGER,  OPTIONAL, INTENT(in) :: MaxItout, MaxItin
    REAL(DP), OPTIONAL, INTENT(in) :: RelTol, AbsTol
    LOGICAL,  OPTIONAL, INTENT(in) :: EmAbPrecond

    IF ( PRESENT( Mout ) ) THEN
      M_outer = Mout
    ELSE
      M_outer = 3
    END IF
    IF ( PRESENT( Min ) ) THEN
      M_inner = Min
    ELSE
      M_inner = 3
    END IF
    IF ( PRESENT( MaxItout ) ) THEN
      MaxIter_outer = MaxItout
    ELSE
      MaxIter_outer = 100
    END IF
    IF ( PRESENT( MaxItin ) ) THEN
      MaxIter_inner = MaxItin
    ELSE
      MaxIter_inner = 100
    END IF
    IF ( PRESENT( RelTol ) ) THEN
      Rtol = RelTol
    ELSE
      Rtol = 1.0d-08
    END IF
    IF ( PRESENT( AbsTol ) ) THEN
      Utol = AbsTol
    ELSE
      Utol = 1.0d-10
    END IF
    IF ( PRESENT( EmAbPrecond ) ) THEN
      UsePreconditionerEmAb = EmAbPrecond
    ELSE
      UsePreconditionerEmAb = .FALSE.
    END IF

    M_FP = MAX(M_outer, M_inner)

  END SUBROUTINE InitializeNeutrinoMatterSolverParameters


  SUBROUTINE FinalizeNeutrinoMatterSolver

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, &
    !$OMP               AMAT, BVEC, TAU, WORK, IPIV, INFO, P1D, P2D, P3D )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, &
    !$ACC         AMAT, BVEC, TAU, WORK, IPIV, INFO, P1D, P2D, P3D )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S )
    DEALLOCATE( AMAT, BVEC, TAU, WORK, IPIV, INFO )
    DEALLOCATE( P1D, P2D, P3D )


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

      W2(iN_E) = WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**2
      W3(iN_E) = WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**3

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveMatterEquations_FP_NestedAA &
    ( dt, J, H_u_1, H_u_2, H_u_3, V_u_1, V_u_2, V_u_3, &
      D, T, Y, E, Gm_dd_11, Gm_dd_22, Gm_dd_33, &
      nIterations_Inner, nIterations_Outer )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                    INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nX_G), INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), &
                                 INTENT(inout) :: J, H_u_1, H_u_2, H_u_3
    REAL(DP), DIMENSION(1:nX_G), INTENT(inout) :: V_u_1, V_u_2, V_u_3
    REAL(DP), DIMENSION(1:nX_G), INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(1:nX_G), INTENT(out)   :: nIterations_Inner, nIterations_Outer


    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: Jold, Jnew
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: H_d_1, H_d_2, H_d_3
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: H1old, H2old, H3old
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: H1new, H2new, H3new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: CJ, CH1, CH2, CH3
    REAL(DP), DIMENSION(       1:nX_G,1:nSpecies) :: Jnorm, Hnorm

    REAL(DP), DIMENSION(1:nX_G) :: Omega, NormVsq
    REAL(DP), DIMENSION(1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(1:nX_G) :: Ef, Efold, S_Ef, C_Ef, Unew_Ef, GVEC_Ef
    REAL(DP), DIMENSION(1:nX_G) :: V_d_1, V1old, S_V1, C_V1, Unew_V1, GVEC_V1
    REAL(DP), DIMENSION(1:nX_G) :: V_d_2, V2old, S_V2, C_V2, Unew_V2, GVEC_V2
    REAL(DP), DIMENSION(1:nX_G) :: V_d_3, V3old, S_V3, C_V3, Unew_V3, GVEC_V3

    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: Phi_0_In_NES
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: Phi_0_Ot_NES
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: Phi_0_In_Pair
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies) :: Phi_0_Ot_Pair

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: J0, Chi, Sig
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: Chi_NES, Eta_NES
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: Chi_Pair, Eta_Pair

    REAL(DP), DIMENSION(1:n_FP_outer,1:M_outer,1:nX_G) :: AMAT_outer, GVEC_outer, FVEC_outer
    REAL(DP), DIMENSION(1:n_FP_outer,          1:nX_G) :: GVECm_outer, FVECm_outer
    REAL(DP), DIMENSION(1:n_FP_outer,          1:nX_G) :: BVEC_outer
    REAL(DP), DIMENSION(             1:M_outer,1:nX_G) :: Alpha_outer

    REAL(DP), DIMENSION(1:n_FP_inner,1:M_inner,1:nX_G) :: AMAT_inner, GVEC_inner, FVEC_inner
    REAL(DP), DIMENSION(1:n_FP_inner,          1:nX_G) :: GVECm_inner, FVECm_inner
    REAL(DP), DIMENSION(1:n_FP_inner,          1:nX_G) :: BVEC_inner
    REAL(DP), DIMENSION(             1:M_inner,1:nX_G) :: Alpha_inner

    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE_OUTER, ITERATE_INNER
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex_outer, UnpackIndex_outer
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex_inner, UnpackIndex_inner

    INTEGER :: k_outer, Mk_outer, nX_P_outer
    INTEGER :: k_inner, Mk_inner, nX_P_inner
    INTEGER :: iN_X


    ITERATE_OUTER(:) = .TRUE.
    ITERATE_INNER(:) = .TRUE.

! --- These pragmas have NOT been updated ---
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE_OUTER, ITERATE_INNER ) &
    !$OMP MAP( alloc: Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP             Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP             Jold_1, Jnew_1, Jnorm_1, &
    !$OMP             Jold_2, Jnew_2, Jnorm_2, &
    !$OMP             Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP             Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP             PackIndex_outer, UnpackIndex_outer, &
    !$OMP             PackIndex_inner, UnpackIndex_inner, &
    !$OMP             AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$OMP             GVECm_outer, FVECm_outer, Alpha_outer, &
    !$OMP             AMAT_inner, BVEC_inner, GVEC_inner, FVEC_inner, &
    !$OMP             GVECm_inner, FVECm_inner, Alpha_inner )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE_OUTER, ITERATE_INNER ) &
    !$ACC CREATE( Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$ACC         GVECm_outer, FVECm_outer, Alpha_outer, &
    !$ACC         AMAT_inner, BVEC_inner, GVEC_inner, FVEC_inner, &
    !$ACC         GVECm_inner, FVECm_inner, Alpha_inner )
#endif


    DO iN_X = 1, nX_G
      V_d_1(iN_X) = Gm_dd_11(iN_X) * V_u_1(iN_X)
      V_d_2(iN_X) = Gm_dd_22(iN_X) * V_u_2(iN_X)
      V_d_3(iN_X) = Gm_dd_33(iN_X) * V_u_3(iN_X)
      H_d_1(:,iN_X,:) = Gm_dd_11(iN_X) * H_u_1(:,iN_X,:)
      H_d_2(:,iN_X,:) = Gm_dd_22(iN_X) * H_u_2(:,iN_X,:)
      H_d_3(:,iN_X,:) = Gm_dd_33(iN_X) * H_u_3(:,iN_X,:)
    END DO

    DO iN_X = 1, nX_G
    ! Compute fluid energy density ***over the matter density***

      ! Ef(iN_X) = D(iN_X) * (E(iN_X) + Half * (V_u_1(iN_X)*V_d_1(iN_X) &
      !                      + V_u_2(iN_X)*V_d_2(iN_X) + V_u_3(iN_X)*V_d_3(iN_X)))

      NormVsq(iN_X) =  V_u_1(iN_X) * V_d_1(iN_X) &
                     + V_u_2(iN_X) * V_d_2(iN_X) &
                     + V_u_3(iN_X) * V_d_3(iN_X)

      Ef(iN_X) = E(iN_X) + Half * NormVsq(iN_X)

    END DO

    CALL ArrayCopy( Y, Ef, Yold, Efold )
    CALL ArrayCopy( V_d_1, V_d_2, V_d_3, V1old, V2old, V3old )

    CALL ArrayCopy( J, Jold )
    CALL ArrayCopy( H_d_1, H_d_2, H_d_3, H1old, H2old, H3old )

    ! IF( UsePreconditionerEmAb )THEN
    !
    !   CALL TimersStart( Timer_Im_EmAb_FP )
    !
    !   CALL SolveMatterEquations_EmAb_FP &
    !          ( dt, iS_1, iS_2, J_1, J_2, &
    !            Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )
    !
    !   CALL TimersStop( Timer_Im_EmAb_FP )
    !
    ! END IF

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( D, T, Y, J0, Chi, Sig, Phi_0_In_NES, Phi_0_Ot_NES, &
             Phi_0_In_Pair, Phi_0_Ot_Pair )

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Initial RHS ---

    CALL InitializeRHS_FP &
           ( J, Jold, CJ, Jnew, H_d_1, H1old, CH1, H1new, &
             H_d_2, H2old, CH2, H2new, H_d_3, H3old, CH3, H3new, &
             D, Y, Yold, C_Y, S_Y, Unew_Y, Ef, Efold, C_Ef, S_Ef, Unew_Ef, &
             V_d_1, V1old, C_V1, S_V1, Unew_V1, V_d_2, V2old, C_V2, S_V2, &
             Unew_V2, V_d_3, V3old, C_V3, S_V3, Unew_V3, &
             Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    k_outer = 0
    DO WHILE( ANY( ITERATE_OUTER(:) ) .AND. k_outer < MaxIter_outer )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      CALL ComputeJNorm &
             ( ITERATE_OUTER, Jnew, Jnorm )

      CALL CreatePackIndex &
             ( ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      ! PRINT*, "kouter = ", k_outer
      ! PRINT*, "D = ", D(:) / (Gram / Centimeter**3)
      ! PRINT*, "T = ", T(:) / (MeV)
      ! PRINT*, "Y = ", Y(:)
      ! PRINT*, "Ef = ", Ef(:)
      ! PRINT*, "E = ", E(:)
      ! PRINT*, "V_d_1 = ", V_d_1(:)
      ! PRINT*, "V_d_2 = ", V_d_2(:)
      ! PRINT*, "V_d_3 = ", V_d_3(:)


      IF ( k_outer > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Im_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( D, T, Y, J0, Chi, Sig, Phi_0_In_NES, Phi_0_Ot_NES, &
                 Phi_0_In_Pair, Phi_0_Ot_Pair, &
                 ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        CALL TimersStop( Timer_Im_ComputeOpacity )

      END IF

      ! START INNER LOOP
      CALL TimersStart( Timer_Im_NestInner )

      k_inner = 0
      DO WHILE( ANY( ITERATE_INNER(:) ) .AND. k_inner < MaxIter_inner )

        k_inner = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

        CALL TimersStart( Timer_Im_ComputeRate )

        CALL ComputeRates_Packed &
               ( Jnew, &
                 Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, Phi_0_Ot_Pair, &
                 Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
                 ITERATE_INNER, nX_P_inner, PackIndex_inner, &
                 UnpackIndex_inner, nX_P_outer )

        CALL TimersStop( Timer_Im_ComputeRate )


        ! --- Compute Omega (Richardson damping coefficient) based on velocity
        DO iN_X = 1, nX_G
          Omega(iN_X) = One / (One + SQRT(NormVsq(iN_X)))
        END DO
        ! --- Right-Hand Side Vectors and Residuals (inner) ---

        CALL ComputeNeutrinoRHS_FP &
               ( ITERATE_INNER, n_FP_inner, FVECm_inner, GVECm_inner, &
                 dt, Omega, V_u_1, V_u_2, V_u_3, &
                 CJ, CH1, CH2, CH3, &
                 Jnew, H1new, H2new, H3new, &
                 J0, Chi, Sig, &
                 Chi_NES, Eta_NES, &
                 Chi_Pair, Eta_Pair, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        ! --- Anderson Acceleration (inner) ---

        CALL TimersStart( Timer_Im_ComputeLS )

        CALL SolveLS_FP &
               ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                 FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner, &
                 AMAT_inner, BVEC_inner, Alpha_inner )

        CALL TimersStop( Timer_Im_ComputeLS )

        ! --- Update Residuals and Solution Vectors (inner) ---

        CALL TimersStart( Timer_Im_UpdateFP )

        CALL UpdateNeutrinoRHS_FP &
               ( ITERATE_INNER, n_FP_inner, &
                 FVECm_inner, GVECm_inner, Jnew, H1new, H2new, H3new )

        ! --- Check Convergence (inner) ---

        CALL CheckConvergenceInner &
               ( ITERATE_INNER, n_FP_inner, k_inner, Rtol, &
                 nIterations_Inner, FVECm_inner, Jnorm )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                 FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Im_UpdateFP )

      END DO

      CALL TimersStop( Timer_Im_NestInner )

      ! PRINT*, "after inner = "
      ! PRINT*, "D = ", D(:) / (Gram / Centimeter**3)
      ! PRINT*, "T = ", T(:) / (MeV)
      ! PRINT*, "Y = ", Y(:)
      ! PRINT*, "Ef = ", Ef(:)
      ! PRINT*, "E = ", E(:)
      ! PRINT*, "V_d_1 = ", V_d_1(:)
      ! PRINT*, "V_d_2 = ", V_d_2(:)
      ! PRINT*, "V_d_3 = ", V_d_3(:)

      ! --- Right-Hand Side Vectors and Residuals (outer) ---

      CALL ComputeMatterRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, FVECm_outer, GVECm_outer, &
               Jnew, H1new, H2new, H3new, &
               C_Y, S_Y, Unew_Y, GVEC_Y, &
               C_Ef, S_Ef, Unew_Ef, GVEC_Ef, &
               V_d_1, C_V1, S_V1, Unew_V1, GVEC_V1, &
               V_d_2, C_V2, S_V2, Unew_V2, GVEC_V2, &
               V_d_3, C_V3, S_V3, Unew_V3, GVEC_V3, &
               Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      ! --- Anderson Acceleration (outer) ---

      CALL TimersStart( Timer_Im_ComputeLS )

      CALL SolveLS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
               FVECm_outer, GVECm_outer, FVEC_outer, GVEC_outer, &
               AMAT_outer, BVEC_outer, Alpha_outer )

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors (outer) ---

      CALL TimersStart( Timer_Im_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, &
             Y, Yold, Unew_Y, &
             Ef, Efold, Unew_Ef, &
             V_d_1, V1old, Unew_V1, &
             V_d_2, V2old, Unew_V2, &
             V_d_3, V3old, Unew_V3, &
             FVECm_outer, GVECm_outer )

      ! --- Update V upper ---

      DO iN_X = 1, nX_G
        V_u_1(iN_X) = V_d_1(iN_X) / Gm_dd_11(iN_X)
        V_u_2(iN_X) = V_d_2(iN_X) / Gm_dd_22(iN_X)
        V_u_3(iN_X) = V_d_3(iN_X) / Gm_dd_33(iN_X)
      END DO

      ! --- Compute E from Ef ---

      DO iN_X = 1, nX_G
        ! E(iN_X) = Ef(iN_X) / D(iN_X) - Half * (V_u_1(iN_X)*V_d_1(iN_X) &
        !                      + V_u_2(iN_X)*V_d_2(iN_X) + V_u_3(iN_X)*V_d_3(iN_X))
        NormVsq(iN_X) =  V_u_1(iN_X) * V_d_1(iN_X) &
                       + V_u_2(iN_X) * V_d_2(iN_X) &
                       + V_u_3(iN_X) * V_d_3(iN_X)
        E(iN_X) = Ef(iN_X) - Half * NormVsq(iN_X)

      END DO

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

       ! PRINT*, "after update = "
       ! PRINT*, "D = ", D(:) / (Gram / Centimeter**3)
       ! PRINT*, "T = ", T(:) / (MeV)
       ! PRINT*, "Y = ", Y(:)
       ! PRINT*, "Ef = ", Ef(:)
       ! PRINT*, "E = ", E(:)
       ! PRINT*, "V_d_1 = ", V_d_1(:)
       ! PRINT*, "V_d_2 = ", V_d_2(:)
       ! PRINT*, "V_d_3 = ", V_d_3(:)

      ! --- Check Convergence (outer) ---

      CALL CheckConvergenceOuter &
             ( ITERATE_OUTER, ITERATE_INNER, n_FP_outer, k_outer, &
               FVECm_outer, Rtol, nIterations_Outer )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    nIterations_Inner &
      = FLOOR( DBLE( nIterations_Inner ) / DBLE( nIterations_Outer ) )

    CALL ArrayCopy( Jnew, H1new, H2new, H3new, J, H_d_1, H_d_2, H_d_3 )

    DO iN_X = 1, nX_G
      H_u_1(:,iN_X,:) = H_d_1(:,iN_X,:) / Gm_dd_11(iN_X)
      H_u_2(:,iN_X,:) = H_d_2(:,iN_X,:) / Gm_dd_22(iN_X)
      H_u_3(:,iN_X,:) = H_d_3(:,iN_X,:) / Gm_dd_33(iN_X)
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE_OUTER, ITERATE_INNER, &
    !$OMP               Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP               Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP               Jold_1, Jnew_1, Jnorm_1, &
    !$OMP               Jold_2, Jnew_2, Jnorm_2, &
    !$OMP               Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP               Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP               Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP               PackIndex_outer, UnpackIndex_outer, &
    !$OMP               PackIndex_inner, UnpackIndex_inner, &
    !$OMP               AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$OMP               GVECm_outer, FVECm_outer, Alpha_outer, &
    !$OMP               AMAT_inner, BVEC_inner, GVEC_inner, FVEC_inner, &
    !$OMP               GVECm_inner, FVECm_inner, Alpha_inner )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE_OUTER, ITERATE_INNER, &
    !$ACC         Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$ACC         GVECm_outer, FVECm_outer, Alpha_outer, &
    !$ACC         AMAT_inner, BVEC_inner, GVEC_inner, FVEC_inner, &
    !$ACC         GVECm_inner, FVECm_inner, Alpha_inner )
#endif

  END SUBROUTINE SolveMatterEquations_FP_NestedAA


  SUBROUTINE ComputeEquilibriumDistributions_Packed &
    ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
      MASK, nX_P, PackIndex, UnpackIndex )

    INTEGER,                                      INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(1:nX_G),        TARGET,   INTENT(in)    :: D, T, Y
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), TARGET,   INTENT(inout) :: J0_1, J0_2

    LOGICAL,  DIMENSION(1:nX_G),        OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                            OPTIONAL, INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(1:nX_G),        OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:),   POINTER :: D_P, T_P, Y_P
    REAL(DP), DIMENSION(:,:), POINTER :: J0_1_P, J0_2_P

    INTEGER :: nX

    IF ( PRESENT( nX_P ) ) THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      D_P => P1D(1:nX,iP1D_D)
      T_P => P1D(1:nX,iP1D_T)
      Y_P => P1D(1:nX,iP1D_Y)

      CALL ArrayPack &
             ( nX, UnpackIndex, &
               D, T, Y, D_P, T_P, Y_P )

      J0_1_P => P2D(:,1:nX,iP2D_J0_1)
      J0_2_P => P2D(:,1:nX,iP2D_J0_2)

    ELSE

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      J0_1_P => J0_1(:,:)
      J0_2_P => J0_2(:,:)

    END IF

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, J0_2_P, iS_1, iS_2 )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               J0_1_P, J0_2_P, J0_1, J0_2 )

    END IF

  END SUBROUTINE ComputeEquilibriumDistributions_Packed


  SUBROUTINE ComputeOpacities_Packed &
    ( D, T, Y, J0, Chi, Sig, Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, Phi_0_Ot_Pair, &
      MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(1:nX_G),                          &
                          TARGET,   INTENT(in)    :: D, T, Y
    REAL(DP), DIMENSION(       1:nE_G,1:nX_G,1:nSpecies), &
                          TARGET,   INTENT(inout) :: J0, Chi, Sig
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies), &
                          TARGET,   INTENT(inout) :: Phi_0_In_NES, Phi_0_Ot_NES
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies), &
                          TARGET,   INTENT(inout) :: Phi_0_In_Pair, Phi_0_Ot_Pair

    LOGICAL,  DIMENSION(1:nX_G), OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                     OPTIONAL, INTENT(in)    :: nX_P, nX_P0
    INTEGER,  DIMENSION(1:nX_G), OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:),     POINTER :: D_P, T_P, Y_P

    ! --- to be changed to handle cases when iSpecies > 2 ---
    REAL(DP), DIMENSION(:,:),   POINTER :: J0_1_P, J0_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Chi_1_P, Chi_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Sig_1_P, Sig_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P
    ! --- to be changed to handle cases when iSpecies > 2 ---


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

      D_P => P1D(1:nX,iP1D_D)
      T_P => P1D(1:nX,iP1D_T)
      Y_P => P1D(1:nX,iP1D_Y)

      CALL ArrayPack &
             ( nX, UnpackIndex, &
               D, T, Y, D_P, T_P, Y_P )

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

      J0_1_P => J0(:,:,iS_1)
      J0_2_P => J0(:,:,iS_2)
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

    ! PRINT*, "D = ", D_P(:) / (Gram / Centimeter**3)
    ! PRINT*, "T = ", T_P(:) / (MeV)
    ! PRINT*, "Y = ", Y_P(:)

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, J0_2_P, iS_1, iS_2 )

    ! --- EmAb ---
    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, Chi_1_P )

    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_2, Chi_2_P )

    ! --- Isoenergetic scattering---
    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, 1, Sig_1_P )

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_2, 1, Sig_2_P )

    ! --- NES Kernels ---

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, iS_2, 1, &
             Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, &
             Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
             P3D(:,:,:,iP3D_WORK1), P3D(:,:,:,iP3D_WORK2) )

    ! --- Pair Kernels ---

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, iS_2, 1, &
             Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
             Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
             P3D(:,:,:,iP3D_WORK1), P3D(:,:,:,iP3D_WORK2) )

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
                 Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
                 Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
                 Phi_0_In_NES(:,:,:,iS_1), Phi_0_Ot_NES(:,:,:,iS_1), &
                 Phi_0_In_NES(:,:,:,iS_2), Phi_0_Ot_NES(:,:,:,iS_2), &
                 Phi_0_In_Pair(:,:,:,iS_1), Phi_0_Ot_Pair(:,:,:,iS_1), &
                 Phi_0_In_Pair(:,:,:,iS_2), Phi_0_Ot_Pair(:,:,:,iS_2) )

      END IF

    END IF

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( J, &
      Phi_0_In_NES, Phi_0_Ot_NES, Phi_0_In_Pair, Phi_0_Ot_Pair, &
      Chi_NES, Eta_NES, Chi_Pair, Eta_Pair, &
      MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies),        TARGET,   INTENT(in)    :: J
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies), TARGET,   INTENT(in)    :: Phi_0_In_NES, Phi_0_Ot_NES
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G,1:nSpecies), TARGET,   INTENT(in)    :: Phi_0_In_Pair, Phi_0_Ot_Pair
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies),        TARGET,   INTENT(inout) :: Chi_NES
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies),        TARGET,   INTENT(inout) :: Eta_NES
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies),        TARGET,   INTENT(inout) :: Chi_Pair
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies),        TARGET,   INTENT(inout) :: Eta_Pair

    LOGICAL,  DIMENSION(1:nX_G),               OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                                   OPTIONAL, INTENT(in)    :: nX_P, nX_P0
    INTEGER,  DIMENSION(1:nX_G),               OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

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
             ( nX, UnpackIndex, &
               J(:,:,iS_1), J(:,:,iS_2), J_1_P, J_2_P )

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
                 Phi_0_In_NES(:,:,:,iS_1), Phi_0_Ot_NES(:,:,:,iS_1), &
                 Phi_0_In_NES(:,:,:,iS_2), Phi_0_Ot_NES(:,:,:,iS_2), &
                 Phi_0_In_Pair(:,:,:,iS_1), Phi_0_Ot_Pair(:,:,:,iS_1), &
                 Phi_0_In_Pair(:,:,:,iS_2), Phi_0_Ot_Pair(:,:,:,iS_2), &
                 Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
                 Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P )

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
               Chi_NES_1_P, Chi_NES_2_P, Eta_NES_1_P, Eta_NES_2_P, &
               Chi_Pair_1_P, Chi_Pair_2_P, Eta_Pair_1_P, Eta_Pair_2_P, &
               Chi_NES(:,:,iS_1), Chi_NES(:,:,iS_2), Eta_NES(:,:,iS_1), Eta_NES(:,:,iS_2), &
               Chi_Pair(:,:,iS_1), Chi_Pair(:,:,iS_2), Eta_Pair(:,:,iS_1), Eta_Pair(:,:,iS_2) )

    END IF

  END SUBROUTINE ComputeRates_Packed


  SUBROUTINE UpdateTemperature_Packed &
    ( D, E, Y, T, MASK, nX_P, PackIndex, UnpackIndex )

    REAL(DP), DIMENSION(1:nX_G), TARGET,   INTENT(in)    :: D, E, Y
    REAL(DP), DIMENSION(1:nX_G), TARGET,   INTENT(inout) :: T

    LOGICAL,  DIMENSION(1:nX_G), OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                     OPTIONAL, INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(1:nX_G), OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:), POINTER :: D_P, E_P, Y_P, T_P

    INTEGER  :: nX

    IF ( PRESENT( nX_P ) ) THEN
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

      CALL ArrayPack &
             ( nX, UnpackIndex, D, Y, T, E, &
               D_P, Y_P, T_P, E_P )
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
    ( J, Jold, CJ, Jnew, &
      H_d_1, H1old, CH1, H1new, &
      H_d_2, H2old, CH2, H2new, &
      H_d_3, H3old, CH3, H3new, &
      D, Y, Yold, C_Y, S_Y, U_Y, &
      Ef, Efold, C_Ef, S_Ef, U_Ef, &
      V_d_1, V_d_1old, C_V_d_1, S_V_d_1, U_V_d_1, &
      V_d_2, V_d_2old, C_V_d_2, S_V_d_2, U_V_d_2, &
      V_d_3, V_d_3old, C_V_d_3, S_V_d_3, U_V_d_3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)  :: J, Jold
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(out) :: CJ, Jnew
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)  :: H_d_1, H1old
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(out) :: CH1, H1new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)  :: H_d_2, H2old
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(out) :: CH2, H2new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)  :: H_d_3, H3old
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(out) :: CH3, H3new


    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: D
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: Y, Yold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: Ef, Efold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_Ef, S_Ef, U_Ef
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: V_d_1, V_d_1old
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_V_d_1, S_V_d_1, U_V_d_1
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: V_d_2, V_d_2old
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_V_d_2, S_V_d_2, U_V_d_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: V_d_3, V_d_3old
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_V_d_3, S_V_d_3, U_V_d_3

    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: vHold, V1Jold, V2Jold, V3Jold
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies,3) :: vKJold

    REAL(DP) :: V_u_1, V_u_2, V_u_3, k_dd(3,3)

    INTEGER  :: iN_E, iN_X


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      V_u_1 = V_d_1old(iN_X)/Gm_dd_11(iN_X)
      V_u_2 = V_d_2old(iN_X)/Gm_dd_22(iN_X)
      V_u_3 = V_d_3old(iN_X)/Gm_dd_33(iN_X)

      vHold(:,iN_X,:) =  V_u_1 * H1old(:,iN_X,:) &
                       + V_u_2 * H2old(:,iN_X,:) &
                       + V_u_3 * H3old(:,iN_X,:)

      V1Jold(:,iN_X,:) = V_d_1old(iN_X)*Jold(:,iN_X,:)
      V2Jold(:,iN_X,:) = V_d_2old(iN_X)*Jold(:,iN_X,:)
      V3Jold(:,iN_X,:) = V_d_3old(iN_X)*Jold(:,iN_X,:)

      DO iN_E = 1, nE_G

        ! PRINT*, "iS = ", iS_1
        ! PRINT*, "Jold = ", Jold(iN_E,iN_X,iS_1)
        ! PRINT*, "H1old = ", H1old(iN_E,iN_X,iS_1)
        ! PRINT*, "H2old = ", H2old(iN_E,iN_X,iS_1)
        ! PRINT*, "H3old = ", H3old(iN_E,iN_X,iS_1)
        !
        ! PRINT*, "G11 = ", Gm_dd_11(iN_X)
        ! PRINT*, "G22 = ", Gm_dd_22(iN_X)
        ! PRINT*, "G33 = ", Gm_dd_33(iN_X)

        k_dd = EddingtonTensorComponents_dd &
                    ( Jold(iN_E,iN_X,iS_1), &
                     H1old(iN_E,iN_X,iS_1), H2old(iN_E,iN_X,iS_1), H3old(iN_E,iN_X,iS_1), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vKJold(iN_E,iN_X,iS_1,1) = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
        vKJold(iN_E,iN_X,iS_1,2) = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
        vKJold(iN_E,iN_X,iS_1,3) = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)
        vKJold(iN_E,iN_X,iS_1,:) = vKJold(iN_E,iN_X,iS_1,:) * Jold(iN_E,iN_X,iS_1)
        ! PRINT*, "iS = ", iS_2
        ! PRINT*, "Jold = ", Jold(iN_E,iN_X,iS_2)
        ! PRINT*, "H1old = ", H1old(iN_E,iN_X,iS_2)
        ! PRINT*, "H2old = ", H2old(iN_E,iN_X,iS_2)
        ! PRINT*, "H3old = ", H3old(iN_E,iN_X,iS_2)
        !
        ! PRINT*, "G11 = ", Gm_dd_11(iN_X)
        ! PRINT*, "G22 = ", Gm_dd_22(iN_X)
        ! PRINT*, "G33 = ", Gm_dd_33(iN_X)

        k_dd = EddingtonTensorComponents_dd &
                    ( Jold(iN_E,iN_X,iS_2), &
                     H1old(iN_E,iN_X,iS_2), H2old(iN_E,iN_X,iS_2), H3old(iN_E,iN_X,iS_2), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vKJold(iN_E,iN_X,iS_2,1) = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
        vKJold(iN_E,iN_X,iS_2,2) = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
        vKJold(iN_E,iN_X,iS_2,3) = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)
        vKJold(iN_E,iN_X,iS_2,:) = vKJold(iN_E,iN_X,iS_2,:) * Jold(iN_E,iN_X,iS_2)

      END DO

       CJ(:,iN_X,:) =  Jold(:,iN_X,:) + vHold(:,iN_X,:)
      CH1(:,iN_X,:) = H1old(:,iN_X,:) + vKJold(:,iN_X,:,1)
      CH2(:,iN_X,:) = H2old(:,iN_X,:) + vKJold(:,iN_X,:,2)
      CH3(:,iN_X,:) = H3old(:,iN_X,:) + vKJold(:,iN_X,:,3)

    END DO



    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,  Jold(:,:,iS_1), nE_G, W2_S, 1, Zero, C_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One,  Jold(:,:,iS_2), nE_G, W2_S, 1,  One, C_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vHold(:,:,iS_1), nE_G, W2_S, 1,  One, C_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, vHold(:,:,iS_2), nE_G, W2_S, 1,  One, C_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,  Jold(:,:,iS_1), nE_G, W3_S, 1, Zero, C_Ef, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,  Jold(:,:,iS_2), nE_G, W3_S, 1,  One, C_Ef, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +Two, vHold(:,:,iS_1), nE_G, W3_S, 1,  One, C_Ef, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +Two, vHold(:,:,iS_2), nE_G, W3_S, 1,  One, C_Ef, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V1Jold(:,:,iS_1), nE_G, W3_S, 1, Zero, C_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V1Jold(:,:,iS_2), nE_G, W3_S, 1,  One, C_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H1old(:,:,iS_1), nE_G, W3_S, 1,  One, C_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H1old(:,:,iS_2), nE_G, W3_S, 1,  One, C_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJold(:,:,iS_1,1), nE_G, W3_S, 1,  One, C_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJold(:,:,iS_2,1), nE_G, W3_S, 1,  One, C_V_d_1, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V2Jold(:,:,iS_1), nE_G, W3_S, 1, Zero, C_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V2Jold(:,:,iS_2), nE_G, W3_S, 1,  One, C_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H2old(:,:,iS_1), nE_G, W3_S, 1,  One, C_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H2old(:,:,iS_2), nE_G, W3_S, 1,  One, C_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJold(:,:,iS_1,2), nE_G, W3_S, 1,  One, C_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJold(:,:,iS_2,2), nE_G, W3_S, 1,  One, C_V_d_2, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V3Jold(:,:,iS_1), nE_G, W3_S, 1, Zero, C_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V3Jold(:,:,iS_2), nE_G, W3_S, 1,  One, C_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H3old(:,:,iS_1), nE_G, W3_S, 1,  One, C_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H3old(:,:,iS_2), nE_G, W3_S, 1,  One, C_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJold(:,:,iS_1,3), nE_G, W3_S, 1,  One, C_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJold(:,:,iS_2,3), nE_G, W3_S, 1,  One, C_V_d_3, 1 )


#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      S_Y(iN_X) = One / ( D(iN_X) * Yold(iN_X) / AtomicMassUnit )
      ! Ef may be negative, may need to change scaling
      ! A potential scaling is the lower bound in the table
      S_Ef(iN_X) = One / ( D(iN_X) * Efold(iN_X) )
      S_V_d_1(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_2(iN_X) = One / ( D(iN_X) * SpeedOfLight )
      S_V_d_3(iN_X) = One / ( D(iN_X) * SpeedOfLight )

      ! Initial guess is the old matter state
      U_Y(iN_X)     = Y(iN_X)     / Yold(iN_X)
      U_Ef(iN_X)    = Ef(iN_X)    / Efold(iN_X)
      U_V_d_1(iN_X) = V_d_1(iN_X) / SpeedOfLight
      U_V_d_2(iN_X) = V_d_2(iN_X) / SpeedOfLight
      U_V_d_3(iN_X) = V_d_3(iN_X) / SpeedOfLight

      ! U_V_d_1(iN_X) = V_d_1(iN_X) / V_d_1old(iN_X) ! --- Initial Guess
      ! U_V_d_2(iN_X) = V_d_2(iN_X) / V_d_2old(iN_X) ! --- Initial Guess
      ! U_V_d_3(iN_X) = V_d_3(iN_X) / V_d_3old(iN_X) ! --- Initial Guess

      ! Now the old matter state is included in C terms
      C_Y(iN_X)     = Yold(iN_X)     / Yold(iN_X)   + C_Y(iN_X)     * S_Y(iN_X)
      C_Ef(iN_X)    = Efold(iN_X)    / Efold(iN_X)  + C_Ef(iN_X)    * S_Ef(iN_X)
      C_V_d_1(iN_X) = V_d_1old(iN_X) / SpeedOfLight + C_V_d_1(iN_X) * S_V_d_1(iN_X)
      C_V_d_2(iN_X) = V_d_2old(iN_X) / SpeedOfLight + C_V_d_2(iN_X) * S_V_d_2(iN_X)
      C_V_d_3(iN_X) = V_d_3old(iN_X) / SpeedOfLight + C_V_d_3(iN_X) * S_V_d_3(iN_X)

    END DO

    CALL ArrayCopy( J, H_d_1, H_d_2, H_d_3, Jnew, H1new, H2new, H3new )

  END SUBROUTINE InitializeRHS_FP


  SUBROUTINE ComputeJNorm &
    ( MASK, J, Jnorm )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: J
    REAL(DP), DIMENSION(1:nX_G,1:nSpecies),        INTENT(inout) :: Jnorm

    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN
        !!! --- TO DO: Change to weighted norm ---
        Jnorm(iN_X,iS_1) = SQRT( SUM( J(:,iN_X,iS_1)**2 ) )
        Jnorm(iN_X,iS_2) = SQRT( SUM( J(:,iN_X,iS_2)**2 ) )
      END IF
    END DO

  END SUBROUTINE ComputeJNorm



  SUBROUTINE ComputeMatterRHS_FP &
    ( MASK, n_FPVar, Fm, Gm, J, H1, H2, H3, &
     C_Y, S_Y, U_Y, G_Y,&
     C_Ef, S_Ef, U_Ef, G_Ef, &
     V_d_1, C_V_d_1, S_V_d_1, U_V_d_1, G_V_d_1, &
     V_d_2, C_V_d_2, S_V_d_2, U_V_d_2, G_V_d_2, &
     V_d_3, C_V_d_3, S_V_d_3, U_V_d_3, G_V_d_3, &
     Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FPVar
    REAL(DP), DIMENSION(1:n_FPVar,1:nX_G), INTENT(inout) :: Fm, Gm

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: J
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: H1, H2, H3

    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: C_Ef, S_Ef, U_Ef
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: V_d_1, C_V_d_1, S_V_d_1, U_V_d_1
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: V_d_2, C_V_d_2, S_V_d_2, U_V_d_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: V_d_3, C_V_d_3, S_V_d_3, U_V_d_3
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out)   :: G_Y, G_Ef, G_V_d_1, G_V_d_2, G_V_d_3

    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies) :: vH, V1J, V2J, V3J
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies,3) :: vKJ

    REAL(DP) :: V_u_1, V_u_2, V_u_3, k_dd(3,3)

    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      V_u_1 = V_d_1(iN_X)/Gm_dd_11(iN_X)
      V_u_2 = V_d_2(iN_X)/Gm_dd_22(iN_X)
      V_u_3 = V_d_3(iN_X)/Gm_dd_33(iN_X)

      vH(:,iN_X,:) =  V_u_1 * H1(:,iN_X,:) &
                    + V_u_2 * H2(:,iN_X,:) &
                    + V_u_3 * H3(:,iN_X,:)

      V1J(:,iN_X,:) = V_d_1(iN_X)*J(:,iN_X,:)
      V2J(:,iN_X,:) = V_d_2(iN_X)*J(:,iN_X,:)
      V3J(:,iN_X,:) = V_d_3(iN_X)*J(:,iN_X,:)

      DO iN_E = 1, nE_G

        k_dd = EddingtonTensorComponents_dd &
                    ( J(iN_E,iN_X,iS_1), &
                     H1(iN_E,iN_X,iS_1), H2(iN_E,iN_X,iS_1), H3(iN_E,iN_X,iS_1), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vKJ(iN_E,iN_X,iS_1,1) = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
        vKJ(iN_E,iN_X,iS_1,2) = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
        vKJ(iN_E,iN_X,iS_1,3) = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)
        vKJ(iN_E,iN_X,iS_1,:) = vKJ(iN_E,iN_X,iS_1,:) * J(iN_E,iN_X,iS_1)

        k_dd = EddingtonTensorComponents_dd &
                    ( J(iN_E,iN_X,iS_2), &
                     H1(iN_E,iN_X,iS_2), H2(iN_E,iN_X,iS_2), H3(iN_E,iN_X,iS_2), &
                     Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )

        vKJ(iN_E,iN_X,iS_2,1) = V_u_1 * k_dd(1,1) + V_u_2 * k_dd(1,2) + V_u_3 * k_dd(1,3)
        vKJ(iN_E,iN_X,iS_2,2) = V_u_1 * k_dd(2,1) + V_u_2 * k_dd(2,2) + V_u_3 * k_dd(2,3)
        vKJ(iN_E,iN_X,iS_2,3) = V_u_1 * k_dd(3,1) + V_u_2 * k_dd(3,2) + V_u_3 * k_dd(3,3)
        vKJ(iN_E,iN_X,iS_2,:) = vKJ(iN_E,iN_X,iS_2,:) * J(iN_E,iN_X,iS_2)
      END DO

    END DO


    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J(:,:,iS_1), nE_G, W2_S, 1, Zero, G_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, J(:,:,iS_2), nE_G, W2_S, 1,  One, G_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vH(:,:,iS_1), nE_G, W2_S, 1,  One, G_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, vH(:,:,iS_2), nE_G, W2_S, 1,  One, G_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J(:,:,iS_1), nE_G, W3_S, 1, Zero, G_Ef, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J(:,:,iS_2), nE_G, W3_S, 1,  One, G_Ef, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +Two, vH(:,:,iS_1), nE_G, W3_S, 1,  One, G_Ef, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +Two, vH(:,:,iS_2), nE_G, W3_S, 1,  One, G_Ef, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V1J(:,:,iS_1), nE_G, W3_S, 1, Zero, G_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V1J(:,:,iS_2), nE_G, W3_S, 1,  One, G_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H1(:,:,iS_1), nE_G, W3_S, 1,  One, G_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H1(:,:,iS_2), nE_G, W3_S, 1,  One, G_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJ(:,:,iS_1,1), nE_G, W3_S, 1,  One, G_V_d_1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJ(:,:,iS_2,1), nE_G, W3_S, 1,  One, G_V_d_1, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V2J(:,:,iS_1), nE_G, W3_S, 1, Zero, G_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V2J(:,:,iS_2), nE_G, W3_S, 1,  One, G_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H2(:,:,iS_1), nE_G, W3_S, 1,  One, G_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H2(:,:,iS_2), nE_G, W3_S, 1,  One, G_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJ(:,:,iS_1,2), nE_G, W3_S, 1,  One, G_V_d_2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJ(:,:,iS_2,2), nE_G, W3_S, 1,  One, G_V_d_2, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V3J(:,:,iS_1), nE_G, W3_S, 1, Zero, G_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,   V3J(:,:,iS_2), nE_G, W3_S, 1,  One, G_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H3(:,:,iS_1), nE_G, W3_S, 1,  One, G_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One,    H3(:,:,iS_2), nE_G, W3_S, 1,  One, G_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJ(:,:,iS_1,3), nE_G, W3_S, 1,  One, G_V_d_3, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, vKJ(:,:,iS_2,3), nE_G, W3_S, 1,  One, G_V_d_3, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_Y(iN_X)      = C_Y(iN_X)     - G_Y(iN_X)     * S_Y(iN_X)
        G_Ef(iN_X)     = C_Ef(iN_X)    - G_Ef(iN_X)    * S_Ef(iN_X)
        G_V_d_1(iN_X)  = C_V_d_1(iN_X) - G_V_d_1(iN_X) * S_V_d_1(iN_X)
        G_V_d_2(iN_X)  = C_V_d_2(iN_X) - G_V_d_2(iN_X) * S_V_d_2(iN_X)
        G_V_d_3(iN_X)  = C_V_d_3(iN_X) - G_V_d_3(iN_X) * S_V_d_3(iN_X)

        Gm( iY,iN_X) = G_Y(iN_X)
        Gm(iEf,iN_X) = G_Ef(iN_X)
        Gm(iV1,iN_X) = G_V_d_1(iN_X)
        Gm(iV2,iN_X) = G_V_d_2(iN_X)
        Gm(iV3,iN_X) = G_V_d_3(iN_X)


        Fm( iY,iN_X) = G_Y(iN_X)     - U_Y(iN_X)
        Fm(iEf,iN_X) = G_Ef(iN_X)    - U_Ef(iN_X)
        Fm(iV1,iN_X) = G_V_d_1(iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = G_V_d_2(iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = G_V_d_3(iN_X) - U_V_d_3(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_FP

  SUBROUTINE ComputeNeutrinoRHS_FP &
    ( MASK, n_FPVar, Fm, Gm, &
      dt, Omega, V_u_1, V_u_2, V_u_3, &
      CJ, CH1, CH2, CH3, &
      Jnew, H1new, H2new, H3new, &
      J0, Chi, Sig, &
      Chi_NES, Eta_NES, &
      Chi_Pair, Eta_Pair, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    LOGICAL,  DIMENSION(1:nX_G),                   INTENT(in)    :: MASK
    INTEGER,                                       INTENT(in)    :: n_FPVar
    REAL(DP), DIMENSION(1:n_FPVar,1:nX_G),         INTENT(inout) :: Fm, Gm
    REAL(DP),                                      INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)    :: Omega
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)    :: V_u_1
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)    :: V_u_2
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)    :: V_u_3
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: CJ
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: CH1
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: CH2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: CH3
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Jnew
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: H1new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: H2new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: H3new
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: J0
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Chi
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Sig
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Chi_NES
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Eta_NES
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Chi_Pair
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(in)    :: Eta_Pair
    REAL(DP), DIMENSION(1:nX_G),                   INTENT(in)    :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: Eta, Eta_T, Chi_T, Kappa
    INTEGER  :: iN_E, iN_X
    REAL(DP) :: vH, vKJ(3), k_dd(3,3)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T, Kappa )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T, Kappa )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T, Kappa )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN
          vH =  V_u_1(iN_X) * H1new(iN_E,iN_X,iS_1) &
              + V_u_2(iN_X) * H2new(iN_E,iN_X,iS_1) &
              + V_u_3(iN_X) * H3new(iN_E,iN_X,iS_1)
          k_dd = EddingtonTensorComponents_dd &
                      ( Jnew(iN_E,iN_X,iS_1), &
                       H1new(iN_E,iN_X,iS_1), H2new(iN_E,iN_X,iS_1), H3new(iN_E,iN_X,iS_1), &
                       Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )
          vKJ(1) = V_u_1(iN_X) * k_dd(1,1) + V_u_2(iN_X) * k_dd(1,2) + V_u_3(iN_X) * k_dd(1,3)
          vKJ(2) = V_u_1(iN_X) * k_dd(2,1) + V_u_2(iN_X) * k_dd(2,2) + V_u_3(iN_X) * k_dd(2,3)
          vKJ(3) = V_u_1(iN_X) * k_dd(3,1) + V_u_2(iN_X) * k_dd(3,2) + V_u_3(iN_X) * k_dd(3,3)
          vKJ = vKJ * Jnew(iN_E,iN_X,iS_1)

          Eta = Chi(iN_E,iN_X,iS_1) * J0(iN_E,iN_X,iS_1)
          Eta_T = Eta + Eta_NES(iN_E,iN_X,iS_1) + Eta_Pair(iN_E,iN_X,iS_1)
          Chi_T = Chi(iN_E,iN_X,iS_1) + Chi_NES(iN_E,iN_X,iS_1) + Chi_Pair(iN_E,iN_X,iS_1)
          Kappa = Chi_T + Sig(iN_E,iN_X,iS_1)

          Gm(OS_JNuE+iN_E,iN_X) = (One - Omega(iN_X)) * Jnew(iN_E,iN_X,iS_1) &
                                  + Omega(iN_X) / ( One + dt * Chi_T ) &
                                  * ( CJ(iN_E,iN_X,iS_1) + dt * Eta_T - vH)
          Fm(OS_JNuE+iN_E,iN_X) = Gm(OS_JNuE+iN_E,iN_X) - Jnew(iN_E,iN_X,iS_1)

          Gm(OS_H1NuE+iN_E,iN_X) = (One - Omega(iN_X)) * H1new(iN_E,iN_X,iS_1) &
                                   + Omega(iN_X) / ( One + dt * Kappa ) &
                                   * ( CH1(iN_E,iN_X,iS_1) - vKJ(1))
          Fm(OS_H1NuE+iN_E,iN_X) = Gm(OS_H1NuE+iN_E,iN_X) - H1new(iN_E,iN_X,iS_1)
          Gm(OS_H2NuE+iN_E,iN_X) = (One - Omega(iN_X)) * H2new(iN_E,iN_X,iS_1) &
                                   + Omega(iN_X) / ( One + dt * Kappa ) &
                                   * ( CH2(iN_E,iN_X,iS_1) - vKJ(2))
          Fm(OS_H2NuE+iN_E,iN_X) = Gm(OS_H2NuE+iN_E,iN_X) - H2new(iN_E,iN_X,iS_1)
          Gm(OS_H3NuE+iN_E,iN_X) = (One - Omega(iN_X)) * H3new(iN_E,iN_X,iS_1) &
                                   + Omega(iN_X) / ( One + dt * Kappa ) &
                                   * ( CH3(iN_E,iN_X,iS_1) - vKJ(3))
          Fm(OS_H3NuE+iN_E,iN_X) = Gm(OS_H3NuE+iN_E,iN_X) - H3new(iN_E,iN_X,iS_1)


          vH =  V_u_1(iN_X) * H1new(iN_E,iN_X,iS_2) &
              + V_u_2(iN_X) * H2new(iN_E,iN_X,iS_2) &
              + V_u_3(iN_X) * H3new(iN_E,iN_X,iS_2)
          k_dd = EddingtonTensorComponents_dd &
                      ( Jnew(iN_E,iN_X,iS_2), &
                       H1new(iN_E,iN_X,iS_2), H2new(iN_E,iN_X,iS_2), H3new(iN_E,iN_X,iS_2), &
                       Gm_dd_11(iN_X), Gm_dd_22(iN_X), Gm_dd_33(iN_X) )
          vKJ(1) = V_u_1(iN_X) * k_dd(1,1) + V_u_2(iN_X) * k_dd(1,2) + V_u_3(iN_X) * k_dd(1,3)
          vKJ(2) = V_u_1(iN_X) * k_dd(2,1) + V_u_2(iN_X) * k_dd(2,2) + V_u_3(iN_X) * k_dd(2,3)
          vKJ(3) = V_u_1(iN_X) * k_dd(3,1) + V_u_2(iN_X) * k_dd(3,2) + V_u_3(iN_X) * k_dd(3,3)
          vKJ = vKJ * Jnew(iN_E,iN_X,iS_2)

          Eta = Chi(iN_E,iN_X,iS_2) * J0(iN_E,iN_X,iS_2)
          Eta_T = Eta + Eta_NES(iN_E,iN_X,iS_2) + Eta_Pair(iN_E,iN_X,iS_2)
          Chi_T = Chi(iN_E,iN_X,iS_2) + Chi_NES(iN_E,iN_X,iS_2) + Chi_Pair(iN_E,iN_X,iS_2)
          Kappa = Chi_T + Sig(iN_E,iN_X,iS_2)

          Gm(OS_JNuE_Bar+iN_E,iN_X) = (One - Omega(iN_X)) * Jnew(iN_E,iN_X,iS_2) &
                                      + Omega(iN_X) / ( One + dt * Chi_T ) &
                                      * ( CJ(iN_E,iN_X,iS_2) + dt * Eta_T - vH)
          Fm(OS_JNuE_Bar+iN_E,iN_X) = Gm(OS_JNuE_Bar+iN_E,iN_X) - Jnew(iN_E,iN_X,iS_2)

          Gm(OS_H1NuE_Bar+iN_E,iN_X) = (One - Omega(iN_X)) * H1new(iN_E,iN_X,iS_2) &
                                       + Omega(iN_X) / ( One + dt * Kappa ) &
                                       * ( CH1(iN_E,iN_X,iS_2) - vKJ(1))
          Fm(OS_H1NuE_Bar+iN_E,iN_X) = Gm(OS_H1NuE_Bar+iN_E,iN_X) - H1new(iN_E,iN_X,iS_2)
          Gm(OS_H2NuE_Bar+iN_E,iN_X) = (One - Omega(iN_X)) * H2new(iN_E,iN_X,iS_2) &
                                       + Omega(iN_X) / ( One + dt * Kappa ) &
                                       * ( CH2(iN_E,iN_X,iS_2) - vKJ(2))
          Fm(OS_H2NuE_Bar+iN_E,iN_X) = Gm(OS_H2NuE_Bar+iN_E,iN_X) - H2new(iN_E,iN_X,iS_2)
          Gm(OS_H3NuE_Bar+iN_E,iN_X) = (One - Omega(iN_X)) * H3new(iN_E,iN_X,iS_2) &
                                       + Omega(iN_X) / ( One + dt * Kappa ) &
                                       * ( CH3(iN_E,iN_X,iS_2) - vKJ(3))
          Fm(OS_H3NuE_Bar+iN_E,iN_X) = Gm(OS_H3NuE_Bar+iN_E,iN_X) - H3new(iN_E,iN_X,iS_2)

        END IF
      END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_FP



  SUBROUTINE UpdateMatterRHS_FP &
    ( MASK, n_FPVar, Y, Yold, U_Y, Ef, Efold, U_Ef, V_d_1, V_d_1old, U_V_d_1, &
      V_d_2, V_d_2old, U_V_d_2, V_d_3, V_d_3old, U_V_d_3, Fm, Gm )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FPVar

    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: Yold, Efold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: V_d_1old, V_d_2old, V_d_3old

    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: Y, Ef, V_d_1, V_d_2, V_d_3

    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: U_Y, U_Ef
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: U_V_d_1, U_V_d_2, U_V_d_3

    REAL(DP), DIMENSION(1:n_FPVar,1:nX_G), INTENT(inout) :: Fm, Gm

    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        Fm( iY,iN_X) = Gm( iY,iN_X) - U_Y(iN_X)
        Fm(iEf,iN_X) = Gm(iEf,iN_X) - U_Ef(iN_X)
        Fm(iV1,iN_X) = Gm(iV1,iN_X) - U_V_d_1(iN_X)
        Fm(iV2,iN_X) = Gm(iV2,iN_X) - U_V_d_2(iN_X)
        Fm(iV3,iN_X) = Gm(iV3,iN_X) - U_V_d_3(iN_X)

        U_Y(iN_X)     = Gm( iY,iN_X)
        U_Ef(iN_X)    = Gm(iEf,iN_X)
        U_V_d_1(iN_X) = Gm(iV1,iN_X)
        U_V_d_2(iN_X) = Gm(iV2,iN_X)
        U_V_d_3(iN_X) = Gm(iV3,iN_X)

        Y(iN_X)     = U_Y(iN_X)     * Yold(iN_X)
        Ef(iN_X)    = U_Ef(iN_X)    * Efold(iN_X)
        V_d_1(iN_X) = U_V_d_1(iN_X) * SpeedOfLight
        V_d_2(iN_X) = U_V_d_2(iN_X) * SpeedOfLight
        V_d_3(iN_X) = U_V_d_3(iN_X) * SpeedOfLight

        ! V_d_1(iN_X) = U_V_d_1(iN_X) * V_d_1old(iN_X)
        ! V_d_2(iN_X) = U_V_d_2(iN_X) * V_d_2old(iN_X)
        ! V_d_3(iN_X) = U_V_d_3(iN_X) * V_d_3old(iN_X)

      END IF
    END DO

  END SUBROUTINE UpdateMatterRHS_FP


  SUBROUTINE UpdateNeutrinoRHS_FP &
    ( MASK, n_FPVar, Fm, Gm, Jnew, H1new, H2new, H3new )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FPVar
    REAL(DP), DIMENSION(1:n_FPVar,1:nX_G), INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(1:nE_G,1:nX_G,1:nSpecies), INTENT(inout) :: Jnew, H1new, H2new, H3new

    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Fm(OS_JNuE     +iN_E,iN_X) = Gm(OS_JNuE     +iN_E,iN_X) -  Jnew(iN_E,iN_X,iS_1)
          Fm(OS_H1NuE    +iN_E,iN_X) = Gm(OS_H1NuE    +iN_E,iN_X) - H1new(iN_E,iN_X,iS_1)
          Fm(OS_H2NuE    +iN_E,iN_X) = Gm(OS_H2NuE    +iN_E,iN_X) - H2new(iN_E,iN_X,iS_1)
          Fm(OS_H3NuE    +iN_E,iN_X) = Gm(OS_H3NuE    +iN_E,iN_X) - H3new(iN_E,iN_X,iS_1)
          Fm(OS_JNuE_Bar +iN_E,iN_X) = Gm(OS_JNuE_Bar +iN_E,iN_X) -  Jnew(iN_E,iN_X,iS_2)
          Fm(OS_H1NuE_Bar+iN_E,iN_X) = Gm(OS_H1NuE_Bar+iN_E,iN_X) - H1new(iN_E,iN_X,iS_2)
          Fm(OS_H2NuE_Bar+iN_E,iN_X) = Gm(OS_H2NuE_Bar+iN_E,iN_X) - H2new(iN_E,iN_X,iS_2)
          Fm(OS_H3NuE_Bar+iN_E,iN_X) = Gm(OS_H3NuE_Bar+iN_E,iN_X) - H3new(iN_E,iN_X,iS_2)

           Jnew(iN_E,iN_X,iS_1) = Gm(OS_JNuE     +iN_E,iN_X)
          H1new(iN_E,iN_X,iS_1) = Gm(OS_H1NuE    +iN_E,iN_X)
          H2new(iN_E,iN_X,iS_1) = Gm(OS_H2NuE    +iN_E,iN_X)
          H3new(iN_E,iN_X,iS_1) = Gm(OS_H3NuE    +iN_E,iN_X)
           Jnew(iN_E,iN_X,iS_2) = Gm(OS_JNuE_Bar +iN_E,iN_X)
          H1new(iN_E,iN_X,iS_2) = Gm(OS_H1NuE_Bar+iN_E,iN_X)
          H2new(iN_E,iN_X,iS_2) = Gm(OS_H2NuE_Bar+iN_E,iN_X)
          H3new(iN_E,iN_X,iS_2) = Gm(OS_H3NuE_Bar+iN_E,iN_X)


        END IF
      END DO
    END DO

  END SUBROUTINE UpdateNeutrinoRHS_FP


  SUBROUTINE SolveLS_FP &
    ( MASK, n_FP, M, Mk, Fm, Gm, F, G, A, B, Alpha )

    LOGICAL,  DIMENSION(1:nX_G),            INTENT(in)    :: MASK
    INTEGER,                                INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(1:n_FP,1:nX_G),     INTENT(inout) :: Fm, Gm, B
    REAL(DP), DIMENSION(1:n_FP,1:M,1:nX_G), INTENT(inout) :: F, G, A
    REAL(DP), DIMENSION(1:M,1:nX_G),        INTENT(inout) :: Alpha

    REAL(DP) :: AA11, AA12, AA21, AA22, AB1, AB2, DET_AA, SUM1
    INTEGER  :: iN_X, iFP, iM

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
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

    IF ( Mk > 1 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
      !$OMP DO SIMD COLLAPSE(2)
#endif
      DO iN_X = 1, nX_G
        DO iFP = 1, n_FP
          IF ( MASK(iN_X) ) THEN
            B(iFP,iN_X) = - Fm(iFP,iN_X)
          END IF
        END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
      !$OMP DO SIMD COLLAPSE(3)
#endif
      DO iN_X = 1, nX_G
        DO iM = 1, Mk-1
          DO iFP = 1, n_FP
            IF ( MASK(iN_X) ) THEN
              A(iFP,iM,iN_X) = F(iFP,iM,iN_X) - Fm(iFP,iN_X)
            END IF
          END DO
        END DO
      END DO

      IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
        !$OMP TARGET TEAMS DISTRIBUTE &
        !$OMP PRIVATE( AA11, AB1 )
#elif defined(THORNADO_OACC)
        !$ACC PARALLEL LOOP GANG &
        !$ACC PRIVATE( AA11, AB1 )
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

              AA11 = AA11 + A(iFP,1,iN_X) * A(iFP,1,iN_X)
              AB1  = AB1  + A(iFP,1,iN_X) * B(iFP,iN_X)

            END DO

            B(1,iN_X) = AB1 / AA11

          END IF
        END DO

      ELSE IF ( Mk == 3 ) THEN

        IF ( n_FP == 2 ) THEN

#if defined(THORNADO_OMP_OL)
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
          !$OMP PRIVATE( AA11, AA12, AA21, AA22, AB1, AB2, DET_AA )
#elif defined(THORNADO_OACC)
          !$ACC PARALLEL LOOP GANG VECTOR &
          !$ACC PRIVATE( AA11, AA12, AA21, AA22, AB1, AB2, DET_AA )
#elif defined(THORNADO_OMP)
          !$OMP PARALLEL DO SIMD &
          !$OMP PRIVATE( AA11, AA12, AA21, AA22, AB1, AB2, DET_AA )
#endif
          DO iN_X = 1, nX_G
            IF ( MASK(iN_X) ) THEN

              AA11 = A(1,1,iN_X)
              AA21 = A(2,1,iN_X)
              AA12 = A(1,2,iN_X)
              AA22 = A(2,2,iN_X)

              AB1 = B(1,iN_X)
              AB2 = B(2,iN_X)

              DET_AA = AA11*AA22 - AA21*AA12

              B(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
              B(2,iN_X) = ( - AA21 * AB1 + AA11 * AB2 ) / DET_AA

            END IF
          END DO

        ELSE

#if defined(THORNADO_OMP_OL)
          !$OMP TARGET TEAMS DISTRIBUTE &
          !$OMP PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
#elif defined(THORNADO_OACC)
          !$ACC PARALLEL LOOP GANG &
          !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA )
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

                AA11 = AA11 + A(iFP,1,iN_X) * A(iFP,1,iN_X)
                AA12 = AA12 + A(iFP,1,iN_X) * A(iFP,2,iN_X)
                AA22 = AA22 + A(iFP,2,iN_X) * A(iFP,2,iN_X)

                AB1  = AB1  + A(iFP,1,iN_X) * B(iFP,iN_X)
                AB2  = AB2  + A(iFP,2,iN_X) * B(iFP,iN_X)

              END DO

              DET_AA = AA11*AA22 - AA12*AA12

              B(1,iN_X) = ( + AA22 * AB1 - AA12 * AB2 ) / DET_AA
              B(2,iN_X) = ( - AA12 * AB1 + AA11 * AB2 ) / DET_AA

            END IF
          END DO

        END IF

      ELSE IF ( Mk > 3 ) THEN

        DO iN_X = 1, nX_G
          IF ( MASK(iN_X) ) THEN

            CALL LinearLeastSquares &
              ( 'N', n_FP, Mk-1, 1, A(:,:,iN_X), n_FP, &
                B(:,iN_X), n_FP, TAU(1,iN_X), WORK(1,iN_X), LWORK, INFO(iN_X) )

          END IF
        END DO

      END IF

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( SUM1 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( SUM1 )
#endif
      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          SUM1 = Zero
          DO iM = 1, Mk-1
            Alpha(iM,iN_X) = B(iM,iN_X)
            SUM1 = SUM1 + B(iM,iN_X)
          END DO
          Alpha(Mk,iN_X) = One - SUM1

        END IF
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( SUM1 )
#elif defined(THORNADO_OMP)
      !$OMP DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( SUM1 )
#endif
      DO iN_X = 1, nX_G
        DO iFP = 1, n_FP
          IF ( MASK(iN_X) ) THEN

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



  SUBROUTINE ShiftRHS_FP &
    ( MASK, n_FPVar, M, Mk, F, G )

    LOGICAL,  DIMENSION(1:nX_G),            INTENT(in)    :: MASK
    INTEGER,                                INTENT(in)    :: n_FPVar, M, Mk
    REAL(DP), DIMENSION(1:n_FPVar,1:M,1:nX_G), INTENT(inout) :: F, G

    REAL(DP) :: FTMP(1:n_FPVar,1:M), GTMP(1:n_FPVar,1:M)
    INTEGER  :: iN_X, iFP, iM

    IF ( Mk == M ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( FTMP, GTMP )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( FTMP, GTMP )
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
          DO iM = 1, Mk-1
            DO iFP = 1, n_FPVar
              FTMP(iFP,iM) = F(iFP,iM+1,iN_X)
              GTMP(iFP,iM) = G(iFP,iM+1,iN_X)
            END DO
          END DO

#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM = 1, Mk-1
            DO iFP = 1, n_FPVar
              F(iFP,iM,iN_X) = FTMP(iFP,iM)
              G(iFP,iM,iN_X) = GTMP(iFP,iM)
            END DO
          END DO

        END IF
      END DO

    END IF

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE CheckConvergenceInner &
    ( MASK, n_FPVar, k_inner, Rtol, nIterations_Inner, Fm, Jnorm )

    LOGICAL,  DIMENSION(1:nX_G),         INTENT(inout)  :: MASK
    INTEGER,                                INTENT(in)  :: n_FPVar, k_inner
    REAL(DP),                               INTENT(in)  :: Rtol
    INTEGER,  DIMENSION(1:nX_G),         INTENT(inout)  :: nIterations_Inner
    REAL(DP), DIMENSION(1:n_FPVar,1:nX_G),  INTENT(in)  :: Fm
    REAL(DP), DIMENSION(1:nX_G,1:nSpecies), INTENT(in)  :: Jnorm

    REAL(DP) :: Fnorm_1, Fnorm_2
    LOGICAL  :: CONVERGED
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_1, Fnorm_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_1, Fnorm_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_1, Fnorm_2 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        Fnorm_1 = SQRT( SUM( Fm(OS_JNuE    +1:OS_JNuE     + 4*nE_G,iN_X)**2 ) )
        Fnorm_2 = SQRT( SUM( Fm(OS_JNuE_Bar+1:OS_JNuE_Bar + 4*nE_G,iN_X)**2 ) )

        CONVERGED = Fnorm_1 <= Rtol * Jnorm(iN_X,iS_1) .AND. &
                    Fnorm_2 <= Rtol * Jnorm(iN_X,iS_2)

        IF ( CONVERGED ) THEN
          MASK(iN_X) = .FALSE.
          nIterations_Inner(iN_X) = nIterations_Inner(iN_X) + k_inner
        END IF

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( MASK )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( MASK )
#endif

  END SUBROUTINE CheckConvergenceInner


  SUBROUTINE CheckConvergenceOuter &
    ( MASK_OUTER, MASK_INNER, n_FPVar, k_outer, Fm, Rtol, nIterations_Outer )

    LOGICAL,  DIMENSION(1:nX_G),           INTENT(inout) :: MASK_OUTER, MASK_INNER
    INTEGER,                               INTENT(in)    :: n_FPVar, k_outer
    REAL(DP), DIMENSION(1:n_FPVar,1:nX_G), INTENT(in)    :: Fm
    REAL(DP),                              INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(1:nX_G),           INTENT(inout) :: nIterations_Outer

    LOGICAL  :: CONVERGED
    REAL(DP) :: Fnorm_Y, Fnorm_Ef, Fnorm_V
    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK_OUTER(iN_X) ) THEN

        Fnorm_Y  = ABS( Fm(iY,iN_X) )
        Fnorm_Ef = ABS( Fm(iEf,iN_X) )
        !!! --- Check if this norm makes sense ---
        Fnorm_V  = SQRT( Fm(iV1,iN_X)**2 + Fm(iV2,iN_X)**2 + Fm(iV3,iN_X)**2 )
        CONVERGED = Fnorm_Y  <= Rtol .AND. &
                    Fnorm_Ef <= Rtol .AND. &
                    Fnorm_V  <= Rtol

        IF ( CONVERGED ) THEN
          MASK_OUTER(iN_X) = .FALSE.
          nIterations_Outer(iN_X) = k_outer
        END IF

        MASK_INNER(iN_X) = MASK_OUTER(iN_X)

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( MASK_OUTER, MASK_INNER )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( MASK_OUTER, MASK_INNER )
#endif

  END SUBROUTINE CheckConvergenceOuter



END MODULE TwoMoment_NeutrinoMatterSolverModule_OrderV
