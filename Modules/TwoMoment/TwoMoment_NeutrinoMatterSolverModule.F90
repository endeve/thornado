#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_IMPLICIT
#endif
MODULE TwoMoment_NeutrinoMatterSolverModule

  USE KindModule, ONLY: &
    DP, Zero, One, FourPi
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
    iNuE, iNuE_Bar
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions, &
    ComputeEquilibriumDistributionAndDerivatives, &
    ComputeNeutrinoOpacities_NES, &
    ComputeNeutrinoOpacities_Pair, &
    ComputeNeutrinoOpacityRates_NES, &
    ComputeNeutrinoOpacityRates_Pair

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: SolveMatterEquations_EmAb
  PUBLIC :: SolveMatterEquations_EmAb_NuE
  PUBLIC :: SolveMatterEquations_EmAb_FP
  PUBLIC :: SolveMatterEquations_FP_Coupled
  PUBLIC :: SolveMatterEquations_FP_NestedAA
  PUBLIC :: SolveMatterEquations_FP_NestedNewton

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

  LOGICAL, PARAMETER :: SolveMatter = .TRUE.
  LOGICAL, PARAMETER :: UsePreconditionerEmAb = .FALSE.
  LOGICAL, PARAMETER :: UsePreconditionerPair = .FALSE.
  LOGICAL, PARAMETER :: UsePreconditionerPairLagAllButJ0 = .FALSE.

  INTEGER,  PARAMETER :: M_FP = 3
  INTEGER,  PARAMETER :: M_outer = 2
  INTEGER,  PARAMETER :: M_inner = 3

  INTEGER  :: nE_G, nX_G, nZ(4)
  INTEGER  :: n_FP, n_FP_inner, n_FP_outer
  INTEGER  :: iE_B, iE_E
  INTEGER  :: iX_B(3), iX_E(3)

  REAL(DP), ALLOCATABLE :: E_N(:)        ! --- Energy Grid
  REAL(DP), ALLOCATABLE :: W2_N(:)       ! --- Ingegration Weights (E^2)
  REAL(DP), ALLOCATABLE :: W3_N(:)       ! --- Integration Weights (E^3)
  REAL(DP), ALLOCATABLE :: W2_S(:)       ! --- Ingegration Weights (E^2) Scaled by (hc)^3
  REAL(DP), ALLOCATABLE :: W3_S(:)       ! --- Ingegration Weights (E^3) Scaled by (hc)^3

  REAL(DP), ALLOCATABLE :: AMAT(:,:,:)
  REAL(DP), ALLOCATABLE :: BVEC(:,:)
  REAL(DP), ALLOCATABLE :: TAU (:,:)
  REAL(DP), ALLOCATABLE :: WORK(:,:)
  INTEGER,  ALLOCATABLE :: IPIV(:,:)
  INTEGER,  ALLOCATABLE :: INFO(:)
  INTEGER               :: LWORK

  ! --- Temporary arrays for scatter/gather (packing)

  REAL(DP), ALLOCATABLE :: ATMP_1(:,:,:)
  REAL(DP), ALLOCATABLE :: BTMP_1(:,:,:)
  REAL(DP), ALLOCATABLE :: ATMP_2(:,:,:)
  REAL(DP), ALLOCATABLE :: BTMP_2(:,:,:)

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
  INTEGER, PARAMETER :: iP2D_S_1        = 15
  INTEGER, PARAMETER :: iP2D_S_2        = 16

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

    n_FP       = 2 + 2*nE_G
    n_FP_inner = 2*nE_G
    n_FP_outer = 2

    ALLOCATE( E_N (nE_G) )
    ALLOCATE( W2_N(nE_G) )
    ALLOCATE( W3_N(nE_G) )
    ALLOCATE( W2_S(nE_G) )
    ALLOCATE( W3_S(nE_G) )

    CALL ComputePointsAndWeightsE( E_N, W2_N, W3_N )

    W2_S(:) = W2_N(:) / ( PlanckConstant * SpeedOfLight )**3
    W3_S(:) = W3_N(:) / ( PlanckConstant * SpeedOfLight )**3

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

    IF ( M_FP > 3 ) THEN

      CALL LinearLeastSquares_LWORK &
             ( 'N', n_FP, MAX(M_FP,M_outer,M_inner)-1, 1, &
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

#if defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_NEWTON)
    ALLOCATE( ATMP_1(nE_G,nE_G,nX_G) )
    ALLOCATE( BTMP_1(nE_G,nE_G,nX_G) )
    ALLOCATE( ATMP_2(nE_G,nE_G,nX_G) )
    ALLOCATE( BTMP_2(nE_G,nE_G,nX_G) )
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: ATMP_1, BTMP_1, ATMP_2, BTMP_2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( ATMP_1, BTMP_1, ATMP_2, BTMP_2 )
#endif
#endif

  END SUBROUTINE InitializeNeutrinoMatterSolver


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

#if defined(NEUTRINO_MATTER_SOLVER_FIXED_POINT_NESTED_NEWTON)
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ATMP_1, BTMP_1, ATMP_2, BTMP_2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ATMP_1, BTMP_1, ATMP_2, BTMP_2 )
#endif
    DEALLOCATE( ATMP_1, BTMP_1, ATMP_2, BTMP_2 )
#endif

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

      W2(iN_E) = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**2
      W3(iN_E) = FourPi * WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**3

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveMatterEquations_EmAb &
    ( J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

    ! --- Electron Neutrinos (1) and Electron Antineutrinos (2) ---

    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: D, T, Y, E

    CALL ComputeEquilibriumDistributions_Packed &
           ( iNuE, iNuE_Bar, D, T, Y, J0_1, J0_2 )

  END SUBROUTINE SolveMatterEquations_EmAb


  SUBROUTINE SolveMatterEquations_EmAb_NuE &
    ( dt, J, Chi, J0, D, T, Y, E, nIterations )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: J
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out)   :: J0
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(1:nX_G),        INTENT(out)   :: nIterations

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: MaxIter = 20
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Gam, GamJ, GamJ0
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: dGamJ0dY, dGamJ0dE
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: dJ0dE, dJ0dY

    REAL(DP), DIMENSION(1:nX_G) :: Yold, C_Y, F_Y, U_Y, dU_Y
    REAL(DP), DIMENSION(1:nX_G) :: Eold, C_E, F_E, U_E, dU_E

    REAL(DP), DIMENSION(1:2,1:2,1:nX_G) :: FJAC
    REAL(DP), DIMENSION(1:nX_G) :: FNRM0

    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex, UnpackIndex

    INTEGER  :: k, nX_P

    ITERATE(:) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE ) &
    !$OMP MAP( alloc: Gam, GamJ, GamJ0, dGamJ0dY, dGamJ0dE, &
    !$OMP             dJ0dE, dJ0dY, &
    !$OMP             Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$OMP             Eold, C_E, F_E, U_E, dU_E, &
    !$OMP             PackIndex, UnpackIndex, &
    !$OMP             FJAC, FNRM0 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( Gam, GamJ, GamJ0, dGamJ0dY, dGamJ0dE, &
    !$ACC         dJ0dE, dJ0dY, &
    !$ACC         Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$ACC         Eold, C_E, F_E, U_E, dU_E, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         FJAC, FNRM0 )
#endif

    ! --- Equilibrium Distribution and Derivatives ---

    CALL ComputeEquilibriumDistributionsAndDerivatives_Packed &
           ( iNuE, D, T, Y, J0, dJ0dY, dJ0dE )

    ! --- Initial RHS ---

    CALL InitializeRHS_EmAb_NuE &
           ( dt, J, Chi, J0, dJ0dY, dJ0dE, D, Y, E, &
             Gam, GamJ0, GamJ, dGamJ0dY, dGamJ0dE, &
             Yold, Eold, U_Y, U_E, C_Y, C_E, F_Y, F_E, FNRM0 )

    k = 0
    DO WHILE( ANY( ITERATE(:) ) .AND. k < MaxIter )

      k  = k + 1

      CALL CreatePackIndex &
             ( ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Jacobian ---

      CALL SolveLS_EmAb_NuE &
             ( ITERATE, D, C_Y, C_E, F_Y, F_E, dGamJ0dY, dGamJ0dE, &
               FJAC, dU_Y, dU_E, U_Y, U_E, Y, E )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Equilibrium Distribution and Derivatives ---

      CALL ComputeEquilibriumDistributionsAndDerivatives_Packed &
             ( iNuE, D, T, Y, J0, dJ0dY, dJ0dE, &
               ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Update Residuals and Solution Vectors ---

      CALL ComputeMatterRHS_EmAb_NuE &
             ( ITERATE, D, U_Y, U_E, C_Y, C_E, Gam, J0, dJ0dY, dJ0dE, &
               GamJ0, dGamJ0dY, dGamJ0dE, F_Y, F_E )

      ! --- Check for Convergence ---

      CALL CheckConvergence_EmAb_NuE &
             ( ITERATE, k, Rtol, Utol, nIterations, &
               F_Y, U_Y, dU_Y, F_E, U_E, dU_E, FNRM0 )

    END DO

    CALL ComputeNumberDensity_EmAb_NuE &
           ( dt, J, Chi, J0 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE, &
    !$OMP               Gam, GamJ, GamJ0, dGamJ0dY, dGamJ0dE, &
    !$OMP               dJ0dE, dJ0dY, &
    !$OMP               Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$OMP               Eold, C_E, F_E, U_E, dU_E, &
    !$OMP               PackIndex, UnpackIndex, &
    !$OMP               FJAC, FNRM0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE, &
    !$ACC         Gam, GamJ, GamJ0, dGamJ0dY, dGamJ0dE, &
    !$ACC         dJ0dE, dJ0dY, &
    !$ACC         Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$ACC         Eold, C_E, F_E, U_E, dU_E, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         FJAC, FNRM0 )
#endif

  END SUBROUTINE SolveMatterEquations_EmAb_NuE


  SUBROUTINE SolveMatterEquations_EmAb_FP &
    ( dt, iS_1, iS_2, J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, &
      D, T, Y, E, nIterations_Out, TOL )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                           INTENT(in)    :: dt
    INTEGER,                            INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(1:nX_G),        INTENT(out), OPTIONAL :: nIterations_Out
    REAL(DP),                           INTENT(in),  OPTIONAL :: TOL

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = iE
    INTEGER,  PARAMETER :: MaxIter = 100

    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E

    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_1, Jnew_1
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_2, Jnew_2
    REAL(DP), DIMENSION(       1:nX_G) :: Jnorm_1, Jnorm_2

    REAL(DP), DIMENSION(1:n_FP,1:M_FP,1:nX_G) :: GVEC, FVEC
    REAL(DP), DIMENSION(1:n_FP,       1:nX_G) :: GVECm, FVECm
    REAL(DP), DIMENSION(       1:M_FP,1:nX_G) :: Alpha

    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex, UnpackIndex

    INTEGER,  DIMENSION(1:nX_G) :: nIterations

    INTEGER  :: k, Mk, iN_X, nX_P
    INTEGER  :: OS_2

    REAL(DP) :: Rtol

    IF(PRESENT(TOL)) THEN
      Rtol = TOL
    ELSE
      Rtol = 1.0d-08
    END IF

    OS_2 = OS_1 + nE_G
    ITERATE(:) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE ) &
    !$OMP MAP( alloc: Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP             Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP             Jold_1, Jnew_1, Jnorm_1, &
    !$OMP             Jold_2, Jnew_2, Jnorm_2, &
    !$OMP             PackIndex, UnpackIndex, &
    !$OMP             GVEC, FVEC, GVECm, FVECm, Alpha, nIterations )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha, nIterations )
#endif

    CALL ArrayCopy( Y, E, Yold, Eold )

    CALL ArrayCopy( J_1, J_2, Jold_1, Jold_2 )

    CALL TimersStart( Timer_Im_ComputeOpacity )

    ! --- Compute Opacity Kernels ---

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_Packed &
           ( iS_1, iS_2, D, T, Y, J0_1, J0_2 )

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Initial RHS ---

    CALL InitializeRHS_FP &
           ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
             D, Y, Yold, C_Y, S_Y, Unew_Y, E, Eold, C_E, S_E, Unew_E )

    CALL ComputeJNorm &
           ( ITERATE, Jold_1, Jold_2, Jnorm_1, Jnorm_2 )

    ! --- Update Neutrino Densities ---

    CALL ComputeNumberDensity_EmAb_FP &
           ( ITERATE, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
             Chi_1, Chi_2, J0_1, J0_2 )

    k = 0
    DO WHILE( ANY( ITERATE(:) ) .AND. k < MaxIter )

      k  = k + 1
      Mk = MIN( M_FP, k )

      CALL CreatePackIndex &
             ( ITERATE, nX_P, PackIndex, UnpackIndex )

      IF ( k > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Im_ComputeOpacity )

        CALL ComputeEquilibriumDistributions_Packed &
               ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
                 ITERATE, nX_P, PackIndex, UnpackIndex )

        CALL TimersStop( Timer_Im_ComputeOpacity )

      END IF

      ! --- Right-Hand Side Vectors and Residuals ---

      CALL ComputeMatterRHS_FP &
             ( ITERATE, n_FP, iY, iE, FVECm, GVECm, Jnew_1, Jnew_2, &
               C_Y, S_Y, Unew_Y, GVEC_Y, C_E, S_E, Unew_E, GVEC_E )

      CALL ComputeNeutrinoRHS_EmAb_FP &
             ( ITERATE, n_FP, OS_1, OS_2, FVECm, GVECm, &
               dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
               Chi_1, Chi_2, J0_1, J0_2 )

      ! --- Anderson Acceleration ---

      CALL TimersStart( Timer_Im_ComputeLS )

      CALL SolveLS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVECm, GVECm, FVEC, GVEC, &
               AMAT, BVEC, Alpha )

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors---

      CALL TimersStart( Timer_Im_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE, n_FP, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm, GVECm )

      CALL UpdateNeutrinoRHS_FP &
             ( ITERATE, n_FP, OS_1, OS_2, &
               FVECm, GVECm, Jnew_1, Jnew_2 )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Check Convergence ---

      CALL CheckConvergenceCoupled &
             ( ITERATE, n_FP, k, iY, iE, OS_1, OS_2, Rtol, &
               nIterations, FVECm, Jnorm_1, Jnorm_2 )

      ! --- Shift History Arrays ---

      CALL ShiftRHS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVEC, GVEC )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    CALL ArrayCopy( Jnew_1, Jnew_2, J_1, J_2 )

    IF(PRESENT(nIterations_Out)) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iN_X = 1, nX_G
        nIterations_Out(iN_X) = nIterations(iN_X)
      END DO
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE, &
    !$OMP               Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP               Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP               Jold_1, Jnew_1, Jnorm_1, &
    !$OMP               Jold_2, Jnew_2, Jnorm_2, &
    !$OMP               PackIndex, UnpackIndex, &
    !$OMP               GVEC, FVEC, GVECm, FVECm, Alpha, nIterations )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE, &
    !$ACC         Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha, nIterations )
#endif

  END SUBROUTINE SolveMatterEquations_EmAb_FP


  SUBROUTINE SolveMatterEquations_FP_Coupled &
    ( dt, iS_1, iS_2, J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      D, T, Y, E, nIterations )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                           INTENT(in)    :: dt
    INTEGER,                            INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Eta_Pair_1, Eta_Pair_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(1:nX_G),        INTENT(out)   :: nIterations

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = iE
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-8
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E

    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_1, Jnew_1
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_2, Jnew_2
    REAL(DP), DIMENSION(       1:nX_G) :: Jnorm_1, Jnorm_2

    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2

    REAL(DP), DIMENSION(1:n_FP,1:M_FP,1:nX_G) :: GVEC, FVEC
    REAL(DP), DIMENSION(1:n_FP,       1:nX_G) :: GVECm, FVECm
    REAL(DP), DIMENSION(       1:M_FP,1:nX_G) :: Alpha

    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex, UnpackIndex

    INTEGER  :: k, Mk, nX_P
    INTEGER  :: OS_2

    OS_2 = OS_1 + nE_G
    ITERATE(:) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE ) &
    !$OMP MAP( alloc: Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP             Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP             Jold_1, Jnew_1, Jnorm_1, &
    !$OMP             Jold_2, Jnew_2, Jnorm_2, &
    !$OMP             Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP             Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP             PackIndex, UnpackIndex, &
    !$OMP             GVEC, FVEC, GVECm, FVECm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha )
#endif

    CALL ArrayCopy( Y, E, Yold, Eold )

    CALL ArrayCopy( J_1, J_2, Jold_1, Jold_2 )

    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J_1, J_2, &
               Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Initial RHS ---

    CALL InitializeRHS_FP &
           ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
             D, Y, Yold, C_Y, S_Y, Unew_Y, E, Eold, C_E, S_E, Unew_E )

    CALL ComputeJNorm &
           ( ITERATE, Jold_1, Jold_2, Jnorm_1, Jnorm_2 )

    ! --- Compute Neutrino Rates ---

    CALL TimersStart( Timer_Im_ComputeRate )

    CALL ComputeRates_Packed &
           ( Jnew_1, Jnew_2, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
             Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
             Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    CALL TimersStop( Timer_Im_ComputeRate )

    ! --- Update Neutrino Densities ---

    CALL ComputeNumberDensity_FP &
           ( ITERATE, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
             Chi_1, Chi_2, J0_1, J0_2, &
             Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
             Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    k = 0
    DO WHILE( ANY( ITERATE(:) ) .AND. k < MaxIter )

      k  = k + 1
      Mk = MIN( M_FP, k )

      CALL CreatePackIndex &
             ( ITERATE, nX_P, PackIndex, UnpackIndex )

      IF ( k > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Im_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 ITERATE, nX_P, PackIndex, UnpackIndex, nX_P )

        CALL TimersStop( Timer_Im_ComputeOpacity )

      END IF

      ! --- NES Emissivities and Opacities ---

      CALL TimersStart( Timer_Im_ComputeRate )

      CALL ComputeRates_Packed &
             ( Jnew_1, Jnew_2, &
               Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
               Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
               Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
               ITERATE, nX_P, PackIndex, UnpackIndex, nX_P )

      CALL TimersStop( Timer_Im_ComputeRate )

      ! --- Right-Hand Side Vectors and Residuals ---

      CALL ComputeMatterRHS_FP &
             ( ITERATE, n_FP, iY, iE, FVECm, GVECm, Jnew_1, Jnew_2, &
               C_Y, S_Y, Unew_Y, GVEC_Y, C_E, S_E, Unew_E, GVEC_E )

      CALL ComputeNeutrinoRHS_FP &
             ( ITERATE, n_FP, OS_1, OS_2, FVECm, GVECm, &
               dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
               Chi_1, Chi_2, J0_1, J0_2, &
               Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
               Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

      ! --- Anderson Acceleration ---

      CALL TimersStart( Timer_Im_ComputeLS )

      CALL SolveLS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVECm, GVECm, FVEC, GVEC, &
               AMAT, BVEC, Alpha )

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors---

      CALL TimersStart( Timer_Im_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE, n_FP, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm, GVECm )

      CALL UpdateNeutrinoRHS_FP &
             ( ITERATE, n_FP, OS_1, OS_2, &
               FVECm, GVECm, Jnew_1, Jnew_2 )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Check Convergence ---

      CALL CheckConvergenceCoupled &
             ( ITERATE, n_FP, k, iY, iE, OS_1, OS_2, Rtol, &
               nIterations, FVECm, Jnorm_1, Jnorm_2 )

      ! --- Shift History Arrays ---

      CALL ShiftRHS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVEC, GVEC )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    CALL ArrayCopy( Jnew_1, Jnew_2, J_1, J_2 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE, &
    !$OMP               Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP               Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP               Jold_1, Jnew_1, Jnorm_1, &
    !$OMP               Jold_2, Jnew_2, Jnorm_2, &
    !$OMP               Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP               Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP               Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP               PackIndex, UnpackIndex, &
    !$OMP               GVEC, FVEC, GVECm, FVECm, Alpha )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE, &
    !$ACC         Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         GVEC, FVEC, GVECm, FVECm, Alpha )
#endif

  END SUBROUTINE SolveMatterEquations_FP_Coupled


  SUBROUTINE SolveMatterEquations_FP_NestedAA &
    ( dt, iS_1, iS_2, J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      D, T, Y, E, nIterations_Inner, nIterations_Outer )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                           INTENT(in)    :: dt
    INTEGER,                            INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Eta_Pair_1, Eta_Pair_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(1:nX_G),        INTENT(out)   :: nIterations_Inner, nIterations_Outer

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = 0
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E

    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_1, Jnew_1
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_2, Jnew_2
    REAL(DP), DIMENSION(       1:nX_G) :: Jnorm_1, Jnorm_2

    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2

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
    INTEGER :: OS_2

    OS_2 = OS_1 + nE_G

    ITERATE_OUTER(:) = .TRUE.
    ITERATE_INNER(:) = .TRUE.

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

    CALL ArrayCopy( Y, E, Yold, Eold )

    CALL ArrayCopy( J_1, J_2, Jold_1, Jold_2 )

    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J_1, J_2, &
               Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Initial RHS ---

    CALL InitializeRHS_FP &
           ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
             D, Y, Yold, C_Y, S_Y, Unew_Y, E, Eold, C_E, S_E, Unew_E )

    k_outer = 0
    DO WHILE( ANY( ITERATE_OUTER(:) ) .AND. k_outer < MaxIter )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      CALL ComputeJNorm &
             ( ITERATE_OUTER, Jnew_1, Jnew_2, Jnorm_1, Jnorm_2 )

      CALL CreatePackIndex &
             ( ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      IF ( k_outer > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Im_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        CALL TimersStop( Timer_Im_ComputeOpacity )

      END IF

      ! START INNER LOOP
      CALL TimersStart( Timer_Im_NestInner )

      k_inner = 0
      DO WHILE( ANY( ITERATE_INNER(:) ) .AND. k_inner < MaxIter )

        k_inner = k_inner + 1
        Mk_inner = MIN( M_inner, k_inner )

        CALL CreatePackIndex &
               ( ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

        CALL TimersStart( Timer_Im_ComputeRate )

        CALL ComputeRates_Packed &
               ( Jnew_1, Jnew_2, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
                 Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
                 ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner, nX_P_outer )

        CALL TimersStop( Timer_Im_ComputeRate )

        ! --- Right-Hand Side Vectors and Residuals (inner) ---

        CALL ComputeNeutrinoRHS_FP &
               ( ITERATE_INNER, n_FP_inner, OS_1, OS_2, FVECm_inner, GVECm_inner, &
                 dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
                 Chi_1, Chi_2, J0_1, J0_2, &
                 Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
                 Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

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
               ( ITERATE_INNER, n_FP_inner, OS_1, OS_2, &
                 FVECm_inner, GVECm_inner, Jnew_1, Jnew_2 )

        ! --- Check Convergence (inner) ---

        CALL CheckConvergenceInner &
               ( ITERATE_INNER, n_FP_inner, k_inner, OS_1, OS_2, Rtol, &
                 nIterations_Inner, FVECm_inner, Jnorm_1, Jnorm_2 )

        ! --- Shift History Arrays (inner) ---

        CALL ShiftRHS_FP &
               ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, FVEC_inner, GVEC_inner )

        CALL TimersStop( Timer_Im_UpdateFP )

      END DO

      CALL TimersStop( Timer_Im_NestInner )

      ! --- Right-Hand Side Vectors and Residuals (outer) ---

      CALL ComputeMatterRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, iY, iE, FVECm_outer, GVECm_outer, &
               Jnew_1, Jnew_2, C_Y, S_Y, Unew_Y, GVEC_Y, C_E, S_E, Unew_E, GVEC_E )

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
             ( ITERATE_OUTER, n_FP_outer, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm_outer, GVECm_outer )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      ! --- Check Convergence (outer) ---

      CALL CheckConvergenceOuter &
             ( ITERATE_OUTER, ITERATE_INNER, n_FP_outer, iY, iE, k_outer, &
               FVECm_outer, Rtol, nIterations_Outer )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    nIterations_Inner &
      = FLOOR( DBLE( nIterations_Inner ) / DBLE( nIterations_Outer ) )

    CALL ArrayCopy( Jnew_1, Jnew_2, J_1, J_2 )

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


  SUBROUTINE SolveMatterEquations_FP_NestedNewton &
    ( dt, iS_1, iS_2, J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      D, T, Y, E, nIterations_Inner, nIterations_Outer )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                           INTENT(in)    :: dt
    INTEGER,                            INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Eta_Pair_1, Eta_Pair_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(1:nX_G),        INTENT(out)   :: nIterations_Inner, nIterations_Outer

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = 0
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E

    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_1, Jnew_1, S_1, SJ_1
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_2, Jnew_2, S_2, SJ_2
    REAL(DP), DIMENSION(       1:nX_G) :: Jnorm_1, Jnorm_2

    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2

    REAL(DP), DIMENSION(1:n_FP_outer,1:M_outer,1:nX_G) :: AMAT_outer, GVEC_outer, FVEC_outer
    REAL(DP), DIMENSION(1:n_FP_outer,          1:nX_G) :: GVECm_outer, FVECm_outer
    REAL(DP), DIMENSION(1:n_FP_outer,          1:nX_G) :: BVEC_outer
    REAL(DP), DIMENSION(             1:M_outer,1:nX_G) :: Alpha_outer

    REAL(DP), DIMENSION(1:n_FP_inner,1:n_FP_inner,1:nX_G) :: GJAC_inner
    REAL(DP), DIMENSION(1:n_FP_inner,             1:nX_G) :: GVECm_inner, FVECm_inner

    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE_OUTER, ITERATE_INNER
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex_outer, UnpackIndex_outer
    INTEGER,  DIMENSION(1:nX_G) :: PackIndex_inner, UnpackIndex_inner

    INTEGER :: k_outer, Mk_outer, nX_P_outer
    INTEGER :: k_inner, Mk_inner, nX_P_inner
    INTEGER :: OS_2

    OS_2 = OS_1 + nE_G

    ITERATE_OUTER(:) = .TRUE.
    ITERATE_INNER(:) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE_OUTER, ITERATE_INNER ) &
    !$OMP MAP( alloc: Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP             Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP             Jold_1, Jnew_1, Jnorm_1, S_1, SJ_1, &
    !$OMP             Jold_2, Jnew_2, Jnorm_2, S_2, SJ_2, &
    !$OMP             Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP             Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP             Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP             PackIndex_outer, UnpackIndex_outer, &
    !$OMP             PackIndex_inner, UnpackIndex_inner, &
    !$OMP             AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$OMP             GVECm_outer, FVECm_outer, Alpha_outer, &
    !$OMP             GJAC_inner, GVECm_inner, FVECm_inner )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE_OUTER, ITERATE_INNER ) &
    !$ACC CREATE( Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, S_1, SJ_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, S_2, SJ_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$ACC         GVECm_outer, FVECm_outer, Alpha_outer, &
    !$ACC         GJAC_inner, GVECm_inner, FVECm_inner )
#endif

    CALL ArrayCopy( Y, E, Yold, Eold )

    CALL ArrayCopy( J_1, J_2, Jold_1, Jold_2 )

    IF( UsePreconditionerEmAb )THEN

      CALL TimersStart( Timer_Im_EmAb_FP )

      CALL SolveMatterEquations_EmAb_FP &
             ( dt, iS_1, iS_2, J_1, J_2, &
               Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

      CALL TimersStop( Timer_Im_EmAb_FP )

    END IF

    ! --- Compute Opacity Kernels ---

    CALL TimersStart( Timer_Im_ComputeOpacity )

    CALL ComputeOpacities_Packed &
           ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
             Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
             Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

    CALL TimersStop( Timer_Im_ComputeOpacity )

    ! --- Initial RHS ---

    CALL InitializeRHS_Newton &
           ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
             dt, Chi_1, Chi_2, S_1, S_2, &
             D, Y, Yold, C_Y, S_Y, Unew_Y, E, Eold, C_E, S_E, Unew_E )

    k_outer = 0
    DO WHILE( ANY( ITERATE_OUTER(:) ) .AND. k_outer < MaxIter )

      k_outer  = k_outer + 1
      Mk_outer = MIN( M_outer, k_outer )

      CALL ComputeJNorm &
             ( ITERATE_OUTER, Jnew_1, Jnew_2, Jnorm_1, Jnorm_2 )

      CALL CreatePackIndex &
             ( ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      IF ( k_outer > 1 ) THEN

        ! --- Recompute Opacity Kernels ---

        CALL TimersStart( Timer_Im_ComputeOpacity )

        CALL ComputeOpacities_Packed &
               ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

        CALL TimersStop( Timer_Im_ComputeOpacity )

      END IF

      ! START INNER LOOP
      CALL TimersStart( Timer_Im_NestInner )

      k_inner = 0
      DO WHILE( ANY( ITERATE_INNER(:) ) .AND. k_inner < MaxIter )

        k_inner = k_inner + 1

        CALL CreatePackIndex &
               ( ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

        ! --- Compute Neutrino Rates ---

        CALL TimersStart( Timer_Im_ComputeRate )

        CALL ComputeRates_Packed &
               ( Jnew_1, Jnew_2, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
                 Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
                 ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner, nX_P_outer )

        CALL TimersStop( Timer_Im_ComputeRate )

        CALL ComputeNeutrinoRHS_Newton &
               ( ITERATE_INNER, nX_P_inner, UnpackIndex_inner, n_FP_inner, OS_1, OS_2, &
                 dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
                 S_1, S_2, Chi_1, Chi_2, J0_1, J0_2, &
                 Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
                 Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
                 GVECm_inner )

        ! --- Build and Solve Linear System (inner) ---

        CALL TimersStart( Timer_Im_ComputeLS )

        CALL BuildJacobian_Newton &
               ( ITERATE_INNER, nX_P_inner, UnpackIndex_inner, n_FP_inner, OS_1, OS_2, &
                 dt, Jnew_1, Jnew_2, S_1, S_2, &
                 Chi_NES_1, Chi_NES_2, Chi_Pair_1, Chi_Pair_2, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 GJAC_inner )

        CALL SolveLS_Newton &
               ( ITERATE_INNER, nX_P_inner, PackIndex_inner, n_FP_inner, OS_1, OS_2, &
                 Jnew_1, Jnew_2, FVECm_inner, GVECm_inner, GJAC_inner )

        CALL TimersStop( Timer_Im_ComputeLS )

        ! --- Update Residuals and Solution Vectors (inner) ---

        CALL TimersStart( Timer_Im_UpdateFP )

        ! --- Check Convergence (inner) ---

        CALL CheckConvergenceInner &
               ( ITERATE_INNER, n_FP_inner, k_inner, OS_1, OS_2, Rtol, &
                 nIterations_Inner, FVECm_inner, Jnorm_1, Jnorm_2 )

        CALL TimersStop( Timer_Im_UpdateFP )

      END DO

      CALL TimersStop( Timer_Im_NestInner )

      ! --- Right-Hand Side Vectors and Residuals (outer) ---

      CALL ComputeMatterRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, iY, iE, FVECm_outer, GVECm_outer, &
               Jnew_1, Jnew_2, C_Y, S_Y, Unew_Y, GVEC_Y, C_E, S_E, Unew_E, GVEC_E )

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
             ( ITERATE_OUTER, n_FP_outer, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm_outer, GVECm_outer )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      ! --- Check Convergence (outer) ---

      CALL CheckConvergenceOuter &
             ( ITERATE_OUTER, ITERATE_INNER, n_FP_outer, iY, iE, k_outer, &
               FVECm_outer, Rtol, nIterations_Outer )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    CALL ArrayCopy( Jnew_1, Jnew_2, J_1, J_2 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE_OUTER, ITERATE_INNER, &
    !$OMP               Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$OMP               Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$OMP               Jold_1, Jnew_1, Jnorm_1, S_1, SJ_1, &
    !$OMP               Jold_2, Jnew_2, Jnorm_2, S_2, SJ_2, &
    !$OMP               Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$OMP               Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$OMP               Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$OMP               Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$OMP               PackIndex_outer, UnpackIndex_outer, &
    !$OMP               PackIndex_inner, UnpackIndex_inner, &
    !$OMP               AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$OMP               GVECm_outer, FVECm_outer, Alpha_outer, &
    !$OMP               GJAC_inner, GVECm_inner, FVECm_inner )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE_OUTER, ITERATE_INNER, &
    !$ACC         Yold, S_Y, C_Y, Unew_Y, GVEC_Y, &
    !$ACC         Eold, S_E, C_E, Unew_E, GVEC_E, &
    !$ACC         Jold_1, Jnew_1, Jnorm_1, S_1, SJ_1, &
    !$ACC         Jold_2, Jnew_2, Jnorm_2, S_2, SJ_2, &
    !$ACC         Phi_0_In_NES_1, Phi_0_Ot_NES_1, &
    !$ACC         Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
    !$ACC         Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, &
    !$ACC         Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
    !$ACC         PackIndex_outer, UnpackIndex_outer, &
    !$ACC         PackIndex_inner, UnpackIndex_inner, &
    !$ACC         AMAT_outer, BVEC_outer, GVEC_outer, FVEC_outer, &
    !$ACC         GVECm_outer, FVECm_outer, Alpha_outer, &
    !$ACC         GJAC_inner, GVECm_inner, FVECm_inner )
#endif

  END SUBROUTINE SolveMatterEquations_FP_NestedNewton


  SUBROUTINE ComputeEquilibriumDistributionsAndDerivatives_Packed &
    ( iSpecies, D, T, Y, J0, dJ0dY, dJ0dE, &
      MASK, nX_P, PackIndex, UnpackIndex )

    INTEGER,                                      INTENT(in)    :: iSpecies
    REAL(DP), DIMENSION(1:nX_G),        TARGET,   INTENT(in)    :: D, T, Y
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), TARGET,   INTENT(inout) :: J0, dJ0dY, dJ0dE

    LOGICAL,  DIMENSION(1:nX_G),        OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                            OPTIONAL, INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(1:nX_G),        OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:),   POINTER :: D_P, T_P, Y_P
    REAL(DP), DIMENSION(:,:), POINTER :: J0_P, dJ0dY_P, dJ0dE_P

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

      J0_P    => P2D(:,1:nX,iP2D_J0_1)
      dJ0dY_P => P2D(:,1:nX,iP2D_J0_2)
      dJ0dE_P => P2D(:,1:nX,iP2D_S_1)

    ELSE

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      J0_P    => J0(:,:)
      dJ0dY_P => dJ0dY(:,:)
      dJ0dE_P => dJ0dE(:,:)

    END IF

    ! --- Equilibrium Distributions ---

    !CALL ComputeEquilibriumDistributionAndDerivatives_Points &
    !       ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_P, dJ0dY_P, dJ0dE_P, iSpecies )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               J0_P, dJ0dY_P, dJ0dE_P, J0, dJ0dY, dJ0dE )

    END IF

  END SUBROUTINE ComputeEquilibriumDistributionsAndDerivatives_Packed


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

    !CALL ComputeEquilibriumDistributions_Points &
    !       ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, J0_2_P, iS_1, iS_2 )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               J0_1_P, J0_2_P, J0_1, J0_2 )

    END IF

  END SUBROUTINE ComputeEquilibriumDistributions_Packed


  SUBROUTINE ComputeOpacities_Packed &
    ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
      Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
      Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
      MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    INTEGER,                                             INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(1:nX_G),               TARGET,   INTENT(in)    :: D, T, Y
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET,   INTENT(inout) :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(inout) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(inout) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(inout) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(inout) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2

    LOGICAL,  DIMENSION(1:nX_G),               OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                                   OPTIONAL, INTENT(in)    :: nX_P, nX_P0
    INTEGER,  DIMENSION(1:nX_G),               OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:),     POINTER :: D_P, T_P, Y_P
    REAL(DP), DIMENSION(:,:),   POINTER :: J0_1_P, J0_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P

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

      J0_1_P => P2D(:,1:nX,iP2D_J0_1)
      J0_2_P => P2D(:,1:nX,iP2D_J0_2)

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

      J0_1_P => J0_1(:,:)
      J0_2_P => J0_2(:,:)

      Phi_0_In_NES_1_P  => Phi_0_In_NES_1 (:,:,:)
      Phi_0_Ot_NES_1_P  => Phi_0_Ot_NES_1 (:,:,:)
      Phi_0_In_NES_2_P  => Phi_0_In_NES_2 (:,:,:)
      Phi_0_Ot_NES_2_P  => Phi_0_Ot_NES_2 (:,:,:)
      Phi_0_In_Pair_1_P => Phi_0_In_Pair_1(:,:,:)
      Phi_0_Ot_Pair_1_P => Phi_0_Ot_Pair_1(:,:,:)
      Phi_0_In_Pair_2_P => Phi_0_In_Pair_2(:,:,:)
      Phi_0_Ot_Pair_2_P => Phi_0_Ot_Pair_2(:,:,:)

    END IF

    ! --- Equilibrium Distributions ---

    !CALL ComputeEquilibriumDistributions_Points &
    !       ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, J0_2_P, iS_1, iS_2 )

    ! --- NES Kernels ---

    !CALL ComputeNeutrinoOpacities_NES_Points &
    !       ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, iS_2, 1, &
    !         Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, &
    !         Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
    !         P3D(:,:,:,iP3D_WORK1), P3D(:,:,:,iP3D_WORK2) )

    ! --- Pair Kernels ---

    !CALL ComputeNeutrinoOpacities_Pair_Points &
    !       ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, iS_2, 1, &
    !         Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
    !         Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
    !         P3D(:,:,:,iP3D_WORK1), P3D(:,:,:,iP3D_WORK2) )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               J0_1_P, J0_2_P, J0_1, J0_2 )

      IF ( nX < nX0 ) THEN

        CALL ArrayUnpack &
               ( nX, MASK, PackIndex, &
                 Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
                 Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

      END IF

    END IF

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( J_1, J_2, &
      Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
      Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET,   INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(in)    :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(in)    :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(in)    :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET,   INTENT(in)    :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET,   INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET,   INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET,   INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET,   INTENT(inout) :: Eta_Pair_1, Eta_Pair_2

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
               J_1, J_2, J_1_P, J_2_P )

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
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
                 Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P )

      END IF

    ELSE

      Chi_NES_1_P  => Chi_NES_1 (:,:)
      Chi_NES_2_P  => Chi_NES_2 (:,:)
      Eta_NES_1_P  => Eta_NES_1 (:,:)
      Eta_NES_2_P  => Eta_NES_2 (:,:)
      Chi_Pair_1_P => Chi_Pair_1(:,:)
      Chi_Pair_2_P => Chi_Pair_2(:,:)
      Eta_Pair_1_P => Eta_Pair_1(:,:)
      Eta_Pair_2_P => Eta_Pair_2(:,:)

      J_1_P => J_1(:,:)
      J_2_P => J_2(:,:)

      Phi_0_In_NES_1_P  => Phi_0_In_NES_1 (:,:,:)
      Phi_0_Ot_NES_1_P  => Phi_0_Ot_NES_1 (:,:,:)
      Phi_0_In_NES_2_P  => Phi_0_In_NES_2 (:,:,:)
      Phi_0_Ot_NES_2_P  => Phi_0_Ot_NES_2 (:,:,:)
      Phi_0_In_Pair_1_P => Phi_0_In_Pair_1(:,:,:)
      Phi_0_Ot_Pair_1_P => Phi_0_Ot_Pair_1(:,:,:)
      Phi_0_In_Pair_2_P => Phi_0_In_Pair_2(:,:,:)
      Phi_0_Ot_Pair_2_P => Phi_0_Ot_Pair_2(:,:,:)

    END IF

    ! --- NES Emissivities and Opacities ---

    !CALL ComputeNeutrinoOpacitiesRates_NES_Points &
    !       ( 1, nE_G, 1, nX, W2_N, J_1_P, &
    !         Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P, &
    !         Eta_NES_1_P, Chi_NES_1_P )

    !CALL ComputeNeutrinoOpacitiesRates_NES_Points &
    !       ( 1, nE_G, 1, nX, W2_N, J_2_P, &
    !         Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P, &
    !         Eta_NES_2_P, Chi_NES_2_P )

    ! --- Pair Emissivities and Opacities ---

    !CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
    !       ( 1, nE_G, 1, nX, W2_N, J_2_P, &
    !         Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P, &
    !         Eta_Pair_1_P, Chi_Pair_1_P )

    !CALL ComputeNeutrinoOpacitiesRates_Pair_Points &
    !       ( 1, nE_G, 1, nX, W2_N, J_1_P, &
    !         Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P, &
    !         Eta_Pair_2_P, Chi_Pair_2_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX, MASK, PackIndex, &
               Chi_NES_1_P, Chi_NES_2_P, Eta_NES_1_P, Eta_NES_2_P, &
               Chi_Pair_1_P, Chi_Pair_2_P, Eta_Pair_1_P, Eta_Pair_2_P, &
               Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
               Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

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
    ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
      D, Y, Yold, C_Y, S_Y, U_Y, E, Eold, C_E, S_E, U_E )

    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)  :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)  :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: D
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: Y, Yold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: E, Eold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_E, S_E, U_E

    INTEGER  :: iN_E, iN_X

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_1, nE_G, W2_S, 1, Zero, C_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, Jold_2, nE_G, W2_S, 1,  One, C_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_1, nE_G, W3_S, 1, Zero, C_E, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, Jold_2, nE_G, W3_S, 1,  One, C_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G

      S_Y(iN_X) = One / ( D(iN_X) * Yold(iN_X) / AtomicMassUnit )
      S_E(iN_X) = One / ( D(iN_X) * Eold(iN_X) )

      C_Y(iN_X) = C_Y(iN_X) * S_Y(iN_X)
      C_E(iN_X) = C_E(iN_X) * S_E(iN_X)

      U_Y(iN_X) = Y(iN_X) / Yold(iN_X) ! --- Initial Guess
      U_E(iN_X) = E(iN_X) / Eold(iN_X) ! --- Initial Guess

    END DO

    CALL ArrayCopy( J_1, J_2, Jnew_1, Jnew_2 )

  END SUBROUTINE InitializeRHS_FP


  SUBROUTINE InitializeRHS_Newton &
    ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
      dt, Chi_1, Chi_2, S_1, S_2, &
      D, Y, Yold, C_Y, S_Y, U_Y, E, Eold, C_E, S_E, U_E )

    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)  :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)  :: J_1, J_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out) :: Jnew_1, Jnew_2
    REAL(DP),                           INTENT(in)  :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)  :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out) :: S_1, S_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: D
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: Y, Yold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: E, Eold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: C_E, S_E, U_E

    INTEGER  :: iN_E, iN_X

    CALL InitializeRHS_FP &
           ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
             D, Y, Yold, C_Y, S_Y, U_Y, E, Eold, C_E, S_E, U_E )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        S_1(iN_E,iN_X) = One / ( One + dt * Chi_1(iN_E,iN_X) )
        S_2(iN_E,iN_X) = One / ( One + dt * Chi_2(iN_E,iN_X) )

      END DO
    END DO

  END SUBROUTINE InitializeRHS_Newton


  SUBROUTINE InitializeRHS_EmAb_NuE &
    ( dt, J, Chi, J0, dJ0dY, dJ0dE, D, Y, E, &
      Gam, GamJ, GamJ0, dGamJ0dY, dGamJ0dE, &
      Yold, Eold, U_Y, U_E, C_Y, C_E, F_Y, F_E, FNRM0 )

    REAL(DP),                           INTENT(in)  :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)  :: J, Chi, J0, dJ0dY, dJ0dE
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)  :: D, Y, E
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(out) :: Gam, GamJ, GamJ0, dGamJ0dY, dGamJ0dE
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out) :: Yold, Eold, U_Y, U_E, C_Y, C_E, F_Y, F_E, FNRM0

    REAL(DP) :: N_B
    INTEGER  :: iN_E, iN_X

    ! --- Auxiliary Variables ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        Gam     (iN_E,iN_X) = dt * Chi(iN_E,iN_X) / ( One + dt * Chi(iN_E,iN_X) )
        GamJ    (iN_E,iN_X) = Gam(iN_E,iN_X) * J    (iN_E,iN_X)
        GamJ0   (iN_E,iN_X) = Gam(iN_E,iN_X) * J0   (iN_E,iN_X)
        dGamJ0dY(iN_E,iN_X) = Gam(iN_E,iN_X) * dJ0dY(iN_E,iN_X)
        dGamJ0dE(iN_E,iN_X) = Gam(iN_E,iN_X) * dJ0dE(iN_E,iN_X)

      END DO
    END DO

    ! --- Old States (Constant) ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, GamJ , nE_G, W2_S, 1, Zero, C_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, GamJ , nE_G, W3_S, 1, Zero, C_E, 1 )

    ! --- Electron Fraction Equation ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, GamJ0, nE_G, W2_S, 1, Zero, F_Y, 1 )

    ! --- Internal Energy Equation ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, GamJ0, nE_G, W3_S, 1, Zero, F_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( N_B )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( N_B )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( N_B )
#endif
    DO iN_X = 1, nX_G

      N_B = D(iN_X) / AtomicMassUnit

      ! --- Initial Guess ---

      Yold(iN_X) = Y(iN_X)
      Eold(iN_X) = E(iN_X)

      U_Y(iN_X) = Yold(iN_X)
      U_E(iN_X) = Eold(iN_X)

      ! --- Scale Equations and Save Initial Evaluation ---

      ! --- Electron Fraction Equation ---

      C_Y(iN_X) =   C_Y(iN_X) + N_B     * U_Y(iN_X)
      F_Y(iN_X) = ( F_Y(iN_X) + N_B     * U_Y(iN_X) - C_Y(iN_X) ) / C_Y(iN_X)

      ! --- Internal Energy Equation ---

      C_E(iN_X) =   C_E(iN_X) + D(iN_X) * U_E(iN_X)
      F_E(iN_X) = ( F_E(iN_X) + D(iN_X) * U_E(iN_X) - C_E(iN_X) ) / C_E(iN_X)

      FNRM0(iN_X) = SQRT( F_Y(iN_X)**2 + F_E(iN_X)**2 )

    END DO

  END SUBROUTINE InitializeRHS_EmAb_NuE


  SUBROUTINE ComputeJNorm &
    ( MASK, J_1, J_2, Jnorm_1, Jnorm_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: Jnorm_1, Jnorm_2

    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN
        Jnorm_1(iN_X) = SQRT( SUM( J_1(:,iN_X)**2 ) )
        Jnorm_2(iN_X) = SQRT( SUM( J_2(:,iN_X)**2 ) )
      END IF
    END DO

  END SUBROUTINE ComputeJNorm


  SUBROUTINE ComputeNumberDensity_FP &
    ( MASK, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Eta_Pair_1, Eta_Pair_2

    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
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

  END SUBROUTINE ComputeNumberDensity_FP


  SUBROUTINE ComputeNumberDensity_EmAb_FP &
    ( MASK, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J0_1, J0_2

    REAL(DP) :: Eta
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Eta = Chi_1(iN_E,iN_X) * J0_1(iN_E,iN_X)
          Jnew_1(iN_E,iN_X) = ( Jold_1(iN_E,iN_X) + dt * Eta ) / ( One + dt * Chi_1(iN_E,iN_X) )

          Eta = Chi_2(iN_E,iN_X) * J0_2(iN_E,iN_X)
          Jnew_2(iN_E,iN_X) = ( Jold_2(iN_E,iN_X) + dt * Eta ) / ( One + dt * Chi_2(iN_E,iN_X) )

        END IF
      END DO
    END DO

  END SUBROUTINE ComputeNumberDensity_EmAb_FP


  SUBROUTINE ComputeNumberDensity_EmAb_NuE &
    ( dt, J, Chi, J0 )

    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: J
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J0

    REAL(DP) :: Eta
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        Eta = Chi(iN_E,iN_X) * J0(iN_E,iN_X)
        J(iN_E,iN_X) = ( J(iN_E,iN_X) + dt * Eta ) / ( One + dt * Chi(iN_E,iN_X) )

      END DO
    END DO

  END SUBROUTINE ComputeNumberDensity_EmAb_NuE


  SUBROUTINE ComputeMatterRHS_FP &
    ( MASK, n_FP, iY, iE, Fm, Gm, J_1, J_2, C_Y, S_Y, U_Y, G_Y, C_E, S_E, U_E, G_E )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, iY, iE
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: C_E, S_E, U_E
    REAL(DP), DIMENSION(1:nX_G),        INTENT(out)   :: G_Y, G_E

    INTEGER  :: iN_X

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J_1, nE_G, W2_S, 1, Zero, G_Y, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, -One, J_2, nE_G, W2_S, 1,  One, G_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J_1, nE_G, W3_S, 1, Zero, G_E, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, +One, J_2, nE_G, W3_S, 1,  One, G_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        G_Y(iN_X) = One + C_Y(iN_X) - G_Y(iN_X) * S_Y(iN_X)
        G_E(iN_X) = One + C_E(iN_X) - G_E(iN_X) * S_E(iN_X)

        Gm(iY,iN_X) = G_Y(iN_X)
        Gm(iE,iN_X) = G_E(iN_X)

        Fm(iY,iN_X) = G_Y(iN_X) - U_Y(iN_X)
        Fm(iE,iN_X) = G_E(iN_X) - U_E(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_FP


  SUBROUTINE ComputeMatterRHS_EmAb_NuE &
    ( MASK, D, U_Y, U_E, C_Y, C_E, Gam, J0, dJ0dY, dJ0dE, &
      GamJ0, dGamJ0dY, dGamJ0dE, F_Y, F_E )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: D, U_Y, U_E, C_Y, C_E
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Gam, J0, dJ0dY, dJ0dE
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: GamJ0, dGamJ0dY, dGamJ0dE
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: F_Y, F_E

    REAL(DP) :: N_B
    INTEGER  :: iN_E, iN_X

    ! --- Auxiliary Variables ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          GamJ0   (iN_E,iN_X) = Gam(iN_E,iN_X) * J0   (iN_E,iN_X)
          dGamJ0dY(iN_E,iN_X) = Gam(iN_E,iN_X) * dJ0dY(iN_E,iN_X)
          dGamJ0dE(iN_E,iN_X) = Gam(iN_E,iN_X) * dJ0dE(iN_E,iN_X)

        END IF
      END DO
    END DO

    ! --- Electron Fraction Equation ---

    CALL MatrixVectorMultiply &
           ( 'T', nE_G, nX_G, One, GamJ0, nE_G, W2_S, 1, Zero, F_Y, 1 )

    ! --- Internal Energy Equation ---

    CALL MatrixVectorMultiply &
           ( 'T', nE_G, nX_G, One, GamJ0, nE_G, W3_S, 1, Zero, F_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( N_B )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( N_B )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( N_B )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        N_B = D(iN_X) / AtomicMassUnit

        F_Y(iN_X) = ( F_Y(iN_X) + N_B     * U_Y(iN_X) - C_Y(iN_X) ) / C_Y(iN_X)
        F_E(iN_X) = ( F_E(iN_X) + D(iN_X) * U_E(iN_X) - C_E(iN_X) ) / C_E(iN_X)

      END IF
    END DO

  END SUBROUTINE ComputeMatterRHS_EmAb_NuE


  SUBROUTINE ComputeNeutrinoRHS_FP &
    ( MASK, n_FP, OS_1, OS_2, Fm, Gm, &
      dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(inout) :: Fm, Gm
    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Eta_Pair_1, Eta_Pair_2

    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta, Eta_T, Chi_T )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
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

  END SUBROUTINE ComputeNeutrinoRHS_FP


  SUBROUTINE ComputeNeutrinoRHS_Newton &
    ( MASK, nX_P, UnpackIndex, n_FP, OS_1, OS_2, &
      dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      S_1, S_2, Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      G )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(1:nX_G),        INTENT(in)    :: UnpackIndex
    INTEGER,                            INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: S_1, S_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J0_1, J0_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Eta_Pair_1, Eta_Pair_2
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(inout) :: G

    REAL(DP) :: EtaT_1, EtaT_2
    INTEGER  :: iN_E, iN_X, iX_P

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iN_X, EtaT_1, EtaT_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iN_X, EtaT_1, EtaT_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_X, EtaT_1, EtaT_2 )
#endif
    DO iX_P = 1, nX_P
      DO iN_E = 1, nE_G

        iN_X = UnpackIndex(iX_P)

        EtaT_1 = Eta_NES_1(iN_E,iN_X) + Eta_Pair_1(iN_E,iN_X) + Chi_1(iN_E,iN_X) * J0_1(iN_E,iN_X)

        G(OS_1+iN_E,iX_P) &
          = ( One + dt * ( Chi_NES_1(iN_E,iN_X) + Chi_Pair_1(iN_E,iN_X) ) * S_1(iN_E,iN_X) ) &
              * Jnew_1(iN_E,iN_X) - ( Jold_1(iN_E,iN_X) + dt * EtaT_1 ) * S_1(iN_E,iN_X)

        EtaT_2 = Eta_NES_2(iN_E,iN_X) + Eta_Pair_2(iN_E,iN_X) + Chi_2(iN_E,iN_X) * J0_2(iN_E,iN_X)

        G(OS_2+iN_E,iX_P) &
          = ( One + dt * ( Chi_NES_2(iN_E,iN_X) + Chi_Pair_2(iN_E,iN_X) ) * S_2(iN_E,iN_X) ) &
              * Jnew_2(iN_E,iN_X) - ( Jold_2(iN_E,iN_X) + dt * EtaT_2 ) * S_2(iN_E,iN_X)

      END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_Newton


  SUBROUTINE ComputeNeutrinoRHS_EmAb_FP &
    ( MASK, n_FP, OS_1, OS_2, Fm, Gm, &
      dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(inout) :: Fm, Gm
    REAL(DP),                           INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(in)    :: J0_1, J0_2

    REAL(DP) :: Eta
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Eta = Chi_1(iN_E,iN_X) * J0_1(iN_E,iN_X)
          Gm(OS_1+iN_E,iN_X) = ( Jold_1(iN_E,iN_X) + dt * Eta ) / ( One + dt * Chi_1(iN_E,iN_X) )
          Fm(OS_1+iN_E,iN_X) = Gm(OS_1+iN_E,iN_X) - Jnew_1(iN_E,iN_X)

          Eta = Chi_2(iN_E,iN_X) * J0_2(iN_E,iN_X)
          Gm(OS_2+iN_E,iN_X) = ( Jold_2(iN_E,iN_X) + dt * Eta ) / ( One + dt * Chi_2(iN_E,iN_X) )
          Fm(OS_2+iN_E,iN_X) = Gm(OS_2+iN_E,iN_X) - Jnew_2(iN_E,iN_X)

        END IF
      END DO
    END DO

  END SUBROUTINE ComputeNeutrinoRHS_EmAb_FP


  SUBROUTINE UpdateMatterRHS_FP &
    ( MASK, n_FP, iY, iE, Yold, Eold, Y, E, U_Y, U_E, Fm, Gm )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, iY, iE
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: Yold, Eold
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: Y, E
    REAL(DP), DIMENSION(1:nX_G),        INTENT(inout) :: U_Y, U_E
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(inout) :: Fm, Gm

    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        Fm(iY,iN_X) = Gm(iY,iN_X) - U_Y(iN_X)
        Fm(iE,iN_X) = Gm(iE,iN_X) - U_E(iN_X)

        U_Y(iN_X) = Gm(iY,iN_X)
        U_E(iN_X) = Gm(iE,iN_X)

        Y(iN_X) = U_Y(iN_X) * Yold(iN_X)
        E(iN_X) = U_E(iN_X) * Eold(iN_X)

      END IF
    END DO

  END SUBROUTINE UpdateMatterRHS_FP


  SUBROUTINE UpdateNeutrinoRHS_FP &
    ( MASK, n_FP, OS_1, OS_2, Fm, Gm, Jnew_1, Jnew_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(in)    :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(1:nE_G,1:nX_G), INTENT(inout) :: Jnew_1, Jnew_2

    INTEGER  :: iN_E, iN_X, iFP

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Fm(OS_1+iN_E,iN_X) = Gm(OS_1+iN_E,iN_X) - Jnew_1(iN_E,iN_X)
          Fm(OS_2+iN_E,iN_X) = Gm(OS_2+iN_E,iN_X) - Jnew_2(iN_E,iN_X)

          Jnew_1(iN_E,iN_X)  = Gm(OS_1+iN_E,iN_X)
          Jnew_2(iN_E,iN_X)  = Gm(OS_2+iN_E,iN_X)

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
    !$OMP PARALLEL DO COLLAPSE(2)
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
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( AA11, AB1 )
#endif
        DO iN_X = 1, nX_G
          IF ( MASK(iN_X) ) THEN

            AA11 = Zero
            AB1 = Zero

#if defined(THORNADO_OMP_OL)
            !$OMP PARALLEL DO &
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
          !$OMP PARALLEL DO &
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

#if defined(THORNADO_OMP_OL)
              !$OMP PARALLEL DO &
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
      !$OMP PARALLEL DO &
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


  SUBROUTINE BuildJacobian_Newton &
    ( MASK, nX_P, UnpackIndex, n_FP, OS_1, OS_2, &
      dt, Jnew_1, Jnew_2, S_1, S_2, &
      Chi_NES_1, Chi_NES_2, Chi_Pair_1, Chi_Pair_2, &
      Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
      Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
      GJAC )

    LOGICAL,  DIMENSION(1:nX_G),                       INTENT(in)    :: MASK
    INTEGER,                                           INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(1:nX_G),                       INTENT(in)    :: UnpackIndex
    INTEGER,                                           INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP),                                          INTENT(in)    :: dt
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),                INTENT(in)    :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET, INTENT(in)    :: S_1, S_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET, INTENT(in)    :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        TARGET, INTENT(in)    :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET, INTENT(in)    :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET, INTENT(in)    :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET, INTENT(in)    :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(1:nE_G,1:nE_G,1:nX_G), TARGET, INTENT(in)    :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2
    REAL(DP), DIMENSION(1:n_FP,1:n_FP,1:nX_G),         INTENT(inout) :: GJAC

    REAL(DP), DIMENSION(:,:),   POINTER :: S_1_P, S_2_P, SJ_1_P, SJ_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Chi_NES_1_P, Chi_NES_2_P
    REAL(DP), DIMENSION(:,:),   POINTER :: Chi_Pair_1_P, Chi_Pair_2_P

    REAL(DP) :: GJACT(1:n_FP,1:n_FP)
    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_E1, iN_E2, iN_X, iX_P

    IF ( nX_P < nX_G ) THEN

      ! --- Pack Arrays ---

      Chi_NES_1_P  => P2D(:,1:nX_P,iP2D_Chi_NES_1)
      Chi_NES_2_P  => P2D(:,1:nX_P,iP2D_Chi_NES_2)
      Chi_Pair_1_P => P2D(:,1:nX_P,iP2D_Chi_Pair_1)
      Chi_Pair_2_P => P2D(:,1:nX_P,iP2D_Chi_Pair_2)

      S_1_P => P2D(:,1:nX_P,iP2D_S_1)
      S_2_P => P2D(:,1:nX_P,iP2D_S_2)

      CALL ArrayPack &
             ( nX_P, UnpackIndex, &
               S_1, S_2, S_1_P, S_2_P )

      Phi_0_In_NES_1_P  => P3D(:,:,1:nX_P,iP3D_Phi_0_In_NES_1)
      Phi_0_Ot_NES_1_P  => P3D(:,:,1:nX_P,iP3D_Phi_0_Ot_NES_1)
      Phi_0_In_NES_2_P  => P3D(:,:,1:nX_P,iP3D_Phi_0_In_NES_2)
      Phi_0_Ot_NES_2_P  => P3D(:,:,1:nX_P,iP3D_Phi_0_Ot_NES_2)
      Phi_0_In_Pair_1_P => P3D(:,:,1:nX_P,iP3D_Phi_0_In_Pair_1)
      Phi_0_Ot_Pair_1_P => P3D(:,:,1:nX_P,iP3D_Phi_0_Ot_Pair_1)
      Phi_0_In_Pair_2_P => P3D(:,:,1:nX_P,iP3D_Phi_0_In_Pair_2)
      Phi_0_Ot_Pair_2_P => P3D(:,:,1:nX_P,iP3D_Phi_0_Ot_Pair_2)

    ELSE

      Chi_NES_1_P  => Chi_NES_1 (:,:)
      Chi_NES_2_P  => Chi_NES_2 (:,:)
      Chi_Pair_1_P => Chi_Pair_1(:,:)
      Chi_Pair_2_P => Chi_Pair_2(:,:)

      S_1_P => S_1 (:,:)
      S_2_P => S_2 (:,:)

      Phi_0_In_NES_1_P  => Phi_0_In_NES_1 (:,:,:)
      Phi_0_Ot_NES_1_P  => Phi_0_Ot_NES_1 (:,:,:)
      Phi_0_In_NES_2_P  => Phi_0_In_NES_2 (:,:,:)
      Phi_0_Ot_NES_2_P  => Phi_0_Ot_NES_2 (:,:,:)
      Phi_0_In_Pair_1_P => Phi_0_In_Pair_1(:,:,:)
      Phi_0_Ot_Pair_1_P => Phi_0_Ot_Pair_1(:,:,:)
      Phi_0_In_Pair_2_P => Phi_0_In_Pair_2(:,:,:)
      Phi_0_Ot_Pair_2_P => Phi_0_Ot_Pair_2(:,:,:)

    END IF

    SJ_1_P => P2D(:,1:nX_P,iP2D_J_1)
    SJ_2_P => P2D(:,1:nX_P,iP2D_J_2)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iN_X )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iN_X )
#endif
    DO iX_P = 1, nX_P
      DO iN_E = 1, nE_G

        iN_X = UnpackIndex(iX_P)

        SJ_1_P(iN_E,iX_P) = S_1_P(iN_E,iX_P) * Jnew_1(iN_E,iN_X)
        SJ_2_P(iN_E,iX_P) = S_2_P(iN_E,iX_P) * Jnew_2(iN_E,iN_X)

      END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nE_G, nE_G*nX_P, 1, One, W2_N, nE_G, &
             SJ_1_P, 1, Zero, ATMP_1, nE_G )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nE_G, nE_G*nX_P, 1, One, W2_N, nE_G, &
             S_1_P, 1, Zero, BTMP_1, nE_G )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nE_G, nE_G*nX_P, 1, One, W2_N, nE_G, &
             SJ_2_P, 1, Zero, ATMP_2, nE_G )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nE_G, nE_G*nX_P, 1, One, W2_N, nE_G, &
             S_2_P, 1, Zero, BTMP_2, nE_G )

    ! --- Jacobian matrix ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( GJACT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE( GJACT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( GJACT )
#endif
    DO iX_P = 1, nX_P

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iN_E2 = 1, nE_G
        DO iN_E1 = 1, nE_G

          GJACT(OS_1+iN_E1,OS_1+iN_E2) &
            = - dt * (   ATMP_1(iN_E1,iN_E2,iX_P) &
                         * (   Phi_0_Ot_NES_1_P(iN_E1,iN_E2,iX_P) &
                             - Phi_0_In_NES_1_P(iN_E1,iN_E2,iX_P) ) &
                       + BTMP_1(iN_E1,iN_E2,iX_P) &
                         *     Phi_0_In_NES_1_P(iN_E1,iN_E2,iX_P) )

          GJACT(OS_2+iN_E1,OS_1+iN_E2) &
            = + dt * (   ATMP_1(iN_E1,iN_E2,iX_P) &
                         * (   Phi_0_Ot_Pair_1_P(iN_E1,iN_E2,iX_P) &
                             - Phi_0_In_Pair_1_P(iN_E1,iN_E2,iX_P) ) &
                       + BTMP_1(iN_E1,iN_E2,iX_P) &
                         *     Phi_0_In_Pair_1_P(iN_E1,iN_E2,iX_P) )

        END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iN_E2 = 1, nE_G
        DO iN_E1 = 1, nE_G

          GJACT(OS_1+iN_E1,OS_2+iN_E2) &
            = + dt * (   ATMP_2(iN_E1,iN_E2,iX_P) &
                         * (   Phi_0_Ot_Pair_2_P(iN_E1,iN_E2,iX_P) &
                             - Phi_0_In_Pair_2_P(iN_E1,iN_E2,iX_P) ) &
                       + BTMP_2(iN_E1,iN_E2,iX_P) &
                         *     Phi_0_In_Pair_2_P(iN_E1,iN_E2,iX_P) )

          GJACT(OS_2+iN_E1,OS_2+iN_E2) &
            = - dt * (   ATMP_2(iN_E1,iN_E2,iX_P) &
                         * (   Phi_0_Ot_NES_2_P(iN_E1,iN_E2,iX_P) &
                             - Phi_0_In_NES_2_P(iN_E1,iN_E2,iX_P) ) &
                       + BTMP_2(iN_E1,iN_E2,iX_P) &
                         *     Phi_0_In_NES_2_P(iN_E1,iN_E2,iX_P) )

        END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iN_E2 = 1, n_FP
        DO iN_E1 = 1, n_FP
          GJAC(iN_E1,iN_E2,iX_P) = GJACT(iN_E2,iN_E1)
        END DO
      END DO

    END DO

    ! --- diagonal terms ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX_P = 1, nX_P
      DO iN_E = 1, nE_G

        GJAC(OS_1+iN_E,OS_1+iN_E,iX_P) &
          = GJAC(OS_1+iN_E,OS_1+iN_E,iX_P) &
            + One + dt * ( Chi_NES_1_P(iN_E,iX_P) + Chi_Pair_1_P(iN_E,iX_P) ) * S_1_P(iN_E,iX_P)

        GJAC(OS_2+iN_E,OS_2+iN_E,iX_P) &
          = GJAC(OS_2+iN_E,OS_2+iN_E,iX_P) &
            + One + dt * ( Chi_NES_2_P(iN_E,iX_P) + Chi_Pair_2_P(iN_E,iX_P) ) * S_2_P(iN_E,iX_P)

      END DO
    END DO

  END SUBROUTINE BuildJacobian_Newton


  SUBROUTINE SolveLS_Newton &
    ( MASK, nX_P, PackIndex, n_FP, OS_1, OS_2, &
      Jnew_1, Jnew_2, F, G, GJAC )

    LOGICAL,  DIMENSION(1:nX_G),               INTENT(in)    :: MASK
    INTEGER,                                   INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(1:nX_G),               INTENT(in)    :: PackIndex
    INTEGER,                                   INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),        INTENT(inout) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(1:n_FP,1:nX_G),        INTENT(inout) :: F, G
    REAL(DP), DIMENSION(1:n_FP,1:n_FP,1:nX_G), INTENT(inout) :: GJAC

    INTEGER  :: iN_E, iN_X

    CALL LinearSolveBatched &
           ( 'N', n_FP, 1, GJAC, n_FP, IPIV, G, n_FP, INFO, nX_P )

    CALL ArrayUnpack &
           ( nX_P, MASK, PackIndex, G, F )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        IF ( MASK(iN_X) ) THEN

          Jnew_1(iN_E,iN_X) = Jnew_1(iN_E,iN_X) - F(OS_1+iN_E,iN_X)
          Jnew_2(iN_E,iN_X) = Jnew_2(iN_E,iN_X) - F(OS_2+iN_E,iN_X)

        END IF
      END DO
    END DO

  END SUBROUTINE SolveLS_Newton


  SUBROUTINE SolveLS_EmAb_NuE &
    ( MASK, D, C_Y, C_E, F_Y, F_E, dGamJ0dY, dGamJ0dE, &
      FJAC, dU_Y, dU_E, U_Y, U_E, Y, E )

    LOGICAL,  DIMENSION(1:nX_G),         INTENT(in)    :: MASK
    REAL(DP), DIMENSION(1:nX_G),         INTENT(in)    :: D, C_Y, C_E, F_Y, F_E
    REAL(DP), DIMENSION(1:nE_G,1:nX_G),  INTENT(in)    :: dGamJ0dY, dGamJ0dE
    REAL(DP), DIMENSION(1:2,1:2,1:nX_G), INTENT(inout) :: FJAC
    REAL(DP), DIMENSION(1:nX_G),         INTENT(inout) :: dU_Y, dU_E, U_Y, U_E, Y, E

    REAL(DP) :: N_B, DJAC
    INTEGER  :: iN_X

    ! --- Build Jacobian ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, dGamJ0dY, nE_G, W2_S, 1, Zero, FJAC(1,1,1), 4 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, dGamJ0dY, nE_G, W3_S, 1, Zero, FJAC(2,1,1), 4 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, dGamJ0dE, nE_G, W2_S, 1, Zero, FJAC(1,2,1), 4 )
    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, dGamJ0dE, nE_G, W3_S, 1, Zero, FJAC(2,2,1), 4 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( DJAC, N_B )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( DJAC, N_B )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( DJAC, N_B )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        N_B = D(iN_X) / AtomicMassUnit

        ! --- Scale Jacobian ---

        FJAC(1,1,iN_X) = ( FJAC(1,1,iN_X) + N_B ) / C_Y(iN_X)
        FJAC(1,2,iN_X) = FJAC(1,2,iN_X) / C_Y(iN_X)

        FJAC(2,1,iN_X) = FJAC(2,1,iN_X) / C_E(iN_X)
        FJAC(2,2,iN_X) = ( FJAC(2,2,iN_X) + D  (iN_X) ) / C_E(iN_X)

        ! --- Determinant of Jacobian ---

        DJAC = FJAC(1,1,iN_X) * FJAC(2,2,iN_X) - FJAC(2,1,iN_X) * FJAC(1,2,iN_X)

        ! --- Correction ---

        dU_Y(iN_X) = - ( + FJAC(2,2,iN_X) * F_Y(iN_X) - FJAC(1,2,iN_X) * F_E(iN_X) ) / DJAC
        dU_E(iN_X) = - ( - FJAC(2,1,iN_X) * F_Y(iN_X) + FJAC(1,1,iN_X) * F_E(iN_X) ) / DJAC

        ! --- Apply Correction ---

        U_Y(iN_X) = U_Y(iN_X) + dU_Y(iN_X)
        U_E(iN_X) = U_E(iN_X) + dU_E(iN_X)

        Y(iN_X) = U_Y(iN_X)
        E(iN_X) = U_E(iN_X)

      END IF
    END DO

  END SUBROUTINE SolveLS_EmAb_NuE


  SUBROUTINE ShiftRHS_FP &
    ( MASK, n_FP, M, Mk, F, G )

    LOGICAL,  DIMENSION(1:nX_G),            INTENT(in)    :: MASK
    INTEGER,                                INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(1:n_FP,1:M,1:nX_G), INTENT(inout) :: F, G

    REAL(DP) :: FTMP(1:n_FP,1:M), GTMP(1:n_FP,1:M)
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
          !$OMP PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM = 1, Mk-1
            DO iFP = 1, n_FP
              FTMP(iFP,iM) = F(iFP,iM+1,iN_X)
              GTMP(iFP,iM) = G(iFP,iM+1,iN_X)
            END DO
          END DO

#if defined(THORNADO_OMP_OL)
          !$OMP PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
          !$ACC LOOP VECTOR COLLAPSE(2)
#endif
          DO iM = 1, Mk-1
            DO iFP = 1, n_FP
              F(iFP,iM,iN_X) = FTMP(iFP,iM)
              G(iFP,iM,iN_X) = GTMP(iFP,iM)
            END DO
          END DO

        END IF
      END DO

    END IF

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE CheckConvergenceInner &
    ( MASK, n_FP, k_inner, OS_1, OS_2, Rtol, nIterations_Inner, Fm, Jnorm_1, Jnorm_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(inout) :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, k_inner, OS_1, OS_2
    REAL(DP),                           INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(1:nX_G),        INTENT(inout) :: nIterations_Inner
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(in)    :: Fm
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: Jnorm_1, Jnorm_2

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
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, Fnorm_1, Fnorm_2 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        Fnorm_1 = SQRT( SUM( Fm(OS_1+1:OS_1+nE_G,iN_X)**2 ) )
        Fnorm_2 = SQRT( SUM( Fm(OS_2+1:OS_2+nE_G,iN_X)**2 ) )

        CONVERGED = Fnorm_1 <= Rtol * Jnorm_1(iN_X) .AND. &
                    Fnorm_2 <= Rtol * Jnorm_2(iN_X)

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
    ( MASK_OUTER, MASK_INNER, n_FP, iY, iE, k_outer, Fm, Rtol, nIterations_Outer )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(inout) :: MASK_OUTER, MASK_INNER
    INTEGER,                            INTENT(in)    :: n_FP, iY, iE, k_outer
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(in)    :: Fm
    REAL(DP),                           INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(1:nX_G),        INTENT(inout) :: nIterations_Outer

    LOGICAL  :: CONVERGED
    REAL(DP) :: Fnorm_Y, Fnorm_E
    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK_OUTER(iN_X) ) THEN

        Fnorm_Y = ABS( Fm(iY,iN_X) )
        Fnorm_E = ABS( Fm(iE,iN_X) )

        CONVERGED = Fnorm_Y <= Rtol .AND. &
                    Fnorm_E <= Rtol

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


  SUBROUTINE CheckConvergenceCoupled &
    ( MASK, n_FP, k, iY, iE, OS_1, OS_2, Rtol, nIterations, Fm, Jnorm_1, Jnorm_2 )

    LOGICAL,  DIMENSION(1:nX_G),        INTENT(inout) :: MASK
    INTEGER,                            INTENT(in)    :: n_FP, k, iY, iE, OS_1, OS_2
    REAL(DP),                           INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(1:nX_G),        INTENT(inout) :: nIterations
    REAL(DP), DIMENSION(1:n_FP,1:nX_G), INTENT(in)    :: Fm
    REAL(DP), DIMENSION(1:nX_G),        INTENT(in)    :: Jnorm_1, Jnorm_2

    REAL(DP) :: Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2
    LOGICAL  :: CONVERGED
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2 )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        Fnorm_Y = ABS( Fm(iY,iN_X) )
        Fnorm_E = ABS( Fm(iE,iN_X) )

        Fnorm_1 = SQRT( SUM( Fm(OS_1+1:OS_1+nE_G,iN_X)**2 ) )
        Fnorm_2 = SQRT( SUM( Fm(OS_2+1:OS_2+nE_G,iN_X)**2 ) )

        CONVERGED = Fnorm_Y <= Rtol .AND. &
                    Fnorm_E <= Rtol .AND. &
                    Fnorm_1 <= Rtol * Jnorm_1(iN_X) .AND. &
                    Fnorm_2 <= Rtol * Jnorm_2(iN_X)

        IF ( CONVERGED ) THEN
          MASK(iN_X) = .FALSE.
          nIterations(iN_X) = k
        END IF

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( MASK )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( MASK )
#endif

  END SUBROUTINE CheckConvergenceCoupled


  SUBROUTINE CheckConvergence_EmAb_NuE &
    ( MASK, k, Rtol, Utol, nIterations, &
      F_Y, U_Y, dU_Y, F_E, U_E, dU_E, FNRM0 )

    LOGICAL,  DIMENSION(1:nX_G), INTENT(inout) :: MASK
    INTEGER,                     INTENT(in)    :: k
    REAL(DP),                    INTENT(in)    :: Rtol, Utol
    INTEGER,  DIMENSION(1:nX_G), INTENT(inout) :: nIterations
    REAL(DP), DIMENSION(1:nX_G), INTENT(in)    :: F_Y, U_Y, dU_Y, F_E, U_E, dU_E, FNRM0

    LOGICAL  :: CONVERGED
    REAL(DP) :: FERR, UERR
    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, FERR, UERR )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, FERR, UERR )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( CONVERGED, FERR, UERR )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN

        FERR = SQRT( F_Y(iN_X)**2 + F_E(iN_X)**2 )

        UERR = SQRT( (dU_Y(iN_X)/U_Y(iN_X))**2 + (dU_E(iN_X)/U_E(iN_X))**2 )

        CONVERGED = FERR <= Rtol * FNRM0(iN_X) .OR. &
                    UERR <= Utol

        IF ( CONVERGED ) THEN
          MASK(iN_X) = .FALSE.
          nIterations(iN_X) = k
        END IF

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( MASK )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( MASK )
#endif

  END SUBROUTINE CheckConvergence_EmAb_NuE


END MODULE TwoMoment_NeutrinoMatterSolverModule
