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
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeEquilibriumDistributions_Points, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacities_NES_Point, &
    ComputeNeutrinoOpacities_Pair_Points, &
    ComputeNeutrinoOpacities_Pair_Point, &
    ComputeNeutrinoOpacitiesRates_NES_Points, &
    ComputeNeutrinoOpacitiesRates_Pair_Points, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeNeutrinoMatterSolver
  PUBLIC :: FinalizeNeutrinoMatterSolver
  PUBLIC :: SolveMatterEquations_EmAb
  PUBLIC :: SolveMatterEquations_EmAb_NuE
  PUBLIC :: SolveMatterEquations_EmAb_FP
  PUBLIC :: SolveMatterEquations_FP_Coupled
  PUBLIC :: SolveMatterEquations_FP_NestedAA

  ! --- Units Only for Displaying to Screen ---

  REAL(DP), PARAMETER :: Unit_D = Gram / Centimeter**3
  REAL(DP), PARAMETER :: Unit_T = MeV

  LOGICAL, PARAMETER :: SolveMatter = .TRUE.
  LOGICAL, PARAMETER :: UsePreconditionerEmAb = .TRUE.
  LOGICAL, PARAMETER :: UsePreconditionerPair = .FALSE.
  LOGICAL, PARAMETER :: UsePreconditionerPairLagAllButJ0 = .FALSE.

  INTEGER,  PARAMETER :: M_FP = 3
  INTEGER,  PARAMETER :: M_outer = 3
  INTEGER,  PARAMETER :: M_inner = 3
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
  INTEGER,  ALLOCATABLE :: INFO(:)
  INTEGER               :: LWORK

  ! --- Temporary arrays for scatter/gather (packing)

  REAL(DP), ALLOCATABLE, TARGET :: P1D_1(:)
  REAL(DP), ALLOCATABLE, TARGET :: P1D_2(:)
  REAL(DP), ALLOCATABLE, TARGET :: P1D_3(:)
  REAL(DP), ALLOCATABLE, TARGET :: P1D_4(:)

  REAL(DP), ALLOCATABLE, TARGET :: P2D_1 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_2 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_3 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_4 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_5 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_6 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_7 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_8 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_9 (:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P2D_10(:,:)

  REAL(DP), ALLOCATABLE, TARGET :: P3D_1(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_2(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_3(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_4(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_5(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_6(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_7(:,:,:)
  REAL(DP), ALLOCATABLE, TARGET :: P3D_8(:,:,:)

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

    W2_S(:) = WFactor_FP * W2_N(:)
    W3_S(:) = WFactor_FP * W3_N(:)

    ALLOCATE( AMAT(n_FP,M_FP,nX_G) )
    ALLOCATE( BVEC(n_FP,nX_G) )
    ALLOCATE( TAU (n_FP,nX_G) )
    ALLOCATE( INFO(nX_G) )

    ALLOCATE( P1D_1(nX_G) )
    ALLOCATE( P1D_2(nX_G) )
    ALLOCATE( P1D_3(nX_G) )
    ALLOCATE( P1D_4(nX_G) )

    ALLOCATE( P2D_1 (nE_G,nX_G) )
    ALLOCATE( P2D_2 (nE_G,nX_G) )
    ALLOCATE( P2D_3 (nE_G,nX_G) )
    ALLOCATE( P2D_4 (nE_G,nX_G) )
    ALLOCATE( P2D_5 (nE_G,nX_G) )
    ALLOCATE( P2D_6 (nE_G,nX_G) )
    ALLOCATE( P2D_7 (nE_G,nX_G) )
    ALLOCATE( P2D_8 (nE_G,nX_G) )
    ALLOCATE( P2D_9 (nE_G,nX_G) )
    ALLOCATE( P2D_10(nE_G,nX_G) )

    ALLOCATE( P3D_1(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_2(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_3(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_4(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_5(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_6(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_7(nE_G,nE_G,nX_G) )
    ALLOCATE( P3D_8(nE_G,nE_G,nX_G) )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$OMP MAP( alloc: AMAT, BVEC, TAU, INFO, &
    !$OMP             P1D_1, P1D_2, P1D_3, P1D_4, &
    !$OMP             P2D_1, P2D_2, P2D_3, P2D_4, P2D_5, P2D_6, P2D_7, P2D_8, P2D_9, P2D_10, &
    !$OMP             P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( E_N, W2_N, W3_N, W2_S, W3_S ) &
    !$ACC CREATE( AMAT, BVEC, TAU, INFO, &
    !$ACC         P1D_1, P1D_2, P1D_3, P1D_4, &
    !$ACC         P2D_1, P2D_2, P2D_3, P2D_4, P2D_5, P2D_6, P2D_7, P2D_8, P2D_9, P2D_10, &
    !$ACC         P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8 )
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

  END SUBROUTINE InitializeNeutrinoMatterSolver


  SUBROUTINE FinalizeNeutrinoMatterSolver

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: E_N, W2_N, W3_N, W2_S, W3_S, &
    !$OMP               AMAT, BVEC, TAU, WORK, INFO, &
    !$OMP               P1D_1, P1D_2, P1D_3, P1D_4, &
    !$OMP               P2D_1, P2D_2, P2D_3, P2D_4, P2D_5, P2D_6, P2D_7, P2D_8, P2D_9, P2D_10, &
    !$OMP               P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( E_N, W2_N, W3_N, W2_S, W3_S, &
    !$ACC         AMAT, BVEC, TAU, WORK, INFO, &
    !$ACC         P1D_1, P1D_2, P1D_3, P1D_4, &
    !$ACC         P2D_1, P2D_2, P2D_3, P2D_4, P2D_5, P2D_6, P2D_7, P2D_8, P2D_9, P2D_10, &
    !$ACC         P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8 )
#endif

    DEALLOCATE( E_N, W2_N, W3_N, W2_S, W3_S )
    DEALLOCATE( AMAT, BVEC, TAU, WORK, INFO )
    DEALLOCATE( P1D_1, P1D_2, P1D_3, P1D_4 )
    DEALLOCATE( P2D_1, P2D_2, P2D_3, P2D_4, P2D_5, P2D_6, P2D_7, P2D_8, P2D_9, P2D_10 )
    DEALLOCATE( P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8 )

  END SUBROUTINE FinalizeNeutrinoMatterSolver


  SUBROUTINE ComputePointsAndWeightsE( E, W2, W3 )

    REAL(DP), INTENT(out) :: E (:)
    REAL(DP), INTENT(out) :: W2(:)
    REAL(DP), INTENT(out) :: W3(:)

    INTEGER  :: iN_E, iE, iNodeE

    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      E (iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

      W2(iN_E) = WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**2
      W3(iN_E) = WeightsE(iNodeE) * MeshE % Width(iE) * E(iN_E)**3

    END DO

  END SUBROUTINE ComputePointsAndWeightsE


  SUBROUTINE SolveMatterEquations_EmAb &
    ( J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, D, T, Y, E )

    ! --- Electron Neutrinos (1) and Electron Antineutrinos (2) ---

    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(:),   INTENT(in)    :: D, T, Y, E

    CALL ComputeEquilibriumDistributions_Packed &
           ( iNuE, iNuE_Bar, D, T, Y, J0_1, J0_2 )

  END SUBROUTINE SolveMatterEquations_EmAb


  SUBROUTINE SolveMatterEquations_EmAb_NuE &
    ( dt, J, Chi, J0, D, T, Y, E, nIterations )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                 INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: J
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi
    REAL(DP), DIMENSION(:,:), INTENT(out)   :: J0
    REAL(DP), DIMENSION(:),   INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(:),   INTENT(out)   :: nIterations

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: MaxIter = 20
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Theta, Theta_J, Theta_J0
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Theta_1, Theta_2

    REAL(DP), DIMENSION(1:nX_G) :: N_B
    REAL(DP), DIMENSION(1:nX_G) :: Yold, C_Y, F_Y, U_Y, dU_Y
    REAL(DP), DIMENSION(1:nX_G) :: Eold, C_E, F_E, U_E, dU_E

    REAL(DP), DIMENSION(1:nX_G) :: Mnu, dMnudT, dMnudY
    REAL(DP), DIMENSION(1:nX_G) :: dEdD, dEdT, dEdY

    INTEGER,  DIMENSION(1:nX_G) :: PackIndex, UnpackIndex
    REAL(DP), DIMENSION(1:nX_G) :: D_P, T_P, Y_P, E_P
    REAL(DP), DIMENSION(1:nX_G) :: Mnu_P, dMnudT_P, dMnudY_P
    REAL(DP), DIMENSION(1:nX_G) :: dEdD_P, dEdT_P, dEdY_P

    REAL(DP), DIMENSION(1:nX_G) :: FJAC11, FJAC21, FJAC12, FJAC22
    REAL(DP), DIMENSION(1:nX_G) :: FNRM0
    LOGICAL,  DIMENSION(1:nX_G) :: ITERATE

    REAL(DP) :: kT, Eta
    REAL(DP) :: dJ0dT_Y, dJ0dY_T, dJ0dE_Y, dJ0dY_E
    REAL(DP) :: DJAC, FERR, UERR
    INTEGER  :: k, iN_X, iN_E, iX_P, nX_P

    ITERATE(:) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: ITERATE ) &
    !$OMP MAP( alloc: Theta, Theta_J, Theta_J0, Theta_1, Theta_2, N_B, &
    !$OMP             Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$OMP             Eold, C_E, F_E, U_E, dU_E, &
    !$OMP             Mnu, dMnudT, dMnudY, dEdD, dEdT, dEdY, &
    !$OMP             PackIndex, UnpackIndex, &
    !$OMP             D_P, T_P, Y_P, E_P, &
    !$OMP             Mnu_P, dMnudT_P, dMnudY_P, dEdD_P, dEdT_P, dEdY_P, &
    !$OMP             FJAC11, FJAC21, FJAC12, FJAC22, FNRM0 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( ITERATE ) &
    !$ACC CREATE( Theta, Theta_J, Theta_J0, Theta_1, Theta_2, N_B, &
    !$ACC         Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$ACC         Eold, C_E, F_E, U_E, dU_E, &
    !$ACC         Mnu, dMnudT, dMnudY, dEdD, dEdT, dEdY, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         D_P, T_P, Y_P, E_P, &
    !$ACC         Mnu_P, dMnudT_P, dMnudY_P, dEdD_P, dEdT_P, dEdY_P, &
    !$ACC         FJAC11, FJAC21, FJAC12, FJAC22, FNRM0 )
#endif

    ! --- Neutrino Chemical Potential and Derivatives ---

    CALL ComputeNeutrinoChemicalPotentials &
           ( D, T, Y, Mnu, dMnudT, dMnudY, iSpecies = iNuE )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( kT ) &
    !$ACC PRESENT( E_N, T, Mnu, Chi, Theta, J0, Theta_J0, J, Theta_J )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( kT )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G

        ! --- Equilibrium Distribution ---

        kT = BoltzmannConstant * T(iN_X)
        J0(iN_E,iN_X) = FermiDirac( E_N(iN_E), Mnu(iN_X), kT )

        ! --- Auxiliary Variables ---

        Theta   (iN_E,iN_X) = dt * Chi(iN_E,iN_X) / ( One + dt * Chi(iN_E,iN_X) )
        Theta_J0(iN_E,iN_X) = Theta(iN_E,iN_X) * J0(iN_E,iN_X)
        Theta_J (iN_E,iN_X) = Theta(iN_E,iN_X) * J(iN_E,iN_X)

      END DO
    END DO

    ! --- Old States (Constant) ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, Theta_J , nE_G, W2_S, 1, Zero, C_Y, 1 )

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, Theta_J , nE_G, W3_S, 1, Zero, C_E, 1 )

    ! --- Electron Fraction Equation ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, Theta_J0, nE_G, W2_S, 1, Zero, F_Y, 1 )

    ! --- Internal Energy Equation ---

    CALL MatrixVectorMultiply &
      ( 'T', nE_G, nX_G, One, Theta_J0, nE_G, W3_S, 1, Zero, F_E, 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( D, N_B, Y, Yold, C_Y, F_Y, U_Y, E, Eold, C_E, F_E, U_E, FNRM0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      N_B(iN_X) = D(iN_X) / AtomicMassUnit

      ! --- Initial Guess ---

      Yold(iN_X) = Y(iN_X)
      Eold(iN_X) = E(iN_X)

      U_Y(iN_X) = Yold(iN_X)
      U_E(iN_X) = Eold(iN_X)

      ! --- Scale Equations and Save Initial Evaluation ---

      ! --- Electron Fraction Equation ---

      C_Y(iN_X) = C_Y(iN_X) + N_B(iN_X) * U_Y(iN_X)
      F_Y(iN_X) = F_Y(iN_X) + N_B(iN_X) * U_Y(iN_X) - C_Y(iN_X)
      F_Y(iN_X) = F_Y(iN_X) / C_Y(iN_X)

      ! --- Internal Energy Equation ---

      C_E(iN_X) = C_E(iN_X) + D  (iN_X) * U_E(iN_X)
      F_E(iN_X) = F_E(iN_X) + D  (iN_X) * U_E(iN_X) - C_E(iN_X)
      F_E(iN_X) = F_E(iN_X) / C_E(iN_X)

      FNRM0(iN_X) = SQRT( F_Y(iN_X)**2 + F_E(iN_X)**2 )

    END DO

    k = 0
    DO WHILE( ANY( ITERATE(:) ) .AND. k < MaxIter )

      k  = k + 1

      CALL CreatePackIndex &
             ( nX_G, ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Internal Energy Derivatives ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iN_X )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iN_X ) &
      !$ACC PRESENT( UnpackIndex, D, T, Y, D_P, T_P, Y_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( iN_X )
#endif
      DO iX_P = 1, nX_P
        iN_X = UnpackIndex(iX_P)
        D_P(iX_P) = D(iN_X)
        T_P(iX_P) = T(iN_X)
        Y_P(iX_P) = Y(iN_X)
      END DO

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), E_P(1:nX_P), &
               dEdD_P(1:nX_P), dEdT_P(1:nX_P), dEdY_P(1:nX_P) )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( iX_P )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( iX_P ) &
      !$ACC PRESENT( ITERATE, PackIndex, dEdT, dEdY, dEdT_P, dEdY_P )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( iX_P )
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN
          iX_P = PackIndex(iN_X)
          dEdT(iN_X) = dEdT_P(iX_P)
          dEdY(iN_X) = dEdY_P(iX_P)
        END IF
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( kT, dJ0dT_Y, dJ0dY_T, dJ0dE_Y, dJ0dY_E )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( kT, dJ0dT_Y, dJ0dY_T, dJ0dE_Y, dJ0dY_E ) &
      !$ACC PRESENT( ITERATE, T, E_N, Mnu, dMnudT, dMnudY, &
      !$ACC          dEdT, dEdY, Theta, Theta_1, Theta_2 )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( kT, dJ0dT_Y, dJ0dY_T, dJ0dE_Y, dJ0dY_E )
#endif
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G
          IF ( ITERATE(iN_X) ) THEN

            kT = BoltzmannConstant * T(iN_X)

            ! --- Derivative of J0 wrt. T (Constant Y) ---

            dJ0dT_Y = dFermiDiracdT( E_N(iN_E), Mnu(iN_X), kT, dMnudT(iN_X), T(iN_X) )

            ! --- Derivative of J0 wrt. Y (Constant T) ---

            dJ0dY_T = dFermiDiracdY( E_N(iN_E), Mnu(iN_X), kT, dMnudY(iN_X), T(iN_X) )

            ! --- Derivative of J0 wrt. E (Constant Y) ---

            dJ0dE_Y = dJ0dT_Y / dEdT(iN_X)

            ! --- Derivative of J0 wrt. Y (Constant E) ---

            dJ0dY_E = dJ0dY_T - dJ0dE_Y * dEdY(iN_X)

            ! --- Auxiliary Variables ---

            Theta_1(iN_E,iN_X) = Theta(iN_E,iN_X) * dJ0dY_E
            Theta_2(iN_E,iN_X) = Theta(iN_E,iN_X) * dJ0dE_Y

          END IF
        END DO
      END DO

      ! --- Jacobian ---

      CALL MatrixVectorMultiply &
        ( 'T', nE_G, nX_G, One, Theta_1, nE_G, W2_S, 1, Zero, FJAC11, 1 )
      CALL MatrixVectorMultiply &
        ( 'T', nE_G, nX_G, One, Theta_1, nE_G, W3_S, 1, Zero, FJAC21, 1 )
      CALL MatrixVectorMultiply &
        ( 'T', nE_G, nX_G, One, Theta_2, nE_G, W2_S, 1, Zero, FJAC12, 1 )
      CALL MatrixVectorMultiply &
        ( 'T', nE_G, nX_G, One, Theta_2, nE_G, W3_S, 1, Zero, FJAC22, 1 )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( DJAC )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( DJAC ) &
      !$ACC PRESENT( ITERATE, D, N_B, FJAC11, FJAC21, FJAC12, FJAC22, &
      !$ACC          C_Y, F_Y, U_Y, dU_Y, C_E, F_E, U_E, dU_E )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( DJAC )
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          ! --- Scale Jacobian ---

          FJAC11(iN_X) = ( FJAC11(iN_X) + N_B(iN_X) ) / C_Y(iN_X)
          FJAC12(iN_X) = FJAC12(iN_X) / C_Y(iN_X)

          FJAC21(iN_X) = FJAC21(iN_X) / C_E(iN_X)
          FJAC22(iN_X) = ( FJAC22(iN_X) + D  (iN_X) ) / C_E(iN_X)

          ! --- Determinant of Jacobian ---

          DJAC = FJAC11(iN_X) * FJAC22(iN_X) - FJAC21(iN_X) * FJAC12(iN_X)

          ! --- Correction ---

          dU_Y(iN_X) = - ( + FJAC22(iN_X) * F_Y(iN_X) - FJAC12(iN_X) * F_E(iN_X) ) / DJAC
          dU_E(iN_X) = - ( - FJAC21(iN_X) * F_Y(iN_X) + FJAC11(iN_X) * F_E(iN_X) ) / DJAC

          ! --- Apply Correction ---

          U_Y(iN_X) = U_Y(iN_X) + dU_Y(iN_X)
          U_E(iN_X) = U_E(iN_X) + dU_E(iN_X)

          Y(iN_X) = U_Y(iN_X)
          E(iN_X) = U_E(iN_X)

          !F_Y(iN_X) = N_B(iN_X) * U_Y(iN_X) - C_Y(iN_X)
          !F_E(iN_X) = D  (iN_X) * U_E(iN_X) - C_E(iN_X)

        END IF
      END DO

      ! --- Update Temperature ---

      CALL ArrayPack &
             ( nX_G, nX_P, UnpackIndex, Y, E, Y_P, E_P )

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D_P(1:nX_P), E_P(1:nX_P), Y_P(1:nX_P), T_P(1:nX_P) )

      ! --- Neutrino Chemical Potential and Derivatives ---

      CALL ComputeNeutrinoChemicalPotentials &
             ( D_P(1:nX_P), T_P(1:nX_P), Y_P(1:nX_P), Mnu_P(1:nX_P), &
               dMnudT_P(1:nX_P), dMnudY_P(1:nX_P), iSpecies = iNuE )

      CALL ArrayUnpack &
             ( nX_G, nX_P, ITERATE, PackIndex, &
               T_P, Mnu_P, dMnudT_P, dMnudY_P, T, Mnu, dMnudT, dMnudY )

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( kT )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( kT ) &
      !$ACC PRESENT( ITERATE, E_N, T, Mnu, J0, Theta_J0, Theta )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( kT )
#endif
      DO iN_X = 1, nX_G
        DO iN_E = 1, nE_G
          IF ( ITERATE(iN_X) ) THEN

            ! --- Equilibrium Distribution ---

            kT = BoltzmannConstant * T(iN_X)
            J0(iN_E,iN_X) = FermiDirac( E_N(iN_E), Mnu(iN_X), kT )

            ! --- Auxiliary Variables ---

            Theta_J0(iN_E,iN_X) = Theta(iN_E,iN_X) * J0(iN_E,iN_X)

          END IF
        END DO
      END DO

      ! --- Electron Fraction Equation ---

      CALL MatrixVectorMultiply &
        ( 'T', nE_G, nX_G, One, Theta_J0, nE_G, W2_S, 1, Zero, F_Y, 1 )

      ! --- Internal Energy Equation ---

      CALL MatrixVectorMultiply &
        ( 'T', nE_G, nX_G, One, Theta_J0, nE_G, W3_S, 1, Zero, F_E, 1 )

      ! --- Check for Convergence ---

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( FERR, UERR )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( FERR, UERR ) &
      !$ACC PRESENT( ITERATE, nIterations, FNRM0, D, N_B, &
      !$ACC          C_Y, F_Y, U_Y, dU_Y, C_E, F_E, U_E, dU_E )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( FERR, UERR )
#endif
      DO iN_X = 1, nX_G
        IF ( ITERATE(iN_X) ) THEN

          F_Y(iN_X) = F_Y(iN_X) + N_B(iN_X) * U_Y(iN_X) - C_Y(iN_X)
          F_Y(iN_X) = F_Y(iN_X) / C_Y(iN_X)

          F_E(iN_X) = F_E(iN_X) + D  (iN_X) * U_E(iN_X) - C_E(iN_X)
          F_E(iN_X) = F_E(iN_X) / C_E(iN_X)

          FERR = SQRT( F_Y(iN_X)**2 + F_E(iN_X)**2 )

          UERR = SQRT( (dU_Y(iN_X)/U_Y(iN_X))**2 + (dU_E(iN_X)/U_E(iN_X))**2 )

          ITERATE(iN_X) = FERR > Rtol * FNRM0(iN_X) &
                    .AND. UERR > Utol

          IF ( .NOT. ITERATE(iN_X) ) nIterations(iN_X) = k

        END IF
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET UPDATE FROM( ITERATE )
#elif defined(THORNADO_OACC)
      !$ACC UPDATE HOST( ITERATE )
#endif

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta ) &
    !$ACC PRESENT( Chi, J0, J )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#endif
    DO iN_X = 1, nX_G
      DO iN_E = 1, nE_G
        Eta = Chi(iN_E,iN_X) * J0(iN_E,iN_X)
        J(iN_E,iN_X) = ( J(iN_E,iN_X) + dt * Eta ) / ( One + dt * Chi(iN_E,iN_X) )
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: ITERATE, &
    !$OMP               Theta, Theta_J, Theta_J0, Theta_1, Theta_2, N_B, &
    !$OMP               Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$OMP               Eold, C_E, F_E, U_E, dU_E, &
    !$OMP               Mnu, dMnudT, dMnudY, dEdD, dEdT, dEdY, &
    !$OMP               PackIndex, UnpackIndex, &
    !$OMP               D_P, T_P, Y_P, E_P, &
    !$OMP               Mnu_P, dMnudT_P, dMnudY_P, dEdD_P, dEdT_P, dEdY_P, &
    !$OMP               FJAC11, FJAC21, FJAC12, FJAC22, FNRM0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( ITERATE, &
    !$ACC         Theta, Theta_J, Theta_J0, Theta_1, Theta_2, N_B, &
    !$ACC         Yold, C_Y, F_Y, U_Y, dU_Y, &
    !$ACC         Eold, C_E, F_E, U_E, dU_E, &
    !$ACC         Mnu, dMnudT, dMnudY, dEdD, dEdT, dEdY, &
    !$ACC         PackIndex, UnpackIndex, &
    !$ACC         D_P, T_P, Y_P, E_P, &
    !$ACC         Mnu_P, dMnudT_P, dMnudY_P, dEdD_P, dEdT_P, dEdY_P, &
    !$ACC         FJAC11, FJAC21, FJAC12, FJAC22, FNRM0 )
#endif

  END SUBROUTINE SolveMatterEquations_EmAb_NuE


  SUBROUTINE SolveMatterEquations_EmAb_FP &
    ( dt, iS_1, iS_2, J_1, J_2, Chi_1, Chi_2, J0_1, J0_2, &
      D, T, Y, E, nIterations_Out, TOL )

    ! --- Neutrino (1) and Antineutrino (2) ---

    REAL(DP),                 INTENT(in)    :: dt
    INTEGER,                  INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(:),   INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(:),   INTENT(out), OPTIONAL :: nIterations_Out(:)
    REAL(DP),                 INTENT(in),  OPTIONAL :: TOL

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = iE
    INTEGER,  PARAMETER :: MaxIter = 100

    ! --- Local Variables ---

    REAL(DP), DIMENSION(       1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(       1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_1, Jnew_1
    REAL(DP), DIMENSION(1:nE_G,1:nX_G) :: Jold_2, Jnew_2


    REAL(DP), DIMENSION(1:n_FP,1:M_FP,1:nX_G) :: GVEC, FVEC
    REAL(DP), DIMENSION(1:n_FP,       1:nX_G) :: GVECm, FVECm
    REAL(DP), DIMENSION(       1:M_FP,1:nX_G) :: Alpha

    REAL(DP), DIMENSION(1:nX_G) :: Jnorm_1, Jnorm_2

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

    CALL ArrayCopy( nX_G, Y, E, Yold, Eold )

    CALL ArrayCopy( nE_G, nX_G, J_1, J_2, Jold_1, Jold_2 )

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
             ( nX_G, ITERATE, nX_P, PackIndex, UnpackIndex )

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

      CALL StoreRHS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVECm, GVECm, FVEC, GVEC )

      ! --- Anderson Acceleration ---

      CALL TimersStart( Timer_Im_ComputeLS )

      IF ( Mk > 1 ) THEN

        CALL BuildLS_FP &
               ( ITERATE, n_FP, M_FP, Mk, FVECm, FVEC, AMAT, BVEC )

        CALL SolveLS_FP &
               ( ITERATE, n_FP, M_FP, Mk, AMAT, BVEC, Alpha )

        CALL ApplyAA_FP &
               ( ITERATE, n_FP, M_FP, Mk, Alpha, GVEC, GVECm )

      END IF

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors---

      CALL TimersStart( Timer_Im_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm, GVECm )

      CALL UpdateNeutrinoRHS_FP &
             ( ITERATE, n_FP, OS_1, OS_2, &
               FVECm, GVECm, Jnew_1, Jnew_2 )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Check Convergence ---

      CALL CheckConvergenceCoupled_FP &
             ( ITERATE, n_FP, k, iY, iE, OS_1, OS_2, Rtol, &
               nIterations, FVECm, Jnorm_1, Jnorm_2 )

      ! --- Shift History Arrays ---

      CALL ShiftRHS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVEC, GVEC )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    CALL ArrayCopy( nE_G, nX_G, Jnew_1, Jnew_2, J_1, J_2 )

    IF(PRESENT(nIterations_Out)) THEN
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRESENT( nIterations_Out, nIterations )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
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

    REAL(DP),                 INTENT(in)    :: dt
    INTEGER,                  INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Eta_Pair_1, Eta_Pair_2
    REAL(DP), DIMENSION(:),   INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(:),   INTENT(out)   :: nIterations

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

    REAL(DP), DIMENSION(1:n_FP,1:M_FP,1:nX_G) :: GVEC, FVEC
    REAL(DP), DIMENSION(1:n_FP,       1:nX_G) :: GVECm, FVECm
    REAL(DP), DIMENSION(       1:M_FP,1:nX_G) :: Alpha

    REAL(DP), DIMENSION(1:nX_G) :: Jnorm_1, Jnorm_2

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

    CALL ArrayCopy( nX_G, Y, E, Yold, Eold )

    CALL ArrayCopy( nE_G, nX_G, J_1, J_2, Jold_1, Jold_2 )

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
             ( nX_G, ITERATE, nX_P, PackIndex, UnpackIndex )

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

      CALL StoreRHS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVECm, GVECm, FVEC, GVEC )

      ! --- Anderson Acceleration ---

      CALL TimersStart( Timer_Im_ComputeLS )

      IF ( Mk > 1 ) THEN

        CALL BuildLS_FP &
               ( ITERATE, n_FP, M_FP, Mk, FVECm, FVEC, AMAT, BVEC )

        CALL SolveLS_FP &
               ( ITERATE, n_FP, M_FP, Mk, AMAT, BVEC, Alpha )

        CALL ApplyAA_FP &
               ( ITERATE, n_FP, M_FP, Mk, Alpha, GVEC, GVECm )

      END IF

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors---

      CALL TimersStart( Timer_Im_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm, GVECm )

      CALL UpdateNeutrinoRHS_FP &
             ( ITERATE, n_FP, OS_1, OS_2, &
               FVECm, GVECm, Jnew_1, Jnew_2 )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE, nX_P, PackIndex, UnpackIndex )

      ! --- Check Convergence ---

      CALL CheckConvergenceCoupled_FP &
             ( ITERATE, n_FP, k, iY, iE, OS_1, OS_2, Rtol, &
               nIterations, FVECm, Jnorm_1, Jnorm_2 )

      ! --- Shift History Arrays ---

      CALL ShiftRHS_FP &
             ( ITERATE, n_FP, M_FP, Mk, FVEC, GVEC )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    CALL ArrayCopy( nE_G, nX_G, Jnew_1, Jnew_2, J_1, J_2 )

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

    REAL(DP),                 INTENT(in)    :: dt
    INTEGER,                  INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: J_1, J_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(out)   :: J0_1, J0_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Eta_Pair_1, Eta_Pair_2
    REAL(DP), DIMENSION(:),   INTENT(inout) :: D, T, Y, E
    INTEGER,  DIMENSION(:),   INTENT(out)   :: nIterations_Inner, nIterations_Outer

    ! --- Solver Parameters ---

    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iE = 2
    INTEGER,  PARAMETER :: OS_1 = 0
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-08
    REAL(DP), PARAMETER :: Utol = 1.0d-10

    ! --- Local Variables ---

    REAL(DP), DIMENSION(              1:nX_G) :: Yold, S_Y, C_Y, Unew_Y, GVEC_Y
    REAL(DP), DIMENSION(              1:nX_G) :: Eold, S_E, C_E, Unew_E, GVEC_E
    REAL(DP), DIMENSION(1:nE_G,       1:nX_G) :: Jold_1, Jold_2, Jnew_1, Jnew_2
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

    REAL(DP), DIMENSION(1:nX_G) :: Jnorm_1, Jnorm_2

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

    CALL ArrayCopy( nX_G, Y, E, Yold, Eold )

    CALL ArrayCopy( nE_G, nX_G, J_1, J_2, Jold_1, Jold_2 )

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
             ( nX_G, ITERATE_OUTER, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

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
               ( nX_G, ITERATE_INNER, nX_P_inner, PackIndex_inner, UnpackIndex_inner )

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

        CALL StoreRHS_FP &
               ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                 FVECm_inner, GVECm_inner, FVEC_inner, GVEC_inner )

        ! --- Anderson Acceleration (inner) ---

        CALL TimersStart( Timer_Im_ComputeLS )

        IF ( Mk_inner > 1 ) THEN

          CALL BuildLS_FP &
                 ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                   FVECm_inner, FVEC_inner, AMAT_inner, BVEC_inner )

          CALL SolveLS_FP &
                 ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                   AMAT_inner, BVEC_inner, Alpha_inner )

          CALL ApplyAA_FP &
                 ( ITERATE_INNER, n_FP_inner, M_inner, Mk_inner, &
                   Alpha_inner, GVEC_inner, GVECm_inner )

        END IF

        CALL TimersStop( Timer_Im_ComputeLS )

        ! --- Update Residuals and Solution Vectors (inner) ---

        CALL TimersStart( Timer_Im_UpdateFP )

        CALL UpdateNeutrinoRHS_FP &
               ( ITERATE_INNER, n_FP_inner, OS_1, OS_2, &
                 FVECm_inner, GVECm_inner, Jnew_1, Jnew_2 )

        ! --- Check Convergence (inner) ---

        CALL CheckConvergenceInner_FP &
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

      CALL StoreRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
               FVECm_outer, GVECm_outer, FVEC_outer, GVEC_outer )

      ! --- Anderson Acceleration (outer) ---

      CALL TimersStart( Timer_Im_ComputeLS )

      IF ( Mk_outer > 1 ) THEN

        CALL BuildLS_FP &
               ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
                 FVECm_outer, FVEC_outer, AMAT_outer, BVEC_outer )

        CALL SolveLS_FP &
               ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
                 AMAT_outer, BVEC_outer, Alpha_outer )

        CALL ApplyAA_FP &
               ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, &
                 Alpha_outer, GVEC_outer, GVECm_outer )

      END IF

      CALL TimersStop( Timer_Im_ComputeLS )

      ! --- Update Residuals and Solution Vectors (outer) ---

      CALL TimersStart( Timer_Im_UpdateFP )

      CALL UpdateMatterRHS_FP &
             ( ITERATE_OUTER, iY, iE, Yold, Eold, Y, E, &
               Unew_Y, Unew_E, FVECm_outer, GVECm_outer )

      ! --- Update Temperature ---

      CALL UpdateTemperature_Packed &
             ( D, E, Y, T, &
               ITERATE_outer, nX_P_outer, PackIndex_outer, UnpackIndex_outer )

      ! --- Check Convergence (outer) ---

      CALL CheckConvergenceOuter_FP &
             ( ITERATE_OUTER, ITERATE_INNER, n_FP_outer, iY, iE, k_outer, &
               FVECm_outer, Rtol, nIterations_Outer )

      ! --- Shift History Arrays (outer) ---

      CALL ShiftRHS_FP &
             ( ITERATE_OUTER, n_FP_outer, M_outer, Mk_outer, FVEC_outer, GVEC_outer )

      CALL TimersStop( Timer_Im_UpdateFP )

    END DO

    CALL ArrayCopy( nE_G, nX_G, Jnew_1, Jnew_2, J_1, J_2 )

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

    INTEGER,                              INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(:),     TARGET,   INTENT(in)    :: D, T, Y
    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(inout) :: J0_1, J0_2
    LOGICAL,  DIMENSION(:),     OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                    OPTIONAL, INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(:),     OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:),     POINTER,  CONTIGUOUS    :: D_P, T_P, Y_P
    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: J0_1_P, J0_2_P

    INTEGER :: nX

    IF ( PRESENT( nX_P ) ) THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      CALL ArrayPack &
             ( nX_G, nX, UnpackIndex, &
               D, T, Y, &
               P1D_1, P1D_2, P1D_3 )

      D_P => P1D_1(1:nX)
      T_P => P1D_2(1:nX)
      Y_P => P1D_3(1:nX)

      J0_1_P => P2D_1(:,1:nX)
      J0_2_P => P2D_2(:,1:nX)

    ELSE

      D_P => D(:)
      T_P => T(:)
      Y_P => Y(:)

      J0_1_P => J0_1(:,:)
      J0_2_P => J0_2(:,:)

    END IF

    ! --- Equilibrium Distributions ---

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, iS_1 )

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_2_P, iS_2 )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nE_G, nX_G, nX, MASK, PackIndex, &
               P2D_1, P2D_2, &
               J0_1, J0_2 )

    END IF

    D_P => NULL()
    T_P => NULL()
    Y_P => NULL()

    J0_1_P => NULL()
    J0_2_P => NULL()

  END SUBROUTINE ComputeEquilibriumDistributions_Packed


  SUBROUTINE ComputeOpacities_Packed &
    ( iS_1, iS_2, D, T, Y, J0_1, J0_2, &
      Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
      Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
      MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    INTEGER,                              INTENT(in)    :: iS_1, iS_2
    REAL(DP), DIMENSION(:),     TARGET,   INTENT(in)    :: D, T, Y
    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(inout) :: J0_1, J0_2
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(inout) :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(inout) :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(inout) :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(inout) :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2
    LOGICAL,  DIMENSION(:),     OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                    OPTIONAL, INTENT(in)    :: nX_P, nX_P0
    INTEGER,  DIMENSION(:),     OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:),     POINTER,  CONTIGUOUS    :: D_P, T_P, Y_P
    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: J0_1_P, J0_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P

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

      CALL ArrayPack &
             ( nX_G, nX, UnpackIndex, &
               D, T, Y, &
               P1D_1, P1D_2, P1D_3 )

      D_P => P1D_1(1:nX)
      T_P => P1D_2(1:nX)
      Y_P => P1D_3(1:nX)

      J0_1_P => P2D_1(:,1:nX)
      J0_2_P => P2D_2(:,1:nX)

      Phi_0_In_NES_1_P  => P3D_1(:,:,1:nX)
      Phi_0_Ot_NES_1_P  => P3D_2(:,:,1:nX)
      Phi_0_In_NES_2_P  => P3D_3(:,:,1:nX)
      Phi_0_Ot_NES_2_P  => P3D_4(:,:,1:nX)
      Phi_0_In_Pair_1_P => P3D_5(:,:,1:nX)
      Phi_0_Ot_Pair_1_P => P3D_6(:,:,1:nX)
      Phi_0_In_Pair_2_P => P3D_7(:,:,1:nX)
      Phi_0_Ot_Pair_2_P => P3D_8(:,:,1:nX)

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

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_1_P, iS_1 )

    CALL ComputeEquilibriumDistributions_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, J0_2_P, iS_2 )

    ! --- NES Kernels ---

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, 1, &
             Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P )

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_2, 1, &
             Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P )

    ! --- Pair Kernels ---

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_1, 1, &
             Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P )

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nE_G, 1, nX, E_N, D_P, T_P, Y_P, iS_2, 1, &
             Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nE_G, nX_G, nX, MASK, PackIndex, &
               P2D_1, P2D_2, &
               J0_1, J0_2 )

      IF ( nX < nX0 ) THEN

        CALL ArrayUnpack &
               ( nE_G, nE_G, nX_G, nX, MASK, PackIndex, &
                 P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2 )

      END IF

    END IF

    D_P => NULL()
    T_P => NULL()
    Y_P => NULL()

    J0_1_P => NULL()
    J0_2_P => NULL()

    Phi_0_In_NES_1_P  => NULL()
    Phi_0_Ot_NES_1_P  => NULL()
    Phi_0_In_NES_2_P  => NULL()
    Phi_0_Ot_NES_2_P  => NULL()
    Phi_0_In_Pair_1_P => NULL()
    Phi_0_Ot_Pair_1_P => NULL()
    Phi_0_In_Pair_2_P => NULL()
    Phi_0_Ot_Pair_2_P => NULL()

  END SUBROUTINE ComputeOpacities_Packed


  SUBROUTINE ComputeRates_Packed &
    ( J_1, J_2, &
      Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
      Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2, &
      MASK, nX_P, PackIndex, UnpackIndex, nX_P0 )

    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(in)    :: Phi_0_In_NES_1, Phi_0_Ot_NES_1
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(in)    :: Phi_0_In_NES_2, Phi_0_Ot_NES_2
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(in)    :: Phi_0_In_Pair_1, Phi_0_Ot_Pair_1
    REAL(DP), DIMENSION(:,:,:), TARGET,   INTENT(in)    :: Phi_0_In_Pair_2, Phi_0_Ot_Pair_2
    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(inout) :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(inout) :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(inout) :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(:,:),   TARGET,   INTENT(inout) :: Eta_Pair_1, Eta_Pair_2

    LOGICAL,  DIMENSION(:),     OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                    OPTIONAL, INTENT(in)    :: nX_P, nX_P0
    INTEGER,  DIMENSION(:),     OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: J_1_P, J_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_NES_1_P, Phi_0_Ot_NES_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_NES_2_P, Phi_0_Ot_NES_2_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_Pair_1_P, Phi_0_Ot_Pair_1_P
    REAL(DP), DIMENSION(:,:,:), POINTER,  CONTIGUOUS    :: Phi_0_In_Pair_2_P, Phi_0_Ot_Pair_2_P
    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: Chi_NES_1_P, Chi_NES_2_P
    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: Eta_NES_1_P, Eta_NES_2_P
    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: Chi_Pair_1_P, Chi_Pair_2_P
    REAL(DP), DIMENSION(:,:),   POINTER,  CONTIGUOUS    :: Eta_Pair_1_P, Eta_Pair_2_P

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

      Chi_NES_1_P  => P2D_1(:,1:nX)
      Chi_NES_2_P  => P2D_2(:,1:nX)
      Eta_NES_1_P  => P2D_3(:,1:nX)
      Eta_NES_2_P  => P2D_4(:,1:nX)
      Chi_Pair_1_P => P2D_5(:,1:nX)
      Chi_Pair_2_P => P2D_6(:,1:nX)
      Eta_Pair_1_P => P2D_7(:,1:nX)
      Eta_Pair_2_P => P2D_8(:,1:nX)

      CALL ArrayPack &
             ( nE_G, nX_G, nX, UnpackIndex, &
               J_1, J_2, &
               P2D_9, P2D_10 )

      J_1_P => P2D_9 (:,1:nX)
      J_2_P => P2D_10(:,1:nX)

      IF ( nX < nX0 ) THEN

        CALL ArrayPack &
               ( nE_G, nE_G, nX_G, nX, UnpackIndex, &
                 Phi_0_In_NES_1, Phi_0_Ot_NES_1, Phi_0_In_NES_2, Phi_0_Ot_NES_2, &
                 Phi_0_In_Pair_1, Phi_0_Ot_Pair_1, Phi_0_In_Pair_2, Phi_0_Ot_Pair_2, &
                 P3D_1, P3D_2, P3D_3, P3D_4, P3D_5, P3D_6, P3D_7, P3D_8 )

      END IF

      Phi_0_In_NES_1_P  => P3D_1(:,:,1:nX)
      Phi_0_Ot_NES_1_P  => P3D_2(:,:,1:nX)
      Phi_0_In_NES_2_P  => P3D_3(:,:,1:nX)
      Phi_0_Ot_NES_2_P  => P3D_4(:,:,1:nX)
      Phi_0_In_Pair_1_P => P3D_5(:,:,1:nX)
      Phi_0_Ot_Pair_1_P => P3D_6(:,:,1:nX)
      Phi_0_In_Pair_2_P => P3D_7(:,:,1:nX)
      Phi_0_Ot_Pair_2_P => P3D_8(:,:,1:nX)

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
             ( nE_G, nX_G, nX, MASK, PackIndex, &
               P2D_1, P2D_2, P2D_3, P2D_4, P2D_5, P2D_6, P2D_7, P2D_8, &
               Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
               Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    END IF

    Chi_NES_1_P  => NULL()
    Chi_NES_2_P  => NULL()
    Eta_NES_1_P  => NULL()
    Eta_NES_2_P  => NULL()
    Chi_Pair_1_P => NULL()
    Chi_Pair_2_P => NULL()
    Eta_Pair_1_P => NULL()
    Eta_Pair_2_P => NULL()

    J_1_P => NULL()
    J_2_P => NULL()

    Phi_0_In_NES_1_P  => NULL()
    Phi_0_Ot_NES_1_P  => NULL()
    Phi_0_In_NES_2_P  => NULL()
    Phi_0_Ot_NES_2_P  => NULL()
    Phi_0_In_Pair_1_P => NULL()
    Phi_0_Ot_Pair_1_P => NULL()
    Phi_0_In_Pair_2_P => NULL()
    Phi_0_Ot_Pair_2_P => NULL()

  END SUBROUTINE ComputeRates_Packed


  SUBROUTINE UpdateTemperature_Packed &
    ( D, E, Y, T, MASK, nX_P, PackIndex, UnpackIndex )

    REAL(DP), DIMENSION(:), TARGET,   INTENT(in)    :: D, E, Y
    REAL(DP), DIMENSION(:), TARGET,   INTENT(inout) :: T

    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(in)    :: MASK
    INTEGER,                OPTIONAL, INTENT(in)    :: nX_P
    INTEGER,  DIMENSION(:), OPTIONAL, INTENT(in)    :: PackIndex, UnpackIndex

    REAL(DP), DIMENSION(:), POINTER,  CONTIGUOUS    :: D_P, E_P, Y_P, T_P

    INTEGER  :: nX

    IF ( PRESENT( nX_P ) ) THEN
      nX = nX_P
    ELSE
      nX = nX_G
    END IF

    IF ( nX < nX_G ) THEN

      ! --- Pack Arrays ---

      CALL ArrayPack &
             ( nX_G, nX, UnpackIndex, &
               D, E, Y, T, &
               P1D_1, P1D_2, P1D_3, P1D_4 )

      D_P => P1D_1(1:nX)
      E_P => P1D_2(1:nX)
      Y_P => P1D_3(1:nX)
      T_P => P1D_4(1:nX)

    ELSE

      D_P => D(:)
      E_P => E(:)
      Y_P => Y(:)
      T_P => T(:)

    END IF

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D_P, E_P, Y_P, T_P )

    IF ( nX < nX_G ) THEN

      ! --- Unpack Results ---

      CALL ArrayUnpack &
             ( nX_G, nX, MASK, PackIndex, &
               P1D_4, &
               T )

    END IF

    D_P => NULL()
    E_P => NULL()
    Y_P => NULL()
    T_P => NULL()

  END SUBROUTINE UpdateTemperature_Packed


  SUBROUTINE InitializeRHS_FP &
    ( Jold_1, Jold_2, J_1, J_2, Jnew_1, Jnew_2, &
      D, Y, Yold, C_Y, S_Y, U_Y, E, Eold, C_E, S_E, U_E )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: Jold_1, Jold_2
    REAL(DP), DIMENSION(:,:), INTENT(in)  :: J_1, J_2
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(:),   INTENT(in)  :: D
    REAL(DP), DIMENSION(:),   INTENT(in)  :: Y, Yold
    REAL(DP), DIMENSION(:),   INTENT(out) :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, Eold
    REAL(DP), DIMENSION(:),   INTENT(out) :: C_E, S_E, U_E

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
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( Y, Yold, S_Y, C_Y, U_Y, &
    !$ACC          E, Eold, S_E, C_E, U_E, D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G

      S_Y(iN_X) = One / ( D(iN_X) * Yold(iN_X) / AtomicMassUnit )
      S_E(iN_X) = One / ( D(iN_X) * Eold(iN_X) )

      C_Y(iN_X) = C_Y(iN_X) * S_Y(iN_X)
      C_E(iN_X) = C_E(iN_X) * S_E(iN_X)

      U_Y(iN_X) = Y(iN_X) / Yold(iN_X) ! --- Initial Guess
      U_E(iN_X) = E(iN_X) / Eold(iN_X) ! --- Initial Guess

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

  END SUBROUTINE InitializeRHS_FP


  SUBROUTINE ComputeMatterRHS_FP &
    ( MASK, n_FP, iY, iE, Fm, Gm, J_1, J_2, C_Y, S_Y, U_Y, G_Y, C_E, S_E, U_E, G_E )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    INTEGER,                  INTENT(in)    :: n_FP, iY, iE
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(:),   INTENT(in)    :: C_Y, S_Y, U_Y
    REAL(DP), DIMENSION(:),   INTENT(in)    :: C_E, S_E, U_E
    REAL(DP), DIMENSION(:),   INTENT(out)   :: G_Y, G_E

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
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, C_Y, S_Y, U_Y, G_Y, C_E, S_E, U_E, G_E, Fm, Gm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
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


  SUBROUTINE UpdateMatterRHS_FP &
    ( MASK, iY, iE, Yold, Eold, Y, E, U_Y, U_E, Fm, Gm )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    INTEGER,                  INTENT(in)    :: iY, iE
    REAL(DP), DIMENSION(:),   INTENT(in)    :: Yold, Eold
    REAL(DP), DIMENSION(:),   INTENT(inout) :: Y, E
    REAL(DP), DIMENSION(:),   INTENT(inout) :: U_Y, U_E
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Fm, Gm

    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, Fm, Gm, Y, E, Yold, Eold, U_Y, U_E )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
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


  SUBROUTINE ComputeNumberDensity_FP &
    ( MASK, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    REAL(DP),                 INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J0_1, J0_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Eta_Pair_1, Eta_Pair_2

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

  END SUBROUTINE ComputeNumberDensity_FP


  SUBROUTINE ComputeNumberDensity_EmAb_FP &
    ( MASK, dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2 )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    REAL(DP),                 INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J0_1, J0_2

    REAL(DP) :: Eta
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta ) &
    !$ACC PRESENT( MASK, &
    !$ACC          Jold_1, Jnew_1, J0_1, Chi_1, &
    !$ACC          Jold_2, Jnew_2, J0_2, Chi_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta, Eta_T, Chi_T )
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


  SUBROUTINE ComputeNeutrinoRHS_FP &
    ( MASK, n_FP, OS_1, OS_2, Fm, Gm, &
      dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2, &
      Chi_NES_1, Chi_NES_2, Eta_NES_1, Eta_NES_2, &
      Chi_Pair_1, Chi_Pair_2, Eta_Pair_1, Eta_Pair_2 )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    INTEGER,                  INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Fm, Gm
    REAL(DP),                 INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J0_1, J0_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_NES_1, Chi_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Eta_NES_1, Eta_NES_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_Pair_1, Chi_Pair_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Eta_Pair_1, Eta_Pair_2

    REAL(DP) :: Eta, Eta_T, Chi_T
    INTEGER  :: iN_E, iN_X, iFP

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

  END SUBROUTINE ComputeNeutrinoRHS_FP


  SUBROUTINE ComputeNeutrinoRHS_EmAb_FP &
    ( MASK, n_FP, OS_1, OS_2, Fm, Gm, &
      dt, Jold_1, Jold_2, Jnew_1, Jnew_2, &
      Chi_1, Chi_2, J0_1, J0_2 )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    INTEGER,                  INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Fm, Gm
    REAL(DP),                 INTENT(in)    :: dt
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Jold_1, Jold_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Jnew_1, Jnew_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Chi_1, Chi_2
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J0_1, J0_2

    REAL(DP) :: Eta
    INTEGER  :: iN_E, iN_X, iFP

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Eta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( Eta ) &
    !$ACC PRESENT( MASK, Fm, Gm, &
    !$ACC          Jold_1, Jnew_1, J0_1, Chi_1, &
    !$ACC          Jold_2, Jnew_2, J0_2, Chi_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
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


  SUBROUTINE UpdateNeutrinoRHS_FP &
    ( MASK, n_FP, OS_1, OS_2, Fm, Gm, Jnew_1, Jnew_2 )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    INTEGER,                  INTENT(in)    :: n_FP, OS_1, OS_2
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:), INTENT(inout) :: Jnew_1, Jnew_2

    INTEGER  :: iN_E, iN_X, iFP

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, Fm, Gm, Jnew_1, Jnew_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2)
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


  SUBROUTINE StoreRHS_FP &
    ( MASK, n_FP, M, Mk, Fm, Gm, F, G )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Fm, Gm
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F, G

    INTEGER  :: iN_X, iFP

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

  END SUBROUTINE StoreRHS_FP


  SUBROUTINE BuildLS_FP &
    ( MASK, n_FP, M, Mk, Fm, F, A, B )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: Fm
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: F
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: A
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: B

    INTEGER  :: iN_X, iFP, iM

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRESENT( MASK, Fm, B )
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
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( MASK, F, Fm, A )
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

  END SUBROUTINE BuildLS_FP


  SUBROUTINE SolveLS_FP &
    ( MASK, n_FP, M, Mk, A, B, Alpha )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: A
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: B
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Alpha

    REAL(DP) :: AA11, AA12, AA21, AA22, AB1, AB2, DET_AA, SUM1
    INTEGER  :: iN_X, iFP, iM

    IF ( Mk == 2 ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE &
      !$OMP PRIVATE( AA11, AB1 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC PRIVATE( AA11, AB1 ) &
      !$ACC PRESENT( MASK, A, B )
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
        !$ACC PRIVATE( AA11, AA12, AA21, AA22, AB1, AB2, DET_AA ) &
        !$ACC PRESENT( MASK, A, B )
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
        !$ACC PRIVATE( AA11, AA12, AA22, AB1, AB2, DET_AA ) &
        !$ACC PRESENT( MASK, A, B )
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

    ELSE

      DO iN_X = 1, nX_G
        IF ( MASK(iN_X) ) THEN

          CALL LinearLeastSquares &
            ( 'N', n_FP, Mk-1, 1, A(1,1,iN_X), n_FP, &
              B(1,iN_X), n_FP, TAU(1,iN_X), WORK(1,iN_X), LWORK, INFO(iN_X) )

        END IF
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( SUM1 ) &
    !$ACC PRESENT( MASK, Alpha, B )
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

  END SUBROUTINE SolveLS_FP


  SUBROUTINE ApplyAA_FP &
    ( MASK, n_FP, M, Mk, Alpha, G, Gm )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:),   INTENT(in)    :: Alpha
    REAL(DP), DIMENSION(:,:,:), INTENT(in)    :: G
    REAL(DP), DIMENSION(:,:),   INTENT(inout) :: Gm

    REAL(DP) :: SUM1
    INTEGER  :: iN_X, iFP, iM

    !CALL MatrixMatrixMultiplyBatched &
    !       ( 'N', 'N', n_FP, 1, Mk, Zero, G, n_FP, n_FP*M, &
    !         Alpha, M, M, Gm, n_FP, n_FP, nX_G )

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
          DO iM = 1, Mk
            SUM1 = SUM1 + G(iFP,iM,iN_X) * Alpha(iM,iN_X)
          END DO
          Gm(iFP,iN_X) = SUM1

        END IF
      END DO
    END DO

  END SUBROUTINE ApplyAA_FP


  SUBROUTINE ShiftRHS_FP &
    ( MASK, n_FP, M, Mk, F, G )

    LOGICAL,  DIMENSION(:),     INTENT(in)    :: MASK
    INTEGER,                    INTENT(in)    :: n_FP, M, Mk
    REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: F, G

    REAL(DP) :: FTMP(1:n_FP,1:M), GTMP(1:n_FP,1:M)
    INTEGER  :: iN_X, iFP, iM

    IF ( Mk == M ) THEN

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
          DO iM = 1, Mk-1
            DO iFP = 1, n_FP
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
            DO iFP = 1, n_FP
              F(iFP,iM,iN_X) = FTMP(iFP,iM)
              G(iFP,iM,iN_X) = GTMP(iFP,iM)
            END DO
          END DO

        END IF
      END DO

    END IF

  END SUBROUTINE ShiftRHS_FP


  SUBROUTINE ComputeJNorm &
    ( MASK, J_1, J_2, Jnorm_1, Jnorm_2 )

    LOGICAL,  DIMENSION(:),   INTENT(in)    :: MASK
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: J_1, J_2
    REAL(DP), DIMENSION(:),   INTENT(inout) :: Jnorm_1, Jnorm_2

    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRESENT( MASK, Jnorm_1, Jnorm_2, J_1, J_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD
#endif
    DO iN_X = 1, nX_G
      IF ( MASK(iN_X) ) THEN
        Jnorm_1(iN_X) = SQRT( SUM( J_1(:,iN_X)**2 ) )
        Jnorm_2(iN_X) = SQRT( SUM( J_2(:,iN_X)**2 ) )
      END IF
    END DO

  END SUBROUTINE ComputeJNorm


  SUBROUTINE CheckConvergenceInner_FP &
    ( MASK, n_FP, k_inner, OS_1, OS_2, Rtol, nIterations_Inner, Fm, Jnorm_1, Jnorm_2 )

    LOGICAL,  DIMENSION(:),   INTENT(inout) :: MASK
    INTEGER,                  INTENT(in)    :: n_FP, k_inner, OS_1, OS_2
    REAL(DP),                 INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(:),   INTENT(inout) :: nIterations_Inner
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Fm
    REAL(DP), DIMENSION(:),   INTENT(in)    :: Jnorm_1, Jnorm_2

    REAL(DP) :: Fnorm_1, Fnorm_2
    LOGICAL  :: CONVERGED
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_1, Fnorm_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_1, Fnorm_2 ) &
    !$ACC PRESENT( MASK, nIterations_Inner, Fm, Jnorm_1, Jnorm_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
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

  END SUBROUTINE CheckConvergenceInner_FP


  SUBROUTINE CheckConvergenceOuter_FP &
    ( MASK_OUTER, MASK_INNER, n_FP, iY, iE, k_outer, Fm, Rtol, nIterations_Outer )

    LOGICAL,  DIMENSION(:),   INTENT(inout) :: MASK_OUTER, MASK_INNER
    INTEGER,                  INTENT(in)    :: n_FP, iY, iE, k_outer
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Fm
    REAL(DP),                 INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(:),   INTENT(inout) :: nIterations_Outer

    REAL(DP) :: Fnorm_Y, Fnorm_E
    INTEGER  :: iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( Fnorm_Y, Fnorm_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( Fnorm_Y, Fnorm_E ) &
    !$ACC PRESENT( MASK_OUTER, MASK_INNER, nIterations_Outer, Fm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( AERR_Y, AERR_E, RERR_Y, RERR_E )
#endif
    DO iN_X = 1, nX_G
      IF ( MASK_OUTER(iN_X) ) THEN

        Fnorm_Y = ABS( Fm(iY,iN_X) )
        Fnorm_E = ABS( Fm(iE,iN_X) )

        MASK_OUTER(iN_X) = Fnorm_Y > Rtol &
                      .OR. Fnorm_E > Rtol

        MASK_INNER(iN_X) = MASK_OUTER(iN_X)

        IF ( .NOT. MASK_OUTER(iN_X) ) nIterations_Outer(iN_X) = k_outer

      END IF
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( MASK_OUTER, MASK_INNER )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( MASK_OUTER, MASK_INNER )
#endif

  END SUBROUTINE CheckConvergenceOuter_FP


  SUBROUTINE CheckConvergenceCoupled_FP &
    ( MASK, n_FP, k, iY, iE, OS_1, OS_2, Rtol, nIterations, Fm, Jnorm_1, Jnorm_2 )

    LOGICAL,  DIMENSION(:),   INTENT(inout) :: MASK
    INTEGER,                  INTENT(in)    :: n_FP, k, iY, iE, OS_1, OS_2
    REAL(DP),                 INTENT(in)    :: Rtol
    INTEGER,  DIMENSION(:),   INTENT(inout) :: nIterations
    REAL(DP), DIMENSION(:,:), INTENT(in)    :: Fm
    REAL(DP), DIMENSION(:),   INTENT(in)    :: Jnorm_1, Jnorm_2

    REAL(DP) :: Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2
    LOGICAL  :: CONVERGED
    INTEGER  :: iN_E, iN_X

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( CONVERGED, Fnorm_Y, Fnorm_E, Fnorm_1, Fnorm_2 ) &
    !$ACC PRESENT( MASK, nIterations, Fm, Jnorm_1, Jnorm_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD &
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

  END SUBROUTINE CheckConvergenceCoupled_FP


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
