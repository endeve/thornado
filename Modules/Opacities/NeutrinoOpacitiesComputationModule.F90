#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_OPACITY
#endif

MODULE NeutrinoOpacitiesComputationModule

  USE KindModule, ONLY: &
    DP, Zero, One, SqrtTiny, Pi, TwoPi, FourPi
  USE PhysicalConstantsModule, ONLY: &
    AvogadroConstantMKS, &
    SpeedOfLightCGS, &
    hbarMeVs
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    PlanckConstant, &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFE, nDOFX, &
    nE, nNodesE
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModule, ONLY: &
    Weights_Q
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InterpMat_E
  USE ReferenceElementModule_Lagrange, ONLY: &
    InterpMat
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeProtonMassFraction_TABLE, &
    ComputeNeutronMassFraction_TABLE, &
    ComputeHeavyMassFraction_TABLE, &
    ComputeHeavyMassNumber_TABLE
  USE OpacityModule_TABLE, ONLY: &
#ifdef MICROPHYSICS_WEAKLIB
    OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem, &
    EmAb_T, Iso_T, NES_T, Pair_T, Brem_T, &
    NES_AT, Pair_AT, Brem_AT, &
    use_EC_table,             &
    OS_EmAb_EC_spec_AT, EmAb_EC_spec_AT, &
    OS_EmAb_EC_rate, EmAb_EC_rate_T, & 
#endif
    LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    Ds_T, Ts_T, Ys_T, &
    C1, C2
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar, LeptonNumber

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
  PRIVATE

  PUBLIC :: ComputeNeutrinoOpacities_EC
  PUBLIC :: ComputeNeutrinoOpacities_EC_Vector
  PUBLIC :: ComputeNeutrinoOpacities_ES
  PUBLIC :: ComputeNeutrinoOpacities_ES_Vector
  PUBLIC :: ComputeNeutrinoOpacities_NES
  PUBLIC :: ComputeNeutrinoOpacityRates_NES
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_NES
  PUBLIC :: ComputeNeutrinoOpacities_Pair
  PUBLIC :: ComputeNeutrinoOpacityRates_Pair
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_Pair
  PUBLIC :: ComputeEquilibriumDistributions_Point
  PUBLIC :: ComputeEquilibriumDistributions
  PUBLIC :: ComputeEquilibriumDistributions_DG
  PUBLIC :: ComputeEquilibriumDistributionAndDerivatives
  PUBLIC :: ComputeNeutrinoOpacities_Brem
  PUBLIC :: ComputeNeutrinoOpacityRates_Brem
  PUBLIC :: ComputeNeutrinoOpacityRates_LinearCorrections_Brem
  PUBLIC :: FermiDirac
  PUBLIC :: dFermiDiracdT
  PUBLIC :: dFermiDiracdY

  REAL(DP), PARAMETER :: Log1d100 = LOG( 1.0d100 )
  REAL(DP), PARAMETER :: UnitD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitT    = Kelvin
  REAL(DP), PARAMETER :: UnitY    = One
  REAL(DP), PARAMETER :: UnitE    = MeV
  REAL(DP), PARAMETER :: UnitMn   = MeV
  REAL(DP), PARAMETER :: UnitMp   = MeV
  REAL(DP), PARAMETER :: UnitMe   = MeV
  REAL(DP), PARAMETER :: UnitEta  = One
  REAL(DP), PARAMETER :: UnitEC   = One / Centimeter
  REAL(DP), PARAMETER :: UnitES   = One / ( Centimeter * MeV**2 )
  REAL(DP), PARAMETER :: UnitNES  = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitPair = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitBrem = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: f0_Max   = One - EPSILON( One )

  REAL(DP), PARAMETER :: C_A      = -1.26d0/2.0d0 ! C_A from HR98
  REAL(DP), PARAMETER :: G_F      = 1.166d-11 ! G_F from HR98 in [MeV**-2]
  REAL(DP), PARAMETER :: nb_14    = AvogadroConstantMKS * 1d14
  REAL(DP), PARAMETER :: Brem_const = C_A**2 * G_F**2 * nb_14 &
                                    * hbarMeVs**2 * SpeedOfLightCGS**2 &
                                    / TwoPi**3 * UnitBrem

  REAL(dp), PARAMETER :: Alpha_Brem(3) = [ 1.0d0, 1.0d0, 28.d0/3.d0 ]

  INTERFACE ComputeEquilibriumDistributions_DG
    MODULE PROCEDURE ComputeEquilibriumDistributions_DG_E
    MODULE PROCEDURE ComputeEquilibriumDistributions_DG_Z
  END INTERFACE ComputeEquilibriumDistributions_DG

  INTERFACE FermiDirac
    MODULE PROCEDURE FermiDirac_Scalar
    MODULE PROCEDURE FermiDirac_Vector
  END INTERFACE

  INTERFACE dFermiDiracdT
    MODULE PROCEDURE dFermiDiracdT_Scalar
    MODULE PROCEDURE dFermiDiracdT_Vector
  END INTERFACE

  INTERFACE dFermiDiracdY
    MODULE PROCEDURE dFermiDiracdY_Scalar
    MODULE PROCEDURE dFermiDiracdY_Vector
  END INTERFACE

CONTAINS


  SUBROUTINE ComputeEquilibriumDistributions_Point( E, D, T, Y, f0, iSpecies )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distributions (Single E,D,T,Y) ---

    REAL(DP), INTENT(in)  :: E
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: f0
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: Me, Mp, Mn
    REAL(DP) :: Mnu, kT

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    kT = BoltzmannConstant * T

    IF ( iSpecies > iNuE_Bar ) THEN
      Mnu = Zero
    ELSE
      Mnu = ( ( Me + Mp ) - Mn ) * LeptonNumber(iSpecies)
    END IF

    f0 = FermiDirac( E, Mnu, kT )

  END SUBROUTINE ComputeEquilibriumDistributions_Point


  SUBROUTINE ComputeEquilibriumDistributions &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0 )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    REAL(DP), INTENT(out) :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: Me(iX_B:iX_E), Mp(iX_B:iX_E), Mn(iX_B:iX_E)
    REAL(DP) :: Mnu, kT
    INTEGER  :: iE, iS, iX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Me, Mp, Mn, f0 ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Me, Mp, Mn, f0 ) &
    !$ACC COPYIN( E, D, T, Y )
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Mnu, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Mnu, kT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Mnu, kT )
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E

      kT = BoltzmannConstant * T(iX)

      IF ( iS > iNuE_Bar ) THEN
        Mnu = Zero
      ELSE
        Mnu = ( ( Me(iX) + Mp(iX) ) - Mn(iX) ) * LeptonNumber(iS)
      END IF

      f0(iE,iS,iX) &
        = One / ( EXP( MIN( MAX( ( E(iE) - Mnu ) / kT, - Log1d100 ), + Log1d100 ) ) + One )

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Me, Mp, Mn, E, D, T, Y ) &
    !$OMP MAP( from: f0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Me, Mp, Mn, E, D, T, Y ) &
    !$ACC COPYOUT( f0 )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions


  SUBROUTINE ComputeEquilibriumDistributions_DG_E &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0 )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in) , TARGET :: E(iE_B:iE_E)
    REAL(DP), INTENT(in) , TARGET :: D(iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: T(iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: Y(iX_B:iX_E)
    REAL(DP), INTENT(out), TARGET :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP), POINTER :: E_Q(:,:), f0_Q(:,:,:,:)

    REAL(DP) :: f0_K(1:(iE_E-iE_B+1)/nDOFE,1:(iS_E-iS_B+1),1:(iX_E-iX_B+1))
    REAL(DP) :: N_K, V_K, f0_Min, f0_P, Min_K, Max_K, Theta
    INTEGER  :: iE, iS, iX, iNodeE, iP_E, nE, nS, nX

    nE = ( iE_E - iE_B + 1 ) / nDOFE
    nS = ( iS_E - iS_B + 1 )
    nX = ( iX_E - iX_B + 1 )

    ! --- Permute Data ---

    E_Q (1:nDOFE,1:nE)           => E (:)
    f0_Q(1:nDOFE,1:nE,1:nS,1:nX) => f0(:,:,:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: f0, f0_K ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( f0, f0_K ) &
    !$ACC COPYIN( E, D, T, Y )
#endif

    CALL ComputeEquilibriumDistributions &
           ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0 )

    ! --- Cell Average of Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( V_K, N_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( V_K, N_K ) &
    !$ACC PRESENT( WeightsE )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( V_K, N_K )
#endif
    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      V_K = Zero
      N_K = Zero

      DO iNodeE = 1, nDOFE
        V_K = V_K + WeightsE(iNodeE) * E_Q(iNodeE,iE)**2
        N_K = N_K + WeightsE(iNodeE) * E_Q(iNodeE,iE)**2 * f0_Q(iNodeE,iE,iS,iX)
      END DO

      f0_K(iE,iS,iX) = N_K / V_K

      IF( f0_K(iE,iS,iX) > f0_Max )THEN
        f0_K(iE,iS,iX) = f0_Max
        DO iNodeE = 1, nDOFE
          f0_Q(iNodeE,iE,iS,iX) = f0_Max
        END DO
      END IF

    END DO
    END DO
    END DO

    ! --- Limit Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta, f0_P )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( f0_Min, Min_K, Max_K, Theta, f0_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta, f0_P )
#endif
    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      IF( iE < nE ) THEN
        f0_Min = f0_K(iE+1,iS,iX)
      ELSE
        ! --- Estimate Cell Average in Outer Ghost Element ---
        f0_Min = f0_K(nE,iS,iX)**2 / f0_K(nE-1,iS,iX)
      END IF

      Max_K = f0_Max
      Min_K = f0_Min
      DO iP_E = 1, nDOFE+2
        f0_P = Zero
        DO iNodeE = 1, nDOFE
          f0_P = f0_P + InterpMat_E(iP_E,iNodeE) * f0_Q(iNodeE,iE,iS,iX)
        END DO
        Max_K = MAX( Max_K, f0_P )
        Min_K = MIN( Min_K, f0_P )
      END DO

      IF( Min_K < f0_Min .OR. Max_K > f0_Max )THEN

        Theta &
          = MIN( One, &
                 ABS((f0_Min-f0_K(iE,iS,iX))/(Min_K-f0_K(iE,iS,iX)+SqrtTiny)), &
                 ABS((f0_Max-f0_K(iE,iS,iX))/(Max_K-f0_K(iE,iS,iX)+SqrtTiny)) )

        DO iNodeE = 1, nDOFE
          f0_Q(iNodeE,iE,iS,iX) &
            = ( One - Theta ) * f0_K(iE,iS,iX) &
              +       Theta   * f0_Q(iNodeE,iE,iS,iX)
        END DO

      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: f0_K, E, D, T, Y ) &
    !$OMP MAP( from: f0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( f0_K, E, D, T, Y ) &
    !$ACC COPYOUT( f0 )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions_DG_E


  SUBROUTINE ComputeEquilibriumDistributions_DG_Z &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, SqrtGm, f0 )

    INTEGER,  INTENT(in)          :: iE_B, iE_E
    INTEGER,  INTENT(in)          :: iS_B, iS_E
    INTEGER,  INTENT(in)          :: iX_B, iX_E
    REAL(DP), INTENT(in) , TARGET :: E (iE_B:iE_E)
    REAL(DP), INTENT(in) , TARGET :: D                     (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: T                     (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: Y                     (iX_B:iX_E)
    REAL(DP), INTENT(in) , TARGET :: SqrtGm                (iX_B:iX_E)
    REAL(DP), INTENT(out), TARGET :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP), PARAMETER :: f0_Min = SqrtTiny

    INTEGER  :: nE, nX, nS
    INTEGER  :: iE, iX, iS, iP
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    REAL(DP) :: kT, Mnu, Exponent, N_K, V_K, f0_P
    REAL(DP) :: Min_K, Max_K, Theta
    REAL(DP), POINTER :: E_Q(:,:), D_Q(:,:), T_Q(:,:), Y_Q(:,:)
    REAL(DP), POINTER :: Me_Q(:,:), Mp_Q(:,:), Mn_Q(:,:), SqrtGm_Q(:,:)
    REAL(DP), TARGET  :: Me(iX_B:iX_E), Mp(iX_B:iX_E), Mn(iX_B:iX_E)
    REAL(DP) :: &
      Tau_Q(1:nDOFE*nDOFX, &
            1:(iE_E-iE_B+1)/nDOFE, &
            1:(iX_E-iX_B+1)/nDOFX)
    REAL(DP) :: &
      f0_K (1:(iE_E-iE_B+1)/nDOFE, &
            1:(iS_E-iS_B+1)      , &
            1:(iX_E-iX_B+1)/nDOFX)
    REAL(DP) :: &
      f0_Q (1:nDOFE*nDOFX, &
            1:(iE_E-iE_B+1)/nDOFE, &
            1:(iS_E-iS_B+1)      , &
            1:(iX_E-iX_B+1)/nDOFX)

    nE = (iE_E-iE_B+1)/nDOFE ! --- Number of energy  elements
    nX = (iX_E-iE_B+1)/nDOFX ! --- Number of spatial elements
    nS = (iS_E-iS_B+1)       ! --- Number of species

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    E_Q (1:nDOFE,1:nE) => E (:)
    D_Q (1:nDOFX,1:nX) => D (:)
    T_Q (1:nDOFX,1:nX) => T (:)
    Y_Q (1:nDOFX,1:nX) => Y (:)
    Me_Q(1:nDOFX,1:nX) => Me(:)
    Mp_Q(1:nDOFX,1:nX) => Mp(:)
    Mn_Q(1:nDOFX,1:nX) => Mn(:)

    SqrtGm_Q(1:nDOFX,1:nX) => SqrtGm(:)

    ! --- Compute Fermi-Dirac Distribution in Elements ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        kT = BoltzmannConstant * T_Q(iNodeX,iX)

        IF( iS > iNuE_Bar )THEN
          Mnu = Zero
        ELSE
          Mnu = (   Me_Q(iNodeX,iX) &
                  + Mp_Q(iNodeX,iX) &
                  - Mn_Q(iNodeX,iX) ) * LeptonNumber(iS)
        END IF

        Exponent &
          = MIN( MAX( ( E_Q(iNodeE,iE) - Mnu ) / kT, - Log1d100 ), Log1d100 )

        f0_Q(iNodeZ,iE,iS,iX) = One / ( EXP( Exponent ) + One )

      END DO
      END DO

    END DO
    END DO
    END DO

    ! --- Energy-Position Space Volume Jacobian ---

    DO iX = 1, nX
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        Tau_Q(iNodeZ,iE,iX) = E_Q(iNodeE,iE)**2 * SqrtGm_Q(iNodeX,iX)

      END DO
      END DO

    END DO
    END DO

    ! --- Cell Average of Equilibrium Distribution ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      V_K = Zero
      N_K = Zero

      DO iNodeZ = 1, nDOFE * nDOFX

        V_K = V_K + Weights_Q(iNodeZ) * Tau_Q(iNodeZ,iE,iX)
        N_K = N_K + Weights_Q(iNodeZ) * Tau_Q(iNodeZ,iE,iX) &
                      * f0_Q(iNodeZ,iE,iS,iX)

      END DO

      f0_K(iE,iS,iX) = MAX( MIN( N_K / V_K, f0_Max ), f0_Min )

    END DO
    END DO
    END DO

    ! --- Apply Bound-Enforcing Limiter to Equilibrium Distribution ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      Max_K = f0_Max
      Min_K = f0_Min
      DO iP = 1, SIZE( InterpMat, 1 )
        f0_P = Zero
        DO iNodeZ = 1, nDOFE * nDOFX
          f0_P = f0_P + InterpMat(iP,iNodeZ) * f0_Q(iNodeZ,iE,iS,iX)
        END DO
        Max_K = MAX( Max_K, f0_P )
        Min_K = MIN( Min_K, f0_P )
      END DO

      IF( Min_K < f0_Min .OR. Max_K > f0_Max )THEN

        Theta &
          = MIN( One, &
                 ABS((f0_Min-f0_K(iE,iS,iX))/(Min_K-f0_K(iE,iS,iX)+SqrtTiny)), &
                 ABS((f0_Max-f0_K(iE,iS,iX))/(Max_K-f0_K(iE,iS,iX)+SqrtTiny)) )

        Theta = ( One - 1.0d-2 ) * Theta

        DO iNodeZ = 1, nDOFE * nDOFX

          f0_Q(iNodeZ,iE,iS,iX) &
            = ( One - Theta ) * f0_K(iE,iS,iX) &
              +       Theta   * f0_Q(iNodeZ,iE,iS,iX)

        END DO

      END IF

    END DO
    END DO
    END DO

    ! --- Permute for Output ---

    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        f0((iE-1)*nDOFE+iNodeE,iS,(iX-1)*nDOFX+iNodeX) &
          = f0_Q(iNodeZ,iE,iS,iX)

      END DO
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeEquilibriumDistributions_DG_Z


  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0, df0dY, df0dU )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    REAL(DP), INTENT(out) :: f0   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: df0dY(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: df0dU(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: Me(iX_B:iX_E), dMedT(iX_B:iX_E), dMedY(iX_B:iX_E)
    REAL(DP) :: Mp(iX_B:iX_E), dMpdT(iX_B:iX_E), dMpdY(iX_B:iX_E)
    REAL(DP) :: Mn(iX_B:iX_E), dMndT(iX_B:iX_E), dMndY(iX_B:iX_E)
    REAL(DP) :: U (iX_B:iX_E), dUdT (iX_B:iX_E), dUdY (iX_B:iX_E), dUdD(iX_B:iX_E)
    REAL(DP) :: Mnu          , dMnudT          , dMnudY

    REAL(DP) :: kT, df0dT_Y, df0dY_T
    INTEGER  :: iE, iS, iX

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Me,  dMedT,  dMedY, &
    !$OMP             Mp,  dMpdT,  dMpdY, &
    !$OMP             Mn,  dMndT,  dMndY, &
    !$OMP             U,   dUdT,   dUdY, dUdD, &
    !$OMP             f0,  df0dY,  df0dU ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Me,  dMedT,  dMedY, &
    !$ACC         Mp,  dMpdT,  dMpdY, &
    !$ACC         Mn,  dMndT,  dMndY, &
    !$ACC         U,   dUdT,   dUdY, dUdD, &
    !$ACC         f0,  df0dY,  df0dU ) &
    !$ACC COPYIN( E, D, T, Y )
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me, dUdD, dMedT, dMedY )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp, dUdD, dMpdT, dMpdY )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn, dUdD, dMndT, dMndY )

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Y, U,  dUdD, dUdT,  dUdY  )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( kT, Mnu, dMnudT, dMnudY, df0dT_Y, df0dY_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( kT, Mnu, dMnudT, dMnudY, df0dT_Y, df0dY_T )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( kT, Mnu, dMnudT, dMnudY, df0dT_Y, df0dY_T )
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E

      kT = BoltzmannConstant * T(iX)

      IF ( iS > iNuE_Bar ) THEN
        Mnu    = Zero
        dMnudT = Zero
        dMnudY = Zero
      ELSE
        Mnu    = ( ( Me   (iX) + Mp   (iX) ) - Mn   (iX) ) * LeptonNumber(iS)
        dMnudT = ( ( dMedT(iX) + dMpdT(iX) ) - dMndT(iX) ) * LeptonNumber(iS)
        dMnudY = ( ( dMedY(iX) + dMpdY(iX) ) - dMndY(iX) ) * LeptonNumber(iS)
      END IF

      f0(iE,iS,iX) = FermiDirac   ( E(iE), Mnu, kT )
      df0dT_Y      = dFermiDiracdT( E(iE), Mnu, kT, dMnudT, T(iX) ) ! Constant T
      df0dY_T      = dFermiDiracdY( E(iE), Mnu, kT, dMnudY, T(iX) ) ! Constant Y

      df0dU(iE,iS,iX) = df0dT_Y / dUdT(iX)
      df0dY(iE,iS,iX) = df0dY_T - df0dU(iE,iS,iX) * dUdY(iX)

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Me,  dMedT,  dMedY, &
    !$OMP               Mp,  dMpdT,  dMpdY, &
    !$OMP               Mn,  dMndT,  dMndY, &
    !$OMP               U,   dUdT,   dUdY, dUdD, &
    !$OMP               E, D, T, Y ) &
    !$OMP MAP( from: f0, df0dU, df0dY )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Me,  dMedT,  dMedY, &
    !$ACC         Mp,  dMpdT,  dMpdY, &
    !$ACC         Mn,  dMndT,  dMndY, &
    !$ACC         U,   dUdT,   dUdY, dUdD, &
    !$ACC         E, D, T, Y ) &
    !$ACC COPYOUT( f0, df0dU, df0dY )
#endif

  END SUBROUTINE ComputeEquilibriumDistributionAndDerivatives


  SUBROUTINE ComputeNeutrinoOpacities_EC &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: f0(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: opEC(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: LogE_P
    REAL(DP) :: LogD_P, LogT_P
    REAL(DP) :: D_P, T_P, Y_P
    INTEGER  :: iE, iS, iX

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA       &
    !$OMP MAP( to:    E, D, T, Y ) &
    !$OMP MAP( alloc: opEC )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA              &
    !$ACC COPYIN  ( E, D, T, Y )  &
    !$ACC CREATE  ( opEC )
#endif

!do EmAb on nucleons first
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( LogE_P, LogD_P, LogT_P, Y_P ) 
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( LogE_P, LogD_P, LogT_P, Y_P ) 
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E

      IF ( iS > iNuE_Bar ) THEN

        opEC(iE,iS,iX) = Zero

      ELSE

        LogE_P = LOG10( E(iE) / UnitE )

        LogD_P = LOG10( D(iX) / UnitD )
        LogT_P = LOG10( T(iX) / UnitT )
        Y_P    =        Y(iX) / UnitY

        CALL LogInterpolateSingleVariable_4D_Custom_Point &
               ( LogE_P , LogD_P , LogT_P , Y_P , &
                 LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                 OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(iE,iS,iX) )

        opEC(iE,iS,iX) = opEC(iE,iS,iX) * UnitEC

      END IF

    END DO
    END DO
    END DO

!now add contribution from EC on heavy nuclei
  IF(use_EC_table .gt. 0) THEN

  EC_Table: BLOCK

  INTEGER :: iE1, iNodeE1

  REAL(dp) :: rate, loctot
  REAL(dp) :: tmev
  REAL(dp) :: spec_nodes(iE_B:iE_E)

  REAL(dp) :: Xnuc
  REAL(dp), PARAMETER :: coeff = 2.0d0 * (Pi*SpeedOfLightCGS)**2 &
                               * hbarMeVs**3
  REAL(dp), PARAMETER :: dmnp  = 1.29333d+00   ! mn - mp
  REAL(dp), PARAMETER :: kmev  = 8.61733d-11   ! Boltzmann's constant [MeV K^{-1}]
  REAL(dp) :: Me(iX_B:iX_E), Mn(iX_B:iX_E), Mp(iX_B:iX_E)
  REAL(dp) :: Xh(iX_B:iX_E), Ah(iX_B:iX_E)
  REAL(dp) :: EC_rate(iX_B:iX_E)
  REAL(dp) :: opECT(iE_B:iE_E,iX_B:iX_E)

  REAL(dp) :: E_node
  REAL(dp) :: CenterE(nE), WidthE(nE), NodesE(nE * nNodesE)

  INTEGER  :: iD, iT, iY
  REAL(dp) :: dD, dT, dY

  CenterE(:) = MeshE % Center(:)
  WidthE(:)  = MeshE % Width(:)
  NodesE(:)  = MeshE % Nodes(:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Me, Mn, Mp, Xh, Ah, EC_rate, opECT, spec_nodes )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Me, Mn, Mp, Xh, Ah, EC_rate, opECT, spec_nodes )
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    ! --- Compute heavy mass fraction and number ---

    CALL ComputeHeavyMassFraction_TABLE &
           ( D, T, Y, Xh )

    CALL ComputeHeavyMassNumber_TABLE &
           ( D, T, Y, Ah )

    ! --- Interpolate EC rate
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( LogD_P, LogT_P, D_P, T_P, Y_P ) 
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR                &
    !$ACC PRIVATE( LogD_P, LogT_P, D_P, T_P, Y_P ) 
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO                              &
    !$OMP PRIVATE( D_P, T_P, Y_P )
#endif
    DO iX = iX_B, iX_E

        !LogD_P = LOG10( D(iX) / UnitD )
        !LogT_P = LOG10( T(iX) / UnitT )
        !Y_P    =        Y(iX) / UnitY
        D_P = D(iX) / UnitD
        T_P = T(iX) / UnitT
        Y_P = Y(iX) / UnitY

        CALL LogInterpolateSingleVariable_3D_Custom_Point &
               ( D_P,  T_P,  Y_P,  &
                 Ds_T, Ts_T, Ys_T, &
                 OS_EmAb_EC_rate(1), EmAb_EC_rate_T, EC_rate(iX)) 

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE                                 &
    !$OMP PRIVATE( D_P, T_P, Y_P,                 &
    !$OMP          loctot, spec_nodes, tmev, Xnuc )               &
    !$OMP MAP    ( to: CenterE, WidthE, NodesE )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG                                      &
    !$ACC PRIVATE( iE, iX, D_P, T_P, Y_P,         &
    !$ACC          loctot, spec_nodes, tmev, Xnuc )               &
    !$ACC COPYIN ( CenterE, WidthE, NodesE, iD, iT, iY)           &
    !$ACC PRESENT( Me, Mn, Mp, Xh, Ah, EC_rate, OS_EmAb_EC_rate,  &
    !$ACC          OS_EmAb_EC_spec_AT, EmAb_EC_spec_AT, WeightsE, &
    !$ACC          f0) 
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO                                             &
    !$OMP PRIVATE( D_P, T_P, Y_P, tmev, E_node,   &
    !$OMP          loctot, spec_nodes, iE1, iNodeE1, Xnuc )
#endif
    DO iX = iX_B, iX_E

      IF((D(iX) / UnitD) > 1.0d13 .or. Ah(iX) < 40.0d0 ) THEN

        DO iE = iE_B, iE_E
          opECT(iE,iX) = 0.0d0
        END DO

      ELSE

        spec_nodes = 0.0d0
        loctot     = 0.0d0

        !LogD_P     = LOG10( D(iX) / UnitD )
        !LogT_P     = LOG10( T(iX) / UnitT )
        D_P     = D(iX) / UnitD
        T_P     = T(iX) / UnitT
        Y_P     = Y(iX) / UnitY

        tmev    = T(iX) / UnitT * kmev
        Xnuc    = Xh(iX) / (Ah(iX) / AvogadroConstantMKS) * D(iX) / UnitD 

        if(iX == 1) then
          !CALL GetIndexAndDelta_Log( LogD_P, LogDs_T, iD, dD )
          !CALL GetIndexAndDelta_Log( LogT_P, LogTs_T, iT, dT )
          CALL GetIndexAndDelta_Log( D_P, Ds_T, iD, dD )
          CALL GetIndexAndDelta_Log( T_P, Ts_T, iT, dT )
          CALL GetIndexAndDelta_Lin( Y_P, Ys_T, iY, dY )
          print *, iD, iT, iY
          print *, 'rho(81)', Ds_T(iD) !D(iX) / UnitD
          print *, 'T(28)',   Ts_T(iT) !T(iX) / UnitT
          print *, 'T(29)',   Ts_T(iT+1) !T(iX) / UnitT
          print *, 'T input', T(iX) / UnitT
          print *, 'Ye(41)',  Ys_T(iY) !!Y(iX) / UnitY
          print *, 'tmev', tmev
          print *, 'Xnuc', Xnuc
          print *, 'Xh', Xh(iX)
          print *, 'Ah', Ah(iX)
          print *, 'Mn', Mn(iX) / UnitMn
          print *, 'Mp', Mp(iX) / UnitMp
          print *, 'Me', Me(iX) / UnitMe
          print *, 'coeff', coeff
          print *, 'EC_rate', EC_rate(iX)
          print *, 'EC_rate_T', 10.0d0**EmAb_EC_rate_T(iD,iT+1,iY)
          print *, 'OS_EmAb_EC_rate', OS_EmAb_EC_rate(1)
          print *, 'rho(76)', LogDs_T(iD)
          print *, 'log10(rho)', log10(D(iX) / UnitD)
          print *, 'rho(77)', LogDs_T(iD+1)
          print *, 'T(30)', LogTs_T(iT)
          print *, 'log10(T)', log10(T(iX) / UnitT)
          print *, 'T(31)', LogTs_T(iT+1)
          print *, 'Ye(40)', Ys_T(iY)
          print *, 'Ye', Y(iX) / UnitY
          print *, 'Ye(41)', Ys_T(iY+1)
        endif

#if defined(THORNADO_OMP_OL)
    !$OMP PARALLEL DO SIMD            &
    !$OMP REDUCTION(+:loctot)         &
    !$OMP PRIVATE( iE1, iNodeE1 )
#elif defined(THORNADO_OACC)
    !$ACC LOOP VECTOR                 &
    !$ACC REDUCTION(+:loctot)         &
    !$ACC PRIVATE( iE1, iNodeE1 ) 
#endif
        DO iE = iE_B, iE_E

          iE1     = MOD( (iE-1) / nNodesE, nE      ) + 1
          iNodeE1 = MOD( (iE-1)          , nNodesE ) + 1

          !CALL LogInterpolateSingleVariable_3D_Custom_Point &
          !       ( LogD_P,  LogT_P,  Y_P,  &
          !         LogDs_T, LogTs_T, Ys_T, &
          !         OS_EmAb_EC_spec_AT, EmAb_EC_spec_AT(:,:,:,iE), spec_nodes(iE)) 
          CALL LogInterpolateSingleVariable_3D_Custom_Point &
                 ( D_P,  T_P,  Y_P,  &
                   Ds_T, Ts_T, Ys_T, &
                   OS_EmAb_EC_spec_AT, EmAb_EC_spec_AT(:,:,:,iE), spec_nodes(iE)) 

          if(spec_nodes(iE) < 1.0d-80) spec_nodes(iE) = Zero

          loctot = loctot + WidthE(iE1)/UnitE * WeightsE(iNodeE1) &
                          * spec_nodes(iE)

          if(iX==1) then
            if(iE==1) print *, 'EmAb_EC_spec_AT(81,29,41)'
            print *, iE, 10.0d0**EmAb_EC_spec_AT(81,28,41,iE) * EC_rate(iX)
          endif
          if(iX==1) then
            if(iE==1) print *, 'EmAb_EC_spec_AT(81+1,29+1,41+1)'
            print *, iE, 10.0d0**EmAb_EC_spec_AT(81+1,28+1,41+1,iE) * EC_rate(iX)
          endif
          if(iX==1) then
            if(iE==1) print *, 'spec_rate_nodes'
            print *, iE, spec_nodes(iE) * EC_rate(iX)
          endif

        END DO

        if(ix==1) then
        print *, 'loctot thornado', loctot, loctot * EC_rate(iX) * 262481508.18084475d0 / (D(iX) / UnitD)
        print *, 'coeff', Xnuc * coeff, 262481508.18084475d0

        endif

#if defined(THORNADO_OMP_OL)
    !$OMP PARALLEL DO SIMD                    &
    !$OMP PRIVATE( iE1, iNodeE1, E_node )
#elif defined(THORNADO_OACC)
    !$ACC LOOP VECTOR                         &
    !$ACC PRIVATE( iE1, iNodeE1, E_node ) 
#endif
        DO iE = iE_B, iE_E
          iE1     = MOD( (iE-1) / nNodesE, nE      ) + 1
          iNodeE1 = MOD( (iE-1)          , nNodesE ) + 1
          E_node  = NodeCoordinate( CenterE(iE1), WidthE(iE1), NodesE(iNodeE1) ) / UnitE

          spec_nodes(iE) = spec_nodes(iE) / loctot * EC_rate(iX) / (E_node)**2
          !spec_nodes(iE) = spec_nodes(iE) * EC_rate(iX) / (E_node)**2
          !spec_nodes(iE) = spec_nodes(iE) / loctot / (E_node)**2

          !opECT(iE,iX) = Xnuc * coeff * spec_nodes(iE) * (1.0d0 + &
          !EXP((E_node + dmnp + Mn(iX)/UnitMn - Mp(iX)/UnitMp - Me(iX)/UnitMe)/tmev))

          !opECT(iE,iX) = Xnuc * coeff * spec_nodes(iE) * (1.0d0 - f0(iE,1,iX)) / f0(iE,1,iX)
          !opECT(iE,iX) = Xnuc * coeff * spec_nodes(iE) / f0(iE,1,iX)
          opECT(iE,iX) = Xnuc * coeff * spec_nodes(iE) !* SpeedOfLightCGS !* 4.0d0 / (PlanckConstant*SpeedOfLightCGS)**3

          !opECT(iE,iX) = Xnuc * coeff * (10.0d0**EmAb_EC_spec_AT(95,33,40,iE) - OS_EmAb_EC_spec_AT) / f0(iE,1,iX)

          if(iX==1) then
          if(iE==1) print *, 'E_node, spec_nodes'
          print *, iE, E_node, spec_nodes(iE)
          endif
          if(iX==1) then
          if(iE==1) print *, 'E, opECT'
          print *, iE, E_node, opECT(iE,iX) !, f0(iE,1,iX)
          endif
          if(iX==1) then
          if(iE==1) print *, '(1-f0)/f0, 1+exp((enu+dmnp+mn-mp-me)/tmev)'
          print *, iE, (1.0d0 - f0(iE,1,iX)) / f0(iE,1,iX), ( &
          EXP((E_node + dmnp + Mn(iX)/UnitMn - Mp(iX)/UnitMp - Me(iX)/UnitMe)/tmev))
   if(iE==1) print *, 'loctot', loctot 
   if(iE==1) print *, 'EC_rate', EC_rate(iX)
          endif
        END DO

      ENDIF
    END DO

!now add Chi from EC table on heavy nuclei to Chi from EmAb on nucleons
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) 
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) 
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) 
#endif
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E

        !opEC(iE,1,iX) = opEC(iE,1,iX) + opECT(iE,iX)
        opEC(iE,1,iX) = opECT(iE,iX) * UnitEC

    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA                                          &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P,                &
    !$OMP               Me, Mn, Mp, Xh, Ah, EC_rate, opECT, spec_nodes )     
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA                                                 &
    !$ACC DELETE  ( LogE_P, LogD_P, LogT_P, Y_P,                    &
    !$ACC           Me, Mn, Mp, Xh, Ah, EC_rate, opECT, spec_nodes) 
#endif

  END BLOCK EC_Table

  ENDIF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA          &
    !$OMP MAP (release: E, D. T, Y) &
    !$OMP MAP (from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA                 &
    !$ACC DELETE (E, D, T, Y )      &
    !$ACC COPYOUT ( opEC )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( opEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E
      opEC(iE,iS,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC

  
  SUBROUTINE ComputeNeutrinoOpacities_EC_Vector &
    ( iP_B, iP_E, iS_B, iS_E, E, D, T, Y, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: D(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: T(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: Y(iP_B:iP_E)
    REAL(DP), INTENT(out) :: opEC(iP_B:iP_E,iS_B:iS_E)

    REAL(DP) :: LogE_P(iP_B:iP_E)
    REAL(DP) :: LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E)
    INTEGER  :: iP, iS

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P, opEC ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P, opEC ) &
    !$ACC COPYIN( E, D, T, Y )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Y_P   (iP) =        Y(iP) / UnitY
      LogE_P(iP) = LOG10( E(iP) / UnitE )
    END DO

    DO iS = iS_B, iS_E

      CALL LogInterpolateSingleVariable_4D_Custom &
             ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(:,iS) )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS = iS_B, iS_E
    DO iP = iP_B, iP_E
      opEC(iP,iS) = opEC(iP,iS) * UnitEC
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Y ) &
    !$OMP MAP( from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Y ) &
    !$ACC COPYOUT( opEC )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP MAP( from: opEC )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC COPYOUT( opEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS = iS_B, iS_E
    DO iP = iP_B, iP_E
      opEC(iP,iS) = Zero
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Vector


  SUBROUTINE ComputeNeutrinoOpacities_ES &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iMoment, opES )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: LogE_P(iE_B:iE_E)
    REAL(DP) :: LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    INTEGER  :: iX, iE
    REAL(DP) :: Phi0(iE_B:iE_E,iX_B:iX_E), Phi1(iE_B:iE_E,iX_B:iX_E)

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P, opES, Phi0, Phi1 ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P, opES, Phi0, Phi1 ) &
    !$ACC COPYIN( E, D, T, Y )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    !CALL LogInterpolateSingleVariable_1D3D_Custom &
    !       ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
    !         OS_Iso(1,iMoment), Iso_T(:,:,:,:,iMoment,1), opES )
    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(1,1), Iso_T(:,:,:,:,1,1), Phi0 )
    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(1,2), Iso_T(:,:,:,:,2,1), Phi1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E
      !opES(iE,iX) = opES(iE,iX) * UnitES
      opES(iE,iX) = (Phi0(iE,iX) - Phi1(iE,iX) / 3.0d0) * UnitES
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Y, Phi0, Phi1 ) &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Y, Phi0, Phi1 ) &
    !$ACC COPYOUT( opES )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC COPYOUT( opES )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E
      opES(iE,iX) = Zero
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES

  
  SUBROUTINE ComputeNeutrinoOpacities_ES_Vector &
    ( iP_B, iP_E, E, D, T, Y, iMoment, opES )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    REAL(DP), INTENT(in)  :: E(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: D(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: T(iP_B:iP_E)
    REAL(DP), INTENT(in)  :: Y(iP_B:iP_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(iP_B:iP_E)

    REAL(DP) :: LogE_P(iP_B:iP_E)
    REAL(DP) :: LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E)
    INTEGER  :: iP

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P, opES ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P, opES ) &
    !$ACC COPYIN( E, D, T, Y )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Y_P(iP) = Y(iP) / UnitY
      LogE_P(iP) = LOG10( E(iP) / UnitE )
    END DO

    CALL LogInterpolateSingleVariable_4D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(1,iMoment), Iso_T(:,:,:,:,iMoment,1), opES )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      opES(iP) = opES(iP) * UnitES
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Y ) &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P, E, D, T, Y ) &
    !$ACC COPYOUT( opES )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: opES )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( opES )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      opES(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Vector


  SUBROUTINE ComputeNeutrinoOpacities_NES &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, iMoment, H_I, H_II )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: H_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: H_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    INTEGER  :: iX, iE1, iE2, iH_I, iH_II

#ifdef MICROPHYSICS_WEAKLIB

    iH_I  = ( iMoment - 1 ) * 2 + 1
    iH_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogT_P, LogEta_P, H_I, H_II ) &
    !$OMP MAP( to: D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogT_P, LogEta_P, H_I, H_II ) &
    !$ACC COPYIN( D, T, Y )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate HI  ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH_I), NES_AT(:,:,:,:,iH_I,1), H_I )

    ! --- Interpolate HII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH_II), NES_AT(:,:,:,:,iH_II,1), H_II )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogT_P, LogEta_P, D, T, Y ) &
    !$OMP MAP( from: H_I, H_II )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogT_P, LogEta_P, D, T, Y ) &
    !$ACC COPYOUT( H_I, H_II )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: H_I, H_II )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( H_I, H_II )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      H_I (iE1,iE2,iX) = Zero
      H_II(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES


  SUBROUTINE ComputeNeutrinoOpacityRates_NES &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, J, J0, H_I, H_II, Eta, Chi )

    ! --- Neutrino-Electron Scattering Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2  (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J0  (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_Out, Phi_In
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iE1, iE2, iS, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, DetBal, Phi_In, Phi_Out ) &
    !$OMP MAP( to: H_I, H_II, W2, J, J0 ) &
    !$OMP MAP( from: Eta, Chi )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( SUM1, SUM2, DetBal, Phi_In, Phi_Out ) &
    !$ACC COPYIN( H_I, H_II, W2, J, J0 ) &
    !$ACC COPYOUT( Eta, Chi ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, DetBal, Phi_In, Phi_Out )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iS,iX) * ( One - J0(iE1,iS,iX) ) ) &
                 / ( J0(iE1,iS,iX) * ( One - J0(iE2,iS,iX) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_Out = (   C1(iS) * H_I (iE1,iE2,iX) &
                      + C2(iS) * H_II(iE1,iE2,iX) ) * UnitNES
          Phi_In  = Phi_Out * DetBal
        ELSE
          Phi_In  = (   C1(iS) * H_I (iE2,iE1,iX) &
                      + C2(iS) * H_II(iE2,iE1,iX) ) * UnitNES
          Phi_Out = Phi_In / DetBal
        END IF

        SUM1 = SUM1 + Phi_In  * W2(iE1) * J(iE1,iS,iX)
        SUM2 = SUM2 + Phi_Out * W2(iE1) * ( One - J(iE1,iS,iX) )

      END DO

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_NES


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NES &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, H_1, H_2, H_3, J0, &
      H_I_1, H_II_1, L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 )

    ! --- Neutrino-Electron Scattering Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J0     (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_I_1  (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_II_1 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_1 (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_2 (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_In_3 (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Out_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Out, Phi_1_In
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_In, Phi_1_Out ) &
    !$OMP MAP( to: H_I_1, H_II_1, W2, H_1, H_2, H_3, J0 ) &
    !$OMP MAP( from: L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_In, Phi_1_Out ) &
    !$ACC COPYIN( H_I_1, H_II_1, W2, H_1, H_2, H_3, J0 ) &
    !$ACC COPYOUT( L_In_1, L_In_2, L_In_3, L_Out_1, L_Out_2, L_Out_3 ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_In, Phi_1_Out )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iS,iX) * ( One - J0(iE1,iS,iX) ) ) &
                 / ( J0(iE1,iS,iX) * ( One - J0(iE2,iS,iX) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_1_Out = (   C1(iS) * H_I_1 (iE1,iE2,iX) &
                        + C2(iS) * H_II_1(iE1,iE2,iX) ) * UnitNES
          Phi_1_In  = Phi_1_Out * DetBal
        ELSE
          Phi_1_In  = (   C1(iS) * H_I_1 (iE2,iE1,iX) &
                        + C2(iS) * H_II_1(iE2,iE1,iX) ) * UnitNES
          Phi_1_Out = Phi_1_In / DetBal
        END IF

        SUM1 = SUM1 + Phi_1_In  * W2(iE1) * H_1(iE1,iS,iX)
        SUM2 = SUM2 + Phi_1_Out * W2(iE1) * H_1(iE1,iS,iX)
        SUM3 = SUM3 + Phi_1_In  * W2(iE1) * H_2(iE1,iS,iX)
        SUM4 = SUM4 + Phi_1_Out * W2(iE1) * H_2(iE1,iS,iX)
        SUM5 = SUM5 + Phi_1_In  * W2(iE1) * H_3(iE1,iS,iX)
        SUM6 = SUM6 + Phi_1_Out * W2(iE1) * H_3(iE1,iS,iX)

      END DO

      L_In_1 (iE2,iS,iX) = SUM1
      L_Out_1(iE2,iS,iX) = SUM2
      L_In_2 (iE2,iS,iX) = SUM3
      L_Out_2(iE2,iS,iX) = SUM4
      L_In_3 (iE2,iS,iX) = SUM5
      L_Out_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_NES


  SUBROUTINE ComputeNeutrinoOpacities_Pair &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, iMoment, J_I, J_II )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: J_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: J_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    INTEGER  :: iX, iE1, iE2, iJ_I, iJ_II

#ifdef MICROPHYSICS_WEAKLIB

    iJ_I  = ( iMoment - 1 ) * 2 + 1
    iJ_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: LogT_P, LogEta_P, J_I, J_II ) &
    !$OMP MAP( to: D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( LogT_P, LogEta_P, J_I, J_II ) &
    !$ACC COPYIN( D, T, Y )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      LogT_P  (iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) &
                            / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate JI  ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_I), Pair_AT(:,:,:,:,iJ_I,1), J_I )

    ! --- Interpolate JII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_II), Pair_AT(:,:,:,:,iJ_II,1), J_II )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogT_P, LogEta_P, D, T, Y ) &
    !$OMP MAP( from: J_I, J_II )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( LogT_P, LogEta_P, D, T, Y ) &
    !$ACC COPYOUT( J_I, J_II )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: J_I, J_II)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( J_I, J_II)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      J_I (iE1,iE2,iX) = Zero
      J_II(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair


  SUBROUTINE ComputeNeutrinoOpacityRates_Pair &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, J, J0, J_I, J_II, Eta, Chi )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2  (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J   (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J0  (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_I (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_II(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE1, iE2, iS, iS_A

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Pro, Phi_0_Ann ) &
    !$OMP MAP( to: J_I, J_II, W2, J, J0 ) &
    !$OMP MAP( from: Eta, Chi )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Pro, Phi_0_Ann ) &
    !$ACC COPYIN( J_I, J_II, W2, J, J0 ) &
    !$ACC COPYOUT( Eta, Chi ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Pro, Phi_0_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iS,iX) * J0(iE1,iS_A,iX) ) &
                 / ( ( One - J0(iE2,iS,iX) ) * ( One - J0(iE1,iS_A,iX) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_0_Ann = (   C1(iS) * J_I (iE1,iE2,iX) &
                        + C2(iS) * J_II(iE1,iE2,iX) ) * UnitPair
        ELSE
          Phi_0_Ann = (   C1(iS) * J_II(iE2,iE1,iX) &
                        + C2(iS) * J_I (iE2,iE1,iX) ) * UnitPair
        END IF
        Phi_0_Pro = Phi_0_Ann * DetBal

        SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iS_A,iX) )
        SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iS_A,iX)

      END DO

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_Pair


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Pair &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, H_1, H_2, H_3, J0, &
      J_I_1, J_II_1, L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )

    ! --- e^+ e^- Pair Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J0     (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_I_1  (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J_II_1 (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Ann, Phi_1_Pro
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iS_A, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann ) &
    !$OMP MAP( to: J_I_1, J_II_1, W2, H_1, H_2, H_3, J0 ) &
    !$OMP MAP( from: L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_Pro, Phi_1_Ann ) &
    !$ACC COPYIN( J_I_1, J_II_1, W2, H_1, H_2, H_3, J0 ) &
    !$ACC COPYOUT( L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 ) &
    !$ACC PRESENT( C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iS,iX) * J0(iE1,iS_A,iX) ) &
                 / ( ( One - J0(iE2,iS,iX) ) * ( One - J0(iE1,iS_A,iX) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_1_Ann = (   C1(iS) * J_I_1 (iE1,iE2,iX) &
                        + C2(iS) * J_II_1(iE1,iE2,iX) ) * UnitPair
        ELSE
          Phi_1_Ann = (   C1(iS) * J_II_1(iE2,iE1,iX) &
                        + C2(iS) * J_I_1 (iE2,iE1,iX) ) * UnitPair
        END IF
        Phi_1_Pro = Phi_1_Ann * DetBal

        SUM1 = SUM1 + Phi_1_Pro * W2(iE1) * H_1(iE1,iS_A,iX)
        SUM2 = SUM2 + Phi_1_Ann * W2(iE1) * H_1(iE1,iS_A,iX)
        SUM3 = SUM3 + Phi_1_Pro * W2(iE1) * H_2(iE1,iS_A,iX)
        SUM4 = SUM4 + Phi_1_Ann * W2(iE1) * H_2(iE1,iS_A,iX)
        SUM5 = SUM5 + Phi_1_Pro * W2(iE1) * H_3(iE1,iS_A,iX)
        SUM6 = SUM6 + Phi_1_Ann * W2(iE1) * H_3(iE1,iS_A,iX)

      END DO

      L_Pro_1(iE2,iS,iX) = SUM1
      L_Ann_1(iE2,iS,iX) = SUM2
      L_Pro_2(iE2,iS,iX) = SUM3
      L_Ann_2(iE2,iS,iX) = SUM4
      L_Pro_3(iE2,iS,iX) = SUM5
      L_Ann_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Pair


  SUBROUTINE ComputeNeutrinoOpacities_Brem &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, S_sigma )

    ! --- Brem Opacities (Multiple D,T) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    REAL(DP), INTENT(out) :: S_Sigma(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    INTEGER  :: iX, iE1, iE2
    REAL(DP) :: Xp(iX_B:iX_E), Xn(iX_B:iX_E) !Proton and neutron mass fractions
    REAL(DP) :: LogT_P(iX_B:iX_E)
    REAL(DP) :: LogDX_P(3,iX_B:iX_E)

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Xp, Xn, LogT_P, LogDX_P, S_sigma ) &
    !$OMP MAP( to: D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Xp, Xn, LogT_P, LogDX_P, S_sigma ) &
    !$ACC COPYIN( D, T, Y )
#endif

    ! --- Compute proton and neutron fractions ---

    CALL ComputeProtonMassFraction_TABLE &
           ( D, T, Y, Xp )

    CALL ComputeNeutronMassFraction_TABLE &
           ( D, T, Y, Xn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      LogDX_P(1,iX) = LOG10( D(iX) * Xp(iX) / UnitD )
      LogDX_P(2,iX) = LOG10( D(iX) * Xn(iX) / UnitD )
      LogDX_P(3,iX) = LOG10( D(iX) * SQRT(ABS(Xp(iX)*Xn(iX))) / UnitD )

      LogT_P(iX)    = LOG10( T(iX) / UnitT )    

    END DO

    CALL SumLogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogDX_P, LogT_P, LogDs_T, LogTs_T, Alpha_Brem, &
             OS_Brem(1,1), Brem_AT(:,:,:,:,1,1), S_sigma )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Xp, Xn, LogT_P, LogDX_P, D, T, Y ) &
    !$OMP MAP( from: S_sigma )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Xp, Xn, LogT_P, LogDX_P, D, T, Y ) &
    !$ACC COPYOUT( S_sigma )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: S_sigma )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( S_sigma )
#endif
    DO iX  = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      S_sigma(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Brem

  
  SUBROUTINE ComputeNeutrinoOpacityRates_Brem &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, J, J0, S_sigma, Eta, Chi )

    ! --- Brem Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: J      (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J0     (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: S_Sigma(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Eta    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Chi    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE1, iE2, iS, iS_A

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro ) &
    !$OMP MAP( to: S_sigma, W2, J, J0 ) &
    !$OMP MAP( from: Eta, Chi )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro ) &
    !$ACC COPYIN( S_sigma, W2, J, J0 ) &
    !$ACC COPYOUT( Eta, Chi )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal = ( J0(iE2,iS,iX) * J0(iE1,iS_A,iX) ) &
                 / ( ( One - J0(iE2,iS,iX) ) * ( One - J0(iE1,iS_A,iX) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_0_Ann = S_Sigma(iE1,iE2,iX) * 3.0d0 * Brem_const
        ELSE
          Phi_0_Ann = S_Sigma(iE2,iE1,iX) * 3.0d0 * Brem_const
        END IF
        Phi_0_Pro = Phi_0_Ann * DetBal

        SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iS_A,iX) )
        SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iS_A,iX)

      END DO

      Eta(iE2,iS,iX) = SUM1
      Chi(iE2,iS,iX) = SUM1 + SUM2

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_Brem


  SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Brem &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, H_1, H_2, H_3, J0, &
      S_Sigma, L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )

    ! --- N N Brem Rates (Linear Corrections) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (iE_B:iE_E)
    REAL(DP), INTENT(in)  :: H_1    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_2    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: H_3    (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: J0     (iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(in)  :: S_Sigma(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Pro_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_1(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_2(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: L_Ann_3(iE_B:iE_E,iS_B:iS_E,iX_B:iX_E)

    REAL(DP) :: DetBal, Phi_1_Ann, Phi_1_Pro
    REAL(DP) :: SUM1, SUM2, SUM3, SUM4, SUM5, SUM6
    INTEGER  :: iE1, iE2, iS, iS_A, iX

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann ) &
    !$OMP MAP( to: S_Sigma, W2, H_1, H_2, H_3, J0 ) &
    !$OMP MAP( from: L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$ACC          Phi_1_Pro, Phi_1_Ann ) &
    !$ACC COPYIN( S_Sigma, W2, H_1, H_2, H_3, J0 ) &
    !$ACC COPYOUT( L_Pro_1, L_Pro_2, L_Pro_3, L_Ann_1, L_Ann_2, L_Ann_3 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, DetBal, &
    !$OMP          Phi_1_Pro, Phi_1_Ann )
#endif
    DO iX  = iX_B, iX_E
    DO iS  = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero
      SUM3 = Zero
      SUM4 = Zero
      SUM5 = Zero
      SUM6 = Zero

      DO iE1 = iE_B, iE_E

        DetBal = ( J0(iE2,iS,iX) * J0(iE1,iS_A,iX) ) &
                 / ( ( One - J0(iE2,iS,iX) ) * ( One - J0(iE1,iS_A,iX) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_1_Ann = - S_Sigma(iE1,iE2,iX) * Brem_const
        ELSE
          Phi_1_Ann = - S_Sigma(iE2,iE1,iX) * Brem_const
        END IF
        Phi_1_Pro = Phi_1_Ann * DetBal

        SUM1 = SUM1 + Phi_1_Pro * W2(iE1) * H_1(iE1,iS_A,iX)
        SUM2 = SUM2 + Phi_1_Ann * W2(iE1) * H_1(iE1,iS_A,iX)
        SUM3 = SUM3 + Phi_1_Pro * W2(iE1) * H_2(iE1,iS_A,iX)
        SUM4 = SUM4 + Phi_1_Ann * W2(iE1) * H_2(iE1,iS_A,iX)
        SUM5 = SUM5 + Phi_1_Pro * W2(iE1) * H_3(iE1,iS_A,iX)
        SUM6 = SUM6 + Phi_1_Ann * W2(iE1) * H_3(iE1,iS_A,iX)

      END DO

      L_Pro_1(iE2,iS,iX) = SUM1
      L_Ann_1(iE2,iS,iX) = SUM2
      L_Pro_2(iE2,iS,iX) = SUM3
      L_Ann_2(iE2,iS,iX) = SUM4
      L_Pro_3(iE2,iS,iX) = SUM5
      L_Ann_3(iE2,iS,iX) = SUM6

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_LinearCorrections_Brem

  
  FUNCTION FermiDirac_Scalar( E, Mu, kT ) RESULT( FermiDirac )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT
    REAL(DP) :: FermiDirac

    REAL(DP) :: Exponent

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac_Scalar


  FUNCTION FermiDirac_Vector( E, Mu, kT ) RESULT( FermiDirac )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT
    REAL(DP) :: FermiDirac(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac_Vector


  FUNCTION dFermiDiracdT_Scalar( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudT, T
    REAL(DP) :: dFermiDiracdT

    REAL(DP) :: Exponent, FD

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdT &
      = ( FD * EXP( Exponent ) ) * FD * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT_Scalar


  FUNCTION dFermiDiracdT_Vector( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudT, T
    REAL(DP) :: dFermiDiracdT(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E)), FD(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdT &
      = ( FD * EXP( Exponent ) ) * FD * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT_Vector


  FUNCTION dFermiDiracdY_Scalar( E, Mu, kT, dMudY, T ) RESULT( dFermiDiracdY )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY

    REAL(DP) :: Exponent, FD

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Scalar


  FUNCTION dFermiDiracdY_Vector( E, Mu, kT, dMudY, T ) RESULT( dFermiDiracdY )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E)), FD(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Vector


END MODULE NeutrinoOpacitiesComputationModule
