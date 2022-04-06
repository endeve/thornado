#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_OPACITY
#endif

MODULE NeutrinoOpacitiesComputationModule

  USE KindModule, ONLY: &
    DP, Zero, One, SqrtTiny, TwoPi
  USE PhysicalConstantsModule, ONLY: &
    AvogadroConstantMKS, &
    SpeedOfLightCGS
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFE, nDOFX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InterpMat_E
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeProtonMassFraction_TABLE, &
    ComputeNeutronMassFraction_TABLE
  USE OpacityModule_TABLE, ONLY: &
#ifdef MICROPHYSICS_WEAKLIB
    OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem, &
    EmAb_T, Iso_T, NES_T, Pair_T, Brem_T, &
    NES_AT, Pair_AT, Brem_AT, &
#endif
    LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    C1, C2
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar, LeptonNumber

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules ---

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_4D_Custom, &
    LogInterpolateSingleVariable_4D_Custom_Point, &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned, &
    SumLogInterpolateSingleVariable_2D2D_Custom_Aligned

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
  PUBLIC :: ComputeNeutrinoOpacities_Pair
  PUBLIC :: ComputeNeutrinoOpacityRates_Pair
  PUBLIC :: ComputeEquilibriumDistributions_Point
  PUBLIC :: ComputeEquilibriumDistributions
  PUBLIC :: ComputeEquilibriumDistributions_DG
  PUBLIC :: ComputeEquilibriumDistributionAndDerivatives
  PUBLIC :: ComputeNeutrinoOpacities_Brem
  PUBLIC :: ComputeNeutrinoOpacityRates_Brem
  PUBLIC :: FermiDirac
  PUBLIC :: dFermiDiracdT
  PUBLIC :: dFermiDiracdY

  REAL(DP), PARAMETER :: Log1d100 = LOG( 1.0d100 )
  REAL(DP), PARAMETER :: UnitD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitT    = Kelvin
  REAL(DP), PARAMETER :: UnitY    = One
  REAL(DP), PARAMETER :: UnitE    = MeV
  REAL(DP), PARAMETER :: UnitEta  = One
  REAL(DP), PARAMETER :: UnitEC   = One / Centimeter
  REAL(DP), PARAMETER :: UnitES   = One / ( Centimeter * MeV**2 )
  REAL(DP), PARAMETER :: UnitNES  = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitPair = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitBrem = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: f0_Max   = One - EPSILON( One )

  REAL(DP), PARAMETER :: hbarMeVs = 6.582119569d-22
  REAL(DP), PARAMETER :: C_A      = -1.26d0/2.0d0 ! C_A from HR98
  REAL(DP), PARAMETER :: G_F      = 1.166d-11 ! G_F from HR98 in [MeV**-2]
  REAL(DP), PARAMETER :: nb_14    = AvogadroConstantMKS * 1d14
  REAL(DP), PARAMETER :: Brem_const = C_A**2 * G_F**2 * nb_14 &
                                    * hbarMeVs**2 * SpeedOfLightCGS**2 &
                                    / TwoPi**3 * UnitBrem

  REAL(dp), PARAMETER :: Alpha_Brem(3) = [ 1.0d0, 1.0d0, 28.d0/3.d0 ]

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
    REAL(DP), INTENT(in)  :: E(iE_B:)
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    REAL(DP), INTENT(out) :: f0(iE_B:,iS_B:,iX_B:)

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

      f0(iE,iS,iX) = FermiDirac( E(iE), Mnu, kT )

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


  SUBROUTINE ComputeEquilibriumDistributions_DG &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0 )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in) , TARGET, CONTIGUOUS :: E(iE_B:)
    REAL(DP), INTENT(in) , TARGET, CONTIGUOUS :: D(iX_B:)
    REAL(DP), INTENT(in) , TARGET, CONTIGUOUS :: T(iX_B:)
    REAL(DP), INTENT(in) , TARGET, CONTIGUOUS :: Y(iX_B:)
    REAL(DP), INTENT(out), TARGET, CONTIGUOUS :: f0(iE_B:,iS_B:,iX_B:)

    REAL(DP), POINTER :: E_Q(:,:), f0_Q(:,:,:,:)

    REAL(DP) :: f0_K(          1:(iE_E-iE_B+1)/nDOFE+1,1:(iS_E-iS_B+1),1:(iX_E-iX_B+1))
    REAL(DP) :: f0_P(1:nDOFE+2,1:(iE_E-iE_B+1)/nDOFE  ,1:(iS_E-iS_B+1),1:(iX_E-iX_B+1))
    REAL(DP) :: N_K, V_K, f0_Min, Min_K, Max_K, Theta
    INTEGER  :: iE, iS, iX, iNodeE, nE, nS, nX

    nE = ( iE_E - iE_B + 1 ) / nDOFE
    nS = ( iS_E - iS_B + 1 )
    nX = ( iX_E - iX_B + 1 )

    ! --- Permute Data ---

    E_Q (1:nDOFE,1:nE)           => E (:)
    f0_Q(1:nDOFE,1:nE,1:nS,1:nX) => f0(:,:,:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: f0, f0_K, f0_P ) &
    !$OMP MAP( to: E, D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( f0, f0_K, f0_P ) &
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

    ! --- Estimate Cell Average in Outer Ghost Element ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX = 1, nX
    DO iS = 1, nS

      f0_K(nE+1,iS,iX) = f0_K(nE,iS,iX)**2 / f0_K(nE-1,iS,iX)

    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFE+2, nE*nS*nX, nDOFE, One, InterpMat_E, nDOFE+2, &
             f0_Q, nDOFE, Zero, f0_P, nDOFE+2 )

    ! --- Limit Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( f0_Min, Min_K, Max_K, Theta )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta )
#endif
    DO iX = 1, nX
    DO iS = 1, nS
    DO iE = 1, nE

      f0_Min = f0_K(iE+1,iS,iX)

      Max_K = f0_Max
      Min_K = f0_Min
      DO iNodeE = 1, nDOFE+2
        Max_K = MAX( Max_K, f0_P(iNodeE,iE,iS,iX) )
        Min_K = MIN( Min_K, f0_P(iNodeE,iE,iS,iX) )
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
    !$OMP MAP( release: f0_K, f0_P, E, D, T, Y ) &
    !$OMP MAP( from: f0 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( f0_K, f0_P, E, D, T, Y ) &
    !$ACC COPYOUT( f0 )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions_DG


  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, f0, df0dY, df0dU )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:)
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    REAL(DP), INTENT(out) :: f0   (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(out) :: df0dY(iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(out) :: df0dU(iE_B:,iS_B:,iX_B:)

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
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, E, D, T, Y, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:)
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    REAL(DP), INTENT(out) :: opEC(iE_B:,iS_B:,iX_B:)

    REAL(DP) :: LogE_P(iE_B:iE_E)
    REAL(DP) :: LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    INTEGER  :: iE, iS, iX

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
    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P   (iX) =        Y(iX) / UnitY
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE = iE_B, iE_E

      IF ( iS > iNuE_Bar ) THEN

        opEC(iE,iS,iX) = Zero

      ELSE

        CALL LogInterpolateSingleVariable_4D_Custom_Point &
               ( LogE_P(iE), LogD_P(iX), LogT_P(iX), Y_P(iX), &
                 LogEs_T   , LogDs_T   , LogTs_T   , Ys_T   , &
                 OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(iE,iS,iX) )

        opEC(iE,iS,iX) = opEC(iE,iS,iX) * UnitEC

      END IF

    END DO
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
    REAL(DP), INTENT(in)  :: E(iP_B:)
    REAL(DP), INTENT(in)  :: D(iP_B:)
    REAL(DP), INTENT(in)  :: T(iP_B:)
    REAL(DP), INTENT(in)  :: Y(iP_B:)
    REAL(DP), INTENT(out) :: opEC(iP_B:,iS_B:)

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
    REAL(DP), INTENT(in)  :: E(iE_B:)
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(iE_B:,iX_B:)

    REAL(DP) :: LogE_P(iE_B:iE_E)
    REAL(DP) :: LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    INTEGER  :: iX, iE

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

    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(1,iMoment), Iso_T(:,:,:,:,iMoment,1), opES )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E
      opES(iE,iX) = opES(iE,iX) * UnitES
    END DO
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
    REAL(DP), INTENT(in)  :: E(iP_B:)
    REAL(DP), INTENT(in)  :: D(iP_B:)
    REAL(DP), INTENT(in)  :: T(iP_B:)
    REAL(DP), INTENT(in)  :: Y(iP_B:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(iP_B:)

    REAL(DP) :: LogE_P(iP_B:iP_E)
    REAL(DP) :: LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E)
    INTEGER  :: iP, iS

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
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: H_I (iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(out) :: H_II(iE_B:,iE_B:,iX_B:)

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
    REAL(DP), INTENT(in)  :: W2  (iE_B:)
    REAL(DP), INTENT(in)  :: J   (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(in)  :: J0  (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(in)  :: H_I (iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(in)  :: H_II(iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(out) :: Eta (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(out) :: Chi (iE_B:,iS_B:,iX_B:)

    REAL(DP) :: DetBal, Phi_Out, Phi_In
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iE, iE1, iE2, iS, iX

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
          Phi_Out = ( C1(iS) * H_I(iE1,iE2,iX) + C2(iS) * H_II(iE1,iE2,iX) ) * UnitNES
          Phi_In  = Phi_Out * DetBal
        ELSE
          Phi_In  = ( C1(iS) * H_I(iE2,iE1,iX) + C2(iS) * H_II(iE2,iE1,iX) ) * UnitNES
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


  SUBROUTINE ComputeNeutrinoOpacities_Pair &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, iMoment, J_I, J_II )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: J_I (iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(out) :: J_II(iE_B:,iE_B:,iX_B:)

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
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
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
    DO iX = iX_B, iX_E
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
    REAL(DP), INTENT(in)  :: W2  (iE_B:)
    REAL(DP), INTENT(in)  :: J   (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(in)  :: J0  (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(in)  :: J_I (iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(in)  :: J_II(iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(out) :: Eta (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(out) :: Chi (iE_B:,iS_B:,iX_B:)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, iS, iS_A

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
          Phi_0_Ann = ( C1(iS) * J_I (iE1,iE2,iX) + C2(iS) * J_II(iE1,iE2,iX) ) * UnitPair
        ELSE
          Phi_0_Ann = ( C1(iS) * J_II(iE2,iE1,iX) + C2(iS) * J_I (iE2,iE1,iX) ) * UnitPair
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


  SUBROUTINE ComputeNeutrinoOpacities_Brem &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, S_Sigma )

    ! --- Brem Opacities (Multiple D,T) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(iX_B:)
    REAL(DP), INTENT(in)  :: T(iX_B:)
    REAL(DP), INTENT(in)  :: Y(iX_B:)
    REAL(DP), INTENT(out) :: S_Sigma(iE_B:,iE_B:,iX_B:)

    INTEGER  :: iX, iE1, iE2
    REAL(DP) :: Xp(iX_B:iX_E), Xn(iX_B:iX_E) !Proton and neutron mass fractions
    REAL(DP) :: LogT_P(iX_B:iX_E)
    REAL(DP) :: LogDX_P(3,iX_B:iX_E)

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( alloc: Xp, Xn, LogT_P, LogDX_P, S_Sigma ) &
    !$OMP MAP( to: D, T, Y )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC CREATE( Xp, Xn, LogT_P, LogDX_P, S_Sigma ) &
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

      LogDX_P(1,iX) = LOG10(D(iX) * Xp(iX) / UnitD)
      LogDX_P(2,iX) = LOG10(D(iX) * Xn(iX) / UnitD)
      LogDX_P(3,iX) = LOG10(D(iX) * SQRT(ABS(Xp(iX)*Xn(iX))) / UnitD)

      LogT_P(iX)   = LOG10(T(iX) / UnitT)    

    END DO

    CALL SumLogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogDX_P, LogT_P, LogDs_T, LogTs_T, Alpha_Brem, &
             OS_Brem(1,1), Brem_AT(:,:,:,:,1,1), S_Sigma )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Xp, Xn, LogT_P, LogDX_P, D, T, Y ) &
    !$OMP MAP( from: S_Sigma )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( Xp, Xn, LogT_P, LogDX_P, D, T, Y ) &
    !$ACC COPYOUT( S_Sigma )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP MAP( from: S_Sigma )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC COPYOUT( S_Sigma )
#endif
    DO iX = iX_B, iX_E
    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      S_Sigma(iE1,iE2,iX) = Zero
    END DO
    END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Brem

  
  SUBROUTINE ComputeNeutrinoOpacityRates_Brem &
    ( iE_B, iE_E, iS_B, iS_E, iX_B, iX_E, W2, J, J0, S_Sigma, Eta, Chi )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (iE_B:)
    REAL(DP), INTENT(in)  :: J      (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(in)  :: J0     (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(in)  :: S_Sigma(iE_B:,iE_B:,iX_B:)
    REAL(DP), INTENT(out) :: Eta    (iE_B:,iS_B:,iX_B:)
    REAL(DP), INTENT(out) :: Chi    (iE_B:,iS_B:,iX_B:)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, iS, iS_A

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro ) &
    !$OMP MAP( to: S_Sigma, W2, J, J0 ) &
    !$OMP MAP( from: Eta, Chi )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro ) &
    !$ACC COPYIN( S_Sigma, W2, J, J0 ) &
    !$ACC COPYOUT( Eta, Chi )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro )
#endif
    DO iX = iX_B, iX_E
    DO iS = iS_B, iS_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iS,iX) * J0(iE1,iS_A,iX) ) &
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
