#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_OPACITY
#endif

MODULE NeutrinoOpacitiesComputationModule

  USE KindModule, ONLY: &
    DP, Zero, One, SqrtTiny
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
  USE DeviceModule, ONLY: &
    QueryOnGpu
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
    iNuE, iNuE_Bar

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules ---

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_4D_Custom, &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned, &
    SumLogInterpolateSingleVariable_2D2D_Custom_Aligned

  ! ----------------------------------------------

#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

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
                                    * hbarMeVs**2 * SpeedOfLightCGS**2 * UnitBrem

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


  SUBROUTINE ComputeEquilibriumDistributions_Point &
    ( E, D, T, Y, f0, iSpecies )
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

    REAL(DP) :: Me, Mp, Mn, Mnu, kT

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    kT = BoltzmannConstant * T

    IF ( iSpecies == iNuE ) THEN
      Mnu = ( Me + Mp ) - Mn
    ELSE IF ( iSpecies == iNuE_Bar ) THEN
      Mnu = Mn - ( Me + Mp )
    ELSE
      Mnu = Zero
    END IF

    f0 = FermiDirac( E, Mnu, kT )

  END SUBROUTINE ComputeEquilibriumDistributions_Point


  SUBROUTINE ComputeEquilibriumDistributions &
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, E, D, T, Y, f0 )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: f0(:,:,:)

    INTEGER  :: iX, iE, iS
    REAL(DP) :: Me(iX_B:iX_E), Mp(iX_B:iX_E), Mn(iX_B:iX_E)
    REAL(DP) :: Mnu, kT
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( f0 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeEquilibriumDistributions] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions]   Y missing'
      IF ( .not. QueryOnGPU( f0 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions]   f0 missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: Me, Mp, Mn )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me, Mp, Mn )
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
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( Mnu, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( Mnu, kT ) &
    !$ACC PRESENT( Me, Mp, Mn, E, T, f0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Mnu, kT )
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E

      kT = BoltzmannConstant * T(iX)

      IF ( iS == iNuE ) THEN
        Mnu = ( Me(iX) + Mp(iX) ) - Mn(iX)
      ELSE IF ( iS == iNuE_Bar ) THEN
        Mnu = Mn(iX) - ( Me(iX) + Mp(iX) )
      ELSE
        Mnu = Zero
      END IF

      f0(iE,iX,iS) = FermiDirac( E(iE), Mnu, kT )

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me, Mp, Mn )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me, Mp, Mn )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions


  SUBROUTINE ComputeEquilibriumDistributions_DG &
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, E, D, T, Y, f0 )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in) , TARGET, CONTIGUOUS :: E(:)
    REAL(DP), INTENT(in) , TARGET :: D(:)
    REAL(DP), INTENT(in) , TARGET :: T(:)
    REAL(DP), INTENT(in) , TARGET :: Y(:)
    REAL(DP), INTENT(out), TARGET, CONTIGUOUS :: f0(:,:,:)

    INTEGER  :: iX, iE, iS, iNodeE, nE, nX, nS
    REAL(DP) :: N_K, V_K, f0_Min, Min_K, Max_K, Theta
    REAL(DP), POINTER :: E_Q(:,:), f0_Q(:,:,:,:)
    REAL(DP) :: f0_K(          1:(iE_E-iE_B+1)/nDOFE+1,1:(iX_E-iX_B+1),iS_B:iS_E)
    REAL(DP) :: f0_P(1:nDOFE+2,1:(iE_E-iE_B+1)/nDOFE  ,1:(iX_E-iX_B+1),iS_B:iS_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( f0 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeEquilibriumDistributions_DG] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_DG]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_DG]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_DG]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_DG]   Y missing'
      IF ( .not. QueryOnGPU( f0 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_DG]   f0 missing'
    END IF
#endif

    CALL ComputeEquilibriumDistributions &
           ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, E, D, T, Y, f0 )

    nE = ( iE_E - iE_B + 1 ) / nDOFE
    nX = ( iX_E - iX_B + 1 )
    nS = ( iS_E - iS_B + 1 )

    ! --- Permute Data ---

    E_Q(1:nDOFE,1:nE) => E(:)
    f0_Q(1:nDOFE,1:nE,1:nX,iS_B:iS_E) => f0(:,:,:)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: E_Q, f0_K, f0_Q, f0_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( E_Q, f0_K, f0_Q, f0_P )
#endif

    ! --- Cell Average of Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( V_K, N_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( V_K, N_K ) &
    !$ACC PRESENT( WeightsE, E_Q, f0_K, f0_Q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( V_K, N_K )
#endif
    DO iS = iS_B, iS_E
    DO iX = 1, nX
    DO iE = 1, nE

      V_K = Zero
      N_K = Zero

      DO iNodeE = 1, nDOFE
        V_K = V_K + WeightsE(iNodeE) * E_Q(iNodeE,iE)**2
        N_K = N_K + WeightsE(iNodeE) * E_Q(iNodeE,iE)**2 * f0_Q(iNodeE,iE,iX,iS)
      END DO

      f0_K(iE,iX,iS) = N_K / V_K

      IF( f0_K(iE,iX,iS) > f0_Max )THEN
        f0_K(iE,iX,iS) = f0_Max
        DO iNodeE = 1, nDOFE
          f0_Q(iNodeE,iE,iX,iS) = f0_Max
        END DO
      END IF

    END DO
    END DO
    END DO

    ! --- Estimate Cell Average in Outer Ghost Element ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( f0_K )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS = iS_B, iS_E
    DO iX = 1, nX

      f0_K(nE+1,iX,iS) = f0_K(nE,iX,iS)**2 / f0_K(nE-1,iX,iS)

    END DO
    END DO

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFE+2, nX*nE*nS, nDOFE, One, InterpMat_E, nDOFE+2, &
             f0_Q, nDOFE, Zero, f0_P, nDOFE+2 )

    ! --- Limit Equilibrium Distributions ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( f0_Min, Min_K, Max_K, Theta ) &
    !$ACC PRESENT( f0_K, f0_Q, f0_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( f0_Min, Min_K, Max_K, Theta )
#endif
    DO iS = iS_B, iS_E
    DO iX = 1, nX
    DO iE = 1, nE

      f0_Min = f0_K(iE+1,iX,iS)

      Max_K = f0_Max
      Min_K = f0_Min
      DO iNodeE = 1, nDOFE+2
        Max_K = MAX( Max_K, f0_P(iNodeE,iE,iX,iS) )
        Min_K = MIN( Min_K, f0_P(iNodeE,iE,iX,iS) )
      END DO

      IF( Min_K < f0_Min .OR. Max_K > f0_Max )THEN

        Theta &
          = MIN( One, &
                 ABS((f0_Min-f0_K(iE,iX,iS))/(Min_K-f0_K(iE,iX,iS)+SqrtTiny)), &
                 ABS((f0_Max-f0_K(iE,iX,iS))/(Max_K-f0_K(iE,iX,iS)+SqrtTiny)) )

        DO iNodeE = 1, nDOFE
          f0_Q(iNodeE,iE,iX,iS) &
            = ( One - Theta ) * f0_K(iE,iX,iS) &
              +       Theta   * f0_Q(iNodeE,iE,iX,iS)
        END DO

      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: E_Q, f0_K, f0_Q, f0_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( E_Q, f0_K, f0_Q, f0_P )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions_DG


  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives &
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, E, D, T, Y, f0, df0dY, df0dU )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: f0   (:,:,:)
    REAL(DP), INTENT(out) :: df0dY(:,:,:)
    REAL(DP), INTENT(out) :: df0dU(:,:,:)

    REAL(DP), DIMENSION(iX_B:iX_E) :: Me,  dMedT , dMedY
    REAL(DP), DIMENSION(iX_B:iX_E) :: Mp,  dMpdT , dMpdY
    REAL(DP), DIMENSION(iX_B:iX_E) :: Mn,  dMndT , dMndY
    REAL(DP), DIMENSION(iX_B:iX_E) :: U,   dUdT,   dUdY, dUdD
    REAL(DP)                       :: Mnu, dMnudT, dMnudY

    REAL(DP) :: kT, df0dT_Y, df0dY_T
    INTEGER  :: iX, iE, iS
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( f0, df0dY, df0dU )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   Y missing'
      IF ( .not. QueryOnGPU( f0 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   f0 missing'
      IF ( .not. QueryOnGPU( df0dY ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   df0dY missing'
      IF ( .not. QueryOnGPU( df0dU ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives]   df0dU missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: Me,  dMedT,  dMedY, &
    !$OMP             Mp,  dMpdT,  dMpdY, &
    !$OMP             Mn,  dMndT,  dMndY, &
    !$OMP             U,   dUdT,   dUdY, dUdD )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me,  dMedT,  dMedY, &
    !$ACC         Mp,  dMpdT,  dMpdY, &
    !$ACC         Mn,  dMndT,  dMndY, &
    !$ACC         U,   dUdT,   dUdY, dUdD )
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
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT, Mnu, dMnudT, dMnudY, df0dT_Y, df0dY_T )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT, Mnu, dMnudT, dMnudY, df0dT_Y, df0dY_T ) &
    !$ACC PRESENT( f0, df0dY, df0dU, E, T, dUdT, dUdY )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( kT, Mnu, dMnudT, dMnudY, df0dT_Y, df0dY_T )
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E

      kT = BoltzmannConstant * T(iX)

      IF ( iS == iNuE ) THEN
        Mnu    = ( Me   (iX) + Mp   (iX) ) - Mn   (iX)
        dMnudT = ( dMedT(iX) + dMpdT(iX) ) - dMndT(iX)
        dMnudY = ( dMedY(iX) + dMpdY(iX) ) - dMndY(iX)
      ELSE IF ( iS == iNuE_Bar ) THEN
        Mnu    = Mn   (iX) - ( Me   (iX) + Mp   (iX) )
        dMnudT = dMndT(iX) - ( dMedT(iX) + dMpdT(iX) )
        dMnudY = dMndY(iX) - ( dMedY(iX) + dMpdY(iX) )
      ELSE
        Mnu    = Zero
        dMnudT = Zero
        dMnudY = Zero
      END IF

      f0(iE,iX,iS) = FermiDirac   ( E(iE), Mnu, kT )
      df0dT_Y      = dFermiDiracdT( E(iE), Mnu, kT, dMnudT, T(iX) ) ! Constant T
      df0dY_T      = dFermiDiracdY( E(iE), Mnu, kT, dMnudY, T(iX) ) ! Constant Y

      df0dU(iE,iX,iS) = df0dT_Y / dUdT(iX)
      df0dY(iE,iX,iS) = df0dY_T - df0dU(iE,iX,iS) * dUdY(iX)

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me,  dMedT,  dMedY, &
    !$OMP               Mp,  dMpdT,  dMpdY, &
    !$OMP               Mn,  dMndT,  dMndY, &
    !$OMP               U,   dUdT,   dUdY, dUdD )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me,  dMedT,  dMedY, &
    !$ACC         Mp,  dMpdT,  dMpdY, &
    !$ACC         Mn,  dMndT,  dMndY, &
    !$ACC         U,   dUdT,   dUdY, dUdD )
#endif

  END SUBROUTINE ComputeEquilibriumDistributionAndDerivatives


  SUBROUTINE ComputeNeutrinoOpacities_EC &
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, E, D, T, Y, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: opEC(:,:,:)

    INTEGER  :: iX, iE, iS
    REAL(DP) :: LogE_P(iE_B:iE_E), LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opEC )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_EC] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC]   Y missing'
      IF ( .not. QueryOnGPU( opEC ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC]   opEC missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, LogE_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    DO iS = iS_B, iS_E

      CALL LogInterpolateSingleVariable_1D3D_Custom &
             ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(:,:,iS), &
               GPU_Option = do_gpu )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opEC )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E
      opEC(iE,iX,iS) = opEC(iE,iX,iS) * UnitEC
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opEC = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC

  
  SUBROUTINE ComputeNeutrinoOpacities_EC_Vector &
    ( iP_B, iP_E, iS_B, iS_E, E, D, T, Y, opEC )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: opEC(:,:)

    INTEGER  :: iP, iS
    REAL(DP) :: LogE_P(iP_B:iP_E), LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opEC )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Vector] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Vector]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Vector]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Vector]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Vector]   Y missing'
      IF ( .not. QueryOnGPU( opEC ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Vector]   opEC missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P, E, LogE_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Y_P(iP) = Y(iP) / UnitY
      LogE_P(iP) = LOG10( E(iP) / UnitE )
    END DO

    DO iS = iS_B, iS_E

      CALL LogInterpolateSingleVariable_4D_Custom &
             ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OS_EmAb(iS), EmAb_T(:,:,:,:,iS), opEC(:,iS), &
               GPU_Option = do_gpu )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opEC )
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
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opEC = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Vector


  SUBROUTINE ComputeNeutrinoOpacities_ES &
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, E, D, T, Y, iMoment, opES )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(:,:,:)

    INTEGER  :: iX, iE, iS
    REAL(DP) :: LogE_P(iE_B:iE_E), LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opES )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_ES] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES]   Y missing'
      IF ( .not. QueryOnGPU( opES ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES]   opES missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, LogE_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    DO iS = iS_B, iS_E

      CALL LogInterpolateSingleVariable_1D3D_Custom &
             ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OS_Iso(iS,iMoment), Iso_T(:,:,:,:,iMoment,iS), opES(:,:,iS), &
               GPU_Option = do_gpu )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opES )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(3)
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE = iE_B, iE_E
      opES(iE,iX,iS) = opES(iE,iX,iS) * UnitES
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opES = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES

  
  SUBROUTINE ComputeNeutrinoOpacities_ES_Vector &
    ( iP_B, iP_E, iS_B, iS_E, E, D, T, Y, iMoment, opES )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES(:,:)

    INTEGER  :: iP, iS
    REAL(DP) :: LogE_P(iP_B:iP_E), LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opES )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Vector] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Vector]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Vector]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Vector]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Vector]   Y missing'
      IF ( .not. QueryOnGPU( opES ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Vector]   opES missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P, E, LogE_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Y_P(iP) = Y(iP) / UnitY
      LogE_P(iP) = LOG10( E(iP) / UnitE )
    END DO

    DO iS = iS_B, iS_E

      CALL LogInterpolateSingleVariable_4D_Custom &
             ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OS_Iso(iS,iMoment), Iso_T(:,:,:,:,iMoment,iS), opES(:,iS), &
               GPU_Option = do_gpu )

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opES )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(2)
#endif
    DO iS = iS_B, iS_E
    DO iP = iP_B, iP_E
      opES(iP,iS) = opES(iP,iS) * UnitES
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opES = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Vector


  SUBROUTINE ComputeNeutrinoOpacities_NES &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, iMoment, H_I, H_II )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: H_I (:,:,:)
    REAL(DP), INTENT(out) :: H_II(:,:,:)

    INTEGER  :: iX, iE1, iE2, iH_I, iH_II
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( D, T, Y ) &
       .AND. QueryOnGPU( H_I, H_II )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) &
        '[ComputeNeutrinoOpacities_NES] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES]   Y missing'
      IF ( .not. QueryOnGPU( H_I ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES]   H_I missing'
      IF ( .not. QueryOnGPU( H_II ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES]   H_II missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

    iH_I  = ( iMoment - 1 ) * 2 + 1
    iH_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogT_P, LogEta_P )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
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
             OS_NES(1,iH_I), NES_AT(:,:,:,:,iH_I,1), H_I, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

    ! --- Interpolate HII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH_II), NES_AT(:,:,:,:,iH_II,1), H_II, &
             GPU_Option = do_gpu, ASYNC_Option = 1  )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogT_P, LogEta_P )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( H_I, H_II )
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
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, W2, J, J0, H_I, H_II, Eta, Chi )

    ! --- Neutrino-Electron Scattering Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: W2  (:)
    REAL(DP), INTENT(in)  :: J   (:,:,:)
    REAL(DP), INTENT(in)  :: J0  (:,:,:)
    REAL(DP), INTENT(in)  :: H_I (:,:,:)
    REAL(DP), INTENT(in)  :: H_II(:,:,:)
    REAL(DP), INTENT(out) :: Eta (:,:,:)
    REAL(DP), INTENT(out) :: Chi (:,:,:)

    REAL(DP) :: DetBal, Phi_Out, Phi_In
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, iS
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( W2 ) &
       .AND. QueryOnGPU( J, J0, Eta, Chi, H_I, H_II )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES] Data not present on device'
      IF ( .not. QueryOnGPU( W2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   W2 missing'
      IF ( .not. QueryOnGPU( J ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   J missing'
      IF ( .not. QueryOnGPU( J0 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   J0 missing'
      IF ( .not. QueryOnGPU( Eta ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   Eta missing'
      IF ( .not. QueryOnGPU( Chi ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   Chi missing'
      IF ( .not. QueryOnGPU( H_I ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   H_I missing'
      IF ( .not. QueryOnGPU( H_II ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_NES]   H_II missing'
    END IF
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, DetBal, Phi_In, Phi_Out )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( SUM1, SUM2, DetBal, Phi_In, Phi_Out ) &
    !$ACC PRESENT( H_I, H_II, Eta, Chi, W2, J, J0, C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( SUM1, SUM2, DetBal, Phi_In, Phi_Out )
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE2 = iE_B, iE_E

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iX,iS) * ( One - J0(iE1,iX,iS) ) ) &
                 / ( J0(iE1,iX,iS) * ( One - J0(iE2,iX,iS) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_Out = ( C1(iS) * H_I(iE1,iE2,iX) + C2(iS) * H_II(iE1,iE2,iX) ) * UnitNES
          Phi_In  = Phi_Out * DetBal
        ELSE
          Phi_In  = ( C1(iS) * H_I(iE2,iE1,iX) + C2(iS) * H_II(iE2,iE1,iX) ) * UnitNES
          Phi_Out = Phi_In / DetBal
        END IF

        SUM1 = SUM1 + Phi_In  * W2(iE1) * J(iE1,iX,iS)
        SUM2 = SUM2 + Phi_Out * W2(iE1) * ( One - J(iE1,iX,iS) )

      END DO

      Eta(iE2,iX,iS) = SUM1
      Chi(iE2,iX,iS) = SUM1 + SUM2

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_NES


  SUBROUTINE ComputeNeutrinoOpacities_Pair &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, iMoment, J_I, J_II )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: J_I (:,:,:)
    REAL(DP), INTENT(out) :: J_II(:,:,:)

    INTEGER  :: iX, iE1, iE2, iJ_I, iJ_II
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( D, T, Y ) &
       .AND. QueryOnGPU( J_I, J_II )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_Pair] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair]   Y missing'
      IF ( .not. QueryOnGPU( J_I ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair]   J_I missing'
      IF ( .not. QueryOnGPU( J_II ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair]   J_II missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

    iJ_I  = ( iMoment - 1 ) * 2 + 1
    iJ_II = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogT_P, LogEta_P )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
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
             OS_Pair(1,iJ_I), Pair_AT(:,:,:,:,iJ_I,1), J_I, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

    ! --- Interpolate JII ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ_II), Pair_AT(:,:,:,:,iJ_II,1), J_II, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogT_P, LogEta_P )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( J_I, J_II)
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
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, W2, J, J0, J_I, J_II, Eta, Chi )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: W2 (:)
    REAL(DP), INTENT(in)  :: J   (:,:,:)
    REAL(DP), INTENT(in)  :: J0  (:,:,:)
    REAL(DP), INTENT(in)  :: J_I (:,:,:)
    REAL(DP), INTENT(in)  :: J_II(:,:,:)
    REAL(DP), INTENT(out) :: Eta (:,:,:)
    REAL(DP), INTENT(out) :: Chi (:,:,:)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, iS, iS_A
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( W2 ) &
       .AND. QueryOnGPU( J, J0, Eta, Chi, J_I, J_II )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair] Data not present on device'
      IF ( .not. QueryOnGPU( W2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair]   W2 missing'
      IF ( .not. QueryOnGPU( J ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair]   J missing'
      IF ( .not. QueryOnGPU( Eta ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair]   Eta missing'
      IF ( .not. QueryOnGPU( Chi ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair]   Chi missing'
      IF ( .not. QueryOnGPU( Phi_0_Pro ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair]   Phi_0_Pro missing'
      IF ( .not. QueryOnGPU( Phi_0_Ann ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Pair]   Phi_0_Ann missing'
    END IF
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Pro, Phi_0_Ann )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Pro, Phi_0_Ann ) &
    !$ACC PRESENT( J_I, J_II, Eta, Chi, W2, J, J0, C1, C2 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Pro, Phi_0_Ann )
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iX,iS) * J0(iE1,iX,iS_A) ) &
                 / ( ( One - J0(iE2,iX,iS) ) * ( One - J0(iE1,iX,iS_A) ) )

        IF ( iE1 <= iE2 ) THEN
          Phi_0_Ann = ( C1(iS) * J_I (iE1,iE2,iX) + C2(iS) * J_II(iE1,iE2,iX) ) * UnitPair
        ELSE
          Phi_0_Ann = ( C1(iS) * J_II(iE2,iE1,iX) + C2(iS) * J_I (iE2,iE1,iX) ) * UnitPair
        END IF
        Phi_0_Pro = Phi_0_Ann * DetBal

        SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iX,iS) )
        SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iX,iS)

      END DO

      Eta(iE2,iX,iS) = SUM1
      Chi(iE2,iX,iS) = SUM1 + SUM2

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeNeutrinoOpacityRates_Pair


  SUBROUTINE ComputeNeutrinoOpacities_Brem &
    ( iE_B, iE_E, iX_B, iX_E, D, T, Y, S_Sigma )

    ! --- Brem Opacities (Multiple D,T) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: S_Sigma(:,:,:)

    INTEGER  :: iX, iE1, iE2
    REAL(DP) :: LogT_P(iX_B:iX_E)
    REAL(DP) :: LogD_X(3,iX_B:iX_E)
    REAL(DP) :: Xp(iX_B:iX_E), Xn(iX_B:iX_E) !Proton and neutron mass fractions
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( D, T, Y ) &
       .AND. QueryOnGPU( S_Sigma )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_Brem] Data not present on device'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Brem]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Brem]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Brem]   Y missing'
      IF ( .not. QueryOnGPU( S_Sigma ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Brem]   S_Sigma missing'
    END IF
#endif


#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: Xp, Xn, LogT_P, LogD_X )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Xp, Xn, LogT_P, LogD_X )
#endif

    ! --- Compute proton and neutron fractions ---

    CALL ComputeProtonMassFraction_TABLE &
           ( D, T, Y, Xp )

    CALL ComputeNeutronMassFraction_TABLE &
           ( D, T, Y, Xn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, T, Xp, Xn, LogT_P, LogD_X )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iX = iX_B, iX_E

      LogD_X(1,iX) = LOG10(D(iX) * Xp(iX) / UnitD)
      LogD_X(2,iX) = LOG10(D(iX) * Xn(iX) / UnitD)
      LogD_X(3,iX) = LOG10(D(iX) * SQRT(ABS(Xp(iX)*Xn(iX))) / UnitD)

      LogT_P(iX)   = LOG10(T(iX) / UnitT)    

    END DO

    CALL SumLogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogD_X, LogT_P, LogDs_T, LogTs_T, Alpha_Brem, &
             OS_Brem(1,1), Brem_AT(:,:,:,:,1,1), S_Sigma, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Xp, Xn, LogT_P, LogD_X )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Xp, Xn, LogT_P, LogD_X )

    !$ACC WAIT(1)
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( S_Sigma )
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
    ( iE_B, iE_E, iX_B, iX_E, iS_B, iS_E, W2, J, J0, S_Sigma, Eta, Chi )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    INTEGER,  INTENT(in)  :: iS_B, iS_E
    REAL(DP), INTENT(in)  :: W2     (:)
    REAL(DP), INTENT(in)  :: J      (:,:,:)
    REAL(DP), INTENT(in)  :: J0     (:,:,:)
    REAL(DP), INTENT(in)  :: S_Sigma(:,:,:)
    REAL(DP), INTENT(out) :: Eta    (:,:,:)
    REAL(DP), INTENT(out) :: Chi    (:,:,:)

    REAL(DP) :: DetBal, Phi_0_Ann, Phi_0_Pro
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, iS, iS_A
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( W2 ) &
       .AND. QueryOnGPU( J, J0, Eta, Chi, S_Sigma )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacityRates_Brem] Data not present on device'
      IF ( .not. QueryOnGPU( W2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Brem]   W2 missing'
      IF ( .not. QueryOnGPU( J ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Brem]   J missing'
      IF ( .not. QueryOnGPU( Eta ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Brem]   Eta missing'
      IF ( .not. QueryOnGPU( Chi ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Brem]   Chi missing'
      IF ( .not. QueryOnGPU( S_Sigma ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacityRates_Brem]   S_Sigma missing'
    END IF
#endif

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro ) &
    !$ACC PRESENT( S_Sigma, Eta, Chi, W2, J )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iS_A, SUM1, SUM2, DetBal, Phi_0_Ann, Phi_0_Pro )
#endif
    DO iS = iS_B, iS_E
    DO iX = iX_B, iX_E
    DO iE2 = iE_B, iE_E

      ! Get index for corresponding anti-neutrino
      iS_A = iS + 2*MOD(iS,2) - 1

      SUM1 = Zero
      SUM2 = Zero

      DO iE1 = iE_B, iE_E

        DetBal =   ( J0(iE2,iX,iS) * J0(iE1,iX,iS_A) ) &
                 / ( ( One - J0(iE2,iX,iS) ) * ( One - J0(iE1,iX,iS_A) ) )

        Phi_0_Ann = S_Sigma(iE1,iE2,iX) * 3.0d0 * Brem_const
        Phi_0_Pro = Phi_0_Ann * DetBal

        SUM1 = SUM1 + Phi_0_Pro * W2(iE1) * ( One - J(iE1,iX,iS) )
        SUM2 = SUM2 + Phi_0_Ann * W2(iE1) * J(iE1,iX,iS)

      END DO

      Eta(iE2,iX,iS) = SUM1
      Chi(iE2,iX,iS) = SUM1 + SUM2

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
